import numpy as np
import awkward as ak
import correctionlib._core as correction
from dataclasses import dataclass
from typing import Optional

import JecConfigReader
import ColumnarHelper
import PrintHelper
import CollectJetMet


# ============================================================
# Small helpers & context
# ============================================================

def phi_mpi_pi(dphi):
    """Wrap delta-phi into [-pi, pi]."""
    return np.arctan2(np.sin(dphi), np.cos(dphi))


def representative_run_number(year: str, runs):
    """Use a representative run for run-based residuals (as in C++)."""
    if year == "2016Pre":
        return ak.ones_like(runs, dtype="float64") * 272007
    if year == "2016Post":
        return ak.ones_like(runs, dtype="float64") * 278769
    if year == "2017":
        return ak.ones_like(runs, dtype="float64") * 297020
    if year == "2018":
        return ak.ones_like(runs, dtype="float64") * 315252
    if year == "2022Pre":
        return ak.ones_like(runs, dtype="float64") * 355065.0
    if year == "2022Post":
        return ak.ones_like(runs, dtype="float64") * 359022.0
    if year == "2023Pre":
        return ak.ones_like(runs, dtype="float64") * 367080.0
    if year == "2023Post":
        return ak.ones_like(runs, dtype="float64") * 369803.0
    if year == "2024":
        return ak.ones_like(runs, dtype="float64") * 379412.0
    if year == "2025":
        return ak.ones_like(runs, dtype="float64") * 392159.0
    if year == "2026":
        return ak.ones_like(runs, dtype="float64") * 401837.0
    return ak.values_astype(runs, "float64")


def has_phi_dependent_L2(year: str) -> bool:
    return year in ("2023Post", "2024", "2025", "2026")


def requires_run_based_residual(year: str) -> bool:
    return year in ("2016Pre", "2016Post", "2017", "2018", "2022Pre", "2022Post","2023Pre", "2023Post", "2024", "2025", "2026")


@dataclass
class SystematicOptions:
    """
    Python analogue of the C++ SystematicOptions:

      struct SystematicOptions {
        std::string jerVar;          // "nom" | "up" | "down"
        std::optional<JerBin> jerRegion;
        std::optional<Correction::Ref> jesSystRef;
        std::string jesSystVar;      // "Up" | "Down"
      };

    Here we store names/tags instead of Correction::Ref directly.
    """
    syst_kind: str = "Nominal"  # "Nominal", "JES", "JER"
    jerVar: str = "nom"
    jerRegion: Optional[JecConfigReader.JerBin] = None
    jesSystTag: Optional[str] = None
    jesSystVar: Optional[str] = None


@dataclass
class JecContext:
    """
    Lightweight context for JEC/MET operations, mirroring the spirit of:

      struct JecContext {
        std::string year;
        bool isData;
        std::optional<double> runNumber;
        bool isDebug;
        JesHandles jes;
        JerHandles jer;
      };
    """
    year: str
    isData: bool
    refsAK4: JecConfigReader.CorrectionRefs
    print_evt: Optional[int] = None


# ============================================================
# JES nominal (vectorized)
# ============================================================

def jes_nominal(pt_raw, eta, phi, area, rho, year: str, isData: bool, runs,
                refs: JecConfigReader.CorrectionRefs):
    """
    Starting from raw (pT undo), apply L1, L2, (residuals on data).
    Returns:
      (pt_corrected, pt_after_l1)
    """
    pt = pt_raw

    # L1FastJet(area, eta, pt, rho)
    c1 = ColumnarHelper.eval_numeric(refs.cL1, area, eta, pt, rho)
    pt = pt * c1
    pt_after_l1 = pt

    # L2Relative
    if has_phi_dependent_L2(year):
        c2 = ColumnarHelper.eval_numeric(refs.cL2, eta, phi, pt)
    else:
        c2 = ColumnarHelper.eval_numeric(refs.cL2, eta, pt)
    pt = pt * c2

    # Data residuals
    if isData:
        if requires_run_based_residual(year):
            run_for_res = representative_run_number(year, runs)
            cR = ColumnarHelper.eval_numeric(refs.cL2L3Res, run_for_res, eta, pt)
        else:
            cR = ColumnarHelper.eval_numeric(refs.cL2L3Res, eta, pt)
        pt = pt * cR

    return pt, pt_after_l1


# ============================================================
# JES uncertainty shift (vectorized)
# ============================================================

def jes_syst_shift(pt_nom, eta, refs_cs: correction.CorrectionSet,
                   syst_name: str, var: str, print_evt: int):
    """
    Apply JES uncertainty shift on a nominal pT:

      syst node returns a fractional scale s(η, pT),
      we apply: pT → pT * (1 ± s).
    """
    syst = refs_cs[syst_name]
    scale = ColumnarHelper.eval_numeric(syst, eta, pt_nom)
    PrintHelper.debug_values("[Syst] Input jet pT:", pt_nom,
                             indent=8, event_index=print_evt)
    if var == "Up":
        shift = (1.0 + scale)
    else:
        shift = (1.0 - scale)

    out = shift * pt_nom
    PrintHelper.debug_values("[Syst] Shifted jet pT:", out,
                             indent=8, event_index=print_evt)
    return out


# ============================================================
# EventID helper for JER smearing
# ============================================================

def _event_id32_from_events(events, like, *, salt: Optional[int] = None):
    """
    Build a positive 32-bit EventID from 64-bit NanoAOD `event` only.
    Optionally XOR a constant `salt` (e.g. per-sample) to avoid collisions
    when stitching datasets. Broadcast to `like` (jet-shaped) arrays.
    """
    ev = ak.to_numpy(ak.values_astype(events, np.uint64))  # 1D: per event
    if salt is not None:
        ev ^= np.uint64(salt)
    ev = (ev & np.uint64(0x7FFFFFFF)).astype(np.int64)     # keep as positive int32 (stored in int64)
    ev_ids = ak.Array(ev)                                  # per-event
    ev_ids_b, _ = ak.broadcast_arrays(ev_ids, like)        # per-jet shape
    return ak.values_astype(ev_ids_b, np.int64)


def _correction_inputs(corr):
    """Return correctionlib input metadata as ``(name, type)`` pairs."""
    return [(inp.name, str(inp.type)) for inp in corr.inputs]


def _is_pt_input(name: str) -> bool:
    return name.lower() in {"jetpt", "pt"}


def _numeric_jer_args(corr, eta, pt):
    args = []
    for name, input_type in _correction_inputs(corr):
        if input_type == "string":
            continue
        if _is_pt_input(name):
            args.append(pt)
        else:
            args.append(eta)
    return args


def _jer_sf_from_scale_factor(refs, year, eta, pt, cat: str):
    inputs = _correction_inputs(refs.cJerSF)
    numeric_args = _numeric_jer_args(refs.cJerSF, eta, pt)
    string_positions = [i for i, (_, input_type) in enumerate(inputs) if input_type == "string"]

    if not string_positions:
        return ColumnarHelper.eval_numeric(refs.cJerSF, *numeric_args)

    if string_positions == [0]:
        return ColumnarHelper.eval_with_string_first(refs.cJerSF, *numeric_args, cat=cat)
    if string_positions == [len(inputs) - 1]:
        return ColumnarHelper.eval_with_string_last(refs.cJerSF, *numeric_args, cat=cat)

    raise RuntimeError(
        "Unsupported JER scale-factor signature for "
        f"{refs.cJerSF.name}: {inputs}"
    )


def _jer_sf_nominal(refs, year, eta, pt):
    return _jer_sf_from_scale_factor(refs, year, eta, pt, cat="nom")


def _jer_sf_uncertainty(refs, eta, pt):
    return ColumnarHelper.eval_numeric(refs.cJerSFUnc, *_numeric_jer_args(refs.cJerSFUnc, eta, pt))


def _jer_sf(refs, year, eta, pt, cat: str):
    if getattr(refs, "cJerSFUnc", None) is None:
        return _jer_sf_from_scale_factor(refs, year, eta, pt, cat=cat)

    sf = _jer_sf_nominal(refs, year, eta, pt)
    if cat == "up":
        return sf*(1 + _jer_sf_uncertainty(refs, eta, pt))
    if cat == "down":
        return sf*(1 - _jer_sf_uncertainty(refs, eta, pt))
    return sf


# ============================================================
# JER smearing (vectorized)
# ============================================================

def jer_smear(pt_in, eta, phi, rho, year: str, refs: JecConfigReader.CorrectionRefs,
              var: str, events, genIdx, gen_pt, gen_eta, gen_phi,
              mindr: float, region: Optional[JecConfigReader.JerBin],
              print_evt: int):
    """
    Apply JER smearing using the dedicated jer_smear.json.gz node.
    """
    PrintHelper.debug_values("[JER] Input jet pT:", pt_in,
                             indent=8, event_index=print_evt)

    # 1) Resolution
    reso = ColumnarHelper.eval_numeric(refs.cReso, eta, pt_in, rho)

    # 2) SF (optional region gating)
    if region is not None:
        aeta   = abs(eta)
        in_reg = (aeta >= region.etaMin) & (aeta < region.etaMax) & \
                 (pt_in >= region.ptMin) & (pt_in < region.ptMax)
        sf_nom = _jer_sf(refs, year, eta, pt_in, cat="nom")
        sf_var = _jer_sf(refs, year, eta, pt_in, cat=var)
        sf     = ak.where(in_reg, sf_var, sf_nom)
        PrintHelper.debug_values("[JER] sf :", sf, indent=8, event_index=print_evt)
    else:
        sf = _jer_sf(refs, year, eta, pt_in, cat=var)
        PrintHelper.debug_values("[JER] sf :", sf, indent=8, event_index=print_evt)

    # 3) Gen match
    gpt, valid_idx = PrintHelper.take_by_index(gen_pt,  genIdx, 0.0)
    geta, _        = PrintHelper.take_by_index(gen_eta, genIdx, 0.0)
    gphi, _        = PrintHelper.take_by_index(gen_phi, genIdx, 0.0)

    dphi     = phi_mpi_pi(phi - gphi)
    dRmatch  = np.sqrt((eta - geta) * (eta - geta) + dphi * dphi)
    is_match = valid_idx & (dRmatch < mindr) & (abs(pt_in - gpt) < 3.0 * reso * pt_in)

    PrintHelper.debug_values("[JER] reso :", reso, indent=8, event_index=print_evt)
    PrintHelper.debug_values("[JER] is_match :", is_match, indent=8, event_index=print_evt)
    PrintHelper.debug_values("[JER] dRmatch :", dRmatch, indent=8, event_index=print_evt)

    genpt_for_smear = ak.where(is_match, gpt, ak.zeros_like(pt_in) - 1.0)

    # 4) Deterministic EventID from event only (add salt=... if stitching)
    event_id = _event_id32_from_events(events, pt_in)

    # ---------- broadcast → flatten ----------
    (pt_b, eta_b, gpt_b, rho_b, eid_b, reso_b, sf_b) = ak.broadcast_arrays(
        pt_in, eta, genpt_for_smear, rho, event_id, reso, sf
    )
    counts = ak.num(pt_b)

    def _flat(x):  # -> np.ndarray (1D)
        return ak.to_numpy(ak.flatten(x))

    JetPt   = _flat(pt_b)
    JetEta  = _flat(eta_b)
    GenPt   = _flat(gpt_b)
    Rho     = _flat(rho_b)
    EventID = _flat(eid_b).astype(np.int64, copy=False)
    JER     = _flat(reso_b)
    JERSF   = _flat(sf_b)

    PrintHelper.debug_values("[JER] EventID (per-event, before flatten):",
                         event_id, indent=8, event_index=print_evt)

    # ---------- evaluate JERSmear ----------
    corr = refs.cJerSmear
    try:
        # Preferred (vectorized, by-name)
        scale_flat = corr.evaluate(
            JetPt=JetPt,
            JetEta=JetEta,
            GenPt=GenPt,
            Rho=Rho,
            EventID=EventID,
            JER=JER,
            JERSF=JERSF,
        )
    except Exception:
        # Fallback: fast Python loop (guaranteed compatible)
        scale_flat = np.fromiter(
            (corr.evaluate(float(a), float(b), float(c), float(d), int(e), float(f), float(g))
             for a, b, c, d, e, f, g in zip(JetPt, JetEta, GenPt, Rho, EventID, JER, JERSF)),
            dtype=np.float64,
            count=len(JetPt),
        )

    # ---------- unflatten back to jagged ----------
    scale = ak.unflatten(scale_flat, counts)
    out   = pt_in * scale

    PrintHelper.debug_values("[JER] JER scale :", scale,
                             indent=8, event_index=print_evt)
    return out


# ============================================================
# Type-1 MET correction tools
# ============================================================

def _apply_l1_to_met_jets(jets: CollectJetMet.JetForMetColumns, rho, ctx: JecContext):
    """Apply L1FastJet to combined jets, returning pt after L1."""
    pt_raw_mu = jets.rawPt_mu_sub
    l1_corr   = ColumnarHelper.eval_numeric(ctx.refsAK4.cL1,
                                            jets.area, jets.eta, pt_raw_mu, rho)
    pt_l1 = pt_raw_mu * l1_corr
    PrintHelper.debug_values("[MET] Combined jet pt after L1FastJet:",
                             pt_l1, indent=8, event_index=ctx.print_evt)
    return pt_l1


def _apply_l2_to_met_jets(jets: CollectJetMet.JetForMetColumns, pt_l1, ctx: JecContext):
    """Apply L2Relative to combined jets."""
    if has_phi_dependent_L2(ctx.year):
        l2_corr = ColumnarHelper.eval_numeric(ctx.refsAK4.cL2, jets.eta, jets.phi, pt_l1)
    else:
        l2_corr = ColumnarHelper.eval_numeric(ctx.refsAK4.cL2, jets.eta, pt_l1)
    pt_tmp = pt_l1 * l2_corr
    PrintHelper.debug_values("[MET] Combined jet pt after L2Relative:",
                             pt_tmp, indent=8, event_index=ctx.print_evt)
    return pt_tmp


def _apply_residuals_to_met_jets(jets: CollectJetMet.JetForMetColumns, pt_in,
                                 arrs, ctx: JecContext):
    """Apply residuals (data only) to combined jets."""
    pt_tmp = pt_in

    if ctx.isData and hasattr(ctx.refsAK4, "cL2L3Res") and ctx.refsAK4.cL2L3Res is not None:
        if requires_run_based_residual(ctx.year):
            run_for_res = representative_run_number(ctx.year, arrs["run"])
            res_corr = ColumnarHelper.eval_numeric(ctx.refsAK4.cL2L3Res,
                                                   run_for_res, jets.eta, pt_tmp)
        else:
            res_corr = ColumnarHelper.eval_numeric(ctx.refsAK4.cL2L3Res,
                                                   jets.eta, pt_tmp)
        pt_tmp = pt_tmp * res_corr
        PrintHelper.debug_values("[MET] Combined jet pt after residuals:",
                                 pt_tmp, indent=8, event_index=ctx.print_evt)

    return pt_tmp


def _apply_jes_syst_to_met_jets(jets: CollectJetMet.JetForMetColumns, pt_in,
                                ctx: JecContext, syst: SystematicOptions):
    """Apply JES systematics (MC only) to combined jets."""
    pt_tmp = pt_in

    if (not ctx.isData) and syst.syst_kind == "JES" and syst.jesSystTag:
        cnode = ctx.refsAK4.cs[syst.jesSystTag]
        scale = ColumnarHelper.eval_numeric(cnode, jets.eta, pt_tmp)
        shift = (1.0 + scale) if syst.jesSystVar == "Up" else (1.0 - scale)
        pt_tmp = pt_tmp * shift
        PrintHelper.debug_values("[MET] Combined jet pt after JES systematic:",
                                 pt_tmp, indent=8, event_index=ctx.print_evt)

    return pt_tmp


def _apply_jer_to_met_jets(jets: CollectJetMet.JetForMetColumns, pt_in,
                           arrs, ctx: JecContext, syst: SystematicOptions):
    """
    Apply JER smearing to combined jets for MET propagation.
    Data: returns pt_in unchanged, with a debug print.
    """
    if ctx.isData:
        PrintHelper.debug_values("[MET] Combined jet pt for MET propagation:",
                                 pt_in, indent=8, event_index=ctx.print_evt)
        return pt_in

    gen_pt  = arrs["GenJet_pt"]
    gen_eta = arrs["GenJet_eta"]
    gen_phi = arrs["GenJet_phi"]

    # ΔR matching for CorrT1 jets
    ct1_eta = arrs["CorrT1METJet_eta"]
    ct1_phi = arrs["CorrT1METJet_phi"]

    dphi_ct1_gen = phi_mpi_pi(ct1_phi[:, :, None] - gen_phi[:, None, :])
    deta_ct1_gen = ct1_eta[:, :, None] - gen_eta[:, None, :]
    dr2_ct1_gen  = dphi_ct1_gen * dphi_ct1_gen + deta_ct1_gen * deta_ct1_gen
    argmin_gen   = ak.argmin(dr2_ct1_gen, axis=2, keepdims=False)
    min_dr2      = ak.min(dr2_ct1_gen, axis=2, initial=1e9)
    ct1_gen_idx  = ak.where(min_dr2 < (0.2 * 0.2), argmin_gen, -1)
    ct1_gen_idx  = ak.values_astype(ct1_gen_idx, np.int32)

    ak4_gen_idx = ak.values_astype(
        arrs.get("Jet_genJetIdx", ak.zeros_like(arrs["Jet_pt"]) - 1),
        np.int32,
    )

    comb_gen_idx = ak.concatenate([ct1_gen_idx, ak4_gen_idx], axis=1)

    jer_var_eff = syst.jerVar if syst.syst_kind == "JER" else "nom"

    rho = arrs["Rho_fixedGridRhoFastjetAll"]
    pt_corr = jer_smear(
        pt_in=pt_in, eta=jets.eta, phi=jets.phi, rho=rho,
        year=ctx.year, refs=ctx.refsAK4,
        var=jer_var_eff, events=arrs["event"],
        genIdx=comb_gen_idx, gen_pt=gen_pt,
        gen_eta=gen_eta, gen_phi=gen_phi,
        mindr=0.2, region=syst.jerRegion,
        print_evt=ctx.print_evt,
    )

    PrintHelper.debug_values("[MET] Combined jet pt after JER smearing:",
                             pt_corr, indent=8, event_index=ctx.print_evt)
    return pt_corr


def _recompute_type1_met(met_raw: CollectJetMet.MetInputs,
                         jets: CollectJetMet.JetForMetColumns,
                         pt_l1, pt_corr,
                         ctx: JecContext):
    """
    Standard Type-1 MET recomputation:

      - Select jets with pT>15, |η|<5.2, EM frac < 0.9
      - Compute ΔpT = pT_corr - pT_L1
      - Shift MET by -Σ ΔpT · (cosφ, sinφ)
    """
    pass_for_met = (pt_corr > 15.0) & (abs(jets.eta) < 5.2) & \
                   ((jets.chEmEf + jets.neEmEf) < 0.9)

    dpt    = ak.where(pass_for_met, pt_corr - pt_l1, 0.0)
    sum_px = ak.sum(dpt * np.cos(jets.phi), axis=1)
    sum_py = ak.sum(dpt * np.sin(jets.phi), axis=1)

    met_px = met_raw.px - sum_px
    met_py = met_raw.py - sum_py
    met_pt_corr = np.hypot(met_px, met_py)

    PrintHelper.debug_values("[MET] Type-1 corrected MET pt (CorrT1+AK4):",
                             met_pt_corr, indent=4, event_index=ctx.print_evt)
    return met_pt_corr


# ============================================================
# Public MET interface (same signature as before)
# ============================================================

def get_corrected_met(
    *,
    arrs,
    year: str,
    isData: bool,
    refsAK4: JecConfigReader.CorrectionRefs,
    syst_kind: str,                 # "Nominal", "JES", "JER"
    jes_ak4_tag: Optional[str] = None,
    jes_var: Optional[str] = None,  # "Up"/"Down" when syst_kind == "JES"
    jer_bin: Optional[JecConfigReader.JerBin] = None,
    jer_var: Optional[str] = None,  # "up"/"down"/"nom"
    print_evt: Optional[int] = None,
):
    """
    Compute Type-1 corrected MET (CorrT1 + AK4) with JES/JER propagation.

    Returns:
        met_pt_corr : awkward.Array[float]  # per event
    """
    # Context + systematic options (POD-like)
    ctx = JecContext(year=year, isData=isData, refsAK4=refsAK4,
                     print_evt=print_evt)
    syst = SystematicOptions(
        syst_kind=syst_kind,
        jerVar=jer_var or "nom",
        jerRegion=jer_bin,
        jesSystTag=jes_ak4_tag,
        jesSystVar=jes_var,
    )

    # 1) Raw MET
    met_raw = CollectJetMet.build_met_inputs(arrs, year, print_evt)

    # 2) Build CorrT1 + AK4 jet collections
    jets_ct1 = CollectJetMet.build_corrt1_jets_for_met(arrs, print_evt)
    jets_ak4 = CollectJetMet.build_ak4_jets_for_met(arrs, print_evt)
    jets_all = CollectJetMet.combine_met_jets(jets_ct1, jets_ak4)

    # 3) Apply L1, L2, residuals, JES syst
    rho       = arrs["Rho_fixedGridRhoFastjetAll"]
    pt_l1     = _apply_l1_to_met_jets(jets_all, rho, ctx)
    pt_afterL2= _apply_l2_to_met_jets(jets_all, pt_l1, ctx)
    pt_res    = _apply_residuals_to_met_jets(jets_all, pt_afterL2, arrs, ctx)
    pt_jes    = _apply_jes_syst_to_met_jets(jets_all, pt_res, ctx, syst)

    # 4) Apply JER (MC only)
    pt_corr_met = _apply_jer_to_met_jets(jets_all, pt_jes, arrs, ctx, syst)

    # 5) Type-1 recompute
    met_pt_corr = _recompute_type1_met(met_raw, jets_all, pt_l1, pt_corr_met, ctx)
    return met_pt_corr

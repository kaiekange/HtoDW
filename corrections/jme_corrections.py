"""JME (JEC + JER) for the custom NanoAOD, via the official jerc-application-tutorial.

WHY THIS EXISTS
---------------
Our custom NanoAOD stores `Jet_pt` already corrected by whatever JEC payload the
GlobalTag `130X_mcRun3_2022_realistic_v5` carried, and **not smeared at all**:
`grep -n smear jetsAK4_Puppi_cff.py` returns nothing, so the standard nano jet chain
(jetPuppiCorrFactorsNano -> updatedJetsPuppi -> updatedJetsPuppiWithUserData ->
jetPuppiTable) does L1FastJet+L2Relative+L3Absolute+L2L3Residual and stops. For MC
L3Absolute = L2L3Residual = 1.0, so the level structure is identical to what the
current payload applies -- the residual nano-vs-rederived difference is a *payload
version* difference, which was verified to be deterministic (spread inside fine
(raw pt, |eta|) bins is 0.1-1%, bin means drift smoothly 1.014 -> 0.983).

So the analysis quantity is NOT `Jet_pt`. It is:

    pt_raw  = Jet_pt * (1 - Jet_rawFactor)          # undo the GT payload
    pt_jes  = jes_nominal(pt_raw, ...)              # Summer22_22Sep2023_V4_MC
    pt_corr = jer_smear(pt_jes, ..., var="nom")     # Summer22_22Sep2023_JRV2_MC

which is exactly the official `ForPythonUser/applyJecAndJvmLite.py` chain.

WHICH pt EACH CONSUMER USES  (read this before touching anything downstream)
---------------------------------------------------------------------------
There are two jet pt's in play and mixing them silently is the main hazard here:

  * nano `Jet_pt`  -- what the *producer* saw. Any closure check against the stored
    producer quantities must keep using this, because it validates the
    producer's own arithmetic against the stored `n_b_jets_tight_dR04`. Feeding it
    corrected pt would make it fail by construction.
  * `pt_corr` from here -- the analysis quantity. Everything physics-facing uses it:
    jet_ht / num_jets* / jet_leadpt, the b-tag SF pt argument, the efficiency-map pt
    binning and its pt>20 selection, and MET (via get_corrected_met, which redoes its
    own internal jet chain -- do not hand-roll a second one).

TWO SEPARATE QUESTIONS: does JME move a quantity, and can nano rebuild it?
-------------------------------------------------------------------------
JME rescales jet **pt only** -- eta, phi and the b-tag discriminator are untouched.
So any stored quantity whose definition contains no pt is JME-invariant and its stored
value stays correct no matter what we do here.

Separately, nano's jetTable stores `finalJetsPuppi` = pt > 15 GeV, while every jet
loop in HDsPhiPiMuProducer.cc cuts only on userInt("tightId")!=0 (== nano `Jet_jetId`
bit 2) and, in the candidate loop, |eta|<2.5 -- **no pt floor anywhere**. So the
producer counted jets nano does not contain. Measured on Signal_DsW output_0: mean
1.16 jets and 14.0 GeV of HT missing per event, i.e. 9.6% of the stored `jet_ht`
(median 7.5%, p90 23.8%), from jets averaging 12.1 GeV.

Crossing the two:

  JME-invariant, so nothing to do (their definitions have no pt):
      num_jets, jet_maxbscore, Ds_nearestjet_dR, Ds_nearestjet_bscore,
      n_b_jets_{loose,medium,tight}_dR0{3,4}
      -- several of these do NOT reproduce from nano (Ds_nearestjet_dR mismatches
      34-48% of events, 65-82% of that from the missing sub-15 GeV jets), but that
      is irrelevant here: we never need to recompute them.

  JME-dependent and rebuildable from nano (closure 0.03-0.6%, residuals all in
  events whose jets sit below nano's threshold):
      num_jets_pt20/30/40/50, jet_ht_pt20, jet_leadpt, n_b_jets_*_pt20_*,
      and MET via get_corrected_met

  JME-dependent and NOT rebuildable -- exactly one:
      jet_ht

`jet_ht` is a DNN feature and it sums pt with no threshold, so 9.6% of it lives in
jets nano dropped. Its pt>20 counterpart `jet_ht_pt20` is already stored, closes, and
discards only a soft component that is largely pileup anyway. Swapping the feature set
onto `jet_ht_pt20` is a retraining decision and is deliberately NOT made here.

The b-veto jet set is safe: `n_b_jets_tight_dR04` rebuilds from nano with no tolerance
at all on 34,234/34,235 events across signal and all three W+jets jet bins, and the
single mismatch is not in an event with a missing sub-15 GeV jet. A jet under 15 GeV
passing the tight PNet-B WP (0.6734) essentially does not occur.

ENVIRONMENT
-----------
correctionlib >= 2.6 (el9). el8 CMSSW_13_0_13 ships 2.1.0/2.2.2, which cannot parse
these schema_version=2 JSONs.

The tutorial modules are an unmodified copy under ../external/jerc-application-tutorial
(upstream commit in its UPSTREAM file), so a resync is a plain diff against upstream.
They resolve their config files by paths relative to the
CWD -- `JecConfigReader` line 93 hardcodes "../../jer_smear.json.gz", so the chdir is
not avoidable by setting CONFIG_PATHS. `_in_tutorial_dir()` below contains that.
"""

import contextlib
import os
import sys

import awkward as ak
import numpy as np

TUTORIAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "external",
                            "jerc-application-tutorial", "ApplyOnNanoAOD", "ForPythonUser")
TUTORIAL_DIR = os.path.normpath(TUTORIAL_DIR)

if TUTORIAL_DIR not in sys.path:
    sys.path.insert(0, TUTORIAL_DIR)

import JecApplication      # noqa: E402
import JecConfigReader     # noqa: E402
import JvmApplication      # noqa: E402
import JvmConfigReader     # noqa: E402

# 2022 C+D == Summer22 preEE. Key into JecConfigAK4.json.
YEAR_PREEE = "2022Pre"

# AK4 R=0.4 -> gen-match cone R/2, as in applyJecAndJvmLite.py.
JER_MATCH_DR = 0.2

# The audit settled on the Total (single-nuisance) JES scheme; the JER counterpart is
# the one-bin JerUncertaintySetTotal covering |eta|<6, all pt. The JES Total name is the
# same in every era; the JER Total label is era-specific (2022Pre: CMS_res_j_2022,
# 2022Post: CMS_res_j_2022EE), so it is taken as the set's single entry, not by name.
JES_TOTAL_NAME = "CMS_scale_j_Total"


@contextlib.contextmanager
def _in_tutorial_dir():
    """Run inside ForPythonUser/ -- its config paths are relative to the CWD."""
    prev = os.getcwd()
    os.chdir(TUTORIAL_DIR)
    try:
        yield
    finally:
        os.chdir(prev)


class JmeRefs:
    """Everything the correction calls need, loaded once."""

    def __init__(self, year=YEAR_PREEE, is_data=False, era=None):
        self.year = year
        self.is_data = is_data
        with _in_tutorial_dir():
            self.cfg = JecConfigReader.CorrCache.get_cfg("AK4")
            self.refs = JecConfigReader.CorrCache.get_refs("AK4", year, is_data, era)

            if is_data:
                self.jes_total_tag = None
                self.jer_total_bin = None
                return

            # get_syst_tag_names yields (fullTag, baseName) pairs.
            syst_tags = JecConfigReader.get_syst_tag_names(self.cfg, year)
            self.jes_total_tag = next(
                tag for tag, base in syst_tags["JesUncertaintySetTotal"]
                if base == JES_TOTAL_NAME
            )

            # get_jer_sets indexes the *full* config by year -- pass self.cfg, not cfg[year].
            jer_sets = JecConfigReader.get_jer_sets(self.cfg, year)
            total_bins = jer_sets["JerUncertaintySetTotal"]
            if len(total_bins) != 1:
                raise RuntimeError(f"{year}: expected exactly one JerUncertaintySetTotal bin,"
                                   f" got {[b.label for b in total_bins]}")
            self.jer_total_bin = total_bins[0]


def _genjets_safe(arrs):
    """GenJet_{pt,eta,phi} and Jet_genJetIdx made safe for the tutorial's gen matching.

    PrintHelper.take_by_index replaces an invalid index by 0 and then does arr[idx_safe],
    which RAISES on an event with zero GenJets (found in W+jets 0J output_110.root: an event
    with 5 jets, 0 GenJets, all genJetIdx=-1). Only empty GenJet lists are padded, with one
    dummy entry at eta=999 so its dR to any reco jet is enormous: AK4 jets keep genIdx=-1 and
    the MET path's CorrT1 dR matching (argmin over GenJets, dR<0.2) can never pick it. Events
    with at least one GenJet are untouched, so results there are bit-identical. As a belt-and-
    braces guard, a genJetIdx beyond the stored GenJet count is forced to -1 (nano writes -1
    for a match that was not stored, so this should never trigger).
    """
    n_gen = ak.num(arrs["GenJet_pt"])

    def pad(x, dummy):
        return ak.fill_none(ak.pad_none(x, 1, axis=1), dummy)

    gen_pt, gen_eta, gen_phi = (pad(arrs["GenJet_pt"], 0.0), pad(arrs["GenJet_eta"], 999.0),
                                pad(arrs["GenJet_phi"], 0.0))
    idx = arrs["Jet_genJetIdx"]
    n_gen_b, _ = ak.broadcast_arrays(n_gen, idx)
    idx = ak.where(idx < n_gen_b, idx, -1)
    return gen_pt, gen_eta, gen_phi, idx


def _jet_inputs(arrs):
    """Unpack the jet columns the JEC chain needs, broadcasting rho to jet shape."""
    pt_raw = arrs["Jet_pt"] * (1.0 - arrs["Jet_rawFactor"])
    rho_b, _ = ak.broadcast_arrays(arrs["Rho_fixedGridRhoFastjetAll"], arrs["Jet_pt"])
    return pt_raw, rho_b


def corrected_jet_pt(arrs, refs, variation="nominal"):
    """Jet pt after the full official chain, for one systematic variation.

    variation: "nominal" | "jes_up" | "jes_down" | "jer_up" | "jer_down"

    Order: JES nominal -> JES shift (if any) -> JER smear. The JES shift goes in
    *before* smearing, and the JER variation is selected by `var` inside jer_smear
    rather than applied afterwards.

    That ordering is the official one, from the primary source: applyJecAndJvm.C runs
    one full pass per variation and its step 4 is commented "JER smearing (MC only,
    after JES)". (applyJecAndJvmLite.py shows the same sequence but on a single dummy
    jet, so on its own it could not settle this.)

    Consequence worth knowing: jer_smear.json's stochastic branch draws from
    `hashprng` with entropy inputs [JetPt, JetEta, Rho, EventID] -- JetPt is in the
    hash, so shifting JES *redraws* the random number for gen-unmatched jets (20% of
    jets here). Measured on Signal_DsW output_0: the event-level <jet_ht_pt20> JES-Up
    ratio is 1.0479 with the official ordering vs 1.0412 with the draw held fixed, so
    the redraw inflates that systematic by 0.68pp out of ~4% (medians 1.0271 vs
    1.0260, i.e. it lives in the tails). We follow the official ordering anyway; this
    note exists so the 0.7pp is a known quantity rather than an unexplained one.
    """
    if variation not in ("nominal", "jes_up", "jes_down", "jer_up", "jer_down"):
        raise ValueError(f"unknown variation {variation!r}")

    pt_raw, rho_b = _jet_inputs(arrs)
    eta, phi = arrs["Jet_eta"], arrs["Jet_phi"]
    # data: L2L3Residual is looked up at a representative run per year, but the tutorial
    # still wants the run column broadcast to jet shape
    runs = ak.broadcast_arrays(arrs["run"], arrs["Jet_pt"])[0] if refs.is_data else None

    with _in_tutorial_dir():
        pt_jes, _ = JecApplication.jes_nominal(
            pt_raw=pt_raw, eta=eta, phi=phi, area=arrs["Jet_area"], rho=rho_b,
            year=refs.year, isData=refs.is_data, runs=runs, refs=refs.refs,
        )

        # JES systematics and JER smearing are MC-only (JecApplication applies the JES shift
        # only when not isData, and jer_smear needs GenJets, which data does not have)
        if refs.is_data:
            if variation != "nominal":
                raise ValueError(f"data supports only the nominal jet pt, got {variation!r}")
            return pt_jes

        if variation.startswith("jes_"):
            pt_jes = JecApplication.jes_syst_shift(
                pt_jes, eta, refs.refs.cs, refs.jes_total_tag,
                "Up" if variation == "jes_up" else "Down", None,
            )

        if variation == "jer_up":
            jer_var, region = "up", refs.jer_total_bin
        elif variation == "jer_down":
            jer_var, region = "down", refs.jer_total_bin
        else:
            jer_var, region = "nom", None

        gen_pt, gen_eta, gen_phi, gen_idx = _genjets_safe(arrs)
        pt_corr = JecApplication.jer_smear(
            pt_in=pt_jes, eta=eta, phi=phi, rho=rho_b,
            year=refs.year, refs=refs.refs, var=jer_var,
            events=arrs["event"],
            genIdx=gen_idx,
            gen_pt=gen_pt, gen_eta=gen_eta,
            gen_phi=gen_phi,
            mindr=JER_MATCH_DR, region=region, print_evt=None,
        )
    return pt_corr


def jet_veto_flag(arrs, jet_pt, year=YEAR_PREEE):
    """Per-event True if any jet lands in the jet veto map, via the official routine.

    JvmApplication.apply_jvm_flag selects jets before the map lookup -- |eta|<=5.191,
    |phi|<=pi, jetId>=6 (TightLepVeto), pt>=15, chEmEF+neEmEF<=0.9 -- identical to the
    C++ VetoChecker::passMinimalSelection (JvmApplication.h/.cpp).

    `jet_pt` must be the CORRECTED pt of the variation being evaluated. Primary source:
    applyJecAndJvm.C applies JES and JER in place to nanoT.Jet_pt (ApplicationHelper.hpp
    `nanoT.Jet_pt[idx] *= s` / `*= sJer`) and only then runs the veto map as step 6, once
    per pass -- so the pt>=15 cut sees the corrected pt, and a JES/JER variation can move
    a jet across it. (applyJecAndJvmLite.py passes nano pt, but that is a one-jet demo.)
    The tag/key come from the tutorial's JvmConfig.json: for 2022Pre that is
    Summer22_23Sep2023_RunCD_V1 / "jetvetomap".
    """
    with _in_tutorial_dir():
        jvm_ref, jvm_key = JvmConfigReader.CorrCache.get_jvm(year)
        return JvmApplication.apply_jvm_flag(
            ak4_eta=arrs["Jet_eta"], ak4_phi=arrs["Jet_phi"], ak4_pt=jet_pt,
            ak4_id=arrs["Jet_jetId"], ak4_chEm=arrs["Jet_chEmEF"],
            ak4_neEm=arrs["Jet_neEmEF"], jvm_ref=jvm_ref, jvm_key=jvm_key,
        )


def corrected_met_pt(arrs, refs, variation="nominal"):
    """Type-1 corrected MET pt, JES/JER propagated, via the official routine.

    get_corrected_met rebuilds its own CorrT1+AK4 jet chain internally -- pass it the
    raw nano columns, not the pt from corrected_jet_pt().
    """
    # CollectJetMet calls arrs.get(...) for optional columns, so it needs a real dict;
    # an awkward record array raises AttributeError on .get.
    if not isinstance(arrs, dict):
        arrs = {f: arrs[f] for f in arrs.fields}
    else:
        arrs = dict(arrs)
    # same zero-GenJet crash guard as corrected_jet_pt (see _genjets_safe); data has no GenJets
    if not refs.is_data:
        (arrs["GenJet_pt"], arrs["GenJet_eta"], arrs["GenJet_phi"],
         arrs["Jet_genJetIdx"]) = _genjets_safe(arrs)
    kind, jes_tag, jes_var, jer_bin, jer_var = "Nominal", None, None, None, "nom"
    if variation == "jes_up":
        kind, jes_tag, jes_var = "JES", refs.jes_total_tag, "Up"
    elif variation == "jes_down":
        kind, jes_tag, jes_var = "JES", refs.jes_total_tag, "Down"
    elif variation == "jer_up":
        kind, jer_bin, jer_var = "JER", refs.jer_total_bin, "up"
    elif variation == "jer_down":
        kind, jer_bin, jer_var = "JER", refs.jer_total_bin, "down"
    elif variation != "nominal":
        raise ValueError(f"unknown variation {variation!r}")

    with _in_tutorial_dir():
        return JecApplication.get_corrected_met(
            arrs=arrs, year=refs.year, isData=refs.is_data, refsAK4=refs.refs,
            syst_kind=kind, jes_ak4_tag=jes_tag, jes_var=jes_var,
            jer_bin=jer_bin, jer_var=jer_var, print_evt=None,
        )


def corrected_met_pt_phi(arrs, refs, variation="nominal"):
    """(pt, phi) of the Type-1 MET, same arithmetic as corrected_met_pt().

    get_corrected_met returns pt only, but W_tmass needs the phi too. This runs the very
    same official steps as JecApplication.get_corrected_met (steps 1-4 called verbatim)
    and replaces only the last line of _recompute_type1_met, np.hypot(met_px, met_py), by
    keeping met_px/met_py. Callers should assert pt == corrected_met_pt() on real input.
    """
    if not isinstance(arrs, dict):
        arrs = {f: arrs[f] for f in arrs.fields}
    else:
        arrs = dict(arrs)
    if not refs.is_data:
        (arrs["GenJet_pt"], arrs["GenJet_eta"], arrs["GenJet_phi"],
         arrs["Jet_genJetIdx"]) = _genjets_safe(arrs)
    kind, jes_tag, jes_var, jer_bin, jer_var = "Nominal", None, None, None, "nom"
    if refs.is_data and variation != "nominal":
        raise ValueError(f"data supports only the nominal MET, got {variation!r}")
    if variation == "jes_up":
        kind, jes_tag, jes_var = "JES", refs.jes_total_tag, "Up"
    elif variation == "jes_down":
        kind, jes_tag, jes_var = "JES", refs.jes_total_tag, "Down"
    elif variation == "jer_up":
        kind, jer_bin, jer_var = "JER", refs.jer_total_bin, "up"
    elif variation == "jer_down":
        kind, jer_bin, jer_var = "JER", refs.jer_total_bin, "down"
    elif variation != "nominal":
        raise ValueError(f"unknown variation {variation!r}")

    JA = JecApplication
    with _in_tutorial_dir():
        ctx = JA.JecContext(year=refs.year, isData=refs.is_data, refsAK4=refs.refs, print_evt=None)
        syst = JA.SystematicOptions(syst_kind=kind, jerVar=jer_var or "nom", jerRegion=jer_bin,
                                    jesSystTag=jes_tag, jesSystVar=jes_var)
        met_raw = JA.CollectJetMet.build_met_inputs(arrs, refs.year, None)
        jets_all = JA.CollectJetMet.combine_met_jets(
            JA.CollectJetMet.build_corrt1_jets_for_met(arrs, None),
            JA.CollectJetMet.build_ak4_jets_for_met(arrs, None))
        rho = arrs["Rho_fixedGridRhoFastjetAll"]
        pt_l1 = JA._apply_l1_to_met_jets(jets_all, rho, ctx)
        pt_l2 = JA._apply_l2_to_met_jets(jets_all, pt_l1, ctx)
        pt_res = JA._apply_residuals_to_met_jets(jets_all, pt_l2, arrs, ctx)
        pt_jes = JA._apply_jes_syst_to_met_jets(jets_all, pt_res, ctx, syst)
        pt_corr = JA._apply_jer_to_met_jets(jets_all, pt_jes, arrs, ctx, syst)
        # _recompute_type1_met, verbatim except that px/py are kept
        jets = jets_all
        pass_for_met = (pt_corr > 15.0) & (abs(jets.eta) < 5.2) & \
                       ((jets.chEmEf + jets.neEmEf) < 0.9)
        dpt = ak.where(pass_for_met, pt_corr - pt_l1, 0.0)
        met_px = met_raw.px - ak.sum(dpt * np.cos(jets.phi), axis=1)
        met_py = met_raw.py - ak.sum(dpt * np.sin(jets.phi), axis=1)
    return np.hypot(met_px, met_py), np.arctan2(met_py, met_px)


# Branches corrected_jet_pt() needs; corrected_met_pt() needs these plus MET_BRANCHES.
JET_BRANCHES = [
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Jet_genJetIdx",
    "GenJet_pt", "GenJet_eta", "GenJet_phi",
    "Rho_fixedGridRhoFastjetAll", "event",
]

JVM_BRANCHES = ["Jet_eta", "Jet_phi", "Jet_jetId", "Jet_chEmEF", "Jet_neEmEF"]

MET_BRANCHES = [
    "RawPuppiMET_pt", "RawPuppiMET_phi", "RawMET_pt", "RawMET_phi", "MET_pt",
    "CorrT1METJet_rawPt", "CorrT1METJet_eta", "CorrT1METJet_phi", "CorrT1METJet_area",
    "CorrT1METJet_muonSubtrFactor",
    "Jet_muonSubtrFactor", "Jet_chEmEF", "Jet_neEmEF",
    "FatJet_pt", "FatJet_eta", "FatJet_phi",
]

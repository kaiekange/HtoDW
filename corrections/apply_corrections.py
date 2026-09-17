#!/usr/bin/env python3
"""POG corrections for ONE custom-nano output file, written as a friend tree.

Usage: apply_corrections.py <.../NanoMVA/MC/<era>/<sample dir>/output_N.root>
       apply_corrections.py <.../NanoMVA/Data/<era>/Muon/<task>/<timestamp>/0000/nano_N.root> [output]
Output:                     .../NanoMVA/MC/<era>/Corrections/<sample dir>/output_N.root
                            (or [output], e.g. for a side-by-side comparison run)
                            tree "Events", one row per nano Events row, same order.

The era and sample are read from the input path, so a job cannot be handed a JSON set,
WP, JME year or efficiency map that does not belong to the file.

Join: row i of the friend is row i of the nano file (the nano is one row per event).
run / luminosityBlock / event are copied in so that every reader can -- and should --
assert they match before using a column.

Physics choices are the ones settled for this analysis (memory
project_corrections_pilot_implementation). The MuonScaRe port lives in
muo_scalesmearing.py.

VARIATIONS. Suffixes: _nom; muon momentum _muScaleUp/_muScaleDown/_muResUp/_muResDown;
JME _jesUp/_jesDown/_jerUp/_jerDown (Total). Each variation re-derives only what it moves:
  muon : mu_pt, pass_mu, w_mu (id x iso x trig, nominal SFs at that pt), W_tmass
         (with nominal MET), Dsmu_{pt,eta,phi,mass,Ht}, Dsstmu_{pt,eta,phi,mass,Ht}
         (mu_px/py/pz follow from the stored direction times mu_pt; mu_IsoR0x_RelIso is
         left as stored, since the POG iso WP and its SF are defined on the uncorrected pt)
  JME  : jet_ht_pt20, num_jets_pt30, jet_leadpt, met_pt, met_phi, W_tmass (with nominal
         muon), jvm_veto, n_bjet, w_btag (central SF at that JME variation)
SF systematics at the nominal kinematics: w_muid/w_muiso/w_mutrig _nom/_up/_down,
w_btag_up/_down, w_pu_nom/_up/_down. (W_tmass_nom, from nominal muon and MET, is written
once.)

COLUMN SEMANTICS
  pass_mu      selected muon: corrected pt > 26, |eta| < 2.4, PFIsoTight (Muon_pfIsoId>=4)
  pass_muveto  no other muon with tight ID, |eta|<2.4, Loose PF iso, corrected pt >= 15
               (nominal muon correction only)
  jvm_veto     True = event VETOED by the jet veto map (official apply_jvm_flag)
  n_bjet       jets with tight ID, |eta|<2.5, dR(Ds)>0.4, PNet B > era tight WP,
               corrected pt >= 20; the b-veto is n_bjet == 0
  w_mu*        NaN where pass_mu (of that variation) is False
  w_btag*      NaN where n_bjet (of that variation) > 0
NaN, not 1, so a weight used without its cut poisons the sum instead of passing silently.

NOT CORRECTED (stored nano values stay as they are): met_sig and met_sumEt -- the JME
tutorial's Type-1 propagation defines no corrected significance or sumEt. jet_ht and
num_jets cannot be rebuilt (nano stores jets with pt > 15 only); jet_ht_pt20 replaces jet_ht.
The producer's candidate choice and its Dsmu_mass < 125 cut used the uncorrected muon.

Runs inside cmssw-el9 + CMSSW_15_0_5 (correctionlib 2.6.4).
"""
import os
import shutil
import sys
import tempfile

import awkward as ak
import correctionlib
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_btag_efficiency_map as bem  # noqa: E402  (ERAS: jme_year, wp_tight; OUT_DIR)
import muo_scalesmearing as muo  # noqa: E402  (MuonScaRe port, MUO keys)
import jme_corrections as jme  # noqa: E402

POG_DIR = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
CAMPAIGN = {
    "2022": "Run3-22CDSep23-Summer22-NanoAODv12",
    "2022EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
}
PU_KEY = {
    "2022": "Collisions2022_355100_357900_eraBCD_GoldenJson",
    "2022EE": "Collisions2022_359022_362760_eraEFG_GoldenJson",
}
# sample directory -> efficiency-map group (build_btag_efficiency_map.GROUPS minus era).
# Every sample uses its OWN map. The three W+jets jet bins are separate samples: measured
# on 2022 they have clearly different efficiencies in the same (|eta|, pt) cell (b jets
# 34% / 40% / 57% overall for 0J / 1J / 2J), so a pooled map would impose the 2J-dominated
# efficiency on 0J and 1J. The bins are combined only afterwards, each with its own norm.
EFF_GROUP = {
    "Signal_DsW_H_phipimu": "DsW",
    "Signal_DsstW_H_phipimu": "DsstW",
    "WtoLNu2Jets_0J_H_phipimu_GenTag": "WtoLNu2Jets_0J",
    "WtoLNu2Jets_1J_H_phipimu_GenTag": "WtoLNu2Jets_1J",
    "WtoLNu2Jets_2J_H_phipimu_GenTag": "WtoLNu2Jets_2J",
}

MU_PRECUT = 24.0       # HDsPhiPiMuProducer: best muon = highest-pt tight muon, pt >= 24
MU_PT_MIN = 26.0       # trigger-SF floor, on corrected pt
MU_ETA_MAX = 2.4
VETO_MU_PT_MIN = 15.0  # ID-SF floor
BTAG_JET_PT_MIN = 20.0
MASS_MU = 0.1056583755
NULL = -777.0

# DATA. Recognised by the path. Differences from MC, all of them because the quantity does
# not exist for data rather than by choice: JEC gains L2L3Residual (era key below), there is
# no JER smearing and no JES variation (JecApplication applies neither on data), the muon
# gets the scale correction from the data branch of MuonScaRe and no resolution smearing,
# and no weight is defined (no PU, no scale factors, no b-tag weight). The selection flags
# are computed exactly as for MC, because the same cuts are applied to data.
DATA_ERA_KEY = {"2022": "Era2022PreAll", "2022EE": "Era2022PostAll"}

MU_VARS = ("nom", "muScaleUp", "muScaleDown", "muResUp", "muResDown")
JME_VARS = {"nom": "nominal", "jesUp": "jes_up", "jesDown": "jes_down",
            "jerUp": "jer_up", "jerDown": "jer_down"}

MC_ONLY_BRANCHES = {"Pileup_nTrueInt", "Jet_hadronFlavour", "Jet_genJetIdx",
                    "GenJet_pt", "GenJet_eta", "GenJet_phi"}

BRANCHES = sorted(set([
    "run", "luminosityBlock", "event", "Pileup_nTrueInt",
    "mu_pt", "mu_eta", "mu_phi", "mu_charge", "mu_px", "mu_py", "mu_pz",
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_nTrackerLayers",
    "Muon_pfIsoId",
    "met_pt", "met_phi", "W_tmass",
    "Dsmu_px", "Dsmu_py", "Dsmu_pz", "Dsmu_mass", "Dsmu_pt", "Dsmu_eta", "Dsmu_phi", "Dsmu_Ht",
    "Dsstmu_px", "Dsstmu_py", "Dsstmu_pz", "Dsstmu_mass", "Dsstmu_pt", "Dsstmu_eta", "Dsstmu_phi",
    "Dsstmu_Ht",
    "Ds_eta", "Ds_phi", "Jet_jetId", "Jet_btagPNetB", "Jet_hadronFlavour",
] + jme.JET_BRANCHES + jme.JVM_BRANCHES + jme.MET_BRANCHES))


def parse_path(path):
    """(era, sample, output path, is_data) from the input path.

    MC  : .../NanoMVA/MC/<era>/<sample dir>/output_N.root
    data: .../NanoMVA/Data/<era>/Muon/<task>/<timestamp>/0000/nano_N.root  (CRAB layout)
          -> .../NanoMVA/Data/<era>/Corrections/<task>/nano_N.root
    """
    parts = os.path.abspath(path).split(os.sep)
    if "Data" in parts:
        i = len(parts) - 1 - parts[::-1].index("Data")
        era = parts[i + 1]
        if era not in CAMPAIGN or len(parts) < i + 6:
            sys.exit(f"cannot derive era/task from {path} (era={era!r})")
        task, fname = parts[i + 3], parts[-1]
        return era, task, os.sep.join(parts[:i + 2] + ["Corrections", task, fname]), True
    era, sample, fname = parts[-3], parts[-2], parts[-1]
    if parts[-4] != "MC" or era not in CAMPAIGN or sample not in EFF_GROUP:
        sys.exit(f"cannot derive era/sample from {path} (era={era!r}, sample={sample!r})")
    return era, sample, os.sep.join(parts[:-2] + ["Corrections", sample, fname]), False


def np_(a, b, dtype=float):
    return ak.to_numpy(a[b]).astype(dtype)


def locate_selected_muon(a):
    """Index in Muon_* of the producer's best muon; asserts it reproduces mu_pt/eta/charge."""
    sel = (a["Muon_pt"] >= MU_PRECUT) & a["Muon_tightId"]
    idx = ak.argmax(ak.where(sel, a["Muon_pt"], -1.0), axis=1, keepdims=True)
    pt, eta, ch = (ak.to_numpy(ak.flatten(a[b][idx])) for b in ("Muon_pt", "Muon_eta", "Muon_charge"))
    ok = ((np.abs(pt - np_(a, "mu_pt")) < 1e-3) & (np.abs(eta - np_(a, "mu_eta")) < 1e-3)
          & (ch == np_(a, "mu_charge", int)))
    if not ok.all():
        sys.exit(f"selected-muon re-id FAILED on {int((~ok).sum())}/{len(ok)} rows")
    return ak.to_numpy(ak.flatten(idx))


def muon_momenta(a, idx, cset, is_data=False):
    """Corrected pt of the selected muon for every MU_VARS entry (official MuonScaRe port,
    same calling convention as the corrections pilot)."""
    n = len(a)
    out = {v: np.empty(n) for v in MU_VARS}
    pt_raw = np_(a, "mu_pt")
    eta = np.clip(np_(a, "mu_eta"), -2.399, 2.399)
    phi = np_(a, "mu_phi")
    ch = np_(a, "mu_charge", int)
    nL = ak.to_numpy(ak.flatten(a["Muon_nTrackerLayers"][idx[:, None]])).astype(float)
    evt, lumi = ak.to_numpy(a["event"]), ak.to_numpy(a["luminosityBlock"])
    for i in range(n):
        p_scale = muo.muo_pt_scale(cset, is_data, pt_raw[i], eta[i], phi[i], ch[i])
        if is_data:
            out["nom"][i] = p_scale          # no resolution smearing on data
            continue
        p = muo.muo_pt_resol(cset, p_scale, eta[i], phi[i], nL[i], evt[i], lumi[i])
        out["nom"][i] = p
        try:
            out["muScaleUp"][i] = muo.muo_pt_scale_var(cset, p, eta[i], phi[i], ch[i], "up")
            out["muScaleDown"][i] = muo.muo_pt_scale_var(cset, p, eta[i], phi[i], ch[i], "dn")
        except (ZeroDivisionError, ValueError, OverflowError):
            out["muScaleUp"][i] = out["muScaleDown"][i] = p
        try:
            out["muResUp"][i] = muo.muo_pt_resol_var(cset, p_scale, p, eta[i], "up")
            out["muResDown"][i] = muo.muo_pt_resol_var(cset, p_scale, p, eta[i], "dn")
        except (ZeroDivisionError, ValueError, OverflowError):
            out["muResUp"][i] = out["muResDown"][i] = p
    return out


def prompt_muon_veto(a, idx, cset, is_data=False):
    """True = passes (no other prompt-like muon). Nominal muon correction."""
    n = len(a)
    passed = np.ones(n, dtype=bool)
    evt, lumi = ak.to_numpy(a["event"]), ak.to_numpy(a["luminosityBlock"])
    cand = ((a["Muon_tightId"] == 1) & (abs(a["Muon_eta"]) < MU_ETA_MAX)
            & (a["Muon_pfIsoId"] >= 2) & (a["Muon_pt"] >= VETO_MU_PT_MIN)
            & (ak.local_index(a["Muon_pt"]) != idx))
    for i in np.where(ak.to_numpy(ak.any(cand, axis=1)))[0]:
        for j in np.where(ak.to_numpy(cand[i]))[0]:
            eta_j = float(np.clip(a["Muon_eta"][i][j], -2.399, 2.399))
            phi_j = float(a["Muon_phi"][i][j])
            pt_j = muo.muo_pt_scale(cset, is_data, float(a["Muon_pt"][i][j]), eta_j, phi_j,
                                    int(a["Muon_charge"][i][j]))
            if not is_data:
                pt_j = muo.muo_pt_resol(cset, pt_j, eta_j, phi_j,
                                        float(a["Muon_nTrackerLayers"][i][j]), evt[i], lumi[i])
            if pt_j >= VETO_MU_PT_MIN:
                passed[i] = False
                break
    return passed


def rescale_muon_in_sum(a, prefix, k):
    """(pt, eta, phi, mass, Ht) of <prefix> = X + mu after scaling the muon 3-momentum by k.

    Ht: the producer's <prefix>_Ht is a scalar sum of track pt's plus mu_tlv.Pt(), so only the
    muon term moves: Ht + (k - 1) * mu_pt.

    Uses only the stored sum and muon vectors: E_X = E_sum - E_mu, so how the producer built
    X does not matter. Rows where the stored sum is the NULL sentinel stay NULL.
    E_sum is built from the stored px/py/pz, NOT from <prefix>_p: the float32 branches
    disagree by one LSB (|p|/p - 1 ~ 5e-8), which E^2 - p^2 amplifies to a 1e-3 GeV mass
    error at p ~ 950 GeV (measured on 2022EE DsW output_0).
    """
    px, py, pz = (np_(a, f"{prefix}_{c}") for c in ("px", "py", "pz"))
    e_sum = np.sqrt(px ** 2 + py ** 2 + pz ** 2 + np_(a, f"{prefix}_mass") ** 2)
    mpx, mpy, mpz = (np_(a, f"mu_{c}") for c in ("px", "py", "pz"))
    mp = np.sqrt(mpx ** 2 + mpy ** 2 + mpz ** 2)
    e_new = e_sum - np.sqrt(mp ** 2 + MASS_MU ** 2) + np.sqrt((k * mp) ** 2 + MASS_MU ** 2)
    px, py, pz = px + (k - 1) * mpx, py + (k - 1) * mpy, pz + (k - 1) * mpz
    pt = np.hypot(px, py)
    eta = np.arcsinh(np.divide(pz, pt, out=np.zeros_like(pz), where=pt > 0))
    phi = np.arctan2(py, px)
    mass = np.sqrt(np.maximum(e_new ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
    ht = np_(a, f"{prefix}_Ht") + (k - 1) * np_(a, "mu_pt")
    null = np_(a, f"{prefix}_mass") == NULL
    return tuple(np.where(null, NULL, x) for x in (pt, eta, phi, mass, ht))


def w_tmass(mu_pt, mu_phi, met_pt, met_phi):
    """HDsPhiPiMuProducer.cc: sqrt(2 pt_mu MET (1 - cos(phi_mu - phi_MET)))."""
    return np.sqrt(2.0 * mu_pt * met_pt * (1.0 - np.cos(mu_phi - met_phi)))


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(f"usage: {sys.argv[0]} <nano file> [output]")
    in_path = sys.argv[1]
    era, sample, out_path, is_data = parse_path(in_path)
    if len(sys.argv) == 3:
        out_path = os.path.abspath(sys.argv[2])
    campaign = CAMPAIGN[era]
    wp_tight = bem.ERAS[era]["wp_tight"]
    jme_year = bem.ERAS[era]["jme_year"]
    mu_vars = ("nom",) if is_data else MU_VARS
    jme_vars = {"nom": "nominal"} if is_data else JME_VARS
    eff_json = None if is_data else f"{bem.OUT_DIR}/{EFF_GROUP[sample]}_{era}.json"
    print(f"input {in_path}\n {'data' if is_data else 'MC'} era={era} sample={sample}"
          f" JME year={jme_year} WP_T={wp_tight}\n eff map {eff_json}\n output {out_path}")
    if eff_json and not os.path.exists(eff_json):
        sys.exit(f"missing efficiency map {eff_json}")

    branches = [b for b in BRANCHES if not (is_data and b in MC_ONLY_BRANCHES)]
    a = uproot.open(in_path)["Events"].arrays(branches)
    n = len(a)
    print(f" rows: {n}")
    if n == 0:
        sys.exit("empty Events tree -- nothing to correct")
    cols = {b: ak.to_numpy(a[b]) for b in ("run", "luminosityBlock", "event")}

    # ---- closure: the W_tmass formula reproduces the stored value ----
    mu_phi = np_(a, "mu_phi")
    wt = w_tmass(np_(a, "mu_pt"), mu_phi, np_(a, "met_pt"), np_(a, "met_phi"))
    bad = np.abs(wt - np_(a, "W_tmass")) > 1e-3 * np.maximum(np_(a, "W_tmass"), 1.0)
    if bad.any():
        sys.exit(f"W_tmass closure FAILED on {int(bad.sum())}/{n} rows")

    # ---- LUM (JSON flow=clamp, nTrueInt passed as is). MC only. ----
    if not is_data:
        pu = correctionlib.CorrectionSet.from_file(
            f"{POG_DIR}/LUM/{campaign}/2024-01-31/puWeights.json.gz")[PU_KEY[era]]
        nti = np_(a, "Pileup_nTrueInt")
        for suffix, v in (("nom", "nominal"), ("up", "up"), ("down", "down")):
            cols[f"w_pu_{suffix}"] = pu.evaluate(nti, v)

    # ---- MUO ----
    cset_ss = correctionlib.CorrectionSet.from_file(
        f"{POG_DIR}/MUO/{campaign}/2026-06-18/muon_scalesmearing.json.gz")
    cset_z = correctionlib.CorrectionSet.from_file(f"{POG_DIR}/MUO/{campaign}/2026-06-18/muon_Z.json.gz")
    idx = locate_selected_muon(a)
    mu_pt = muon_momenta(a, idx, cset_ss, is_data)
    mu_eta = np_(a, "mu_eta")  # SIGNED: the trigger SF is asymmetric in eta
    iso_tight = ak.to_numpy(ak.flatten(a["Muon_pfIsoId"][idx[:, None]])) >= 4
    pt_raw = np_(a, "mu_pt")

    def mu_sf(key, pt, var, sel):
        # flow="error" keys: evaluated only on selected rows, no clamping
        w = np.full(n, np.nan)
        w[sel] = cset_z[key].evaluate(mu_eta[sel], pt[sel], var)
        return w

    for v in mu_vars:
        sel = (np.abs(mu_eta) < MU_ETA_MAX) & (mu_pt[v] > MU_PT_MIN) & iso_tight
        cols[f"mu_pt_{v}"] = mu_pt[v]
        cols[f"pass_mu_{v}"] = sel
        if not is_data:
            cols[f"w_mu_{v}"] = (mu_sf(muo.MUO_ID_KEY, mu_pt[v], "nominal", sel)
                                 * mu_sf(muo.MUO_ISO_KEY, mu_pt[v], "nominal", sel)
                                 * mu_sf(muo.MUO_TRIG_KEY, mu_pt[v], "nominal", sel))
        k = mu_pt[v] / pt_raw
        for prefix in ("Dsmu", "Dsstmu"):
            (cols[f"{prefix}_pt_{v}"], cols[f"{prefix}_eta_{v}"], cols[f"{prefix}_phi_{v}"],
             cols[f"{prefix}_mass_{v}"], cols[f"{prefix}_Ht_{v}"]) = rescale_muon_in_sum(a, prefix, k)
    sel_nom = cols["pass_mu_nom"]
    if not is_data:
        for label, key in (("muid", muo.MUO_ID_KEY), ("muiso", muo.MUO_ISO_KEY),
                           ("mutrig", muo.MUO_TRIG_KEY)):
            for suffix, v in (("nom", "nominal"), ("up", "systup"), ("down", "systdown")):
                cols[f"w_{label}_{suffix}"] = mu_sf(key, mu_pt["nom"], v, sel_nom)
    cols["pass_muveto"] = prompt_muon_veto(a, idx, cset_ss, is_data)

    # closure of the rescaling at k = 1, against every stored quantity it reproduces
    for prefix in ("Dsmu", "Dsstmu"):
        ok = np_(a, f"{prefix}_mass") != NULL
        if not ok.any():
            continue
        got = dict(zip(("pt", "eta", "phi", "mass", "Ht"), rescale_muon_in_sum(a, prefix, np.ones(n))))
        for q, x in got.items():
            d = np.abs(x - np_(a, f"{prefix}_{q}"))[ok]
            if q == "phi":
                d = np.minimum(d, 2 * np.pi - d)
            if (d > 1e-3).any():
                sys.exit(f"{prefix} k=1 closure FAILED on {q}: max |d| {d.max():.3g}")

    # ---- JME + BTV ----
    refs = jme.JmeRefs(year=jme_year, is_data=is_data,
                       era=DATA_ERA_KEY[era] if is_data else None)
    btv = correctionlib.CorrectionSet.from_file(f"{POG_DIR}/BTV/{campaign}/2025-08-20/btagging.json.gz")
    wp_json = btv["particleNet_wp_values"].evaluate("T")
    if abs(wp_json - wp_tight) > 1e-6:
        sys.exit(f"tight WP mismatch: JSON {wp_json} vs {wp_tight}")
    if not is_data:
        c_comb, c_light = btv["particleNet_comb"], btv["particleNet_light"]
        eff_map = correctionlib.CorrectionSet.from_file(eff_json)["particleNet_efficiency"]

    tight = (a["Jet_jetId"] & 2) == 2
    dphi = abs(a["Jet_phi"] - a["Ds_phi"])
    dphi = ak.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    dR = np.sqrt((a["Jet_eta"] - a["Ds_eta"]) ** 2 + dphi ** 2)
    examined = (abs(a["Jet_eta"]) < 2.5) & (a["Jet_jetId"] >= 2) & (dR > 0.4)
    tagged = examined & (a["Jet_btagPNetB"] > wp_tight)

    def w_btag(pt_j, veto_pass, systs):
        """Pi (1-SF*eps)/(1-eps) over the examined pt>=20 jets (all untagged in veto-passing
        events); SF and eps at the corrected pt of this JME variation."""
        scored = examined & (pt_j >= BTAG_JET_PT_MIN)
        counts = ak.to_numpy(ak.num(pt_j[scored], axis=1))
        pt_s = ak.to_numpy(ak.flatten(pt_j[scored])).astype(float)
        abseta = np.abs(ak.to_numpy(ak.flatten(a["Jet_eta"][scored]))).astype(float)
        flav = ak.to_numpy(ak.flatten(a["Jet_hadronFlavour"][scored])).astype(int)
        eps = eff_map.evaluate("T", flav, abseta, pt_s)
        heavy = (flav == 4) | (flav == 5)
        out = {}
        for s in systs:
            sf = np.empty_like(pt_s)
            sf[heavy] = c_comb.evaluate(s, "T", flav[heavy], abseta[heavy], pt_s[heavy])
            sf[~heavy] = c_light.evaluate(s, "T", flav[~heavy], abseta[~heavy], pt_s[~heavy])
            per_jet = np.where(eps < 1.0, (1.0 - sf * eps) / (1.0 - eps), 1.0)
            w = ak.to_numpy(ak.prod(ak.unflatten(per_jet, counts), axis=1))
            out[s] = np.where(veto_pass, w, np.nan)
        return out

    met_nom = None
    for v, jv in jme_vars.items():
        pt_j = jme.corrected_jet_pt(a, refs, jv)
        cols[f"jvm_veto_{v}"] = ak.to_numpy(jme.jet_veto_flag(a, pt_j, year=jme_year))
        cols[f"jet_ht_pt20_{v}"] = ak.to_numpy(ak.sum(ak.where(tight & (pt_j > 20), pt_j, 0.0), axis=1))
        cols[f"num_jets_pt30_{v}"] = ak.to_numpy(ak.sum(tight & (pt_j > 30), axis=1))
        cols[f"jet_leadpt_{v}"] = ak.to_numpy(ak.fill_none(ak.max(pt_j[tight], axis=1), NULL))
        met_pt, met_phi = (ak.to_numpy(x) for x in jme.corrected_met_pt_phi(a, refs, jv))
        if v == "nom":
            ref = ak.to_numpy(jme.corrected_met_pt(a, refs, jv))
            if not np.array_equal(met_pt, ref):
                sys.exit(f"MET (pt,phi) wrapper disagrees with get_corrected_met:"
                         f" max |d| {np.max(np.abs(met_pt - ref)):.3g}")
            met_nom = (met_pt, met_phi)
        else:
            cols[f"W_tmass_{v}"] = w_tmass(mu_pt["nom"], mu_phi, met_pt, met_phi)
        cols[f"met_pt_{v}"], cols[f"met_phi_{v}"] = met_pt, met_phi
        nb = ak.to_numpy(ak.sum(tagged & (pt_j >= BTAG_JET_PT_MIN), axis=1))
        cols[f"n_bjet_{v}"] = nb
        if not is_data:
            wb = w_btag(pt_j, nb == 0, ("central", "up", "down") if v == "nom" else ("central",))
            cols[f"w_btag_{v}"] = wb["central"]
            if v == "nom":
                cols["w_btag_up"], cols["w_btag_down"] = wb["up"], wb["down"]
    for v in mu_vars:
        cols[f"W_tmass_{v}"] = w_tmass(mu_pt[v], mu_phi, *met_nom)

    # ---- summary ----
    print(" summary (unweighted; finite values only):")
    for k in sorted(cols):
        x = np.asarray(cols[k])
        if k in ("run", "luminosityBlock", "event"):
            continue
        if x.dtype == bool:
            print(f"   {k:26s} frac True {x.mean():.4f}")
        else:
            f = np.isfinite(x)
            print(f"   {k:26s} mean {np.mean(x[f]):12.5g}  min {np.min(x[f]):12.5g}"
                  f"  max {np.max(x[f]):12.5g}  NaN {int((~f).sum())}")

    # ---- write (local temp, then copy), then read back the keys ----
    for k, x in cols.items():
        if len(x) != n:
            sys.exit(f"column {k} has {len(x)} rows, expected {n}")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        local = os.path.join(tmp, "friend.root")
        with uproot.recreate(local) as f:
            f["Events"] = {k: np.asarray(x) for k, x in cols.items()}
        shutil.copyfile(local, out_path)
    with uproot.open(out_path) as f:
        t = f["Events"]
        back = t.arrays(["run", "luminosityBlock", "event"], library="np")
        if t.num_entries != n or any(not np.array_equal(back[b], cols[b]) for b in back):
            sys.exit(f"read-back of {out_path} does not match the input keys")
    print(f"Wrote {out_path} ({n} rows, {len(cols)} columns)")


if __name__ == "__main__":
    main()

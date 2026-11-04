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
project_corrections_pilot_implementation). One module per POG, each calling the POG's own
example code where there is one (external/, unmodified):
  lum.py  pileup weights
  muo.py  muonscarekit MuonScaRe.cc (via ROOT) + muon SFs, selected muon, prompt-muon veto
  btv.py  b-jet count + btag_weight.cc (C++ port of BTV examples compute_btag_event_weight_multi_wp)
  jme.py  jet veto map + nano jet/MET columns
This file reads the nano, calls them, runs the closures and writes the friend.

NO JME CORRECTION (2026-10-21, user decision). In the Z->mumu CR the Type-1 + Summer22 V4
JEC + MC JER chain made the MC MET harder than data while the nano MET agreed. Jets and MET
are therefore the nano values: Jet_pt (MiniAOD updatedJetsPuppi, i.e. the production GT's
JEC, no JER smearing in MC) and met_pt/met_phi (slimmedMETsPuppi). The column NAMES are
unchanged (jet_ht_pt20_nom, met_pt_nom, n_bjet_nom, ...) so make_samples_nano.py reads them
as before; their *_corr names in the merged samples now carry only the muon correction.
The jet veto map is kept (a map, not an energy correction), evaluated on nano Jet_pt.
The b-tag efficiency maps must be built on nano Jet_pt too (build_btag_efficiency_map.py).
The JME version is in git history (commit ea09e9c and earlier).

VARIATIONS. Suffixes: _nom; muon momentum _muScaleUp/_muScaleDown/_muResUp/_muResDown
(no JES/JER variations any more). Each muon variation re-derives:
  mu_pt, pass_mu, w_mu (id x iso x trig, nominal SFs at that pt), W_tmass (with the nano
  MET), Dsmu_{pt,eta,phi,mass,Ht}, Dsstmu_{pt,eta,phi,mass,Ht}
  (mu_px/py/pz follow from the stored direction times mu_pt; mu_IsoR0x_RelIso is
  left as stored, since the POG iso WP and its SF are defined on the uncorrected pt)
Written once (_nom, nano jets/MET): jet_ht_pt20, num_jets_pt30, jet_leadpt, met_pt, met_phi,
jvm_veto, n_bjet, w_btag.
SF systematics at the nominal kinematics: w_muid/w_muiso/w_mutrig _nom/_up/_down,
w_btag_up/_down, w_btag_any_nom/_up/_down, w_pu_nom/_up/_down.

COLUMN SEMANTICS
  pass_mu      selected muon: corrected pt > 26, |eta| < 2.4, PFIsoTight (Muon_pfIsoId>=4)
  pass_muveto  no other muon with tight ID, |eta|<2.4, Loose PF iso, corrected pt >= 15
               (nominal muon correction only)
  jvm_veto     True = event VETOED by the jet veto map (official apply_jvm_flag)
  n_bjet       jets with tight ID, |eta|<2.5, dR(Ds)>0.4, PNet B > era tight WP,
               nano pt >= 20; the b-veto is n_bjet == 0
  w_mu*        NaN where pass_mu (of that variation) is False
  w_btag*      NaN where n_bjet (of that variation) > 0
  w_btag_any_* (nom/up/down only) valid at any n_bjet: tagged jets give SF, untagged (1-SF*eps)/(1-eps)
NaN, not 1, so a weight used without its cut poisons the sum instead of passing silently.

NOT WRITTEN: met_sig and met_sumEt (use the nano values). jet_ht and num_jets cannot be
rebuilt (nano stores jets with pt > 15 only); jet_ht_pt20 replaces jet_ht.
The producer's candidate choice and its Dsmu_mass < 125 cut used the uncorrected muon.

Runs inside cmssw-el9 + CMSSW_15_0_5 (correctionlib 2.6.4, ROOT for MuonScaRe.cc).
"""
import os
import shutil
import sys
import tempfile

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_btag_efficiency_map as bem  # noqa: E402  (ERAS: jme_year, wp_tight; OUT_DIR)
import btv  # noqa: E402
import jme  # noqa: E402
import lum  # noqa: E402
import muo  # noqa: E402

CAMPAIGN = {
    "2022": "Run3-22CDSep23-Summer22-NanoAODv12",
    "2022EE": "Run3-22EFGSep23-Summer22EE-NanoAODv12",
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
    "WtoLNu2Jets_0J_H_phipimu_GenTag_SS": "WtoLNu2Jets_0J_SS",
    "WtoLNu2Jets_1J_H_phipimu_GenTag_SS": "WtoLNu2Jets_1J_SS",
    "WtoLNu2Jets_2J_H_phipimu_GenTag_SS": "WtoLNu2Jets_2J_SS",
    "TTto2L2Nu_H_phipimu_GenTag": "TTto2L2Nu",
    "TTtoLNu2Q_H_phipimu_GenTag": "TTtoLNu2Q",
    # ttbar _ext1 shares the nominal map (one pooled map, see build_btag_efficiency_map.py)
    "TTto2L2Nu_ext1_H_phipimu_GenTag": "TTto2L2Nu",
    "TTtoLNu2Q_ext1_H_phipimu_GenTag": "TTtoLNu2Q",
    "TTto2L2Nu_H_phipimu_GenTag_SS": "TTto2L2Nu_SS",
    "TTtoLNu2Q_H_phipimu_GenTag_SS": "TTtoLNu2Q_SS",
    "TTto2L2Nu_ext1_H_phipimu_GenTag_SS": "TTto2L2Nu_SS",
    "TTtoLNu2Q_ext1_H_phipimu_GenTag_SS": "TTtoLNu2Q_SS",
    "DYto2L2Jets_0J_H_phipimu_GenTag": "DYto2L2Jets_0J",
    "DYto2L2Jets_1J_H_phipimu_GenTag": "DYto2L2Jets_1J",
    "DYto2L2Jets_2J_H_phipimu_GenTag": "DYto2L2Jets_2J",
    "DYto2L2Jets_0J_H_phipimu_GenTag_SS": "DYto2L2Jets_0J_SS",
    "DYto2L2Jets_1J_H_phipimu_GenTag_SS": "DYto2L2Jets_1J_SS",
    "DYto2L2Jets_2J_H_phipimu_GenTag_SS": "DYto2L2Jets_2J_SS",
}

MASS_MU = 0.1056583755
NULL = -777.0

# DATA. Recognised by the path. Differences from MC, all of them because the quantity does
# not exist for data rather than by choice: the muon gets the scale correction from the data
# branch of MuonScaRe and no resolution smearing, and no weight is defined (no PU, no scale
# factors, no b-tag weight). The selection flags are computed exactly as for MC, because the
# same cuts are applied to data. Jets/MET are the nano values for data as for MC.

MC_ONLY_BRANCHES = {"Pileup_nTrueInt"} | set(btv.MC_BRANCHES)

BRANCHES = sorted(set([
    "run", "luminosityBlock", "event", "Pileup_nTrueInt",
    "mu_pt", "mu_eta", "mu_phi", "mu_charge", "mu_px", "mu_py", "mu_pz",
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_nTrackerLayers",
    "Muon_pfIsoId",
    "met_pt", "met_phi", "W_tmass",
    "Dsmu_px", "Dsmu_py", "Dsmu_pz", "Dsmu_mass", "Dsmu_pt", "Dsmu_eta", "Dsmu_phi", "Dsmu_Ht",
    "Dsstmu_px", "Dsstmu_py", "Dsstmu_pz", "Dsstmu_mass", "Dsstmu_pt", "Dsstmu_eta", "Dsstmu_phi",
    "Dsstmu_Ht",
] + btv.BRANCHES + btv.MC_BRANCHES + jme.BRANCHES))


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
    mu_vars = ("nom",) if is_data else muo.MU_VARS
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

    # ---- LUM (MC only) ----
    if not is_data:
        for suffix, w in lum.pu_weights(campaign, era, np_(a, "Pileup_nTrueInt")).items():
            cols[f"w_pu_{suffix}"] = w

    # ---- MUO ----
    muo.load_scalesmearing(campaign)
    cset_z = muo.sf_set(campaign)
    idx = muo.locate_selected_muon(a)
    mu_pt = muo.muon_momenta(a, idx, is_data)
    mu_eta = np_(a, "mu_eta")  # SIGNED, as in the MUO example
    iso_tight = ak.to_numpy(ak.flatten(a["Muon_pfIsoId"][idx[:, None]])) >= 4
    pt_raw = np_(a, "mu_pt")
    for v in mu_vars:
        sel = (np.abs(mu_eta) < muo.MU_ETA_MAX) & (mu_pt[v] > muo.MU_PT_MIN) & iso_tight
        cols[f"mu_pt_{v}"] = mu_pt[v]
        cols[f"pass_mu_{v}"] = sel
        if not is_data:
            cols[f"w_mu_{v}"] = (muo.sf(cset_z, muo.MUO_ID_KEY, mu_eta, mu_pt[v], "nominal", sel)
                                 * muo.sf(cset_z, muo.MUO_ISO_KEY, mu_eta, mu_pt[v], "nominal", sel)
                                 * muo.sf(cset_z, muo.MUO_TRIG_KEY, mu_eta, mu_pt[v], "nominal", sel))
        k = mu_pt[v] / pt_raw
        for prefix in ("Dsmu", "Dsstmu"):
            (cols[f"{prefix}_pt_{v}"], cols[f"{prefix}_eta_{v}"], cols[f"{prefix}_phi_{v}"],
             cols[f"{prefix}_mass_{v}"], cols[f"{prefix}_Ht_{v}"]) = rescale_muon_in_sum(a, prefix, k)
    sel_nom = cols["pass_mu_nom"]
    if not is_data:
        for label, key in (("muid", muo.MUO_ID_KEY), ("muiso", muo.MUO_ISO_KEY),
                           ("mutrig", muo.MUO_TRIG_KEY)):
            for suffix, v in (("nom", "nominal"), ("up", "systup"), ("down", "systdown")):
                cols[f"w_{label}_{suffix}"] = muo.sf(cset_z, key, mu_eta, mu_pt["nom"], v, sel_nom)
    cols["pass_muveto"] = muo.prompt_muon_veto(a, idx, is_data)

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

    # ---- JME (jet veto map, nano jets/MET) ----
    cols.update(jme.columns(a, jme_year, NULL))

    # ---- BTV ----
    btv_set = btv.sf_set(campaign, wp_tight)
    nb = btv.n_bjet(a, wp_tight)
    cols["n_bjet_nom"] = nb
    if not is_data:
        cols.update(btv.weights(a, campaign, btv_set, eff_json, nb))
    for v in mu_vars:
        cols[f"W_tmass_{v}"] = w_tmass(mu_pt[v], mu_phi, cols["met_pt_nom"], cols["met_phi_nom"])

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

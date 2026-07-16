"""Score every Samples/ file with a trained BDT and write a copy with score_dsmu added.

All 36 samples are scored (2 signals + 34 backgrounds), not just the ones used in
training: the diagnostics need the per-process survival rates and the mass spectra
of the backgrounds the MVA never saw.

Rows whose features are non-finite get score_dsmu = -777.

Usage:  python3 predict.py <model_name>
"""

import sys
import os
import numpy as np
import uproot
import joblib

model_name  = sys.argv[1]
model_dir   = f"/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/{model_name}"
SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUT_DIR     = f"{model_dir}/Prediction"
TREE        = "Hphipimu/Events"

os.makedirs(OUT_DIR, exist_ok=True)

# Must stay identical to BDT_phipimu.py.
FEATURE_NAMES = [
    # lifetime / vertex -- prompt-from-H vs displaced-from-B (the real-Ds background)
    "Ds_ipsig",
    "Ds_FDsig",
    "dira",
    "Ds_prob",
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "mu_ipsig",
    # 2D lifetime/vertex: transverse resolution beats longitudinal.
    # FDxysig separates DY (fake Ds from prompt tracks -> no flight distance);
    # ipxysig separates TT (real Ds from a displaced B -> momentum misses the PV).
    "Ds_FDxysig",
    "Ds_ipxysig",
    # Are the Ds tracks from the primary vertex? (DY's fakes are prompt tracks)
    "Kp_PVassoc",
    "Km_PVassoc",
    "pi_PVassoc",
    # isolation -- the components, not just RelIso. sumChargedHadronPt is the single
    # strongest variable against every background (AUC 0.91 vs TT): a Ds buried in a
    # b jet is surrounded by charged hadrons. RelIso divides by Ds_pt, which dilutes it.
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "Ds_IsoR04_RelIso",
    "Ds_IsoR03_RelIso",
    "mu_IsoR04_RelIso",
    # MET / event. num_jets_pt20..pt60 are near-duplicates (AUC all ~0.80), so only
    # one threshold is kept; num_jets (no threshold) adds the soft-jet activity.
    "met_pt",
    "met_sig",
    "met_sumEt",
    "W_tmass",
    "num_jets",
    "num_jets_pt30",
    # Ds + mu kinematics -- hard from a 125 GeV H vs soft inside a jet
    "Ds_pt",
    "Dsmu_pt",
    "Dsmu_eta",
    "mu_pt",
    "mu_eta",
    "dR_Ds_mu",
    "dR_Kp_mu",
    "dR_Km_mu",
    "dR_pi_mu",
    "dR_phi_mu",
    # Ds+mu mass: H-mass constrained for signal, continuum for background
    "Dsmu_mass",
]

SAMPLES = [
    "HtoDsW", "HtoDsstW",
    "WtoLNu_2Jets", "TTtoLNu2Q", "WGtoLNuG_1Jets", "TQbartoLNu", "TbarQtoLNu",
    "WWtoLNu2Q", "VBFtoLNu", "TbarWplustoLNu2Q", "TWminustoLNu2Q", "WZtoLNu2Q",
    "WWGtoLNu2QG", "TTZ_ZtoQQ_TTtoLNu2Q", "WZto3LNu", "WZtoL3Nu",
    "TTto2L2Nu", "TbarWplusto2L2Nu", "TWminusto2L2Nu", "WWto2L2Nu",
    "ggWWto2L2Nu_LL", "ggWWto2L2Nu_LT", "ggWWto2L2Nu_TL", "ggWWto2L2Nu_TT",
    "ZZto2L2Nu", "DYto2L_MLL4to10", "DYto2L_MLL10to50", "DYto2L_MLL50",
    "DYGto2LG_MLL4to50_PTG10to100", "DYGto2LG_MLL4to50_PTG100to200",
    "DYGto2LG_MLL4to50_PTG200", "DYGto2LG_MLL50", "VBFto2L_MLL50",
    "WZto2L2Q", "ZZto2L2Q", "ZGto2LG",
]

dira_idx = FEATURE_NAMES.index("dira")

print(f"Loading model from {model_dir}...")
model = joblib.load(f"{model_dir}/bdt_model.pkl")

for name in SAMPLES:
    with uproot.open(f"{SAMPLES_DIR}/{name}.root")[TREE] as t:
        data = t.arrays(library="np")

    X = np.column_stack([data[b].astype(np.float64) for b in FEATURE_NAMES])
    n = len(X)

    mask = np.isfinite(X).all(axis=1) & (X[:, dira_idx] > -1)
    scores = np.full(n, -777.0, dtype=np.float32)
    scores[mask] = model.predict_proba(X[mask])[:, 1].astype(np.float32)

    data["score_dsmu"] = scores
    out_file = f"{OUT_DIR}/{name}.root"
    with uproot.recreate(out_file) as f:
        f.mktree(TREE, {k: v.dtype for k, v in data.items()})
        f[TREE].extend(data)

    print(f"  {name:32s} {n:>8,} events, {int((~mask).sum()):>4} non-finite -> {out_file}")

print("Done.")

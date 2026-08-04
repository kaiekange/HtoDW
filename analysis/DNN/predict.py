"""Score every current Samples/ file with the trained DNN and write a copy with
score_dsmu added.

Both OS and SS variants are scored, and data as well as MC: the diagnostics need
the score distribution in the SS control region and in data sidebands, not just in
the samples the MVA was trained on.

The event selection (1 tight muon + b-jet veto) is deliberately NOT applied here --
every row gets a score and the selection is applied downstream, so a cut can be
revisited without re-scoring. Rows with non-finite features get score_dsmu = -777.

The preprocessing chain below (sentinel fix -> log transform -> scaler) must match
DNN.py exactly, in that order, or the scores are meaningless.

Usage:  python3 predict.py
"""

import os
import numpy as np
import uproot
import torch
import torch.nn as nn

MODEL_DIR   = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_wjets_test"
SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUT_DIR     = f"{MODEL_DIR}/Prediction"
TREE        = "Events"

os.makedirs(OUT_DIR, exist_ok=True)

# Must stay identical to DNN.py (and to BDT/BDT.py + BDT/predict.py).
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
    "Ds_FDxysig",
    "Ds_ipxysig",
    # Are the Ds tracks from the primary vertex?
    "Kp_PVassoc",
    "Km_PVassoc",
    "pi_PVassoc",
    # isolation -- the components, not just RelIso.
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "Ds_IsoR04_RelIso",
    "Ds_IsoR03_RelIso",
    "mu_IsoR04_RelIso",
    # MET / event
    "met_pt",
    "met_sig",
    "met_sumEt",
    "W_tmass",
    "num_jets",
    "num_jets_pt30",
    # Ds + mu kinematics
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
    "Dsmu_mass",
    # jet/b-tag environment
    "Ds_nearestjet_bscore",
    "Ds_nearestjet_dR",
    "jet_leadpt",
    "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR04",
]

LOG_FEATURES = [
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "Ds_ipsig",
    "mu_ipsig",
    "Ds_FDsig",
    "Ds_FDxysig",
    "Ds_ipxysig",
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
]

# Must stay identical to DNN.py -- see the comment there for what each sentinel means.
SENTINEL_FIX = {
    "Ds_nearestjet_bscore": -0.1,
    "Ds_nearestjet_dR": 5.0,
    "jet_leadpt": 0.0,
}

# Current Samples/ contents. The old inclusive WtoLNu2Jets{,_SS} are deliberately
# left out -- superseded by the 0J/1J/2J jet-binned samples.
SAMPLES = [
    "HtoDsW", "HtoDsW_SS",
    "WtoLNu2Jets_0J", "WtoLNu2Jets_1J", "WtoLNu2Jets_2J",
    "WtoLNu2Jets_0J_SS", "WtoLNu2Jets_1J_SS", "WtoLNu2Jets_2J_SS",
    "Data_2022", "Data_2022_SS",
]

dira_idx = FEATURE_NAMES.index("dira")
log_idx  = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]


class phipimuNet(nn.Module):
    def __init__(self, n_in, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64),  nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 32),   nn.BatchNorm1d(32),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, 16),   nn.BatchNorm1d(16),  nn.ReLU(),
            nn.Linear(16, 1),
        )

    def forward(self, x):
        return self.net(x)


print(f"Loading model from {MODEL_DIR}...")
model = phipimuNet(len(FEATURE_NAMES))
model.load_state_dict(torch.load(f"{MODEL_DIR}/phipimu_model.pt", map_location="cpu"))
model.eval()

scaler_mean  = np.load(f"{MODEL_DIR}/scaler_mean.npy")
scaler_scale = np.load(f"{MODEL_DIR}/scaler_scale.npy")

for name in SAMPLES:
    path = f"{SAMPLES_DIR}/{name}.root"
    if not os.path.exists(path):
        print(f"  {name:24s} MISSING -- skipped")
        continue

    with uproot.open(path)[TREE] as t:
        data = t.arrays(library="np")

    X = np.column_stack([data[b].astype(np.float64) for b in FEATURE_NAMES])
    n = len(X)

    # Sentinels first, exactly as in training -- before the finite mask, the log
    # transform and the scaler. Skipping this would push -777 through the scaler
    # and produce garbage scores for those rows.
    for feat, fill in SENTINEL_FIX.items():
        i = FEATURE_NAMES.index(feat)
        X[X[:, i] < -1, i] = fill

    mask = np.isfinite(X).all(axis=1) & (X[:, dira_idx] > -1)
    scores = np.full(n, -777.0, dtype=np.float32)

    Xp = X[mask]
    Xp[:, log_idx]  = np.log1p(Xp[:, log_idx])
    Xp[:, dira_idx] = np.log(1.01 - Xp[:, dira_idx])
    Xp = (Xp - scaler_mean) / scaler_scale

    with torch.no_grad():
        s = torch.sigmoid(model(torch.tensor(Xp, dtype=torch.float32))).numpy().flatten()
    scores[mask] = s.astype(np.float32)

    data["score_dsmu"] = scores
    out_file = f"{OUT_DIR}/{name}.root"
    with uproot.recreate(out_file) as f:
        f.mktree(TREE, {k: v.dtype for k, v in data.items()})
        f[TREE].extend(data)

    print(f"  {name:24s} {n:>9,} events, {int((~mask).sum()):>5} non-finite -> {out_file}")

print("Done.")

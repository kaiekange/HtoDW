"""Apply a trained DsmuNet model to a ROOT file and write a new file with a 'score_dsmu' branch added.

Usage:
    python3 predict.py <model_name> <input_root> <tree_name> [output_root]

Outputs:
    <output_root>  (default: input file with _scored.root suffix)
      copy of the input tree with an extra float32 branch 'score_dsmu' in [0,1].
      Entries that fail the finite filter get score_dsmu = -777.
"""

import sys
import os
import numpy as np
import uproot
import torch
import torch.nn as nn

# ========= config =========
model_name = sys.argv[1]
input_file = sys.argv[2]
tree_name  = sys.argv[3]
out_root   = sys.argv[4] if len(sys.argv) > 4 else input_file.replace(".root", "_scored.root")

model_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/{model_name}"

FEATURE_NAMES = [
    # MET / event
    "met_sig",
    "met_pt",
    "met_sumEt",
    "W_tmass",
    "num_jets_pt30",
    # muon
    "mu_eta",
    "mu_pt",
    # K+
    "Kp_eta",
    "Kp_pt",
    "Kp_PVassoc",
    # K-
    "Km_eta",
    "Km_pt",
    "Km_PVassoc",
    # pi
    "pi_eta",
    "pi_pt",
    "pi_PVassoc",
    # phi
    "phi_eta",
    "phi_pt",
    # Ds
    "Ds_eta",
    "Ds_pt",
    # Ds+mu system
    "Dsmu_eta",
    "Dsmu_pt",
    # Ds fit quality
    "Ds_prob",
    "dira",
    # dR (KKpi subsystem)
    "dR_Kp_Km",
    "dR_Kp_pi",
    "dR_Kp_Ds",
    "dR_Km_pi",
    "dR_Km_Ds",
    "dR_pi_Ds",
    "dR_phi_pi",
    "dR_phi_Ds",
    # dR (muon)
    "dR_Kp_mu",
    "dR_Km_mu",
    "dR_pi_mu",
    "dR_phi_mu",
    "dR_Ds_mu",
    # IP significances (3D)
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "Ds_ipsig",
    "mu_ipsig",
    # Ds flight distance (3D)
    "Ds_FDsig",
    # isolation
    "mu_IsoR04_RelIso",
    "Ds_IsoR04_RelIso",
]

LOG_FEATURES = [
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "Ds_ipsig",
    "mu_ipsig",
    "Ds_FDsig",
]

log_idx  = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]
dira_idx = FEATURE_NAMES.index("dira")


# ========= model definition (must match DNN_phipimu.py) =========
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


# ========= load artifacts =========
print(f"Loading model from {model_dir}...")
scaler_mean  = np.load(f"{model_dir}/scaler_mean.npy")
scaler_scale = np.load(f"{model_dir}/scaler_scale.npy")

model = phipimuNet(len(FEATURE_NAMES))
model.load_state_dict(torch.load(f"{model_dir}/phipimu_model.pt", map_location="cpu"))
model.eval()

# ========= load data =========
print(f"Reading {input_file} [{tree_name}]...")
with uproot.open(input_file) as f:
    all_data = f[tree_name].arrays(library="np")

X = np.column_stack([all_data[b] for b in FEATURE_NAMES])
n_entries = len(X)
print(f"Loaded {n_entries} entries")

# ========= preprocessing (must match DNN_Dsmu.py exactly) =========
mask = np.isfinite(X).all(axis=1)
mask &= X[:, dira_idx] > -1
n_dropped = int((~mask).sum())
if n_dropped:
    print(f"Dropping {n_dropped} non-finite entries (score_dsmu = -777)")
X_pass = X[mask].copy()

X_pass[:, log_idx]  = np.log1p(X_pass[:, log_idx])
X_pass[:, dira_idx] = np.log(1.01 - X_pass[:, dira_idx])
X_pass = (X_pass - scaler_mean) / scaler_scale

# ========= inference =========
print("Running inference...")
with torch.no_grad():
    scores_pass = torch.sigmoid(
        model(torch.tensor(X_pass, dtype=torch.float32))
    ).numpy().flatten()

scores_full = np.full(n_entries, -777.0, dtype=np.float32)
scores_full[mask] = scores_pass

passing = scores_full[mask]
print(f"Score stats: min={passing.min():.4f}  max={passing.max():.4f}  mean={passing.mean():.4f}")

# ========= write output ROOT file =========
print(f"Writing {out_root}...")
out_branches = dict(all_data)
out_branches["score_dsmu"] = scores_full

branch_types = {k: v.dtype for k, v in out_branches.items()}
with uproot.recreate(out_root) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print(f"Done. {n_entries} entries with 'score_dsmu' branch → {out_root}")

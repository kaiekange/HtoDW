"""Apply a trained GammaNet model to a ROOT file and write a new file with
a 'score_gamma' branch added.

Usage:
    python3 predict.py <model_name> <input_root> <tree_name> [output_root]

Entries that fail the finite filter get score_gamma = -777.
"""

import sys
import numpy as np
import uproot
import torch
import torch.nn as nn

model_name = sys.argv[1]
input_file = sys.argv[2]
tree_name  = sys.argv[3]
out_root   = sys.argv[4] if len(sys.argv) > 4 else input_file.replace(".root", "_scored.root")

model_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/{model_name}"

FEATURE_NAMES = [
    "gamma_pt",
    "dR_Ds_gamma",
    "dR_Kp_gamma",
    "dR_Km_gamma",
    "dR_pi_gamma",
    "dR_gamma_mu",
    "dR_Dsst_mu",
    "Dsst_pt",
    "Dsst_p",
    "Dsst_IsoR04_RelIso",
    "Dsst_IsoR04_sumPhotonEt",
    "Dsstmu_pt",
    "Dsstmu_eta",
]
ALL_FEATURES = FEATURE_NAMES

LOG_FEATURES = []
log_idx = [ALL_FEATURES.index(f) for f in LOG_FEATURES]


class GammaNet(nn.Module):
    def __init__(self, n_in, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, 32),  nn.BatchNorm1d(32),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, 16),   nn.BatchNorm1d(16),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(16, 8),    nn.BatchNorm1d(8),   nn.ReLU(),
            nn.Linear(8, 1),
        )
    def forward(self, x):
        return self.net(x)


# ========= load artifacts =========
print(f"Loading model from {model_dir}...")
scaler_mean  = np.load(f"{model_dir}/scaler_mean.npy")
scaler_scale = np.load(f"{model_dir}/scaler_scale.npy")

model = GammaNet(len(ALL_FEATURES))
model.load_state_dict(torch.load(f"{model_dir}/gamma_model.pt", map_location="cpu"))
model.eval()

# ========= load data =========
print(f"Reading {input_file} [{tree_name}]...")
with uproot.open(input_file)[tree_name] as t:
    all_data = t.arrays(library="np")

X  = np.column_stack([all_data[b] for b in FEATURE_NAMES])
n_entries = len(X)
print(f"Loaded {n_entries} entries")

# ========= preprocessing =========
mask = np.isfinite(X).all(axis=1)
n_dropped = int((~mask).sum())
if n_dropped:
    print(f"Dropping {n_dropped} non-finite entries (score_gamma = -777)")
X_pass = X[mask].copy()

X_pass[:, log_idx] = np.log1p(X_pass[:, log_idx])
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

# ========= write output =========
print(f"Writing {out_root}...")
out_branches = dict(all_data)
out_branches["score_gamma"] = scores_full

branch_types = {k: v.dtype for k, v in out_branches.items()}
with uproot.recreate(out_root) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print(f"Done. {n_entries} entries with 'score_gamma' branch → {out_root}")

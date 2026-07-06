"""Apply a trained XGBoost BDT model to a ROOT file and write a new file with 'score_dsmu' added.

Usage:
    python3 predict.py <model_name> <input_root> <tree_name> [output_root]
"""

import sys
import os
import numpy as np
import uproot
import joblib

model_name = sys.argv[1]
input_file = sys.argv[2]
tree_name  = sys.argv[3]
out_root   = sys.argv[4] if len(sys.argv) > 4 else input_file.replace(".root", "_scored.root")

model_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/{model_name}"

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

dira_idx = FEATURE_NAMES.index("dira")

print(f"Loading model from {model_dir}...")
model = joblib.load(f"{model_dir}/bdt_model.pkl")

print(f"Reading {input_file} [{tree_name}]...")
with uproot.open(input_file) as f:
    all_data = f[tree_name].arrays(library="np")

X = np.column_stack([all_data[b] for b in FEATURE_NAMES])
n_entries = len(X)
print(f"Loaded {n_entries} entries")

mask = np.isfinite(X).all(axis=1)
mask &= X[:, dira_idx] > -1
n_dropped = int((~mask).sum())
if n_dropped:
    print(f"Dropping {n_dropped} non-finite entries (score_dsmu = -777)")

print("Running inference...")
scores_pass = model.predict_proba(X[mask])[:, 1].astype(np.float32)

scores_full = np.full(n_entries, -777.0, dtype=np.float32)
scores_full[mask] = scores_pass

passing = scores_full[mask]
print(f"Score stats: min={passing.min():.4f}  max={passing.max():.4f}  mean={passing.mean():.4f}")

print(f"Writing {out_root}...")
out_branches = dict(all_data)
out_branches["score_dsmu"] = scores_full

branch_types = {k: v.dtype for k, v in out_branches.items()}
with uproot.recreate(out_root) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print(f"Done. {n_entries} entries with 'score_dsmu' branch → {out_root}")

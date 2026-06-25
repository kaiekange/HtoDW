"""Redraw all DNN_gamma figures from saved artifacts — no retraining needed.

Usage:
    python3 plot.py <model_name>
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import roc_curve, roc_auc_score

import mplhep as hep
plt.style.use(hep.style.CMS)

model_name = sys.argv[1]
out_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/{model_name}"

ALL_FEATURES = [
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
    "Dsstmu_p",
    "dM",
]

# ========= Load artifacts =========
scores       = np.load(f"{out_dir}/scores.npy")
y_test       = np.load(f"{out_dir}/y_test.npy")
train_losses = np.load(f"{out_dir}/train_losses.npy")
val_losses   = np.load(f"{out_dir}/val_losses.npy")
val_aucs     = np.load(f"{out_dir}/val_aucs.npy")
importances  = np.load(f"{out_dir}/feature_importances.npy")
X_sig        = np.load(f"{out_dir}/X_sig.npy")
X_bkg        = np.load(f"{out_dir}/X_bkg.npy")

auc = roc_auc_score(y_test, scores)
print(f"Test AUC = {auc:.4f}")

# ========= Variable comparison =========
print("Drawing variable_comparison.pdf...")
N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS

with PdfPages(f"{out_dir}/variable_comparison.pdf") as pdf:
    for page_start in range(0, len(ALL_FEATURES), PER_PAGE):
        page_feats = ALL_FEATURES[page_start:page_start + PER_PAGE]
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(25, 20))
        axes_flat = axes.flatten()

        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            fi = ALL_FEATURES.index(feat)
            s_vals = X_sig[:, fi]
            b_vals = X_bkg[:, fi]
            all_vals = np.concatenate([s_vals, b_vals])
            lo, hi = np.percentile(all_vals, [1, 99])
            if lo == hi:
                lo, hi = all_vals.min(), all_vals.max()
            bins = np.linspace(lo, hi, 50)
            ax.hist(b_vals, bins=bins, density=True, alpha=0.6, color="steelblue", label="Fake gamma")
            ax.hist(s_vals, bins=bins, density=True, alpha=0.6, color="firebrick",  label="True gamma")
            ax.set_xlabel(feat)
            ax.set_ylabel("a.u.")
            ax.legend(fontsize=13)
            if i == 0:
                hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

        for j in range(len(page_feats), len(axes_flat)):
            axes_flat[j].set_visible(False)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print("Saved variable_comparison.pdf")

# ========= Training results =========
print("Drawing results.png...")
fig, axes = plt.subplots(2, 3, figsize=(21, 12))
axes = axes.flatten()

axes[0].plot(train_losses, label="Train Loss")
axes[0].plot(val_losses, label="Val Loss")
axes[0].set_xlabel("Epoch"); axes[0].set_ylabel("Loss"); axes[0].legend()
hep.cms.label("Simulation", ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

axes[1].plot(val_aucs, color="steelblue")
axes[1].set_xlabel("Epoch"); axes[1].set_ylabel("Val AUC")
hep.cms.label("Simulation", ax=axes[1], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

fpr, tpr, _ = roc_curve(y_test, scores)
axes[2].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[2].plot([0, 1], [0, 1], "k--")
axes[2].set_xlabel("False Positive Rate"); axes[2].set_ylabel("True Positive Rate"); axes[2].legend()
hep.cms.label("Simulation", ax=axes[2], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

sig_scores = scores[y_test == 1]
bkg_scores = scores[y_test == 0]
bins = np.linspace(0, 1, 51)

axes[3].hist(bkg_scores, bins=bins, alpha=0.6, color="steelblue", label="Fake gamma", density=True)
axes[3].hist(sig_scores, bins=bins, alpha=0.6, color="firebrick",  label="True gamma", density=True)
axes[3].set_xlabel("DNN Score"); axes[3].set_ylabel("a.u."); axes[3].legend()
hep.cms.label("Simulation", ax=axes[3], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

axes[4].hist(bkg_scores, bins=bins, alpha=0.6, color="steelblue", label="Fake gamma", density=True)
axes[4].hist(sig_scores, bins=bins, alpha=0.6, color="firebrick",  label="True gamma", density=True)
axes[4].set_xlabel("DNN Score"); axes[4].set_ylabel("a.u."); axes[4].set_yscale("log"); axes[4].legend()
hep.cms.label("Simulation", ax=axes[4], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

axes[5].set_visible(False)

plt.tight_layout()
plt.savefig(f"{out_dir}/results.png", dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= Feature importance =========
print("Drawing feature_importance.png...")
sorted_idx = np.argsort(importances)
colors = ["firebrick" if importances[i] > 0 else "steelblue" for i in sorted_idx]

fig, ax = plt.subplots(figsize=(8, len(ALL_FEATURES) * 0.28))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color=colors)
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([ALL_FEATURES[i] for i in sorted_idx], fontsize=13)
ax.set_xlabel("AUC drop when feature permuted")
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

print("Done.")

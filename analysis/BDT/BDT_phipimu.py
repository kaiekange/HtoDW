"""
XGBoost BDT: true Ds+mu vs combinatorial (W+jets), phipimu channel.
Signal:     Signal.root (output of Prepare_signal.cc) — DsW+DsstW, Ds_match&&mu_match, one per event
Background: WJets_Background.root (output of Prepare.cc) — W+jets, one per event

Run both Prepare macros first:
  root -l -b -q Prepare_signal.cc
  root -l -b -q Prepare.cc

Train/val/test split by event_id to avoid leakage.
"""
import sys
import numpy as np
import uproot
import xgboost as xgb
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages

import mplhep as hep
plt.style.use(hep.style.CMS)

from sklearn.metrics import roc_curve, roc_auc_score
import joblib

# ========= 1. 配置 =========
SIG_FILE   = "/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/Signal.root"
WJETS_FILE = "/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/WJets_Background.root"

TREE = "Events"

model_name = sys.argv[1] if len(sys.argv) > 1 else "v1"
out_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/{model_name}"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)

# ========= 2. 特征定义 =========
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

LOAD_BRANCHES = FEATURE_NAMES + ["event_id"]

# ========= 3. 加载数据 =========
def load_file(fname):
    with uproot.open(fname)[TREE] as t:
        return t.arrays(LOAD_BRANCHES, library="np")

def to_matrix(arrs):
    return np.column_stack([arrs[f] for f in FEATURE_NAMES])

print("Loading signal...")
sig    = load_file(SIG_FILE)
X_sig  = to_matrix(sig)
ev_sig = sig["event_id"].astype(np.int64)

print("Loading W+jets background...")
bkg    = load_file(WJETS_FILE)
X_bkg  = to_matrix(bkg)
ev_bkg = bkg["event_id"].astype(np.int64)

print(f"Signal: {len(X_sig)}  Background: {len(X_bkg)}")

X      = np.vstack([X_sig, X_bkg])
y      = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
ev_all = np.concatenate([ev_sig, ev_bkg])

dira_idx = FEATURE_NAMES.index("dira")

# ========= 4. 有限值 + dira 过滤 =========
mask = np.isfinite(X).all(axis=1)
mask &= X[:, dira_idx] > -1
X, y, ev_all = X[mask], y[mask], ev_all[mask]
print(f"After finite filter: {X.shape}")

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 5. 按 event_id 划分 train/val/test =========
unique_ev = np.unique(ev_all)
rng2 = np.random.default_rng(SEED)
rng2.shuffle(unique_ev)
n      = len(unique_ev)
n_test = int(0.2 * n)
n_val  = int(0.1 * n)
ev_test  = set(unique_ev[:n_test])
ev_val   = set(unique_ev[n_test:n_test + n_val])
ev_train = set(unique_ev[n_test + n_val:])

train_m = np.array([e in ev_train for e in ev_all])
val_m   = np.array([e in ev_val   for e in ev_all])
test_m  = np.array([e in ev_test  for e in ev_all])

X_train, y_train = X[train_m], y[train_m]
X_val,   y_val   = X[val_m],   y[val_m]
X_test,  y_test  = X[test_m],  y[test_m]
print(f"Train: {len(X_train)} | Val: {len(X_val)} | Test: {len(X_test)}")

# ========= 6. 变量分布对比图 =========
print("Drawing variable comparison plots...")
X_sig_plot = X[y == 1]
X_bkg_plot = X[y == 0]
N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS

with PdfPages(f"{out_dir}/variable_comparison.pdf") as pdf:
    for page_start in range(0, len(FEATURE_NAMES), PER_PAGE):
        page_feats = FEATURE_NAMES[page_start:page_start + PER_PAGE]
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(25, 20))
        axes_flat = axes.flatten()
        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            fi = FEATURE_NAMES.index(feat)
            s_vals = X_sig_plot[:, fi]
            b_vals = X_bkg_plot[:, fi]
            all_vals = np.concatenate([s_vals, b_vals])
            lo, hi = np.percentile(all_vals, [1, 99])
            if lo == hi:
                lo, hi = all_vals.min(), all_vals.max()
            bins = np.linspace(lo, hi, 50)
            ax.hist(b_vals, bins=bins, density=True, alpha=0.6, color="steelblue", label="W+jets")
            ax.hist(s_vals, bins=bins, density=True, alpha=0.6, color="firebrick",  label="Signal")
            ax.set_xlabel(feat)
            ax.set_ylabel("a.u.")
            ax.legend(fontsize=13)
            if i == 0:
                hep.cms.label("Preliminary", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
        for j in range(len(page_feats), len(axes_flat)):
            axes_flat[j].set_visible(False)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
print("Saved variable_comparison.pdf")

# ========= 7. 类别权重 =========
n_sig_tr = int(np.sum(y_train == 1))
n_bkg_tr = int(np.sum(y_train == 0))
scale_pos_weight = n_bkg_tr / n_sig_tr
print(f"Train signal: {n_sig_tr} | Train background: {n_bkg_tr} | scale_pos_weight: {scale_pos_weight:.3f}")

# ========= 8. 训练 =========
print("\nStart training...")
model = xgb.XGBClassifier(
    n_estimators=2000,
    max_depth=5,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    min_child_weight=10,
    scale_pos_weight=scale_pos_weight,
    eval_metric="auc",
    early_stopping_rounds=30,
    random_state=SEED,
    n_jobs=4,
    verbosity=1,
)

model.fit(
    X_train, y_train,
    eval_set=[(X_val, y_val)],
    verbose=10,
)

best_iter = model.best_iteration
print(f"\nBest iteration: {best_iter}")

# ========= 9. 测试 =========
scores = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_test, scores)
print(f"Test AUC: {auc:.4f}")

np.save(f"{out_dir}/scores.npy",  scores)
np.save(f"{out_dir}/y_test.npy",  y_test)

val_results = model.evals_result()
val_aucs = val_results["validation_0"]["auc"]
np.save(f"{out_dir}/val_aucs.npy", np.array(val_aucs))

# ========= 10. 结果图 =========
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

axes[0].plot(val_aucs, color="steelblue")
axes[0].set_xlabel("Boosting round"); axes[0].set_ylabel("Val AUC")
axes[0].axvline(best_iter, color="firebrick", ls="--", label=f"best={best_iter}")
axes[0].legend()
hep.cms.label("Preliminary", ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

fpr, tpr, _ = roc_curve(y_test, scores)
axes[1].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[1].plot([0, 1], [0, 1], "k--")
axes[1].set_xlabel("False Positive Rate"); axes[1].set_ylabel("True Positive Rate"); axes[1].legend()
hep.cms.label("Preliminary", ax=axes[1], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

bins = np.linspace(0, 1, 51)
axes[2].hist(scores[y_test == 0], bins=bins, alpha=0.6, color="steelblue", label="W+jets", density=True)
axes[2].hist(scores[y_test == 1], bins=bins, alpha=0.6, color="firebrick", label="Signal", density=True)
axes[2].set_xlabel("BDT Score"); axes[2].set_ylabel("a.u."); axes[2].legend()
axes[2].set_yscale("log")
hep.cms.label("Preliminary", ax=axes[2], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

plt.tight_layout()
plt.savefig(f"{out_dir}/results.png", dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= 11. 特征重要性 =========
importances = model.feature_importances_
np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)
fig, ax = plt.subplots(figsize=(8, len(FEATURE_NAMES) * 0.35))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color="steelblue")
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([FEATURE_NAMES[i] for i in sorted_idx], fontsize=13)
ax.set_xlabel("XGBoost feature importance (gain)")
hep.cms.label("Preliminary", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 12. 保存模型 =========
joblib.dump(model, f"{out_dir}/bdt_model.pkl")
print("Done.")

"""
XGBoost BDT: signal vs W+jets / TT / DY, phipimu channel.

Reads the one-row-per-event Samples tuples directly (no Prepare macros needed);
each carries a Lumi_reweight branch used as the training weight.

Feature set (26) deliberately excludes every decay-kinematics variable: those
reconstruct phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the
resonance peaks used for the final fit. What is left (lifetime/vertex,
isolation/MET, Ds+mu kinematics) leaks < 0.04 into all three fit variables.
Dsmu_mass is a feature here, not a fit variable -- it has no resonant peak, and
keeping it out of the fit is what lets Ds_pt stay in the MVA.
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
SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
TREE = "Hphipimu/Events"

SIGNALS = ["HtoDsW", "HtoDsstW"]

# Four physically distinct enemies (used by v2/v3):
#   W+jets     combinatorial fakes from tracks inside jets
#   TT         REAL Ds from b->B->Ds (displaced, in a b jet) + high MET
#   DY         combinatorial fakes from PROMPT tracks + real Ds from jet charm
# The three training versions differ only in how the backgrounds are mixed:
#   v1  W+jets only                     -- baseline
#   v2  four backgrounds, PROPORTIONAL  -- subsample counts ~ physical yield (81:7.5:7.2:4.2)
#   v3  four backgrounds, FLAT          -- equal subsample count per sample
# Flat suppresses TT/DY 10-17x better than W+jets-only but lets W+jets through
# more; proportional restores W+jets emphasis at the cost of thin TT/DY stats.
FOUR = ["WtoLNu_2Jets", "TTtoLNu2Q", "TTto2L2Nu", "DYto2L_MLL50"]
CONFIG = {
    "v1": (["WtoLNu_2Jets"], "flat"),
    "v2": (FOUR, "proportional"),
    "v3": (FOUR, "flat"),
}

model_name = sys.argv[1] if len(sys.argv) > 1 else "v3"
BACKGROUNDS, SAMPLE_MODE = CONFIG[model_name]
out_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/{model_name}"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)

# ========= 2. 特征定义 =========
# Must stay identical to DNN_phipimu/DNN_phipimu.py.
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

SELECTION = ["num_os_mu_tight", "num_ss_mu_tight"]
LOAD_BRANCHES = FEATURE_NAMES + SELECTION + ["Lumi_reweight"]

# ========= 3. 加载数据 =========
def load(name, truth_filter):
    branches = LOAD_BRANCHES + (["Ds_match", "mu_match"] if truth_filter else [])
    with uproot.open(f"{SAMPLES_DIR}/{name}.root")[TREE] as t:
        a = t.arrays(branches, library="np")
    # exactly one tight muon in the event (the best muon is os to the pi)
    m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    X = np.column_stack([a[f][m] for f in FEATURE_NAMES])
    w = a["Lumi_reweight"][m].astype(np.float64)
    return X, w

print("Loading signal...")
sig = [load(s, True) for s in SIGNALS]
X_sig = np.vstack([x for x, _ in sig])
w_sig = np.concatenate([w for _, w in sig])

print(f"Loading {len(BACKGROUNDS)} backgrounds...")
bkg = [load(b, False) for b in BACKGROUNDS]

# Mix the backgrounds by SUBSAMPLING, not by down-weighting: every kept event
# stays at weight 1 (so N_eff = N), and the training mixture is set purely by how
# many events are drawn from each sample.
#   flat          equal count per sample (each background weighs the same)
#   proportional  count ~ physical weighted yield; capped by availability, so the
#                 sample that runs out first (W+jets: 81% of the yield but only
#                 38.7k events) fixes the total and the rest scale down from it.
counts = [len(x) for x, _ in bkg]
if SAMPLE_MODE == "flat":
    target = [min(counts)] * len(bkg)
elif SAMPLE_MODE == "proportional":
    frac = np.array([wt.sum() for _, wt in bkg]); frac /= frac.sum()
    K = min(counts[i] / frac[i] for i in range(len(bkg)))
    target = [min(counts[i], int(round(K * frac[i]))) for i in range(len(bkg))]
else:
    raise ValueError(f"unknown SAMPLE_MODE {SAMPLE_MODE}")

sub = []
for name, (x, wt), t, n_all in zip(BACKGROUNDS, bkg, target, counts):
    if n_all > t:
        keep = rng.choice(n_all, t, replace=False)
        x, wt = x[keep], wt[keep]
        print(f"  {name:16s} {t:>8,} events (subsampled from {n_all:,})")
    else:
        print(f"  {name:16s} {n_all:>8,} events (all)")
    sub.append((x, wt))
bkg = sub

X_bkg = np.vstack([x for x, _ in bkg])
w_bkg = np.concatenate([w for _, w in bkg])
# which background sample each row came from (-1 = signal)
src = np.concatenate([np.full(len(X_sig), -1)] +
                     [np.full(len(x), i) for i, (x, _) in enumerate(bkg)])

print(f"Signal: {len(X_sig):,} (weighted {w_sig.sum():.3f})")
print(f"Background: {len(X_bkg):,} ({SAMPLE_MODE})")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
w = np.concatenate([w_sig, w_bkg])

# ========= 4. 有限值 + dira 过滤 =========
dira_idx = FEATURE_NAMES.index("dira")
mask = np.isfinite(X).all(axis=1) & (X[:, dira_idx] > -1)
X, y, w, src = X[mask], y[mask], w[mask], src[mask]
print(f"After finite filter: {X.shape}")

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 5. 训练权重 =========
# Background: every row weighs the same here; the sample mix is set upstream by
# SAMPLE_MODE (flat or proportional), not by the weights. Lumi_reweight is
# deliberately NOT used: weighting by luminosity would hand the gradient straight
# back to W+jets (weight 288 vs 5 for TT, so it carries ~97% of sum(w^2)), which
# is the failure that left TT/DY unsuppressed in v1.
#
# Measured on the full 34-background weighted residual: proportional (v2) beats
# flat (v3) at every working point. Flat suppresses each background harder in
# isolation, but W+jets is 77% of the weighted yield, so what flat lets through
# there outweighs what it saves on TT/DY.
#
# Signal: by Lumi_reweight, which sets the DsW/DsstW relative contribution.
w_train = np.empty_like(w)
sig_m = y == 1
w_train[sig_m] = w[sig_m] / w[sig_m].sum()
w_train[~sig_m] = 1.0 / (~sig_m).sum()
# equal total weight for the two classes
for cls in (0, 1):
    m = y == cls
    w_train[m] *= len(y) / 2 / w_train[m].sum()

# ========= 6. 划分 train/val/test（每事例一行，直接按行随机划分）=========
idx = rng.permutation(len(y))
n_test, n_val = int(0.2 * len(y)), int(0.1 * len(y))
test_i, val_i, train_i = idx[:n_test], idx[n_test:n_test + n_val], idx[n_test + n_val:]

X_train, y_train, sw_train = X[train_i], y[train_i], w_train[train_i]
X_val,   y_val,   sw_val   = X[val_i],   y[val_i],   w_train[val_i]
X_test,  y_test,  sw_test  = X[test_i],  y[test_i],  w_train[test_i]
print(f"Train: {len(X_train):,} | Val: {len(X_val):,} | Test: {len(X_test):,}")

# ========= 7. 变量分布对比图 =========
print("Drawing variable comparison plots...")
N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS
X_sig_plot, w_sig_plot = X[y == 1], w[y == 1]
X_bkg_plot, w_bkg_plot = X[y == 0], w[y == 0]

with PdfPages(f"{out_dir}/variable_comparison.pdf") as pdf:
    for page_start in range(0, len(FEATURE_NAMES), PER_PAGE):
        page_feats = FEATURE_NAMES[page_start:page_start + PER_PAGE]
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(25, 20))
        axes_flat = axes.flatten()
        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            fi = FEATURE_NAMES.index(feat)
            s_vals, b_vals = X_sig_plot[:, fi], X_bkg_plot[:, fi]
            lo, hi = np.percentile(np.concatenate([s_vals, b_vals]), [1, 99])
            if lo == hi:
                lo, hi = min(s_vals.min(), b_vals.min()), max(s_vals.max(), b_vals.max())
            bins = np.linspace(lo, hi, 50)
            ax.hist(b_vals, bins=bins, weights=w_bkg_plot, density=True, alpha=0.6,
                    color="steelblue", label="Background (4)")
            ax.hist(s_vals, bins=bins, weights=w_sig_plot, density=True, alpha=0.6,
                    color="firebrick", label="Signal")
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

# ========= 8. 训练 =========
print("\nStart training...")
model = xgb.XGBClassifier(
    n_estimators=2000,
    max_depth=5,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    min_child_weight=10,
    eval_metric="auc",
    early_stopping_rounds=30,
    random_state=SEED,
    n_jobs=4,
    verbosity=1,
)

model.fit(
    X_train, y_train,
    sample_weight=sw_train,
    eval_set=[(X_val, y_val)],
    sample_weight_eval_set=[sw_val],
    verbose=10,
)

best_iter = model.best_iteration
print(f"\nBest iteration: {best_iter}")

# ========= 9. 测试 =========
scores = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_test, scores, sample_weight=sw_test)
print(f"Test AUC (weighted): {auc:.4f}")

np.save(f"{out_dir}/scores.npy", scores)
np.save(f"{out_dir}/y_test.npy", y_test)
np.save(f"{out_dir}/w_test.npy", w[test_i])          # physical Lumi_reweight, for yields

val_aucs = model.evals_result()["validation_0"]["auc"]
np.save(f"{out_dir}/val_aucs.npy", np.array(val_aucs))

# ========= 10. 结果图 =========
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

axes[0].plot(val_aucs, color="steelblue")
axes[0].set_xlabel("Boosting round"); axes[0].set_ylabel("Val AUC")
axes[0].axvline(best_iter, color="firebrick", ls="--", label=f"best={best_iter}")
axes[0].legend()
hep.cms.label("Simulation", ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

fpr, tpr, _ = roc_curve(y_test, scores, sample_weight=sw_test)
axes[1].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[1].plot([0, 1], [0, 1], "k--")
axes[1].set_xlabel("False Positive Rate"); axes[1].set_ylabel("True Positive Rate"); axes[1].legend()
hep.cms.label("Simulation", ax=axes[1], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

bins = np.linspace(0, 1, 51)
axes[2].hist(scores[y_test == 0], bins=bins, weights=sw_test[y_test == 0], alpha=0.6,
             color="steelblue", label="Background (4)", density=True)
axes[2].hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
             color="firebrick", label="Signal", density=True)
axes[2].set_xlabel("BDT Score"); axes[2].set_ylabel("a.u."); axes[2].legend()
axes[2].set_yscale("log")
hep.cms.label("Simulation", ax=axes[2], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

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
hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 12. 保存模型 =========
joblib.dump(model, f"{out_dir}/bdt_model.pkl")
print("Done.")

"""
XGBoost BDT: H->DsW signal vs W+jets, phipimu channel -- TEST BUILD.

Single-signal, single-background test run, put together before DsstW's tuplizer
output and the other background processes (TT, DY) are ready in the new Samples/
format (see project_tuplizer_dedup_and_qc / project_wjets_background_incomplete
memories). Once DsstW + the other backgrounds are back, this needs revisiting --
in particular the "single background source" weighting below stops applying once
there's more than one physically distinct background to mix.

Reads the one-row-per-event Samples tuples directly (flat `Events` tree -- these
files, produced by highest_pt_dedup.py, carry no embedded Lumi_reweight branch,
unlike the older 34-background Samples/ files; see project_tuplizer_dedup_and_qc).

Feature set excludes every decay-kinematics variable: those reconstruct
phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the resonance peaks
used for the final fit (measured on the original 37-feature set; the 18 jet/b-tag
features appended below look mass-blind by construction -- dR's, b-scores, HT,
jet counts -- but have NOT been re-screened with the same R^2 test yet. Do that
before this stops being a test build.)
Dsmu_mass is a feature here, not a fit variable -- it has no resonant peak, and
keeping it out of the fit is what lets Ds_pt stay in the MVA.
"""
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
TREE = "Events"

SIGNAL = "HtoDsW"   # DsstW tuplizer output not ready yet -- single signal for this test

# WtoLNu2Jets is split into three jet-multiplicity bins of the SAME physical
# process (not physically distinct backgrounds like the old TT/DY mix), so they
# are weight-merged here using the full statistics of all three files -- not
# subsampled/mixed the way genuinely different processes would be. Weights are
# the equivalent-luminosity reweighting factors (external -- these files have no
# embedded Lumi_reweight branch).
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}

out_dir = "/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/dsw_wjets_test"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)

# ========= 2. 特征定义 =========
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
    # --- jet/b-tag environment, added 2026-08-03. Five of the eighteen new
    # branches earn a slot; the rest are either redundant with these or scored
    # ~0 importance against W+jets (they were designed to fight TT, which this
    # test build does not contain -- revisit when TT is in the mix).
    # The Ds's own local environment: real Ds from H is isolated, a fake built
    # from tracks inside a jet is not. Highest-ranked of the new variables.
    "Ds_nearestjet_bscore",
    "Ds_nearestjet_dR",
    "jet_leadpt",
    # Soft b jets that the tight veto below does not catch: tight is a subset of
    # medium is a subset of loose, so requiring 0 tight b jets leaves these two
    # free to vary -- they are not degenerate with the cut.
    "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR04",
]

# Exactly one tight muon in the event (the best muon itself, no extra of either
# charge). All samples loaded below are OS files, so this single expression
# applies unmirrored to signal and all three WtoLNu2Jets bins -- see
# H_phipimu.cc and project_data2022_and_sr_definition memory.
#
# Plus a b-jet veto. This one is a CUT, not a feature: its job is to remove the
# process the background model does not contain at all (tt, with real b jets),
# not to discriminate signal from W+jets. Measured on data sidebands (the only
# place the missing tt is visible), the tight WP keeps 99.3% of signal while
# removing 25% of data; medium removes more background (34%) but costs 3.1% of
# signal -- tight was chosen to protect the scarce signal statistics, leaving
# the residual tt to be modelled properly later. dR04 matches the AK4 jet cone
# radius, so the exclusion cone lines up with "the jet the Ds itself would be
# clustered into".
B_VETO = "n_b_jets_tight_dR04"
SELECTION = ["num_os_mu_tight", "num_ss_mu_tight", B_VETO]
LOAD_BRANCHES = list(dict.fromkeys(FEATURE_NAMES + SELECTION))

# ========= 3. 加载数据 =========
def load(name, truth_filter):
    branches = LOAD_BRANCHES + (["Ds_match", "mu_match"] if truth_filter else [])
    with uproot.open(f"{SAMPLES_DIR}/{name}.root")[TREE] as t:
        a = t.arrays(branches, library="np")
    n_in = len(a[B_VETO])
    m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= a[B_VETO] == 0
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    X = np.column_stack([a[f][m] for f in FEATURE_NAMES])
    print(f"  {name:18s} {n_in:>8,} -> {m.sum():>8,} rows after selection")
    return X

print("Loading signal (HtoDsW)...")
X_sig = load(SIGNAL, truth_filter=True)
# Single signal process in this test build -- no embedded Lumi_reweight branch
# either way (see docstring), and with only one process its physical value
# wouldn't affect training (any positive constant collapses to the same
# normalized w_train below). Uniform weight, NOT a physical yield.
w_sig = np.ones(len(X_sig))
print(f"Signal: {len(X_sig):,} events")

print(f"Loading WtoLNu2Jets ({'+'.join(WJETS_WEIGHTS)}, lumi-reweighted)...")
bkg_X, bkg_w = [], []
for name, wgt in WJETS_WEIGHTS.items():
    x = load(name, truth_filter=False)
    bkg_X.append(x)
    bkg_w.append(np.full(len(x), wgt))
X_bkg = np.vstack(bkg_X)
w_bkg = np.concatenate(bkg_w)

print(f"Signal: {len(X_sig):,} (uniform weight -- single process, see above)")
print(f"Background: {len(X_bkg):,} (WtoLNu2Jets, full stats, lumi-reweighted merge of 3 jet bins)")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
w = np.concatenate([w_sig, w_bkg])

# ========= 4. 有限值 + dira 过滤 =========
dira_idx = FEATURE_NAMES.index("dira")
mask = np.isfinite(X).all(axis=1) & (X[:, dira_idx] > -1)
X, y, w = X[mask], y[mask], w[mask]
print(f"After finite filter: {X.shape}")

# Verify the lumi-reweighting actually took effect (not just carried through as
# dead weight): weighted mean should match the physical mixture (~1.095), not the
# raw-event-count mixture (~1.553, dominated by 2J simply having more events) --
# see project_wjets_background_incomplete memory for where 1.553 came from.
jet30_idx = FEATURE_NAMES.index("num_jets_pt30")
check = np.average(X[y == 0][:, jet30_idx], weights=w[y == 0])
print(f"[check] weighted <num_jets_pt30> for background = {check:.3f} (expect ~1.095, not ~1.553)")

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 5. 训练权重 =========
# Signal: uniform (single process, see above).
# Background: weighted BY the physical Lumi_reweight (not subsampled) so the
# 0J/1J/2J jet-multiplicity mixture matches the true cross-section ratio
# (2.33 : 0.54 : 0.27) instead of these files' raw event-count ratio, which is
# dominated by 2J having more triggered/skimmed events than 0J.
w_train = np.empty_like(w)
sig_m = y == 1
w_train[sig_m] = w[sig_m] / w[sig_m].sum()
w_train[~sig_m] = w[~sig_m] / w[~sig_m].sum()
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
                    color="steelblue", label="WtoLNu2Jets")
            ax.hist(s_vals, bins=bins, weights=w_sig_plot, density=True, alpha=0.6,
                    color="firebrick", label="HtoDsW")
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
# NOT a physical yield for signal (uniform 1.0 placeholder, see step 3) -- only
# meaningful for background-internal ratios (e.g. jet-bin composition) until a
# real signal Lumi_reweight is available again.
np.save(f"{out_dir}/w_test.npy", w[test_i])

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
             color="steelblue", label="WtoLNu2Jets", density=True)
axes[2].hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
             color="firebrick", label="HtoDsW", density=True)
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
fig, ax = plt.subplots(figsize=(8, len(FEATURE_NAMES) * 0.3))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color="steelblue")
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([FEATURE_NAMES[i] for i in sorted_idx], fontsize=11)
ax.set_xlabel("XGBoost feature importance (gain)")
hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 12. 保存模型 =========
joblib.dump(model, f"{out_dir}/bdt_model.pkl")
print("Done.")

"""
XGBoost BDT: H->Ds(*)W signal (OS) vs one of two background choices, phipimu channel.

    python3 BDT.py ssdata           # background = SS data (Data_SS_merged, even event)
    python3 BDT.py wjets            # background = WtoLNu2Jets MC (0J+1J+2J, lumi-reweighted)
    python3 BDT.py ssdata 101       # same config, different seed -> ..._f44_seed101

**Everything except the algorithm is now shared with DNN.py via ../DNN/common.py**:
the 41 features and their order, the preprocessing chain, the selection, the signal
mixture, the background table, the class weighting, the 70/10/20 split and the seed.
That is the point of this rewrite (2026-08-14). Before it, BDT differed from the DNN
in three things at once -- 42 features (its own fork: BDT-only {Kp,Km,pi}_PVassoc,
DNN-only {Kp,Km,pi}_ipxysig + jet_ht), W+jets MC as background, and HtoDsW as the
only signal -- so "BDT vs DNN" measured nothing in particular.

The preprocessing is common.preprocess() unchanged, log transform included. Trees are
invariant under a monotone transform so the log costs XGBoost nothing, and feeding
both algorithms bit-identical inputs is what makes the head-to-head clean; a
`log=False` branch would only be one more thing that can drift. One cost is logged
rather than fixed: SENTINEL_FIX collapses Ds_nearestjet_bscore's -777 (no jet) and
~-10 (tagger unavailable) onto a single -0.1, and a tree could have split those apart.
Changing it would alter two things at once, so it stays a follow-up.

Hyperparameters and eval_metric are deliberately UNCHANGED from the 42-feature build,
so this run isolates the feature/background/signal alignment. eval_metric stays "auc"
rather than following the DNN to a tail metric: it is natively weight-aware, and the
choice barely matters here -- on the previous build 490 of the 1085 boosting rounds
sat within 1e-4 of the best val AUC, i.e. where boosting stops is nearly free. Test
rejection at the actual working points is printed below for comparability instead.

Feature set excludes every decay-kinematics variable: those reconstruct
phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the resonance peaks used
for the final fit. Dsmu_mass is a feature here, not a fit variable -- it has no
resonant peak, and keeping it out of the fit is what lets Ds_pt stay in the MVA. See
common.py for the per-feature rationale, including the jet/b-tag block that was
re-screened against SS on 2026-08-10.

Reads NanoMVA/Samples/ (custom NanoAOD + POG corrections, since the 2026-09-17
migration -- see ../DNN/common.py).
"""
import os
import sys

import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import mplhep as hep
plt.style.use(hep.style.CMS)

from sklearn.metrics import roc_curve, roc_auc_score
import joblib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNN"))
import common
from common import SIGNALS

# ========= 1. 配置 =========
# Feature set, preprocessing, selection and the background table all live in
# common.py -- see its docstring for why nothing here may re-declare them.
if len(sys.argv) not in (2, 3):
    sys.exit(f"usage: python3 {os.path.basename(__file__)} "
             f"<{'|'.join(common.TRAININGS)}> [seed]")
BKG = sys.argv[1]
cfg = common.resolve(BKG)          # validates before anything is created on disk
feats = cfg["features"]
SEED = int(sys.argv[2]) if len(sys.argv) == 3 else common.SEED_DEFAULT

out_dir = common.out_dir_for(BKG, tool="BDT", seed=SEED)
os.makedirs(out_dir, exist_ok=True)
print(f"Background: {BKG} ({cfg['label']}), seed={SEED} -> {out_dir}")
print(f"Feature set: {cfg['feature_set']} ({len(feats.names)} features)")

rng = np.random.default_rng(SEED)

# A stage-2 config pre-filters its training rows on a stage-1 score cut. Resolve the
# threshold once and stamp it -- it is the definition of the training population.
STAGE1_CUT = None
if cfg.get("stage1_cut"):
    _thr, _nsig = common.stage1_threshold(cfg["stage1_cut"])
    STAGE1_CUT = dict(cfg["stage1_cut"], thr=_thr)
    with open(f"{out_dir}/stage1_filter.txt", "w") as _f:
        for _k, _v in {**STAGE1_CUT, "n_sig_rows": _nsig}.items():
            _f.write(f"{_k} {_v}\n")
    print(f"Stage-1 pre-filter: {cfg['stage1_cut']['key']} score_dsmu > {_thr:.6f}")

# Every panel's CMS stamp follows the background: real data gets "Preliminary" and
# data=True, simulation gets "Simulation" and data=False. This file used to hardcode
# "Simulation" in three places, which would have labelled every SS-data panel wrong.
CMS_TEXT = "Preliminary" if cfg["is_data"] else "Simulation"
BKG_LABEL = cfg["label"]


def cms_label(ax, fontsize=13):
    hep.cms.label(CMS_TEXT, ax=ax, rlabel="2022 (13.6 TeV)",
                  data=cfg["is_data"], fontsize=fontsize)


# ========= 2. 加载数据 =========
print(f"Loading signal ({'+'.join(SIGNALS)})...")
# weight is an ABSOLUTE per-row column (sigma*L/sum(w) per (channel, era) piece) as
# of the NanoMVA migration -- it already encodes the correct DsW:DsstW and
# preEE:postEE mixture. SIGNALS is a one-element list (Signal_merged, both
# channels/eras inside); the loop is kept so this reads identically to the
# background loop below.
sig_X, sig_w = [], []
for name in SIGNALS:
    x, ex = common.load(name, feats=feats, truth_filter=True, stage1=STAGE1_CUT,
                        extra_branches=("weight",))
    sig_X.append(x)
    sig_w.append(ex["weight"])
X_sig = np.vstack(sig_X)
w_sig = np.concatenate(sig_w)

parity_note = "" if cfg["parity"] is None else (
    f", event parity=={cfg['parity']} only -- the other parity is held out "
    f"for a later non-circular sculpting check")
print(f"Loading background ({cfg['label']}: {'+'.join(cfg['samples'])}{parity_note})...")
# weighted=True reads the per-event weight from the file (W+jets: sigma*L/sum(w)
# carrying gen_weight's sign); weighted=False uses the flat per-sample value from the
# table (data, which carries no luminosity weight). The W+jets branch used to multiply
# by three stale per-file constants instead -- wrong, because the amcatnloFXFX negative
# fraction varies 3.4x across the jet bins, so no per-file constant survives a selection.
bkg_X, bkg_w = [], []
for name, wgt in cfg["samples"].items():
    load_kw = dict(feats=feats, ss=cfg["ss"], parity=cfg["parity"],
                   gen_mode_filter=cfg.get("gen_mode_filter"), stage1=STAGE1_CUT)
    if cfg["weighted"]:
        x, ex = common.load(name, extra_branches=("weight",), **load_kw)
        bkg_w.append(ex["weight"])
    else:
        x, _ = common.load(name, **load_kw)
        bkg_w.append(np.full(len(x), wgt))
    bkg_X.append(x)
X_bkg = np.vstack(bkg_X)
w_bkg = np.concatenate(bkg_w)

print(f"Signal: {len(X_sig):,} (both channels at their physical proportion)")
print(f"Background: {len(X_bkg):,} ({cfg['label']})")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
w = np.concatenate([w_sig, w_bkg])

# ========= 3. 预处理：sentinel 规整 + 有限值过滤 + 特征变换 =========
# The same common.preprocess DNN.py trains with and predict.py scores with.
X, mask = common.preprocess(X, feats)
X, y, w = X[mask], y[mask], w[mask]
print(f"After finite filter: {X.shape}")

# Verify the lumi-reweighting actually took effect: the weighted mean should match the
# physical jet-bin mixture, not the raw event-count one. Only meaningful when the
# samples carry unequal weights -- the SS-data config is a single sample at weight 1.
if cfg["weighted"]:
    jet30_idx = feats.names.index("num_jets_pt30_corr")
    check = np.average(X[y == 0][:, jet30_idx], weights=w[y == 0])
    print(f"[check] weighted <num_jets_pt30_corr> for background = {check:.3f} "
          f"(compare against the raw-concat mean printed above)")

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 4. 训练权重 =========
# The mixture WITHIN each class is physical (DsW:DsstW by production size, jet bins
# by cross section), then the two classes are scaled to equal total weight, so
# signal:background in the loss is 1:1 by construction.
#
# |w|, not w: an amcatnloFXFX background (WtoLNu2Jets_GenTag here) carries negative
# gen weights -- 27% of the post-stage-1 modes-1-5 rows. XGBoost's hist DMatrix
# rejects negative sample_weight outright ("Weights must be positive values"), and a
# per-event P(class|x) fit has no use for the NLO subtraction term anyway. Training
# on |w| keeps every row (N_eff 3.6k -> 13k) at the cost of a slightly conservative
# background shape. The SIGNED weight is what gets saved for downstream yields/FOM
# (w_test.npy below) -- only the fit sees |w|. Signal weights are positive (relative
# luminosity), so abs is a no-op there. DNN.py keeps signed w: torch BCE does not
# validate, so its wjets model silently trained on ~15% negatives.
aw = np.abs(w)
w_train = np.empty_like(aw)
sig_m = y == 1
w_train[sig_m] = aw[sig_m] / aw[sig_m].sum()
w_train[~sig_m] = aw[~sig_m] / aw[~sig_m].sum()
for cls in (0, 1):
    m = y == cls
    w_train[m] *= len(y) / 2 / w_train[m].sum()

# ========= 5. 划分 train/val/test（每事例一行，直接按行随机划分）=========
# Same fractions and same rng stream as DNN.py.
idx = rng.permutation(len(y))
n_test, n_val = int(0.2 * len(y)), int(0.1 * len(y))
test_i, val_i, train_i = idx[:n_test], idx[n_test:n_test + n_val], idx[n_test + n_val:]

X_train, y_train, sw_train = X[train_i], y[train_i], w_train[train_i]
X_val,   y_val,   sw_val   = X[val_i],   y[val_i],   w_train[val_i]
X_test,  y_test,  sw_test  = X[test_i],  y[test_i],  w_train[test_i]
print(f"Train: {len(X_train):,} | Val: {len(X_val):,} | Test: {len(X_test):,}")

# ========= 6. 变量分布对比图 =========
print("Drawing variable comparison plots...")
N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS
X_sig_plot, w_sig_plot = X[y == 1], w[y == 1]
X_bkg_plot, w_bkg_plot = X[y == 0], w[y == 0]

with PdfPages(common.fig_path("training", "variable_comparison.pdf",
                              key=BKG, tool="BDT", seed=SEED)) as pdf:
    for page_start in range(0, len(feats.names), PER_PAGE):
        page_feats = feats.names[page_start:page_start + PER_PAGE]
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(25, 20))
        axes_flat = axes.flatten()
        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            fi = feats.names.index(feat)
            s_vals, b_vals = X_sig_plot[:, fi], X_bkg_plot[:, fi]
            lo, hi = np.percentile(np.concatenate([s_vals, b_vals]), [1, 99])
            if lo == hi:
                lo, hi = min(s_vals.min(), b_vals.min()), max(s_vals.max(), b_vals.max())
            bins = np.linspace(lo, hi, 50)
            ax.hist(b_vals, bins=bins, weights=w_bkg_plot, density=True, alpha=0.6,
                    color="steelblue", label=BKG_LABEL)
            ax.hist(s_vals, bins=bins, weights=w_sig_plot, density=True, alpha=0.6,
                    color="firebrick", label="Signal (DsW+DsstW)")
            ax.set_xlabel(feat)
            ax.set_ylabel("a.u.")
            ax.legend(fontsize=13)
            if i == 0:
                cms_label(ax)
        for j in range(len(page_feats), len(axes_flat)):
            axes_flat[j].set_visible(False)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
print("Saved variable_comparison.pdf")

# ========= 7. 训练 =========
print("\nStart training...")
hp = common.BDT_PARAMS[cfg["feature_set"]]
print(f"XGBoost hyperparameters ({cfg['feature_set']}): {hp}")
model = xgb.XGBClassifier(
    **hp,
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

# ========= 8. 测试 =========
scores = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_test, scores, sample_weight=sw_test)
print(f"Test AUC (weighted): {auc:.4f}")

# AUC is saturated (~0.996-0.998) and is a bulk measure; the analysis works at
# eff_sig 50-90%. Print rejection there too so BDT and DNN can be compared on the
# quantity that actually matters. The real head-to-head is still OS sideband data
# (compare_models.py) -- this is the in-sample cross-check, measured against the
# training background rather than against unseen data.
print("\nRejection at matched signal efficiency (test set, vs the training background):")
ss_, sb_ = np.sort(scores[y_test == 1]), scores[y_test == 0]
print(f"{'eff_sig':>8} {'thr':>10} {'eps_bkg':>10} {'N_bkg':>8} {'rejection':>11}")
for e in (0.90, 0.70, 0.50):
    thr = np.quantile(ss_, 1 - e)
    n = int((sb_ > thr).sum())
    eb = n / len(sb_)
    print(f"{100*e:7.0f}% {thr:10.5f} {eb:10.5f} {n:8,} "
          f"{1/eb if eb > 0 else np.inf:11.1f}")

np.save(f"{out_dir}/scores.npy", scores)
np.save(f"{out_dir}/y_test.npy", y_test)
# Relative weights only, not absolute yields: within signal the true DsW:DsstW mixture
# and within background the true jet-bin mixture, but the two classes are on
# different, unrelated scales -- do not compute a cross-class significance from this.
np.save(f"{out_dir}/w_test.npy", w[test_i])

val_aucs = model.evals_result()["validation_0"]["auc"]
np.save(f"{out_dir}/val_aucs.npy", np.array(val_aucs))

# ========= 9. 结果图 =========
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

axes[0].plot(val_aucs, color="steelblue")
axes[0].set_xlabel("Boosting round"); axes[0].set_ylabel("Val AUC")
axes[0].axvline(best_iter, color="firebrick", ls="--", label=f"best={best_iter}")
axes[0].legend()
cms_label(axes[0])

fpr, tpr, _ = roc_curve(y_test, scores, sample_weight=sw_test)
axes[1].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[1].plot([0, 1], [0, 1], "k--")
axes[1].set_xlabel("False Positive Rate"); axes[1].set_ylabel("True Positive Rate"); axes[1].legend()
cms_label(axes[1])

bins = np.linspace(0, 1, 51)
axes[2].hist(scores[y_test == 0], bins=bins, weights=sw_test[y_test == 0], alpha=0.6,
             color="steelblue", label=BKG_LABEL, density=True)
axes[2].hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
             color="firebrick", label="Signal (DsW+DsstW)", density=True)
axes[2].set_xlabel("BDT Score"); axes[2].set_ylabel("a.u."); axes[2].legend()
axes[2].set_yscale("log")
cms_label(axes[2])

plt.tight_layout()
plt.savefig(common.fig_path("training", "results.png", key=BKG, tool="BDT", seed=SEED),
            dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= 10. 特征重要性 =========
# XGBoost's native gain importance. NOT the same quantity as DNN.py's permutation
# importance (AUC drop), so the two rankings are indicative, not directly comparable.
importances = model.feature_importances_
np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)
fig, ax = plt.subplots(figsize=(8, len(feats.names) * 0.3))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color="steelblue")
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([feats.names[i] for i in sorted_idx], fontsize=11)
ax.set_xlabel("XGBoost feature importance (gain)")
cms_label(ax)
plt.tight_layout()
plt.savefig(common.fig_path("training", "feature_importance.png",
                            key=BKG, tool="BDT", seed=SEED), dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 11. 保存模型 =========
joblib.dump(model, f"{out_dir}/bdt_model.pkl")
# Stamp the exact feature list (and order) this model was trained on. predict.py
# refuses to score against a dir whose stamp disagrees with common.py -- without it,
# a feature added or reordered later would silently produce wrong scores. This guard
# is what makes wiring BDT into common.py safe.
common.save_feature_names(out_dir, feats)
print(f"Done. ({len(feats.names)} features, background={BKG}, seed={SEED})")

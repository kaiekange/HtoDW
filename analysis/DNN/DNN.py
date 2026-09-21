"""
DNN: H->Ds(*)W signal (OS) vs one of two background choices, phipimu channel.

    python3 DNN.py ssdata           # background = SS data (Data_SS_merged, even event)
    python3 DNN.py wjets            # background = WtoLNu2Jets MC (0J+1J+2J, lumi-reweighted)
    python3 DNN.py ssdata 101       # same config, different seed -> ..._f44_seed101

The optional seed exists to measure the run-to-run noise floor. Any claim that a
model change improved background rejection has to clear that floor, and until it
was measured (2026-08-14) there was no number for it. The seed drives BOTH the
numpy rng (which fixes the train/val/test split) and torch's init, so the spread it
produces is the total run-to-run variation -- deliberately including split variance.
Do not try to decompose it afterwards. Seed 42 (common.SEED_DEFAULT) writes to the
canonical dir; anything else gets a _seed{N} suffix so the canonical models stay put.

Everything except the background -- feature set, seed, split fractions, network,
class balancing, early stopping -- is identical between the two, so the two trained
models differ in exactly one thing and can be compared head to head. That was NOT
true before 2026-08-10: the W+jets model was 40 features and the SS model 43, which
silently invalidated the first comparison. The feature list now lives in common.py
so it cannot drift again.

**W+jets background** was the original choice. It is known to be missing tt/DY
entirely (project_wjets_background_incomplete) -- the data/MC gap widened with the
score cut, meaning the real background the classifier never saw was exactly the
hardest part.

**SS data background** contains the genuine process mixture (tt, DY, QCD, ...) with
no need to add each one to Samples/ separately, at the cost of the SS-purity caveat
in project_mva_sculpting_finding (SS is not purely combinatorial after tight cuts --
some of it is real Ds with a wrong-charge pi). Signal stays OS-selected, so with
this background the two training classes differ structurally in charge
configuration, not just physics. Checked before adopting it (see
check_os_vs_ss_features.py, OS sideband vs SS post-veto, no signal involved): the
largest single-feature Cohen's d across all 43 features is 0.10 (Ds_pt), KS <=0.041
-- the charge-config artifact is small at the univariate level. That does not rule
out the classifier finding a multivariate combination encoding OS-vs-SS instead of
signal-vs-background; the real check is scoring OS sideband data with the trained
model and confirming it does not come out looking signal-like.

Not kept in lockstep with BDT/BDT.py: as of 2026-08-04 the two use different feature
sets and are not directly comparable. BDT deliberately does not import common.py.

Reads NanoMVA/Samples/ (custom NanoAOD + POG corrections, since the 2026-09-17
migration -- see common.py). Signal is one file (Signal_merged, both channels and
eras inside), and every row carries a per-event absolute `weight` column already.
"""
import sys
import os

import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import mplhep as hep
plt.style.use(hep.style.CMS)

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

import common
from common import SIGNALS

# ========= 1. 配置 =========
# Feature set, preprocessing, network and the background table all live in
# common.py -- see its docstring for why nothing here may re-declare them.
if len(sys.argv) not in (2, 3):
    sys.exit(f"usage: python3 {os.path.basename(__file__)} "
             f"<{'|'.join(common.TRAININGS)}> [seed]")
BKG = sys.argv[1]
cfg = common.resolve(BKG)          # validates before anything is created on disk
feats = cfg["features"]
SEED = int(sys.argv[2]) if len(sys.argv) == 3 else common.SEED_DEFAULT

out_dir = common.out_dir_for(BKG, tool="DNN", seed=SEED)
os.makedirs(out_dir, exist_ok=True)
print(f"Background: {BKG} ({cfg['label']}), seed={SEED} -> {out_dir}")
print(f"Feature set: {cfg['feature_set']} ({len(feats.names)} features)")

rng = np.random.default_rng(SEED)
torch.manual_seed(SEED)

# A stage-2 config pre-filters its training rows on a stage-1 score cut. Resolve the
# threshold once here and stamp it -- it is the definition of the training
# population, not just provenance, and it moves if stage-1 is ever re-scored.
STAGE1_CUT = None
if cfg.get("stage1_cut"):
    _thr, _nsig = common.stage1_threshold(cfg["stage1_cut"])
    STAGE1_CUT = dict(cfg["stage1_cut"], thr=_thr)
    with open(f"{out_dir}/stage1_filter.txt", "w") as _f:
        for _k, _v in {**STAGE1_CUT, "n_sig_rows": _nsig}.items():
            _f.write(f"{_k} {_v}\n")
    print(f"Stage-1 pre-filter: {cfg['stage1_cut']['key']} score_dsmu > {_thr:.6f}")

# Every panel's CMS stamp follows the background: real data gets "Preliminary" and
# data=True, simulation gets "Simulation" and data=False. Hardcoding it is what made
# the pre-refactor SS run label every plot "Simulation" while training on real data.
CMS_TEXT = "Preliminary" if cfg["is_data"] else "Simulation"
BKG_LABEL = cfg["label"]


def cms_label(ax, fontsize=13):
    hep.cms.label(CMS_TEXT, ax=ax, rlabel="2022 (13.6 TeV)",
                  data=cfg["is_data"], fontsize=fontsize)


print(f"Loading signal ({'+'.join(SIGNALS)})...")
# weight is an ABSOLUTE per-row column (sigma*L/sum(w) per (channel, era) piece) as
# of the NanoMVA migration, not the old relative SIGNAL_ERA_WEIGHT -- it already
# encodes the correct DsW:DsstW and preEE:postEE mixture. SIGNALS is a one-element
# list (Signal_merged, both channels/eras inside), the loop is kept so this reads
# identically to the background loop below.
sig_X, sig_w = [], []
for name in SIGNALS:
    x, ex = common.load(name, feats=feats, truth_filter=True, stage1=STAGE1_CUT,
                        extra_branches=("weight",))
    sig_X.append(x)
    sig_w.append(ex["weight"])
X_sig = np.vstack(sig_X)
w_sig = np.concatenate(sig_w)

# Background: one loop over the configured samples. The SS config is a single
# sample at weight 1; the W+jets config is the three jet bins at their
# equivalent-luminosity weights. Same code path, so nothing but the table differs.
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

print(f"Signal: {len(X_sig):,} (both channels at their physical proportion, see above)")
print(f"Background: {len(X_bkg):,} ({cfg['label']})")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
w = np.concatenate([w_sig, w_bkg])

# ========= 4. 预处理：sentinel 规整 + 有限值过滤 + 特征变换 =========
# One shared implementation (common.preprocess) so predict.py cannot drift from it.
X, mask = common.preprocess(X, feats)
X, y, w = X[mask], y[mask], w[mask]
print(f"After finite filter: {X.shape}")

# Pileup control models blank a subset of the inputs to a constant (see
# common.neutralize). The feature COUNT is unchanged, so the scaler, the input layer
# and the feature_names.txt stamp stay valid and predict.py needs no special case
# beyond applying the identical blanking -- which it does, from the same config.
_blanked = common.neutralize(X, cfg.get("neutralize", ()), feats)
if _blanked:
    print(f"Neutralised {len(_blanked)} feature(s) to a constant: {', '.join(_blanked)}")

# Verify the lumi-reweighting actually took effect: the weighted mean should match
# the physical jet-bin mixture, not the raw event-count one (which is dominated by
# 2J simply having more events). Only meaningful when the samples carry unequal
# weights -- the SS-data config is a single sample at weight 1.
if cfg["weighted"]:
    jet30_idx = feats.names.index("num_jets_pt30_corr")
    check = np.average(X[y == 0][:, jet30_idx], weights=w[y == 0])
    print(f"[check] weighted <num_jets_pt30_corr> for background = {check:.3f} "
          f"(compare against the raw-concat mean printed above)")

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 6. 训练权重 =========
# Signal: DsW and DsstW at equal per-event weight, which is their physical mixture
# because the two productions were sized in proportion to cross section (see SIGNALS).
# Background: weighted BY the physical Lumi_reweight (not subsampled) so the
# 0J/1J/2J jet-multiplicity mixture matches the true cross-section ratio
# (2.33 : 0.54 : 0.27) instead of these files' raw event-count ratio.
# The two classes are then scaled to equal total weight below, so the signal:background
# ratio in the loss is 1:1 by construction -- only the mixture *within* each class is
# physical.
w_train = np.empty_like(w)
sig_m = y == 1
w_train[sig_m] = w[sig_m] / w[sig_m].sum()
w_train[~sig_m] = w[~sig_m] / w[~sig_m].sum()
# equal total weight for the two classes
for cls in (0, 1):
    m = y == cls
    w_train[m] *= len(y) / 2 / w_train[m].sum()

# ========= 7. 划分 train/val/test（每事例一行，直接按行随机划分）=========
idx = rng.permutation(len(y))
n_test, n_val = int(0.2 * len(y)), int(0.1 * len(y))
test_i, val_i, train_i = idx[:n_test], idx[n_test:n_test + n_val], idx[n_test + n_val:]

X_train, y_train, sw_train = X[train_i], y[train_i], w_train[train_i]
X_val,   y_val,   sw_val   = X[val_i],   y[val_i],   w_train[val_i]
X_test,  y_test,  sw_test  = X[test_i],  y[test_i],  w_train[test_i]
print(f"Train: {len(X_train):,} | Val: {len(X_val):,} | Test: {len(X_test):,}")

# ========= 8. 变量分布对比图 =========
print("Drawing variable comparison plots...")
N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS
X_sig_plot, w_sig_plot = X[y == 1], w[y == 1]
X_bkg_plot, w_bkg_plot = X[y == 0], w[y == 0]

with PdfPages(common.fig_path("training", "variable_comparison.pdf",
                              key=BKG, tool="DNN", seed=SEED)) as pdf:
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

# ========= 9. 标准化 =========
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val   = scaler.transform(X_val)
X_test  = scaler.transform(X_test)
np.save(f"{out_dir}/scaler_mean.npy",  scaler.mean_)
np.save(f"{out_dir}/scaler_scale.npy", scaler.scale_)

# ========= 10. 转 tensor =========
def to_tensor(X, y=None):
    Xt = torch.tensor(X, dtype=torch.float32)
    if y is not None:
        return Xt, torch.tensor(y.reshape(-1, 1), dtype=torch.float32)
    return Xt

X_train_t, y_train_t = to_tensor(X_train, y_train)
X_val_t,   y_val_t   = to_tensor(X_val,   y_val)
X_test_t             = to_tensor(X_test)

# ========= 11. 模型 =========
# Defined in common.py so predict.py rebuilds the identical architecture.
model = common.phipimuNet(X_train.shape[1])
print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

# ========= 12. 训练 =========
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
criterion = nn.BCEWithLogitsLoss()
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

EPOCHS, PATIENCE, BATCH_SIZE = 200, 30, 1024

# Sample by the physical luminosity weight (already class-balanced above), so the
# network sees the true background composition rather than a flat MC mixture.
sampler = torch.utils.data.WeightedRandomSampler(
    torch.tensor(sw_train, dtype=torch.float64), len(y_train), replacement=True
)
dataset = torch.utils.data.TensorDataset(X_train_t, y_train_t)
loader  = torch.utils.data.DataLoader(dataset, batch_size=BATCH_SIZE, sampler=sampler)

best_val_auc, patience_counter, best_state = -1.0, 0, None
train_losses, val_losses, val_aucs = [], [], []

print("\nStart training...")
for epoch in range(EPOCHS):
    model.train()
    epoch_loss = 0.0
    for xb, yb in loader:
        pred = model(xb)
        loss = criterion(pred, yb)
        optimizer.zero_grad(); loss.backward(); optimizer.step()
        epoch_loss += loss.item() * len(xb)
    epoch_loss /= len(X_train_t)

    model.eval()
    with torch.no_grad():
        val_pred   = model(X_val_t)
        val_loss   = criterion(val_pred, y_val_t).item()
        val_scores = torch.sigmoid(val_pred).numpy().flatten()
    val_auc = roc_auc_score(y_val, val_scores, sample_weight=sw_val)

    train_losses.append(epoch_loss)
    val_losses.append(val_loss)
    val_aucs.append(val_auc)
    scheduler.step(val_loss)

    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1:03d} | Train {epoch_loss:.4f} | Val {val_loss:.4f} | AUC {val_auc:.4f}")

    if val_auc > best_val_auc:
        best_val_auc = val_auc
        patience_counter = 0
        best_state = {k: v.clone() for k, v in model.state_dict().items()}
    else:
        patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"Early stopping at epoch {epoch+1}.")
            break

model.load_state_dict(best_state)
print(f"\nBest val AUC: {best_val_auc:.4f}")

# ========= 13. 测试 =========
model.eval()
with torch.no_grad():
    scores = torch.sigmoid(model(X_test_t)).numpy().flatten()
auc = roc_auc_score(y_test, scores, sample_weight=sw_test)
print(f"Test AUC (weighted): {auc:.4f}")

np.save(f"{out_dir}/train_losses.npy", np.array(train_losses))
np.save(f"{out_dir}/val_losses.npy",   np.array(val_losses))
np.save(f"{out_dir}/val_aucs.npy",     np.array(val_aucs))
np.save(f"{out_dir}/scores.npy",       scores)
np.save(f"{out_dir}/y_test.npy",       y_test)
# Relative weights only, not absolute yields: within signal it is the true DsW:DsstW
# mixture (both 1.0 by construction) and within background the true jet-bin mixture,
# but signal and background are on different, unrelated scales -- do not compute a
# cross-class significance from this.
np.save(f"{out_dir}/w_test.npy",       w[test_i])

# ========= 14. 结果图 =========
fig, axes = plt.subplots(2, 3, figsize=(21, 12))
axes = axes.flatten()

axes[0].plot(train_losses, label="Train Loss"); axes[0].plot(val_losses, label="Val Loss")
axes[0].set_xlabel("Epoch"); axes[0].set_ylabel("Loss"); axes[0].legend()
cms_label(axes[0])

axes[1].plot(val_aucs, color="steelblue")
axes[1].set_xlabel("Epoch"); axes[1].set_ylabel("Val AUC")
cms_label(axes[1])

fpr, tpr, _ = roc_curve(y_test, scores, sample_weight=sw_test)
axes[2].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[2].plot([0, 1], [0, 1], "k--")
axes[2].set_xlabel("False Positive Rate"); axes[2].set_ylabel("True Positive Rate"); axes[2].legend()
cms_label(axes[2])

bins = np.linspace(0, 1, 51)
for ax, logy in ((axes[3], False), (axes[4], True)):
    ax.hist(scores[y_test == 0], bins=bins, weights=sw_test[y_test == 0], alpha=0.6,
            color="steelblue", label=BKG_LABEL, density=True)
    ax.hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
            color="firebrick", label="Signal (DsW+DsstW)", density=True)
    ax.set_xlabel("DNN Score"); ax.set_ylabel("a.u."); ax.legend()
    if logy:
        ax.set_yscale("log")
    cms_label(ax)

axes[5].set_visible(False)

plt.tight_layout()
plt.savefig(common.fig_path("training", "results.png", key=BKG, tool="DNN", seed=SEED),
            dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= 15. 特征重要性 =========
print("Computing feature importance...")
importances = np.zeros(len(feats.names))
model.eval()
for i in range(len(feats.names)):
    X_perm = X_test.copy()
    X_perm[:, i] = rng.permutation(X_perm[:, i])
    with torch.no_grad():
        perm_scores = torch.sigmoid(model(to_tensor(X_perm))).numpy().flatten()
    importances[i] = auc - roc_auc_score(y_test, perm_scores, sample_weight=sw_test)

np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)
colors = ["firebrick" if importances[i] > 0 else "steelblue" for i in sorted_idx]
fig, ax = plt.subplots(figsize=(8, len(feats.names) * 0.3))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color=colors)
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([feats.names[i] for i in sorted_idx], fontsize=11)
ax.set_xlabel("AUC drop when feature permuted")
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
cms_label(ax)
plt.tight_layout()
plt.savefig(common.fig_path("training", "feature_importance.png",
                            key=BKG, tool="DNN", seed=SEED), dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 16. 保存模型 =========
torch.save(model.state_dict(), f"{out_dir}/phipimu_model.pt")
# Stamp the exact feature list (and order) this model was trained on. predict.py
# refuses to score against a dir whose stamp disagrees with common.py -- without it,
# a feature added or reordered later would silently produce wrong scores.
common.save_feature_names(out_dir, feats)
print(f"Done. ({len(feats.names)} features, background={BKG})")

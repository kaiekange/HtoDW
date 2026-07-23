"""
DNN: signal vs W+jets / TT / DY, phipimu channel.

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
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages

import mplhep as hep
plt.style.use(hep.style.CMS)

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

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
out_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/{model_name}"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)

# ========= 2. 特征定义 =========
# Must stay identical to BDT/BDT_phipimu.py.
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


# ========= 5. 特征变换 =========
log_idx = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]
X[:, log_idx]  = np.log1p(X[:, log_idx])
X[:, dira_idx] = np.log(1.01 - X[:, dira_idx])

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 6. 训练权重 =========
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

model = phipimuNet(X_train.shape[1])
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
np.save(f"{out_dir}/w_test.npy",       w[test_i])     # physical Lumi_reweight, for yields

# ========= 14. 结果图 =========
fig, axes = plt.subplots(2, 3, figsize=(21, 12))
axes = axes.flatten()

axes[0].plot(train_losses, label="Train Loss"); axes[0].plot(val_losses, label="Val Loss")
axes[0].set_xlabel("Epoch"); axes[0].set_ylabel("Loss"); axes[0].legend()
hep.cms.label("Simulation", ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

axes[1].plot(val_aucs, color="steelblue")
axes[1].set_xlabel("Epoch"); axes[1].set_ylabel("Val AUC")
hep.cms.label("Simulation", ax=axes[1], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

fpr, tpr, _ = roc_curve(y_test, scores, sample_weight=sw_test)
axes[2].plot(fpr, tpr, label=f"AUC = {auc:.4f}"); axes[2].plot([0, 1], [0, 1], "k--")
axes[2].set_xlabel("False Positive Rate"); axes[2].set_ylabel("True Positive Rate"); axes[2].legend()
hep.cms.label("Simulation", ax=axes[2], rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

bins = np.linspace(0, 1, 51)
for ax, logy in ((axes[3], False), (axes[4], True)):
    ax.hist(scores[y_test == 0], bins=bins, weights=sw_test[y_test == 0], alpha=0.6,
            color="steelblue", label="Background (4)", density=True)
    ax.hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
            color="firebrick", label="Signal", density=True)
    ax.set_xlabel("DNN Score"); ax.set_ylabel("a.u."); ax.legend()
    if logy:
        ax.set_yscale("log")
    hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)

axes[5].set_visible(False)

plt.tight_layout()
plt.savefig(f"{out_dir}/results.png", dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= 15. 特征重要性 =========
print("Computing feature importance...")
importances = np.zeros(len(FEATURE_NAMES))
model.eval()
for i in range(len(FEATURE_NAMES)):
    X_perm = X_test.copy()
    X_perm[:, i] = rng.permutation(X_perm[:, i])
    with torch.no_grad():
        perm_scores = torch.sigmoid(model(to_tensor(X_perm))).numpy().flatten()
    importances[i] = auc - roc_auc_score(y_test, perm_scores, sample_weight=sw_test)

np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)
colors = ["firebrick" if importances[i] > 0 else "steelblue" for i in sorted_idx]
fig, ax = plt.subplots(figsize=(8, len(FEATURE_NAMES) * 0.35))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color=colors)
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([FEATURE_NAMES[i] for i in sorted_idx], fontsize=13)
ax.set_xlabel("AUC drop when feature permuted")
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=13)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 16. 保存模型 =========
torch.save(model.state_dict(), f"{out_dir}/phipimu_model.pt")
print("Done.")

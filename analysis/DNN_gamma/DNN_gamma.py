"""
DNN: true gamma vs fake gamma
Signal:  DsstW _best rows with Ds_match && mu_match && gamma_match
Background: DsstW _best rows with Ds_match && mu_match && has_gamma && !gamma_match
          + DsW   _best rows with Ds_match && mu_match && has_gamma

One gamma randomly sampled per event. Train/val/test split by event_id.
Input files are the output of best_cand.py (one Ds per event, all gammas kept).
"""
import sys
import numpy as np
import uproot
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

import mplhep as hep
plt.style.use(hep.style.CMS)

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

# ========= 1. 配置 =========
DSSTW_FILE = "/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/v1/Prediction/Signal_DsstW_H_KKpimu_best.root"
DSW_FILE   = "/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/v1/Prediction/Signal_H_KKpimu_best.root"
DSSTW_TREE = "HKKpimuDsstW/Events"
DSW_TREE   = "HKKpimuDsW/Events"

model_name = sys.argv[1] if len(sys.argv) > 1 else "v1"
out_dir = f"/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/{model_name}"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)

# ========= 2. 特征定义 =========
# Focus: Ds-gamma relationship
FEATURE_NAMES = [
    # gamma kinematics
    "gamma_pt",
    # angular distances: gamma vs Ds and daughters
    "dR_Ds_gamma",
    "dR_Kp_gamma",
    "dR_Km_gamma",
    "dR_pi_gamma",
    "dR_gamma_mu",
    "dR_Dsst_mu",
    # Ds* system (all change with gamma choice)
    "Dsst_pt",
    "Dsst_p",
    "Dsst_IsoR04_RelIso",
    "Dsst_IsoR04_sumPhotonEt",
    # Ds*+mu system
    "Dsstmu_pt",
    "Dsstmu_eta",
]

ALL_FEATURES = FEATURE_NAMES

LOAD_BRANCHES = FEATURE_NAMES + ["event_id", "Ds_match", "mu_match",
                                  "gamma_match", "has_gamma"]

LOG_FEATURES = []

# ========= 3. 数据加载函数 =========
def load_rows(fname, tree_name, sig_mask_fn):
    """Load rows satisfying sig_mask_fn(arrs) -> bool array."""
    with uproot.open(fname)[tree_name] as t:
        arrs = t.arrays(LOAD_BRANCHES, library="np")
    mask = sig_mask_fn(arrs)
    return {b: arrs[b][mask] for b in LOAD_BRANCHES}

def one_per_event(arrs):
    """Randomly pick one row per event_id. Returns index array."""
    ev_to_rows = defaultdict(list)
    for i, ev in enumerate(arrs["event_id"]):
        ev_to_rows[int(ev)].append(i)
    chosen = []
    for rows in ev_to_rows.values():
        chosen.append(rng.choice(rows))
    return np.array(chosen, dtype=int)

def to_matrix(arrs, indices=None):
    cols = []
    for f in ALL_FEATURES:
        v = arrs[f] if indices is None else arrs[f][indices]
        cols.append(v)
    return np.column_stack(cols)

# ========= 4. 加载信号和本底 =========
print("Loading DsstW signal (gamma_match)...")
sig_arrs = load_rows(
    DSSTW_FILE, DSSTW_TREE,
    lambda a: a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) & a["gamma_match"].astype(bool)
)
sig_idx  = one_per_event(sig_arrs)
X_sig    = to_matrix(sig_arrs, sig_idx)
ev_sig   = sig_arrs["event_id"][sig_idx]

print("Loading DsstW background (Ds_match && mu_match && has_gamma && !gamma_match)...")
bkg1_arrs = load_rows(
    DSSTW_FILE, DSSTW_TREE,
    lambda a: (a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) &
               a["has_gamma"].astype(bool) & ~a["gamma_match"].astype(bool))
)
bkg1_idx = one_per_event(bkg1_arrs)
X_bkg1   = to_matrix(bkg1_arrs, bkg1_idx)
ev_bkg1  = bkg1_arrs["event_id"][bkg1_idx]

print("Loading DsW background (Ds_match && mu_match && has_gamma)...")
bkg2_arrs = load_rows(
    DSW_FILE, DSW_TREE,
    lambda a: (a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) &
               a["has_gamma"].astype(bool))
)
bkg2_idx = one_per_event(bkg2_arrs)
X_bkg2   = to_matrix(bkg2_arrs, bkg2_idx)
ev_bkg2  = bkg2_arrs["event_id"][bkg2_idx]

X_bkg = np.vstack([X_bkg1, X_bkg2])
ev_bkg = np.concatenate([ev_bkg1, ev_bkg2])

print(f"Signal:        {len(X_sig)}")
print(f"Bkg (DsstW):   {len(X_bkg1)}")
print(f"Bkg (DsW):     {len(X_bkg2)}")
print(f"Bkg total:     {len(X_bkg)}")

# ========= 4b. bkg1 vs bkg2 诊断图 =========
print("Drawing bkg1 vs bkg2 diagnostic plots...")
m1 = np.isfinite(X_bkg1).all(axis=1)
m2 = np.isfinite(X_bkg2).all(axis=1)
Xb1_diag, Xb2_diag = X_bkg1[m1], X_bkg2[m2]

with PdfPages(f"{out_dir}/bkg_comparison.pdf") as pdf:
    n_feats = len(ALL_FEATURES)
    N_COLS_D, N_ROWS_D = 4, 4
    PER_PAGE_D = N_COLS_D * N_ROWS_D
    for page_start in range(0, n_feats, PER_PAGE_D):
        page_feats = ALL_FEATURES[page_start:page_start + PER_PAGE_D]
        fig, axes = plt.subplots(N_ROWS_D, N_COLS_D, figsize=(20, 20))
        axes_flat = axes.flatten()
        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            fi = ALL_FEATURES.index(feat)
            b1_vals = Xb1_diag[:, fi]
            b2_vals = Xb2_diag[:, fi]
            all_vals = np.concatenate([b1_vals, b2_vals])
            lo, hi = np.percentile(all_vals, [1, 99])
            if lo == hi:
                lo, hi = all_vals.min(), all_vals.max()
            bins = np.linspace(lo, hi, 50)
            ax.hist(b1_vals, bins=bins, density=True, alpha=0.6,
                    color="darkorange", label="DsstW fake γ (bkg1)")
            ax.hist(b2_vals, bins=bins, density=True, alpha=0.6,
                    color="steelblue",  label="DsW any γ (bkg2)")
            ax.set_xlabel(feat)
            ax.set_ylabel("a.u.")
            ax.legend(fontsize=10)
            if i == 0:
                hep.cms.label("Simulation", ax=ax, rlabel="2022 (13.6 TeV)",
                              data=False, fontsize=10)
        for j in range(len(page_feats), len(axes_flat)):
            axes_flat[j].set_visible(False)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
print("Saved bkg_comparison.pdf")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
ev_all = np.concatenate([ev_sig, ev_bkg])

# ========= 5. 有限值过滤 =========
mask = np.isfinite(X).all(axis=1)
X, y, ev_all = X[mask], y[mask], ev_all[mask]
print(f"After finite filter: {X.shape}")

# ========= 6. 特征变换 =========
log_idx = [ALL_FEATURES.index(f) for f in LOG_FEATURES]
X[:, log_idx] = np.log1p(X[:, log_idx])

np.save(f"{out_dir}/X_sig.npy", X[y == 1])
np.save(f"{out_dir}/X_bkg.npy", X[y == 0])

# ========= 7. 按 event_id 划分 train/val/test =========
unique_ev = np.unique(ev_all)
rng2 = np.random.default_rng(SEED)
rng2.shuffle(unique_ev)
n = len(unique_ev)
n_test = int(0.2 * n)
n_val  = int(0.1 * n)
ev_test = set(unique_ev[:n_test])
ev_val  = set(unique_ev[n_test:n_test + n_val])
ev_train = set(unique_ev[n_test + n_val:])

train_mask = np.array([ev in ev_train for ev in ev_all])
val_mask   = np.array([ev in ev_val   for ev in ev_all])
test_mask  = np.array([ev in ev_test  for ev in ev_all])

X_train, y_train = X[train_mask], y[train_mask]
X_val,   y_val   = X[val_mask],   y[val_mask]
X_test,  y_test  = X[test_mask],  y[test_mask]
print(f"Train: {len(X_train)} | Val: {len(X_val)} | Test: {len(X_test)}")

# ========= 8. 变量分布对比图 =========
print("Drawing variable comparison plots...")
X_sig_plot = X[y == 1]
X_bkg_plot = X[y == 0]
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
            s_vals = X_sig_plot[:, fi]
            b_vals = X_bkg_plot[:, fi]
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
X_test_t              = to_tensor(X_test)

# ========= 11. 模型 =========
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

model = GammaNet(X_train.shape[1])
print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

# ========= 12. 训练 =========
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
criterion = nn.BCEWithLogitsLoss()
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

EPOCHS, PATIENCE, BATCH_SIZE = 200, 30, 1024

n_sig_tr = int(np.sum(y_train == 1))
n_bkg_tr = int(np.sum(y_train == 0))
print(f"Train signal: {n_sig_tr} | Train background: {n_bkg_tr}")
sample_weights = np.where(y_train == 1, 1.0 / n_sig_tr, 1.0 / n_bkg_tr)
sampler = torch.utils.data.WeightedRandomSampler(
    torch.tensor(sample_weights, dtype=torch.float64), len(y_train), replacement=True
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
        val_pred = model(X_val_t)
        val_loss = criterion(val_pred, y_val_t).item()
        val_scores = torch.sigmoid(val_pred).numpy().flatten()
    val_auc = roc_auc_score(y_val, val_scores)

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
auc = roc_auc_score(y_test, scores)
print(f"Test AUC: {auc:.4f}")

np.save(f"{out_dir}/train_losses.npy", np.array(train_losses))
np.save(f"{out_dir}/val_losses.npy",   np.array(val_losses))
np.save(f"{out_dir}/val_aucs.npy",     np.array(val_aucs))
np.save(f"{out_dir}/scores.npy",       scores)
np.save(f"{out_dir}/y_test.npy",       y_test)

# ========= 14. 结果图 =========
fig, axes = plt.subplots(2, 3, figsize=(21, 12))
axes = axes.flatten()

axes[0].plot(train_losses, label="Train Loss"); axes[0].plot(val_losses, label="Val Loss")
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

# ========= 15. 特征重要性 =========
print("Computing feature importance...")
importances = np.zeros(len(ALL_FEATURES))
model.eval()
for i in range(len(ALL_FEATURES)):
    X_perm = X_test.copy()
    X_perm[:, i] = rng.permutation(X_perm[:, i])
    with torch.no_grad():
        perm_scores = torch.sigmoid(model(to_tensor(X_perm))).numpy().flatten()
    importances[i] = auc - roc_auc_score(y_test, perm_scores)

np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)
colors = ["firebrick" if importances[i] > 0 else "steelblue" for i in sorted_idx]
fig, ax = plt.subplots(figsize=(8, len(ALL_FEATURES) * 0.35))
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

# ========= 16. 保存模型 =========
torch.save(model.state_dict(), f"{out_dir}/gamma_model.pt")
print("Done.")

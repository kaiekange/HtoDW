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

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

# ========= 1. 参数 =========
sig_file = "/eos/user/k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/Signal.root"
bkg_file = "/eos/user/k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/Background.root"
tree_name = "Events"
model_name = sys.argv[1]

out_dir = f"/eos/user/k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/DNN/{model_name}"
os.makedirs(out_dir, exist_ok=True)

FEATURE_NAMES = [
    # MET
    "met_sig",
    "met_pt",
    "met_sumEt",
    # jets
    "num_jets_pt30",
    # muon
    "mu_eta",
    "mu_pt",
    "mu_p",
    # K+
    "Kp_eta",
    "Kp_pt",
    "Kp_p",
    "Kp_PVassoc",
    # K-
    "Km_eta",
    "Km_pt",
    "Km_p",
    "Km_PVassoc",
    # pi
    "pi_eta",
    "pi_pt",
    "pi_p",
    "pi_PVassoc",
    # KK system
    "KK_eta",
    "KK_pt",
    "KK_p",
    # K+pi system
    "Kppi_eta",
    "Kppi_pt",
    "Kppi_p",
    # K-pi system
    "Kmpi_eta",
    "Kmpi_pt",
    "Kmpi_p",
    # Ds
    "Ds_eta",
    "Ds_pt",
    "Ds_p",
    # vis (Ds + mu)
    "vis_eta",
    "vis_pt",
    "vis_p",
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
    # dR (muon)
    "dR_Kp_mu",
    "dR_Km_mu",
    "dR_pi_mu",
    "dR_Ds_mu",
    # IP significances
    "Kp_ipsig",
    "Kp_ipxysig",
    "Km_ipsig",
    "Km_ipxysig",
    "pi_ipsig",
    "pi_ipxysig",
    "Ds_ipsig",
    "Ds_ipxysig",
    "mu_ipsig",
    "mu_ipxysig",
    # Ds flight distance
    "Ds_FDsig",
    "Ds_FDxysig",
    # isolation
    "mu_IsoR04_RelIso",
    "Ds_IsoR04_RelIso",
]

# log 变换的特征（长尾分布）
LOG_FEATURES = [
    "Kp_ipsig",
    "Kp_ipxysig",
    "Km_ipsig",
    "Km_ipxysig",
    "pi_ipsig",
    "pi_ipxysig",
    "Ds_ipsig",
    "Ds_ipxysig",
    "mu_ipsig",
    "mu_ipxysig",
    "Ds_FDsig",
    "Ds_FDxysig",
]


# ========= 2. 读取数据 =========
print("Loading ROOT files...")
sig_tree = uproot.open(sig_file)[tree_name]
bkg_tree = uproot.open(bkg_file)[tree_name]

sig_data = sig_tree.arrays(FEATURE_NAMES, library="np")
bkg_data = bkg_tree.arrays(FEATURE_NAMES, library="np")

X_sig = np.column_stack([sig_data[b] for b in FEATURE_NAMES])
X_bkg = np.column_stack([bkg_data[b] for b in FEATURE_NAMES])
X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])

print(f"Signal: {len(X_sig)} | Background: {len(X_bkg)}")
print(f"Data shape: {X.shape}")
print(f"Signal fraction: {np.mean(y):.4f}")

# ========= 3. 清理数据 =========
mask = np.isfinite(X).all(axis=1)
X = X[mask]
y = y[mask]
print(f"After finite filter: {X.shape}")

# ========= 4. Log 变换 =========
log_idx = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]
X[:, log_idx] = np.log1p(X[:, log_idx])

# dira 单独处理：log(1.01 - dira)
dira_idx = FEATURE_NAMES.index("dira")
X[:, dira_idx] = np.log(1.01 - X[:, dira_idx])

# ========= 4b. 变量分布对比 =========
print("Drawing variable comparison plots...")
X_sig_plot = X[y == 1]
X_bkg_plot = X[y == 0]

N_COLS, N_ROWS = 5, 4
PER_PAGE = N_COLS * N_ROWS

with PdfPages(f"{out_dir}/variable_comparison.pdf") as pdf:
    for page_start in range(0, len(FEATURE_NAMES), PER_PAGE):
        page_feats = FEATURE_NAMES[page_start : page_start + PER_PAGE]
        fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(25, 20))
        axes_flat = axes.flatten()

        for i, feat in enumerate(page_feats):
            ax = axes_flat[i]
            idx = FEATURE_NAMES.index(feat)
            s_vals = X_sig_plot[:, idx]
            b_vals = X_bkg_plot[:, idx]

            all_vals = np.concatenate([s_vals, b_vals])
            lo, hi = np.percentile(all_vals, [1, 99])
            if lo == hi:
                lo, hi = all_vals.min(), all_vals.max()
            bins = np.linspace(lo, hi, 50)

            ax.hist(b_vals, bins=bins, density=True, alpha=0.6, color="steelblue", label="Background")
            ax.hist(s_vals, bins=bins, density=True, alpha=0.6, color="firebrick", label="Signal")
            ax.set_xlabel(feat)
            ax.set_ylabel("a.u.")
            ax.legend(fontsize=12)
            if i == 0:
                hep.cms.label("Preliminary", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=12)

        for j in range(len(page_feats), len(axes_flat)):
            axes_flat[j].set_visible(False)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print("Saved variable_comparison.pdf")

np.save(f"{out_dir}/X_sig.npy", X_sig_plot)
np.save(f"{out_dir}/X_bkg.npy", X_bkg_plot)

# ========= 5. 划分数据 =========
# stratify=y preserves the signal fraction in every split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1, random_state=42, stratify=y_train)

print(f"Train: {X_train.shape[0]} | Val: {X_val.shape[0]} | Test: {X_test.shape[0]}")

# ========= 6. 标准化 =========
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)
X_test = scaler.transform(X_test)


# ========= 7. 转 tensor =========
def to_tensor(X, y=None):
    X_t = torch.tensor(X, dtype=torch.float32)
    if y is not None:
        y_t = torch.tensor(y.reshape(-1, 1), dtype=torch.float32)
        return X_t, y_t
    return X_t


X_train_t, y_train_t = to_tensor(X_train, y_train)
X_val_t, y_val_t = to_tensor(X_val, y_val)
X_test_t = to_tensor(X_test)


# ========= 8. 定义模型 =========
class DsNet(nn.Module):
    def __init__(self, n_features, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(32, 16),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Linear(16, 1),
            # no Sigmoid — BCEWithLogitsLoss expects raw logits
        )

    def forward(self, x):
        return self.net(x)


model = DsNet(X_train.shape[1])
print(f"\nModel parameters: {sum(p.numel() for p in model.parameters()):,}")

# ========= 9. 训练配置 =========
# AdamW: weight decay is decoupled from the adaptive lr update (unlike Adam+L2)
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
# balanced sampling handles class imbalance → standard BCE (no pos_weight needed)
criterion = nn.BCEWithLogitsLoss()
# scheduler tracks val_loss (smooth signal); early stopping tracks val_AUC (decoupled)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

EPOCHS = 200
PATIENCE = 30
BATCH_SIZE = 1024

# WeightedRandomSampler: signal and background contribute equally to every batch
n_sig = int(np.sum(y_train == 1))
n_bkg = int(np.sum(y_train == 0))
print(f"Train signal: {n_sig} | Train background: {n_bkg}")
sample_weights = np.where(y_train == 1, 1.0 / n_sig, 1.0 / n_bkg)
sampler = torch.utils.data.WeightedRandomSampler(
    weights=torch.tensor(sample_weights, dtype=torch.float64),
    num_samples=len(y_train),
    replacement=True,
)

best_val_auc = -1.0
patience_counter = 0
best_state = None

train_losses, val_losses, val_aucs = [], [], []

# ========= 10. 训练循环 =========
print("\nStart training...")
dataset = torch.utils.data.TensorDataset(X_train_t, y_train_t)
loader = torch.utils.data.DataLoader(dataset, batch_size=BATCH_SIZE, sampler=sampler)

for epoch in range(EPOCHS):
    # --- 训练 ---
    model.train()
    epoch_loss = 0.0
    for xb, yb in loader:
        pred = model(xb)
        loss = criterion(pred, yb)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item() * len(xb)
    epoch_loss /= len(X_train_t)

    # --- 验证 ---
    model.eval()
    with torch.no_grad():
        val_pred = model(X_val_t)
        val_loss = criterion(val_pred, y_val_t).item()
        val_scores_ep = torch.sigmoid(val_pred).numpy().flatten()
    val_auc = roc_auc_score(y_val, val_scores_ep)

    train_losses.append(epoch_loss)
    val_losses.append(val_loss)
    val_aucs.append(val_auc)

    # scheduler steps on val_loss (smooth per-epoch metric)
    scheduler.step(val_loss)

    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1:03d} | Train Loss = {epoch_loss:.4f} | Val Loss = {val_loss:.4f} | Val AUC = {val_auc:.4f}")

    # early stopping on val_AUC (decoupled from scheduler's val_loss metric)
    if val_auc > best_val_auc:
        best_val_auc = val_auc
        patience_counter = 0
        best_state = {k: v.clone() for k, v in model.state_dict().items()}
    else:
        patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"\nEarly stopping at epoch {epoch+1}.")
            break

model.load_state_dict(best_state)
print(f"\nBest val AUC: {best_val_auc:.4f}")

np.save(f"{out_dir}/train_losses.npy", np.array(train_losses))
np.save(f"{out_dir}/val_losses.npy", np.array(val_losses))
np.save(f"{out_dir}/val_aucs.npy", np.array(val_aucs))

# ========= 11. 测试 =========
model.eval()
with torch.no_grad():
    scores = torch.sigmoid(model(X_test_t)).numpy().flatten()

auc = roc_auc_score(y_test, scores)
print(f"Test AUC = {auc:.4f}")

# ========= 12. 训练结果图 =========
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

# Loss curves
axes[0].plot(train_losses, label="Train Loss")
axes[0].plot(val_losses, label="Val Loss")
axes[0].set_xlabel("Epoch")
axes[0].set_ylabel("Loss")
axes[0].legend()
hep.cms.label("Preliminary", ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, fontsize=12)

# Val AUC curve
axes[1].plot(val_aucs, color="steelblue")
axes[1].set_xlabel("Epoch")
axes[1].set_ylabel("Val AUC")
hep.cms.label("Preliminary", ax=axes[1], rlabel="2022 (13.6 TeV)", data=False, fontsize=12)

# ROC curve
fpr, tpr, _ = roc_curve(y_test, scores)
axes[2].plot(fpr, tpr, label=f"AUC = {auc:.4f}")
axes[2].plot([0, 1], [0, 1], "k--")
axes[2].set_xlabel("False Positive Rate")
axes[2].set_ylabel("True Positive Rate")
axes[2].legend()
hep.cms.label("Preliminary", ax=axes[2], rlabel="2022 (13.6 TeV)", data=False, fontsize=12)

# Score distribution
sig_scores = scores[y_test == 1]
bkg_scores = scores[y_test == 0]
axes[3].hist(bkg_scores, bins=50, alpha=0.6, color="steelblue", label="Background", density=True)
axes[3].hist(sig_scores, bins=50, alpha=0.6, color="firebrick", label="Signal", density=True)
axes[3].set_xlabel("DNN Score")
axes[3].set_ylabel("a.u.")
axes[3].legend()
hep.cms.label("Preliminary", ax=axes[3], rlabel="2022 (13.6 TeV)", data=False, fontsize=12)

plt.tight_layout()
plt.savefig(f"{out_dir}/results.png", dpi=150)
plt.close(fig)
print("Saved results.png")

# ========= 13. 特征重要性 =========
print("Computing permutation feature importance...")
rng = np.random.default_rng(seed=42)
importances = np.zeros(len(FEATURE_NAMES))

model.eval()
for i in range(len(FEATURE_NAMES)):
    X_perm = X_test.copy()
    X_perm[:, i] = rng.permutation(X_perm[:, i])
    with torch.no_grad():
        perm_scores = torch.sigmoid(model(to_tensor(X_perm))).numpy().flatten()
    importances[i] = auc - roc_auc_score(y_test, perm_scores)

np.save(f"{out_dir}/feature_importances.npy", importances)
sorted_idx = np.argsort(importances)  # ascending → most important at top of barh
colors = ["firebrick" if importances[i] > 0 else "steelblue" for i in sorted_idx]

fig, ax = plt.subplots(figsize=(8, len(FEATURE_NAMES) * 0.28))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color=colors)
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([FEATURE_NAMES[i] for i in sorted_idx], fontsize=8)
ax.set_xlabel("AUC drop when feature permuted")
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
hep.cms.label("Preliminary", ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=12)
plt.tight_layout()
plt.savefig(f"{out_dir}/feature_importance.png", dpi=150)
plt.close(fig)
print("Saved feature_importance.png")

# ========= 14. 保存 =========
torch.save(model.state_dict(), f"{out_dir}/ds_model.pt")
np.save(f"{out_dir}/scaler_mean.npy", scaler.mean_)
np.save(f"{out_dir}/scaler_scale.npy", scaler.scale_)
np.save(f"{out_dir}/scores.npy", scores)
np.save(f"{out_dir}/y_test.npy", y_test)

print("Done.")

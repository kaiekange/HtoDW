import sys
import numpy as np
import uproot
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import os

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

# ========= 1. 参数 =========
root_file = f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/MVAflatten.root"
tree_name = "Events"
model_name = sys.argv[1]

branches = [
    "Kp_eta",
    "Kp_phi",
    "Kp_pt",
    "Kp_p",
    "Km_eta",
    "Km_phi",
    "Km_pt",
    "Km_p",
    "pi_eta",
    "pi_phi",
    "pi_pt",
    "pi_p",
    "KK_eta",
    "KK_phi",
    "KK_pt",
    "KK_p",
    "Kppi_eta",
    "Kppi_phi",
    "Kppi_pt",
    "Kppi_p",
    "Kmpi_eta",
    "Kmpi_phi",
    "Kmpi_pt",
    "Kmpi_p",
    "Ds_eta",
    "Ds_phi",
    "Ds_pt",
    "Ds_p",
    "Ds_prob",
    "dira",
    "dR_Kp_Km",
    "dR_Kp_pi",
    "dR_Km_pi",
    "dR_Kp_Ds",
    "dR_Km_Ds",
    "dR_pi_Ds",
    "Kp_ipsig",
    "Kp_ipxysig",
    "Km_ipsig",
    "Km_ipxysig",
    "pi_ipsig",
    "pi_ipxysig",
    "Ds_ipsig",
    "Ds_ipxysig",
    "Ds_FDsig",
    "Ds_FDxysig",
    "isSig",
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
    "Ds_FDsig",
    "Ds_FDxysig",
]

FEATURE_NAMES = branches[:-1]


# ========= 2. 读取数据 =========
print("Loading ROOT file...")
file = uproot.open(root_file)
tree = file[tree_name]
data = tree.arrays(branches, library="np")

X = np.column_stack([data[b] for b in FEATURE_NAMES])
y = data["isSig"]

print(f"Data shape: {X.shape}")
print(f"Signal fraction: {np.mean(y):.4f}")

# ========= 3. 清理数据 =========
mask = np.isfinite(X).all(axis=1)
X = X[mask]
y = y[mask]
print(f"After finite filter: {X.shape}")

# ========= 4. Log 变换 =========
log_idx = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]
# 取绝对值后 log(1+x)，避免负值和零
X[:, log_idx] = np.log1p(np.abs(X[:, log_idx]))

# dira 单独处理：log(1.01 - dira)
dira_idx = FEATURE_NAMES.index("dira")
X[:, dira_idx] = np.log(1.01 - X[:, dira_idx])

# ========= 5. 划分数据 =========
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1, random_state=42)

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
            nn.Linear(n_features, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            # nn.Linear(n_features, 128),
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
            nn.Linear(32, 1),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.net(x)


model = DsNet(X_train.shape[1])
print(f"\nModel parameters: {sum(p.numel() for p in model.parameters()):,}")

# ========= 9. 训练配置 =========
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
criterion = nn.BCELoss()
# scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5, verbose=True)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)

# 早停参数
EPOCHS = 60
PATIENCE = 15
BATCH_SIZE = 1024

best_val_loss = float("inf")
patience_counter = 0
best_state = None

train_losses, val_losses = [], []

# ========= 10. 训练循环 =========
print("\nStart training...")
dataset = torch.utils.data.TensorDataset(X_train_t, y_train_t)
loader = torch.utils.data.DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)

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

    train_losses.append(epoch_loss)
    val_losses.append(val_loss)
    scheduler.step(val_loss)

    if (epoch + 1) % 10 == 0:
        print(f"Epoch {epoch+1:03d} | Train Loss = {epoch_loss:.4f} | Val Loss = {val_loss:.4f}")

    # --- 早停 ---
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        patience_counter = 0
        best_state = {k: v.clone() for k, v in model.state_dict().items()}
    else:
        patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"\nEarly stopping at epoch {epoch+1}.")
            break

# 加载最优权重
model.load_state_dict(best_state)
print(f"\nBest val loss: {best_val_loss:.4f}")

# ========= 11. 测试 =========
model.eval()
with torch.no_grad():
    scores = model(X_test_t).numpy().flatten()

auc = roc_auc_score(y_test, scores)
print(f"Test AUC = {auc:.4f}")

# ========= 12. 画图 =========
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Loss 曲线
axes[0].plot(train_losses, label="Train")
axes[0].plot(val_losses, label="Val")
axes[0].set_xlabel("Epoch")
axes[0].set_ylabel("Loss")
axes[0].set_title("Loss Curve")
axes[0].legend()

# ROC 曲线
fpr, tpr, _ = roc_curve(y_test, scores)
axes[1].plot(fpr, tpr, label=f"AUC = {auc:.4f}")
axes[1].plot([0, 1], [0, 1], "k--")
axes[1].set_xlabel("False Positive Rate")
axes[1].set_ylabel("True Positive Rate")
axes[1].set_title("ROC Curve")
axes[1].legend()

# Score 分布
sig_scores = scores[y_test == 1]
bkg_scores = scores[y_test == 0]
axes[2].hist(bkg_scores, bins=50, alpha=0.6, label="Background", density=True)
axes[2].hist(sig_scores, bins=50, alpha=0.6, label="Signal", density=True)
axes[2].set_xlabel("Score")
axes[2].set_ylabel("Density")
axes[2].set_title("Score Distribution")
axes[2].legend()

plt.tight_layout()
plt.savefig(f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/results.png", dpi=150)
print("Saved results.png")

# ========= 13. 保存 =========
torch.save(model.state_dict(), f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/ds_model.pt")
np.save(f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/scaler_mean.npy", scaler.mean_)
np.save(f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/scaler_scale.npy", scaler.scale_)
np.save(f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/scores.npy", scores)
np.save(f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/{model_name}/y_test.npy", y_test)

print("Done.")

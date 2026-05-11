import sys
import numpy as np
import uproot
import torch
import torch.nn as nn
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

# ========= 1. 参数 =========
root_file = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/MVAflatten.root"
tree_name = "Events"

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


# ========= 2. 读取数据 =========
print("Loading ROOT file...")
file = uproot.open(root_file)
tree = file[tree_name]

data = tree.arrays(branches, library="np")

X = np.column_stack([data[b] for b in branches[:-1]])
y = data["isSig"]

print("Data shape:", X.shape)
print("Signal fraction:", np.mean(y))

# ========= 3. 清理数据 =========
mask = np.isfinite(X).all(axis=1)
X = X[mask]
y = y[mask]

# ========= 4. 划分数据 =========
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# ========= 5. 标准化 =========
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# ========= 6. 转 tensor =========
X_train_t = torch.tensor(X_train, dtype=torch.float32)
y_train_t = torch.tensor(y_train.reshape(-1, 1), dtype=torch.float32)

X_test_t = torch.tensor(X_test, dtype=torch.float32)


# ========= 7. 定义模型 =========
class DsNet(nn.Module):
    def __init__(self, n_features):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(n_features, 128), nn.ReLU(), nn.Dropout(0.2), nn.Linear(128, 64), nn.ReLU(), nn.Linear(64, 1), nn.Sigmoid())

    def forward(self, x):
        return self.net(x)


model = DsNet(X_train.shape[1])

# ========= 8. 训练配置 =========
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
criterion = nn.BCELoss()

# ========= 9. 训练 =========
print("Start training...")
for epoch in range(20):
    model.train()

    pred = model(X_train_t)
    loss = criterion(pred, y_train_t)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    print(f"Epoch {epoch:02d} | Loss = {loss.item():.4f}")

# ========= 10. 测试 =========
model.eval()
with torch.no_grad():
    scores = model(X_test_t).numpy().flatten()

auc = roc_auc_score(y_test, scores)
print("AUC =", auc)

# ========= 11. 画 ROC =========
fpr, tpr, _ = roc_curve(y_test, scores)

plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.title("ROC Curve")
plt.savefig("roc.png")
print("Saved roc.png")

# ========= 12. 保存模型 =========
torch.save(model.state_dict(), "ds_model.pt")
np.save("scaler_mean.npy", scaler.mean_)
np.save("scaler_scale.npy", scaler.scale_)

print("Saved model and scaler")

# ========= 13. 保存score =========
np.save("scores.npy", scores)
np.save("y_test.npy", y_test)

print("Done.")

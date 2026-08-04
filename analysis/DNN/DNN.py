"""
DNN: H->Ds(*)W signal vs W+jets, phipimu channel.

Both signal channels (DsW + DsstW) against a single background process. The other
backgrounds (TT, DY, ...) are still missing from the new Samples/ format, and they
are known to be the ones this classifier finds hardest -- see the
project_wjets_background_incomplete / project_mva_sculpting_finding notes: the
data/MC gap widens with the score cut, so the rejection quoted here is optimistic.

No longer kept in lockstep with BDT/BDT.py: as of 2026-08-04 the two use different
feature sets (BDT is still on the 42-feature test build) and are not directly
comparable. Feature set here is 40: the 42 minus n_b_jets_{loose,medium}_dR04 and
the three *_PVassoc, plus the three per-track *_ipxysig.

Reads the one-row-per-event Samples tuples directly (flat `Events` tree -- these
files, produced by highest_pt_dedup.py, carry no embedded Lumi_reweight branch,
unlike the older 34-background Samples/ files).

Feature set excludes every decay-kinematics variable: those reconstruct
phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the resonance peaks
used for the final fit. Dsmu_mass is a feature here, not a fit variable -- it has
no resonant peak, and keeping it out of the fit is what lets Ds_pt stay in the MVA.
"""
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
TREE = "Events"

# Both signal channels, each event at weight 1. The two productions were sized in
# proportion to their cross sections (500 vs ~992 PAT files, N_processed 1.54M vs
# 2.86M = 1:1.86), so N_gen is already proportional to sigma x BR and the
# equivalent-luminosity weight per event is the same for both -- concatenating them
# reproduces the physical DsW:DsstW mixture, and their selection efficiencies are
# near identical anyway (9.32% vs 9.49%), so it survives the selection too.
# If the two channels ever need genuinely unequal weights, put them here.
SIGNALS = {
    "HtoDsW":   1.0,
    "HtoDsstW": 1.0,
}

# WtoLNu2Jets is split into three jet-multiplicity bins of the SAME physical
# process, so they are weight-merged using the full statistics of all three
# files. Weights are the equivalent-luminosity reweighting factors (external --
# these files have no embedded Lumi_reweight branch).
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}

out_dir = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets"
os.makedirs(out_dir, exist_ok=True)

SEED = 42
rng = np.random.default_rng(SEED)
torch.manual_seed(SEED)

# ========= 2. 特征定义 =========
# Must stay identical to BDT/BDT.py.
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
    # Per-track 2D IP significance. Added 2026-08-04: the same "transverse beats
    # longitudinal" argument that motivated Ds_ipxysig applies to the tracks, and it
    # had simply never been carried down to them. Measured on this sample they beat
    # their 3D counterparts (AUC 0.725/0.723/0.780 vs 0.710/0.709/0.774) and are only
    # weakly correlated with them (r = 0.14/0.19/0.27) -- the 3D versions are diluted
    # by the much worse z resolution, so the two carry largely complementary
    # information rather than being two views of the same quantity.
    # mu_ipxysig deliberately left out: AUC 0.504, i.e. nothing -- the muon comes from
    # a prompt W and has no displacement to measure (its 3D version, mu_ipsig, is
    # equally useless at 0.503 but is kept for continuity with the earlier feature set).
    "Kp_ipxysig",
    "Km_ipxysig",
    "pi_ipxysig",
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
    "Ds_nearestjet_bscore",
    "Ds_nearestjet_dR",
    "jet_leadpt",
    # n_b_jets_{loose,medium}_dR04 were dropped 2026-08-04: kept in the first test build
    # on the argument that tight subset medium subset loose leaves them free to vary
    # after the tight veto -- true, but they then ranked dead last (41st and 42nd of 42,
    # importance 0.0016/0.0020). "Not degenerate with the cut" did not imply "useful".
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
    "Kp_ipxysig",
    "Km_ipxysig",
    "pi_ipxysig",
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
]

# Sentinel handling -- needed here but NOT in the BDT: trees just split around a
# -777, but StandardScaler would let it wreck the mean/variance of the whole
# feature. Each sentinel is mapped to the physically sensible edge of its range
# ("no jet nearby" / "tagger gave nothing"), keeping the input bounded.
#   Ds_nearestjet_bscore  -777 (no jet) and ~-10 (tagger unavailable, ~14% of rows)
#   Ds_nearestjet_dR      -777 (no jet)  -> far away
#   jet_leadpt            -777 (no jet)  -> zero pt
SENTINEL_FIX = {
    "Ds_nearestjet_bscore": -0.1,   # just below the [0,1] tagger range
    "Ds_nearestjet_dR": 5.0,        # beyond any real dR in the sample
    "jet_leadpt": 0.0,
}

# Exactly one tight muon in the event (the best muon itself, no extra of either
# charge), plus a b-jet veto. The veto is a CUT, not a feature: its job is to
# remove the process the background model does not contain at all (tt, with real
# b jets). Tight WP keeps 99.3% of signal while removing 25% of data; dR04
# matches the AK4 jet cone radius. See BDT/BDT.py for the full reasoning.
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

print(f"Loading signal ({'+'.join(SIGNALS)})...")
sig_X, sig_w = [], []
for name, wgt in SIGNALS.items():
    x = load(name, truth_filter=True)
    sig_X.append(x)
    sig_w.append(np.full(len(x), wgt))
X_sig = np.vstack(sig_X)
w_sig = np.concatenate(sig_w)

print(f"Loading WtoLNu2Jets ({'+'.join(WJETS_WEIGHTS)}, lumi-reweighted)...")
bkg_X, bkg_w = [], []
for name, wgt in WJETS_WEIGHTS.items():
    x = load(name, truth_filter=False)
    bkg_X.append(x)
    bkg_w.append(np.full(len(x), wgt))
X_bkg = np.vstack(bkg_X)
w_bkg = np.concatenate(bkg_w)

print(f"Signal: {len(X_sig):,} (both channels at their physical proportion, see above)")
print(f"Background: {len(X_bkg):,} (WtoLNu2Jets, full stats, lumi-reweighted merge of 3 jet bins)")

X = np.vstack([X_sig, X_bkg])
y = np.concatenate([np.ones(len(X_sig)), np.zeros(len(X_bkg))])
w = np.concatenate([w_sig, w_bkg])

# ========= 4. sentinel 规整 + 有限值 + dira 过滤 =========
for feat, fill in SENTINEL_FIX.items():
    i = FEATURE_NAMES.index(feat)
    bad = X[:, i] < -1
    if bad.any():
        print(f"  sentinel: {feat:24s} {bad.sum():>8,} rows ({100*bad.mean():5.2f}%) -> {fill}")
        X[bad, i] = fill

dira_idx = FEATURE_NAMES.index("dira")
mask = np.isfinite(X).all(axis=1) & (X[:, dira_idx] > -1)
X, y, w = X[mask], y[mask], w[mask]
print(f"After finite filter: {X.shape}")

# Verify the lumi-reweighting actually took effect: weighted mean should match the
# physical jet-bin mixture, not the raw event-count one (which is dominated by 2J
# simply having more events).
jet30_idx = FEATURE_NAMES.index("num_jets_pt30")
check = np.average(X[y == 0][:, jet30_idx], weights=w[y == 0])
print(f"[check] weighted <num_jets_pt30> for background = {check:.3f} (raw-concat would be ~1.55)")

# ========= 5. 特征变换 =========
log_idx = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]
X[:, log_idx]  = np.log1p(X[:, log_idx])
X[:, dira_idx] = np.log(1.01 - X[:, dira_idx])

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
            color="steelblue", label="WtoLNu2Jets", density=True)
    ax.hist(scores[y_test == 1], bins=bins, weights=sw_test[y_test == 1], alpha=0.6,
            color="firebrick", label="HtoDsW", density=True)
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
fig, ax = plt.subplots(figsize=(8, len(FEATURE_NAMES) * 0.3))
ax.barh(range(len(sorted_idx)), importances[sorted_idx], color=colors)
ax.set_yticks(range(len(sorted_idx)))
ax.set_yticklabels([FEATURE_NAMES[i] for i in sorted_idx], fontsize=11)
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

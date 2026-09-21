"""OS vs SS: is the charge selection itself moving the MVA input variables?

    python3 check_os_vs_ss_features.py            # 2022 data (default, both eras)

NanoMVA migration (2026-09-18): re-pointed at NanoMVA/Samples/ (Data_OS_merged,
Data_SS_merged; weight column is now `weight`, not `Lumi_reweight`), both eras
(era=1-only restriction dropped). The `wjets` mode is DROPPED, not just disabled by
default: there is no W+jets SS file in NanoMVA/Samples/ yet (only Data_OS_merged,
Data_SS_merged, Signal_merged, WtoLNu2Jets_merged exist) -- see
draw_dsmass_os_ss.py's module docstring for the same constraint. Re-add it once a
W+jets SS production exists.

Signal is selected OS (pi opposite-sign to mu). Training the classifier against an
SS background means the two training classes differ structurally in charge
configuration, not just in physics content. Any feature built from the Ds/pi
4-vector relative to the muon (dR_pi_mu, dR_Ds_mu, Dsmu_mass, Dsmu_pt/eta, ...) could
encode "OS vs SS" rather than "signal-like vs background-like" -- and a classifier
leaning on that would score real OS background as signal-like too, since real OS
background also sits on the OS side of the artifact.

This script uses NO signal. It compares two background-only samples that share the
b-veto and differ only in charge configuration, and reports per feature a
standardized mean difference (Cohen's d) and a KS statistic. Large values flag
features where the charge selection, not physics, drives the distribution.

Measured 2026-08-10 on the pre-migration data (43 features, preEE-only, old
tuplizer): max |d| = 0.10 (Ds_pt), max KS = 0.041. Re-measure on this run --
different production, different feature set (41), both eras now.

The feature list comes from common.py, so this can never drift from what is trained.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_os_vs_ss_features.py
"""

import os
import sys

import numpy as np
import uproot
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HALFWIDTH = 1.96835, 0.04
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]

if len(sys.argv) > 1:
    sys.exit("usage: python3 check_os_vs_ss_features.py  (data only -- "
             "no W+jets SS file exists in NanoMVA/Samples/ yet)")

FEATURES = common.FEATURE_NAMES


def load(name, ss):
    br = list(dict.fromkeys(FEATURES + common.SELECTION + ["Ds_mass"]))
    a = uproot.open(f"{common.SAMPLES_DIR}/{name}.root")["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= a[common.B_VETO] == 0
    # OS data is SR-blinded, so the same sideband requirement goes on both sides to
    # compare like with like.
    m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    out = {k: v[m] for k, v in a.items()}
    for feat, fill in common.SENTINEL_FIX.items():
        out[feat][out[feat] < -1] = fill
    return out


def load_group(ss):
    d = load("Data_SS_merged" if ss else "Data_OS_merged", ss)
    return d, np.ones(len(d["Ds_mass"]))


print("Loading OS and SS (data, both eras, sideband-only)...")
os_, os_w = load_group(ss=False)
ss_, ss_w = load_group(ss=True)
print(f"OS: {len(os_['Ds_mass']):,} | SS: {len(ss_['Ds_mass']):,}")

print(f"\n{'feature':32s} {'mean_OS':>10s} {'mean_SS':>10s} {'cohen_d':>9s} {'KS':>7s}")
rows = []
for f in FEATURES:
    x, y = os_[f], ss_[f]
    # a few ip-significance rows are inf/nan; common.preprocess drops them the same way
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    mx, my, sx, sy = np.mean(x), np.mean(y), np.std(x), np.std(y)
    pooled = np.sqrt((sx**2 + sy**2) / 2)
    d = (my - mx) / pooled if pooled > 0 else np.nan
    rows.append((f, mx, my, d, ks_2samp(x, y, method="asymp").statistic))

rows.sort(key=lambda r: -abs(r[3]))
for f, mx, my, d, ks in rows:
    flag = "  <-- LARGE" if abs(d) > 0.2 or ks > 0.1 else ""
    print(f"{f:32s} {mx:10.4f} {my:10.4f} {d:9.3f} {ks:7.3f}{flag}")

# ---------- plot: features sorted most-divergent first, plus Ds_mass ----------
ordered = [r[0] for r in rows]
dmap = {r[0]: r[3] for r in rows}
panels = ordered + ["Ds_mass"]
N_COLS = 5
N_ROWS = int(np.ceil(len(panels) / N_COLS))
fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(5.5 * N_COLS, 4.6 * N_ROWS))
axes_flat = axes.flatten()

for i, f in enumerate(panels):
    ax = axes_flat[i]
    x, y = os_[f], ss_[f]
    fx, fy = np.isfinite(x), np.isfinite(y)
    if f == "Ds_mass":
        bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)
    elif f.startswith("n_b_jets") or f.startswith("num_jets"):
        bins = np.arange(-0.5, 12.5, 1.0)
    else:
        lo, hi = np.percentile(np.concatenate([x[fx], y[fy]]), [0.5, 99.5])
        bins = np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)
    ax.hist(x[fx], bins=bins, weights=os_w[fx], density=True, histtype="step", lw=2,
            color="steelblue", label="OS (control region)")
    ax.hist(y[fy], bins=bins, weights=ss_w[fy], density=True, histtype="step", lw=2,
            ls="--", color="firebrick", label="SS (control region)")
    if f == "Ds_mass":
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.set_title("Ds_mass (fit variable, not an input)", fontsize=12)
    else:
        ax.set_title(f, fontsize=12)
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=9)
    if i == 0:
        hep.cms.label(ax=ax, rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", data=True,
                      text="Preliminary", fontsize=10)

for j in range(len(panels), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.tight_layout()
out = common.fig_path("features", "os_vs_ss_features.png")
plt.savefig(out, dpi=100)
print(f"\nSaved {out}")

# ---------- focused figure: an explicit fixed feature list ----------
# dR_pi_mu/dR_Km_mu/dR_phi_mu/dR_Kp_mu move together with dR_Ds_mu (same geometric
# quantity from different legs of the same decay, see printed table) -- keep only
# dR_Ds_mu as the representative, no auto-fill from the ranking.
top = ["Ds_pt", "Ds_IsoR04_sumPhotonEt", "Ds_IsoR04_sumNeutralHadronEt", "dira",
       "Ds_nearestjet_dR", "Dsmu_mass_corr", "dR_Ds_mu"]
N_COLS2 = 4
N_ROWS2 = int(np.ceil(len(top) / N_COLS2))
fig2, axes2 = plt.subplots(N_ROWS2, N_COLS2, figsize=(5.5 * N_COLS2, 4.6 * N_ROWS2))
axes2_flat = axes2.flatten()

for i, f in enumerate(top):
    ax = axes2_flat[i]
    x, y = os_[f], ss_[f]
    fx, fy = np.isfinite(x), np.isfinite(y)
    if f.startswith("n_b_jets") or f.startswith("num_jets"):
        bins = np.arange(-0.5, 12.5, 1.0)
    else:
        lo, hi = np.percentile(np.concatenate([x[fx], y[fy]]), [0.5, 99.5])
        bins = np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)
    ax.hist(x[fx], bins=bins, weights=os_w[fx], density=True, histtype="step", lw=2,
            color="steelblue", label="OS (control region)")
    ax.hist(y[fy], bins=bins, weights=ss_w[fy], density=True, histtype="step", lw=2,
            ls="--", color="firebrick", label="SS (control region)")
    ax.set_title(f, fontsize=12)
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=9)
    if i == 0:
        hep.cms.label(ax=ax, rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", data=True,
                      text="Preliminary", fontsize=10)

for j in range(len(top), len(axes2_flat)):
    axes2_flat[j].set_visible(False)

plt.tight_layout()
out2 = common.fig_path("features", "os_vs_ss_features_selected.png")
plt.savefig(out2, dpi=110)
print(f"Saved {out2}")

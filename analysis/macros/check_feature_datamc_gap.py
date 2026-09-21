"""Rank the 41 DNN input features (common.FEATURE_NAMES) by shape disagreement
between W+jets MC and OS data, in the Ds mass sideband (OS data is SR-blind by
construction; W+jets is cut the same way so the comparison is on equal footing --
same convention as draw_features.py). This is the input-level picture behind the
tail-cut study in check_feature_tail.py, but does not need a trained model or its
Prediction/ dir -- it reads NanoMVA/Samples/ directly.

Metric: total variation distance (TVD) between the two normalized shapes,
TVD = 0.5 * sum_bins |p_data - p_wjets|, bounded in [0, 1]. 0 = identical shape,
1 = disjoint support. W+jets is weighted (absolute sigma*L/sum(w)); OS data is
unweighted counts. Sentinel rows (-777 "no jet", and ~-10 "tagger unavailable" for
Ds_nearestjet_bscore) are excluded from the shape comparison per feature, same
convention as draw_features.py -- only common.SENTINEL_FIX features carry one.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_feature_datamc_gap.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HALF = 1.96835, 0.040
FEATURES = common.FEATURE_NAMES
SEL = [*common.SELECTION]


def ossel(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
             & (a[common.B_VETO] == 0))


def sssel(a):
    return ((a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
             & (a[common.B_VETO] == 0))


def load(name, extra, sel=ossel):
    a = uproot.open(f"{common.SAMPLES_DIR}/{name}.root")["Events"].arrays(
        FEATURES + SEL + ["Ds_mass"] + extra, library="np")
    m = sel(a) & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALF)   # sideband, both samples
    return {k: v[m] for k, v in a.items()}


print("Loading (OS + b-veto + Ds mass sideband)...")
d = load("Data_OS_merged", [])
w = load("WtoLNu2Jets_merged", ["weight"])
w_w = w["weight"]
print(f"  OS data {len(d['Ds_mass']):,} rows | W+jets {len(w['Ds_mass']):,} rows, "
      f"{w_w.sum():,.0f} weighted")

print("Loading SS data (b-veto + Ds mass sideband)...")
ss = load("Data_SS_merged", [], sel=sssel)
print(f"  SS data {len(ss['Ds_mass']):,} rows")


def bins_for(v, pooled):
    if v.startswith("n_b_jets") or v.startswith("num_jets"):
        return np.arange(-0.5, min(pooled.max(), 12) + 1.5, 1.0)
    lo, hi = np.percentile(pooled, [0.5, 99.5])
    return np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)


rows = []
for v in FEATURES:
    dv, wv = d[v], w[v]
    thr = -1.0 if v in common.SENTINEL_FIX else -700.0
    dok, wok = dv > thr, wv > thr
    pooled = np.concatenate([dv[dok], wv[wok]])
    bins = bins_for(v, pooled)
    width = np.diff(bins)

    hd, _ = np.histogram(dv[dok], bins=bins, density=True)
    hw, _ = np.histogram(wv[wok], bins=bins, weights=w_w[wok], density=True)
    tvd = 0.5 * np.sum(np.abs(hd - hw) * width)

    mu_d, mu_w = dv[dok].mean(), np.average(wv[wok], weights=w_w[wok])
    var_w = np.average((wv[wok] - mu_w) ** 2, weights=w_w[wok])
    sig_w = np.sqrt(var_w) if var_w > 0 else 0.0   # guards a near-zero-variance feature
    shift = (mu_d - mu_w) / sig_w if sig_w > 0 else np.nan

    n_eff_w = w_w[wok].sum() ** 2 / (w_w[wok] ** 2).sum()
    rows.append((v, tvd, shift, dok.mean(), n_eff_w))

rows.sort(key=lambda r: -r[1])

print(f"\n=== ranked by TVD (0=identical shape, 1=disjoint), OS sideband, "
      f"{len(FEATURES)} DNN features ===")
print(f"{'feature':>32s} {'TVD':>7s} {'mean shift (sigma_wj)':>22s} "
      f"{'data frac. valid':>17s} {'W+jets N_eff':>13s}")
for v, tvd, shift, fvalid, neff in rows:
    print(f"{v:>32s} {tvd:7.3f} {shift:22.2f} {100*fvalid:16.1f}% {neff:13,.0f}")

# ---------- bar chart, all features ranked ----------
fig, ax = plt.subplots(figsize=(10, 12))
names = [r[0] for r in rows]
tvds = [r[1] for r in rows]
ypos = np.arange(len(names))[::-1]
ax.barh(ypos, tvds, color="steelblue")
ax.set_yticks(ypos)
ax.set_yticklabels(names, fontsize=9)
ax.set_xlabel("shape difference, OS data vs W+jets MC (sideband)")
ax.set_title("DNN feature shape mismatch, ranked", fontsize=14)
hep.cms.label(ax=ax, rlabel="2022+2022EE (13.6 TeV)", data=True, text="Preliminary",
              fontsize=10)
plt.tight_layout()
out1 = common.fig_path("features", "feature_datamc_tvd_ranked.png", key=None)
plt.savefig(out1, dpi=120)
plt.close(fig)
print(f"\nSaved {out1}")

def draw_shape_grid(features, out_name, n_cols):
    n_rows = int(np.ceil(len(features) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 4.3 * n_rows))
    axes_flat = np.atleast_1d(axes).flatten()
    for i, v in enumerate(features):
        ax = axes_flat[i]
        dv, wv = d[v], w[v]
        thr = -1.0 if v in common.SENTINEL_FIX else -700.0
        dok, wok = dv > thr, wv > thr
        bins = bins_for(v, np.concatenate([dv[dok], wv[wok]]))
        ax.hist(wv[wok], bins=bins, weights=w_w[wok], density=True, alpha=0.5,
                color="seagreen", label="W+jets MC")
        ax.hist(dv[dok], bins=bins, density=True, histtype="step", lw=2, color="black",
                label="OS data")
        ax.set_title(v, fontsize=11)
        ax.set_ylabel("a.u.")
        ax.legend(fontsize=8)
    for j in range(len(features), len(axes_flat)):
        axes_flat[j].set_visible(False)
    plt.tight_layout()
    out = common.fig_path("features", out_name, key=None)
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"Saved {out}")


# ---------- shape overlay for the top offenders ----------
draw_shape_grid([r[0] for r in rows[:8]], "feature_datamc_tvd_top8_shapes.png", 4)

# ---------- shape overlay, jet- and muon-related features ----------
JET_FEATURES = ["num_jets_pt30_corr", "jet_leadpt_corr", "Ds_nearestjet_bscore",
                "Ds_nearestjet_dR", "n_b_jets_loose_dR04", "n_b_jets_medium_dR04"]
MU_FEATURES = ["mu_pt_corr", "mu_IsoR04_RelIso"]
JET_MU = JET_FEATURES + MU_FEATURES
draw_shape_grid(JET_MU, "feature_datamc_jet_mu_shapes.png", 4)


def draw_os_ss_grid(features, out_name, n_cols):
    n_rows = int(np.ceil(len(features) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 4.3 * n_rows))
    axes_flat = np.atleast_1d(axes).flatten()
    for i, v in enumerate(features):
        ax = axes_flat[i]
        dv, sv = d[v], ss[v]
        thr = -1.0 if v in common.SENTINEL_FIX else -700.0
        dok, sok = dv > thr, sv > thr
        bins = bins_for(v, np.concatenate([dv[dok], sv[sok]]))
        ax.hist(dv[dok], bins=bins, density=True, histtype="step", lw=2,
                color="black", label="OS data")
        ax.hist(sv[sok], bins=bins, density=True, histtype="step", lw=2, ls="--",
                color="steelblue", label="SS data")
        ax.set_title(v, fontsize=11)
        ax.set_ylabel("a.u.")
        ax.legend(fontsize=8)
    for j in range(len(features), len(axes_flat)):
        axes_flat[j].set_visible(False)
    plt.tight_layout()
    out = common.fig_path("features", out_name, key=None)
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"Saved {out}")


draw_os_ss_grid(JET_MU, "feature_os_ss_jet_mu_shapes.png", 4)

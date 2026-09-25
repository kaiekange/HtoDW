"""All 43 MVA input features: signal vs W+jets MC vs OS data vs SS data.

This is the pre-MVA input-level picture -- it reads Samples/ directly and needs no
DNN score, so it stays valid across retrainings. Two questions it answers at a
glance:

  1. **Separation** -- where does signal actually differ from background, and is the
     feature doing anything at all?
  2. **Data/MC modelling** -- does W+jets MC reproduce real data? A feature the MC
     gets wrong is a liability in a classifier no matter how discriminating it looks,
     and W+jets is known to be an incomplete background model
     (project_wjets_background_incomplete), so gaps here are expected and worth
     seeing rather than assuming.

The feature list comes from common.py, so this figure can never show a different set
than the one actually trained on.

Values are plotted RAW (no log transform / sentinel remap) because the point is
physical interpretation; the ranges are percentile-clipped and sentinel rows
(-777 etc.) are excluded from the plotted range, with their rate reported in the
title. common.preprocess is deliberately NOT applied.

All three BACKGROUND samples are restricted to the Ds mass sideband so they are on
equal footing: OS data is SR-blinded by construction, and cutting the other two the
same way removes an artificial difference. Signal is left unrestricted -- it is real
Ds and lives in the SR by definition. The features are designed to be uncorrelated
with M(KKpi) (that is why decay kinematics are excluded), so this does not bias the
signal/background comparison.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_features.py
"""

import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HALFWIDTH = 1.96835, 0.04

FEATURES = common.FEATURE_NAMES


def load(name, *, era=None, ss=False, sig=False, sideband=False):
    br = list(dict.fromkeys(FEATURES + common.SELECTION + ["Ds_mass", "Lumi_reweight"]
                            + (["era"] if era is not None else [])
                            + (["Ds_match", "mu_match"] if sig else [])))
    a = uproot.open(f"{common.SAMPLES_DIR}/{name}.root")["Events"].arrays(br, library="np")
    # OS/SS is a property of the FILE (producer's sameSignPi), so it is an ASSERTION now,
    # not a selection; the muon cuts come from common.selection_mask. The old
    # num_os/ss_mu_tight ==1/==0 pair conflated the two and additionally imposed an
    # uncalibrated soft-in-jet-muon tag. See common.check_charge_config / selection_mask.
    common.check_charge_config(a, ss, f"load({name!r})")
    m = common.selection_mask(a)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if sideband:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


print("Loading...")
# Signal: both eras, as before -- but now WEIGHTED. This used to concatenate the four
# per-era signal files at equal weight, silently dropping the postEE factor of ~3.34.
sig_parts = [load(n, sig=True) for n in common.SIGNALS]
sig = {k: np.concatenate([p[k] for p in sig_parts]) for k in sig_parts[0]}
sig_w = sig["Lumi_reweight"]

# W+jets and data stay preEE-only (era=1), which is what this macro always plotted --
# WJETS_WEIGHTS only ever held preEE keys and the data literal was Data_2022. Dropping
# `era=1` is now all it takes to include postEE, but that changes every panel, so it is
# left as a deliberate choice rather than a side effect of the merge.
wj = load("WtoLNu2Jets", era=1, sideband=True)
wj_w = wj["Lumi_reweight"]

os_ = load("Data", era=1, sideband=True)
ss_ = load("Data_SS", era=1, ss=True, sideband=True)

print(f"signal {len(sig['Ds_mass']):,} | W+jets {len(wj['Ds_mass']):,} "
      f"| OS data {len(os_['Ds_mass']):,} | SS data {len(ss_['Ds_mass']):,}")

N_COLS = 5
N_ROWS = int(np.ceil(len(FEATURES) / N_COLS))
fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(5.5 * N_COLS, 4.6 * N_ROWS))
axes_flat = axes.flatten()

for i, v in enumerate(FEATURES):
    ax = axes_flat[i]
    s, b, o, x = sig[v], wj[v], os_[v], ss_[v]
    # Sentinel exclusion must be per-feature. Only the three SENTINEL_FIX variables
    # carry a below-range filler (-777 "no jet", and ~-10 "tagger unavailable" for
    # the bscore); for everything else the only sentinel is the -777 null. Applying
    # a blanket z > -1 would silently chop the negative half off eta variables.
    thr = -1.0 if v in common.SENTINEL_FIX else -700.0
    ok = lambda z: z > thr
    pool = np.concatenate([s[ok(s)], b[ok(b)], o[ok(o)], x[ok(x)]])
    if v.startswith("n_b_jets") or v.startswith("num_jets"):
        bins = np.arange(-0.5, min(pool.max(), 12) + 1.5, 1.0)
    else:
        lo, hi = np.percentile(pool, [0.5, 99.5])
        bins = np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)

    ax.hist(b, bins=bins, weights=wj_w, density=True, alpha=0.45, color="seagreen",
            label="W+jets MC")
    ax.hist(o, bins=bins, density=True, histtype="step", lw=2, ls="--", color="black",
            label="OS data")
    ax.hist(x, bins=bins, density=True, histtype="step", lw=2, color="steelblue",
            label="SS data")
    ax.hist(s, bins=bins, weights=sig_w, density=True, alpha=0.4, color="firebrick",
            label="signal")

    frac_sentinel = 100 * np.mean(~ok(b))
    title = v if frac_sentinel < 0.01 else f"{v}\n({frac_sentinel:.1f}% sentinel, off-range)"
    ax.set_title(title, fontsize=12)
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=9)
    if i == 0:
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=10)

for j in range(len(FEATURES), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.tight_layout()
out = common.fig_path("features", "features_sig_wjets_data.png")
plt.savefig(out, dpi=100)
print(f"Saved {out}")

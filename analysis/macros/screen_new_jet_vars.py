"""Screen the 18 jet/b-tag variables added to the tuplizer on 2026-07-23 before
committing them to the MVA feature set.

Re-run 2026-08-10 for the planned switch to SS data as the MVA training background
(instead of W+jets MC): the tight b-veto (n_b_jets_tight_dR04==0) is now applied
before every comparison -- it already strips most tt-like events, so separation
seen without it (the 2026-08-03 version of this script) may not survive. The "Data"
reference is now SS (Data_2022_SS) instead of OS SR-blinded, since SS is the
candidate training background and needs no blind window.

Three things decide whether a variable earns its slot:
  1. separation      -- does signal differ from background, AFTER the b-veto?
  2. data/MC modelling -- does the background MC reproduce SS data? A variable the MC
     gets wrong is a liability in a classifier no matter how discriminating it looks.
  3. mass sculpting  -- correlation with Ds_mass / phi_mass, measured on SS data
     itself (the actual candidate training background), not on W+jets MC. The
     existing 37-feature set was chosen to keep this near zero so the MVA does not
     carve a fake peak into the fit variable.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 screen_new_jet_vars.py
"""

import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

# The b-veto and the lumi weights come from common.py so this screening always runs
# against the accepted baseline. NEW_VARS below stays local on purpose: its job is to
# evaluate variables that are NOT yet in the feature set, so it must not be coupled
# to common.FEATURE_NAMES -- some of the 18 have since been accepted, some rejected,
# and the point of keeping the full list is to preserve the comparison that decided.
SAMPLES_DIR = common.SAMPLES_DIR
M_DS = 1.96835
SR_HALFWIDTH = 0.04

NEW_VARS = [
    "Ds_nearestjet_bscore", "Ds_nearestjet_dR",
    "jet_ht", "jet_ht_pt20", "jet_leadpt", "jet_maxbscore",
    "n_b_jets_loose_dR03", "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR03", "n_b_jets_medium_dR04",
    "n_b_jets_tight_dR03", "n_b_jets_tight_dR04",
    "n_b_jets_loose_pt20_dR03", "n_b_jets_loose_pt20_dR04",
    "n_b_jets_medium_pt20_dR03", "n_b_jets_medium_pt20_dR04",
    "n_b_jets_tight_pt20_dR03", "n_b_jets_tight_pt20_dR04",
]


SEL = common.SELECTION
EXTRA = ["Ds_mass", "phi_mass"]


def one_tight_mu(a, ss=False):
    if ss:
        return (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0) & (a[common.B_VETO] == 0)
    return (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0) & (a[common.B_VETO] == 0)


def load(name, truth_filter=False, ss=False, era=None):
    br = (NEW_VARS + SEL + EXTRA + ["Lumi_reweight"] + (["era"] if era else [])
          + (["Ds_match", "mu_match"] if truth_filter else []))
    a = uproot.open(f"{SAMPLES_DIR}/{name}.root")["Events"].arrays(
        list(dict.fromkeys(br)), library="np")
    m = one_tight_mu(a, ss=ss)
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


print("Loading...")
# preEE only (era=1), as before: the old WJETS_WEIGHTS had only preEE keys and the
# data/signal literals carried _2022. Per-event Lumi_reweight replaces the constants.
sig = load("HtoDsW", truth_filter=True, era=1)
data = load("Data_SS", ss=True, era=1)
bkg = load("WtoLNu2Jets", era=1)
bkg_w = bkg["Lumi_reweight"]

print(f"signal {len(sig['Ds_mass']):,} | background {len(bkg['Ds_mass']):,} | data {len(data['Ds_mass']):,}")

# W+jets MC restricted to the SR sideband, purely for the plotted reference curve
# (SS data itself is not blinded, so it needs no such restriction)
bkg_sb = np.abs(bkg["Ds_mass"] - M_DS) >= SR_HALFWIDTH
print(f"W+jets MC in sideband (plot reference only): {bkg_sb.sum():,} "
      f"({100 * bkg_sb.mean():.1f}%)")

# ---------- sentinel audit ----------
print("\n=== sentinel fractions (-777 = not computable, ~-10 = tagger unavailable) ===")
for v in NEW_VARS:
    s777 = 100 * np.mean(bkg[v] == -777)
    sneg = 100 * np.mean((bkg[v] < -1) & (bkg[v] != -777))
    if s777 > 0.01 or sneg > 0.01:
        print(f"  {v:32s} bkg: -777 {s777:5.2f}%   other-negative {sneg:5.2f}%")

# ---------- sculpting check ----------
# Pearson r against the two fit variables, measured on SS data -- the actual
# candidate training background, not W+jets MC. |r| well under ~0.1 is the
# working target.
print("\n=== sculpting risk: corr(feature, fit variable) on SS data ===")
print(f"{'feature':32s} {'r(Ds_mass)':>11s} {'r(phi_mass)':>12s}")
sculpt = {}
for v in NEW_VARS:
    ok = data[v] > -1          # drop sentinels, they would fake a correlation
    if ok.sum() < 100:
        continue
    r_ds = np.corrcoef(data[v][ok], data["Ds_mass"][ok])[0, 1]
    r_phi = np.corrcoef(data[v][ok], data["phi_mass"][ok])[0, 1]
    sculpt[v] = (r_ds, r_phi)
    flag = "  <-- check" if max(abs(r_ds), abs(r_phi)) > 0.05 else ""
    print(f"{v:32s} {r_ds:11.4f} {r_phi:12.4f}{flag}")

# ---------- comparison plot ----------
print("\nDrawing...")
N_COLS, N_ROWS = 5, 4
fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(26, 19))
axes_flat = axes.flatten()

for i, v in enumerate(NEW_VARS):
    ax = axes_flat[i]
    s, b, d = sig[v], bkg[v], data[v]
    bw = bkg_w
    # sentinels are excluded from the plotted range but their rate is in the title
    ok_all = np.concatenate([s[s > -1], b[b > -1], d[d > -1]])
    if v.startswith("n_b_jets"):
        bins = np.arange(-0.5, 5.5, 1.0)
    else:
        lo, hi = np.percentile(ok_all, [0.5, 99.5])
        bins = np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)

    ax.hist(b[bkg_sb], bins=bins, weights=bw[bkg_sb], density=True, histtype="step",
            lw=2, color="steelblue", label="W+jets (sideband)")
    ax.hist(d, bins=bins, density=True, histtype="step", lw=2, ls="--",
            color="black", label="Data SS")
    ax.hist(s, bins=bins, density=True, alpha=0.45, color="firebrick", label="HtoDsW")
    r = sculpt.get(v, (np.nan, np.nan))[0]
    ax.set_title(f"{v}\n" + r"$r(M_{D_s})$=" + f"{r:+.3f}", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=11)
    if i == 0:
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=12)

for j in range(len(NEW_VARS), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.tight_layout()
out = common.fig_path("features", "new_jet_vars_screening.png")
plt.savefig(out, dpi=130)
print(f"Saved {out}")

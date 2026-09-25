"""Dsst_mass (M(Ds+gamma)) for the two signal channels, stacked by gamma truth category.

Both channels use the same event selection: Ds_match & mu_match (Dsmu truth-matched,
see project memory) & has_gamma (Dsst_mass is null otherwise). Rows are per-candidate
(Ds x gamma combination), not deduped to one row per event -- same convention as
draw_wp_mass.py.

DsW has no true gen-level gamma at all (Ds*W is a different channel), so gamma_match
is false by construction there -- the DsW panel is a trivial one-component "stack"
(100% fake gamma), included for a like-for-like comparison against DsstW.

DsstW splits Dsst_match (== Ds_match & gamma_match, see H_phipimu.cc) from
has_gamma & !gamma_match ("fake"/combinatorial gamma).

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_signal_dsst_mass_truth.py
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

SIG_FILE = f"{common.SAMPLES_DIR}/Signal_merged.root"
BINS = np.arange(1.95, 2.26 + 1e-9, 0.005)
BC = 0.5 * (BINS[1:] + BINS[:-1])
WIDTH = BINS[1] - BINS[0]

CHANNELS = [(0, "DsW"), (1, "DsstW")]
COLOR_TRUE, COLOR_FAKE = "firebrick", "seagreen"
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]

a = uproot.open(SIG_FILE)["Events"].arrays(
    ["channel", "Ds_match", "mu_match", "has_gamma", "gamma_match", "Dsst_mass", "weight"],
    library="np")
base = a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) & a["has_gamma"].astype(bool)

fig, axes = plt.subplots(1, 2, figsize=(17.0, 7.0), gridspec_kw={"wspace": 0.24})

for col, (ch, label) in enumerate(CHANNELS):
    ax = axes[col]
    sel = base & (a["channel"] == ch)
    true_m = sel & a["gamma_match"].astype(bool)
    fake_m = sel & ~a["gamma_match"].astype(bool)

    h_true, _ = np.histogram(a["Dsst_mass"][true_m], bins=BINS, weights=a["weight"][true_m])
    h_fake, _ = np.histogram(a["Dsst_mass"][fake_m], bins=BINS, weights=a["weight"][fake_m])
    e_tot = np.sqrt(np.histogram(a["Dsst_mass"][sel], bins=BINS, weights=a["weight"][sel] ** 2)[0])
    h_tot = h_fake + h_true
    n_true, n_fake = true_m.sum(), fake_m.sum()

    ax.bar(BC, h_fake, width=WIDTH, color=COLOR_FAKE, alpha=0.55, edgecolor="none",
           label=f"fake/combinatorial $\\gamma$  N={n_fake:,}  yield={h_fake.sum():.4f}")
    ax.bar(BC, h_true, width=WIDTH, bottom=h_fake, color=COLOR_TRUE, alpha=0.55,
           edgecolor="none", label=f"true gen $\\gamma$  N={n_true:,}  yield={h_true.sum():.4f}")
    ax.bar(BC, 2 * e_tot, width=WIDTH, bottom=h_tot - e_tot, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3, label="MC stat")

    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_ylim(0, h_tot.max() * 1.55)
    ax.set_xlabel(r"$M(D_s^\pm\gamma)$ [GeV]")
    if col == 0:
        ax.set_ylabel(f"Events / {WIDTH*1000:.0f} MeV")
    # panel description goes in the legend title, not ax.set_title(), so it never
    # fights hep.cms.label() for the same row above the axes
    ax.legend(fontsize=8.5, loc="upper right",
              title=f"{label} signal, Dsmu truth-matched  N={sel.sum():,}", title_fontsize=8.5)
    if col == 0:
        hep.cms.label(ax=ax, data=False, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

out = common.fig_path("preselection", "signal_dsst_mass_gamma_truth.png", key=None)
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

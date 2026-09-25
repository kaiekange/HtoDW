"""n_b_jets multiplicity for signal / W+jets MC / data, all 12 variants
({loose,medium,tight} x {default,pt20} x {dR03,dR04}), to pick the b-jet veto
working point.

The veto's job is not signal/background discrimination -- it is removing the
process the background model does not contain (tt, real b jets), which is why
data is the meaningful reference here and not the W+jets MC. Data is the OS
SR-blinded sample, so the same sideband requirement is applied to the background
MC to compare like with like.

RETIRED 2026-09-04 -- THIS CANNOT RUN AGAINST THE CURRENT Samples/.

The b-jet veto (n_b_jets_tight_dR04 == 0) is now applied by make_samples.py BEFORE the
dedup and baked into every file, so every one of the 12 variants this macro plots is
truncated to its veto-passing part: the tight/dR04 one is a single bar at zero by
construction, and the other 11 are the heavily correlated remainder. The figure would
still be drawn, and it would be wrong in a way that looks plausible -- which is why this
exits instead.

The working point it was written to choose is settled (tight, dR04); this file is kept
as the record of that decision. To run it again you need Samples/ built without the
b_veto cut -- drop "b_veto" from the `cuts` tuple in make_samples.py's SAMPLES table and
rebuild into a separate directory.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_nbjets_veto.py
"""

import sys

import numpy as np
import uproot

sys.exit(__doc__[__doc__.index("RETIRED"):__doc__.index("Usage:")].strip())
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

SAMPLES_DIR = common.SAMPLES_DIR
M_DS, SR_HALFWIDTH = 1.96835, 0.04

WPS = ["loose", "medium", "tight"]
PTS = ["", "_pt20"]
DRS = ["dR03", "dR04"]
VARS = [f"n_b_jets_{wp}{pt}_{dr}" for wp in WPS for pt in PTS for dr in DRS]

WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J_2022": 2.3293047885,
    "WtoLNu2Jets_1J_2022": 0.5382864853,
    "WtoLNu2Jets_2J_2022": 0.2748710631,
}
SEL = ["num_os_mu_tight", "num_ss_mu_tight"]


def load(name, truth_filter=False):
    br = VARS + SEL + ["Ds_mass"] + (["Ds_match", "mu_match"] if truth_filter else [])
    a = uproot.open(f"{SAMPLES_DIR}/{name}.root")["Events"].arrays(br, library="np")
    # bveto=False: this macro PICKS the b-veto working point, so applying it here would
    # collapse the very distribution being plotted. Muon cuts still come from common.
    m = common.selection_mask(a, bveto=False)
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


print("Loading...")
sig = load("HtoDsW_2022", truth_filter=True)
data = load("Data_2022")
parts = [load(n) for n in WJETS_WEIGHTS]
bkg = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
bkg_w = np.concatenate([np.full(len(p["Ds_mass"]), w)
                        for p, w in zip(parts, WJETS_WEIGHTS.values())])
sb = np.abs(bkg["Ds_mass"] - M_DS) >= SR_HALFWIDTH      # match the blinded data
bkg = {k: v[sb] for k, v in bkg.items()}
bkg_w = bkg_w[sb]
print(f"signal {len(sig['Ds_mass']):,} | W+jets(sideband) {len(bkg['Ds_mass']):,} "
      f"| data(sideband) {len(data['Ds_mass']):,}")

bins = np.arange(-0.5, 5.5, 1.0)
c = 0.5 * (bins[1:] + bins[:-1])

# One figure per dR working point; within each, rows = b-tag WP, columns = jet pt
# threshold. The two files are drawn identically so they can be flipped between.
print(f"\n{'variable':30s} {'eps_sig(==0)':>13s} {'eps_Wjets':>10s} {'eps_data':>9s} {'FoM':>8s}")
for dr in DRS:
    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    for row, wp in enumerate(WPS):
        for col, pt in enumerate(PTS):
            v = f"n_b_jets_{wp}{pt}_{dr}"
            ax = axes[row][col]
            s, b, d = sig[v], bkg[v], data[v]

            hs, _ = np.histogram(s, bins=bins); hs = hs / hs.sum()
            hb, _ = np.histogram(b, bins=bins, weights=bkg_w); hb = hb / hb.sum()
            hd, _ = np.histogram(d, bins=bins); hd = hd / hd.sum()

            es, eb, ed = (s == 0).mean(), np.average(b == 0, weights=bkg_w), (d == 0).mean()
            fom = es / np.sqrt(ed)
            print(f"{v:30s} {100*es:12.2f}% {100*eb:9.2f}% {100*ed:8.2f}% {fom:8.4f}")

            hep.histplot(hs, bins, ax=ax, histtype="fill", alpha=0.45, color="firebrick",
                         label=r"H$\to$D$_s$W")
            hep.histplot(hb, bins, ax=ax, histtype="step", lw=2, color="steelblue",
                         label="W+jets")
            hep.histplot(hd, bins, ax=ax, yerr=False, histtype="errorbar", color="black",
                         marker="o", markersize=5, label="Data")
            ax.set_xlabel(v.replace("_", " "))
            ax.set_ylabel("Fraction of events")
            ax.set_ylim(0, 1.25)
            ax.set_xticks(range(0, 5))
            ax.legend(fontsize=13, loc="upper right")
            hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=13)

    plt.tight_layout()
    out = common.fig_path("preselection", f"n_b_jets_veto_{dr}.png")
    plt.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  -> saved {out}")

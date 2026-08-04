"""Compare the b-jet-veto counts (n_b_jets_{loose,medium,tight}_dR{03,04}, jets
outside the Ds cone) between the H->DsW signal and the WtoLNu2Jets background,
using the deduped (highest-Ds_pt) Samples/ tuples.

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_jet_vars.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"

WPS = ["loose", "medium", "tight"]
DRS = ["dR03", "dR04"]
BRANCHES = [f"n_b_jets_{wp}_{dr}" for wp in WPS for dr in DRS] + ["event_match"]

sig = uproot.open(f"{SAMPLES_DIR}/HtoDsW.root")["Events"].arrays(BRANCHES, library="np")
bkg = uproot.open(f"{SAMPLES_DIR}/WtoLNu2Jets.root")["Events"].arrays(BRANCHES, library="np")

# Signal: require the truth-matched Ds+mu candidate. Background has doTruth=False
# (event_match is meaningless there), so no truth filter is applied.
sig_sel = sig["event_match"].astype(bool)
print(f"signal: {sig_sel.sum()}/{len(sig_sel)} truth-matched events kept")
print(f"background (WtoLNu2Jets): {len(bkg['event_match'])} events, no truth filter")

bins = np.arange(-0.5, 5.5, 1.0)
fig, axes = plt.subplots(2, 3, figsize=(24, 14), sharey="row")

for row, dr in enumerate(DRS):
    for col, wp in enumerate(WPS):
        branch = f"n_b_jets_{wp}_{dr}"
        s, b = sig[branch][sig_sel], bkg[branch]
        ax = axes[row][col]
        ax.hist(b, bins=bins, density=True, alpha=0.6, color="steelblue",
                 label=f"W+jets (mean {b.mean():.4f})")
        ax.hist(s, bins=bins, density=True, alpha=0.6, color="firebrick",
                 label=f"Signal (mean {s.mean():.4f})")
        ax.set_xlabel(f"N b-jets ({wp} WP, {dr})")
        ax.set_ylabel("a.u.")
        ax.legend(fontsize=15)
        if row == 0 and col == 0:
            hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=16)

plt.tight_layout()
out = f"{OUTDIR}/n_b_jets_signal_vs_wjets.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

# ---- numeric summary ----
print("\n=== summary: P(>=1 b-jet outside the Ds cone) ===")
for dr in DRS:
    for wp in WPS:
        branch = f"n_b_jets_{wp}_{dr}"
        s, b = sig[branch][sig_sel], bkg[branch]
        print(f"  {wp:6s} {dr}:  Signal mean={s.mean():.4f} P(>=1)={100*(s>=1).sum()/s.size:5.2f}%   "
              f"W+jets mean={b.mean():.4f} P(>=1)={100*(b>=1).sum()/b.size:5.2f}%   "
              f"ratio(W+jets/Signal)={b.mean()/s.mean():.2f}")

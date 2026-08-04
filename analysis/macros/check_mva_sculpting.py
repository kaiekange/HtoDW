"""Does cutting on the MVA score reshape the Ds_mass spectrum of the background?

The 42-feature set was built to avoid exactly this -- no decay-kinematics variable
is used, so a linear correlation with Ds_mass stays near zero. But a near-zero
Pearson r does not rule out the MVA *enriching* a real-Ds component that peaks at
M_Ds, which would bias any sideband-based background estimate. This script looks at
the shape directly instead of trusting the correlation coefficient.

Left panels: background M(KKpi) at increasing score cuts, each area-normalized, so
a growing peak shows up as a shape change rather than just fewer events.
Right panels: the same for data (sidebands only -- the SR is blinded).

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_mva_sculpting.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight", "n_b_jets_tight_dR04"]
CUTS = [0.0, 0.5, 0.9, 0.99]
COLORS = ["black", "steelblue", "darkorange", "firebrick"]

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)
c = 0.5 * (bins[1:] + bins[:-1])


def load(model, name):
    p = f"/eos/home-k/kakang/H2DW/Analysis/MVA/{model}/dsw_wjets_test/Prediction/{name}.root"
    a = uproot.open(p)["Events"].arrays(BR, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a["n_b_jets_tight_dR04"] == 0) & (a["score_dsmu"] > -1))
    return {k: v[m] for k, v in a.items()}


fig, axes = plt.subplots(2, 2, figsize=(20, 15))

for row, model in enumerate(["BDT", "DNN"]):
    parts = [load(model, n) for n in WJETS_WEIGHTS]
    bkg = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    bw = np.concatenate([np.full(len(p["Ds_mass"]), w)
                         for p, w in zip(parts, WJETS_WEIGHTS.values())])
    dat = load(model, "Data_2022")

    print(f"\n=== {model} ===")
    for col, (src, label) in enumerate([("mc", "W+jets MC"), ("data", "Data (SR blinded)")]):
        ax = axes[row][col]
        for cut, color in zip(CUTS, COLORS):
            if src == "mc":
                m = bkg["score_dsmu"] > cut
                v, wt = bkg["Ds_mass"][m], bw[m]
                h, _ = np.histogram(v, bins=bins, weights=wt)
                n_eff = wt.sum()
            else:
                m = dat["score_dsmu"] > cut
                v = dat["Ds_mass"][m]
                h, _ = np.histogram(v, bins=bins)
                n_eff = len(v)
            if n_eff < 50:
                continue
            hep.histplot(h / h.sum(), bins, ax=ax, histtype="step", lw=2, color=color,
                         label=f"score > {cut:g}  (N={n_eff:,.0f})")
            if src == "mc":
                frac = np.average(np.abs(bkg["Ds_mass"][m] - M_DS) < SR_HALFWIDTH, weights=wt)
                print(f"  MC   score>{cut:<5g} frac in SR = {100*frac:5.2f}%")

        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15,
                   label="signal region" if cut == CUTS[-1] else None)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        ax.set_ylabel("Fraction of events / bin")
        ax.set_title(f"{model} -- {label}", fontsize=19)
        ax.legend(fontsize=13)
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=(src == "data"),
                      text="Preliminary", fontsize=13)

plt.tight_layout()
out = f"{OUTDIR}/mva_score_sculpting.png"
plt.savefig(out, dpi=120)
print(f"\nSaved {out}")

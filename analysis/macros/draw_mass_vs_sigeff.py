"""M(KKpi) in data and W+jets MC at a series of DNN working points, labelled by the
signal efficiency they correspond to rather than by the raw score.

All the thresholds here sit above score 0.9 (90% signal efficiency is already at
0.917), i.e. this zooms into the tight end where the analysis will actually live.

Both are area-normalized: the absolute data/MC normalization is known to be off by
roughly a factor 2 and to worsen with the cut (the background model is still missing
tt and friends), so overlaying raw yields would only show that gap again. What this
plot is for is the *shape* -- specifically how the D+- and Ds real-charm peaks grow
out of the combinatorial continuum as the cut tightens, in data and in MC.

Data is OS and SR-blinded, so its Ds region is empty by construction; the D+- peak at
1.8696 sits outside the blind window and is the one that can be compared directly.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_mass_vs_sigeff.py
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

PRED = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets/Prediction"
OUTDIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets/figures"
os.makedirs(OUTDIR, exist_ok=True)

M_DS, M_D, SR_HALFWIDTH = 1.96835, 1.8696, 0.04
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}
SIGNALS = ["HtoDsW", "HtoDsstW"]
EFFS = [0.90, 0.70, 0.50, 0.30, 0.10]

BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      "n_b_jets_tight_dR04", "Ds_match", "mu_match"]

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)


def load(name, truth):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(BR, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a["n_b_jets_tight_dR04"] == 0) & (a["score_dsmu"] > -1))
    if truth:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


sig_score = np.concatenate([load(n, True)["score_dsmu"] for n in SIGNALS])
parts = [load(n, False) for n in WJETS_WEIGHTS]
bkg = {k: np.concatenate([p[k] for p in parts]) for k in BR}
bkg_w = np.concatenate([np.full(len(p["Ds_mass"]), w)
                        for p, w in zip(parts, WJETS_WEIGHTS.values())])
dat = load("Data_2022", False)

fig, axes = plt.subplots(1, len(EFFS), figsize=(6 * len(EFFS), 6.5), sharex=True)

print(f"{'eff_sig':>8s} {'score cut':>11s} {'N_Wjets(wt)':>13s} {'N_data':>9s}")
for ax, eff in zip(axes, EFFS):
    thr = np.quantile(sig_score, 1 - eff)

    mb = bkg["score_dsmu"] > thr
    hb, _ = np.histogram(bkg["Ds_mass"][mb], bins=bins, weights=bkg_w[mb])
    md = dat["score_dsmu"] > thr
    hd, _ = np.histogram(dat["Ds_mass"][md], bins=bins)
    print(f"{100*eff:7.0f}% {thr:11.5f} {bkg_w[mb].sum():13,.0f} {md.sum():9,}")

    hep.histplot(hb / hb.sum(), bins, ax=ax, histtype="step", lw=2, color="steelblue",
                 label=f"W+jets ({bkg_w[mb].sum():,.0f})")
    hep.histplot(hd / hd.sum(), bins, ax=ax, histtype="errorbar", color="black",
                 marker="o", markersize=4, yerr=np.sqrt(hd) / hd.sum(),
                 label=f"Data ({md.sum():,})")

    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    ax.set_title(rf"$\varepsilon_{{sig}}$ = {100*eff:.0f}%   (score > {thr:.4f})", fontsize=17)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.legend(fontsize=12, loc="upper left")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)

axes[0].set_ylabel("Fraction of events / bin")
for ax in axes:
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)

plt.tight_layout()
out = f"{OUTDIR}/mass_vs_sigeff.png"
plt.savefig(out, dpi=120)
print(f"\nSaved {out}")

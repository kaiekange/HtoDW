"""OS vs SS version of the MVA sculpting check, on the DNN.

The point this plot has to settle: SS was adopted as a combinatorial-background
template because at inclusive selection it shows no real-Ds peak. That still holds
without an MVA cut, but the numbers say it stops holding once the score cut is
tightened -- Data SS develops a Ds peak the W+jets SS MC does not predict, presumably
because charge-symmetric real charm (tt-like) populates both signs while W+c does not.

Blinding only ever applied to OS: SS excludes real signal by charge conservation, so
there was never anything there to hide. That asymmetry is what makes this measurable
at all -- the Ds region is simply there to look at in SS, while in OS it is cut out.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_mva_sculpting_ss.py
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

OUTDIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets/figures"
os.makedirs(OUTDIR, exist_ok=True)
PRED = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets/Prediction"
M_DS, M_D, SR_HALFWIDTH = 1.96835, 1.8696, 0.04
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight", "n_b_jets_tight_dR04"]
CUTS = [0.0, 0.9, 0.99]
COLORS = ["black", "darkorange", "firebrick"]

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)


def load(name, ss):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(BR, library="np")
    # os/ss counting is defined relative to the pi charge and is flipped between the
    # two plugins, so the "exactly one tight muon" cut mirrors for the SS samples.
    m = ((a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)) if ss \
        else ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0))
    m &= (a["n_b_jets_tight_dR04"] == 0) & (a["score_dsmu"] > -1)
    return {k: v[m] for k, v in a.items()}


def merged(suffix, ss):
    parts = [load(n + suffix, ss) for n in WJETS_WEIGHTS]
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    w = np.concatenate([np.full(len(p["Ds_mass"]), x)
                        for p, x in zip(parts, WJETS_WEIGHTS.values())])
    return d, w


panels = [
    ("W+jets MC, OS", *merged("", False)),
    ("W+jets MC, SS", *merged("_SS", True)),
]
d = load("Data_2022", False);    panels.append(("Data OS (SR blinded)", d, np.ones(len(d["Ds_mass"]))))
d = load("Data_2022_SS", True);  panels.append(("Data SS", d, np.ones(len(d["Ds_mass"]))))

fig, axes = plt.subplots(2, 2, figsize=(20, 15))
for ax, (title, dat, wt) in zip(axes.flatten(), panels):
    is_data = title.startswith("Data")
    for cut, color in zip(CUTS, COLORS):
        m = dat["score_dsmu"] > cut
        h, _ = np.histogram(dat["Ds_mass"][m], bins=bins, weights=wt[m])
        n = wt[m].sum()
        if n < 50:
            continue
        hep.histplot(h / h.sum(), bins, ax=ax, histtype="step", lw=2, color=color,
                     label=f"score > {cut:g}  (N={n:,.0f})")
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.text(M_D, ax.get_ylim()[1] * 0.97, r"D$^\pm$", color="green",
            ha="center", va="top", fontsize=15)
    ax.text(M_DS, ax.get_ylim()[1] * 0.97, r"D$_s^\pm$", color="dimgrey",
            ha="center", va="top", fontsize=15)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel("Fraction of events / bin")
    ax.set_title(title, fontsize=19)
    ax.legend(fontsize=13, loc="upper right")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=is_data, text="Preliminary", fontsize=13)

plt.tight_layout()
out = f"{OUTDIR}/mva_score_sculpting_os_vs_ss.png"
plt.savefig(out, dpi=120)
print(f"Saved {out}")

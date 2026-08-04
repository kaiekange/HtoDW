"""OS vs SS comparison of the luminosity-reweighted sum of the three WtoLNu2Jets
jet-bin MC samples (0J/1J/2J), absolute scale (no shape normalization) -- same
per-sample weights as draw_data_vs_wjets{,_ss}.py, applied separately to the OS and
SS sums so the two curves are on the same physical scale and can be compared directly
(not just their shapes).

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_wjets_os_vs_ss_weighted.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS = 1.96835

WEIGHTS = {"0J": 2.3293047885, "1J": 0.5382864853, "2J": 0.2748710631}

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)   # 10 MeV bins
c = 0.5 * (bins[1:] + bins[:-1])


# Exactly one tight muon in the event (the best muon itself, no extra of either
# charge). os/ss are relative to the pi charge, which is flipped between the OS/SS
# plugins, so the two selections are each other's mirror -- see H_phipimu(_SS).cc
# and draw_data_vs_wjets{,_ss}.py.
SEL_BRANCHES = ["num_os_mu_tight", "num_ss_mu_tight"]


def weighted_sum(suffix, sel):
    counts = np.zeros(len(bins) - 1)
    var = np.zeros(len(bins) - 1)
    for jetbin, w in WEIGHTS.items():
        fname = f"WtoLNu2Jets_{jetbin}{suffix}.root"
        a = uproot.open(f"{SAMPLES_DIR}/{fname}")["Events"].arrays(["Ds_mass"] + SEL_BRANCHES, library="np")
        v = a["Ds_mass"][sel(a)]
        h, _ = np.histogram(v, bins=bins)
        print(f"{fname}: {len(v):,} candidates (1-tight-mu), weight={w}, weighted yield={w * len(v):,.1f}")
        counts += w * h
        var += (w ** 2) * h
    return counts, np.sqrt(var)


print("=== OS ===")
os_counts, os_err = weighted_sum("", lambda a: (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0))
print("=== SS ===")
ss_counts, ss_err = weighted_sum("_SS", lambda a: (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0))

with np.errstate(divide="ignore", invalid="ignore"):
    ratio = os_counts / ss_counts
    ratio_err = ratio * np.sqrt((os_err / os_counts) ** 2 + (ss_err / ss_counts) ** 2)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

ax.errorbar(c, os_counts, yerr=os_err, fmt="o", ms=4, color="steelblue",
            label=f"OS (weighted yield {os_counts.sum():,.0f})")
ax.errorbar(c, ss_counts, yerr=ss_err, fmt="s", ms=4, color="firebrick",
            label=f"SS (weighted yield {ss_counts.sum():,.0f})")
ax.set_ylabel("# events (lumi-reweighted)")
ax.legend(fontsize=16)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False)

rax.errorbar(c, ratio, yerr=ratio_err, fmt="o", ms=4, color="black")
rax.axhline(1.0, color="gray", ls="--", lw=1.2)
rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
rax.set_ylabel("OS / SS")

plt.tight_layout()
out = f"{OUTDIR}/wjets_os_vs_ss_weighted_Ds_mass.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

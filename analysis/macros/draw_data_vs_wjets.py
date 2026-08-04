"""Data (OS, SR-blinded) vs the luminosity-reweighted sum of the three WtoLNu2Jets
jet-bin MC samples (0J/1J/2J, OS), absolute scale (no shape normalization) -- the
per-sample weights are the equivalent-luminosity reweighting factors, so summing the
weighted MC directly tests whether the predicted rate/shape matches data, not just
the shape.

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_data_vs_wjets.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS = 1.96835

MC_WEIGHTS = {
    "WtoLNu2Jets_0J.root": 2.3293047885,
    "WtoLNu2Jets_1J.root": 0.5382864853,
    "WtoLNu2Jets_2J.root": 0.2748710631,
}

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)   # 10 MeV bins
c = 0.5 * (bins[1:] + bins[:-1])

# OS: exactly one tight muon in the event (the best muon itself, no extra of either
# charge). os/ss are relative to the pi charge, not the best muon -- see H_phipimu.cc.
SEL_BRANCHES = ["num_os_mu_tight", "num_ss_mu_tight"]


def os_sel(a):
    return (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)


mc_counts = np.zeros(len(bins) - 1)
mc_var = np.zeros(len(bins) - 1)
for fname, w in MC_WEIGHTS.items():
    a = uproot.open(f"{SAMPLES_DIR}/{fname}")["Events"].arrays(["Ds_mass"] + SEL_BRANCHES, library="np")
    v = a["Ds_mass"][os_sel(a)]
    h, _ = np.histogram(v, bins=bins)
    print(f"{fname}: {len(v):,} candidates (1-tight-mu), weight={w}, weighted yield={w * len(v):,.1f}")
    mc_counts += w * h
    mc_var += (w ** 2) * h
mc_err = np.sqrt(mc_var)

data_a = uproot.open(f"{SAMPLES_DIR}/Data_2022.root")["Events"].arrays(["Ds_mass"] + SEL_BRANCHES, library="np")
data_v = data_a["Ds_mass"][os_sel(data_a)]
data_h, _ = np.histogram(data_v, bins=bins)
data_err = np.sqrt(data_h)
print(f"Data_2022.root (OS, SR blinded, 1-tight-mu): {len(data_v):,} candidates")

with np.errstate(divide="ignore", invalid="ignore"):
    ratio = data_h / mc_counts
    ratio_err = ratio * np.sqrt((data_err / data_h) ** 2 + (mc_err / mc_counts) ** 2)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

hep.histplot((mc_counts, bins), ax=ax, histtype="step", color="firebrick",
             label="WtoLNu2Jets (0J+1J+2J, lumi-reweighted)")
ax.fill_between(c, mc_counts - mc_err, mc_counts + mc_err, step="mid", color="firebrick", alpha=0.2)
ax.errorbar(c, data_h, yerr=data_err, fmt="o", ms=4, color="black",
            label=f"Data 2022 (OS, SR blinded, {len(data_v):,} cand.)")
ax.set_ylabel("# events")
ax.legend(fontsize=15)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary")

rax.errorbar(c, ratio, yerr=ratio_err, fmt="o", ms=4, color="black")
rax.axhline(1.0, color="gray", ls="--", lw=1.2)
rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
rax.set_ylabel("Data / MC")
rax.set_ylim(0, 3)

plt.tight_layout()
out = f"{OUTDIR}/data_vs_wjets_Ds_mass.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

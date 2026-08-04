"""OS vs SS comparison of Ds_mass for 2022 data, raw counts (not normalized), using the
deduped (highest-Ds_pt) Samples/ tuples. OS has the |M-M_Ds|<40 MeV signal region blinded
(dropped entirely), so the OS curve has a gap there; SS is unblinded.

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_data_os_vs_ss.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS = 1.96835

# Exactly one tight muon in the event (the best muon itself, no extra of either
# charge). os/ss are relative to the pi charge, which is flipped between the OS/SS
# plugins, so the two selections are each other's mirror -- see H_phipimu(_SS).cc
# and draw_data_vs_wjets{,_ss}.py.
SEL_BRANCHES = ["num_os_mu_tight", "num_ss_mu_tight"]

os_a = uproot.open(f"{SAMPLES_DIR}/Data_2022.root")["Events"].arrays(["Ds_mass"] + SEL_BRANCHES, library="np")
ss_a = uproot.open(f"{SAMPLES_DIR}/Data_2022_SS.root")["Events"].arrays(["Ds_mass"] + SEL_BRANCHES, library="np")
os_ = os_a["Ds_mass"][(os_a["num_os_mu_tight"] == 1) & (os_a["num_ss_mu_tight"] == 0)]
ss_ = ss_a["Ds_mass"][(ss_a["num_ss_mu_tight"] == 1) & (ss_a["num_os_mu_tight"] == 0)]
print(f"OS candidates: {len(os_):,}   SS candidates: {len(ss_):,}  (both after 1-tight-mu cut)")

# Full range the tuplizer actually keeps (|M-M_Ds|<0.3), not a zoomed sub-range.
bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)   # 10 MeV bins
c = 0.5 * (bins[1:] + bins[:-1])


def counts_ep(x):
    h, _ = np.histogram(x, bins=bins)
    return h, np.sqrt(h)


yo, eo = counts_ep(os_)
ys, es = counts_ep(ss_)

with np.errstate(divide="ignore", invalid="ignore"):
    ratio = yo / ys
    ratio_err = ratio * np.sqrt((eo / yo) ** 2 + (es / ys) ** 2)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
ax.errorbar(c, yo, yerr=eo, fmt="o", ms=4, color="steelblue", label=f"OS ({len(os_):,} cand., SR blinded, 1-tight-mu)")
ax.errorbar(c, ys, yerr=es, fmt="s", ms=4, color="firebrick",  label=f"SS ({len(ss_):,} cand., 1-tight-mu)")
ax.set_ylabel("# events")
ax.legend(fontsize=16)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary")

rax.errorbar(c, ratio, yerr=ratio_err, fmt="o", ms=4, color="black")
rax.axhline(1.0, color="gray", ls="--", lw=1.2)
rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
rax.set_ylabel("OS / SS")

plt.tight_layout()
out = f"{OUTDIR}/data_Ds_mass_os_vs_ss.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

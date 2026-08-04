"""OS vs SS comparison of Ds_mass for the WtoLNu2Jets background, area-normalized,
using the deduped (highest-Ds_pt) Samples/ tuples. The point: the OS distribution
carries a small real-Ds peak on top of the smooth combinatorial continuum, while
the SS distribution has no visible peak -- SS reproduces the combinatorial
background but not the real-Ds component (which is charge-correlated, e.g. W+c).

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_wjets_os_vs_ss.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS = 1.96835

os_ = uproot.open(f"{SAMPLES_DIR}/WtoLNu2Jets.root")["Events"].arrays(["Ds_mass"], library="np")["Ds_mass"]
ss_ = uproot.open(f"{SAMPLES_DIR}/WtoLNu2Jets_SS.root")["Events"].arrays(["Ds_mass"], library="np")["Ds_mass"]
print(f"OS candidates: {len(os_):,}   SS candidates: {len(ss_):,}")

# Full range the tuplizer actually keeps (|M-M_Ds|<0.3), not a zoomed sub-range.
# Drawn as points with Poisson error bars (ep style), area-normalized.
bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)   # 10 MeV bins
c = 0.5 * (bins[1:] + bins[:-1])
bw = bins[1] - bins[0]


def norm_ep(x):
    h, _ = np.histogram(x, bins=bins)
    n = len(x)
    y = h / (n * bw)          # area-normalized density
    yerr = np.sqrt(h) / (n * bw)
    return y, yerr


yo, eo = norm_ep(os_)
ys, es = norm_ep(ss_)

ratio = yo / ys
ratio_err = ratio * np.sqrt((eo / yo) ** 2 + (es / ys) ** 2)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
ax.errorbar(c, yo, yerr=eo, fmt="o", ms=4, color="steelblue", label=f"OS ({len(os_):,} cand.)")
ax.errorbar(c, ys, yerr=es, fmt="s", ms=4, color="firebrick",  label=f"SS ({len(ss_):,} cand.)")
ax.set_ylabel("a.u. (area-normalized)")
ax.legend(fontsize=16)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False)

rax.errorbar(c, ratio, yerr=ratio_err, fmt="o", ms=4, color="black")
rax.axhline(1.0, color="gray", ls="--", lw=1.2)
rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
rax.set_ylabel("OS / SS")

plt.tight_layout()
out = f"{OUTDIR}/wjets_Ds_mass_os_vs_ss.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

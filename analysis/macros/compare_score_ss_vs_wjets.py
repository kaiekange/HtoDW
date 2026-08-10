"""v3 (SS-data-trained) DNN: score distribution AND Ds_mass shape, for signal, SS
data (the class it was trained on), W+jets MC (a real background process it never
saw during training), and OS data (the actual real-world application target).

W+jets was the OLD training background (v2); v3 never saw it. If the score still
separates W+jets from signal well, that's evidence the SS-training approach
generalizes to real background processes beyond SS itself, not just to SS's own
quirks. If W+jets scores much higher (more signal-like) than SS does, that would
suggest the model over-fit to something specific to the SS sample rather than
learning genuine signal-vs-background physics. OS data is included as the ultimate
reference: it is real background, blinded in |M-M_Ds|<40 MeV, but otherwise
unrestricted -- no score cut is applied here, so this is the raw, uncut population.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_score_ss_vs_wjets.py
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

MODEL_DIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_ssdata"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
os.makedirs(OUTDIR, exist_ok=True)

M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D = 1.8696

WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}

BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight", "n_b_jets_tight_dR04"]
BR_SIG = BR + ["Ds_match", "mu_match"]


def load(name, ss=False, sig=False):
    br = BR_SIG if sig else BR
    a = uproot.open(f"{MODEL_DIR}/Prediction/{name}.root")["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a["n_b_jets_tight_dR04"] == 0) & (a["score_dsmu"] > -1)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


print("Loading...")
sig = load("HtoDsW", sig=True)
ss = load("Data_2022_SS", ss=True)
os_ = load("Data_2022", ss=False)  # SR-blinded: no |M-M_Ds|<40 MeV events

wj_parts, wj_w = [], []
for name, wgt in WJETS_WEIGHTS.items():
    d = load(name, ss=False)
    wj_parts.append(d)
    wj_w.append(np.full(len(d["Ds_mass"]), wgt))
wj = {k: np.concatenate([p[k] for p in wj_parts]) for k in wj_parts[0]}
wj_weight = np.concatenate(wj_w)

print(f"Signal HtoDsW (OS, truth-matched): {len(sig['score_dsmu']):,}")
print(f"SS data (trained-on background):   {len(ss['score_dsmu']):,}")
print(f"OS data (SR-blinded, real target): {len(os_['score_dsmu']):,}")
print(f"W+jets MC (never trained on):      {len(wj['score_dsmu']):,} (lumi-reweighted)")

print(f"\n{'sample':12s} {'mean':>7s} {'>0.5':>7s} {'>0.9':>7s} {'>0.99':>7s}")
for name, d, w in [("Signal", sig, None), ("SS data", ss, None),
                    ("OS data", os_, None), ("W+jets MC", wj, wj_weight)]:
    s = d["score_dsmu"]
    mean = np.average(s, weights=w)
    print(f"{name:12s} {mean:7.4f} {100*np.average(s>0.5, weights=w):6.2f}% "
          f"{100*np.average(s>0.9, weights=w):6.2f}% {100*np.average(s>0.99, weights=w):6.2f}%")

# ---------- score distribution ----------
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
bins = np.linspace(0, 1, 51)
for ax, logy in ((axes[0], False), (axes[1], True)):
    ax.hist(wj["score_dsmu"], bins=bins, weights=wj_weight, density=True, alpha=0.5,
            color="seagreen", label="W+jets MC (not trained on)")
    ax.hist(ss["score_dsmu"], bins=bins, density=True, histtype="step", lw=2.5,
            color="steelblue", label="SS data (trained on)")
    ax.hist(os_["score_dsmu"], bins=bins, density=True, histtype="step", lw=2.5,
            ls="--", color="black", label="OS data (SR-blinded)")
    ax.hist(sig["score_dsmu"], bins=bins, density=True, alpha=0.4,
            color="firebrick", label="HtoDsW signal")
    ax.set_xlabel("DNN score (v3)")
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=13)
    if logy:
        ax.set_yscale("log")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=12)

plt.tight_layout()
out = f"{OUTDIR}/score_ss_vs_wjets.png"
plt.savefig(out, dpi=130)
print(f"\nSaved {out}")

# ---------- Ds_mass at increasing score cuts: does a real Ds peak grow? ----------
# Only OS data is blinded; SS and W+jets are NOT, and hiding their Ds window would
# throw away the one thing this plot exists to show -- whether the score cut
# enriches real Ds in the very region the fit will use. So: plot the FULL range for
# SS/W+jets, and normalize every sample to its own SIDEBAND integral. The
# combinatorial continuum then lines up across samples, and any real Ds excess
# appears as a bump rising above it. OS data is simply empty in that window.
# Score thresholds come from the signal's own quantiles; all samples are scored by
# the same v3 model, so one efficiency scale is valid across them.
EFFS_MASS = [1.0, 0.90, 0.70, 0.50]
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)   # 20 MeV bins; blind window falls on bin edges
mcent = 0.5 * (mbins[1:] + mbins[:-1])
sb_bin = np.abs(mcent - M_DS) >= SR_HALFWIDTH     # sideband bins, used for normalization
sr_bin = ~sb_bin


def sideband_normed(mass, weight, cut_mask):
    h, _ = np.histogram(mass[cut_mask], bins=mbins,
                        weights=None if weight is None else weight[cut_mask])
    norm = h[sb_bin].sum()
    return h / norm if norm > 0 else h, h


print(f"\n=== SR/sideband ratio (real-Ds enrichment; OS is blind so it is N/A) ===")
print(f"{'eff_sig':>8s} {'SS data':>10s} {'W+jets MC':>11s}")

fig, axes = plt.subplots(1, len(EFFS_MASS), figsize=(9 * len(EFFS_MASS), 8))
for ax, eff in zip(axes, EFFS_MASS):
    thr = -1.0 if eff >= 1.0 else np.quantile(sig["score_dsmu"], 1 - eff)
    c_wj = wj["score_dsmu"] > thr
    c_ss = ss["score_dsmu"] > thr
    c_os = os_["score_dsmu"] > thr

    hw, raw_w = sideband_normed(wj["Ds_mass"], wj_weight, c_wj)
    hs, raw_s = sideband_normed(ss["Ds_mass"], None, c_ss)
    ho, raw_o = sideband_normed(os_["Ds_mass"], None, c_os)

    r_ss = raw_s[sr_bin].sum() / raw_s[sb_bin].sum()
    r_wj = raw_w[sr_bin].sum() / raw_w[sb_bin].sum()
    lbl_eff = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
    print(f"{lbl_eff:>8s} {r_ss:10.4f} {r_wj:11.4f}")

    hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                 label=f"W+jets MC (N={raw_w.sum():,.0f})")
    hep.histplot(hs, mbins, ax=ax, histtype="step", lw=2.5, color="steelblue",
                 label=f"SS data (N={raw_s.sum():,.0f})")
    hep.histplot(ho, mbins, ax=ax, histtype="step", lw=2.5, ls="--", color="black",
                 label=f"OS data, blind (N={raw_o.sum():,.0f})")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel("Normalized to sideband integral")
    lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}% (score>{thr:.3f})"
    ax.set_title(lbl, fontsize=15)
    ax.legend(fontsize=11, loc="upper right")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=10)

plt.tight_layout()
out2 = f"{OUTDIR}/mass_ss_vs_wjets_vs_os.png"
plt.savefig(out2, dpi=120)
print(f"\nSaved {out2}")

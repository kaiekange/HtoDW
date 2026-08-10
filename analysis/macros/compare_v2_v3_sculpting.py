"""v2 (WtoLNu2Jets-trained) vs v3 (SS-data-trained) DNN: does the score cut reshape
OS sideband data the same way, and does more or less combinatorial background
survive at a matched signal efficiency?

Two things this answers:
  1. Shape: M(KKpi) in OS sideband data at increasing score cuts, v2 (top row) vs
     v3 (bottom row) side by side -- the same check as check_mva_sculpting.py, run
     for both models so the D+- peak enrichment can be compared directly.
  2. Combinatorial yield at matched signal efficiency: pick a DNN score threshold
     per model that keeps a fixed fraction of truth-matched OS signal (50/70/90%),
     then in OS sideband data estimate the D+- peak excess via a local sideband
     subtraction (a flat window on each side of the D+- peak) and report
     combinatorial = total - peak_excess. Lower combinatorial at the same signal
     efficiency is a real improvement; a shift that's all in the peak term is not
     (that's real charm, not something a score cut is supposed to remove).

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_v2_v3_sculpting.py
"""

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
os.makedirs(OUTDIR, exist_ok=True)

M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D = 1.8696
PEAK_HALFWIDTH = 0.02
SB_LO_HI = (0.02, 0.06)  # local sideband strip on each side of the D+- peak

MODELS = {
    "v2 (W+jets-trained)": "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets",
    "v3 (SS-data-trained)": "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_ssdata",
}
COLORS = ["black", "steelblue", "darkorange", "firebrick"]
EFFS = [0.90, 0.70, 0.50]  # signal efficiency working points; "no cut" is added separately

BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight", "n_b_jets_tight_dR04"]
BR_SIG = BR + ["Ds_match", "mu_match"]

bins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)  # 20 MeV/bin -- 60 bins was too noisy at the tight-cut working points (N~300)


def load(model_dir, name, sig=False):
    br = BR_SIG if sig else BR
    a = uproot.open(f"{model_dir}/Prediction/{name}.root")["Events"].arrays(br, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a["n_b_jets_tight_dR04"] == 0) & (a["score_dsmu"] > -1))
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


# ---------- 1. shape plot, v2 vs v3, side by side ----------
# Cuts are matched by SIGNAL EFFICIENCY, not raw score -- the two models' score
# distributions are not on the same scale, so "score > 0.99" means a different
# amount of retained signal in each. A fixed-efficiency cut is the fair comparison.
fig, axes = plt.subplots(2, 1, figsize=(11, 16))
for row, (label, model_dir) in enumerate(MODELS.items()):
    ax = axes[row]
    sig = load(model_dir, "HtoDsW", sig=True)
    dat = load(model_dir, "Data_2022")
    for eff, color in zip([1.0] + EFFS, COLORS):
        thr = -1.0 if eff >= 1.0 else np.quantile(sig["score_dsmu"], 1 - eff)
        m = dat["score_dsmu"] > thr
        v = dat["Ds_mass"][m]
        if len(v) < 50:
            continue
        h, _ = np.histogram(v, bins=bins)
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}% (score>{thr:.3f})"
        hep.histplot(h / h.sum(), bins, ax=ax, histtype="step", lw=2, color=color,
                     label=f"{lbl}  (N={len(v):,})")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=15)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=15)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel("Fraction of events / bin")
    ax.set_title(f"OS sideband data -- {label}", fontsize=17)
    ax.legend(fontsize=12)
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=12)

plt.tight_layout()
out = f"{OUTDIR}/compare_v2_v3_sculpting.png"
plt.savefig(out, dpi=120)
print(f"Saved {out}")

# ---------- 2. combinatorial vs D+- peak yield at matched signal efficiency ----------
print(f"\n{'model':24s} {'eff_sig':>8s} {'cut':>8s} {'N_total':>9s} "
      f"{'N_peak_excess(D+-)':>19s} {'N_combinatorial':>16s}")
for label, model_dir in MODELS.items():
    sig = load(model_dir, "HtoDsW", sig=True)
    dat = load(model_dir, "Data_2022")
    for eff in EFFS:
        thr = np.quantile(sig["score_dsmu"], 1 - eff)
        m = dat["score_dsmu"] > thr
        mass = dat["Ds_mass"][m]
        n_total = len(mass)
        n_peak = np.sum(np.abs(mass - M_D) < PEAK_HALFWIDTH)
        n_sb = np.sum((np.abs(mass - M_D) > SB_LO_HI[0]) & (np.abs(mass - M_D) < SB_LO_HI[1]))
        # two SB strips of the same total width as the peak window (2 * 0.04-wide
        # strips vs 1 * 0.04-wide peak window -- average, not sum)
        n_peak_expected_bkg = n_sb / 2
        n_peak_excess = max(n_peak - n_peak_expected_bkg, 0)
        n_comb = n_total - n_peak_excess
        print(f"{label:24s} {100*eff:7.0f}% {thr:8.4f} {n_total:9,} "
              f"{n_peak_excess:19,.0f} {n_comb:16,.0f}")
    print()

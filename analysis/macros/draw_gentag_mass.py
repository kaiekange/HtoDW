"""M(KKpi) split by what the three tracks truly are, versus the DNN working point.

Reads the H_phipimu_GenTag production (WtoLNu2Jets_2022_GenTag, the three preEE OS jet
bins merged) after scoring, and draws the mass spectrum of each truth category.

Two rows, one column per working point:
  top     every category stacked -- the actual composition, absolute lumi-weighted yield
  bottom  only the "all three tracks from ONE charm hadron" categories, overlaid, so the
          peak POSITION of each is readable (they are buried under combinatorics on top)

W+jets MC is not blinded, so the Ds window is drawn: whether a category peaks inside it
is the whole question (see feedback_blind_sideband_care -- only OS data gets the SR dug
out). Absolute yields, never shape-normalized, per feedback_plot_style.

Thresholds come from truth-matched signal (Ds_match & mu_match), the same convention the
other macros use.

NOTE ON NEGATIVE BINS: Lumi_reweight carries the amcatnloFXFX sign (~15.6% of events have
negative weight), so a bin holding only a couple of raw events can sum negative. Bins that
sparse carry no information -- read the shape, not individual bins.

    python3 draw_gentag_mass.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

hep.style.use("CMS")

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL)}/Prediction"

M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.040
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)          # 20 MeV, SR on bin edges
SEL = ["num_os_mu_tight", "num_ss_mu_tight", "n_b_jets_tight_dR04"]
EFFS = [1.0, 0.90, 0.70, 0.50]

# (label, mode values, colour). Order is the stacking order, smallest/most interesting last.
CATS = [
    ("combinatorial (0)",        [0],    "#B8BDC4"),
    ("partial match (8)",        [8],    "#8C9196"),
    (r"$\geq$2 charm ancestors (7)", [7], "#5A6068"),
    (r"D/D$_s\to\pi\pi\pi$ (10)", [10],  "#C46A9E"),
    ("other 1-ancestor (6,9)",   [6, 9], "#9B7FD4"),
    (r"D$\to$KK$\pi$ (4)",       [4],    "#4C9F70"),
    (r"D$\to\phi\pi$ (3)",       [3],    "#2E86AB"),
    (r"D$\to$K$\pi\pi$ misID (5)", [5],  "#E8871A"),
    (r"true D$_s$ (1+2)",        [1, 2], "#D1495B"),
]
ONE_ANCESTOR = [c for c in CATS if c[0].startswith(("D", "true", "other"))]


def sel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    # NOTE: this also moves the b-veto from the hardcoded n_b_jets_tight_dR04 to
    # common.B_VETO (n_bjet_corr), the corrected branch the rest of the chain uses.
    return common.selection_mask(a)


# ---------- signal thresholds, from truth-matched signal ----------
sig = []
for s in ("HtoDsW_2022", "HtoDsstW_2022"):
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) & (a["score_dsmu"] > -1)
    sig.append(a["score_dsmu"][m])
sig = np.concatenate(sig)
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig, 1 - e))) for e in EFFS}
print(f"truth-matched signal: {len(sig):,}   thresholds: "
      + "  ".join(f"{e:.0%}>{t:.4f}" for e, t in THR.items()))

# ---------- background ----------
b = uproot.open(f"{PRED}/WtoLNu2Jets_2022_GenTag.root")["Events"].arrays(
    ["score_dsmu", "cand_gen_mode", "Ds_mass", "Lumi_reweight"] + SEL, library="np")
k = sel(b) & (b["score_dsmu"] > -1)
sc, mode, mass, w = (b["score_dsmu"][k], b["cand_gen_mode"][k],
                     b["Ds_mass"][k], b["Lumi_reweight"][k])
print(f"scored background rows after selection: {k.sum():,}")


def hist(m, cut):
    h, _ = np.histogram(mass[cut & m], bins=mbins, weights=w[cut & m])
    v, _ = np.histogram(mass[cut & m], bins=mbins, weights=w[cut & m] ** 2)
    return h, np.sqrt(v)


def decorate(ax, col):
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        ax.set_ylabel("Weighted events / 20 MeV")


fig, axes = plt.subplots(2, len(EFFS), figsize=(8.5 * len(EFFS), 15), squeeze=False)
for col, eff in enumerate(EFFS):
    cut = sc > THR[eff]
    lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"

    ax = axes[0][col]
    stack, var_stack, labels, colors = [], [], [], []
    for name, ms, c in CATS:
        h, e = hist(np.isin(mode, ms), cut)
        stack.append(h); var_stack.append(e ** 2)
        labels.append(f"{name}  N={h.sum():,.0f}"); colors.append(c)
    hep.histplot(stack, mbins, ax=ax, histtype="fill", stack=True,
                 color=colors, label=labels, edgecolor="white", linewidth=0.4)
    # Error bar on the STACK TOTAL, not per-category: with categories stacked, a
    # per-layer error bar would sit at the wrong height (it belongs at the top of
    # that layer's own band, not the cumulative one) and stacked Poisson/weighted
    # variances don't add along the bars themselves. The total is the one quantity
    # with an unambiguous height and uncertainty here.
    total = np.sum(stack, axis=0)
    total_err = np.sqrt(np.sum(var_stack, axis=0))
    mcent_top = 0.5 * (mbins[1:] + mbins[:-1])
    ax.errorbar(mcent_top, total, yerr=total_err, fmt="none", ecolor="black",
                elinewidth=1.0, capsize=0, alpha=0.7)
    # pad clears the CMS label, which sits above the axes on the first column.
    ax.set_title(f"all categories, stacked -- {lbl}", fontsize=15, pad=26)
    ax.legend(fontsize=10, loc="upper right")
    decorate(ax, col)
    if col == 0:
        # data=False already prepends "Simulation"; text= must not repeat it.
        hep.cms.label(ax=ax, rlabel="2022 preEE (13.6 TeV)", data=False,
                      text="Preliminary", fontsize=11)

    ax = axes[1][col]
    for name, ms, c in ONE_ANCESTOR:
        h, e = hist(np.isin(mode, ms), cut)
        hep.histplot(h, mbins, ax=ax, yerr=e, histtype="step", color=c, linewidth=2.0,
                     label=f"{name}  N={h.sum():,.0f}")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.set_title(f"one common charm ancestor only -- {lbl}", fontsize=15)
    ax.legend(fontsize=10, loc="upper right")
    decorate(ax, col)

plt.tight_layout()
out = common.fig_path("gentag", "gentag_mass_by_mode.png", key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

# ---------------------------------------------------------------------------
# Second figure: the genuine D -> K+K-pi decays on their own.
#
# Modes 3 (D->phi pi) and 4 (D->KKpi direct) are a REAL D+- reconstructed with the
# correct mass hypotheses -- unlike mode 5, which is D->K pi pi with a pion taken for a
# kaon. So this peak is an in-situ measurement of the KKpi mass resolution in exactly the
# final state the signal uses, and the true-Ds peak is drawn alongside for comparison.
#
# Finer bins than the 20 MeV house binning, which is too coarse to read a width off a peak
# this narrow -- but the SAME full mass range as every other spectrum plot. Binning and
# range are independent choices; narrowing the range would hide whatever this category
# does in the tails and would stop the figure being comparable to gentag_mass_by_mode.
# ---------------------------------------------------------------------------
zbins = np.arange(M_DS - 0.3, M_DS + 0.3001, 0.010)
DCATS = [(r"D$\to\phi\pi$ (3)", [3], "#2E86AB"),
         (r"D$\to$KK$\pi$ direct (4)", [4], "#4C9F70"),
         (r"D$\to$KK$\pi$ total (3+4)", [3, 4], "#1B3A54")]


def zhist(m, cut):
    h, _ = np.histogram(mass[cut & m], bins=zbins, weights=w[cut & m])
    v, _ = np.histogram(mass[cut & m], bins=zbins, weights=w[cut & m] ** 2)
    return h, np.sqrt(v)


def peak_stats(m, cut, centre, half=0.060):
    """Weighted mean and RMS inside +-half of centre. Not a fit: with negative
    amcatnloFXFX weights in sparse bins a fit is not robust, and the RMS is enough to
    say whether two peaks have the same width."""
    s = cut & m & (np.abs(mass - centre) < half)
    if w[s].sum() <= 0:
        return float("nan"), float("nan"), 0.0
    mu = np.average(mass[s], weights=w[s])
    var = np.average((mass[s] - mu) ** 2, weights=w[s])
    return mu, np.sqrt(max(var, 0.0)), w[s].sum()


fig, axes = plt.subplots(1, len(EFFS), figsize=(8.5 * len(EFFS), 7.5), squeeze=False)
for col, eff in enumerate(EFFS):
    cut = sc > THR[eff]
    lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
    ax = axes[0][col]
    for name, ms, c in DCATS:
        h, e = zhist(np.isin(mode, ms), cut)
        hep.histplot(h, zbins, ax=ax, yerr=e, histtype="step", color=c,
                     linewidth=2.2 if len(ms) > 1 else 1.6,
                     label=f"{name}  N={h.sum():,.0f}")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        ax.set_ylabel("Weighted events / 10 MeV")
        hep.cms.label(ax=ax, rlabel="2022 preEE (13.6 TeV)", data=False,
                      text="Preliminary", fontsize=11)
    ax.set_title(rf"genuine D$\to$K$^+$K$^-\pi$ -- {lbl}", fontsize=15, pad=26)
    ax.legend(fontsize=11, loc="upper right")
plt.tight_layout()
out2 = common.fig_path("gentag", "gentag_mass_D_to_KKpi.png", key=MODEL)
plt.savefig(out2, dpi=110)
plt.close(fig)
print(f"Saved {out2}")

print(f"\npeak position / RMS within +-60 MeV of the nominal mass:")
print(f"{'':<28}{'eps_sig':>9}{'<M> [GeV]':>12}{'RMS [MeV]':>12}{'yield':>10}")
for eff in EFFS:
    cut = sc > THR[eff]
    for name, ms, centre in ((r"D->KKpi (3+4)", [3, 4], M_D),
                             (r"true Ds (1+2)", [1, 2], M_DS)):
        mu, rms, y = peak_stats(np.isin(mode, ms), cut, centre)
        e = "no cut" if eff >= 1.0 else f"{eff:.0%}"
        print(f"{name:<28}{e:>9}{mu:>12.4f}{1000 * rms:>12.1f}{y:>10,.0f}")

# Peak position of each single-ancestor category, no cut -- the number the plot is for.
print(f"\n{'category':<26}{'N (no cut)':>12}{'in +-40 MeV':>13}{'peak%':>8}{'<M> [GeV]':>11}")
for name, ms, _ in ONE_ANCESTOR:
    m = np.isin(mode, ms)
    n = w[m].sum()
    pk = np.abs(mass - M_DS) < SR_HALF
    print(f"{name:<26}{n:>12,.0f}{w[m & pk].sum():>13,.0f}"
          f"{100 * w[m & pk].sum() / max(n, 1e-9):>7.1f}%{np.average(mass[m], weights=w[m]):>11.4f}")

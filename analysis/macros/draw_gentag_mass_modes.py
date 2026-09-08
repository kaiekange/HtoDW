"""M(KKpi) of one GenTag truth category at a time, versus the DNN working point.

Reads the merged (preEE+postEE) H_phipimu_GenTag production after ssdata scoring
(WtoLNu2Jets_GenTag) and, for each of four truth groups, draws the lumi-weighted mass
spectrum at three signal-efficiency working points.

Four output PNGs, one per group, each a 1x3 grid (eps_sig = 90 / 70 / 50 %):
  mode12       true Ds        (cand_gen_mode 1 + 2)
  mode34       genuine D->KKpi (3 + 4)
  mode5        D->Kpipi misID reflection (5)
  mode_other   everything else (0, 6, 7, 8, 9, 10)

Style follows compare_incl_by_era.py: absolute lumi-weighted yield (never shape
normalized), 10 MeV bins over M_Ds +- 0.3, seagreen fill, blind window greyed, the
D+- line and the D+-/Ds+- labels. Only the one group is drawn per panel -- no
full-background reference overlay (user's call).

Selection is OS + tight-b-veto, matching every other closure macro.

Thresholds come from truth-matched signal (Ds_match & mu_match), both channels. The
quantile is UNWEIGHTED so these panels share an x-axis with the other figures; the
script prints the Lumi_reweight-weighted thresholds alongside as a check that the
choice does not move the cut (the signal file now spans two eras at a 3.34x relative
weight -- project_samples_merge_era_branch).

NEGATIVE BINS: Lumi_reweight carries the amcatnloFXFX sign, so a sparse bin can sum
negative; a zero line is drawn. Read the shape, not individual bins. The raw row
count is printed next to the weighted N= in every legend for the same reason.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_gentag_mass_modes.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL)}/Prediction"

M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.040
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)      # 10 MeV, SR on bin edges
mcent = 0.5 * (mbins[1:] + mbins[:-1])
EFFS = [0.90, 0.70, 0.50]
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]  # merged file spans both

SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "score_dsmu"]

# (filename tag, legend/title label, cand_gen_mode values, colour). Mode ids are wire
# format; filenames key off the ids, not the physics labels, which get renamed.
GROUPS = [
    ("mode12",     r"true D$_s$  (mode 1+2)",            [1, 2],            "#D1495B"),
    ("mode34",     r"genuine D$\to$KK$\pi$  (mode 3+4)", [3, 4],            "#4C9F70"),
    ("mode5",      r"D$\to$K$\pi\pi$ misID  (mode 5)",   [5],               "#E8871A"),
    ("mode_other", r"other  (mode 0,6,7,8,9,10)",        [0, 6, 7, 8, 9, 10], "#B8BDC4"),
]


def sel(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
            & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1))


# ---------- signal thresholds, from truth-matched signal ----------
s_score, s_w = [], []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match", "Lumi_reweight"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s_score.append(a["score_dsmu"][m])
    s_w.append(a["Lumi_reweight"][m])
s_score = np.concatenate(s_score)
s_w = np.concatenate(s_w)


def wquantile(x, w, q):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    return np.interp(q * w.sum(), c, x)


THR = {e: float(np.quantile(s_score, 1 - e)) for e in EFFS}
THR_W = {e: float(wquantile(s_score, s_w, 1 - e)) for e in EFFS}
print(f"truth-matched signal: {len(s_score):,} rows")
print(f"{'eps_sig':>8}{'unweighted':>14}{'lumi-weighted':>16}{'d(score)':>12}")
for e in EFFS:
    print(f"{e:>7.0%}{THR[e]:>14.5f}{THR_W[e]:>16.5f}{THR_W[e]-THR[e]:>12.5f}")

# ---------- background ----------
b = uproot.open(f"{PRED}/WtoLNu2Jets_GenTag.root")["Events"].arrays(
    ["score_dsmu", "cand_gen_mode", "Ds_mass", "Lumi_reweight"] + SEL, library="np")
k = sel(b)
sc, mode, mass, w = (b["score_dsmu"][k], b["cand_gen_mode"][k],
                     b["Ds_mass"][k], b["Lumi_reweight"][k])
print(f"\nscored background rows after selection: {k.sum():,}")


def hist(rowmask):
    h, _ = np.histogram(mass[rowmask], bins=mbins, weights=w[rowmask])
    v, _ = np.histogram(mass[rowmask], bins=mbins, weights=w[rowmask] ** 2)
    return h, np.sqrt(v)


for tag, label, modes, _c in GROUPS:
    in_grp = np.isin(mode, modes)
    fig, axes = plt.subplots(1, len(EFFS), figsize=(8.5 * len(EFFS), 7.5), squeeze=False)
    for col, eff in enumerate(EFFS):
        ax = axes[0][col]
        rowmask = in_grp & (sc > THR[eff])
        h, herr = hist(rowmask)

        # histtype="fill" silently drops yerr (feedback_mplhep_fill_yerr): draw the
        # fill without it, then the stat error as a separate errorbar call.
        hep.histplot(h, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"W+jets GenTag  N={h.sum():,.0f}  (n={rowmask.sum():,})")
        ax.errorbar(mcent, h, yerr=herr, fmt="none", ecolor="black",
                    elinewidth=1.0, capsize=0, alpha=0.8)
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)

        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        lo, hi = ax.get_ylim()
        yt = hi - 0.03 * (hi - lo)
        ax.text(M_D, yt, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, yt, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel("Weighted events / 10 MeV")
            hep.cms.label(ax=ax, rlabel=f"{LUMI/1000:.2f} fb$^{{-1}}$ (13.6 TeV)",
                          data=False, text="Preliminary", fontsize=11)
        # pad lifts the title clear of the CMS label, which sits above the axes on col 0.
        ax.set_title(rf"{label}  --  $\varepsilon_{{sig}}$={100*eff:.0f}%", fontsize=14, pad=26)
        ax.legend(fontsize=11, loc="upper right")

    plt.tight_layout()
    out = common.fig_path("gentag", f"gentag_mass_{tag}_by_cut.png", key=MODEL)
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# One combined figure, same four groups, both views (after gentag_mass_by_mode.png):
#   top row     the four groups stacked -- the composition and its total
#   bottom row  the four groups overlaid as steps -- each group's own shape/peak
# Columns are the three working points. The shape barely moves across the cuts,
# which is the point of showing them side by side.
# ---------------------------------------------------------------------------
def decorate(ax, col):
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    lo, hi = ax.get_ylim()
    yt = hi - 0.03 * (hi - lo)
    ax.text(M_D, yt, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, yt, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        ax.set_ylabel("Weighted events / 10 MeV")


fig, axes = plt.subplots(2, len(EFFS), figsize=(8.5 * len(EFFS), 14), squeeze=False)
for col, eff in enumerate(EFFS):
    cut = sc > THR[eff]
    grp_rows = {tag: np.isin(mode, ms) & cut for tag, _l, ms, _c in GROUPS}
    hh = {tag: hist(m) for tag, m in grp_rows.items()}

    # ---- top: stacked. Grey "other" on the bottom, true Ds on top.
    ax = axes[0][col]
    order = ["mode_other", "mode5", "mode34", "mode12"]
    G = {tag: (lab, c) for tag, lab, _ms, c in GROUPS}
    stack = [hh[t][0] for t in order]
    var = [hh[t][1] ** 2 for t in order]
    labels = [f"{G[t][0]}  N={hh[t][0].sum():,.0f}  (n={grp_rows[t].sum():,})"
              for t in order]
    hep.histplot(stack, mbins, ax=ax, histtype="fill", stack=True,
                 color=[G[t][1] for t in order], label=labels,
                 edgecolor="white", linewidth=0.4)
    total = np.sum(stack, axis=0)
    total_err = np.sqrt(np.sum(var, axis=0))
    ax.errorbar(mcent, total, yerr=total_err, fmt="none", ecolor="black",
                elinewidth=1.0, capsize=0, alpha=0.7)
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.set_title(rf"all four groups, stacked  --  $\varepsilon_{{sig}}$={100*eff:.0f}%",
                 fontsize=14, pad=26)
    ax.legend(fontsize=10, loc="upper right")
    decorate(ax, col)
    if col == 0:
        hep.cms.label(ax=ax, rlabel=f"{LUMI/1000:.2f} fb$^{{-1}}$ (13.6 TeV)",
                      data=False, text="Preliminary", fontsize=11)

    # ---- bottom: overlaid steps, each group's own shape.
    ax = axes[1][col]
    for tag, lab, ms, c in GROUPS:
        h, e = hh[tag]
        hep.histplot(h, mbins, ax=ax, yerr=e, histtype="step", color=c, linewidth=2.0,
                     label=f"{lab}  N={h.sum():,.0f}")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.set_title(rf"all four groups, overlaid  --  $\varepsilon_{{sig}}$={100*eff:.0f}%",
                 fontsize=14, pad=6)
    ax.legend(fontsize=10, loc="upper right")
    decorate(ax, col)

plt.tight_layout()
out = common.fig_path("gentag", "gentag_mass_groups_stack_overlay.png", key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

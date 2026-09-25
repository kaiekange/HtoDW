"""preEE vs postEE for each GenTag truth group, at a common integrated luminosity.

For the four combined groups (cand_gen_mode 1+2 / 3+4 / 5 / rest) this overlays the
2022 preEE and 2022EE postEE M(KKpi) spectra with BOTH scaled to the preEE
luminosity (7.99 fb-1), so the vertical gap between the two curves is the thing to
read: a group sitting higher in postEE is produced-and-reconstructed more often per
unit luminosity there.

WHAT THE RATIO MEANS. W+jets Lumi_reweight is weight_mode="xsec": it already carries
sigma * L / sum_w per era, and the cross sections are identical for the two eras by
construction (make_samples.py uses the same 51846/8860/3284 pb for both pieces).
sum_w divides out. So after rescaling era 2 by L_pre/L_post the postEE/preEE ratio is
a ratio of *selection x reconstruction efficiency* between the two eras (pileup,
detector conditions), NOT of generated cross section. "More generated" in the loose
sense the question was asked, but the honest quantity is efficiency.

Thresholds are computed once from truth-matched signal (both channels) and reused for
both eras. Signal MC is 100% preEE, so a per-era quantile would silently make the two
curves different cuts (compare_incl_by_era.py makes the same point).

The printed table carries a no-cut row: whether postEE/preEE is flat in the cut
(a production/efficiency offset) or moves with it (the classifier responding to the
era differently) are different stories and only the trend separates them.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_gentag_eras.py
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
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)
mcent = 0.5 * (mbins[1:] + mbins[:-1])
EFFS_FIG = [0.90, 0.70, 0.50]
EFFS_TAB = [1.0, 0.90, 0.70, 0.50]
MIN_NEFF = 30

L_PRE, L_POST = common.ERAS[1]["lumi_pb"], common.ERAS[2]["lumi_pb"]
SF_POST = L_PRE / L_POST                       # bring postEE onto the preEE luminosity
print(f"postEE -> preEE luminosity scale factor: {SF_POST:.5f}")

SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "score_dsmu"]

GROUPS = [
    ("mode12",     r"true D$_s$  (mode 1+2)",            [1, 2]),
    ("mode34",     r"genuine D$\to$KK$\pi$  (mode 3+4)", [3, 4]),
    ("mode5",      r"D$\to$K$\pi\pi$ misID  (mode 5)",   [5]),
    ("mode_other", r"other  (mode 0,6,7,8,9,10)",        [0, 6, 7, 8, 9, 10]),
]


def sel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a) & (a["score_dsmu"] > -1)


# ---------- thresholds from truth-matched signal ----------
s_score = []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s_score.append(a["score_dsmu"][m])
s_score = np.concatenate(s_score)
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(s_score, 1 - e))) for e in EFFS_TAB}
print(f"truth-matched signal: {len(s_score):,} rows")
for e in EFFS_TAB:
    if e < 1.0:
        print(f"  eps_sig={e:.0%}  ->  score > {THR[e]:.5f}")

# ---------- background, split by era ----------
b = uproot.open(f"{PRED}/WtoLNu2Jets_GenTag.root")["Events"].arrays(
    ["score_dsmu", "cand_gen_mode", "Ds_mass", "Lumi_reweight", "era", "jet_bin"] + SEL,
    library="np")
k = sel(b) & (np.abs(b["Ds_mass"] - M_DS) < 0.3)
sc, mode, mass, w, era, jb = (b["score_dsmu"][k], b["cand_gen_mode"][k], b["Ds_mass"][k],
                              b["Lumi_reweight"][k].astype(float), b["era"][k], b["jet_bin"][k])
w[era == 2] *= SF_POST                          # postEE now on the preEE luminosity

is_pre, is_post = era == 1, era == 2
print(f"\nrows after selection (|M-M_Ds|<0.3): preEE {is_pre.sum():,}  postEE {is_post.sum():,}")

# ---------- diagnostic: per era x jet bin, catches a missing merge piece ----------
print(f"\n{'era':>7}{'jet_bin':>9}{'rows':>12}{'weighted (@preEE L)':>22}")
for e, etag in ((1, "preEE"), (2, "postEE")):
    for j in (0, 1, 2):
        s = (era == e) & (jb == j)
        print(f"{etag:>7}{j:>9}{s.sum():>12,}{w[s].sum():>22,.1f}")


def integ(m):
    """Weighted yield, effective entries, and stat error for a row mask."""
    ww = w[m]
    sw, sw2 = ww.sum(), (ww ** 2).sum()
    neff = sw ** 2 / sw2 if sw2 > 0 else 0.0
    return sw, np.sqrt(sw2), neff


# ---------- table ----------
print("\n" + "=" * 92)
print("postEE / preEE at common luminosity  (= ratio of selection x reco efficiency)")
print("=" * 92)
print(f"{'group':<22}{'eps_sig':>8}{'N preEE':>12}{'N postEE':>12}"
      f"{'post/pre':>16}{'N_eff pre':>11}{'N_eff post':>11}")
for tag, label, modes in [("all", "all candidates", None)] + GROUPS:
    ing = np.ones(len(mode), bool) if modes is None else np.isin(mode, modes)
    print("-" * 92)
    for e in EFFS_TAB:
        cut = sc > THR[e]
        npre, epre, nfpre = integ(ing & cut & is_pre)
        npost, epost, nfpost = integ(ing & cut & is_post)
        if npre > 0 and npost > 0:
            r = npost / npre
            rerr = r * np.sqrt((epost / npost) ** 2 + (epre / npre) ** 2)
            rtxt = f"{r:6.3f} +- {rerr:5.3f}"
        else:
            rtxt = "         --     "
        flag = " <-" if min(nfpre, nfpost) < MIN_NEFF else ""
        lbl = "no cut" if e >= 1.0 else f"{e:.0%}"
        print(f"{label:<22}{lbl:>8}{npre:>12,.0f}{npost:>12,.0f}{rtxt:>16}"
              f"{nfpre:>11,.0f}{nfpost:>11,.0f}{flag}")


# ---------- figures: one per group, 1x3 over the working points ----------
def hist(m):
    h, _ = np.histogram(mass[m], bins=mbins, weights=w[m])
    v, _ = np.histogram(mass[m], bins=mbins, weights=w[m] ** 2)
    return h, np.sqrt(v)


RATIO_YMAX = 3.0     # fixed, not autoscaled: feedback_plot_style's autoscale incident

# modes=None -> every candidate, no truth split (the overall shape).
# Two rows: top the overlaid spectra, bottom the per-bin postEE/preEE ratio at the
# common luminosity (so 1.0 = the two eras agree bin for bin).
for tag, label, modes in [("all", "all candidates (no mode selection)", None)] + GROUPS:
    ing = np.ones(len(mode), bool) if modes is None else np.isin(mode, modes)
    fig, axes = plt.subplots(2, len(EFFS_FIG), figsize=(8.5 * len(EFFS_FIG), 10.0),
                             squeeze=False, sharex="col",
                             gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.06))
    for col, e in enumerate(EFFS_FIG):
        ax, rax = axes[0][col], axes[1][col]
        cut = sc > THR[e]
        hpre, epre = hist(ing & cut & is_pre)
        hpost, epost = hist(ing & cut & is_post)
        npre, _, nfpre = integ(ing & cut & is_pre)
        npost, _, nfpost = integ(ing & cut & is_post)

        # ---- top: overlaid spectra ----
        hep.histplot(hpre, mbins, ax=ax, yerr=epre, histtype="step", color="steelblue",
                     linewidth=2.0, label=f"2022 preEE   N={npre:,.0f}")
        hep.histplot(hpost, mbins, ax=ax, yerr=epost, histtype="step", color="firebrick",
                     linewidth=2.0, label=f"2022EE postEE (scaled)   N={npost:,.0f}")
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)
        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        lo, hi = ax.get_ylim()
        yt = hi - 0.03 * (hi - lo)
        ax.text(M_D, yt, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, yt, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)

        rtxt = "post/pre = --"
        if npre > 0 and npost > 0:
            rerr = (npost / npre) * np.sqrt((np.sqrt((epost ** 2).sum()) / npost) ** 2
                                            + (np.sqrt((epre ** 2).sum()) / npre) ** 2)
            rtxt = rf"post/pre = {npost/npre:.2f} $\pm$ {rerr:.2f}"
        warn = "\n(N_eff < 30)" if min(nfpre, nfpost) < MIN_NEFF else ""
        ax.text(0.03, 0.97, rtxt + warn, transform=ax.transAxes, va="top", fontsize=12,
                bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.85))
        if col == 0:
            ax.set_ylabel(r"Weighted events / 10 MeV  (@ 7.99 fb$^{-1}$)")
            hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False,
                          text="Preliminary", fontsize=11)
        ax.set_title(rf"{label}  --  $\varepsilon_{{sig}}$={100*e:.0f}%", fontsize=14, pad=26)
        ax.legend(fontsize=11, loc="upper right")

        # ---- bottom: per-bin ratio ----
        # Only bins where both eras have a positive weighted content: a negative or
        # empty bin (which tight cuts do produce) has no meaningful ratio.
        ok = (hpre > 0) & (hpost > 0)
        r = np.where(ok, hpost / np.where(hpre > 0, hpre, 1), np.nan)
        rerr_bin = np.where(ok, r * np.sqrt((epost / np.where(hpost > 0, hpost, 1)) ** 2
                                            + (epre / np.where(hpre > 0, hpre, 1)) ** 2), np.nan)
        inside = ok & (r <= RATIO_YMAX)
        rax.errorbar(mcent[inside], r[inside], yerr=rerr_bin[inside], fmt="s",
                     color="black", markersize=4.0, elinewidth=1.2)
        above = ok & (r > RATIO_YMAX)
        if above.any():
            rax.plot(mcent[above], np.full(above.sum(), RATIO_YMAX * 0.97), "^",
                     color="black", markersize=7)
        rax.axhline(1.0, color="firebrick", lw=1.4)
        rax.set_ylim(0, RATIO_YMAX)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        rax.axvline(M_D, color="green", ls=":", lw=1.5)
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            rax.set_ylabel("postEE / preEE")
        n_bad = int((~ok).sum())
        notes = [f"{n_bad} bins: no ratio"] if n_bad else []
        if above.any():
            notes.append(f"{int(above.sum())} above axis")
        if notes:
            rax.text(0.98, 0.93, "    ".join(notes), transform=rax.transAxes,
                     fontsize=9, color="firebrick", va="top", ha="right")

    out = common.fig_path("era", f"gentag_era_compare_{tag}.png", key=MODEL)
    plt.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

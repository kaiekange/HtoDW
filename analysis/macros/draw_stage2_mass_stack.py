"""M(KKpi) with vs without stage-2, at matched signal efficiency, W+jets split by GenTag.

Two logical rows, four columns (total signal efficiency 90 / 70 / 50 / 30 %):

  row 1  stage-1 only              score_dsmu > THR1(eps)
  row 2  stage-1 (WP1*=eps90) + stage-2   score_dsmu > THR1(0.90)  AND
                                          score_dsmu_res > THR2(eps / 0.90)

THR2 is a quantile of score_dsmu_res on the truth-matched signal that ALREADY passed
stage-1 at eps90, so row-2 column c keeps exactly the same fraction of signal as
row-1 column c -- the two rows are a like-for-like comparison of the Ds-mass
background at equal signal efficiency. At 90% stage-2 is wide open, so the two rows'
first column is identical by construction (the sanity anchor).

W+jets (WtoLNu2Jets_GenTag, both eras, lumi-weighted, absolute yield) is drawn as a
stack of the four GenTag truth groups. OS data (SR blinded) overlaid as points; the
sub-panel is OS data / total-MC in the sideband bins only.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_stage2_mass_stack.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

S1_KEY = "ssdata"
S2_DIR = common.out_dir_for("wjetsres", tool="BDT", seed=42)
PRED = f"{S2_DIR}/Prediction"                    # holds score_dsmu AND score_dsmu_res
M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.040
WP1_STAR = 0.90
EFFS = [0.90, 0.70, 0.50, 0.30]
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)
mc = 0.5 * (mbins[1:] + mbins[:-1])
SB = np.abs(mc - M_DS) >= SR_HALF
SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO]
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]

GROUPS = [   # (label, modes, colour) -- stack bottom -> top
    (r"other (mode 0,6-10)",        [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",  [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],             "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],             "#c0392b"),
]


def ossel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a) & (a["score_dsmu"] > -1)


def rd(name, extra):
    return uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["Ds_mass", "score_dsmu", "score_dsmu_res"] + SEL + extra, library="np")


# ---- thresholds ----
s1_all, s2_ps1 = [], []
for s in common.SIGNALS:
    a = rd(s, ["Ds_match", "mu_match"])
    m = ossel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s1_all.append(a["score_dsmu"][m])
s1_all = np.concatenate(s1_all)
THR1 = {e: float(np.quantile(s1_all, 1 - e)) for e in EFFS}
THR1_STAR = float(np.quantile(s1_all, 1 - WP1_STAR))

for s in common.SIGNALS:
    a = rd(s, ["Ds_match", "mu_match"])
    m = (ossel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
         & (a["score_dsmu"] > THR1_STAR))
    s2_ps1.append(a["score_dsmu_res"][m])
s2_ps1 = np.concatenate(s2_ps1)
THR2 = {e: (-1.0 if e >= WP1_STAR else float(np.quantile(s2_ps1, 1 - e / WP1_STAR)))
        for e in EFFS}
print(f"signal rows: stage-1 pool {len(s1_all):,}, post-WP1* pool {len(s2_ps1):,}")
for e in EFFS:
    print(f"  eps {e:.0%}:  THR1 {THR1[e]:.4f}   |   WP1* + THR2 {THR2[e]:.4f}")

# ---- W+jets GenTag ----
g = rd("WtoLNu2Jets_GenTag", ["cand_gen_mode", "Lumi_reweight"])
gk = ossel(g)
gm, gmass, g1, g2, gw = (g["cand_gen_mode"][gk], g["Ds_mass"][gk], g["score_dsmu"][gk],
                         g["score_dsmu_res"][gk], g["Lumi_reweight"][gk].astype(float))

# ---- OS data ----
d = rd("Data", [])
dk = ossel(d) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALF)
dmass, d1, d2 = d["Ds_mass"][dk], d["score_dsmu"][dk], d["score_dsmu_res"][dk]


def stack_hists(pass_mc):
    hs, errs = [], []
    for _, modes, _ in GROUPS:
        mm = pass_mc & np.isin(gm, modes)
        h, _ = np.histogram(gmass[mm], bins=mbins, weights=gw[mm])
        e2, _ = np.histogram(gmass[mm], bins=mbins, weights=gw[mm] ** 2)
        hs.append(h); errs.append(np.sqrt(e2))
    return hs, errs


fig = plt.figure(figsize=(7.0 * len(EFFS), 13.0))
outer = gridspec.GridSpec(2, len(EFFS), figure=fig, hspace=0.34, wspace=0.24)
ROWS = [("stage-1 only", None), ("stage-1 (WP1*) + stage-2", None)]

for row in (0, 1):
    for col, e in enumerate(EFFS):
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

        if row == 0:
            pmc = g1 > THR1[e]
            pdt = d1 > THR1[e]
        else:
            pmc = (g1 > THR1_STAR) & (g2 > THR2[e])
            pdt = (d1 > THR1_STAR) & (d2 > THR2[e])

        hs, errs = stack_hists(pmc)
        tot = np.sum(hs, axis=0)
        tot_err = np.sqrt(np.sum(np.square(errs), axis=0))   # sqrt(sum w^2), per bin
        bottom = np.zeros(len(mc))
        for (lab, _, cl), h in zip(GROUPS, hs):
            ax.bar(mc, h, width=0.010, bottom=bottom, color=cl, alpha=0.9,
                   label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
            bottom = bottom + h
        # MC stat uncertainty on the stack total: matters here -- amcatnloFXFX with a
        # 27% negative-weight fraction and N_eff only O(100) in the tight columns.
        ax.bar(mc, 2 * tot_err, width=0.010, bottom=tot - tot_err, color="none",
               hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
               label="W+jets MC stat")
        dh, _ = np.histogram(dmass[pdt], bins=mbins)
        dh_m = np.where(SB, dh, np.nan)
        ax.errorbar(mc, dh_m, yerr=np.sqrt(np.where(SB, dh, 0)), fmt="o", ms=4,
                    color="black", elinewidth=1.2, label=f"OS data (SR blind)  N={dh[SB].sum():,.0f}")
        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
        ax.axvline(M_D, color="green", ls=":", lw=1.4)
        ax.set_xlim(mbins[0], mbins[-1])
        peak = max(bottom.max(), np.nanmax(np.where(SB, dh, np.nan)))
        ax.set_ylim(0, peak * 1.45)
        ax.tick_params(labelbottom=False)
        if col == 0:
            ax.set_ylabel("Events / 10 MeV")
        ax.set_title(rf"{ROWS[row][0]}  --  $\varepsilon_{{sig}}$={e:.0%}", fontsize=14)
        ax.legend(fontsize=8.5, loc="upper right")
        if row == 0 and col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary",
                          rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

        # band = MC stat (sqrt(sum w^2) / tot); point error bar = data stat only.
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(SB & (tot > 0), dh / tot, np.nan)
            re = np.where(SB & (tot > 0), np.sqrt(np.where(SB, dh, 0)) / tot, np.nan)
            mc_band = np.where(SB & (tot > 0), tot_err / tot, np.nan)
        rax.bar(mc, 2 * mc_band, width=0.010, bottom=1 - mc_band, color="none",
                hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(mc, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        rax.set_ylim(0.0, 1.99)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.set_xlim(mbins[0], mbins[-1])
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            rax.set_ylabel("data / MC")

out = common.fig_path("comparisons", "stage2_mass_stack.png", key="wjetsres", tool="BDT")
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

# ---- printed yields: resonant SR yield vs eps, both rows ----
print(f"\n{'eps':>5} {'row':>22} {'modes1-5 SR':>13} {'modes1-2 SR':>13} {'nonres SR':>11} {'OS data SB':>11}")
for e in EFFS:
    for row, lab in ((0, "stage-1 only"), (1, "stage-1+stage-2")):
        if row == 0:
            pmc, pdt = g1 > THR1[e], d1 > THR1[e]
        else:
            pmc = (g1 > THR1_STAR) & (g2 > THR2[e]); pdt = (d1 > THR1_STAR) & (d2 > THR2[e])
        sr = np.abs(gmass - M_DS) < SR_HALF
        y15 = gw[pmc & sr & np.isin(gm, [1, 2, 3, 4, 5])].sum()
        y12 = gw[pmc & sr & np.isin(gm, [1, 2])].sum()
        ynr = gw[pmc & sr & np.isin(gm, [0, 6, 7, 8, 9, 10])].sum()
        dsb = int((d1[pdt] > -np.inf).sum()) if row == 0 else int(pdt.sum())
        print(f"{e:>5.0%} {lab:>22} {y15:>13.1f} {y12:>13.1f} {ynr:>11.1f} {dsb:>11d}")

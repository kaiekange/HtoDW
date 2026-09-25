"""Ds_mass and phi_mass at a single fixed 2D working point (WP1, WP2).

Same drawing conventions as draw_stage2_mass_stack.py (GenTag-mode stack for W+jets,
OS data overlaid as points with the SR blinded, data/MC ratio sub-panel) but for one
WP instead of a 4-column efficiency scan, and two mass variables side by side instead
of one.

WP2 IS UNCONDITIONAL (2026-09-22, matching fom_map_stage2_run3.py): thr2 is the
quantile of score_dsmu_res over ALL truth-matched signal, not over the subset that
already passed thr1, so WP2 names one fixed stage-2 cut whatever WP1 is. eps_total
is therefore NOT WP1*WP2 -- the two cuts are correlated -- and is counted directly
from the AND. Both mass variables use the SAME event selection (score_dsmu > thr1 AND
score_dsmu_res > thr2) -- phi_mass is not itself blinded, but events in the Ds_mass
blind window are excluded from OS data for BOTH panels, because blindness is a
property of the event, not of which variable happens to be on the x-axis.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_wp_mass.py [WP1] [WP2]        # defaults 0.40 0.55
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

WP1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.40
WP2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.55

S2_DIR = common.out_dir_for("wjetsres", tool="BDT", seed=42)
PRED = f"{S2_DIR}/Prediction"
M_DS, M_PHI, SR_HALF = 1.96835, 1.019460, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
SEL = [*common.SELECTION]

GROUPS = [   # (label, modes, colour) -- stack bottom -> top, same as stage2_mass_stack
    (r"other (mode 0,6-10)",        [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",  [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],             "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],             "#c0392b"),
]

VARS = [   # (branch, label, bin edges, blind half-width or None, ref line or None)
    ("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]",
     np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010), SR_HALF, None),
    # range must clear the physical bounds (2m_K ~ 0.98735 lower, M_phi+0.05 = 1.06946
    # upper post-fit cut) with margin, or np.histogram silently drops edge events as
    # overflow -- caught by the N mismatch between this panel and the Ds_mass one,
    # which use the same underlying (dk-filtered) event set and must agree.
    ("phi_mass", r"$M(K^+K^-)$ [GeV]",
     np.arange(0.984, 1.074 + 1e-9, 0.002), None, M_PHI),
]


def ossel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return (common.selection_mask(a) & (a["score_dsmu"] > -1)
            & (a["score_dsmu_res"] > -1))


def rd(name, extra):
    return uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["Ds_mass", "phi_mass", "score_dsmu", "score_dsmu_res"] + SEL + extra,
        library="np")


# ---- thresholds: both unconditional, on all truth-matched signal ----
s1, s2, sw = [], [], []
for s in common.SIGNALS:
    a = rd(s, ["Ds_match", "mu_match", "weight"])
    m = ossel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s1.append(a["score_dsmu"][m])
    s2.append(a["score_dsmu_res"][m])
    sw.append(a["weight"][m])
s1, s2, sw = np.concatenate(s1), np.concatenate(s2), np.concatenate(sw)
thr1 = float(np.quantile(s1, 1 - WP1))
# WP2 >= 1 means "stage-1 only": no stage-2 cut at all. Taking the 0th quantile is NOT
# the same thing -- `> min(score)` would still drop the single lowest-scoring row.
STAGE1_ONLY = WP2 >= 1.0
thr2 = -np.inf if STAGE1_ONLY else float(np.quantile(s2, 1 - WP2))
eps_tot = float(sw[(s1 > thr1) & (s2 > thr2)].sum() / sw.sum())
print(f"WP1={WP1:.3f} (thr1={thr1:.4f})  "
      + ("stage-1 ONLY, no stage-2 cut" if STAGE1_ONLY else
         f"WP2={WP2:.2f} unconditional (thr2={thr2:.4f})")
      + f"  eps_total={eps_tot:.4f}")

# ---- W+jets GenTag: apply the 2D cut once, keep both mass branches ----
g = rd("WtoLNu2Jets_merged", ["cand_gen_mode", "weight"])
gk = ossel(g) & (g["score_dsmu"] > thr1) & (g["score_dsmu_res"] > thr2)
g_ds, g_phi, gm, gw = (g["Ds_mass"][gk], g["phi_mass"][gk],
                       g["cand_gen_mode"][gk], g["weight"][gk].astype(float))

# ---- OS data: same 2D cut, blinded on Ds_mass for BOTH panels (see docstring) ----
d = rd("Data_OS_merged", [])
dk = (ossel(d) & (d["score_dsmu"] > thr1) & (d["score_dsmu_res"] > thr2)
      & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALF))
d_ds, d_phi = d["Ds_mass"][dk], d["phi_mass"][dk]
print(f"W+jets rows passing cut: {gk.sum():,}   OS data rows passing cut "
      f"(Ds SR blinded): {dk.sum():,}")

fig = plt.figure(figsize=(15.0, 7.5))
outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.26)

for col, (branch, xlabel, bins, blind_half, refline) in enumerate(VARS):
    mc_vals = g_ds if branch == "Ds_mass" else g_phi
    d_vals = d_ds if branch == "Ds_mass" else d_phi
    bc = 0.5 * (bins[1:] + bins[:-1])
    width = bins[1] - bins[0]
    blind = np.abs(bc - M_DS) < blind_half if blind_half else np.zeros_like(bc, dtype=bool)
    vis = ~blind

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[col], height_ratios=[3, 1], hspace=0.05)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

    hs, errs = [], []
    for _, modes, _ in GROUPS:
        mm = np.isin(gm, modes)
        h, _ = np.histogram(mc_vals[mm], bins=bins, weights=gw[mm])
        e2, _ = np.histogram(mc_vals[mm], bins=bins, weights=gw[mm] ** 2)
        hs.append(h); errs.append(np.sqrt(e2))
    tot = np.sum(hs, axis=0)
    tot_err = np.sqrt(np.sum(np.square(errs), axis=0))
    bottom = np.zeros(len(bc))
    for (lab, _, cl), h in zip(GROUPS, hs):
        ax.bar(bc, h, width=width, bottom=bottom, color=cl, alpha=0.9,
               label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
        bottom = bottom + h
    ax.bar(bc, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
           label="W+jets MC stat")

    dh, _ = np.histogram(d_vals, bins=bins)
    assert dh.sum() == len(d_vals), (
        f"{branch}: histogram dropped {len(d_vals) - dh.sum()} events as under/overflow "
        f"-- bin range [{bins[0]:.4f},{bins[-1]:.4f}] does not cover the data.")
    dh_m = np.where(vis, dh, np.nan)
    ax.errorbar(bc, dh_m, yerr=np.sqrt(np.where(vis, dh, 0)), fmt="o", ms=4,
                color="black", elinewidth=1.2,
                label=f"OS data{' (SR blind)' if blind_half else ''}  N={dh[vis].sum():,.0f}")
    if blind_half:
        ax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15, zorder=0)
    if refline:
        ax.axvline(refline, color="green", ls=":", lw=1.4)
    ax.set_xlim(bins[0], bins[-1])
    peak = max(bottom.max(), np.nanmax(dh_m)) if np.isfinite(np.nanmax(dh_m)) else bottom.max()
    ax.set_ylim(0, peak * 1.45)
    ax.tick_params(labelbottom=False)
    if col == 0:
        ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
    ax.legend(fontsize=8.5, loc="upper right")
    if col == 0:
        hep.cms.label(ax=ax, data=True, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(vis & (tot > 0), dh / tot, np.nan)
        re = np.where(vis & (tot > 0), np.sqrt(np.where(vis, dh, 0)) / tot, np.nan)
        mc_band = np.where(vis & (tot > 0), tot_err / tot, np.nan)
    rax.bar(bc, 2 * mc_band, width=width, bottom=1 - mc_band, color="none",
            hatch="////", edgecolor="0.35", linewidth=0.0)
    rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
    rax.axhline(1.0, color="red", lw=1.0)
    if blind_half:
        rax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15)
    rax.set_ylim(0.0, 1.99)
    rax.set_yticks([0.5, 1.0, 1.5])
    rax.set_xlim(bins[0], bins[-1])
    rax.set_xlabel(xlabel)
    if col == 0:
        rax.set_ylabel("data / MC")

fig.suptitle((rf"stage-1 only, WP1 $\varepsilon_1$={WP1:.1%}" if STAGE1_ONLY else
              rf"WP1 $\varepsilon_1$={WP1:.0%}  +  WP2 $\varepsilon_2$={WP2:.0%}")
             + rf"  $\to$  $\varepsilon_{{tot}}$={eps_tot:.1%}", fontsize=14)

out = common.fig_path("comparisons", (f"wp_mass_wp1_{WP1:.3f}_stage1only.png" if STAGE1_ONLY else
                        f"wp_mass_wp1_{WP1:.2f}_wp2_{WP2:.2f}.png"),
                      key="wjetsres", tool="BDT")
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

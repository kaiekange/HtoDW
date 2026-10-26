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
property of the event, not of which variable happens to be on the x-axis. Since 1026
the MC in the phi panel drops the Ds-SR rows too (before, only data did, so the phi
data/MC was biased low).

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_wp_mass.py [WP1] [WP2] [--add TT,DY] [--key2=resdy]   # defaults 0.40 0.55

2026-10-26: --key2 picks the stage-2 BDT (its Prediction/ carries both scores; default
resdy). MC weight = skim weight x common.br_reweight (DY also x dy_jetbin_sf), written
in the title. Stack bottom -> top: DY, tt, W+jets non-peak, peaks (user, 1022). Bins as
the tightest column of draw_stage1_mass_by_score_cut.py: Ds_mass 52 MeV (edges on the
blind window, SR one bin), phi_mass 8 MeV around the +-20 MeV phi window.

--add stacks extra background MC on top of the W+jets GenTag groups (TT = TT_merged,
DY = DY_merged; each must have been scored by stage-1 AND stage-2), the ratio uses the
full MC sum, and a low / D+- / high window table is printed. Same option as
draw_stage1_mass_by_score_cut.py.
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

EXTRA_DEF = {   # key -> (Prediction file, legend label, colour); same as draw_stage1_mass_by_score_cut.py
    "TT": ("TT_merged", r"$t\bar{t}$ (2L2Nu + LNu2Q)", "#7f5fbf"),
    "DY": ("DY_merged", r"DY$\to\ell\ell$ (0J+1J+2J)", "#4a90c2"),
}
_argv = sys.argv[1:]
KEY2 = next((a.split("=", 1)[1] for a in _argv if a.startswith("--key2=")), "resdy")
# --ymax=DS,PHI: fixed y-axis tops for the two panels (to compare models side by side)
YMAX = next((tuple(float(x) for x in a.split("=", 1)[1].split(",")) for a in _argv
             if a.startswith("--ymax=")), None)
_argv = [a for a in _argv if not a.startswith(("--key2=", "--ymax="))]
ADD = []
if "--add" in _argv:
    i = _argv.index("--add")
    ADD = [x for x in _argv[i + 1].split(",") if x]
    _argv = _argv[:i] + _argv[i + 2:]
_bad = [k for k in ADD if k not in EXTRA_DEF]
if _bad:
    sys.exit(f"--add: unknown {_bad}; choices {list(EXTRA_DEF)}")
WP1 = float(_argv[0]) if len(_argv) > 0 else 0.40
WP2 = float(_argv[1]) if len(_argv) > 1 else 0.55

S2_DIR = common.out_dir_for(KEY2, tool="BDT", seed=42)
S1_KEY = common.resolve(KEY2)["stage1_cut"]["key"]
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

PHI_WIN, PHI_W, DS_NSIDE = 0.020, 0.008, 5
# Ds_mass: 5 bins per 0.26 GeV sideband (52 MeV), the 80 MeV SR one bin (blind anyway).
DS_EDGES = np.concatenate([np.linspace(M_DS - 0.3, M_DS - SR_HALF, DS_NSIDE + 1),
                           np.linspace(M_DS + SR_HALF, M_DS + 0.3, DS_NSIDE + 1)])
# phi_mass: edges on M_phi +- 20 MeV plus two empty bins either side; the window cut is
# in the skim, and float32 rows a hair outside it are clipped onto it below.
PHI_EDGES = np.linspace(M_PHI - PHI_WIN - 2 * PHI_W, M_PHI + PHI_WIN + 2 * PHI_W,
                        int(round((2 * PHI_WIN + 4 * PHI_W) / PHI_W)) + 1)
VARS = [   # (branch, label, bin edges, blind half-width or None, ref line or None)
    ("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", DS_EDGES, SR_HALF, None),
    ("phi_mass", r"$M(K^+K^-)$ [GeV]", PHI_EDGES, None, M_PHI),
]


def phi_clip(x):
    return np.clip(x, M_PHI - PHI_WIN + 1e-7, M_PHI + PHI_WIN - 1e-7)


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
BRC = ["cand_gen_mode", "cand_gen_common_hf_pdgId", "weight"]
g = rd("WtoLNu2Jets_merged", BRC)
gk = ossel(g) & (g["score_dsmu"] > thr1) & (g["score_dsmu_res"] > thr2)
g = {k: v[gk] for k, v in g.items()}
g_ds, g_phi, gm = g["Ds_mass"], phi_clip(g["phi_mass"]), g["cand_gen_mode"]
gw = g["weight"].astype(float) * common.br_reweight(g)

# ---- extra background MC (--add): same 2D cut, never SR-excluded ----
EXTRA = {}
for key in ADD:
    x = rd(EXTRA_DEF[key][0], BRC + (["jet_bin"] if key == "DY" else []))
    xk = ossel(x) & (x["score_dsmu"] > thr1) & (x["score_dsmu_res"] > thr2)
    x = {k: v[xk] for k, v in x.items()}
    xw = x["weight"].astype(float) * common.br_reweight(x)
    if key == "DY":
        xw = xw * common.dy_jetbin_sf(x)
    EXTRA[key] = dict(Ds_mass=x["Ds_mass"], phi_mass=phi_clip(x["phi_mass"]), w=xw)
    print(f"{key} rows passing cut: {xk.sum():,}  sum(w)={EXTRA[key]['w'].sum():,.1f}")

# ---- OS data: same 2D cut, blinded on Ds_mass for BOTH panels (see docstring) ----
d = rd("Data_OS_merged", [])
dk = (ossel(d) & (d["score_dsmu"] > thr1) & (d["score_dsmu_res"] > thr2)
      & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALF))
d_ds, d_phi = d["Ds_mass"][dk], phi_clip(d["phi_mass"][dk])
print(f"W+jets rows passing cut: {gk.sum():,}   OS data rows passing cut "
      f"(Ds SR blinded): {dk.sum():,}")

fig = plt.figure(figsize=(15.0, 7.5))
outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.26)

for col, (branch, xlabel, bins, blind_half, refline) in enumerate(VARS):
    # phi panel: MC drops Ds-SR rows too, like the blinded data (was data-only before 1026)
    sb = (lambda m: np.ones(len(m), bool)) if branch == "Ds_mass" else (lambda m: np.abs(m - M_DS) >= SR_HALF)
    mc_vals = g_ds if branch == "Ds_mass" else g_phi
    d_vals = d_ds if branch == "Ds_mass" else d_phi
    bc = 0.5 * (bins[1:] + bins[:-1])
    width = np.diff(bins)
    blind = np.abs(bc - M_DS) < blind_half if blind_half else np.zeros_like(bc, dtype=bool)
    vis = ~blind

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[col], height_ratios=[3, 1], hspace=0.05)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

    # stack bottom -> top: DY, tt (--add order reversed so DY is lowest), then W+jets groups
    hs, errs, h_extra = [], [], {}
    stack_extra = sorted(EXTRA, key=lambda k: 0 if k == "DY" else 1)
    for key in stack_extra:
        X = EXTRA[key]
        xm = sb(X["Ds_mass"])
        h, _ = np.histogram(X[branch][xm], bins=bins, weights=X["w"][xm])
        e2, _ = np.histogram(X[branch][xm], bins=bins, weights=X["w"][xm] ** 2)
        hs.append(h); errs.append(np.sqrt(e2)); h_extra[key] = h
    hs_wj = []
    for _, modes, _ in GROUPS:
        mm = np.isin(gm, modes) & sb(g_ds)
        h, _ = np.histogram(mc_vals[mm], bins=bins, weights=gw[mm])
        e2, _ = np.histogram(mc_vals[mm], bins=bins, weights=gw[mm] ** 2)
        hs.append(h); errs.append(np.sqrt(e2)); hs_wj.append(h)
    tot_wj = np.sum(hs_wj, axis=0)
    tot = np.sum(hs, axis=0)
    tot_err = np.sqrt(np.sum(np.square(errs), axis=0))
    # drawn heights per width[0] (the wider Ds SR bin is scaled down); sums/ratios unscaled
    sc = width[0] / width
    bottom = np.zeros(len(bc))
    for (lab, _, cl), h in zip([(EXTRA_DEF[k][1], None, EXTRA_DEF[k][2]) for k in stack_extra] + GROUPS, hs):
        ax.bar(bc, h * sc, width=width, bottom=bottom, color=cl, alpha=0.9,
               label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
        bottom = bottom + h * sc
    ax.bar(bc, 2 * tot_err * sc, width=width, bottom=(tot - tot_err) * sc, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
           label="MC stat")

    dh, _ = np.histogram(d_vals, bins=bins)
    assert dh.sum() == len(d_vals), (
        f"{branch}: histogram dropped {len(d_vals) - dh.sum()} events as under/overflow "
        f"-- bin range [{bins[0]:.4f},{bins[-1]:.4f}] does not cover the data.")
    dh_m = np.where(vis, dh * sc, np.nan)
    ax.errorbar(bc, dh_m, yerr=np.sqrt(np.where(vis, dh, 0)) * sc, fmt="o", ms=4,
                color="black", elinewidth=1.2,
                label=f"OS data{' (SR blind)' if blind_half else ''}  N={dh[vis].sum():,.0f}")
    if blind_half:
        ax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15, zorder=0)
    if refline:
        ax.axvline(refline, color="green", ls=":", lw=1.4)
    ax.set_xlim(bins[0], bins[-1])
    peak = max(bottom.max(), np.nanmax(dh_m)) if np.isfinite(np.nanmax(dh_m)) else bottom.max()
    ax.set_ylim(0, YMAX[col] if YMAX else peak * 1.45)
    ax.tick_params(labelbottom=False)
    if col == 0:
        ax.set_ylabel(f"Events / {width[0]*1000:.0f} MeV")
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

    if branch == "Ds_mass":
        # windows on the UNBINNED mass (bin edges no longer sit on 1.830 / 1.910)
        def win(m, lo, hi):
            return (m >= lo) & (m < hi)
        for wname, lo, hi in (("low <1.830", -np.inf, 1.830), ("D+- 1.83-1.91", 1.830, 1.910),
                              ("high >2.008", M_DS + SR_HALF, np.inf)):
            nd = win(d_ds, lo, hi).sum()
            wwj = gw[win(g_ds, lo, hi)]
            wex = {k: X["w"][win(X["Ds_mass"], lo, hi)] for k, X in EXTRA.items()}
            nw = wwj.sum()
            nall = nw + sum(v.sum() for v in wex.values())
            ew = np.sqrt((wwj ** 2).sum() + sum((v ** 2).sum() for v in wex.values()))
            rr = nd / nall
            shares = "  ".join(f"{k} {v.sum():8,.1f} ({100*v.sum()/nall:4.1f}%)"
                               for k, v in wex.items())
            print(f"  {wname:14}: data {nd:7,.0f}  W+jets {nw:9,.1f}  {shares}"
                  f"  data/Wjets {nd/nw:.3f}"
                  f"  data/allMC {rr:.3f}+-{rr*np.hypot(1/np.sqrt(max(nd,1)), ew/nall):.3f}")

fig.suptitle(f"stage-1 {S1_KEY} + stage-2 {KEY2}: "
             + (rf"stage-1 only, WP1 $\varepsilon_1$={WP1:.1%}" if STAGE1_ONLY else
                rf"WP1 $\varepsilon_1$={WP1:.0%}  +  WP2 $\varepsilon_2$={WP2:.0%}")
             + rf"  $\to$  $\varepsilon_{{tot}}$={eps_tot:.1%}" "\n"
             + "MC x BR rw" + (" (DY x jet-bin SF)" if "DY" in EXTRA else ""), fontsize=13)

out = common.fig_path("comparisons", (f"wp_mass_wp1_{WP1:.3f}_stage1only.png" if STAGE1_ONLY else
                        f"wp_mass_wp1_{WP1:.2f}_wp2_{WP2:.2f}")
                      .replace(".png", "") + ("_" + "_".join(k.lower() for k in EXTRA) if EXTRA else "") + ".png",
                      key=KEY2, tool="BDT")
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

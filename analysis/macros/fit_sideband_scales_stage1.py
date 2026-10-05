"""Ds_mass sideband fit of per-component MC scale factors at stage-1 score cuts.

ONE set of factors, common to every stage-1 cut, fitted to OS data in the Ds-mass
sidebands (20 MeV bins, |M-M_Ds|>40 MeV, same blind as everywhere else) with
    data = W+jets(mode 0,6-10) + W+jets(mode 1+2) + TT            [fixed at MC]
           + r34 * W+jets(mode 3+4) + r5 * W+jets(mode 5) + k * DY  [floated]
Only these three float (user's choice, 2026-10-05): W+jets and TT combinatorial are taken
as correct, and mode 1+2 (true Ds) is kept out of the fit because it sits almost entirely
in the blinded SR, so the sidebands cannot constrain it.
The fit is simultaneous over DISJOINT score slices (each event counted once), with
r34, r5, k shared by all slices: eps90-70, eps70-50, >eps50. Nothing below eps90 is
fitted or drawn (user's choice, 2026-10-05: no fit at "no cut"; with W+jets
combinatorial fixed, that 2.25M-event bulk would dominate and pull r34, r5 negative).
Linear chi2, sigma^2 = data + sum of MC sumw2. Selection/thresholds/colours follow
draw_stage1_mass_by_score_cut.py. Figure: the nested cuts eps90/70/50, every
column drawn with the SAME common factors; ratio panel shows prefit (open) and postfit
(filled) data/MC, chi2 per column is with zero fitted parameters in that column.
Printed after the figure: each slice fitted alone, and the Delta chi2 of sharing each
factor vs all-free (how well one common set actually holds).

    conda activate ds_dnn && python3 fit_sideband_scales_stage1.py
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

MODEL_KEY = "ssdata"
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, M_DP, SR_HALF = 1.96835, 1.8696, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
BINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.020)
CUTS = [(0.90, r"$\varepsilon_{sig}$=90%"),
        (0.70, r"$\varepsilon_{sig}$=70%"), (0.50, r"$\varepsilon_{sig}$=50%")]
# (name, legend, colour, fit parameter); stack bottom -> top
COMPS = [
    ("comb", r"W+jets other (mode 0,6-10)", "#5a9e78", None),
    ("TT",   r"$t\bar{t}$",                 "#7f5fbf", None),
    ("DY",   r"DY$\to\ell\ell$",            "#4a90c2", "k"),
    ("5",    r"W+jets D$\to$K$\pi\pi$ misID (5)",   "#e6c229", "r5"),
    ("34",   r"W+jets genuine D$\to$KK$\pi$ (3+4)", "#e07a3c", "r34"),
    ("12",   r"W+jets true D$_s$ (1+2)",            "#c0392b", None),
]
PARS = ["r34", "r5", "k"]
PAR_TEX = {"r34": "r34 (W+jets mode 3+4)", "r5": "r5 (W+jets mode 5)", "k": "k (DY)"}
WMODES = {"comb": [0, 6, 7, 8, 9, 10], "12": [1, 2], "34": [3, 4], "5": [5]}


def rd(name, extra):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["Ds_mass", "score_dsmu"] + extra + list(common.SELECTION), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    return {x: a[x][m] for x in ["Ds_mass", "score_dsmu"] + extra}


W = rd("WtoLNu2Jets_merged", ["weight", "cand_gen_mode"])
T = rd("TT_merged", ["weight"])
D = rd("DY_merged", ["weight"])
X = rd("Data_OS_merged", [])
X = {k: v[np.abs(X["Ds_mass"] - M_DS) >= SR_HALF] for k, v in X.items()}   # blind


def hist(S, m):
    w = S["weight"][m].astype(float) if "weight" in S else np.ones(m.sum())
    return (np.histogram(S["Ds_mass"][m], BINS, weights=w)[0],
            np.histogram(S["Ds_mass"][m], BINS, weights=w ** 2)[0])


bc = 0.5 * (BINS[1:] + BINS[:-1])
width = BINS[1] - BINS[0]
vis = np.abs(bc - M_DS) >= SR_HALF
FIXED = [g for g, _l, _c, p in COMPS if p is None]


def make_hists(lo, hi):
    """Data and every MC component with lo < score_dsmu <= hi."""
    sel = lambda S: (S["score_dsmu"] > lo) & (S["score_dsmu"] <= hi)
    dh = np.histogram(X["Ds_mass"][sel(X)], BINS)[0].astype(float)
    H = {g: hist(W, sel(W) & np.isin(W["cand_gen_mode"], ms)) for g, ms in WMODES.items()}
    H["TT"] = hist(T, sel(T))
    H["DY"] = hist(D, sel(D))
    return dh, H


def fit(regions, shared):
    """Linear chi2 fit on sideband bins of one or more regions [(dh, H), ...].
    Fixed components are subtracted at their MC value; a parameter in `shared` is one
    number for all regions, otherwise each region gets its own copy."""
    ys, sigs, blocks, names = [], [], [], []
    for dh, H in regions:
        sigs.append(np.sqrt(np.maximum(dh + sum(e2 for _h, e2 in H.values()), 1.0))[vis])
        ys.append((dh - sum(H[g][0] for g in FIXED))[vis])
    for p in PARS:
        t = [sum(H[g][0] for g, _l, _c, pp in COMPS if pp == p)[vis] for _d, H in regions]
        if p in shared or len(regions) == 1:
            blocks.append(np.concatenate(t)); names.append(p)
        else:
            for i in range(len(regions)):
                blocks.append(np.concatenate([t[j] if j == i else 0 * t[j] for j in range(len(t))]))
                names.append(f"{p}[{i}]")
    y, sig = np.concatenate(ys), np.concatenate(sigs)
    Xm = np.column_stack(blocks)
    A = Xm / sig[:, None]
    sol = np.linalg.lstsq(A, y / sig, rcond=None)[0]
    cov = np.linalg.inv(A.T @ A)
    err = np.sqrt(np.diag(cov))
    return dict(S=dict(zip(names, sol)), err=dict(zip(names, err)), rho=cov / np.outer(err, err),
                chi2=(((y - Xm @ sol) / sig) ** 2).sum(), ndf=len(y) - len(sol),
                chi2_pre=(((y - Xm @ np.ones(len(sol))) / sig) ** 2).sum(), nbins=len(y))


fig = plt.figure(figsize=(7.2 * len(CUTS), 7.5))
outer = gridspec.GridSpec(1, len(CUTS), figure=fig, wspace=0.30)
THR = {eff: common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))[0]
       for eff, _l in CUTS if eff is not None}

SLICES = [("eps90-70", THR[0.90], THR[0.70]),
          ("eps70-50", THR[0.70], THR[0.50]), ("top (>eps50)", THR[0.50], 2.0)]
regions = [make_hists(lo, hi) for _n, lo, hi in SLICES]
G = fit(regions, PARS)
S, err = G["S"], [G["err"][p] for p in PARS]
rho_txt = "  ".join(f"rho({PARS[i]},{PARS[j]})={G['rho'][i, j]:+.2f}" for i, j in ((0, 1), (0, 2), (1, 2)))
print(f"COMMON FIT ({', '.join(n for n, _lo, _hi in SLICES)}): "
      + "  ".join(f"{p}={S[p]:.3f}+-{e:.3f}" for p, e in zip(PARS, err))
      + f"  chi2/ndf {G['chi2']:.1f}/{G['ndf']}  {rho_txt}")

for col, (eff, lbl) in enumerate(CUTS):
    dh, H = make_hists(THR[eff], 2.0)
    sig = np.sqrt(np.maximum(dh + sum(e2 for _h, e2 in H.values()), 1.0))[vis]
    y = (dh - sum(H[g][0] for g in FIXED))[vis]
    tmpl = {p: sum(H[g][0] for g, _l, _c, pp in COMPS if pp == p)[vis] for p in PARS}
    chi2 = (((y - sum(S[p] * tmpl[p] for p in PARS)) / sig) ** 2).sum()
    chi2_pre = (((y - sum(tmpl.values())) / sig) ** 2).sum()
    nb = vis.sum()
    print(f"{lbl:28}: common factors -> chi2/nbins post {chi2:.1f}/{nb}  pre {chi2_pre:.1f}/{nb}")

    # ---- draw postfit ----
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[col],
                                             height_ratios=[3, 1], hspace=0.05)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    bottom = np.zeros(len(bc))
    post, post_e2, pre = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
    for g, lab, cl, p in COMPS:
        s = 1.0 if p is None else S[p]
        h, e2 = H[g]
        ax.bar(bc, s * h, width=width, bottom=bottom, color=cl, alpha=0.9, edgecolor="none",
               label=f"{lab}" + (f" ×{s:.2f}" if p else "") + f"  N={s*h.sum():,.0f}")
        bottom += s * h
        post += s * h; post_e2 += s * s * e2; pre += h
    pe = np.sqrt(post_e2)
    ax.bar(bc, 2 * pe, width=width, bottom=post - pe, color="none", hatch="////",
           edgecolor="0.35", linewidth=0.0, zorder=3, label="MC stat (postfit)")
    dm = np.where(vis, dh, np.nan)
    ax.errorbar(bc, dm, yerr=np.sqrt(np.where(vis, dh, 0)), fmt="o", ms=4, color="black",
                label=f"OS data (Ds-SR excl.)  N={np.nansum(dm):,.0f}")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
    ax.axvline(M_DP, color="green", ls=":", lw=1.4)
    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_ylim(0, max(bottom.max(), np.nanmax(dm)) * 1.9)
    ax.tick_params(labelbottom=False)
    ax.legend(fontsize=7.5, loc="upper right")
    ax.text(0.03, 0.97, f"$\\chi^2$/nbins = {chi2:.1f}/{nb}\n(prefit {chi2_pre:.1f}/{nb})",
            transform=ax.transAxes, va="top", ha="left", fontsize=9)
    ax.set_title(f"{lbl} — common-factor postfit", fontsize=13, pad=22)
    if col == 0:
        ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
        hep.cms.label(ax=ax, data=True, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

    with np.errstate(invalid="ignore", divide="ignore"):
        r_post = np.where(vis & (post > 0), dh / post, np.nan)
        r_pre = np.where(vis & (pre > 0), dh / pre, np.nan)
        re_post = np.where(vis & (post > 0), np.sqrt(dh) / post, np.nan)
        re_pre = np.where(vis & (pre > 0), np.sqrt(dh) / pre, np.nan)
        band = np.where(vis & (post > 0), pe / post, np.nan)
    rax.bar(bc, 2 * band, width=width, bottom=1 - band, color="none", hatch="////",
            edgecolor="0.35", linewidth=0.0)
    rax.errorbar(bc - 0.002, r_pre, yerr=re_pre, fmt="o", ms=3.5, mfc="white", color="0.45",
                 elinewidth=1.0, label="prefit")
    rax.errorbar(bc + 0.002, r_post, yerr=re_post, fmt="o", ms=3.5, color="black",
                 elinewidth=1.0, label="postfit")
    n_off = int(np.sum((np.nan_to_num(r_pre, nan=1) > 1.99) | (np.nan_to_num(r_pre, nan=1) < 0)
                       | (np.nan_to_num(r_post, nan=1) > 1.99) | (np.nan_to_num(r_post, nan=1) < 0)))
    if n_off:
        rax.text(0.02, 0.92, f"{n_off} point(s) off scale", transform=rax.transAxes,
                 fontsize=8, va="top")
    rax.axhline(1.0, color="red", lw=1.0)
    rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    rax.set_ylim(0.0, 1.99)
    rax.set_yticks([0.5, 1.0, 1.5])
    rax.set_xlim(BINS[0], BINS[-1])
    rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    rax.legend(fontsize=8, loc="lower right", ncol=2)
    if col == 0:
        rax.set_ylabel("data / MC")

fig.suptitle(r"Sideband fit of MC scale factors, stage-1 (ssdata DNN) cuts — $M(K^+K^-\pi^\pm)$"
             + "\ncommon to all columns, fitted on disjoint slices eps90-70 / 70-50 / >50:  "
             + ",  ".join(f"{PAR_TEX[p]} = {S[p]:.2f} ± {e:.2f}" for p, e in zip(PARS, err))
             + f"  ($\\chi^2$/ndf = {G['chi2']:.1f}/{G['ndf']})",
             fontsize=14, y=1.06)
out = common.fig_path("closure", "stage1_sideband_scale_fit_Ds_mass.png", key=MODEL_KEY)
plt.savefig(out, dpi=120, bbox_inches="tight")
print(f"Saved {out}")


# ---- are the factors the same at every cut? disjoint score slices ----
from scipy.stats import chi2 as chi2_dist
print("\n== disjoint score slices, each fitted alone ==")
for (name, _lo, _hi), reg in zip(SLICES, regions):
    F = fit([reg], PARS)
    print(f"  {name:15} data {reg[0][vis].sum():9,.0f}  "
          + "  ".join(f"{p}={F['S'][p]:.2f}+-{F['err'][p]:.2f}" for p in PARS)
          + f"  chi2/ndf {F['chi2']:.1f}/{F['ndf']}")
for label, regs in (("3 slices above eps90", regions),):
    free = fit(regs, [])
    print(f"== constancy, {label}: all-free chi2/ndf {free['chi2']:.1f}/{free['ndf']} ==")
    for sh in (["r34"], ["r5"], ["k"], PARS):
        F = fit(regs, sh)
        dc, dn = F["chi2"] - free["chi2"], F["ndf"] - free["ndf"]
        print(f"  share {'+'.join(sh):10}: dchi2={dc:6.1f}/{dn}  p={chi2_dist.sf(dc, dn):.3g}   "
              + "  ".join(f"{p}={F['S'][p]:.2f}+-{F['err'][p]:.2f}" for p in sh))

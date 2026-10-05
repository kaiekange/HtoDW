"""DY alone: Ds_mass across STAGE-1-ONLY DNN score cuts, stacked by GenTag group.

Stripped-down copy of draw_stage1_mass_by_score_cut.py (same ossel(), stage-1
thresholds, CUTS, GROUPS, colours, 10 MeV bins) with the data and W+jets parts
removed -- so each column here is exactly the DY slice of that script's --add DY
stack. Absolute sigma*L yields, not normalized (feedback_plot_style). DY is MC:
nothing is blinded; the grey SR band and the D+- line are reference marks only.

DY (amcatnloFXFX) has ~29% negative-weight rows, so per column it prints raw rows,
sum w, N_eff=(sum w)^2/sum w^2, negative-bin counts per group, and a window table
(low <1.830 / D+- 1.83-1.91 / high >2.008) with the low/high ratio -- the tilt.

    conda activate ds_dnn && python3 draw_stage1_dy_mass_by_cut.py
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

MODEL_KEY = "ssdata"
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, M_DP, SR_HALF = 1.96835, 1.8696, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
SEL = [*common.SELECTION]
BINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)

CUTS = [
    (None, "no cut"),
    (0.90, r"$\varepsilon_{sig}$=90%"),
    (0.70, r"$\varepsilon_{sig}$=70%"),
    (0.50, r"$\varepsilon_{sig}$=50%"),
]

GROUPS = [   # stack bottom -> top, same order/colours as draw_stage1_mass_by_score_cut.py
    (r"other (mode 0,6-10)",         [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",   [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],              "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],              "#c0392b"),
]

THR = {eff: common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))[0]
       for eff, _l in CUTS if eff is not None}

a = uproot.open(f"{PRED}/DY_merged.root")["Events"].arrays(
    ["Ds_mass", "score_dsmu", "cand_gen_mode", "weight"] + SEL, library="np")
k = common.selection_mask(a) & (a["score_dsmu"] > -1)
print(f"DY: {len(k):,} rows -> {k.sum():,} after OS+b-veto+score>-1")
ds, sc, mode, w = a["Ds_mass"][k], a["score_dsmu"][k], a["cand_gen_mode"][k], a["weight"][k].astype(float)

n_out = len(ds) - ((ds >= BINS[0]) & (ds <= BINS[-1])).sum()
assert n_out <= 5, f"{n_out} DY rows outside [{BINS[0]:.4f},{BINS[-1]:.4f}] -- would be dropped"
grp_masks = [np.isin(mode, ms) for _l, ms, _c in GROUPS]
assert np.logical_or.reduce(grp_masks).all(), "cand_gen_mode outside the four groups"


def hist(v, ww):
    return (np.histogram(v, bins=BINS, weights=ww)[0],
            np.sqrt(np.histogram(v, bins=BINS, weights=ww ** 2)[0]))


bc = 0.5 * (BINS[1:] + BINS[:-1])
width = BINS[1] - BINS[0]
WINDOWS = [("low <1.830", bc < 1.830), ("D+- 1.83-1.91", (bc >= 1.830) & (bc < 1.910)),
           ("high >2.008", bc > M_DS + SR_HALF)]

fig, axes = plt.subplots(1, len(CUTS), figsize=(7.2 * len(CUTS), 6.5))
for col, ((eff, lbl), ax) in enumerate(zip(CUTS, axes)):
    cut = np.ones(len(ds), dtype=bool) if eff is None else (sc > THR[eff])
    hs, errs = [], []
    for (lab, _ms, _c), gm in zip(GROUPS, grp_masks):
        m = gm & cut
        h, e = hist(ds[m], w[m])
        if (h < 0).any():
            print(f"  [{lbl} {lab}] {(h < 0).sum()} bin(s) negative, min={h.min():.3g}")
        hs.append(h); errs.append(e)
    tot = np.sum(hs, axis=0)
    tot_err = np.sqrt(np.sum(np.square(errs), axis=0))

    bottom = np.zeros(len(bc))
    for (lab, _ms, cl), h in zip(GROUPS, hs):
        ax.bar(bc, h, width=width, bottom=bottom, color=cl, alpha=0.9,
               label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
        bottom = bottom + h
    ax.bar(bc, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3, label="MC stat")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
    ax.axvline(M_DP, color="green", ls=":", lw=1.4)
    ax.set_xlim(BINS[0], BINS[-1])
    peak = max(bottom.max(), (tot + tot_err).max())
    ax.set_ylim(min(0, (tot - tot_err).min() * 1.1), peak * 1.5 if peak > 0 else 1.0)
    ax.text(M_DP, ax.get_ylim()[1] * 0.97, r"D$^\pm$", ha="center", va="top", fontsize=11)
    ax.text(M_DS, ax.get_ylim()[1] * 0.97, r"D$_s^\pm$", ha="center", va="top", fontsize=11)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
        hep.cms.label(ax=ax, data=False, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)
    ax.legend(fontsize=8, loc="upper right", title=r"DY$\to\ell\ell$ (0J+1J+2J)", title_fontsize=9)
    ax.set_title(lbl, fontsize=13, pad=6)

    sw, sw2 = w[cut].sum(), (w[cut] ** 2).sum()
    print(f"{lbl:28}: rows {cut.sum():9,}  sumw {sw:10,.1f}  N_eff {sw**2/sw2:9,.1f}  "
          f"neg-w rows {(w[cut] < 0).mean():.3f}")
    vals = {}
    for wname, wsel in WINDOWS:
        n, e = tot[wsel].sum(), np.sqrt((tot_err[wsel] ** 2).sum())
        vals[wname] = (n, e)
        comb = hs[0][wsel].sum() / n if n else np.nan
        print(f"    {wname:14}: {n:9,.1f} +- {e:7,.1f}  (mode 0,6-10 share {comb:.2f})")
    (nl, el), (nh, eh) = vals["low <1.830"], vals["high >2.008"]
    print(f"    low/high      : {nl/nh:.3f} +- {nl/nh*np.hypot(el/nl, eh/nh):.3f}")

fig.suptitle(r"Stage-1 (ssdata DNN) score cut scan -- DY, $M(K^+K^-\pi^\pm)$", fontsize=15)
out = common.fig_path("closure", "stage1_dy_Ds_mass_by_cut.png", key=MODEL_KEY)
plt.savefig(out, dpi=120, bbox_inches="tight")
print(f"Saved {out}")

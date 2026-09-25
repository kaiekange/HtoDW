"""Ds_mass only, OS data vs W+jets, split preEE (top row) / postEE (bottom row),
scanned across the same stage-1-only DNN score cuts as draw_stage1_mass_by_score_cut.py
(no cut / eps_sig 90/70/50%). W+jets stacked into the four GenTag truth groups.

Same event selection, blinding, threshold and negative-weight conventions as
draw_stage1_mass_by_score_cut.py's Ds_mass panel -- see that script's docstring for
the six points this follows (base score>-1 mask, Ds-SR event exclusion for data,
common.stage1_threshold, group-partition assert, overflow assert). The only addition
here is the `era` branch split (1=preEE, 2=postEE, analysis/eras.py).

Thresholds are computed ONCE from the combined (both-era) truth-matched signal and
reused for both rows -- a per-era quantile would silently make the two rows different
cuts (compare_gentag_eras.py makes the same point). Each row's CMS label carries only
that era's own luminosity (not the combined 34.7 fb^-1).

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_stage1_dsmass_by_era.py
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
M_DS, SR_HALF = 1.96835, 0.040
BINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)
BC = 0.5 * (BINS[1:] + BINS[:-1])
WIDTH = BINS[1] - BINS[0]
XLABEL = r"$M(K^+K^-\pi^\pm)$ [GeV]"
SEL = [*common.SELECTION]  # num_os_mu_tight, num_ss_mu_tight, n_bjet_corr

CUTS = [
    (None, "no cut"),
    (0.90, r"$\varepsilon_{sig}$=90%"),
    (0.70, r"$\varepsilon_{sig}$=70%"),
    (0.50, r"$\varepsilon_{sig}$=50%"),
]

ERA_ROWS = [(1, "preEE"), (2, "postEE")]   # analysis/eras.py ERAS

GROUPS = [   # stack bottom -> top, same order/colours as draw_stage1_mass_by_score_cut.py
    (r"other (mode 0,6-10)",         [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",   [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],              "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],              "#c0392b"),
]


def ossel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a) & (a["score_dsmu"] > -1)


def rd(name, extra):
    return uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["Ds_mass", "score_dsmu", "era"] + SEL + extra, library="np")


# ---- stage-1 thresholds: combined-era signal, same recipe as the other stage-1 macro ----
THR = {}
for eff, _lbl in CUTS:
    if eff is None:
        continue
    thr, _n = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))
    THR[eff] = thr

# ---- W+jets GenTag: base selection (OS + b-veto + score>-1), full mass range ----
w = rd("WtoLNu2Jets_merged", ["cand_gen_mode", "weight"])
wk = ossel(w)
print(f"W+jets: {len(w['Ds_mass']):,} rows -> {wk.sum():,} after OS+b-veto+score>-1")
w_ds = w["Ds_mass"][wk]
w_mode, w_score, w_w = w["cand_gen_mode"][wk], w["score_dsmu"][wk], w["weight"][wk].astype(float)
w_era = w["era"][wk]

n_out = len(w_ds) - ((w_ds >= BINS[0]) & (w_ds <= BINS[-1])).sum()
assert n_out <= max(5, int(1e-5 * len(w_ds))), (
    f"W+jets Ds_mass: {n_out} selected rows fall outside "
    f"[{BINS[0]:.4f},{BINS[-1]:.4f}] -- would be silently dropped as overflow.")
if n_out:
    print(f"  W+jets Ds_mass: {n_out} row(s) just outside the bin range "
          f"(float-precision edge of the post-fit mass cut) -- dropped from the hist")

grp_masks = [np.isin(w_mode, ms) for _l, ms, _c in GROUPS]
covered = np.logical_or.reduce(grp_masks)
assert covered.all(), (
    f"{(~covered).sum()} W+jets rows have cand_gen_mode outside the four groups "
    f"({sorted(np.unique(w_mode[~covered]).tolist())}) -- would silently vanish "
    f"from the stack.")

# ---- OS data: base selection, Ds-mass SR excluded ONCE ----
d = rd("Data_OS_merged", [])
dk = ossel(d) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALF)
print(f"OS data: {len(d['Ds_mass']):,} rows -> {dk.sum():,} after OS+b-veto+score>-1"
      f"+Ds-SR-excluded")
d_ds, d_score, d_era = d["Ds_mass"][dk], d["score_dsmu"][dk], d["era"][dk]

n_out = len(d_ds) - ((d_ds >= BINS[0]) & (d_ds <= BINS[-1])).sum()
assert n_out <= max(5, int(1e-5 * len(d_ds))), (
    f"OS data Ds_mass: {n_out} selected rows fall outside "
    f"[{BINS[0]:.4f},{BINS[-1]:.4f}] -- would be silently dropped as overflow.")
if n_out:
    print(f"  OS data Ds_mass: {n_out} row(s) just outside the bin range "
          f"(float-precision edge of the post-fit mass cut) -- dropped from the hist")

blind_bins = np.abs(BC - M_DS) < SR_HALF
vis = ~blind_bins


def hist(vals, wt=None):
    h, _ = np.histogram(vals, bins=BINS, weights=wt)
    e2, _ = np.histogram(vals, bins=BINS, weights=(wt ** 2) if wt is not None else None)
    return h, np.sqrt(e2)


nrows, ncols = len(ERA_ROWS), len(CUTS)
fig = plt.figure(figsize=(7.2 * ncols, 7.0 * nrows))
outer = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0.30, hspace=0.35)

for row, (era_id, era_lbl) in enumerate(ERA_ROWS):
    era_lumi = common.ERAS[era_id]["lumi_pb"]
    w_era_mask = (w_era == era_id)
    d_era_mask = (d_era == era_id)

    for col, (eff, lbl) in enumerate(CUTS):
        cutmask = w_era_mask if eff is None else (w_era_mask & (w_score > THR[eff]))
        dcutmask = d_era_mask if eff is None else (d_era_mask & (d_score > THR[eff]))

        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

        hs, errs = [], []
        for (glab, _ms, _c), gm in zip(GROUPS, grp_masks):
            m = gm & cutmask
            h, e = hist(w_ds[m], w_w[m])
            neg = (h < 0).sum()
            if neg:
                print(f"  [{era_lbl} col={lbl!r} group={glab!r}] {neg} bin(s) sum "
                      f"negative (amcatnloFXFX sign) -- min={h.min():.3g}")
            hs.append(h); errs.append(e)
        tot = np.sum(hs, axis=0)
        tot_err = np.sqrt(np.sum(np.square(errs), axis=0))

        bottom = np.zeros(len(BC))
        for (glab, _ms, cl), h in zip(GROUPS, hs):
            ax.bar(BC, h, width=WIDTH, bottom=bottom, color=cl, alpha=0.9,
                   label=f"{glab}  N={h.sum():,.0f}", edgecolor="none")
            bottom = bottom + h
        ax.bar(BC, 2 * tot_err, width=WIDTH, bottom=tot - tot_err, color="none",
               hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
               label="W+jets MC stat")
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)

        dh, _ = np.histogram(d_ds[dcutmask], bins=BINS)
        n_drop = dcutmask.sum() - dh.sum()
        assert n_drop <= 5, (
            f"{era_lbl} col={lbl!r}: histogram dropped {n_drop} OS data events "
            f"as under/overflow.")
        dh_m = np.where(vis, dh, np.nan)
        ax.errorbar(BC, dh_m, yerr=np.sqrt(np.where(vis, dh, 0)), fmt="o", ms=4,
                    color="black", elinewidth=1.2,
                    label=f"OS data (Ds-SR excl.)  N={dh.sum():,.0f}")

        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
        ax.set_xlim(BINS[0], BINS[-1])
        peak = np.nanmax([bottom.max(), np.nanmax(dh_m)])
        ax.set_ylim(0, peak * 1.5 if peak > 0 else 1.0)
        ax.tick_params(labelbottom=False)
        if col == 0:
            ax.set_ylabel(f"Events / {WIDTH*1000:.0f} MeV")
        ax.legend(fontsize=7.5, loc="upper right")
        ax.set_title(f"{era_lbl}  --  {lbl}", fontsize=12, pad=6)
        if col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary",
                          rlabel=rf"{era_lumi/1000:.2f} fb$^{{-1}}$ (13.6 TeV)",
                          fontsize=9)

        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(vis & (tot > 0), dh / tot, np.nan)
            re = np.where(vis & (tot > 0), np.sqrt(np.where(vis, dh, 0)) / tot, np.nan)
            mc_band = np.where(vis & (tot > 0), tot_err / tot, np.nan)
        rax.bar(BC, 2 * mc_band, width=WIDTH, bottom=1 - mc_band, color="none",
                hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(BC, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        rax.set_ylim(0.0, 1.99)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.set_xlim(BINS[0], BINS[-1])
        rax.set_xlabel(XLABEL)
        if col == 0:
            rax.set_ylabel("data / MC")

        print(f"{era_lbl} col={lbl!r}: W+jets N={tot.sum():,.1f} (n={cutmask.sum():,})  "
              f"OS data N={dh.sum():,.0f} (n={dcutmask.sum():,})")

fig.suptitle(r"Stage-1 (ssdata DNN) score cut scan -- $M(K^+K^-\pi^\pm)$, preEE vs postEE",
             fontsize=15)

out = common.fig_path("era", "stage1_dsmass_by_cut_preEE_postEE.png", key=MODEL_KEY)
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

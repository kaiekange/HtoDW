"""Ds_mass and phi_mass, OS data vs W+jets, scanned across STAGE-1-ONLY DNN score cuts.

Stage-1 (ssdata) only -- no stage-2 resonant classifier is applied here, unlike
draw_wp_mass.py. Columns are the score_dsmu working point: no cut, then
eps_sig = 90/70/50% (the project-standard scan points, same as
draw_gentag_mass_modes.py / fom_scan_stage1.py). Thresholds come from
common.stage1_threshold(), the single recipe already used to define stage-2's
training population -- reusing it (rather than hand-rolling the quantile) guarantees
this script's 90% column is the same 90% everywhere else in the project.

W+jets is split into the four GenTag truth groups (cand_gen_mode) and STACKED, same
grouping/colours/order as draw_wp_mass.py / draw_gentag_mass_modes.py:
    mode 1+2   true Ds
    mode 3+4   genuine D->KKpi
    mode 5     D->Kpipi misID reflection
    mode 0,6-10  everything else
Two separate output figures (one per mass variable), each a WP1..4-column
main+ratio grid, following draw_wp_mass.py's per-column layout.

Reads from ssdata's OWN Prediction/ dir (Prediction/{Signal,WtoLNu2Jets,Data_OS}_merged
.root) -- predict.py copies every input branch through, so cand_gen_mode, Ds_mass,
phi_mass and weight all live in the same file as score_dsmu; no join needed.

SIX THINGS THIS SCRIPT GETS RIGHT THAT A NAIVE COPY OF ITS TEMPLATES WOULD MISS
(recorded here because two of the templates this was built from -- draw_wp_mass.py
and the uncommitted draw_ds_mass_signal_data_wjets.py -- do NOT all get them right):

1. Base selection includes score_dsmu > -1 (the non-finite-feature sentinel) even in
   the "no cut" column. draw_ds_mass_signal_data_wjets.py's ossel() omits this (it
   reads SAMPLES_DIR, which has no score branch) -- copying it verbatim here would
   make "no cut" and "eps_sig=90%" disagree on which EVENT SET they start from, not
   just on the score threshold (feedback_compare_like_for_like).

2. OS data has its Ds-mass-SR events (|Ds_mass-M_Ds|<40 MeV) excluded ONCE, up front,
   from the event set used for BOTH panels (draw_wp_mass.py's pattern) -- not via a
   per-panel "blind_bins on the x-axis" trick, which only works when the blind
   variable and the plotted variable are the same. On the phi_mass panel the SR cut
   is event-level, not phi-bin-level, so it removes a roughly flat few-percent slice
   of every phi bin, not a contiguous mass window.

3. W+jets is NEVER SR-excluded for the drawn STACK (it is MC, not blinded) -- but for
   the phi_mass RATIO panel specifically, the denominator uses a second, SR-excluded
   W+jets total so the ratio compares like-for-like event sets (same reasoning as
   #2: data's phi spectrum there is missing a flat SR-shaped slice, so the ratio
   needs the same slice missing from W+jets, or every bin reads systematically low
   by the SR event fraction and the panel is not a closure test). The Ds_mass panel
   needs no such fix: its ratio is only ever read at non-blind bins.

4. weight carries the amcatnloFXFX sign (project_genweight_missing memory). A sparse
   group x cut combination (mode 5 or mode 3+4 at eps_sig=50%) can sum negative in a
   bin; per-group negative-bin counts are printed and a zero line is drawn
   (draw_gentag_mass_modes.py convention) rather than clipping, which would silently
   break "stack sums to the total".

5. The four GenTag groups are asserted to partition the selected W+jets rows --
   catches a cand_gen_mode value outside {0..10} silently vanishing from the stack.

6. Both mass windows are asserted to cover 100% of the selected rows (no silent
   under/overflow drop), for OS data and W+jets, before any cut column is applied.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_stage1_mass_by_score_cut.py [--add TT,DY] [--key=ssdata_phi20]

--key=K picks the stage-1 model (default ssdata): its Prediction/ dir, its thresholds and its
figures dir. If K's training config sets phi_win, the same |phi_mass - M_phi| window is
applied to every sample here (ossel), so the plotted population is the one K was trained on.

--add stacks extra background MC samples (each ONE entry, its own sigma*L weight) on top
of the W+jets GenTag groups: TT = TT_merged (TTto2L2Nu + TTtoLNu2Q), DY = DY_merged
(DYto2L2Jets 0J/1J/2J). The ratio denominator becomes the full MC sum, and a window table
(low <1.830 / D+- 1.83-1.91 / high >2.008) prints each sample's share and data/MC with
W+jets alone and with everything. Output names get a _wjets_<tt_dy...> suffix so the
W+jets-only figures survive. (--tt is kept as an alias for --add TT.)

--mc-sr-excl also removes the Ds-mass SR (|Ds_mass-M_Ds|<40 MeV) from every MC sample at
the base selection, so the drawn stack (not only the phi_mass ratio denominator, point 3)
covers the same event set as the OS data points. Output names get a _mcsrexcl suffix.
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

EXTRA_DEF = {   # key -> (Prediction file, legend label, colour)
    "TT": ("TT_merged", r"$t\bar{t}$ (2L2Nu + LNu2Q)", "#7f5fbf"),
    "DY": ("DY_merged", r"DY$\to\ell\ell$ (0J+1J+2J)", "#4a90c2"),
}
ADD = []
if "--tt" in sys.argv:
    ADD = ["TT"]
for i, a in enumerate(sys.argv):
    if a.startswith("--add"):
        val = a.split("=", 1)[1] if "=" in a else sys.argv[i + 1]
        ADD = [x for x in val.split(",") if x]
bad = [k for k in ADD if k not in EXTRA_DEF]
if bad:
    sys.exit(f"--add: unknown {bad}; choices {list(EXTRA_DEF)}")
MC_SR_EXCL = "--mc-sr-excl" in sys.argv
MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata")
PHI_WIN = common.resolve(MODEL_KEY).get("phi_win")
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, M_PHI, SR_HALF = 1.96835, 1.019460, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
SEL = [*common.SELECTION]  # corrected muon cuts + charge branches + n_bjet_corr

CUTS = [
    (None, "no cut"),
    (0.90, r"$\varepsilon_{sig}$=90%"),
    (0.70, r"$\varepsilon_{sig}$=70%"),
    (0.50, r"$\varepsilon_{sig}$=50%"),
]

GROUPS = [   # stack bottom -> top, same order/colours as draw_wp_mass.py
    (r"other (mode 0,6-10)",         [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",   [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],              "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],              "#c0392b"),
]

VARS = [   # (branch, xlabel, bin edges, blind half-width or None, ref line or None)
    ("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]",
     np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010), SR_HALF, None),
    ("phi_mass", r"$M(K^+K^-)$ [GeV]",
     np.arange(0.984, 1.074 + 1e-9, 0.002), None, M_PHI),
]


def ossel(a):
    # common.selection_mask, NOT a hand-written copy. The analysis selection changed on
    # 2026-09-21: the old `num_os_mu_tight == 1 & num_ss_mu_tight == 0` was not a muon
    # selection at all (the producer counts every tight muon with no pt/eta/iso cut, so
    # it acted as an uncalibrated soft-in-jet-muon tag rejecting 7.6% of OS data against
    # 2.0% of signal) and was replaced by the corrected-quantity cuts already baked into
    # NanoMVA/Samples/. Keeping the old copy here would draw this figure on a DIFFERENT
    # event set than the model on the x-axis was trained and scored with -- exactly the
    # like-for-like failure point 1 of this docstring warns about, one level up.
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if PHI_WIN is not None:
        m &= np.abs(a["phi_mass"] - M_PHI) < PHI_WIN
    return m


def rd(name, extra):
    return uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["Ds_mass", "phi_mass", "score_dsmu"] + SEL + extra, library="np")


# ---- stage-1 thresholds: same recipe as stage-2's training population ----
THR = {}
for eff, _lbl in CUTS:
    if eff is None:
        continue
    thr, _n = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))
    THR[eff] = thr

# ---- W+jets GenTag: base selection (OS + b-veto + score>-1), full mass range ----
w = rd("WtoLNu2Jets_merged", ["cand_gen_mode", "weight"])
wk = ossel(w)
if MC_SR_EXCL:
    wk &= np.abs(w["Ds_mass"] - M_DS) >= SR_HALF
print(f"W+jets: {len(w['Ds_mass']):,} rows -> {wk.sum():,} after OS+b-veto+score>-1"
      + ("+Ds-SR-excluded" if MC_SR_EXCL else ""))
w_ds, w_phi = w["Ds_mass"][wk], w["phi_mass"][wk]
w_mode, w_score, w_w = w["cand_gen_mode"][wk], w["score_dsmu"][wk], w["weight"][wk].astype(float)
w_not_sr = np.abs(w_ds - M_DS) >= SR_HALF

for branch, vals in (("Ds_mass", w_ds), ("phi_mass", w_phi)):
    bins = next(b for br, _l, b, _bh, _r in VARS if br == branch)
    n_out = len(vals) - ((vals >= bins[0]) & (vals <= bins[-1])).sum()
    # A handful of rows can sit float-precision-past the post-fit |M-M_Ds|<0.3 cut
    # boundary the tuplizer applied; anything beyond that (a bin-range typo) is loud.
    assert n_out <= max(5, int(1e-5 * len(vals))), (
        f"W+jets {branch}: {n_out} selected rows fall outside "
        f"[{bins[0]:.4f},{bins[-1]:.4f}] -- would be silently dropped as overflow.")
    if n_out:
        print(f"  W+jets {branch}: {n_out} row(s) just outside the bin range "
              f"(float-precision edge of the post-fit mass cut) -- dropped from the hist")

grp_masks = [np.isin(w_mode, ms) for _l, ms, _c in GROUPS]
covered = np.logical_or.reduce(grp_masks)
assert covered.all(), (
    f"{(~covered).sum()} W+jets rows have cand_gen_mode outside the four groups "
    f"({sorted(np.unique(w_mode[~covered]).tolist())}) -- would silently vanish "
    f"from the stack.")

# ---- extra background MC (--add): same base selection, full mass range, never SR-excluded ----
EXTRA = {}
for key in ADD:
    fname, _lab, _col = EXTRA_DEF[key]
    x = rd(fname, ["weight"])
    xk = ossel(x)
    if MC_SR_EXCL:
        xk &= np.abs(x["Ds_mass"] - M_DS) >= SR_HALF
    print(f"{key}: {len(x['Ds_mass']):,} rows -> {xk.sum():,} after OS+b-veto+score>-1")
    EXTRA[key] = dict(Ds_mass=x["Ds_mass"][xk], phi_mass=x["phi_mass"][xk],
                      score=x["score_dsmu"][xk], w=x["weight"][xk].astype(float))
    EXTRA[key]["not_sr"] = np.abs(EXTRA[key]["Ds_mass"] - M_DS) >= SR_HALF

# ---- OS data: base selection, Ds-mass SR excluded ONCE for both panels (point 2) ----
d = rd("Data_OS_merged", [])
dk = ossel(d) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALF)
print(f"OS data: {len(d['Ds_mass']):,} rows -> {dk.sum():,} after OS+b-veto+score>-1"
      f"+Ds-SR-excluded")
d_ds, d_phi, d_score = d["Ds_mass"][dk], d["phi_mass"][dk], d["score_dsmu"][dk]

for branch, vals in (("Ds_mass", d_ds), ("phi_mass", d_phi)):
    bins = next(b for br, _l, b, _bh, _r in VARS if br == branch)
    n_out = len(vals) - ((vals >= bins[0]) & (vals <= bins[-1])).sum()
    assert n_out <= max(5, int(1e-5 * len(vals))), (
        f"OS data {branch}: {n_out} selected rows fall outside "
        f"[{bins[0]:.4f},{bins[-1]:.4f}] -- would be silently dropped as overflow.")
    if n_out:
        print(f"  OS data {branch}: {n_out} row(s) just outside the bin range "
              f"(float-precision edge of the post-fit mass cut) -- dropped from the hist")


def hist(vals, w=None, bins=None):
    h, _ = np.histogram(vals, bins=bins, weights=w)
    e2, _ = np.histogram(vals, bins=bins, weights=(w ** 2) if w is not None else None)
    return h, np.sqrt(e2)


for branch, xlabel, bins, blind_half, refline in VARS:
    bc = 0.5 * (bins[1:] + bins[:-1])
    width = bins[1] - bins[0]
    blind_bins = (np.abs(bc - M_DS) < blind_half) if blind_half else np.zeros_like(bc, dtype=bool)
    vis = ~blind_bins

    w_vals = w_ds if branch == "Ds_mass" else w_phi
    d_vals = d_ds if branch == "Ds_mass" else d_phi

    ncols = len(CUTS)
    fig = plt.figure(figsize=(7.2 * ncols, 7.5))
    outer = gridspec.GridSpec(1, ncols, figure=fig, wspace=0.30)

    for col, (eff, lbl) in enumerate(CUTS):
        cutmask = np.ones(len(w_vals), dtype=bool) if eff is None else (w_score > THR[eff])
        dcutmask = np.ones(len(d_vals), dtype=bool) if eff is None else (d_score > THR[eff])

        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

        # ---- stacked W+jets, full range (never SR-excluded -- point 3) ----
        hs, errs = [], []
        for (lab, _ms, _c), gm in zip(GROUPS, grp_masks):
            m = gm & cutmask
            h, e = hist(w_vals[m], w_w[m], bins)
            neg = (h < 0).sum()
            if neg:
                print(f"  [{branch} col={lbl!r} group={lab!r}] {neg} bin(s) sum negative "
                      f"(amcatnloFXFX sign, point 4) -- min={h.min():.3g}")
            hs.append(h); errs.append(e)
        tot_wj = np.sum(hs, axis=0)
        h_extra, xcuts = {}, {}
        for key, X in EXTRA.items():
            xcuts[key] = np.ones(len(X["w"]), dtype=bool) if eff is None else (X["score"] > THR[eff])
            h, e = hist(X[branch][xcuts[key]], X["w"][xcuts[key]], bins)
            hs.append(h); errs.append(e); h_extra[key] = h
        tot = np.sum(hs, axis=0)
        tot_err = np.sqrt(np.sum(np.square(errs), axis=0))

        bottom = np.zeros(len(bc))
        stack = GROUPS + [(EXTRA_DEF[k][1], None, EXTRA_DEF[k][2]) for k in EXTRA]
        for (lab, _ms, cl), h in zip(stack, hs):
            ax.bar(bc, h, width=width, bottom=bottom, color=cl, alpha=0.9,
                   label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
            bottom = bottom + h
        ax.bar(bc, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
               hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
               label="MC stat")
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)

        # ---- ratio denominator: SR-excluded W+jets total, phi_mass panel only
        # (point 3) -- Ds_mass panel reuses `tot`, read only at vis bins ----
        if blind_half is None:
            m_ratio = cutmask & w_not_sr
            tot_ratio, _ = hist(w_vals[m_ratio], w_w[m_ratio], bins)
            for key, X in EXTRA.items():
                mt = xcuts[key] & X["not_sr"]
                tot_ratio = tot_ratio + hist(X[branch][mt], X["w"][mt], bins)[0]
        else:
            tot_ratio = tot

        dh, _ = np.histogram(d_vals[dcutmask], bins=bins)
        n_drop = dcutmask.sum() - dh.sum()
        assert n_drop <= 5, (
            f"{branch} col={lbl!r}: histogram dropped "
            f"{n_drop} OS data events as under/overflow.")
        dh_m = np.where(vis, dh, np.nan) if blind_half else dh.astype(float)
        ax.errorbar(bc, dh_m, yerr=np.sqrt(np.where(vis, dh, 0)) if blind_half else np.sqrt(dh),
                    fmt="o", ms=4, color="black", elinewidth=1.2,
                    label=f"OS data (Ds-SR excl.)  N={dh.sum():,.0f}")

        if blind_half:
            ax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15, zorder=0)
        if refline:
            ax.axvline(refline, color="green", ls=":", lw=1.4)

        ax.set_xlim(bins[0], bins[-1])
        peak = np.nanmax([bottom.max(), np.nanmax(dh_m)])
        ax.set_ylim(0, peak * 1.5 if peak > 0 else 1.0)
        ax.tick_params(labelbottom=False)
        if col == 0:
            ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_title(f"{lbl}", fontsize=13, pad=6)
        if col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary",
                          rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(vis & (tot_ratio > 0), dh / tot_ratio, np.nan)
            re = np.where(vis & (tot_ratio > 0),
                          np.sqrt(np.where(vis, dh, 0)) / tot_ratio, np.nan)
            mc_band = np.where(vis & (tot_ratio > 0), tot_err / tot_ratio, np.nan)
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

        print(f"{branch} col={lbl!r}: W+jets N={tot_wj.sum():,.1f} (n={cutmask.sum():,})  "
              + "".join(f"{k} N={h.sum():,.1f}  " for k, h in h_extra.items())
              + f"OS data N={dh.sum():,.0f} (n={dcutmask.sum():,})")
        if branch == "Ds_mass" and EXTRA:
            for wname, wsel in (("low <1.830", bc < 1.830),
                                ("D+- 1.83-1.91", (bc >= 1.830) & (bc < 1.910)),
                                ("high >2.008", bc > M_DS + SR_HALF)):
                nd, nw, nall = dh[wsel].sum(), tot_wj[wsel].sum(), tot[wsel].sum()
                ew = np.sqrt((tot_err[wsel] ** 2).sum())
                rr = nd / nall
                shares = "  ".join(f"{k} {h[wsel].sum():9,.1f} ({100*h[wsel].sum()/nall:4.1f}%)"
                                   for k, h in h_extra.items())
                print(f"    {wname:14}: data {nd:9,.0f}  W+jets {nw:11,.1f}  {shares}"
                      f"  data/Wjets {nd/nw:.3f}"
                      f"  data/allMC {rr:.3f}+-{rr*np.hypot(1/np.sqrt(max(nd,1)), ew/nall):.3f}")

    var_tex = r"$M(K^+K^-\pi^\pm)$" if branch == "Ds_mass" else r"$M(K^+K^-)$"
    fig.suptitle(f"Stage-1 (ssdata DNN) score cut scan -- {var_tex}"
                 + ("".join(f" + {k}" for k in EXTRA) and " -- W+jets" + "".join(f" + {k}" for k in EXTRA)),
                 fontsize=15)

    out = common.fig_path("closure", f"stage1_mass_{branch}_by_cut"
                          + ("_wjets_" + "_".join(k.lower() for k in EXTRA) if EXTRA else "")
                          + ("_mcsrexcl" if MC_SR_EXCL else "")
                          + ".png", key=MODEL_KEY)
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

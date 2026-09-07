"""Inclusive W+jets against OS data in M(KKpi), per era, at matched signal efficiency.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_incl_by_era.py            # DNN:ssdata, seed 42
    python3 compare_incl_by_era.py BDT:ssdata

Not folded into check_data_vs_wjets.py: that script is wired to the jet-binned
W+jets set and Data_2022 alone, and produces eight figures. Adding an era axis and
an inclusive/jet-binned axis would change all eight, so this asks the one question
instead -- does the newly produced inclusive sample describe the data, and does the
answer differ between preEE and postEE.

Three things kept honest here, each of which has bitten this analysis before:

* **Thresholds are computed once, from truth-matched signal, and reused for both
  eras.** Signal MC is 100% preEE, so there is no per-era alternative; re-deriving a
  quantile per era would silently make the two rows different cuts.

* **Nothing is ever renormalized to the data.** W+jets carries sigma*L, so its height
  is a prediction, not a template; scaling it to the data sideband would convert the
  result being measured into a flat line at 1. Both figures are at the absolute
  normalization, and "shape" is read off the per-bin ratio panel.

* **Yields are quoted on the sideband only, for both sides.** OS data is SR-blinded
  (verified: zero rows within +-40 MeV in both eras) while W+jets is not, so any
  ratio taken over the full range compares a hole against a peak. The figure is the
  other way round -- W+jets drawn full-range with the blind window greyed -- because
  hiding the MC Ds peak would hide the thing worth looking at.

* **N_eff, not row count, decides whether a number is readable.** The inclusive
  sample carries |w| ~ 9.5-10, an order of magnitude above the jet-binned files, so
  its effective statistics run out much earlier under a tight cut. Cells below
  MIN_NEFF are printed but marked: a weighted histogram can put a convincing-looking
  bump where three events sit.

This is NOT the preEE-vs-postEE era control test that CLAUDE.md lists as blocking
the "the new model fixed the mass tilt" claim. The eras differ in pileup and
detector conditions on their own, so a per-era difference in data/MC closure cannot
be attributed to the model's background-is-77%-postEE / signal-is-100%-preEE
asymmetry. That control still has to compare score distributions between eras.
"""

import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D = 1.8696
EFFS = [1.0, 0.90, 0.70, 0.50]
MIN_NEFF = 30
RATIO_YMAX = 3.0            # see the ratio panel below for why this is not autoscaled

# Shape probe: the low/high double ratio is independent of normalization, so it
# survives any overall mis-normalization. These are the same two windows used
# throughout the 0811/0814 studies, neither touching the blind region.
#
# There is deliberately NO W+jets scale factor here. A provisional SF = 20.79/21.03
# was carried on 2026-08-25 and removed the same day: it was never traced to a source
# and it predates the cross-section update (incl 63990 / 0J 51846 / 1J 8860 /
# 2J 3284 pb), so it was correcting a normalization that no longer exists. W+jets
# carries sigma*L and nothing else -- every absolute yield below is the raw
# prediction.
LOW_MAX, HIGH_MIN = 1.830, 2.008

# Which W+jets description to compare against. The two are ALTERNATIVES for the same
# physical process, never components -- the jet-binned set and the inclusive sample are
# separate MiniAOD productions with their own sumGenWeights_ and their own cross
# sections (incl 63990 = 0J 51846 + 1J 8860 + 2J 3284 pb, exactly). Summing them would
# double-count. Producing both figures is how the two descriptions get compared.
SOURCES = {
    "incl": (["WtoLNu2Jets_incl_{era}"], "W+jets incl", "incl"),
    "bins": (["WtoLNu2Jets_0J_{era}", "WtoLNu2Jets_1J_{era}", "WtoLNu2Jets_2J_{era}"],
             "W+jets 0J+1J+2J", "bins"),
}

ERA_INFO = {
    "2022 (preEE)":    ("Data_2022",   "2022",   "7.99 fb$^{-1}$"),
    "2022EE (postEE)": ("Data_2022EE", "2022EE", "26.68 fb$^{-1}$"),
}

# Lumi_reweight is deliberately absent from BR_SIG: the signal files carry no such
# branch (they need no luminosity weight), and asking for it there is a hard error.
BR = ["score_dsmu", "Ds_mass", "phi_mass", "Lumi_reweight", "event_id",
      "num_os_mu_tight", "num_ss_mu_tight", common.B_VETO]
BR_SIG = [b for b in BR if b != "Lumi_reweight"] + ["Ds_match", "mu_match"]

argv = sys.argv[1:]
SOURCE = "incl"
# OS 还是 SS。SS 两处结构性差异: (1) SS data 没有 blind(SS 结构上排除信号), 所以谱图上
# 它的 SR 不是空的; (2) SS data 有一半(偶数 event_id)参与了 ssdata 训练, 只能用奇数半,
# 否则自证自洁 —— W+jets SS 从未参与训练, 不切 parity。
SIGN = "OS"
# phi 窗半宽(GeV)。None = 用 tuplizer 自带的宽窗 M(KK) < Mphi+0.05。
# 卡紧的理由不是审美: 宽窗里混着非共振 K+K-pi+ 和 K+K*0bar, 它们在质量谱上不成峰而是
# 宽分布, 把连续本底的形状搞复杂; 实测卡到 +-10 MeV 之后 D+- 峰的 sideband 抽取才稳定
# (两组质量窗口给 0.644/0.652, 宽窗时是 0.868/1.160, 差 34%)。
PHI_HW = None
for a in list(argv):
    if a in ("OS", "SS"):
        SIGN = a
        argv.remove(a)
    elif a.startswith("phi="):
        PHI_HW = float(a.split("=")[1])
        argv.remove(a)
for a in list(argv):
    if a in SOURCES:
        SOURCE = a
        argv.remove(a)
WJ_NAMES, WJ_LABEL, WJ_TAG = SOURCES[SOURCE]

SPEC = argv[0] if argv else "ssdata"
parts = SPEC.split(":")
TOOL = parts.pop(0) if parts and parts[0] in common.TOOLS else "DNN"
SEED = int(parts.pop()[4:]) if parts and parts[-1].startswith("seed") else common.SEED_DEFAULT
if len(parts) != 1:
    sys.exit(f"cannot parse model spec '{SPEC}'")
BKG = parts[0]
common.resolve(BKG)
PRED = f"{common.out_dir_for(BKG, tool=TOOL, seed=SEED)}/Prediction"


def load(names, *, sig=False, sideband=False, is_data=False):
    """One sample name or a list of them; a list is concatenated.

    Concatenating is correct for the jet bins because each file already carries its
    own per-event Lumi_reweight (its own sigma and its own sumGenWeights_ divisor), so
    the three are already on a common absolute scale before they are stacked.
    """
    if isinstance(names, str):
        names = [names]
    parts_ = []
    for name in names:
        path = f"{PRED}/{name}.root"
        if not os.path.exists(path):
            sys.exit(f"missing {path} -- run predict.py for this sample first")
        parts_.append(uproot.open(path)["Events"].arrays(BR_SIG if sig else BR,
                                                         library="np"))
    a = parts_[0] if len(parts_) == 1 else {
        k: np.concatenate([p[k] for p in parts_]) for k in parts_[0]}
    if SIGN == "SS" and not sig:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
    if PHI_HW is not None:
        m &= np.abs(a["phi_mass"] - 1.019460) < PHI_HW
    if SIGN == "SS" and is_data:
        m &= (a["event_id"] % 2) == 1          # 训练用了偶数半
    if sig:
        # Truth-matched only: an unmatched signal row is combinatorial background
        # wearing the signal file's name, and would drag the quantiles low.
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if sideband:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    return {k: v[m] for k, v in a.items()}


print(f"Model: {TOOL}/{BKG} seed{SEED}   source: {SOURCE} ({WJ_LABEL})   sign: {SIGN}")

sig_parts = [load(n, sig=True) for n in common.SIGNALS]
sig_score = np.concatenate([p["score_dsmu"] for p in sig_parts])
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig_score, 1.0 - e))) for e in EFFS}
print(f"Signal (truth-matched, both channels): {len(sig_score):,} rows")
for e in EFFS:
    if e < 1.0:
        print(f"  eps_sig={e:.0%}  ->  score > {THR[e]:.6f}")

PHITAG = "" if PHI_HW is None else f"_phi{PHI_HW*1000:.0f}"
SFX = "" if SIGN == "OS" else "_SS"
ERAS = {era: (dname + SFX, [n.format(era=etag) + SFX for n in WJ_NAMES], lumi)
        for era, (dname, etag, lumi) in ERA_INFO.items()}

D = {}
for era, (dname, wnames, _) in ERAS.items():
    # sideband 定义对 OS/SS 一致, 只为让两者可比; SS 并不因 blind 而需要它。
    D[era] = dict(data=load(dname, sideband=True, is_data=True),
                  wj=load(wnames, sideband=True),
                  wj_full=load(wnames),
                  data_full=load(dname, is_data=True))
    print(f"  {era:<17} data {len(D[era]['data']['Ds_mass']):>9,} | "
          f"{WJ_LABEL} {len(D[era]['wj_full']['Ds_mass']):>9,} rows "
          f"({', '.join(wnames)})")


def cell(era, eff):
    """Sideband yields for one (era, working point): data, W+jets, and the tilt."""
    thr = THR[eff]
    d = D[era]["data"]
    w = D[era]["wj"]
    md, mw = d["Ds_mass"][d["score_dsmu"] > thr], w["Ds_mass"][w["score_dsmu"] > thr]
    ww = w["Lumi_reweight"][w["score_dsmu"] > thr]

    n_d = len(md)
    n_w, n_w2 = ww.sum(), (ww ** 2).sum()
    # Effective entries of a weighted sample: (sum w)^2 / sum w^2.  With +-|w| this
    # collapses fast, which is exactly what has to be visible before reading a shape.
    n_eff = n_w ** 2 / n_w2 if n_w2 > 0 else 0.0

    def band(mass, weight):
        lo = weight[mass < LOW_MAX].sum() if weight is not None else (mass < LOW_MAX).sum()
        hi = weight[mass > HIGH_MIN].sum() if weight is not None else (mass > HIGH_MIN).sum()
        return lo, hi

    dl, dh = band(md, None)
    wl, wh = band(mw, ww)
    # Double ratio (data_low/MC_low) / (data_high/MC_high): normalization cancels
    # twice over, so the luminosity weight cannot move it.
    if min(dl, dh, wl, wh) > 0:
        dr = (dl / wl) / (dh / wh)
        # Data Poisson only; the MC term needs sum w^2 per window, kept below.
        wl2 = (ww[mw < LOW_MAX] ** 2).sum()
        wh2 = (ww[mw > HIGH_MIN] ** 2).sum()
        rel = np.sqrt(1 / dl + 1 / dh + wl2 / wl ** 2 + wh2 / wh ** 2)
        dr_err = dr * rel
    else:
        dr = dr_err = float("nan")

    ratio = n_d / n_w if n_w > 0 else float("nan")
    r_err = ratio * np.sqrt(1 / max(n_d, 1) + n_w2 / max(n_w, 1e-9) ** 2)
    return dict(n_d=n_d, n_w=n_w, raw_w=len(mw), n_eff=n_eff,
                ratio=ratio, r_err=r_err, dr=dr, dr_err=dr_err)


print("\n" + "=" * 108)
print("OS sideband only (|M - M_Ds| >= 40 MeV on BOTH sides -- OS data is blinded there)")
print("=" * 108)
hdr = (f"{'era':<17} {'eps_sig':>8} {'N(data)':>10} {'N(W+jets)':>12} {'raw':>8} "
       f"{'N_eff':>8} {'data/MC':>16} {'low/high 双比值':>20}")
for era in ERAS:
    print("-" * 108)
    print(hdr if era == list(ERAS)[0] else "")
    for eff in EFFS:
        c = cell(era, eff)
        flag = " <- N_eff 太低" if c["n_eff"] < MIN_NEFF else ""
        lbl = "no cut" if eff >= 1.0 else f"{eff:.0%}"
        print(f"{era:<17} {lbl:>8} {c['n_d']:>10,} {c['n_w']:>12,.0f} {c['raw_w']:>8,} "
              f"{c['n_eff']:>8,.0f} {c['ratio']:>8.3f} +- {c['r_err']:.3f} "
              f"{c['dr']:>11.3f} +- {c['dr_err']:.3f}{flag}")

# ---------- figure ----------
# Absolute yields, no shape normalization: the W+jets weight makes its height a real
# prediction and the data points are real counts, so the vertical gap between them is
# the result rather than a plotting choice.  W+jets is drawn over the full range with
# the blind window greyed; OS data is simply empty there.
# 10 MeV (2026-08-25, 从 20 改过来). blind 窗仍落在 bin 边界 (0.04/0.01 = 4).
# 理由是实测: 20 MeV 把 1.87 处的 D+- 峰和 MC 的 Ds 峰形状都抹进连续本底看不见,
# 10 MeV 下两者都露出来。代价是 data 连续本底每 bin 误差从 ~7% 到 ~10% (eps_sig=90%),
# 换得两个此前不可见的结构 —— 划算。5 MeV 不值: 误差到 14% 而分辨率不再改善。
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.01)
mcent = 0.5 * (mbins[1:] + mbins[:-1])

fig, axes = plt.subplots(len(ERAS), len(EFFS),
                         figsize=(8.5 * len(EFFS), 7.5 * len(ERAS)), squeeze=False)
for row, (era, (_, _, lumi)) in enumerate(ERAS.items()):
    for col, eff in enumerate(EFFS):
        ax = axes[row][col]
        thr = THR[eff]
        w = D[era]["wj_full"]
        kw = w["score_dsmu"] > thr
        hw, _ = np.histogram(w["Ds_mass"][kw], bins=mbins,
                             weights=w["Lumi_reweight"][kw])
        d = D[era]["data"] if SIGN == "OS" else D[era]["data_full"]
        kd = d["score_dsmu"] > thr
        hd, _ = np.histogram(d["Ds_mass"][kd], bins=mbins)

        hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"{WJ_LABEL} {SIGN} MC  N={hw.sum():,.0f}")
        hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd),
                     color="black", marker="s", markersize=4.5, elinewidth=1.4,
                     label=f"{SIGN} data{' (SR blind)' if SIGN=='OS' else ' (odd parity)'}  N={hd.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel("Events / 10 MeV")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{era} -- {lbl}" + ("" if PHI_HW is None else rf"   $|M_{{KK}}-M_\phi|<${PHI_HW*1000:.0f} MeV"), fontsize=14)
        ax.legend(fontsize=11, loc="upper right")
        if col == 0:
            hep.cms.label(ax=ax, rlabel=f"{lumi} (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path("era", f"compare_{WJ_TAG}_{SIGN}{PHITAG}_by_era_mass.png",
                      key=BKG, tool=TOOL, seed=SEED)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"\nSaved {out}")

# ---------- per-bin data/MC, at the absolute normalization ----------
# Deliberately NOT renormalized to the data sideband.  W+jets carries sigma*L, so
# its height is a prediction rather than a template: rescaling it to the data would
# throw away the very quantity being tested and turn a 30% over-prediction into a
# flat line at 1.  So this panel answers "shape" and "how many" together -- a
# horizontal band away from 1 is a normalization offset, a slope across mass is the
# tilt, and the two are readable off the same picture.
fig, axes = plt.subplots(len(ERAS), len(EFFS),
                         figsize=(8.5 * len(EFFS), 6.0 * len(ERAS)), squeeze=False)
for row, era in enumerate(ERAS):
    for col, eff in enumerate(EFFS):
        ax = axes[row][col]
        thr = THR[eff]
        w, d = D[era]["wj"], D[era]["data"]
        kw, kd = w["score_dsmu"] > thr, d["score_dsmu"] > thr
        ww = w["Lumi_reweight"][kw]
        hw, _ = np.histogram(w["Ds_mass"][kw], bins=mbins, weights=ww)
        hw2, _ = np.histogram(w["Ds_mass"][kw], bins=mbins, weights=ww ** 2)
        hd, _ = np.histogram(d["Ds_mass"][kd], bins=mbins)
        # Only bins where both sides are populated and the MC bin is positive:
        # a negative weighted bin (which the preEE sample does produce at tight
        # cuts) has no meaningful ratio at all.
        ok = (hd > 0) & (hw > 0)
        r = np.where(ok, hd / np.where(hw > 0, hw, 1), np.nan)
        rerr = np.where(ok, r * np.sqrt(1 / np.maximum(hd, 1) + hw2 / np.maximum(hw, 1e-9) ** 2), np.nan)
        # Fixed y-range. Autoscaling is actively misleading here: one outlying bin
        # plus the metre-long error bars of a nearly-exhausted MC stretch the axis
        # over tens of units, and every point that matters -- the ones between 0.4
        # and 1.5 -- collapses onto the reference line. Points outside the window are
        # drawn as triangles at the edge so nothing is silently dropped.
        inside = ok & (r <= RATIO_YMAX)
        ax.errorbar(mcent[inside], r[inside], yerr=rerr[inside], fmt="s",
                    color="black", markersize=4.5, elinewidth=1.4)
        above = ok & (r > RATIO_YMAX)
        if above.any():
            ax.plot(mcent[above], np.full(above.sum(), RATIO_YMAX * 0.97), "^",
                    color="black", markersize=7)
        ax.set_ylim(0, RATIO_YMAX)
        ax.axhline(1.0, color="firebrick", lw=1.4)
        sb = np.abs(mcent - M_DS) >= SR_HALFWIDTH
        n_bad = int((~ok & sb).sum())
        notes = []
        if n_bad:
            notes.append(f"{n_bad} bin 无法作比（MC<=0 或 data=0）")
        if above.any():
            notes.append(f"{int(above.sum())} bin 超出画面（三角）")
        if notes:
            ax.text(0.02, 0.94, "\n".join(notes), transform=ax.transAxes,
                    fontsize=10, color="firebrick", va="top")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel(f"{SIGN} data / {WJ_LABEL}")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{era} -- {lbl}" + ("" if PHI_HW is None else rf"   $|M_{{KK}}-M_\phi|<${PHI_HW*1000:.0f} MeV"), fontsize=14)
plt.tight_layout()
out = common.fig_path("era", f"compare_{WJ_TAG}_{SIGN}{PHITAG}_by_era_ratio.png",
                      key=BKG, tool=TOOL, seed=SEED)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

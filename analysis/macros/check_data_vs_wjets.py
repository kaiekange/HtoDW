"""Does the W+jets background model get worse as the DNN cut tightens?

    python3 check_data_vs_wjets.py wjets ssdata

W+jets MC is an incomplete background -- tt/DY/QCD are simply absent from it -- so it
undershoots data even before any MVA cut. The question that matters is whether the
shortfall is FLAT in the score: a constant deficit can be absorbed into a
normalization, but a deficit that grows with the cut means the classifier is
preferentially rejecting the part of the background the MC does model, leaving the
unmodelled part behind, and the extrapolation into the signal region gets worse
exactly where the analysis lives.

Reports data/MC in the control region (|M-M_Ds| >= 40 MeV) at matched signal
efficiency working points, for both charge configurations and both trained models.

Details that matter:

* Working points are per-model signal-efficiency quantiles, never raw score -- the
  two models' score scales are unrelated.
* SS data is restricted to the odd event_id half, because the ssdata model trained on
  the even half. The W+jets side is scaled by the measured odd fraction so the ratio
  stays a like-for-like comparison rather than picking up a factor 2.
* OS data is blinded in |M-M_Ds| < 40 MeV, so only the control region is quoted for
  it. That is the region the fit extrapolates from anyway.
* The analysis selection (1 tight muon + tight b-veto) is applied throughout; the
  veto alone already removes the largest missing process (tt).

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_data_vs_wjets.py wjets ssdata
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

OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
os.makedirs(OUTDIR, exist_ok=True)
M_DS, SR_HALFWIDTH = 1.96835, 0.04
EFFS = [0.90, 0.70, 0.50]
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      common.B_VETO, "event_id"]

if len(sys.argv) < 2:
    sys.exit(f"usage: python3 {os.path.basename(__file__)} <model> [<model> ...]")
MODELS = sys.argv[1:]
for k in MODELS:
    common.resolve(k)


def load(model_key, name, *, ss=False, sig=False, parity=None):
    br = BR + (["Ds_match", "mu_match"] if sig else [])
    pred = f"{common.resolve(model_key)['out_dir']}/Prediction/{name}.root"
    a = uproot.open(pred)["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if parity is not None:
        m &= (a["event_id"] % 2) == parity
    return {k: v[m] for k, v in a.items()}


def wjets(model_key, ss):
    v, w = [], []
    for n, wgt in common.WJETS_WEIGHTS.items():
        d = load(model_key, f"{n}_SS" if ss else n, ss=ss)
        v.append(d)
        w.append(np.full(len(d["Ds_mass"]), wgt))
    return ({k: np.concatenate([p[k] for p in v]) for k in BR}, np.concatenate(w))


def cr_yield(d, thr, weight=None):
    """Weighted yield in the control region above a score threshold."""
    m = (d["score_dsmu"] > thr) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    return float(m.sum()) if weight is None else float(weight[m].sum())


print("Loading scored samples...")
curves = {}
keep = {}
THR = {}
for k in MODELS:
    sig = {b: np.concatenate([load(k, n, sig=True)[b] for n in common.SIGNALS])
           for b in ["score_dsmu"]}
    os_d = load(k, "Data_2022")
    ss_d_all = load(k, "Data_2022_SS", ss=True)
    ss_d = load(k, "Data_2022_SS", ss=True, parity=1)          # held-out half
    f_odd = len(ss_d["Ds_mass"]) / len(ss_d_all["Ds_mass"])    # ~0.5, measured
    os_m, os_mw = wjets(k, ss=False)
    ss_m, ss_mw = wjets(k, ss=True)
    keep[k] = dict(ss_d=ss_d, ss_m=ss_m, ss_mw=ss_mw * f_odd,
                   os_d=os_d, os_m=os_m, os_mw=os_mw)
    THR[k] = {e: np.quantile(sig["score_dsmu"], 1 - e) for e in EFFS}

    print(f"\n=== {common.resolve(k)['label']} model "
          f"(SS data uses odd event_id only, {100*f_odd:.1f}% of SS; "
          f"W+jets scaled by that fraction) ===")
    print(f"{'eff_sig':>8s} {'cut':>9s} {'OS data':>10s} {'OS W+jets':>11s} {'OS d/MC':>9s} "
          f"{'SS data':>10s} {'SS W+jets':>11s} {'SS d/MC':>9s}")
    rows = []
    for eff in [1.0] + EFFS:
        thr = -1.0 if eff >= 1.0 else np.quantile(sig["score_dsmu"], 1 - eff)
        od, om = cr_yield(os_d, thr), cr_yield(os_m, thr, os_mw)
        sd, sm = cr_yield(ss_d, thr), cr_yield(ss_m, thr, ss_mw) * f_odd
        r_os = od / om if om > 0 else np.nan
        r_ss = sd / sm if sm > 0 else np.nan
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        print(f"{lbl:>8s} {thr:9.4f} {od:10,.0f} {om:11,.0f} {r_os:9.3f} "
              f"{sd:10,.0f} {sm:11,.0f} {r_ss:9.3f}")
        rows.append((eff, r_os, r_ss))
    curves[k] = rows

# ---------- plot ----------
fig, ax = plt.subplots(figsize=(10, 8))
STYLE = {"wjets": dict(color="seagreen", linestyle="-"),
         "ssdata": dict(color="firebrick", linestyle="--")}
for k in MODELS:
    st = STYLE.get(k, dict(color="grey", linestyle="-"))
    xs = [100 * r[0] for r in curves[k]]
    ax.plot(xs, [r[1] for r in curves[k]], marker="o", markersize=9, lw=2, **st,
            label=f"OS -- {common.resolve(k)['label']} model")
    ax.plot(xs, [r[2] for r in curves[k]], marker="^", markersize=9, lw=2,
            alpha=0.55, **st, label=f"SS -- {common.resolve(k)['label']} model")
ax.axhline(1.0, color="black", lw=1, ls=":")
ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
ax.set_ylabel("data / W+jets MC in control region")
ax.invert_xaxis()
ax.legend(fontsize=12)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = f"{OUTDIR}/data_vs_wjets_vs_score.png"
plt.savefig(out, dpi=130)
plt.close(fig)
print(f"\nSaved {out}")


# ---------- Can W+jets MC represent the data after the MVA cut? ----------
# The lumi-weighted W+jets is an ABSOLUTE prediction: cross-section x luminosity, no
# free normalization. So the first requirement is simply that its event count equals
# the data's -- in the OS control region, the SS control region, and the SS signal
# region alike. Renormalizing the MC to the data in the CR would hide exactly the
# failure being looked for, by absorbing the missing events into a fitted factor.
# Shape agreement is only meaningful once the yields already match.
#
# Ratio = 1 in all three regions -> W+jets represents the data.
# Ratio != 1, or ratio drifting between regions or with the cut -> it does not.
#
# SS is the only place the SR can be checked at all: OS data is blinded there.
COLS = [1.0] + EFFS
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SB_BIN = np.abs(mcent - M_DS) >= SR_HALFWIDTH

def yield_ratio(nd, nm, nm2):
    """data/MC with data Poisson + MC weighted error. nm2 = sum of w^2."""
    if nm <= 0 or nd <= 0:
        return np.nan, np.nan
    r = nd / nm
    return r, r * np.sqrt(1 / nd + nm2 / nm ** 2)


fig, axes = plt.subplots(len(MODELS), len(COLS),
                         figsize=(8.5 * len(COLS), 7.5 * len(MODELS)), squeeze=False)
print("\n=== Does the lumi-weighted W+jets reproduce the observed yield? ===")
print("   Absolute prediction, no free normalization. 1.00 = W+jets represents the data.")
print("   OS SR is blinded and cannot be checked at all.\n")
print(f"{'model':8s} {'eff_sig':>8s} {'OS CR':>16s} {'SS CR':>16s} {'SS SR':>16s}")
ratios = {k: [] for k in MODELS}
for row, k in enumerate(MODELS):
    d, m, mw = keep[k]["ss_d"], keep[k]["ss_m"], keep[k]["ss_mw"]
    od, om, omw = keep[k]["os_d"], keep[k]["os_m"], keep[k]["os_mw"]
    for col, eff in enumerate(COLS):
        ax = axes[row][col]
        thr = -1.0 if eff >= 1.0 else THR[k][eff]
        md, mm = d["score_dsmu"] > thr, m["score_dsmu"] > thr
        hd, _ = np.histogram(d["Ds_mass"][md], bins=mbins)
        hm, _ = np.histogram(m["Ds_mass"][mm], bins=mbins, weights=mw[mm])
        hm2, _ = np.histogram(m["Ds_mass"][mm], bins=mbins, weights=mw[mm] ** 2)

        # OS control region (OS SR is blind, so only the CR exists there)
        ocr = (od["score_dsmu"] > thr) & (np.abs(od["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
        mcr = (om["score_dsmu"] > thr) & (np.abs(om["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
        r_os = yield_ratio(float(ocr.sum()), float(omw[mcr].sum()), float((omw[mcr] ** 2).sum()))
        r_scr = yield_ratio(hd[SB_BIN].sum(), hm[SB_BIN].sum(), hm2[SB_BIN].sum())
        r_ssr = yield_ratio(hd[~SB_BIN].sum(), hm[~SB_BIN].sum(), hm2[~SB_BIN].sum())
        ratios[k].append((eff, r_os, r_scr, r_ssr))

        # absolute yields -- the MC is NOT rescaled to the data
        hep.histplot(hm, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"SS W+jets MC (lumi-weighted)  N={hm.sum():,.0f}")
        hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd),
                     color="black", marker="o", markersize=5, elinewidth=1.4,
                     label=f"SS data  N={hd.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel("Events / 20 MeV")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{common.resolve(k)['label']} model -- {lbl}\n"
                     f"data/MC: CR {r_scr[0]:.2f}, SR {r_ssr[0]:.2f}", fontsize=14)
        ax.legend(fontsize=10)
        if row == 0 and col == 0:
            hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=11)

        f2 = lambda r: "n/a" if np.isnan(r[0]) else f"{r[0]:.2f} +- {r[1]:.2f}"
        lbl_t = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        print(f"{k:8s} {lbl_t:>8s} {f2(r_os):>16s} {f2(r_scr):>16s} {f2(r_ssr):>16s}")

plt.tight_layout()
out = f"{OUTDIR}/ss_data_vs_wjets_mass.png"
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"\nSaved {out}")


def bin_ratio(hd, hm, hm2):
    """Per-bin data/MC with error; nan (not 0) where either side is empty -- an
    empty bin is missing information, not a measured ratio of zero, and plotting it
    as 0 would misread as 'MC way overshoots here' instead of 'no data fell in this
    bin' (which is exactly what happens in every OS bin inside the blind window)."""
    r = np.full(len(hd), np.nan)
    e = np.full(len(hd), np.nan)
    ok = (hd > 0) & (hm > 0)
    r[ok] = hd[ok] / hm[ok]
    e[ok] = r[ok] * np.sqrt(1 / hd[ok] + hm2[ok] / hm[ok] ** 2)
    return r, e


# ---------- per-bin data/MC ratio across M(KKpi), OS and SS, at each cut ----------
# The integrated CR/SR numbers above answer "is the total off" -- this answers "is it
# off in a mass-dependent way", which is what would actually distort a fit. OS bins
# inside the blind window are structurally empty and are left as gaps, not zeros.
for k in MODELS:
    d, m, mw = keep[k]["ss_d"], keep[k]["ss_m"], keep[k]["ss_mw"]
    od, om, omw = keep[k]["os_d"], keep[k]["os_m"], keep[k]["os_mw"]
    fig, axes = plt.subplots(2, len(COLS), figsize=(8 * len(COLS), 6 * 2),
                             sharey="row", squeeze=False)
    for col, eff in enumerate(COLS):
        thr = -1.0 if eff >= 1.0 else THR[k][eff]
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"

        for row, (name, dd, mm_, ww) in enumerate((("OS", od, om, omw), ("SS", d, m, mw))):
            ax = axes[row][col]
            sd, sm = dd["score_dsmu"] > thr, mm_["score_dsmu"] > thr
            hd, _ = np.histogram(dd["Ds_mass"][sd], bins=mbins)
            hm, _ = np.histogram(mm_["Ds_mass"][sm], bins=mbins, weights=ww[sm])
            hm2, _ = np.histogram(mm_["Ds_mass"][sm], bins=mbins, weights=ww[sm] ** 2)
            r, e = bin_ratio(hd, hm, hm2)

            ax.errorbar(mcent, r, yerr=e, fmt="o", markersize=5, color="black",
                        elinewidth=1.3, capsize=2)
            ax.axhline(1.0, color="firebrick", lw=1.3, ls=":")
            ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
            # A handful of low-count, large-error bins (esp. at tight cuts) would
            # otherwise stretch the axis so far that every informative point flattens
            # to the bottom. Cap the view; the excluded points are noise, not signal --
            # their error bars already said so.
            ax.set_ylim(0, 6)
            ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
            if col == 0:
                ax.set_ylabel(f"{name} data / W+jets MC")
            ax.set_title(f"{name} -- {lbl}", fontsize=14)
            if row == 0 and col == 0:
                hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True,
                              text="Preliminary", fontsize=10)
    plt.tight_layout()
    out = f"{OUTDIR}/data_wjets_binratio_{k}.png"
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"Saved {out}")

# ---------- fine scan in signal efficiency, one figure per model ----------
# The three coarse working points only show that the ratios move; a dense scan shows
# WHERE they start to break and how fast, which is what picking a working point needs.
# Event counts are printed alongside because the tight end runs out of statistics long
# before the trend does -- past the point where SS SR holds only a handful of events
# the curve is noise, and the error bars are the only honest guide.
SCAN = [1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70,
        0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20]

for k in MODELS:
    sig_sc = np.concatenate([load(k, n, sig=True)["score_dsmu"] for n in common.SIGNALS])
    d, m, mw = keep[k]["ss_d"], keep[k]["ss_m"], keep[k]["ss_mw"]
    od, om, omw = keep[k]["os_d"], keep[k]["os_m"], keep[k]["os_mw"]
    d_sb = np.abs(d["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    m_sb = np.abs(m["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    od_sb = np.abs(od["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    om_sb = np.abs(om["Ds_mass"] - M_DS) >= SR_HALFWIDTH

    print(f"\n=== fine scan, {common.resolve(k)['label']} model ===")
    print(f"{'eff_sig':>8s} {'cut':>9s} {'OS CR':>15s} {'SS CR':>15s} {'SS SR':>15s} "
          f"{'N(SS SR data)':>14s}")
    scan = []
    for eff in SCAN:
        thr = -1.0 if eff >= 1.0 else np.quantile(sig_sc, 1 - eff)
        sd, sm = d["score_dsmu"] > thr, m["score_dsmu"] > thr
        sod, som = od["score_dsmu"] > thr, om["score_dsmu"] > thr
        r_os = yield_ratio(float((sod & od_sb).sum()),
                           float(omw[som & om_sb].sum()),
                           float((omw[som & om_sb] ** 2).sum()))
        r_cr = yield_ratio(float((sd & d_sb).sum()),
                           float(mw[sm & m_sb].sum()),
                           float((mw[sm & m_sb] ** 2).sum()))
        r_sr = yield_ratio(float((sd & ~d_sb).sum()),
                           float(mw[sm & ~m_sb].sum()),
                           float((mw[sm & ~m_sb] ** 2).sum()))
        n_ssr = int((sd & ~d_sb).sum())
        scan.append((eff, r_os, r_cr, r_sr, n_ssr))
        f3 = lambda r: "n/a" if np.isnan(r[0]) else f"{r[0]:.2f}+-{r[1]:.2f}"
        print(f"{100*eff:7.0f}% {thr:9.4f} {f3(r_os):>15s} {f3(r_cr):>15s} "
              f"{f3(r_sr):>15s} {n_ssr:14,}")

    fig, ax = plt.subplots(figsize=(11, 8))
    for name, idx, c, mk in (("OS CR", 1, "steelblue", "o"),
                             ("SS CR", 2, "seagreen", "^"),
                             ("SS SR", 3, "firebrick", "s")):
        xs = [100 * s[0] for s in scan]
        ys = [s[idx][0] for s in scan]
        es = [s[idx][1] for s in scan]
        ax.errorbar(xs, ys, yerr=es, marker=mk, markersize=7, lw=2, capsize=3,
                    color=c, label=name)
    ax.axhline(1.0, color="black", lw=1.5, ls=":")
    ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
    ax.set_ylabel("data / lumi-weighted W+jets")
    ax.set_title(f"{common.resolve(k)['label']} model", fontsize=16)
    ax.invert_xaxis()
    ax.legend(fontsize=13)
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                  fontsize=11)
    plt.tight_layout()
    out = f"{OUTDIR}/wjets_yield_scan_{k}.png"
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out}")

fig, ax = plt.subplots(figsize=(11, 8))
MK = [("OS CR", 1, "o", "-"), ("SS CR", 2, "^", "--"), ("SS SR", 3, "s", ":")]
for k in MODELS:
    base = STYLE.get(k, dict(color="grey"))["color"]
    for name, idx, mk, ls in MK:
        xs = [100 * r[0] for r in ratios[k]]
        ys = [r[idx][0] for r in ratios[k]]
        es = [r[idx][1] for r in ratios[k]]
        ax.errorbar(xs, ys, yerr=es, marker=mk, linestyle=ls, markersize=8, lw=2,
                    capsize=4, color=base, alpha=1.0 if idx == 1 else 0.6,
                    label=f"{name} -- {common.resolve(k)['label']} model")
ax.axhline(1.0, color="black", lw=1.5, ls=":")
ax.text(99, 1.03, "W+jets reproduces data", fontsize=11, color="black")
ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
ax.set_ylabel("data / lumi-weighted W+jets")
ax.set_yscale("log")
ax.invert_xaxis()
ax.legend(fontsize=10, ncol=2)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = f"{OUTDIR}/wjets_yield_closure.png"
plt.savefig(out, dpi=130)
print(f"Saved {out}")


# ---------- OS control region at the tightest scanned working point ----------
# eff_sig=20% is the far end of SCAN above -- included here because it is the tight
# extreme where the closure plots showed the largest deviation and where OS (the
# region that actually matters -- SS is a proxy) has never been shown directly.
EFF_TIGHT = 0.20
for k in MODELS:
    sig_sc = np.concatenate([load(k, n, sig=True)["score_dsmu"] for n in common.SIGNALS])
    od, om, omw = keep[k]["os_d"], keep[k]["os_m"], keep[k]["os_mw"]
    thr = np.quantile(sig_sc, 1 - EFF_TIGHT)
    sd = (od["score_dsmu"] > thr) & (np.abs(od["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    sm = (om["score_dsmu"] > thr) & (np.abs(om["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    hd, _ = np.histogram(od["Ds_mass"][sd], bins=mbins)
    hm, _ = np.histogram(om["Ds_mass"][sm], bins=mbins, weights=omw[sm])
    hm2, _ = np.histogram(om["Ds_mass"][sm], bins=mbins, weights=omw[sm] ** 2)
    r, e = yield_ratio(float(hd.sum()), float(hm.sum()), float(hm2.sum()))

    fig, ax = plt.subplots(figsize=(11, 8))
    hep.histplot(hm, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                 label=f"OS W+jets MC (lumi-weighted)  N={hm.sum():,.1f}")
    hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd), color="black",
                 marker="o", markersize=6, elinewidth=1.5,
                 label=f"OS data  N={hd.sum():,.0f}")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel("Events / 20 MeV")
    ax.set_title(f"{common.resolve(k)['label']} model -- "
                 rf"$\varepsilon_{{sig}}$={100*EFF_TIGHT:.0f}%  (OS CR, score>{thr:.4f})"
                 f"\ndata/MC in CR = {r:.2f} $\\pm$ {e:.2f}", fontsize=15)
    ax.legend(fontsize=12)
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
    plt.tight_layout()
    out = f"{OUTDIR}/os_cr_dist_eff20_{k}.png"
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\n{common.resolve(k)['label']} model, eff_sig=20%, thr={thr:.4f}: "
          f"OS CR data={hd.sum():,.0f}  W+jets={hm.sum():,.1f}  ratio={r:.2f}+-{e:.2f}")
    print(f"Saved {out}")

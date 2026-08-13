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
      common.B_VETO, "event_id", "Ds_pt"]

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


def double_ratio(a, b):
    """(data/MC in region A) / (data/MC in region B), for DISJOINT A and B.

    Why the double ratio is the number that matters for the SR: an overall
    normalization error -- a missing process, a mismodelled selection efficiency,
    anything that scales the whole spectrum -- cancels in it, leaving only the
    mass-dependent part. The absolute SS SR ratio moves partly because its own
    denominator (SS CR) moves, which is not an SR effect.

    A and B being disjoint mass windows of the same samples makes all four counts
    independent Poisson, so the relative errors add in quadrature with no
    correlation term.  a, b are (n_data, n_mc, sum_w2) triplets.
    """
    (nda, nma, n2a), (ndb, nmb, n2b) = a, b
    if min(nda, ndb) <= 0 or min(nma, nmb) <= 0:
        return np.nan, np.nan
    r = (nda / nma) / (ndb / nmb)
    rel2 = 1 / nda + n2a / nma ** 2 + 1 / ndb + n2b / nmb ** 2
    return r, r * np.sqrt(rel2)


def neff(nm, nm2):
    """Effective entry count of a weighted sum -- the honest N behind an MC yield.

    The 0J bin carries weight 2.33, so a weighted yield can be dominated by a few
    entries and its raw row count badly overstates how much it constrains. Every
    ratio quoted at a tight working point has to be read next to this."""
    return nm ** 2 / nm2 if nm2 > 0 else 0.0


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

        # Absolute yields -- the MC is NOT rescaled to the data.
        #
        # Displayed at FULL-SAMPLE equivalent. Internally the SS side runs at half
        # luminosity: the SS data is the odd event_id half (the ssdata model trained
        # on the even one) and the SS MC is scaled down by the same fraction to match,
        # which is correct for every ratio but makes the plotted yields half the
        # physical numbers. A legend reading "SS data N=1,170" next to "OS data
        # N=4,990" elsewhere then invites a comparison that is off by two. Both SS
        # curves are scaled by the SAME factor here, so every ratio, title and printed
        # table below is bit-identical -- this is presentation only.
        s = 1.0 / f_odd
        hep.histplot(hm * s, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"SS W+jets MC (lumi-weighted)  N={hm.sum()*s:,.0f}")
        hep.histplot(hd * s, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd) * s,
                     color="black", marker="o", markersize=5, elinewidth=1.4,
                     label=f"SS data  N={hd.sum()*s:,.0f}")
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

    # Low / high sideband split. The integrated CR ratio averages the two and can sit
    # at 1.0 while the two halves sit either side of it -- which is what a tilt in
    # M(KKpi) looks like, and a tilt is what actually distorts a sideband fit.
    # Low is where partially reconstructed real charm piles up (Ds->phi rho, phi l nu:
    # a neutral is lost so M(KKpi) lands below the true mass), so a low-side-only
    # excess is a specific, checkable physics statement rather than a generic slope.
    d_lo, d_hi = d["Ds_mass"] < M_DS - SR_HALFWIDTH, d["Ds_mass"] > M_DS + SR_HALFWIDTH
    m_lo, m_hi = m["Ds_mass"] < M_DS - SR_HALFWIDTH, m["Ds_mass"] > M_DS + SR_HALFWIDTH
    od_lo, od_hi = od["Ds_mass"] < M_DS - SR_HALFWIDTH, od["Ds_mass"] > M_DS + SR_HALFWIDTH
    om_lo, om_hi = om["Ds_mass"] < M_DS - SR_HALFWIDTH, om["Ds_mass"] > M_DS + SR_HALFWIDTH

    def counts(dmask, mmask, weights):
        """(n_data, n_mc, sum_w2) for one region."""
        return (float(dmask.sum()), float(weights[mmask].sum()),
                float((weights[mmask] ** 2).sum()))

    print(f"\n=== fine scan, {common.resolve(k)['label']} model ===")
    print(f"{'eff_sig':>8s} {'cut':>9s} {'OS CR':>15s} {'SS CR':>15s} {'SS SR':>15s} "
          f"{'SS SR/CR':>15s} {'N(SS SR)':>9s} {'Neff(MC)':>9s}")
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
        # SR relative to its own sideband: the overall normalization cancels.
        c_sr = counts(sd & ~d_sb, sm & ~m_sb, mw)
        c_cr = counts(sd & d_sb, sm & m_sb, mw)
        r_dbl = double_ratio(c_sr, c_cr)
        # Low vs high sideband, both charge configurations.
        lohi = {}
        for nm_, (dm, mm_, ww, lo, hi, mlo, mhi) in (
                ("SS", (sd, sm, mw, d_lo, d_hi, m_lo, m_hi)),
                ("OS", (sod, som, omw, od_lo, od_hi, om_lo, om_hi))):
            c_lo = counts(dm & lo, mm_ & mlo, ww)
            c_hi = counts(dm & hi, mm_ & mhi, ww)
            lohi[nm_] = (yield_ratio(*c_lo), yield_ratio(*c_hi), double_ratio(c_lo, c_hi))
        n_ssr = int((sd & ~d_sb).sum())
        ne_ssr = neff(c_sr[1], c_sr[2])
        scan.append((eff, r_os, r_cr, r_sr, n_ssr, r_dbl, lohi, ne_ssr))
        f3 = lambda r: "n/a" if np.isnan(r[0]) else f"{r[0]:.2f}+-{r[1]:.2f}"
        print(f"{100*eff:7.0f}% {thr:9.4f} {f3(r_os):>15s} {f3(r_cr):>15s} "
              f"{f3(r_sr):>15s} {f3(r_dbl):>15s} {n_ssr:9,} {ne_ssr:9,.0f}")

    # ---------- low vs high sideband, integrated (step 2) ----------
    # Printed separately because the per-bin scatter plot cannot be read by eye at
    # these statistics -- at eps_sig=50% the OS CR holds ~390 events over ~26 bins.
    # Two integrated numbers with errors can be read; 26 scattered points cannot.
    # ---------- where inside the low sideband does the excess sit? ----------
    # The low/high split says there is a tilt; it does not say what shape it has.
    # Two very different things live below the Ds window and they need separating:
    #   D+- core   -- D+- -> K+K-pi+ is a real, charge-correlated charm peak at 1.8696
    #   the rest   -- partially reconstructed Ds (phi rho+, phi l nu: a neutral is
    #                 lost so M(KKpi) lands low) makes a broad shoulder, not a peak
    # A peak-only excess and a broad-shoulder excess call for completely different
    # background models, so this is not cosmetic. Windows are chosen to tile the low
    # side without ever touching the OS blind region [1.929, 2.008].
    MWIN = [("low  <1.830", 1.66835, 1.830),
            ("D+- 1.83-1.91", 1.830, 1.910),
            ("1.91-1.929", 1.910, M_DS - SR_HALFWIDTH),
            ("high >2.008", M_DS + SR_HALFWIDTH, 2.26835)]

    # SS is printed alongside because it is the control that decides the
    # interpretation: SS contains the charge-symmetric background (combinatorics,
    # b->charm) but NOT the charge-correlated real charm (W+c, where the charm
    # hadron's charge is tied to the W's). A structure present in OS and absent in
    # SS is charge-correlated by construction; one present in both is not.
    win_scan = {}
    for cname, dd, mm_, ww in (("OS", od, om, omw), ("SS", d, m, mw)):
        print(f"\n--- {cname}: data/MC in mass windows (blind region never touched), "
              f"{common.resolve(k)['label']} model ---")
        print(f"{'eff_sig':>8s} " + " ".join(f"{n:>18s}" for n, _, _ in MWIN))
        rows = []
        for eff in SCAN:
            thr = -1.0 if eff >= 1.0 else np.quantile(sig_sc, 1 - eff)
            sdd, smm = dd["score_dsmu"] > thr, mm_["score_dsmu"] > thr
            cells, vals, nes = [], [], []
            for _, lo, hi in MWIN:
                wd = sdd & (dd["Ds_mass"] >= lo) & (dd["Ds_mass"] < hi)
                wm = smm & (mm_["Ds_mass"] >= lo) & (mm_["Ds_mass"] < hi)
                c = counts(wd, wm, ww)
                r, e = yield_ratio(*c)
                cells.append("n/a" if np.isnan(r) else f"{r:.3f}+-{e:.3f}")
                vals.append((r, e))
                nes.append(min(c[0], neff(c[1], c[2])))
            rows.append((eff, vals, min(nes)))
            lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
            print(f"{lbl:>8s} " + " ".join(f"{c:>18s}" for c in cells))
        win_scan[cname] = rows

    # ---------- the OS-vs-SS contrast, as a figure ----------
    # The central observation of this study is easier to see than to read off a table:
    # both charge configurations start together at 1.44 and then fan out in OPPOSITE
    # directions as the cut tightens. Same panel layout for both so the mirroring is
    # immediate; the D+- window is drawn too, since its staying flat is part of the
    # result rather than a null to be omitted.
    fig, axs = plt.subplots(1, 2, figsize=(20, 8), sharey=True)
    WCOL = ["steelblue", "seagreen", "darkorange", "firebrick"]
    WMK = ["o", "^", "v", "s"]
    for ax, cname in zip(axs, ("OS", "SS")):
        for j, (wname, _, _) in enumerate(MWIN):
            xs = [100 * e for e, _, _ in win_scan[cname]]
            ys = [v[j][0] for _, v, _ in win_scan[cname]]
            es = [v[j][1] for _, v, _ in win_scan[cname]]
            ax.errorbar(xs, ys, yerr=es, marker=WMK[j], markersize=7, lw=2, capsize=3,
                        color=WCOL[j], label=f"M {wname}")
        # Shade where the thinnest window drops below 30 effective entries. On a log
        # axis the huge error bars out there still LOOK like points, and the SS panel
        # runs out of statistics much earlier than the OS one -- without this the eye
        # reads a trend into the noise.
        lo_x = [100 * e for e, _, ne in win_scan[cname] if ne < 30]
        if lo_x:
            ax.axvspan(min(lo_x) - 2.5, max(lo_x) + 2.5, color="grey", alpha=0.18)
            ax.text(max(lo_x) + 2.0, 5.0, r"$N_{eff}<30$", fontsize=12,
                    color="dimgrey", ha="right", va="top")
        ax.axhline(1.0, color="black", lw=1.5, ls=":")
        ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
        ax.set_yscale("log")
        ax.set_ylim(0.2, 6)
        ax.invert_xaxis()
        ax.set_title(f"{cname}  ({common.resolve(k)['label']} model)", fontsize=16)
        ax.legend(fontsize=12)
    axs[0].set_ylabel("data / lumi-weighted W+jets")
    hep.cms.label(ax=axs[0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                  fontsize=11)
    plt.tight_layout()
    out = f"{OUTDIR}/mass_windows_os_vs_ss_{k}.png"
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\nSaved {out}")

    # ---------- is the mass tilt really about mass, or is it Ds_pt in disguise? ----------
    # Ds_pt is the highest-importance feature AND the worst-mismodelled single
    # variable measured (eps_data/eps_MC 2.52 at eps_MC=1e-2, rising to 4.57 at
    # 5e-4, check_feature_tail.py). M(KKpi) and Ds_pt are not independent for
    # combinatorial candidates, so a cut that leans on Ds_pt can project a pT
    # mismatch onto the mass spectrum and fake a mass-dependent effect.
    #
    # The OS-vs-SS contrast already argues against that -- a pT mismodelling is
    # charge-blind and would have to appear in SS too -- but arguing is not
    # measuring. Redo the mass windows inside narrow Ds_pt bands: if the OS low-mass
    # rise survives at fixed pT, it is a mass effect; if it collapses, it was pT all
    # along and the charge-correlated-charm reading has to be withdrawn.
    PT_BANDS = [(0.0, 12.0), (12.0, 18.0), (18.0, 1e9)]
    # Scanned all the way down, not stopped at 70%: the mass dependence keeps growing
    # past there, so cutting the scan short would hide the part that matters most.
    # N(OS low) is printed on every line so the point where it becomes noise is
    # visible rather than assumed.
    PT_EFFS = [1.0, 0.90, 0.70, 0.50, 0.30]

    print(f"\n--- mass tilt at FIXED Ds_pt (is the tilt really about mass?), "
          f"{common.resolve(k)['label']} model ---")
    print("    low = M<1.830, high = M>2.008; each side is an ABSOLUTE data/MC ratio")
    print("    (lumi-weighted W+jets is a sigma x L prediction -- there is no free")
    print("     normalization, so each window's ratio means something on its own)\n")
    print(f"{'eff_sig':>8s} {'Ds_pt band':>14s} {'OS low/high':>18s} {'SS low/high':>18s} "
          f"{'N(OS low)':>10s}")
    for eff in PT_EFFS:
        thr = -1.0 if eff >= 1.0 else np.quantile(sig_sc, 1 - eff)
        for plo, phi in PT_BANDS:
            cells = []
            n_os_lo = 0
            for cname, dd, mm_, ww in (("OS", od, om, omw), ("SS", d, m, mw)):
                pb_d = (dd["Ds_pt"] >= plo) & (dd["Ds_pt"] < phi)
                pb_m = (mm_["Ds_pt"] >= plo) & (mm_["Ds_pt"] < phi)
                sdd = (dd["score_dsmu"] > thr) & pb_d
                smm = (mm_["score_dsmu"] > thr) & pb_m
                c_lo = counts(sdd & (dd["Ds_mass"] < M_DS - SR_HALFWIDTH),
                              smm & (mm_["Ds_mass"] < M_DS - SR_HALFWIDTH), ww)
                c_hi = counts(sdd & (dd["Ds_mass"] > M_DS + SR_HALFWIDTH),
                              smm & (mm_["Ds_mass"] > M_DS + SR_HALFWIDTH), ww)
                r, e = double_ratio(c_lo, c_hi)
                cells.append("n/a" if np.isnan(r) else f"{r:.3f}+-{e:.3f}")
                if cname == "OS":
                    n_os_lo = int(c_lo[0])
            lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
            band = f"{plo:.0f}-{phi:.0f}" if phi < 1e8 else f">{plo:.0f}"
            print(f"{lbl:>8s} {band:>14s} {cells[0]:>18s} {cells[1]:>18s} {n_os_lo:>10,}")

    # ---------- the D+- peak: a real-charm probe on the OS side ----------
    # OS is where the analysis lives and its signal region is blinded, so the only
    # real-charm peak that can be measured there is D+- -> K+K-pi+- at 1.8696, which
    # sits in the low sideband. Comparing the data and MC D+- EXCESS (not the window
    # yield, which is mostly continuum) says whether W+jets MC has the right amount of
    # genuine, charge-correlated charm -- the thing an MVA trained on SS data has
    # never been shown how to reject.
    #
    # Background under the peak comes from a linear extrapolation of TWO windows both
    # BELOW the peak. The usual symmetric local-sideband subtraction cannot be used
    # here: its upper strip would run to ~1.94 and overlap the OS blind region, which
    # depletes the data side only and inflates the apparent excess. That exact mistake
    # was made once already on this measurement.
    M_D = 1.8696
    PK = (M_D - 0.02, M_D + 0.02)
    # TWO window schemes, deliberately. Scheme "near" reaches up to 1.8196, which is
    # inside the region where partially reconstructed charm (Ds->phi rho+, phi l nu)
    # piles up -- i.e. inside the very structure this study suspects is growing with
    # the cut. If it is, that scheme's extrapolated background is biased high and the
    # excess biased low, by an amount that itself grows with the cut. Scheme "far"
    # backs both windows away from it. The two agreeing within errors is the licence
    # to quote the number at all; disagreeing means the D+- measurement is not usable
    # and must be reported as such rather than quietly averaged.
    SCHEMES = {"near": ((1.7196, 1.7696), (1.7696, 1.8196)),
               "far":  ((1.6796, 1.7196), (1.7196, 1.7696))}

    def peak_excess(mass, sba, sbb, w=None):
        """(excess, error, extrapolated density, measured peak density).

        The two densities are returned so a nonsense extrapolation is visible in the
        printout instead of silently producing a plausible-looking excess.
        """
        def integ(lo, hi):
            m = (mass >= lo) & (mass < hi)
            if w is None:
                n = float(m.sum())
                return n, n            # Poisson: variance = N
            return float(w[m].sum()), float((w[m] ** 2).sum())
        (na, va), (nb, vb) = integ(*sba), integ(*sbb)
        npk, vpk = integ(*PK)
        wa, wb, wpk = sba[1] - sba[0], sbb[1] - sbb[0], PK[1] - PK[0]
        ca, cb = 0.5 * sum(sba), 0.5 * sum(sbb)
        da, db = na / wa, nb / wb
        slope = (db - da) / (cb - ca)
        dpk = db + slope * (M_D - cb)
        bkg = dpk * wpk
        # d(bkg)/d(na), d(bkg)/d(nb) -- the three windows are disjoint, no covariance.
        gb = wpk * (1.0 / wb) * (1.0 + (M_D - cb) / (cb - ca))
        ga = wpk * (1.0 / wa) * (-(M_D - cb) / (cb - ca))
        err = np.sqrt(vpk + ga ** 2 * va + gb ** 2 * vb)
        return npk - bkg, err, dpk, npk / wpk

    # MEASURED 2026-08-11 AND FOUND UNUSABLE -- kept because the failure is the
    # result. The two schemes agree only with no cut (2.42+-0.63 vs 1.89+-0.48) and
    # diverge badly once cut (eps_sig=50%: 0.67+-0.36 vs 0.22+-0.37), while the
    # extrapolated background falls from 89% of the peak-window density to ~40%. A
    # linear extrapolation across ~100 MeV is not reliable once the cut reshapes the
    # low-mass spectrum, which is exactly what this study found it does. Do NOT quote
    # the D+- excess ratio. Use the 1.83-1.91 window ratio from the table above
    # instead -- it needs no extrapolation and answers the same question.
    print(f"\n--- OS D+- excess, data vs W+jets, {common.resolve(k)['label']} model ---")
    print(f"    NOT USABLE -- the two schemes disagree once cut; see comment in source.")
    print(f"    peak {PK[0]:.4f}-{PK[1]:.4f}; background extrapolated from two windows")
    print(f"    BELOW the peak (the symmetric scheme's upper strip would overlap the "
          f"blind region)")
    for sname, (sba, sbb) in SCHEMES.items():
        print(f"      scheme '{sname}': {sba} and {sbb}")
    print(f"\n{'eff_sig':>8s} " + " ".join(
        f"{c:>17s}" for c in ("data exc [near]", "MC exc [near]", "ratio [near]",
                              "ratio [far]", "bkg/peak dens"))
          + f" {'CR ratio':>9s}")
    for eff in SCAN:
        thr = -1.0 if eff >= 1.0 else np.quantile(sig_sc, 1 - eff)
        sod, som = od["score_dsmu"] > thr, om["score_dsmu"] > thr
        res = {}
        for sname, (sba, sbb) in SCHEMES.items():
            xd, xd_e, dpk_d, dm_d = peak_excess(od["Ds_mass"][sod], sba, sbb)
            xm, xm_e, _, _ = peak_excess(om["Ds_mass"][som], sba, sbb, omw[som])
            if xm > 0 and xd > 0:
                rr = xd / xm
                re = rr * np.sqrt((xd_e / xd) ** 2 + (xm_e / xm) ** 2)
                rs = f"{rr:.2f}+-{re:.2f}"
            else:
                rs = "n/a"
            res[sname] = (xd, xd_e, xm, xm_e, rs, dpk_d / dm_d if dm_d > 0 else np.nan)
        xd, xd_e, xm, xm_e, rs_n, frac = res["near"]
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        cr_r = yield_ratio(*counts(sod & od_sb, som & om_sb, omw))[0]
        print(f"{lbl:>8s} {f'{xd:,.0f}+-{xd_e:,.0f}':>17s} "
              f"{f'{xm:,.0f}+-{xm_e:,.0f}':>17s} {rs_n:>17s} {res['far'][4]:>17s} "
              f"{frac:>17.2f} {cr_r:9.2f}")

    print(f"\n--- low vs high sideband, {common.resolve(k)['label']} model ---")
    print(f"{'eff_sig':>8s} " + " ".join(
        f"{c:>16s}" for c in ("OS low", "OS high", "OS low/high",
                              "SS low", "SS high", "SS low/high")))
    for eff, _, _, _, _, _, lohi, _ in scan:
        f3 = lambda r: "n/a" if np.isnan(r[0]) else f"{r[0]:.3f}+-{r[1]:.3f}"
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        print(f"{lbl:>8s} " + " ".join(f"{f3(v):>16s}"
                                       for nm_ in ("OS", "SS") for v in lohi[nm_]))

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

    # ---------- the two normalization-free ratios ----------
    # Left: SR relative to its own sideband -- is there an SR-specific excess, and
    # does it grow? Right: low sideband relative to high -- is there a tilt?
    # Both are double ratios, so the overall data/MC level cancels out of them and
    # what is left is purely the mass dependence. A shaded band marks where the MC
    # side falls below 30 effective entries, past which the points are noise and the
    # curve must not be read for structure.
    fig, axs = plt.subplots(1, 2, figsize=(20, 8))
    xs = [100 * s[0] for s in scan]
    NEFF_FLOOR = 30.0
    lo_x = [x for x, s in zip(xs, scan) if s[7] < NEFF_FLOOR]

    axs[0].errorbar(xs, [s[5][0] for s in scan], yerr=[s[5][1] for s in scan],
                    marker="s", markersize=8, lw=2, capsize=3, color="firebrick",
                    label="SS  SR / CR")
    axs[0].set_ylabel("(data/MC in SR) / (data/MC in CR)")
    axs[0].set_title(f"SR-specific excess -- {common.resolve(k)['label']} model", fontsize=15)

    for nm_, c, mk in (("OS", "steelblue", "o"), ("SS", "seagreen", "^")):
        axs[1].errorbar(xs, [s[6][nm_][2][0] for s in scan],
                        yerr=[s[6][nm_][2][1] for s in scan],
                        marker=mk, markersize=8, lw=2, capsize=3, color=c,
                        label=f"{nm_}  low sideband / high sideband")
    axs[1].set_ylabel("(data/MC below $M_{D_s}$) / (data/MC above)")
    axs[1].set_title(f"mass tilt inside the sideband -- "
                     f"{common.resolve(k)['label']} model", fontsize=15)

    for ax in axs:
        ax.axhline(1.0, color="black", lw=1.5, ls=":")
        if lo_x:
            ax.axvspan(min(lo_x), max(lo_x) + 2.5, color="grey", alpha=0.18)
            ax.text(max(lo_x) + 2.0, ax.get_ylim()[1] * 0.95,
                    f"MC $N_{{eff}}<{NEFF_FLOOR:.0f}$", fontsize=11, color="dimgrey",
                    ha="right", va="top")
        ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
        ax.invert_xaxis()
        ax.legend(fontsize=13)
    hep.cms.label(ax=axs[0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                  fontsize=11)
    plt.tight_layout()
    out = f"{OUTDIR}/wjets_double_ratio_{k}.png"
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

    # ---------- the same thing for SS, where the signal region is visible ----------
    # SS is not blinded, so unlike the OS version above this shows the FULL mass range
    # including the Ds window, and quotes the SR ratio next to the CR one. That is the
    # only place the signal-region closure can be looked at directly at a tight
    # working point.
    #
    # Both SS curves sit at half luminosity: the data is the odd event_id half (the
    # ssdata model trained on the even one) and the MC is scaled by the same measured
    # fraction to match. The ratio is therefore correct; the absolute heights are half
    # the physical numbers, which the legend says explicitly so the N values are not
    # compared against the OS figure's.
    sd_d, mm_d, mmw = keep[k]["ss_d"], keep[k]["ss_m"], keep[k]["ss_mw"]
    ssd = sd_d["score_dsmu"] > thr
    ssm = mm_d["score_dsmu"] > thr
    hd_s, _ = np.histogram(sd_d["Ds_mass"][ssd], bins=mbins)
    hm_s, _ = np.histogram(mm_d["Ds_mass"][ssm], bins=mbins, weights=mmw[ssm])
    hm2_s, _ = np.histogram(mm_d["Ds_mass"][ssm], bins=mbins, weights=mmw[ssm] ** 2)
    r_cr, e_cr = yield_ratio(float(hd_s[SB_BIN].sum()), float(hm_s[SB_BIN].sum()),
                             float(hm2_s[SB_BIN].sum()))
    r_sr, e_sr = yield_ratio(float(hd_s[~SB_BIN].sum()), float(hm_s[~SB_BIN].sum()),
                             float(hm2_s[~SB_BIN].sum()))

    fig, ax = plt.subplots(figsize=(11, 8))
    hep.histplot(hm_s, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                 label=f"SS W+jets MC (lumi-wt. $\\times$ odd frac.)  N={hm_s.sum():,.1f}")
    hep.histplot(hd_s, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd_s),
                 color="black", marker="o", markersize=6, elinewidth=1.5,
                 label=f"SS data (odd event_id half)  N={hd_s.sum():,.0f}")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel("Events / 20 MeV")
    ax.set_title(f"{common.resolve(k)['label']} model -- "
                 rf"$\varepsilon_{{sig}}$={100*EFF_TIGHT:.0f}%  (SS, score>{thr:.4f})"
                 f"\ndata/MC:  CR {r_cr:.2f} $\\pm$ {e_cr:.2f}   "
                 f"SR {r_sr:.2f} $\\pm$ {e_sr:.2f}", fontsize=15)
    ax.legend(fontsize=12)
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
    plt.tight_layout()
    out = f"{OUTDIR}/ss_dist_eff20_{k}.png"
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"{common.resolve(k)['label']} model, eff_sig=20%: "
          f"SS CR data={hd_s[SB_BIN].sum():,.0f} W+jets={hm_s[SB_BIN].sum():,.1f} "
          f"ratio={r_cr:.2f}+-{e_cr:.2f} | "
          f"SS SR data={hd_s[~SB_BIN].sum():,.0f} W+jets={hm_s[~SB_BIN].sum():,.1f} "
          f"ratio={r_sr:.2f}+-{e_sr:.2f}")
    print(f"Saved {out}")

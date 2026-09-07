"""What does the CHARGE-CORRELATED background do under the DNN cut?

    python3 check_os_minus_ss.py

The four-window study showed OS and SS drifting apart in M(KKpi) as the cut tightens,
in opposite directions, with no mechanism identified. Slicing the same two spectra
into more windows only redescribes that. This takes the one step the data still
allows: subtract them.

The pi charge requirement is the only difference between the OS and SS selections.
Background that has no reason to correlate the pi charge with the muon charge --
combinatorics, and charge-symmetric real charm -- populates both equally and cancels
in OS - SS. What survives the subtraction is the charge-correlated part (W+c and
friends, where the charm hadron's charge is tied to the W's). So OS - SS isolates
exactly the component present in OS and absent in SS, which is what the four-window
result pointed at without being able to name.

Two questions, both asked of data and of lumi-weighted W+jets separately and then
compared:
  1. how much charge-correlated background is there, and does the cut enrich it?
  2. what is its M(KKpi) shape, and does MC get that shape right?

Bookkeeping that has to be right or the subtraction is meaningless:

* SS data is the odd-event_id half (the ssdata model trained on the even half), so it
  is scaled to full size by the MEASURED odd fraction, and its variance scales by the
  square of that factor. SS W+jets is used at full weight -- it was never split.
* OS data is blinded in |M-M_Ds| < 40 MeV, so the difference is undefined there. Those
  bins are masked for BOTH data and MC rather than silently compared.
* This is NOT the OS-SS subtraction that would isolate a signal. It isolates a
  charge-correlated BACKGROUND, which is its intended use here.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_os_minus_ss.py
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

MODEL = "ssdata"
PRED = f"{common.resolve(MODEL)['out_dir']}/Prediction"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
EFFS = [1.0, 0.90, 0.70, 0.50]
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      common.B_VETO, "event_id"]

MBINS = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)
MCENT = 0.5 * (MBINS[1:] + MBINS[:-1])
CR_BIN = np.abs(MCENT - M_DS) >= SR_HALFWIDTH      # bins where OS data exists

WINDOWS = [("low  <1.830", 1.66835, 1.830),
           ("D+- 1.83-1.91", 1.830, 1.910),
           ("1.91-1.929", 1.910, M_DS - SR_HALFWIDTH),
           ("high >2.008", M_DS + SR_HALFWIDTH, 2.26835)]


def load(name, *, ss=False, sig=False, parity=None, era=None):
    br = list(dict.fromkeys(BR + ["Lumi_reweight"] + (["era"] if era else [])
                            + (["Ds_match", "mu_match"] if sig else [])))
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if parity is not None:
        m &= (a["event_id"] % 2) == parity
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


print(f"Loading from {PRED} ...")
sig = np.concatenate([load(n, sig=True)["score_dsmu"] for n in common.SIGNALS])
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig, 1 - e))) for e in EFFS}

# preEE only (era=1), as before.
os_d = load("Data", era=1)
ss_d_all = load("Data_SS", ss=True, era=1)
ss_d = load("Data_SS", ss=True, parity=1, era=1)
F_ODD = len(ss_d["Ds_mass"]) / len(ss_d_all["Ds_mass"])
SS_SCALE = 1.0 / F_ODD


def wjets(ss):
    # Per-event Lumi_reweight from the file; preEE only, as before.
    d = load("WtoLNu2Jets_SS" if ss else "WtoLNu2Jets", ss=ss, era=1)
    return d, d["Lumi_reweight"]


os_m, os_mw = wjets(False)
ss_m, ss_mw = wjets(True)
print(f"  OS data {len(os_d['Ds_mass']):>9,} | SS data odd {len(ss_d['Ds_mass']):>9,} "
      f"(x{SS_SCALE:.4f} to full) | OS MC {os_mw.sum():>11,.0f} | SS MC {ss_mw.sum():>11,.0f}")


def hist(d, cut, lo=None, hi=None, w=None):
    """(sum, sum of variances) in a mass window, or binned if lo is None."""
    m = cut if lo is None else cut & (d["Ds_mass"] >= lo) & (d["Ds_mass"] < hi)
    if w is None:
        if lo is None:
            h, _ = np.histogram(d["Ds_mass"][m], bins=MBINS)
            return h.astype(float), h.astype(float)
        n = float(m.sum())
        return n, n
    if lo is None:
        h, _ = np.histogram(d["Ds_mass"][m], bins=MBINS, weights=w[m])
        v, _ = np.histogram(d["Ds_mass"][m], bins=MBINS, weights=w[m] ** 2)
        return h, v
    return float(w[m].sum()), float((w[m] ** 2).sum())


def diff(eff, lo=None, hi=None):
    """(D_data, err, D_MC, err) -- the charge-correlated component, OS minus SS."""
    thr = THR[eff]
    od = hist(os_d, os_d["score_dsmu"] > thr, lo, hi)
    sd = hist(ss_d, ss_d["score_dsmu"] > thr, lo, hi)
    om = hist(os_m, os_m["score_dsmu"] > thr, lo, hi, os_mw)
    sm = hist(ss_m, ss_m["score_dsmu"] > thr, lo, hi, ss_mw)
    # SS data is half the sample; scaling the yield by SS_SCALE scales its variance
    # by SS_SCALE^2. OS and SS are disjoint event sets, so no covariance.
    dd = od[0] - SS_SCALE * sd[0]
    dd_e = np.sqrt(od[1] + SS_SCALE ** 2 * sd[1])
    dm = om[0] - sm[0]
    dm_e = np.sqrt(om[1] + sm[1])
    return dd, dd_e, dm, dm_e


# ---------- 1. how big is it, and does the cut enrich it? ----------
# f_cc = (OS - SS) / OS is the charge-correlated FRACTION of the OS sample. If the
# classifier is picking out genuine charm, this fraction has to rise with the cut --
# and whether MC's rises the same way is the first thing that would show a problem.
print("\n=== the charge-correlated component (OS - SS), control region ===")
print(f"{'eff_sig':>8s} {'D(data)':>16s} {'D(W+jets)':>16s} {'D data/MC':>15s} "
      f"{'f_cc data':>10s} {'f_cc MC':>9s}")
for eff in EFFS:
    thr = THR[eff]
    tot = []
    for lo, hi in ((1.66835, M_DS - SR_HALFWIDTH), (M_DS + SR_HALFWIDTH, 2.26835)):
        tot.append(diff(eff, lo, hi))
    dd = sum(t[0] for t in tot)
    dd_e = np.sqrt(sum(t[1] ** 2 for t in tot))
    dm = sum(t[2] for t in tot)
    dm_e = np.sqrt(sum(t[3] ** 2 for t in tot))
    n_os_d = float(((os_d["score_dsmu"] > thr)
                    & (np.abs(os_d["Ds_mass"] - M_DS) >= SR_HALFWIDTH)).sum())
    sel_m = (os_m["score_dsmu"] > thr) & (np.abs(os_m["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    n_os_m = float(os_mw[sel_m].sum())
    r = dd / dm if dm > 0 else np.nan
    r_e = abs(r) * np.sqrt((dd_e / dd) ** 2 + (dm_e / dm) ** 2) if dm > 0 and dd != 0 else np.nan
    lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
    print(f"{lbl:>8s} {f'{dd:,.0f}+-{dd_e:,.0f}':>16s} {f'{dm:,.0f}+-{dm_e:,.0f}':>16s} "
          f"{f'{r:.3f}+-{r_e:.3f}':>15s} {dd/n_os_d:10.3f} {dm/n_os_m:9.3f}")

# ---------- 2. what is its mass shape, and does MC get it right? ----------
print("\n=== mass shape of the charge-correlated component ===")
print("    D = OS - SS in each window; last two columns are its data/MC ratio\n")
print(f"{'eff_sig':>8s} " + " ".join(f"{n:>17s}" for n, _, _ in WINDOWS))
for eff in EFFS:
    cells = []
    for _, lo, hi in WINDOWS:
        dd, dd_e, dm, dm_e = diff(eff, lo, hi)
        if dm <= 0 or dd == 0:
            cells.append("n/a")
            continue
        r = dd / dm
        r_e = abs(r) * np.sqrt((dd_e / dd) ** 2 + (dm_e / dm) ** 2)
        cells.append(f"{r:.2f}+-{r_e:.2f}")
    lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
    print(f"{lbl:>8s} " + " ".join(f"{c:>17s}" for c in cells))

print("\n    raw yields behind those ratios (data | W+jets):")
print(f"{'eff_sig':>8s} " + " ".join(f"{n:>17s}" for n, _, _ in WINDOWS))
for eff in EFFS:
    cells = []
    for _, lo, hi in WINDOWS:
        dd, dd_e, dm, dm_e = diff(eff, lo, hi)
        cells.append(f"{dd:,.0f}|{dm:,.0f}")
    lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
    print(f"{lbl:>8s} " + " ".join(f"{c:>17s}" for c in cells))

# ---------- plot ----------
fig, axes = plt.subplots(2, len(EFFS), figsize=(7.5 * len(EFFS), 12),
                         gridspec_kw={"height_ratios": [2, 1]})
for col, eff in enumerate(EFFS):
    thr = THR[eff]
    od = hist(os_d, os_d["score_dsmu"] > thr)
    sd = hist(ss_d, ss_d["score_dsmu"] > thr)
    om = hist(os_m, os_m["score_dsmu"] > thr, w=os_mw)
    sm = hist(ss_m, ss_m["score_dsmu"] > thr, w=ss_mw)
    dd = od[0] - SS_SCALE * sd[0]
    dd_e = np.sqrt(od[1] + SS_SCALE ** 2 * sd[1])
    dm = om[0] - sm[0]
    dm_e = np.sqrt(om[1] + sm[1])
    # The blind window has no OS data, so the difference does not exist there --
    # blanked for MC too so the two curves are always compared on the same footing.
    dd, dd_e, dm, dm_e = (np.where(CR_BIN, x, np.nan) for x in (dd, dd_e, dm, dm_e))

    ax = axes[0][col]
    hep.histplot(dm, MBINS, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                 label=f"W+jets MC  N={np.nansum(dm):,.0f}")
    ax.errorbar(MCENT, dd, yerr=dd_e, fmt="o", markersize=5, color="black",
                elinewidth=1.4, label=f"data  N={np.nansum(dd):,.0f}")
    ax.axhline(0.0, color="grey", lw=1, ls=":")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.set_title("OS $-$ SS  " + ("no cut" if eff >= 1.0
                                  else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"), fontsize=15)
    if col == 0:
        ax.set_ylabel("charge-correlated events / 20 MeV")
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=11)
    ax.legend(fontsize=11)

    ax = axes[1][col]
    ok = np.isfinite(dd) & np.isfinite(dm) & (dm > 0) & (dd > 0)
    r = np.full(len(dd), np.nan)
    re = np.full(len(dd), np.nan)
    r[ok] = dd[ok] / dm[ok]
    re[ok] = r[ok] * np.sqrt((dd_e[ok] / dd[ok]) ** 2 + (dm_e[ok] / dm[ok]) ** 2)
    ax.errorbar(MCENT, r, yerr=re, fmt="o", markersize=5, color="black", elinewidth=1.3)
    ax.axhline(1.0, color="firebrick", lw=1.3, ls=":")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.set_ylim(0, 3)
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        ax.set_ylabel("data / MC")
plt.tight_layout()
out = common.fig_path("os_ss", "os_minus_ss.png", key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"\nSaved {out}")

"""DNN score response, OS vs SS, shape-normalized, for data and for W+jets MC.

    python3 draw_score_os_ss.py [model]        # default: ssdata

Two figures, each a unit-normalized score spectrum with its OS/SS ratio underneath:

  1. data      -- control region only (|M-M_Ds| >= 40 MeV); OS is blinded in the
                  signal region, so an OS/SS ratio built there would be meaningless
  2. W+jets MC -- full mass range, it is not blinded

**Every curve is normalized to unit area**, so this asks about SHAPE, not yield. OS
already has ~1.14x more events than SS before any cut, and leaving that in would put a
flat factor 1.14 into the ratio and bury what is being looked for. Normalized, the
ratio reads directly: 1.0 means OS and SS put the same FRACTION of their events at
that score, and a rise towards high score means OS concentrates more of itself in the
signal-like region than SS does.

Why that is worth asking: the classifier was trained against SS data, so it has been
shown what SS-like background looks like and has never been shown the charge-correlated
component that only OS contains. Comparing data's OS/SS shape ratio against the MC's
says whether the MC contains the same charge asymmetry in its score response.

One thing that is not optional: SS data is the odd event_id half. The model trained on
the even half, so those events carry artificially low scores and mixing them in would
distort the very distribution being drawn. Under unit normalization the odd/even scale
factor cancels, so no scaling is applied or needed.

The main panel is log-y: normalized, the lowest score bin holds most of the sample and
the rest spans several decades, so a linear axis would show one spike and nothing
else. The ratio panel is linear.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_score_os_ss.py
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

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
MODEL = ARGS[0] if ARGS else "ssdata"

# --veto-peaks removes BOTH resonance windows, Ds and D+-, from every sample.
#
# The point of it: the OS/SS score asymmetry rising into the tail is naturally read as
# the classifier picking out charge-correlated real charm, which OS has and SS does
# not. If that is the whole story, deleting the two mass peaks should collapse the
# ratio to a flat 1. If it does not, there is charge-correlated NON-peaking background
# as well, and that is a different problem.
#
# It also makes data and MC strictly comparable for the first time in this script: by
# default the W+jets curve spans the full mass range while OS data is missing its Ds
# window to blinding, so the two are not integrated over the same mass. Under
# --veto-peaks every sample gets the identical mass mask.
VETO_PEAKS = "--veto-peaks" in sys.argv

PRED = f"{common.resolve(MODEL)['out_dir']}/Prediction"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D, D_HALFWIDTH = 1.8696, 0.04     # same half-width as the Ds window, deliberately
SBINS = np.linspace(0.0, 1.0, 51)
SCENT = 0.5 * (SBINS[1:] + SBINS[:-1])

BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      common.B_VETO, "event_id"]


def load(name, *, ss=False, parity=None, cr=False):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(BR, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
    if cr or VETO_PEAKS:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    if VETO_PEAKS:
        m &= np.abs(a["Ds_mass"] - M_D) >= D_HALFWIDTH
    if parity is not None:
        m &= (a["event_id"] % 2) == parity
    return {k: v[m] for k, v in a.items()}


print(f"Loading from {PRED} ...")
d_os = load("Data_2022", cr=True)["score_dsmu"]
d_ss = load("Data_2022_SS", ss=True, cr=True, parity=1)["score_dsmu"]


def wjets(ss):
    parts = [load(f"{n}_SS" if ss else n, ss=ss) for n in common.WJETS_WEIGHTS]
    v = np.concatenate([p["score_dsmu"] for p in parts])
    w = np.concatenate([np.full(len(p["score_dsmu"]), wt)
                        for p, wt in zip(parts, common.WJETS_WEIGHTS.values())])
    return v, w


m_os, m_os_w = wjets(False)
m_ss, m_ss_w = wjets(True)
print(f"  data  OS(CR) {len(d_os):>9,} | SS(CR, odd half) {len(d_ss):>9,}")
print(f"  W+jets OS    {m_os_w.sum():>9,.0f} | SS               {m_ss_w.sum():>9,.0f}")


# ---------- which peak, if any, drives the OS/SS rise ----------
# Printed unconditionally because the naive version of this comparison is a trap that
# has already been fallen into once. The default data selection ALREADY has the Ds
# window removed (OS is blinded there, and the same cut is applied to SS to keep the
# two symmetric), while the default W+jets selection spans the full mass range. So
# running --veto-peaks and diffing removes D+- ONLY from data but BOTH peaks from MC,
# and the resulting "data is flat, MC drops" reads as a data/MC difference when it is
# just two different comparisons placed side by side.
#
# Broken out per window, the actual answer is: D+- contributes nothing on either side,
# and the Ds window's contribution is measurable in MC only -- in data the OS
# numerator is blinded there, so any row that keeps the Ds window is comparing a
# blinded OS against an unblinded SS and is meaningless.
def _os_ss(vo, wo, mo, vs, ws, ms, keep, thr):
    ko, ks = keep(mo), keep(ms)
    if wo is None:
        no, to = float((vo[ko] > thr).sum()), float(ko.sum())
        fo, varo = no / to, no / to ** 2
        ns, ts = float((vs[ks] > thr).sum()), float(ks.sum())
        fs, vars_ = ns / ts, ns / ts ** 2
    else:
        so, sso = wo[ko], ws[ks]
        io, is_ = vo[ko] > thr, vs[ks] > thr
        fo, varo = so[io].sum() / so.sum(), (so[io] ** 2).sum() / so.sum() ** 2
        fs, vars_ = sso[is_].sum() / sso.sum(), (sso[is_] ** 2).sum() / sso.sum() ** 2
    r = fo / fs
    return r, r * np.sqrt(varo / fo ** 2 + vars_ / fs ** 2)


REGIONS = [
    ("full range", lambda m: np.ones(len(m), bool), "data row meaningless (OS blinded)"),
    ("Ds removed", lambda m: np.abs(m - M_DS) >= SR_HALFWIDTH, ""),
    ("D+- removed", lambda m: np.abs(m - M_D) >= D_HALFWIDTH,
     "data row meaningless (OS blinded)"),
    ("both removed", lambda m: (np.abs(m - M_DS) >= SR_HALFWIDTH)
     & (np.abs(m - M_D) >= D_HALFWIDTH), ""),
]

if not VETO_PEAKS:
    d_os_m = load("Data_2022")                      # unrestricted mass, for the scan
    d_ss_m = load("Data_2022_SS", ss=True, parity=1)
    m_os_all = [load(n) for n in common.WJETS_WEIGHTS]
    m_ss_all = [load(f"{n}_SS", ss=True) for n in common.WJETS_WEIGHTS]
    mow = np.concatenate([np.full(len(p["Ds_mass"]), w)
                          for p, w in zip(m_os_all, common.WJETS_WEIGHTS.values())])
    msw = np.concatenate([np.full(len(p["Ds_mass"]), w)
                          for p, w in zip(m_ss_all, common.WJETS_WEIGHTS.values())])
    SETS = [("data", d_os_m["score_dsmu"], None, d_os_m["Ds_mass"],
             d_ss_m["score_dsmu"], None, d_ss_m["Ds_mass"]),
            ("W+jets", np.concatenate([p["score_dsmu"] for p in m_os_all]), mow,
             np.concatenate([p["Ds_mass"] for p in m_os_all]),
             np.concatenate([p["score_dsmu"] for p in m_ss_all]), msw,
             np.concatenate([p["Ds_mass"] for p in m_ss_all]))]

    print("\n=== OS/SS above threshold, per mass region -- which peak drives the rise? ===")
    print(f"{'sample':>8s} {'mass region':>14s} {'>0.9':>16s} {'>0.98':>16s}   note")
    for name, vo, wo, mo, vs, ws, ms in SETS:
        for rname, keep, note in REGIONS:
            cells = [f"{r:.3f}+-{e:.3f}" for r, e in
                     (_os_ss(vo, wo, mo, vs, ws, ms, keep, t) for t in (0.9, 0.98))]
            print(f"{name:>8s} {rname:>14s} {cells[0]:>16s} {cells[1]:>16s}   "
                  f"{note if name == 'data' else ''}")

MREG = (r"both peaks removed: $|M-M_{D_s}|>40$, $|M-M_{D^\pm}|>40$ MeV"
        if VETO_PEAKS else None)
FIGS = [
    ("data", "Data -- " + (MREG or "control region"),
     (d_os, None), "Data OS (1-tight-mu + b-veto)",
     (d_ss, None), "Data SS (odd event_id half)"),
    ("wjets", "W+jets MC -- " + (MREG or "full mass range"),
     (m_os, m_os_w), "W+jets MC OS (lumi-weighted)",
     (m_ss, m_ss_w), "W+jets MC SS (lumi-weighted)"),
]


def norm_hist(v, w):
    """Unit-normalized density and its error."""
    if w is None:
        h, _ = np.histogram(v, bins=SBINS)
        h = h.astype(float)
        var, tot = h.copy(), h.sum()
    else:
        h, _ = np.histogram(v, bins=SBINS, weights=w)
        var, _ = np.histogram(v, bins=SBINS, weights=w ** 2)
        tot = w.sum()
    width = np.diff(SBINS)
    return h / tot / width, np.sqrt(var) / tot / width


for stem, title, (vo, wo), lab_os, (vs, ws), lab_ss in FIGS:
    ho, eo = norm_hist(vo, wo)
    hs, es = norm_hist(vs, ws)

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1]})
    ax.errorbar(SCENT, ho, yerr=eo, fmt="o", markersize=4.5, color="steelblue",
                elinewidth=1.2, label=lab_os)
    ax.errorbar(SCENT, hs, yerr=es, fmt="s", markersize=4.5, color="firebrick",
                elinewidth=1.2, label=lab_ss)
    ax.set_yscale("log")
    ax.set_ylabel("normalized to unit area")
    ax.set_title(f"{title} -- {common.resolve(MODEL)['label']} model", fontsize=16)
    ax.legend(fontsize=12, loc="upper center")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                  fontsize=11)

    ok = (ho > 0) & (hs > 0)
    r = np.full(len(ho), np.nan)
    re = np.full(len(ho), np.nan)
    r[ok] = ho[ok] / hs[ok]
    re[ok] = r[ok] * np.sqrt((eo[ok] / ho[ok]) ** 2 + (es[ok] / hs[ok]) ** 2)
    rax.errorbar(SCENT, r, yerr=re, fmt="o", markersize=4.5, color="black",
                 elinewidth=1.2)
    rax.axhline(1.0, color="firebrick", lw=1.3, ls=":")
    rax.set_xlabel(f"DNN score  (trained vs {common.resolve(MODEL)['label']})")
    rax.set_ylabel("OS / SS  (shape)")
    lo_y, hi_y = np.nanpercentile(r, [2, 98])
    pad = 0.35 * (hi_y - lo_y) if hi_y > lo_y else 0.2
    rax.set_ylim(max(0.0, lo_y - pad), hi_y + pad)

    plt.tight_layout()
    out = f"{OUTDIR}/score_os_ss_{stem}_{MODEL}{'_nopeaks' if VETO_PEAKS else ''}.png"
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out}")

    # A few numbers to quote alongside the figure: the shape ratio in the highest
    # score bins, where the per-bin points get noisy and the eye is least reliable.
    def frac(v, w, lo):
        """(fraction above lo, its variance) -- yields divided out, so this is shape."""
        i = v > lo
        if w is None:
            n, tot = float(i.sum()), float(len(v))
            return n / tot, n / tot ** 2
        n, tot = float(w[i].sum()), float(w.sum())
        return n / tot, float((w[i] ** 2).sum()) / tot ** 2

    for lo in (0.9, 0.98):
        fo, vfo = frac(vo, wo, lo)
        fs, vfs = frac(vs, ws, lo)
        r = fo / fs
        re = r * np.sqrt(vfo / fo ** 2 + vfs / fs ** 2)
        print(f"    score>{lo}:  OS {100*fo:.3f}%  SS {100*fs:.3f}%"
              f"   OS/SS = {r:.3f} +- {re:.3f}")

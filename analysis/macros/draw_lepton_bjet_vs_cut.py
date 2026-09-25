"""Lepton and b-jet multiplicities at matched signal efficiency.

    python3 draw_lepton_bjet_vs_cut.py [model]      # default: ssdata

Same layout as compare_models_mass.png -- one column per working point, absolute
yields, linear -- with M(KKpi) replaced by the multiplicity counters. Three figures:

  mult_muons      num_{os,ss}_mu_{loose,medium}
  mult_electrons  num_{os,ss}_e_{cutloose,cuttight,mvawp80,mvawp90}
  mult_bjets      n_b_jets_{tight,medium}_dR04

Rows are the variables, columns the working points.

Two variables are deliberately absent and two selections are deliberately different,
both for the same reason -- a counter that the selection fixes carries no information:

* `num_os_mu_tight` and `num_ss_mu_tight` are not plotted at all. NOTE (2026-09-22):
  the stated reason -- "the analysis selection IS num_os_mu_tight==1 &
  num_ss_mu_tight==0, so they would be single bars" -- EXPIRED on 2026-09-21, when that
  condition was replaced by common.selection_mask (corrected muon cuts). These two
  counters are no longer fixed by the selection and would now vary, so whether they are
  worth a row again is an open question, not a settled one. Loose and medium are
  supersets of tight and do vary.
* The b-jet figure has the **b-veto removed** (`n_b_jets_tight_dR04==0` is the veto
  itself, so with it applied the tight row would be a spike at zero and the medium row
  would be biased by it). Muons and electrons keep the full selection.

Yields are absolute: the lumi-weighted W+jets is a sigma x L prediction with no free
normalization, and OS data is the observed count, so the two are directly comparable
and neither is rescaled. Signal has no meaningful absolute normalization here, so it is
scaled to the data yield of its panel with the factor written into the legend -- it is
there for its SHAPE only.

Data is the OS control region (the signal region is blinded, and the same mass cut is
applied to W+jets and signal so all three are integrated over the same range). Signal
is truth-matched (Ds_match & mu_match).

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_lepton_bjet_vs_cut.py
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

MODEL = sys.argv[1] if len(sys.argv) > 1 else "ssdata"
PRED = f"{common.resolve(MODEL)['out_dir']}/Prediction"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
EFFS = [1.0, 0.90, 0.70, 0.50]
B_VETO = "n_b_jets_tight_dR04"

MU_VARS = [f"num_{c}_mu_{wp}" for wp in ("loose", "medium") for c in ("os", "ss")]
E_VARS = [f"num_{c}_e_{wp}" for wp in ("cutloose", "cuttight", "mvawp80", "mvawp90")
          for c in ("os", "ss")]
B_VARS = ["n_b_jets_tight_dR04", "n_b_jets_medium_dR04"]
BR = (["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight"]
      + MU_VARS + E_VARS + B_VARS)


def load(name, *, drop_bveto, sig=False):
    br = list(dict.fromkeys(BR + (["Ds_match", "mu_match"] if sig else [])))
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(br, library="np")
    m = (common.selection_mask(a, bveto=not drop_bveto) & (a["score_dsmu"] > -1)
         & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH))
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


def gather(drop_bveto):
    sig_parts = [load(n, drop_bveto=drop_bveto, sig=True) for n in common.SIGNALS]
    sig = {k: np.concatenate([p[k] for p in sig_parts]) for k in BR}
    mc_parts = [load(n, drop_bveto=drop_bveto) for n in common.WJETS_WEIGHTS]
    mc = {k: np.concatenate([p[k] for p in mc_parts]) for k in BR}
    mcw = np.concatenate([np.full(len(p["Ds_mass"]), w)
                          for p, w in zip(mc_parts, common.WJETS_WEIGHTS.values())])
    return sig, mc, mcw, load("Data_2022", drop_bveto=drop_bveto)


print(f"Loading from {PRED} ...")
SETS = {d: gather(d) for d in (False, True)}
for d, tag in ((False, "full selection "), (True, "b-veto removed ")):
    s, m, w, dat = SETS[d]
    print(f"  {tag} signal {len(s['Ds_mass']):>8,} | W+jets {w.sum():>10,.0f} "
          f"| OS data {len(dat['Ds_mass']):>9,}")

# Full analysis selection everywhere, b-veto included. That makes
# `n_b_jets_tight_dR04` a single bar at zero -- it IS common.B_VETO, applied before
# scoring -- which is the honest picture of the analysis region rather than a
# convenient one. `n_b_jets_medium_dR04` still varies under the veto because
# medium is a superset of tight, so it is the b-jet counter the analysis actually sees.
FIGS = [("muons", MU_VARS, False), ("electrons", E_VARS, False),
        ("bjets", B_VARS, False)]

for stem, variables, drop in FIGS:
    sig, mc, mcw, dat = SETS[drop]
    thr = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig["score_dsmu"], 1 - e)))
           for e in EFFS}

    fig, axes = plt.subplots(len(variables), len(EFFS),
                             figsize=(6.5 * len(EFFS), 5.0 * len(variables)),
                             squeeze=False)
    print(f"\n=== {stem} ({'b-veto removed' if drop else 'full selection'}) ===")
    for row, var in enumerate(variables):
        hi = int(max(np.percentile(dat[var], 99.95), np.percentile(mc[var], 99.95), 2))
        edges = np.arange(0, hi + 2) - 0.5
        cent = 0.5 * (edges[1:] + edges[:-1])
        print(f"{var:>24s} " + " ".join(
            f"{('no cut' if e >= 1 else f'{100*e:.0f}%'):>26s}" for e in EFFS))
        cells = []
        for col, eff in enumerate(EFFS):
            ax = axes[row][col]
            t = thr[eff]
            sd, sm, ss_ = (dat["score_dsmu"] > t, mc["score_dsmu"] > t,
                           sig["score_dsmu"] > t)
            clip = lambda v: np.clip(v, edges[0], edges[-1] - 1e-6)
            hd, _ = np.histogram(clip(dat[var][sd]), bins=edges)
            hm, _ = np.histogram(clip(mc[var][sm]), bins=edges, weights=mcw[sm])
            hm2, _ = np.histogram(clip(mc[var][sm]), bins=edges, weights=mcw[sm] ** 2)
            hs, _ = np.histogram(clip(sig[var][ss_]), bins=edges)
            # Signal carries no absolute normalization in this build, so it is scaled
            # to the panel's data yield and the factor is stated. Data and W+jets are
            # NOT touched -- their absolute comparison is the point.
            f_sig = hd.sum() / hs.sum() if hs.sum() > 0 else 1.0

            hep.histplot(hm, edges, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                         label=f"W+jets MC  N={hm.sum():,.0f}")
            hep.histplot(hs * f_sig, edges, ax=ax, histtype="step", lw=2.5,
                         color="firebrick",
                         label=f"signal $\\times${f_sig:.3g}  N={hs.sum():,.0f}")
            hep.histplot(hd, edges, ax=ax, histtype="errorbar", yerr=np.sqrt(hd),
                         color="black", marker="o", markersize=5, elinewidth=1.4,
                         label=f"OS data  N={hd.sum():,.0f}")
            ax.set_xticks(np.arange(0, hi + 1))
            ax.set_ylim(0, None)
            ax.set_xlabel(var)
            if col == 0:
                ax.set_ylabel("Events")
            lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
            ax.set_title(lbl, fontsize=14)
            ax.legend(fontsize=9)
            if row == 0 and col == 0:
                hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True,
                              text="Preliminary", fontsize=10)
            # mean +- error, the one number that summarises a multiplicity
            md = dat[var][sd].mean() if sd.sum() else np.nan
            mm_ = (np.average(mc[var][sm], weights=mcw[sm]) if sm.sum() else np.nan)
            ms = sig[var][ss_].mean() if ss_.sum() else np.nan
            cells.append(f"{md:.3f}/{mm_:.3f}/{ms:.3f}")
        print(f"{'mean data/MC/sig':>24s} " + " ".join(f"{c:>26s}" for c in cells))

    plt.tight_layout()
    out = common.fig_path("closure", f"mult_{stem}_{MODEL}.png", key=MODEL)
    plt.savefig(out, dpi=100)
    plt.close(fig)
    print(f"Saved {out}")

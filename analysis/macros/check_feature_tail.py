"""How much of the data/MC drift after a DNN cut comes from single-feature mismodelling?

    python3 check_feature_tail.py

The DNN working points sit at background efficiencies of order 10^-3, i.e. the far
tail of a 43-dimensional distribution. A few-percent data/MC mismatch per feature
compounds there, so the classifier can move the data/MC ratio without any missing
process being involved. This asks how much of that is attributable to individual
features, by cutting on each one ALONE and watching the same quantity.

**The x axis is eps_MC, not eps_sig.** A single feature at eps_sig=50% still keeps
5-20% of the background; the DNN at eps_sig=50% keeps 0.055%. Comparing them at
matched signal efficiency compares two completely different depths into the tail and
answers nothing -- the drift may be a phenomenon that only exists at eps_MC~10^-3,
which no single feature reaches at any signal efficiency. Matching on eps_MC puts
every curve at the same background rejection and makes the comparison meaningful.

Read it as: at the eps_MC where the DNN sits, does any single feature -- cut to the
same background rejection -- reproduce a comparable eps_data/eps_MC? If one does,
that feature's mismodelling is a sufficient explanation and can be handled with a
reweight plus a systematic. If none does, the drift is genuinely multivariate or
comes from composition (a process W+jets MC does not contain), not from one variable.

Everything is measured in the OS control region under the analysis selection, so it
is directly comparable to check_data_vs_wjets.py.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_feature_tail.py
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
MODEL_DIR = common.resolve(MODEL)["out_dir"]
PRED = f"{MODEL_DIR}/Prediction"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
N_FEATURES = 6                      # how many of the top-importance features to scan

# Background efficiencies to compare everything at. The DNN's own working points
# (eps_sig 90/70/50%) land near 1e-2 / 2e-3 / 6e-4, so the grid brackets them.
EPS_MC_GRID = np.array([1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 5e-4])

SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "Ds_mass", "score_dsmu"]


def load(name, *, sig=False, era=None):
    """Scored sample + all features, OS analysis selection, control region only."""
    br = list(dict.fromkeys(SEL + common.FEATURE_NAMES + ["Lumi_reweight"]
                            + (["era"] if era else [])
                            + (["Ds_match", "mu_match"] if sig else [])))
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(br, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
         & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH))
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


print(f"Loading from {PRED} ...")
sig = {k: np.concatenate([load(n, sig=True)[k] for n in common.SIGNALS])
       for k in ["score_dsmu"] + common.FEATURE_NAMES}
# preEE only (era=1), as before. Per-event Lumi_reweight replaces the three constants.
dat = load("Data", era=1)
mc = load("WtoLNu2Jets", era=1)
mc_w = mc["Lumi_reweight"]
N_D, W_MC = len(dat["Ds_mass"]), mc_w.sum()
print(f"  signal {len(sig['score_dsmu']):,} | OS data {N_D:,} | "
      f"W+jets {len(mc_w):,} rows, {W_MC:,.0f} weighted")


def curve(v_sig, v_dat, v_mc, higher_is_signal):
    """eps_data/eps_MC at each target eps_MC, plus eps_sig there.

    The threshold is set on the MC (so eps_MC hits the target exactly), then applied
    to data and signal. Errors: Poisson on the data numerator, sum-of-weights-squared
    on the MC denominator; the MC's own eps enters both the threshold definition and
    the denominator, which is fine -- it is one sample used once.
    """
    s = 1.0 if higher_is_signal else -1.0        # flip so "larger is more signal-like"
    vs, vd, vm = s * v_sig, s * v_dat, s * v_mc
    out = []
    for target in EPS_MC_GRID:
        # Weighted quantile of the MC: the threshold keeping exactly `target` of the
        # weighted background. Unweighted np.quantile would key on the raw jet-bin
        # mixture instead of the physical one.
        order = np.argsort(vm)[::-1]
        cw = np.cumsum(mc_w[order])
        idx = np.searchsorted(cw, target * W_MC)
        if idx >= len(order):
            out.append((target, np.nan, np.nan, np.nan))
            continue
        thr = vm[order[idx]]
        km = vm >= thr
        nm, nm2 = mc_w[km].sum(), (mc_w[km] ** 2).sum()
        nd = int((vd >= thr).sum())
        eps_sig = float((vs >= thr).mean())
        # ACHIEVED eps_MC, not the target. Features with a large spike at one value
        # (sumChargedHadronPt is exactly 0 for a big isolated fraction) saturate: the
        # threshold cannot move past the spike, so every tighter target returns the
        # same selection. That shows up as an identical ratio repeated down the row,
        # which would read as "flat in the tail" when it is really "not in the tail
        # at all". Reporting what was achieved is the only way to see it.
        eps_mc = nm / W_MC
        if nd <= 0 or nm <= 0:
            out.append((target, np.nan, np.nan, eps_sig, eps_mc))
            continue
        r = (nd / N_D) / eps_mc
        e = r * np.sqrt(1 / nd + nm2 / nm ** 2)
        out.append((target, r, e, eps_sig, eps_mc))
    return out


# Which features to scan, and which direction is signal-like. The direction is
# measured rather than assumed: whichever orientation gives the signal a higher
# median than the background is the signal-like one.
imp = np.load(f"{MODEL_DIR}/feature_importances.npy")
top = [common.FEATURE_NAMES[i] for i in np.argsort(imp)[::-1][:N_FEATURES]]
print(f"\nTop-{N_FEATURES} features by permutation importance: {', '.join(top)}")

print("\n=== eps_data / eps_MC at matched BACKGROUND efficiency, OS control region ===")
print("    1.00 = the cut removes data and W+jets MC in the same proportion.")
print("    (pre-cut data/MC is 1.44, so this ratio times 1.44 is the post-cut data/MC)\n")
hdr = f"{'variable':>32s} " + " ".join(f"{f'eps_MC={e:.0e}':>16s}" for e in EPS_MC_GRID)
print(hdr)

curves = {}
for name in ["score_dsmu"] + top:
    hi = np.median(sig[name]) > np.median(mc[name])
    c = curve(sig[name], dat[name], mc[name], hi)
    curves[name] = c
    # "!" marks a point whose achieved eps_MC missed the target by more than 50% --
    # the cut saturated and this is not the requested rejection.
    cells = []
    for tgt, r, e, _, em in c:
        if np.isnan(r):
            cells.append("n/a".rjust(16))
            continue
        flag = "!" if em > 1.5 * tgt else " "
        cells.append(f"{r:.2f}+-{e:.2f}{flag}".rjust(16))
    print(f"{name:>32s} " + " ".join(cells))

print(f"\n{'achieved eps_MC:':>32s}")
for name, c in curves.items():
    print(f"{name:>32s} " + " ".join(f"{em:16.1e}" for _, _, _, _, em in c))
print(f"\n{'eps_sig at the DNN points:':>32s} " + " ".join(
    f"{100*s:15.1f}%" for _, _, _, s, _ in curves['score_dsmu']))

# ---------- plot ----------
fig, ax = plt.subplots(figsize=(12, 8.5))
for i, (name, c) in enumerate(curves.items()):
    # Plotted against the ACHIEVED eps_MC, so a saturated feature stacks its points
    # on top of each other instead of pretending to reach the tail.
    xs = [em for _, r, _, _, em in c if np.isfinite(r)]
    ys = [r for _, r, _, _, _ in c if np.isfinite(r)]
    es = [e for _, r, e, _, _ in c if np.isfinite(r)]
    if name == "score_dsmu":
        ax.errorbar(xs, ys, yerr=es, marker="o", markersize=11, lw=3.5, capsize=4,
                    color="black", zorder=5, label="DNN score (all 43 features)")
    else:
        # Flag saturation in the legend, not only in the log: a feature whose points
        # pile up at one x never reached the tail, and this figure will outlive the
        # printout that says so.
        reached = min(xs) if xs else 1.0
        tag = "" if reached <= 1.5 * EPS_MC_GRID[-1] else \
            rf"  (saturates at $\varepsilon_{{MC}}\approx${reached:.0e})"
        ax.errorbar(xs, ys, yerr=es, marker="s", markersize=6, lw=1.8, capsize=2,
                    alpha=0.85, label=name + tag)
ax.axhline(1.0, color="grey", lw=1.2, ls=":")
ax.set_xscale("log")
ax.invert_xaxis()
ax.set_xlabel(r"background efficiency $\varepsilon_{MC}$  (tighter $\rightarrow$)")
ax.set_ylabel(r"$\varepsilon_{data}\,/\,\varepsilon_{MC}$   (OS control region)")
ax.set_title("Does one feature explain the data/MC drift?", fontsize=16)
ax.legend(fontsize=11, ncol=2)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path("closure", "feature_tail_datamc.png", key=MODEL)
plt.savefig(out, dpi=130)
plt.close(fig)
print(f"\nSaved {out}")

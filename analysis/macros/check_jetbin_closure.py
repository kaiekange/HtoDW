"""Where does the pre-cut W+jets deficit actually live, and what does the DNN do to it?

    python3 check_jetbin_closure.py

Starting point. The lumi-weighted W+jets prediction undershoots data by a factor
1.44 before any MVA cut (with the analysis selection: 1 tight muon + tight b-veto;
without the b-veto it is 2.19). That number has always been quoted as "flat" because
it is flat in M(KKpi) -- which was read as "the shape is right, whole processes are
missing". This script asks whether it is flat in anything else.

It is not. Sliced on reconstructed jet multiplicity the deficit runs
1.11 / 1.14 / 1.79 / 2.76 / 4.76 for 0/1/2/3/>=4 jets. So the 1.44 is a mixing
average over a strongly multiplicity-dependent deficit, and any comparison whose
jet-multiplicity composition changes -- which is exactly what an MVA cut does -- will
move the ratio for reasons that have nothing to do with the classifier.

Two things this separates:

* **Is the "step" from 1.44 to ~1.0 at the first working point a composition
  artifact?** If the ratio inside every slice were unchanged and only the mixture
  moved, yes. Reported here by re-slicing at each working point.

* **Is the high-multiplicity deficit a missing process or the W+jets sample's own
  modelling?** WtoLNu2Jets is 0/1/2-jet merged, so >=3 jets come only from the
  parton shower -- a known place to undershoot. The discriminator is the b-tagged
  region (`n_b_jets_tight_dR04 >= 1`), which is the complement of the analysis
  region and is tt-enriched: tt has two real b jets and four jets at the matrix
  element. If the deficit is tt, it must be much worse there. If it is
  multiplicity modelling, it should look the same in both regions.

Data has no jet binning, so the slicing is on the RECONSTRUCTED multiplicity, not
the generator bins. Migration means every generator bin feeds every reconstructed
slice; the per-slice generator mix is printed alongside so it is visible which bin
a given slice is actually testing.

Control region only (|M-M_Ds| >= 40 MeV), matching check_data_vs_wjets.py so the
numbers are directly comparable to it.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_jetbin_closure.py
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

MODEL = "ssdata"
PRED = f"{common.resolve(MODEL)['out_dir']}/Prediction"
M_DS, SR_HALFWIDTH = 1.96835, 0.04
EFFS = [0.90, 0.70, 0.50]
BR = ["score_dsmu", "Ds_mass", "num_jets_pt30",
      "num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "event_id"]

# Reconstructed jet-multiplicity slices, exhaustive and mutually exclusive.
# Split past 2 on purpose: WtoLNu2Jets has matrix elements only up to 2 jets, so the
# 3 and >=4 slices are the ones that can only be filled by the parton shower.
SLICES = [("0 jet", lambda n: n == 0),
          ("1 jet", lambda n: n == 1),
          ("2 jet", lambda n: n == 2),
          ("3 jet", lambda n: n == 3),
          (">=4 jet", lambda n: n >= 4)]

# The analysis region and its complement. Same everything else, so the difference
# between the two columns is attributable to the b-veto alone.
REGIONS = [("b-veto", lambda nb: nb == 0),
           ("b-tagged", lambda nb: nb >= 1)]


def load(name, *, ss=False, sig=False, parity=None):
    """Load one scored sample with the muon selection only.

    The b-jet cut is deliberately NOT applied here -- both sides of it are wanted,
    so the branch is carried through and masked downstream instead.
    """
    br = BR + (["Ds_match", "mu_match"] if sig else [])
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= a["score_dsmu"] > -1
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if parity is not None:
        m &= (a["event_id"] % 2) == parity
    return {k: v[m] for k, v in a.items()}


def ratio_err(nd, nm, nm2):
    """data/MC with Poisson data error and weighted-MC error. nm2 = sum of w^2."""
    if nd <= 0 or nm <= 0:
        return np.nan, np.nan
    r = nd / nm
    return r, r * np.sqrt(1.0 / nd + nm2 / nm ** 2)


print(f"Loading scored samples from {PRED} ...")
# Signal quantiles define the working points. Signal keeps the analysis b-veto --
# the working point must be the one the analysis would actually use.
sig_parts = [load(n, sig=True) for n in common.SIGNALS]
sig = np.concatenate([p["score_dsmu"][p[common.B_VETO] == 0] for p in sig_parts])
THR = {e: float(np.quantile(sig, 1 - e)) for e in EFFS}

os_d = load("Data_2022")
ss_d_all = load("Data_2022_SS", ss=True)
ss_d = load("Data_2022_SS", ss=True, parity=1)
# Measured, not assumed to be exactly 0.5: the MC must be scaled by the same
# fraction of SS data actually used, or the ratio picks up a spurious factor ~2.
F_ODD = len(ss_d["Ds_mass"]) / len(ss_d_all["Ds_mass"])


def wjets(ss):
    """Per-jet-bin arrays, kept SEPARATE -- the per-slice generator mix needs them."""
    return {n: (load(f"{n}_SS" if ss else n, ss=ss), w * (F_ODD if ss else 1.0))
            for n, w in common.WJETS_WEIGHTS.items()}


os_m, ss_m = wjets(False), wjets(True)
print(f"  OS data {len(os_d['Ds_mass']):>10,} | SS data (odd) {len(ss_d['Ds_mass']):>9,} "
      f"| f_odd = {F_ODD:.4f}")


def base(d, thr, rfn):
    return ((d["score_dsmu"] > thr)
            & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
            & rfn(d[common.B_VETO]))


def table(tag, d, mc, rname, rfn, effs):
    # The ratio alone is misleading where the MC is tiny: 171 in the b-tagged >=4 jet
    # slice and 4.8 in the b-veto one are not comparable statements about how many
    # events are missing. data - MC, and its distribution over the slices, is what
    # says where the unmodelled events actually sit -- and comparing that SHAPE
    # between the b-veto and b-tagged regions is what tests whether the two deficits
    # share an origin (tt leaking through the veto) or are different things.
    print(f"--- {tag}, {rname} region ---")
    print(f"{'eff_sig':>8s} {'slice':>8s} {'data':>10s} {'W+jets':>11s} {'ratio':>15s} "
          f"{'excess':>11s} {'of tot':>7s} {'N_eff':>8s}   gen mix (0J/1J/2J)")
    out = {}
    for eff in effs:
        thr = -1.0 if eff >= 1.0 else THR[eff]
        sd = base(d, thr, rfn)
        smc = {bn: base(bd, thr, rfn) for bn, (bd, _) in mc.items()}
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        row = []
        for sname, sfn in SLICES:
            nd = int((sd & sfn(d["num_jets_pt30"])).sum())
            per_bin, nm, nm2 = [], 0.0, 0.0
            for bn, (bd, bw) in mc.items():
                k = int((smc[bn] & sfn(bd["num_jets_pt30"])).sum())
                per_bin.append(k * bw)
                nm += k * bw
                nm2 += k * bw ** 2
            r, e = ratio_err(nd, nm, nm2)
            ne = nm ** 2 / nm2 if nm2 > 0 else 0.0
            mix = "/".join(f"{100*p/nm:.0f}%" for p in per_bin) if nm > 0 else "-"
            row.append((sname, r, e, nd, nm, ne, nd - nm, mix, nm2))
        tot_exc = sum(x[6] for x in row)
        for sname, r, e, nd, nm, ne, exc, mix, _ in row:
            frac = f"{100*exc/tot_exc:5.1f}%" if tot_exc > 0 else "    -"
            print(f"{lbl:>8s} {sname:>8s} {nd:>10,} {nm:>11,.0f} "
                  f"{f'{r:.3f} +- {e:.3f}':>15s} {exc:>11,.0f} {frac:>7s} "
                  f"{ne:>8,.0f}   {mix}")
        print(f"{'':>8s} {'TOTAL':>8s} {'':>10s} {'':>11s} {'':>15s} {tot_exc:>11,.0f}")
        out[eff] = row
        print()
    return out


print("\n=== data / lumi-weighted W+jets, sliced on reconstructed num_jets_pt30 ===")
print("    Control region only. 1 tight muon throughout; the b cut is the region label.\n")

results = {}
for tag, d, mc in (("OS", os_d, os_m), ("SS", ss_d, ss_m)):
    for rname, rfn in REGIONS:
        # The b-tagged region is a cross-check on the ORIGIN of the deficit, not an
        # analysis region, so it is only shown uncut -- working points there would
        # just be low-statistics noise about a region the analysis never uses.
        effs = [1.0] + EFFS if rname == "b-veto" else [1.0]
        results[(tag, rname)] = table(tag, d, mc, rname, rfn, effs)

# ---------- is the spread between slices significant? ----------
# "The slices agree" has to be a number, not an impression from ratios that look
# close. chi2 of the slices against their own weighted mean.
print("=== consistency of the slices (chi2 vs their weighted mean) ===")
print(f"{'region':>10s} {'charge':>7s} {'eff_sig':>8s} {'mean':>8s} {'chi2/ndf':>10s} {'max pull':>9s}")
for (tag, rname), res in results.items():
    for eff, row in res.items():
        rs = np.array([r for _, r, e, *_ in row if np.isfinite(r)])
        es = np.array([e for _, r, e, *_ in row if np.isfinite(r)])
        if len(rs) < 2:
            continue
        mean = (rs / es ** 2).sum() / (1 / es ** 2).sum()
        chi2 = (((rs - mean) / es) ** 2).sum()
        pull = np.abs((rs - mean) / es).max()
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        print(f"{rname:>10s} {tag:>7s} {lbl:>8s} {mean:8.3f} "
              f"{chi2:7.1f}/{len(rs)-1:<3d} {pull:9.2f}")

# ---------- what does the cut do to the UNMODELLED component? ----------
# The ratio going to 1.0 can mean two very different things, and the ratio alone
# cannot tell them apart:
#   (i)  the classifier throws away almost everything, so both sides shrink to noise
#   (ii) the classifier specifically rejects the events W+jets MC does not contain
# Treat (data - MC) as the unmodelled component and measure ITS efficiency, next to
# the efficiency of the modelled part. If eps(unmodelled) << eps(W+jets), it is (ii):
# the missing background is easier for this classifier to reject than the background
# it can see -- which is what training against SS data (which contains the missing
# processes) was supposed to buy.
print("\n=== efficiency of the modelled vs the UNMODELLED background ===")
print("    unmodelled := data - lumi-weighted W+jets, per working point, OS b-veto CR\n")
print(f"{'region':>7s} {'eff_sig':>8s} {'eps(data)':>11s} {'eps(W+jets)':>12s} "
      f"{'eps(unmodelled)':>20s} {'ratio to W+jets':>16s}")
for tag in ("OS", "SS"):
    res = results[(tag, "b-veto")]
    base_row = res[1.0]
    nd0 = sum(x[3] for x in base_row)
    nm0 = sum(x[4] for x in base_row)
    exc0 = nd0 - nm0
    for eff, row in res.items():
        if eff >= 1.0:
            continue
        nd = sum(x[3] for x in row)
        nm = sum(x[4] for x in row)
        nm2 = sum(x[8] for x in row)
        e_d, e_m = nd / nd0, nm / nm0
        exc = nd - nm
        # Numerator error dominates; the uncut denominator is huge and its error
        # is negligible next to it, so it is not propagated.
        exc_e = np.sqrt(nd + nm2)
        e_x, e_x_err = exc / exc0, exc_e / exc0
        print(f"{tag:>7s} {100*eff:7.0f}% {100*e_d:10.3f}% {100*e_m:11.3f}% "
              f"{f'{100*e_x:.3f} +- {100*e_x_err:.3f}%':>20s} "
              f"{f'{e_x/e_m:+.2f}':>16s}")

# ---------- plot ----------
fig, axes = plt.subplots(1, 3, figsize=(27, 8), sharey=True)
X = np.arange(len(SLICES))
COL = {1.0: "black", 0.90: "steelblue", 0.70: "darkorange", 0.50: "firebrick"}
PANELS = [("OS", "b-veto"), ("SS", "b-veto"), ("OS", "b-tagged")]
for ax, key in zip(axes, PANELS):
    for i, (eff, row) in enumerate(results[key].items()):
        off = (i - 1.5) * 0.12
        ax.errorbar(X + off, [r for _, r, *_ in row], yerr=[e for _, _, e, *_ in row],
                    marker="o", markersize=9, lw=0, elinewidth=2, capsize=4,
                    color=COL.get(eff, "grey"),
                    label="no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%")
    ax.axhline(1.0, color="black", lw=1, ls=":")
    ax.set_xticks(X)
    ax.set_xticklabels([s for s, _ in SLICES], rotation=20)
    ax.set_xlabel(r"reconstructed $N_{jet}$ ($p_T>30$ GeV)")
    ax.set_title(f"{key[0]} -- {key[1]} region", fontsize=16)
    ax.legend(fontsize=12)
axes[0].set_ylabel("data / lumi-weighted W+jets")
axes[0].set_yscale("log")
hep.cms.label(ax=axes[0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = f"{OUTDIR}/jetbin_closure.png"
plt.savefig(out, dpi=130)
plt.close(fig)
print(f"\nSaved {out}")

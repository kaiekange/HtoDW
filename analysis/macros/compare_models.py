"""Head-to-head comparison of the two phipimu DNN background choices.

    python3 compare_models.py wjets ssdata

Both models are trained by the same DNN.py on the same 43 features (common.py) with
the same seed, so the only difference between them is the training background. That
was NOT true of the first comparison on 2026-08-10: the W+jets model was 40 features
and the SS model 43, so "SS background admits 36-54% less background" measured the
feature set as much as the background choice. This script recomputes that number on
a valid basis, with uncertainties.

Outputs three figures plus printed tables:
  1. score distributions per model, over all four samples
  2. M(KKpi) at matched signal efficiency -- does the cut sculpt the fit variable
  3. SR/sideband ratio vs signal efficiency -- the real-charm enrichment systematic

Three things this script is careful about, each of which silently breaks the
comparison if ignored:

* **Never threshold on raw score.** The two models' score scales are unrelated, so
  a shared "score > 0.99" keeps different amounts of signal in each. Every working
  point here is a per-model quantile of the signal score distribution. Two caveats
  on the resulting efficiency: the quantile is taken over HtoDsW+HtoDsstW combined
  (matching what training used), and part of that signal was in each model's own
  training set. Both effects are symmetric across the two models and use the same
  events, so the comparison stays fair -- it is the absolute efficiency that is
  slightly optimistic, not the relative verdict.

* **SS data must use the held-out half.** The ssdata model trained on even event_id;
  scoring it on all of Data_2022_SS is train-set self-validation. Only odd event_id
  rows are used here.

* **W+jets is not a neutral yardstick.** The wjets model saw all of it, so rejection
  quoted on W+jets flatters that model by construction. The primary metric is OS
  sideband data, which neither model ever saw. Every legend marks how each sample
  relates to the model being plotted.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_models.py wjets ssdata
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
M_D = 1.8696
PEAK_HW = 0.02              # D+- counting window
PEAK_SB = (0.03, 0.07)      # local sideband strips either side of it
EFFS = [0.90, 0.70, 0.50]
COLORS = ["black", "steelblue", "darkorange", "firebrick"]

# run_id/lumi_id are carried alongside event_id because event_id is only unique
# WITHIN a run: Data_2022 has 1,083,571 rows but only 1,082,925 distinct event_id,
# i.e. 646 cross-run collisions. Keying identity on event_id alone would quietly
# mis-assign those when computing the two models' overlap below.
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      common.B_VETO, "run_id", "lumi_id", "event_id"]
BR_SIG = BR + ["Ds_match", "mu_match"]

if len(sys.argv) != 3:
    sys.exit(f"usage: python3 {os.path.basename(__file__)} "
             f"<{'|'.join(common.TRAININGS)}> <{'|'.join(common.TRAININGS)}>")
MODELS = sys.argv[1:3]
for k in MODELS:
    common.resolve(k)


def load(model_key, name, *, ss=False, sig=False, sideband=False, parity=None):
    """Load a scored sample from one model's Prediction/ dir, with the analysis cuts."""
    pred = f"{common.resolve(model_key)['out_dir']}/Prediction/{name}.root"
    a = uproot.open(pred)["Events"].arrays(BR_SIG if sig else BR, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if sideband:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALFWIDTH
    if parity is not None:
        m &= (a["event_id"] % 2) == parity
    return {k: v[m] for k, v in a.items()}


def gather(model_key):
    """All four evaluation samples for one model, plus its signal score quantiles."""
    sig_parts = [load(model_key, n, sig=True) for n in common.SIGNALS]
    sig = {k: np.concatenate([p[k] for p in sig_parts]) for k in BR_SIG}
    # SS: odd parity only -- the ssdata model trained on the even half.
    ss = load(model_key, "Data_2022_SS", ss=True, parity=1)
    # OS data is SR-blinded already; the explicit sideband cut makes that intent
    # visible and keeps the yield definitions identical across samples.
    os_ = load(model_key, "Data_2022", sideband=True)
    parts = [load(model_key, n) for n in common.WJETS_WEIGHTS]
    wj = {k: np.concatenate([p[k] for p in parts]) for k in BR}
    wj_w = np.concatenate([np.full(len(p["Ds_mass"]), w)
                           for p, w in zip(parts, common.WJETS_WEIGHTS.values())])
    return dict(sig=sig, ss=ss, os=os_, wj=wj, wj_w=wj_w)


def provenance(model_key, sample):
    """How the sample relates to this model -- belongs in every legend."""
    if sample == "wj":
        return "trained on" if model_key == "wjets" else "unseen"
    if sample == "ss":
        return "held-out half" if model_key == "ssdata" else "unseen"
    if sample == "os":
        return "unseen"
    return "signal"


print("Loading scored samples...")
D = {k: gather(k) for k in MODELS}
for k in MODELS:
    print(f"  {k:8s} signal {len(D[k]['sig']['score_dsmu']):>8,} | "
          f"OS sideband {len(D[k]['os']['score_dsmu']):>8,} | "
          f"SS odd {len(D[k]['ss']['score_dsmu']):>8,} | "
          f"W+jets {len(D[k]['wj']['score_dsmu']):>8,}")

THR = {k: {e: np.quantile(D[k]["sig"]["score_dsmu"], 1 - e) for e in EFFS} for k in MODELS}


# ---------- 1. score distributions ----------
fig, axes = plt.subplots(len(MODELS), 2, figsize=(18, 8 * len(MODELS)))
bins = np.linspace(0, 1, 51)
for row, k in enumerate(MODELS):
    for col, logy in enumerate((False, True)):
        ax = axes[row][col]
        ax.hist(D[k]["wj"]["score_dsmu"], bins=bins, weights=D[k]["wj_w"], density=True,
                alpha=0.5, color="seagreen",
                label=f"W+jets MC [{provenance(k,'wj')}]")
        ax.hist(D[k]["ss"]["score_dsmu"], bins=bins, density=True, histtype="step",
                lw=2.5, color="steelblue", label=f"SS data [{provenance(k,'ss')}]")
        ax.hist(D[k]["os"]["score_dsmu"], bins=bins, density=True, histtype="step",
                lw=2.5, ls="--", color="black", label="OS data sideband [unseen]")
        ax.hist(D[k]["sig"]["score_dsmu"], bins=bins, density=True, alpha=0.4,
                color="firebrick", label="signal")
        ax.set_xlabel(f"DNN score -- trained vs {common.resolve(k)['label']}")
        ax.set_ylabel("a.u.")
        ax.legend(fontsize=12)
        if logy:
            ax.set_yscale("log")
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=11)
plt.tight_layout()
out = f"{OUTDIR}/compare_models_score.png"
plt.savefig(out, dpi=120)
plt.close(fig)
print(f"\nSaved {out}")


# ---------- 2. M(KKpi) at matched signal efficiency ----------
# SS and W+jets are NOT blinded, so their Ds window is shown: whether a real Ds peak
# grows there under the cut is the whole question. Every sample is normalized to its
# own SIDEBAND integral, so the combinatorial continuum lines up and a real-Ds excess
# shows as a bump above it. OS data is simply empty in that window.
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)   # 20 MeV; blind window on bin edges
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SB_BIN = np.abs(mcent - M_DS) >= SR_HALFWIDTH
SR_BIN = ~SB_BIN


def sb_normed(mass, weight, cut):
    h, _ = np.histogram(mass[cut], bins=mbins,
                        weights=None if weight is None else weight[cut])
    n = h[SB_BIN].sum()
    return (h / n if n > 0 else h), h


cols = [1.0] + EFFS
fig, axes = plt.subplots(len(MODELS), len(cols), figsize=(9 * len(cols), 8 * len(MODELS)))
for row, k in enumerate(MODELS):
    for col, eff in enumerate(cols):
        ax = axes[row][col]
        thr = -1.0 if eff >= 1.0 else THR[k][eff]
        hw, rw = sb_normed(D[k]["wj"]["Ds_mass"], D[k]["wj_w"], D[k]["wj"]["score_dsmu"] > thr)
        hs, rs = sb_normed(D[k]["ss"]["Ds_mass"], None, D[k]["ss"]["score_dsmu"] > thr)
        ho, ro = sb_normed(D[k]["os"]["Ds_mass"], None, D[k]["os"]["score_dsmu"] > thr)
        hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"W+jets [{provenance(k,'wj')}] N={rw.sum():,.0f}")
        hep.histplot(hs, mbins, ax=ax, histtype="step", lw=2.5, color="steelblue",
                     label=f"SS [{provenance(k,'ss')}] N={rs.sum():,.0f}")
        hep.histplot(ho, mbins, ax=ax, histtype="step", lw=2.5, ls="--", color="black",
                     label=f"OS sideband N={ro.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        ax.set_ylabel("Normalized to sideband integral")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{common.resolve(k)['label']} -- {lbl}", fontsize=15)
        ax.legend(fontsize=10, loc="upper right")
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=10)
plt.tight_layout()
out = f"{OUTDIR}/compare_models_mass.png"
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")


# ---------- 3. yields on OS sideband data, with uncertainties ----------
def split_yield(mass):
    """(total, D+- excess, combinatorial) with Poisson errors, via local sideband sub.

    Disjoint Poisson pieces: n_pk (peak window), n_sb (the two flanking strips),
    n_other (everything else). The strips are twice the width of the peak window, so
    the expected background under the peak is n_sb/2.
      excess = n_pk - n_sb/2
      comb   = total - excess = n_other + 1.5 * n_sb
    Writing both in terms of disjoint counts is what makes the errors honest: the
    naive sqrt(total) would double-count the peak events.
    """
    d = np.abs(mass - M_D)
    n_pk = int((d < PEAK_HW).sum())
    n_sb = int(((d > PEAK_SB[0]) & (d < PEAK_SB[1])).sum())
    n_tot = len(mass)
    n_other = n_tot - n_pk - n_sb
    exc, exc_e = n_pk - 0.5 * n_sb, np.sqrt(n_pk + 0.25 * n_sb)
    comb, comb_e = n_other + 1.5 * n_sb, np.sqrt(n_other + 2.25 * n_sb)
    return (n_tot, np.sqrt(n_tot)), (exc, exc_e), (comb, comb_e)


def event_keys(d, sel):
    """(run, lumi, event) identity for the selected rows -- event_id alone collides."""
    return set(zip(d["run_id"][sel].tolist(),
                   d["lumi_id"][sel].tolist(),
                   d["event_id"][sel].tolist()))


def ratio_overlapping(sa, sb):
    """Ratio N_b/N_a for two selections of the SAME events, with the correct error.

    The two models cut on the same OS sideband data, so their surviving sets overlap
    heavily and sqrt(1/Na + 1/Nb) would badly overstate the uncertainty. Decompose
    into disjoint Poisson pieces -- A only, B only, common C -- and propagate.
    """
    C = len(sa & sb)
    A, B = len(sa) - C, len(sb) - C
    na, nb = A + C, B + C
    if na == 0:
        return np.nan, np.nan, (A, B, C)
    r = nb / na
    dA, dB, dC = -nb / na**2, 1.0 / na, (A - B) / na**2
    err = np.sqrt(dA**2 * A + dB**2 * B + dC**2 * C)
    return r, err, (A, B, C)


print("\n=== OS sideband data (neither model saw it) -- yields at matched eff_sig ===")
hdr = f"{'eff_sig':>8s} {'model':>8s} {'cut':>8s} {'N_total':>14s} {'D+- excess':>14s} {'combinatorial':>16s}"
print(hdr)
for eff in EFFS:
    kept = {}
    for k in MODELS:
        d = D[k]["os"]
        sel = d["score_dsmu"] > THR[k][eff]
        kept[k] = event_keys(d, sel)
        (nt, nte), (ex, exe), (cb, cbe) = split_yield(d["Ds_mass"][sel])
        print(f"{100*eff:7.0f}% {k:>8s} {THR[k][eff]:8.4f} "
              f"{nt:8,.0f} +-{nte:4.0f} {ex:8,.0f} +-{exe:4.0f} {cb:10,.0f} +-{cbe:4.0f}")
    if len(MODELS) == 2:
        a, b = MODELS
        r, re, (A, B, C) = ratio_overlapping(kept[a], kept[b])
        print(f"{'':>8s} {'ratio':>8s} {b}/{a} = {r:.3f} +- {re:.3f}   "
              f"(overlap {C:,}; {a}-only {A:,}; {b}-only {B:,})")
    print()

print("=== SR/sideband ratio (real-Ds enrichment; OS is blind so it is N/A) ===")
print(f"{'eff_sig':>8s} " + " ".join(f"{k+'/'+s:>16s}" for k in MODELS for s in ("SS", "W+jets")))
for eff in [1.0] + EFFS:
    row = []
    for k in MODELS:
        for s, wt in (("ss", None), ("wj", D[k]["wj_w"])):
            d = D[k][s]
            thr = -1.0 if eff >= 1.0 else THR[k][eff]
            cut = d["score_dsmu"] > thr
            _, h = sb_normed(d["Ds_mass"], wt, cut)
            row.append(h[SR_BIN].sum() / h[SB_BIN].sum() if h[SB_BIN].sum() > 0 else np.nan)
    lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
    print(f"{lbl:>8s} " + " ".join(f"{v:16.4f}" for v in row))

fig, ax = plt.subplots(figsize=(10, 8))
xs = [1.0] + EFFS
for k, ls in zip(MODELS, ["-", "--"]):
    for s, c, wt in (("ss", "steelblue", None), ("wj", "seagreen", D[k]["wj_w"])):
        ys = []
        for eff in xs:
            d = D[k][s]
            thr = -1.0 if eff >= 1.0 else THR[k][eff]
            _, h = sb_normed(d["Ds_mass"], wt, d["score_dsmu"] > thr)
            ys.append(h[SR_BIN].sum() / h[SB_BIN].sum())
        ax.plot([100 * x for x in xs], ys, ls, color=c, marker="o", lw=2,
                label=f"{'SS' if s=='ss' else 'W+jets'} -- {common.resolve(k)['label']} model")
ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
ax.set_ylabel("SR / sideband yield ratio")
ax.invert_xaxis()
ax.legend(fontsize=12)
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = f"{OUTDIR}/compare_models_srsb.png"
plt.savefig(out, dpi=130)
print(f"\nSaved {out}")

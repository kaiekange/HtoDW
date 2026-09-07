"""Head-to-head comparison of two or more trained phipimu classifiers.

    python3 compare_models.py wjets ssdata                      # the two backgrounds
    python3 compare_models.py DNN:ssdata BDT:ssdata             # the two algorithms
    python3 compare_models.py DNN:ssdata:seed101 DNN:ssdata     # the seed noise floor

A model spec is `[TOOL:]background[:seedN]`, resolved through common.out_dir_for.
Bare `ssdata` means `DNN:ssdata` (seed 42), so the original two-argument form still
works. Any number of specs >= 2 may be given; every table and figure gains a row.

Both DNN models are trained by the same DNN.py on the same 43 features (common.py)
with the same seed, so the only difference between them is the training background.
That was NOT true of the first comparison on 2026-08-10: the W+jets model was 40
features and the SS model 43, so "SS background admits 36-54% less background"
measured the feature set as much as the background choice. This script recomputes
that number on a valid basis, with uncertainties.

**Nothing here is meaningful below the seed noise floor.** Two runs of the same
config with different seeds do not select the same events; run the seed form above
and compare any claimed improvement against that spread before believing it.

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
  scoring it on all of Data_SS is train-set self-validation. Only odd event_id
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

# Every figure here overlays N model specs, so none of them belongs to one model:
# they go to the shared comparisons/ directory. Writing into the first spec's
# directory would make the destination depend on argument order.
M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D = 1.8696
PEAK_HW = 0.02              # D+- counting window
PEAK_SB = (0.03, 0.07)      # local sideband strips either side of it
EFFS = [0.90, 0.70, 0.50]
COLORS = ["black", "steelblue", "darkorange", "firebrick"]

# run_id/lumi_id are carried alongside event_id because event_id is only unique
# WITHIN a run: the preEE data slice has ~1.08M rows but fewer distinct event_id,
# i.e. 646 cross-run collisions. Keying identity on event_id alone would quietly
# mis-assign those when computing the two models' overlap below.
BR = ["score_dsmu", "Ds_mass", "num_os_mu_tight", "num_ss_mu_tight",
      common.B_VETO, "run_id", "lumi_id", "event_id"]
BR_SIG = BR + ["Ds_match", "mu_match"]

if len(sys.argv) < 3:
    sys.exit(f"usage: python3 {os.path.basename(__file__)} <spec> <spec> [spec ...]\n"
             f"  spec = [{'|'.join(common.TOOLS)}:]{'|'.join(common.TRAININGS)}[:seedN]\n"
             f"  e.g. ssdata wjets | DNN:ssdata BDT:ssdata | DNN:ssdata:seed101 DNN:ssdata")
MODELS = sys.argv[1:]


def parse_spec(spec):
    """[TOOL:]background[:seedN] -> (tool, bkg, seed). Bare 'ssdata' = DNN, seed 42."""
    parts = spec.split(":")
    tool = "DNN"
    seed = common.SEED_DEFAULT
    if parts and parts[0] in common.TOOLS:
        tool = parts.pop(0)
    if parts and parts[-1].startswith("seed"):
        seed = int(parts.pop()[4:])
    if len(parts) != 1:
        sys.exit(f"cannot parse model spec '{spec}'")
    common.resolve(parts[0])       # fail loudly on a typo'd background
    return tool, parts[0], seed


SPEC = {s: parse_spec(s) for s in MODELS}
LABEL = {s: f"{SPEC[s][0]} / {common.resolve(SPEC[s][1])['label']}"
            + (f" / seed{SPEC[s][2]}" if SPEC[s][2] != common.SEED_DEFAULT else "")
         for s in MODELS}


def load(model_key, name, *, ss=False, sig=False, sideband=False, parity=None, era=None):
    """Load a scored sample from one model's Prediction/ dir, with the analysis cuts.

    Returns None if the file is absent. That is not an error here: the seed-spread
    runs deliberately score only signal + Data (each full Prediction/ is ~2.8 GB
    on EOS), so the SS and W+jets sections below simply have nothing to draw.
    """
    tool, bkg, seed = SPEC[model_key]
    pred = f"{common.out_dir_for(bkg, tool=tool, seed=seed)}/Prediction/{name}.root"
    if not os.path.exists(pred):
        return None
    br = list(dict.fromkeys((BR_SIG if sig else BR) + ["Lumi_reweight", "era"]))
    a = uproot.open(pred)["Events"].arrays(br, library="np")
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
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


def gather(model_key):
    """All four evaluation samples for one model, plus its signal score quantiles.

    Any sample whose Prediction/ file is missing comes back as None; every consumer
    below skips a None rather than assuming all four are present. Signal is the one
    exception -- without it there are no score quantiles and nothing to compare.
    """
    sig_parts = [load(model_key, n, sig=True) for n in common.SIGNALS]
    if any(p is None for p in sig_parts):
        sys.exit(f"{model_key}: signal was not scored -- nothing to set the "
                 f"efficiency quantiles from. Re-run predict.py without --samples, "
                 f"or with the signals included.")
    sig = {k: np.concatenate([p[k] for p in sig_parts]) for k in BR_SIG}
    # SS: odd parity only -- the ssdata model trained on the even half.
    ss = load(model_key, "Data_SS", ss=True, parity=1, era=1)
    # OS data is SR-blinded already; the explicit sideband cut makes that intent
    # visible and keeps the yield definitions identical across samples.
    os_ = load(model_key, "Data", sideband=True, era=1)
    # preEE only (era=1), as before; per-event Lumi_reweight from the file.
    wj = load(model_key, "WtoLNu2Jets", era=1)
    wj_w = None if wj is None else wj["Lumi_reweight"]
    return dict(sig=sig, ss=ss, os=os_, wj=wj, wj_w=wj_w)


def provenance(model_key, sample):
    """Legend tag marking a sample the model was trained on -- empty when it wasn't.

    Only the positive case is annotated: an unlabelled curve is one the model never
    saw, which is the default and does not need saying on every entry. Keyed on the
    parsed background, not the raw spec, so DNN:ssdata:seed101 is still recognised.
    """
    bkg = SPEC[model_key][1]
    if sample == "wj" and bkg == "wjets":
        return " [trained on]"
    if sample == "ss" and bkg == "ssdata":
        return " [trained on half]"
    return ""


def n_of(d, s):
    return f"{len(d[s]['score_dsmu']):>8,}" if d[s] is not None else "   absent"


print("Loading scored samples...")
D = {k: gather(k) for k in MODELS}
W = max(len(k) for k in MODELS)
for k in MODELS:
    print(f"  {k:{W}s} signal {len(D[k]['sig']['score_dsmu']):>8,} | "
          f"OS sideband {n_of(D[k],'os')} | "
          f"SS odd {n_of(D[k],'ss')} | "
          f"W+jets {n_of(D[k],'wj')}")

THR = {k: {e: np.quantile(D[k]["sig"]["score_dsmu"], 1 - e) for e in EFFS} for k in MODELS}


# ---------- 1. score distributions ----------
fig, axes = plt.subplots(len(MODELS), 2, figsize=(18, 8 * len(MODELS)), squeeze=False)
bins = np.linspace(0, 1, 51)
for row, k in enumerate(MODELS):
    for col, logy in enumerate((False, True)):
        ax = axes[row][col]
        if D[k]["wj"] is not None:
            ax.hist(D[k]["wj"]["score_dsmu"], bins=bins, weights=D[k]["wj_w"], density=True,
                    alpha=0.5, color="seagreen",
                    label=f"W+jets MC{provenance(k,'wj')}")
        if D[k]["ss"] is not None:
            ax.hist(D[k]["ss"]["score_dsmu"], bins=bins, density=True, histtype="step",
                    lw=2.5, color="steelblue", label=f"SS data{provenance(k,'ss')}")
        if D[k]["os"] is not None:
            ax.hist(D[k]["os"]["score_dsmu"], bins=bins, density=True, histtype="step",
                    lw=2.5, ls="--", color="black", label="OS data sideband")
        ax.hist(D[k]["sig"]["score_dsmu"], bins=bins, density=True, alpha=0.4,
                color="firebrick", label="signal")
        ax.set_xlabel(f"MVA score -- {LABEL[k]}")
        ax.set_ylabel("a.u.")
        ax.legend(fontsize=12)
        if logy:
            ax.set_yscale("log")
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=11)
plt.tight_layout()
out = common.fig_path("comparisons", "compare_models_score.png")
plt.savefig(out, dpi=120)
plt.close(fig)
print(f"\nSaved {out}")


# ---------- 2. M(KKpi) at matched signal efficiency ----------
# SS and W+jets are NOT blinded, so their Ds window is shown: whether a real Ds peak
# grows there under the cut is the whole question. OS data is simply empty in that
# window.
#
# THIS IS NOT A DATA/MC CLOSURE PLOT and must not be read as one. The three curves in
# a panel are not on a common normalization: W+jets is OS and covers the full mass
# range, SS data is the odd-parity half only (~50% of it), and OS data is
# sideband-only because the Ds window is blinded. So the legend's N values are not
# comparable to each other, and W+jets "towering over" the data points at M_Ds is
# mostly just the blind window. For an actual closure test -- SS against SS, with the
# CR and SR ratios quoted separately -- use check_data_vs_wjets.py.
# (An earlier version of this comment claimed each curve was normalized to its own
# sideband integral. It never was; the code below plots raw counts.)
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 31)   # 20 MeV; blind window on bin edges
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SB_BIN = np.abs(mcent - M_DS) >= SR_HALFWIDTH
SR_BIN = ~SB_BIN


def sb_normed(mass, weight, cut):
    """(normalized hist, raw hist, normalized error, raw error). Error is sqrt(N) for
    unweighted data and sqrt(sum w^2) for the lumi-weighted MC; the sideband integral
    used for the normalization is treated as fixed, which is the usual shape-plot
    convention."""
    m = mass[cut]
    w = None if weight is None else weight[cut]
    h, _ = np.histogram(m, bins=mbins, weights=w)
    if w is None:
        var = h.astype(float)
    else:
        var, _ = np.histogram(m, bins=mbins, weights=w ** 2)
    raw_err = np.sqrt(var)
    n = h[SB_BIN].sum()
    if n <= 0:
        return h, h, np.zeros(len(h)), raw_err
    return h / n, h, raw_err / n, raw_err


# Layout: one row per MODEL, one column per working point, and all three samples
# overlaid within each panel. Every panel title names the classifier it was scored
# with, so no column is ambiguous about which model produced it.
#
# Raw (unnormalized) counts -- W+jets is lumi-weighted so its absolute yield is a
# genuine prediction, and SS/OS are actual observed counts. No per-panel axis is
# shared: the tight columns have orders of magnitude fewer events than "no cut", and
# forcing one y-scale across a row would flatten every tight-cut panel to a line at
# the bottom.
MASS_EFFS = [1.0, 0.90, 0.70, 0.50]

fig, axes = plt.subplots(len(MODELS), len(MASS_EFFS),
                         figsize=(8.5 * len(MASS_EFFS), 7.5 * len(MODELS)), squeeze=False)
for row, k in enumerate(MODELS):
    for col, eff in enumerate(MASS_EFFS):
        ax = axes[row][col]
        thr = -1.0 if eff >= 1.0 else THR[k][eff]
        # MC as a filled histogram; real data as points with Poisson error bars.
        if D[k]["wj"] is not None:
            _, rw, _, ew = sb_normed(D[k]["wj"]["Ds_mass"], D[k]["wj_w"],
                                     D[k]["wj"]["score_dsmu"] > thr)
            hep.histplot(rw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                         label=f"W+jets MC{provenance(k,'wj')}  N={rw.sum():,.1f}")
        if D[k]["ss"] is not None:
            _, rs, _, es = sb_normed(D[k]["ss"]["Ds_mass"], None, D[k]["ss"]["score_dsmu"] > thr)
            hep.histplot(rs, mbins, ax=ax, histtype="errorbar", yerr=es, color="steelblue",
                         marker="o", markersize=5, elinewidth=1.4,
                         label=f"SS data{provenance(k,'ss')}  N={rs.sum():,.0f}")
        if D[k]["os"] is not None:
            _, ro, _, eo = sb_normed(D[k]["os"]["Ds_mass"], None, D[k]["os"]["score_dsmu"] > thr)
            hep.histplot(ro, mbins, ax=ax, histtype="errorbar", yerr=eo, color="black",
                         marker="s", markersize=4.5, elinewidth=1.4,
                         label=f"OS data sideband  N={ro.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel("Events / 20 MeV")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{LABEL[k]} -- {lbl}", fontsize=15)
        ax.legend(fontsize=11, loc="upper right")
        if row == 0 and col == 0:
            hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                          fontsize=11)
plt.tight_layout()
out = common.fig_path("comparisons", "compare_models_mass.png")
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


if all(D[k]["os"] is None for k in MODELS):
    sys.exit("\nNo model has Data scored -- the OS sideband table is the primary "
             "metric and cannot be produced. Nothing further to do.")

print("\n=== OS sideband data (neither model saw it) -- yields at matched eff_sig ===")
print("N_kept below is the metric. Ratios are quoted against the FIRST spec listed;")
print("the error accounts for the two selections sharing most of their events.")
hdr = (f"{'eff_sig':>8s} {'model':>{W}s} {'cut':>8s} {'N_total':>14s} "
       f"{'D+- excess':>14s} {'combinatorial':>16s}")
print(hdr)
for eff in EFFS:
    kept = {}
    for k in MODELS:
        d = D[k]["os"]
        if d is None:
            print(f"{100*eff:7.0f}% {k:>{W}s} {'--':>8s}  (Data not scored)")
            continue
        sel = d["score_dsmu"] > THR[k][eff]
        kept[k] = event_keys(d, sel)
        (nt, nte), (ex, exe), (cb, cbe) = split_yield(d["Ds_mass"][sel])
        print(f"{100*eff:7.0f}% {k:>{W}s} {THR[k][eff]:8.4f} "
              f"{nt:8,.0f} +-{nte:4.0f} {ex:8,.0f} +-{exe:4.0f} {cb:10,.0f} +-{cbe:4.0f}")
    ref = next((k for k in MODELS if k in kept), None)
    for k in MODELS:
        if k == ref or k not in kept:
            continue
        r, re, (A, B, C) = ratio_overlapping(kept[ref], kept[k])
        print(f"{'':>8s} {'ratio':>{W}s} {k}/{ref} = {r:.3f} +- {re:.3f}   "
              f"(overlap {C:,}; {ref}-only {A:,}; {k}-only {B:,})")
    # With >2 specs the spread across them is the number to read (this is how the
    # seed noise floor is measured); with 2 it is just the ratio again.
    ns = [len(kept[k]) for k in MODELS if k in kept]
    if len(ns) > 2:
        print(f"{'':>8s} {'spread':>{W}s} N_kept mean {np.mean(ns):,.1f} "
              f"RMS {np.std(ns, ddof=1):,.1f} ({100*np.std(ns, ddof=1)/np.mean(ns):.2f}%) "
              f"min {min(ns):,} max {max(ns):,}")
    print()

# (model, sample) pairs that actually have data behind them. With seed-spread runs
# both are usually empty, so the whole section is skipped rather than half-drawn.
CURVES = [(k, s) for k in MODELS for s in ("ss", "wj") if D[k][s] is not None]

if not CURVES:
    print("=== SR/sideband ratio: skipped, no SS or W+jets sample was scored ===")
else:
    def srsb(k, s, eff):
        d = D[k][s]
        thr = -1.0 if eff >= 1.0 else THR[k][eff]
        _, h, _, _ = sb_normed(d["Ds_mass"], D[k]["wj_w"] if s == "wj" else None,
                               d["score_dsmu"] > thr)
        return h[SR_BIN].sum() / h[SB_BIN].sum() if h[SB_BIN].sum() > 0 else np.nan

    print("=== SR/sideband ratio (real-Ds enrichment; OS is blind so it is N/A) ===")
    cols = [f"{k}/{'SS' if s == 'ss' else 'W+jets'}" for k, s in CURVES]
    cw = max(16, max(len(c) for c in cols))
    print(f"{'eff_sig':>8s} " + " ".join(f"{c:>{cw}s}" for c in cols))
    for eff in [1.0] + EFFS:
        lbl = "no cut" if eff >= 1.0 else f"{100*eff:.0f}%"
        print(f"{lbl:>8s} " + " ".join(f"{srsb(k, s, eff):{cw}.4f}" for k, s in CURVES))

    fig, ax = plt.subplots(figsize=(10, 8))
    xs = [1.0] + EFFS
    # Each curve gets its own colour AND marker: encoding the model in the linestyle
    # alone made the legend entries hard to tell apart at legend scale. The palette
    # cycles, so an arbitrary number of specs stays readable.
    PALETTE = ["steelblue", "seagreen", "darkorange", "firebrick",
               "purple", "olive", "teal", "crimson"]
    MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*"]
    for i, (k, s) in enumerate(CURVES):
        ys = [srsb(k, s, eff) for eff in xs]
        ax.plot([100 * x for x in xs], ys, lw=2, markersize=9,
                color=PALETTE[i % len(PALETTE)], marker=MARKERS[i % len(MARKERS)],
                linestyle="-" if i % 2 == 0 else "--",
                label=f"{'SS' if s == 'ss' else 'W+jets'} sample -- {LABEL[k]}")
    ax.set_xlabel(r"signal efficiency $\varepsilon_{sig}$ [%]")
    ax.set_ylabel("SR / sideband yield ratio")
    ax.invert_xaxis()
    ax.legend(fontsize=12)
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
    plt.tight_layout()
    out = common.fig_path("comparisons", "compare_models_srsb.png")
    plt.savefig(out, dpi=130)
    print(f"\nSaved {out}")

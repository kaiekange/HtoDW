"""W+jets against data in M(KKpi), preEE and postEE separately, at matched signal eff.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_incl_by_era.py              # DNN:ssdata, jet-binned W+jets, OS
    python3 compare_incl_by_era.py incl         # inclusive W+jets sample instead
    python3 compare_incl_by_era.py SS           # same-sign control region
    python3 compare_incl_by_era.py phi=0.010    # tighten the phi window to +-10 MeV
    python3 compare_incl_by_era.py BDT:ssdata

One figure: rows are (preEE spectrum, preEE data/MC ratio, postEE spectrum, postEE
data/MC ratio), columns are the working points (no cut / 90 / 70 / 50 %). The ratio
sits directly under its spectrum, same x-axis.

Three things kept honest here, each of which has bitten this analysis before:

* Thresholds are computed once, from truth-matched signal, and reused for both eras.
  Signal MC is 100% preEE, so a per-era quantile would silently make the two rows
  different cuts.

* Nothing is renormalized to the data. W+jets carries sigma*L per era (the merged
  file's per-event Lumi_reweight already has the right era luminosity baked in), so
  its height is a prediction; each era's W+jets is compared to that same era's data.

* OS data is SR-blinded (|M-M_Ds| < 40 MeV is empty, verified for both eras), so the
  ratio there is drawn as "no ratio" and sideband yields are quoted on the sideband
  only. W+jets is drawn full-range with the blind window greyed, because hiding the
  MC Ds peak would hide the thing worth looking at.

This is NOT the preEE-vs-postEE era control test CLAUDE.md lists as blocking the
"the new model fixed the mass tilt" claim -- the eras differ in pileup and detector
conditions on their own, so a per-era closure difference cannot be pinned on the
model's background-77%-postEE / signal-100%-preEE asymmetry.
"""
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
RATIO_YMAX = 3.0            # fixed, never autoscaled (feedback_plot_style incident)
LOW_MAX, HIGH_MIN = 1.830, 2.008

# W+jets description: "bins" = the 0J+1J+2J jet-binned production (one merged file,
# jet_bin folded in), "incl" = the separate inclusive MiniAOD. ALTERNATIVES for the
# same process, never summed -- each has its own sigma and sumGenWeights_.
SOURCES = {"bins": ("WtoLNu2Jets",      "W+jets 0J+1J+2J", "bins"),
           "incl": ("WtoLNu2Jets_incl", "W+jets incl",     "incl")}

ERA_INFO = {1: ("2022 (preEE)",    "7.99 fb$^{-1}$"),
            2: ("2022EE (postEE)", "26.68 fb$^{-1}$")}

BR = ["score_dsmu", "Ds_mass", "phi_mass", "Lumi_reweight", "event_id", "era",
      *common.SELECTION]
BR_SIG = [b for b in BR if b not in ("Lumi_reweight", "era")] + ["Ds_match", "mu_match"]

argv = sys.argv[1:]
SOURCE, SIGN, PHI_HW = "bins", "OS", None
for a in list(argv):
    if a in ("OS", "SS"):
        SIGN = a; argv.remove(a)
    elif a.startswith("phi="):
        PHI_HW = float(a.split("=")[1]); argv.remove(a)
    elif a in SOURCES:
        SOURCE = a; argv.remove(a)
WJ_NAME, WJ_LABEL, WJ_TAG = SOURCES[SOURCE]

SPEC = argv[0] if argv else "ssdata"
parts = SPEC.split(":")
TOOL = parts.pop(0) if parts and parts[0] in common.TOOLS else "DNN"
SEED = int(parts.pop()[4:]) if parts and parts[-1].startswith("seed") else common.SEED_DEFAULT
if len(parts) != 1:
    sys.exit(f"cannot parse model spec '{SPEC}'")
BKG = parts[0]
common.resolve(BKG)
PRED = f"{common.out_dir_for(BKG, tool=TOOL, seed=SEED)}/Prediction"


def load(name, *, sig=False, is_data=False, era=None):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        BR_SIG if sig else BR, library="np")
    # OS/SS is a property of the FILE (producer's sameSignPi), so it is an ASSERTION now,
    # not a selection; the muon cuts come from common.selection_mask. The old
    # num_os/ss_mu_tight ==1/==0 pair conflated the two and additionally imposed an
    # uncalibrated soft-in-jet-muon tag. See common.check_charge_config / selection_mask.
    common.check_charge_config(a, SIGN == "SS" and not sig, f"load({name!r})")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if PHI_HW is not None:
        m &= np.abs(a["phi_mass"] - 1.019460) < PHI_HW
    if era is not None:
        m &= a["era"] == era
    if SIGN == "SS" and is_data:
        m &= (a["event_id"] % 2) == 1          # ssdata trained on the even half
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


print(f"Model: {TOOL}/{BKG} seed{SEED}   source: {SOURCE} ({WJ_LABEL})   sign: {SIGN}"
      + ("" if PHI_HW is None else f"   phi +-{PHI_HW*1000:.0f} MeV"))

sig_score = np.concatenate([load(n, sig=True)["score_dsmu"] for n in common.SIGNALS])
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig_score, 1.0 - e))) for e in EFFS}
print(f"Signal (truth-matched, both channels): {len(sig_score):,} rows")
for e in EFFS:
    if e < 1.0:
        print(f"  eps_sig={e:.0%}  ->  score > {THR[e]:.6f}")

DATA_NAME = "Data" + ("_SS" if SIGN == "SS" else "")
WJ_FILE = WJ_NAME + ("_SS" if SIGN == "SS" else "")

# f_odd: SS data is halved by parity; scale W+jets SS by the same measured fraction so
# the comparison stays like-for-like. OS is untouched (f_odd = 1).
D = {}
for era in ERA_INFO:
    data = load(DATA_NAME, is_data=True, era=era)
    wj = load(WJ_FILE, era=era)
    f_odd = 1.0
    if SIGN == "SS":
        all_ss = uproot.open(f"{PRED}/{DATA_NAME}.root")["Events"].arrays(
            [*common.SELECTION, "score_dsmu", "era"], library="np")
        common.check_charge_config(all_ss, True, f"f_odd({DATA_NAME!r})")
        base = (common.selection_mask(all_ss) & (all_ss["score_dsmu"] > -1)
                & (all_ss["era"] == era))
        f_odd = len(data["Ds_mass"]) / max(base.sum(), 1)
    D[era] = dict(data=data, wj=wj, wjw=wj["Lumi_reweight"] * f_odd, f_odd=f_odd)
    print(f"  era {era}: data {len(data['Ds_mass']):>9,} | {WJ_LABEL} "
          f"{len(wj['Ds_mass']):>9,} rows" + (f"  (f_odd={f_odd:.3f})" if SIGN == "SS" else ""))

# 10 MeV bins over M_Ds +- 0.3, blind window on a bin edge (0.04 / 0.01 = 4).
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.01)
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SB = np.abs(mcent - M_DS) >= SR_HALFWIDTH      # sideband bins (OS data lives here only)


def yields(era, eff):
    """Sideband data, W+jets, data/MC, low/high double ratio, N_eff for one cell."""
    thr = THR[eff]
    d, w, ww = D[era]["data"], D[era]["wj"], D[era]["wjw"]
    kd = (d["score_dsmu"] > thr) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    kw = (w["score_dsmu"] > thr) & (np.abs(w["Ds_mass"] - M_DS) >= SR_HALFWIDTH)
    md, mw, w2 = d["Ds_mass"][kd], w["Ds_mass"][kw], ww[kw]
    n_d, n_w, n_w2 = len(md), w2.sum(), (w2 ** 2).sum()
    n_eff = n_w ** 2 / n_w2 if n_w2 > 0 else 0.0
    dl, dh = (md < LOW_MAX).sum(), (md > HIGH_MIN).sum()
    wl, wh = w2[mw < LOW_MAX].sum(), w2[mw > HIGH_MIN].sum()
    if min(dl, dh, wl, wh) > 0:
        dr = (dl / wl) / (dh / wh)
        wl2 = (w2[mw < LOW_MAX] ** 2).sum(); wh2 = (w2[mw > HIGH_MIN] ** 2).sum()
        dr_err = dr * np.sqrt(1 / dl + 1 / dh + wl2 / wl ** 2 + wh2 / wh ** 2)
    else:
        dr = dr_err = float("nan")
    ratio = n_d / n_w if n_w > 0 else float("nan")
    r_err = ratio * np.sqrt(1 / max(n_d, 1) + n_w2 / max(n_w, 1e-9) ** 2)
    return dict(n_d=n_d, n_w=n_w, n_eff=n_eff, ratio=ratio, r_err=r_err,
                dr=dr, dr_err=dr_err)


print("\n" + "=" * 104)
print(f"{SIGN} sideband only (|M - M_Ds| >= 40 MeV)   -- W+jets is an absolute sigma*L prediction")
print("=" * 104)
hdr = (f"{'era':<17}{'eps_sig':>8}{'N(data)':>11}{'N(W+jets)':>12}{'N_eff':>9}"
       f"{'data/MC':>17}{'low/high dbl':>18}")
for era, (ename, _) in ERA_INFO.items():
    print("-" * 104)
    print(hdr if era == 1 else "")
    for e in EFFS:
        c = yields(era, e)
        flag = "  <- N_eff" if c["n_eff"] < MIN_NEFF else ""
        lbl = "no cut" if e >= 1.0 else f"{e:.0%}"
        print(f"{ename:<17}{lbl:>8}{c['n_d']:>11,}{c['n_w']:>12,.0f}{c['n_eff']:>9,.0f}"
              f"{c['ratio']:>10.3f} +-{c['r_err']:6.3f}{c['dr']:>11.3f} +-{c['dr_err']:6.3f}{flag}")

# ---------------------------------------------------------------------------
# Figure: one block per era, each block a spectrum row with the data/MC ratio row
# directly beneath. Nested gridspec so the gap between the two eras is wide while
# spectrum-and-its-ratio stay tight.
# ---------------------------------------------------------------------------
ncol = len(EFFS)
fig = plt.figure(figsize=(8.5 * ncol, 20))
outer = fig.add_gridspec(len(ERA_INFO), 1, hspace=0.20)

for i, (era, (ename, lumi)) in enumerate(ERA_INFO.items()):
    inner = outer[i].subgridspec(2, ncol, height_ratios=[2.4, 1.0],
                                 hspace=0.05, wspace=0.16)
    for col, eff in enumerate(EFFS):
        ax = fig.add_subplot(inner[0, col])
        rax = fig.add_subplot(inner[1, col], sharex=ax)
        ax.tick_params(labelbottom=False)
        thr = THR[eff]
        d, w, ww = D[era]["data"], D[era]["wj"], D[era]["wjw"]
        kd, kw = d["score_dsmu"] > thr, w["score_dsmu"] > thr

        hw, _ = np.histogram(w["Ds_mass"][kw], bins=mbins, weights=ww[kw])
        hw2, _ = np.histogram(w["Ds_mass"][kw], bins=mbins, weights=ww[kw] ** 2)
        hd, _ = np.histogram(d["Ds_mass"][kd], bins=mbins)

        # ---- spectrum ----
        hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"{WJ_LABEL} {SIGN} MC  N={hw.sum():,.0f}")
        # histtype="fill" silently drops yerr (feedback_mplhep_fill_yerr); the MC
        # weighted stat error (sqrt(sum w^2)) goes on as a separate errorbar call.
        ax.errorbar(mcent, hw, yerr=np.sqrt(hw2), fmt="none", ecolor="darkgreen",
                    elinewidth=1.0, capsize=0, alpha=0.7)
        hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd),
                     color="black", marker="s", markersize=4.5, elinewidth=1.4,
                     label=f"{SIGN} data{' (SR blind)' if SIGN == 'OS' else ' (odd parity)'}"
                           f"  N={hd.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        if col == 0:
            ax.set_ylabel("Events / 10 MeV")
            hep.cms.label(ax=ax, rlabel=f"{lumi} (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=11)
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{ename} -- {lbl}"
                     + ("" if PHI_HW is None else rf"   $|M_{{KK}}-M_\phi|<${PHI_HW*1000:.0f} MeV"),
                     fontsize=14, pad=26)
        ax.legend(fontsize=11, loc="upper right")

        # ---- data/MC ratio ----
        ok = (hd > 0) & (hw > 0)
        r = np.where(ok, hd / np.where(hw > 0, hw, 1), np.nan)
        rerr = np.where(ok, r * np.sqrt(1 / np.maximum(hd, 1)
                                        + hw2 / np.maximum(hw, 1e-9) ** 2), np.nan)
        inside = ok & (r <= RATIO_YMAX)
        rax.errorbar(mcent[inside], r[inside], yerr=rerr[inside], fmt="s",
                     color="black", markersize=4.0, elinewidth=1.2)
        above = ok & (r > RATIO_YMAX)
        if above.any():
            rax.plot(mcent[above], np.full(above.sum(), RATIO_YMAX * 0.97), "^",
                     color="black", markersize=7)
        rax.axhline(1.0, color="firebrick", lw=1.4)
        rax.set_ylim(0, RATIO_YMAX)
        rax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        rax.axvline(M_D, color="green", ls=":", lw=1.5)
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            rax.set_ylabel(f"{SIGN} data / MC")
        n_bad = int((~ok & SB).sum())          # sideband bins only; SR is blind for OS
        notes = [f"{n_bad} sideband bins: no ratio"] if n_bad else []
        if above.any():
            notes.append(f"{int(above.sum())} above axis")
        if notes:
            rax.text(0.98, 0.93, "    ".join(notes), transform=rax.transAxes,
                     fontsize=9, color="firebrick", va="top", ha="right")

PHITAG = "" if PHI_HW is None else f"_phi{PHI_HW*1000:.0f}"
out = common.fig_path("era", f"compare_{WJ_TAG}_{SIGN}{PHITAG}_by_era.png",
                      key=BKG, tool=TOOL, seed=SEED)
plt.savefig(out, dpi=110, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

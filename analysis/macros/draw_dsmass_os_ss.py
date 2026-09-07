"""The four basic pre-MVA M(KKpi) comparisons, absolute yields, no classifier involved.

    python3 draw_dsmass_os_ss.py            # analysis selection (1 tight mu + b-veto)
    python3 draw_dsmass_os_ss.py --no-bveto # 1 tight mu only

Four panels, each a spectrum plus its ratio:

  1. Data OS vs Data SS     -- the charge asymmetry actually present in data
  2. W+jets OS vs W+jets SS -- the charge asymmetry the MC predicts
  3. Data OS vs W+jets OS   -- closure, opposite-sign
  4. Data SS vs W+jets SS   -- closure, same-sign

Panels 1 and 2 answer "does the MC reproduce the OS/SS asymmetry"; 3 and 4 answer
"does it reproduce the yield". Together they are the baseline every post-MVA number
has to be read against, and they were never drawn.

Everything is an ABSOLUTE yield. The lumi-weighted W+jets is a sigma x L prediction
with no free normalization, so a ratio away from 1 means something on its own and
nothing is renormalized to anything.

Read from Samples/ directly, not from a model's Prediction/ -- no score is needed here
and going through a model would make a pre-MVA plot depend on which model was trained.
For the same reason SS data is used in FULL: the odd/even parity split exists only to
keep the ssdata model from being validated on its own training half, which cannot
happen when no score cut is applied.

The b-veto is a flag because it matters more than it looks: the OS ratio is 2.19
without it and 1.44 with it (it removes 25% of data and 0.7% of W+jets MC). Quoting
either number without saying which selection produced it has already caused confusion
once.

OS data is blinded in |M-M_Ds| < 40 MeV; those bins are left as gaps, never zeros.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_dsmass_os_ss.py
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

SAMPLES = common.SAMPLES_DIR

M_DS, SR_HALFWIDTH = 1.96835, 0.04
M_D = 1.8696
B_VETO = "--no-bveto" not in sys.argv
SUFFIX = "" if B_VETO else "_nobveto"

# 10 MeV bins: pre-MVA there is enough data to resolve the D+- peak, and the blind
# window edges land exactly on bin boundaries at this width.
MBINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.01)
MCENT = 0.5 * (MBINS[1:] + MBINS[:-1])
BLIND = np.abs(MCENT - M_DS) < SR_HALFWIDTH


def load(name, *, ss, era=None, want_w=False):
    br = ["Ds_mass", "num_os_mu_tight", "num_ss_mu_tight", common.B_VETO]
    if era is not None:
        br.append("era")
    if want_w:
        br.append("Lumi_reweight")
    a = uproot.open(f"{SAMPLES}/{name}.root")["Events"].arrays(br, library="np")
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    if B_VETO:
        m &= a[common.B_VETO] == 0
    if era is not None:
        m &= a["era"] == era
    if want_w:
        return a["Ds_mass"][m], a["Lumi_reweight"][m]
    return a["Ds_mass"][m]


def spectrum(mass, w=None, blind=False):
    """(counts, errors). blind=True blanks the signal window instead of zeroing it --
    an empty blinded bin is missing information, not a measurement of zero."""
    if w is None:
        h, _ = np.histogram(mass, bins=MBINS)
        h = h.astype(float)
        e = np.sqrt(h)
    else:
        h, _ = np.histogram(mass, bins=MBINS, weights=w)
        v, _ = np.histogram(mass, bins=MBINS, weights=w ** 2)
        e = np.sqrt(v)
    if blind:
        h = np.where(BLIND, np.nan, h)
        e = np.where(BLIND, np.nan, e)
    return h, e


print(f"Loading from {SAMPLES} "
      f"(selection: 1 tight muon{' + tight b-veto' if B_VETO else ', NO b-veto'}) ...")
# preEE only (era=1), which is what this macro always plotted: WJETS_WEIGHTS held only
# preEE keys and the data literal was Data_2022. Drop era=1 to include postEE.
data = {c: load("Data" if c == "OS" else "Data_SS", ss=(c == "SS"), era=1)
        for c in ("OS", "SS")}
mc, mcw = {}, {}
for c in ("OS", "SS"):
    # Per-event Lumi_reweight from the file, replacing three per-file constants.
    mc[c], mcw[c] = load("WtoLNu2Jets_SS" if c == "SS" else "WtoLNu2Jets",
                         ss=(c == "SS"), era=1, want_w=True)
for c in ("OS", "SS"):
    print(f"  {c}: data {len(data[c]):>9,} | W+jets {mcw[c].sum():>11,.0f} "
          f"({len(mc[c]):,} rows)")

H = {("data", "OS"): spectrum(data["OS"], blind=True),
     ("data", "SS"): spectrum(data["SS"]),
     ("mc", "OS"): spectrum(mc["OS"], mcw["OS"]),
     ("mc", "SS"): spectrum(mc["SS"], mcw["SS"])}

# (title, numerator, denominator, ratio label)
SEL_TXT = "1-tight-mu" + (" + b-veto" if B_VETO else "")
TAG = {("data", "OS"): f"Data OS (SR blinded, {SEL_TXT})",
       ("data", "SS"): f"Data SS ({SEL_TXT})",
       ("mc", "OS"): f"W+jets MC OS (lumi-weighted, {SEL_TXT})",
       ("mc", "SS"): f"W+jets MC SS (lumi-weighted, {SEL_TXT})"}

# (filename stem, series A, series B, ratio label). A is plotted blue, B red, in
# every figure -- keeping the colour meaning positional rather than per-sample means
# the four figures can be flipped through without re-reading a legend each time.
FIGS = [
    ("data_os_vs_ss", ("data", "OS"), ("data", "SS"), "OS / SS"),
    ("wjets_os_vs_ss", ("mc", "OS"), ("mc", "SS"), "OS / SS"),
    ("os_data_vs_wjets", ("data", "OS"), ("mc", "OS"), "data / MC"),
    ("ss_data_vs_wjets", ("data", "SS"), ("mc", "SS"), "data / MC"),
]

print("\n=== integrated ratios (control region only, |M-M_Ds| >= 40 MeV) ===")
CR = ~BLIND
for stem, ka, kb, rlab in FIGS:
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1]})
    for key, color, mk in ((ka, "steelblue", "o"), (kb, "firebrick", "s")):
        h, e = H[key]
        ax.errorbar(MCENT, h, yerr=e, fmt=mk, markersize=4.5, color=color,
                    elinewidth=1.2, label=f"{TAG[key]}   N={np.nansum(h):,.0f}")
    ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.4)
    # Linear, and anchored at zero: on a log axis a 15% difference between two
    # falling spectra is nearly invisible, which is the whole thing these figures
    # exist to show.
    ax.set_ylim(0, None)
    ax.set_ylabel("# events / 10 MeV")
    ax.legend(fontsize=12, loc="upper right")
    hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                  fontsize=11)

    hn, en = H[ka]
    hd, ed = H[kb]
    ok = np.isfinite(hn) & np.isfinite(hd) & (hn > 0) & (hd > 0)
    r = np.full(len(hn), np.nan)
    re = np.full(len(hn), np.nan)
    r[ok] = hn[ok] / hd[ok]
    re[ok] = r[ok] * np.sqrt((en[ok] / hn[ok]) ** 2 + (ed[ok] / hd[ok]) ** 2)
    rax.errorbar(MCENT, r, yerr=re, fmt="o", markersize=4.5, color="black",
                 elinewidth=1.2)
    rax.axhline(1.0, color="firebrick", lw=1.3, ls=":")
    rax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
    rax.axvline(M_D, color="green", ls=":", lw=1.4)
    rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    rax.set_ylabel(rlab)
    lo, hi = np.nanpercentile(r, [1, 99])
    pad = 0.35 * (hi - lo) if hi > lo else 0.1
    rax.set_ylim(lo - pad, hi + pad)

    # Integrated over the CR only, so both sides are counted on the same bins --
    # otherwise the OS numbers would include the SR for MC and exclude it for data.
    n_a, n_b = np.nansum(hn[CR]), np.nansum(hd[CR])
    e_a = np.sqrt(np.nansum(en[CR] ** 2))
    e_b = np.sqrt(np.nansum(ed[CR] ** 2))
    rr = n_a / n_b
    rre = rr * np.sqrt((e_a / n_a) ** 2 + (e_b / n_b) ** 2)
    print(f"  {stem:20s} {rlab:>10s} = {rr:.4f} +- {rre:.4f}   ({n_a:,.0f} / {n_b:,.0f})")

    plt.tight_layout()
    out = common.fig_path("preselection", f"dsmass_{stem}_pre_mva{SUFFIX}.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"    saved {out}")

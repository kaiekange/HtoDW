"""preEE vs postEE data in M(KKpi), scaled to a common integrated luminosity.

Same idea as compare_gentag_eras.py but for the data itself: postEE is multiplied by
L_pre / L_post = 0.29948 so both spectra represent 7.99 fb-1. If the two eras' data
behaved identically the ratio panel would sit at 1.0; a working-point-dependent
departure means the DNN cut selects the data differently in the two eras.

This pairs with compare_incl_by_era.py: there, each era's data/W+jets closure; here,
data-only across eras. W+jets postEE/preEE is ~1.0 per fb-1 (compare_gentag_eras),
so whatever moves here moves in the data, not the background model.

Tight-b-veto selection. Thresholds from truth-matched signal, both channels,
unweighted quantile, reused for both eras.

  OS (default): sideband only, since OS data is SR-blinded.
  SS          : full mass range (SS is not blinded); odd event_id half only, because
                the even half went into ssdata training. The odd fraction is ~0.5 in
                both eras and cancels in the post/pre ratio, so absolute heights are
                "odd-parity half of 7.99 fb-1" -- the ratio is the number to read.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_data_eras.py           # OS
    python3 compare_data_eras.py SS
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL)}/Prediction"

M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.040
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.01)
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SB = np.abs(mcent - M_DS) >= SR_HALF
EFFS = [1.0, 0.90, 0.70, 0.50]
RATIO_YMAX = 3.0

SIGN = "SS" if "SS" in sys.argv[1:] else "OS"
L_PRE, L_POST = common.ERAS[1]["lumi_pb"], common.ERAS[2]["lumi_pb"]
SF_POST = L_PRE / L_POST
print(f"sign: {SIGN}   postEE -> preEE luminosity scale factor: {SF_POST:.5f}")

SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "score_dsmu"]


def sel(a, sign="OS"):
    if sign == "SS":
        base = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        base = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    return base & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)


# ---------- thresholds from truth-matched signal ----------
s_score = []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s_score.append(a["score_dsmu"][m])
s_score = np.concatenate(s_score)
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(s_score, 1 - e))) for e in EFFS}
print(f"truth-matched signal: {len(s_score):,} rows")
for e in EFFS:
    if e < 1.0:
        print(f"  eps_sig={e:.0%}  ->  score > {THR[e]:.5f}")

# ---------- data, split by era ----------
DATA_FILE = "Data_SS" if SIGN == "SS" else "Data"
extra = ["event_id"] if SIGN == "SS" else []
d = uproot.open(f"{PRED}/{DATA_FILE}.root")["Events"].arrays(
    ["score_dsmu", "Ds_mass", "era"] + SEL + extra, library="np")
k = sel(d, SIGN)
if SIGN == "OS":
    k &= np.abs(d["Ds_mass"] - M_DS) >= SR_HALF            # OS is SR-blinded
    RANGE_NOTE = "sideband"
else:
    k &= (d["event_id"] % 2) == 1                          # held-out half
    RANGE_NOTE = "full range, odd parity"
sc, mass, era = d["score_dsmu"][k], d["Ds_mass"][k], d["era"][k]
is_pre, is_post = era == 1, era == 2
print(f"\n{SIGN} data rows ({RANGE_NOTE}): preEE {is_pre.sum():,}  postEE {is_post.sum():,}")

print("\n" + "=" * 76)
print(f"postEE / preEE {SIGN} data at common luminosity ({RANGE_NOTE})")
print("=" * 76)
print(f"{'eps_sig':>8}  {'N preEE':>11}  {'N postEE':>11}  {'N post scaled':>13}  {'post/pre':>16}")
for e in EFFS:
    cut = sc > THR[e]
    npre = int((cut & is_pre).sum())
    npost = int((cut & is_post).sum())
    nps = npost * SF_POST
    if npre > 0:
        r = nps / npre
        rerr = r * np.sqrt(1 / max(npre, 1) + 1 / max(npost, 1))
        rt = f"{r:.3f} +- {rerr:.3f}"
    else:
        rt = "--"
    lbl = "no cut" if e >= 1.0 else f"{e:.0%}"
    print(f"{lbl:>8}  {npre:>11,}  {npost:>11,}  {nps:>13,.0f}  {rt:>16}")

# ---------- figure: spectra + per-bin ratio ----------
ncol = len(EFFS)
fig, axes = plt.subplots(2, ncol, figsize=(8.5 * ncol, 10.0), squeeze=False,
                         sharex="col",
                         gridspec_kw=dict(height_ratios=[2.4, 1.0], hspace=0.06))
for col, e in enumerate(EFFS):
    ax, rax = axes[0][col], axes[1][col]
    ax.tick_params(labelbottom=False)
    cut = sc > THR[e]

    hpre, _ = np.histogram(mass[cut & is_pre], bins=mbins)
    hpost_raw, _ = np.histogram(mass[cut & is_post], bins=mbins)
    hpost = hpost_raw * SF_POST
    epre, epost = np.sqrt(hpre), np.sqrt(hpost_raw) * SF_POST
    npre, npost = int((cut & is_pre).sum()), int((cut & is_post).sum())

    hep.histplot(hpre, mbins, ax=ax, yerr=epre, histtype="step", color="steelblue",
                 linewidth=2.0, label=f"2022 preEE   N={npre:,}")
    hep.histplot(hpost, mbins, ax=ax, yerr=epost, histtype="step", color="firebrick",
                 linewidth=2.0, label=f"2022EE postEE (scaled)   N={npost*SF_POST:,.0f}")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    ax.axvline(M_D, color="green", ls=":", lw=1.5)
    top = ax.get_ylim()[1]
    ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
    ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
    if npre > 0:
        r = npost * SF_POST / npre
        rerr = r * np.sqrt(1 / max(npre, 1) + 1 / max(npost, 1))
        ax.text(0.03, 0.97, rf"post/pre = {r:.2f} $\pm$ {rerr:.2f}", transform=ax.transAxes,
                va="top", fontsize=12,
                bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.85))
    if col == 0:
        pty = " (odd parity)" if SIGN == "SS" else ""
        ax.set_ylabel(rf"{SIGN} data{pty} / 10 MeV  (@ 7.99 fb$^{{-1}}$)")
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary",
                      fontsize=11)
    lbl = "no cut" if e >= 1.0 else rf"$\varepsilon_{{sig}}$={100*e:.0f}%"
    ax.set_title(f"{SIGN} data, preEE vs postEE -- {lbl}", fontsize=14, pad=26)
    ax.legend(fontsize=11, loc="upper right")

    ok = (hpre > 0) & (hpost > 0)
    r = np.where(ok, hpost / np.where(hpre > 0, hpre, 1), np.nan)
    # both eras Poisson: rel err^2 = 1/N_pre + 1/N_post_raw
    rerr = np.where(ok, r * np.sqrt(1 / np.maximum(hpre, 1)
                                    + 1 / np.maximum(hpost_raw, 1)), np.nan)
    inside = ok & (r <= RATIO_YMAX)
    rax.errorbar(mcent[inside], r[inside], yerr=rerr[inside], fmt="s", color="black",
                 markersize=4.0, elinewidth=1.2)
    above = ok & (r > RATIO_YMAX)
    if above.any():
        rax.plot(mcent[above], np.full(above.sum(), RATIO_YMAX * 0.97), "^",
                 color="black", markersize=7)
    rax.axhline(1.0, color="firebrick", lw=1.4)
    rax.set_ylim(0, RATIO_YMAX)
    rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
    rax.axvline(M_D, color="green", ls=":", lw=1.5)
    rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    if col == 0:
        rax.set_ylabel("postEE / preEE")
    mask_bad = (~ok & SB) if SIGN == "OS" else ~ok
    n_bad = int(mask_bad.sum())
    notes = [f"{n_bad} bins: no ratio"] if n_bad else []
    if above.any():
        notes.append(f"{int(above.sum())} above axis")
    if notes:
        rax.text(0.98, 0.93, "    ".join(notes), transform=rax.transAxes,
                 fontsize=9, color="firebrick", va="top", ha="right")

out = common.fig_path("era", f"compare_data_{SIGN}_by_era.png", key=MODEL)
plt.savefig(out, dpi=110, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

"""data - W+jets on M(KKpi), OS vs SS: is the data/W+jets difference charge-symmetric?

    python3 draw_os_minus_ss.py [key]     (stage-1 training key, default ssdata)

Columns: no cut and stage-1 (<key> f55) eps_sig = 90/70/50%. Stage-1 only -- stage-2 was
trained on W+jets OS, so after a stage-2 cut OS W+jets is in-sample and SS is not.
  row 1: data - W+jets, OS (black) and SS (blue) overlaid
  row 2: (data - W+jets)_OS - (data - W+jets)_SS; nonzero = a charge-correlated component the
         W+jets MC does not describe (missing process, or W+jets' own charge-correlated part wrong)
The model's own training sample (ssdata: SS data, osdata: OS data, oswjets: OS W+jets) is
used on its odd-event half x2 only (the even half trained it). OS data is blinded in the Ds
window, so that window is dropped from every curve to keep OS and SS on the same footing.
ttbar/DY exist only OS and are NOT subtracted.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
import common  # noqa: E402

M_DS, M_D, SR_HALFWIDTH = 1.96835, 1.8696, 0.04
KEY = sys.argv[1] if len(sys.argv) > 1 else "ssdata"
_cfg = common.resolve(KEY)
PRED = common.out_dir_for(KEY, "DNN") + "/Prediction"
# the sample(s) this model trained on: only their held-out odd half is usable
HELD_OUT = set(_cfg["samples"]) if _cfg["parity"] is not None else set()
EFFS = [1.0, 0.90, 0.70, 0.50]
mbins = np.linspace(M_DS - 0.3, M_DS + 0.3, 61)   # 10 MeV; blind window on bin edges
mcent = 0.5 * (mbins[1:] + mbins[:-1])
SR_BIN = np.abs(mcent - M_DS) < SR_HALFWIDTH
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "event", "mu_charge", "pi_charge",
                           *common.SELECTION]))


def get(name, ss):
    with uproot.open(f"{PRED}/{name}.root")[common.TREE] as t:
        a = t.arrays(WANT, library="np")
    common.check_charge_config(a, ss, name)
    m = common.selection_mask(a)
    w = a["weight"].astype(float)
    scale = 1.0
    if name in HELD_OUT:
        m &= (a["event"] % 2) == 1
        scale = 2.0
    return dict(mass=a["Ds_mass"][m], score=a["score_dsmu"][m], w=w[m] * scale,
                data=name.startswith("Data"), scale=scale)


def hist(s, thr):
    k = s["score"] > thr
    h, _ = np.histogram(s["mass"][k], bins=mbins, weights=s["w"][k])
    if s["data"]:   # Poisson on the raw count, times the odd-half scale
        n, _ = np.histogram(s["mass"][k], bins=mbins)
        e = np.sqrt(n) * s["scale"]
    else:
        e = np.sqrt(np.histogram(s["mass"][k], bins=mbins, weights=s["w"][k] ** 2)[0])
    h, e = h.astype(float), e.astype(float)
    h[SR_BIN] = np.nan
    e[SR_BIN] = np.nan
    return h, e


S = {"dOS": get("Data_OS_merged", False), "dSS": get("Data_SS_merged", True),
     "wOS": get("WtoLNu2Jets_merged", False), "wSS": get("WtoLNu2Jets_SS_merged", True)}

plt.style.use(hep.style.CMS)
fig, axes = plt.subplots(2, len(EFFS), figsize=(8.5 * len(EFFS), 14), squeeze=False)
for col, eff in enumerate(EFFS):
    thr = -np.inf if eff >= 1.0 else common.stage1_threshold(
        {"key": KEY, "tool": "DNN", "eff_sig": eff})[0]
    H = {k: hist(s, thr) for k, s in S.items()}
    lbl = "no cut" if eff >= 1.0 else rf"stage-1 {KEY} $\varepsilon_{{sig}}$={100 * eff:.0f}%"
    R = {}
    ax = axes[0][col]
    for key_d, key_w, color, marker, name in (("dOS", "wOS", "black", "s", "OS"),
                                              ("dSS", "wSS", "steelblue", "o", "SS")):
        if S[key_d]["scale"] != 1.0 or S[key_w]["scale"] != 1.0:
            name += " (odd half x2)"
        r = H[key_d][0] - H[key_w][0]
        er = np.hypot(H[key_d][1], H[key_w][1])
        R[key_d] = (r, er)
        hep.histplot(r, mbins, ax=ax, histtype="errorbar", yerr=er, color=color, marker=marker,
                     markersize=4.5, elinewidth=1.4,
                     label=f"{name}: data $-$ W+jets  N={np.nansum(r):,.0f}")
    ax.axhline(0, color="grey", lw=1)
    ax.set_ylabel("(data $-$ W+jets) Events / 10 MeV")
    ax.set_title(lbl, fontsize=15)

    ax = axes[1][col]
    dd = R["dOS"][0] - R["dSS"][0]
    edd = np.hypot(R["dOS"][1], R["dSS"][1])
    hep.histplot(dd, mbins, ax=ax, histtype="errorbar", yerr=edd, color="darkred", marker="D",
                 markersize=4.5, elinewidth=1.4,
                 label=f"(data$-$W+jets)$_{{OS}}$ $-$ (data$-$W+jets)$_{{SS}}$  N={np.nansum(dd):,.0f}")
    ax.axhline(0, color="grey", lw=1)
    ax.set_ylabel("double difference / 10 MeV")

    for row in (0, 1):
        ax = axes[row][col]
        ax.axvspan(M_DS - SR_HALFWIDTH, M_DS + SR_HALFWIDTH, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.25 * (hi - lo))   # room for legend and peak labels
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97 if top > 0 else top, r"D$^\pm$", color="green", ha="center",
                va="top", fontsize=13)
        ax.text(M_DS, top * 0.97 if top > 0 else top, r"D$_s^\pm$", color="dimgrey",
                ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        ax.legend(fontsize=11, loc="upper right")
    print(f"{lbl}: data-W+jets OS {np.nansum(R['dOS'][0]):,.0f} SS {np.nansum(R['dSS'][0]):,.0f}  double diff {np.nansum(dd):,.0f}")

hep.cms.label(ax=axes[0][0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path("os_ss", "data_minus_wjets_os_ss_mass.png", key=KEY, tool="DNN")
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

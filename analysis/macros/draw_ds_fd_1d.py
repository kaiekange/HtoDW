"""Ds flight distance 1D: Ds_FDxy, Ds_FDxysig, Ds_FD, Ds_FDsig; signal, OS data, W+jets overlaid.

One figure per variable, four panels: no cut and stage-1 eps_sig = 90/70/50%
(score_dsmu > common.stage1_threshold). Rows from the stage-1 model's Prediction/ files,
base selection common.selection_mask & score_dsmu > -1 (as draw_stage1_mass_by_score_cut.py).
- OS data: markers, absolute counts; the Ds SR (|Ds_mass - M_Ds| < 0.04) is blinded by the skim.
- W+jets: filled, absolute (skim weight x common.br_reweight), with the same Ds SR removed so it
  covers the same phase space as the data; hatched band = MC stat.
- Signal: line (skim weight, all rows), scaled to the data integral in each panel -- shape only.
Linear y. Same binning in all panels (0 to the 98% quantile of no-cut data+signal; the last
bin holds the overflow).

    python3 draw_ds_fd_1d.py [--key=ssdata_f45]
Output: <SHARED_FIG_DIR>/features/ds_fd_1d_<var>.png
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

MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata_f45")
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
VARS = [("Ds_FDxy", "Ds flight distance xy [cm]"),
        ("Ds_FDxysig", "Ds flight distance xy significance"),
        ("Ds_FD", "Ds flight distance 3D [cm]"),
        ("Ds_FDsig", "Ds flight distance 3D significance")]
EFFS = [None, 0.90, 0.70, 0.50]
NBINS = 40
M_DS, SR_HALF = 1.96835, 0.040
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
BR_COLS = ["cand_gen_mode", "cand_gen_common_hf_pdgId"]

THR = {e: common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=e))[0]
       for e in EFFS if e is not None}
print("thresholds:", THR)


def read(name, wjets=False, data=False):
    cols = [v for v, _ in VARS] + ["Ds_mass", "score_dsmu", "weight", *common.SELECTION]
    if wjets:
        cols += BR_COLS
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(list(dict.fromkeys(cols)), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    w = np.ones(len(m)) if data else a["weight"].astype(float)
    if wjets:
        w = w * common.br_reweight(a)
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HALF
    out = {v: a[v][m] for v, _ in VARS}
    out["score"], out["w"] = a["score_dsmu"][m], w[m]
    print(f"{name}: {len(m):,} rows -> {m.sum():,} selected")
    return out


SIG = read("Signal_merged")
DAT = read("Data_OS_merged", data=True)
WJ = read("WtoLNu2Jets_merged", wjets=True)


def hist(x, w, edges):
    x = np.clip(x, edges[0], edges[-1] - 1e-9 * (edges[-1] - edges[0]))   # overflow -> last bin
    return np.histogram(x, edges, weights=w)[0], np.histogram(x, edges, weights=w * w)[0]


for var, xlabel in VARS:
    hi = np.quantile(np.concatenate([SIG[var], DAT[var]]), 0.98)
    edges = np.linspace(0.0, hi, NBINS + 1)
    fig, axes = plt.subplots(1, len(EFFS), figsize=(11 * len(EFFS), 9))
    for ax, eff in zip(axes, EFFS):
        sel = {k: (np.ones(len(s["score"]), bool) if eff is None else s["score"] > THR[eff])
               for k, s in (("sig", SIG), ("dat", DAT), ("wj", WJ))}
        hd, _ = hist(DAT[var][sel["dat"]], DAT["w"][sel["dat"]], edges)
        hw, hw2 = hist(WJ[var][sel["wj"]], WJ["w"][sel["wj"]], edges)
        hs, _ = hist(SIG[var][sel["sig"]], SIG["w"][sel["sig"]], edges)
        hs = hs * (hd.sum() / hs.sum())
        hep.histplot(hw, edges, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"W+jets MC (N={hw.sum():,.0f}, Ds SR removed)")
        ew = np.sqrt(hw2)
        ax.bar(edges[:-1], 2 * ew, bottom=hw - ew, width=np.diff(edges), align="edge",
               fill=False, hatch="////", lw=0, edgecolor="grey", label="MC stat.")
        hep.histplot(hs, edges, ax=ax, histtype="step", color="red", lw=2,
                     label="Signal (scaled to data)")
        ctr = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(ctr, hd, yerr=np.sqrt(hd), fmt="s", color="black", ms=5,
                    label=f"Data OS (N={hd.sum():,.0f}, Ds SR blinded)")
        cut_txt = "no cut" if eff is None else f"{MODEL_KEY} $\\epsilon_{{sig}}$={round(100 * eff)}%"
        ax.set_title(cut_txt, fontsize=20, pad=30)
        ax.set_xlabel(xlabel + " (last bin = overflow)")
        ax.set_ylabel(f"Candidates / {edges[1] - edges[0]:.3g}")
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0, 1.35 * max(hd.max(), hw.max(), hs.max()))
        ax.legend(fontsize=14, loc="upper right")
        hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6, fontsize=14)
    out = common.fig_path("features", f"ds_fd_1d_{var}.png")
    plt.tight_layout()
    plt.savefig(out, dpi=80, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

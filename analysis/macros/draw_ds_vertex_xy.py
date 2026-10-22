"""Ds decay vertex (kinematic-fit SV) in the transverse plane: Ds_vy vs Ds_vx, 2D colz.

Signal, OS data (Ds SR blinded by the skim) and W+jets, read from the stage-1 model's
Prediction/ files (same rows as NanoMVA/Samples, plus score_dsmu). Base selection as in
draw_stage1_mass_by_score_cut.py: common.selection_mask & score_dsmu > -1. Cuts: no cut and
score_dsmu > common.stage1_threshold at eps_sig = 90/70/50%. Signal = MC x skim weight (all
rows, no truth match); W+jets = MC x skim weight x common.br_reweight. log z, |x|,|y| < 2 cm.
--rel-pv: Ds vertex - primary vertex (primvtx_vx/vy) instead of absolute coordinates.

    python3 draw_ds_vertex_xy.py [--rel-pv] [--key=ssdata_f45]
Output: <SHARED_FIG_DIR>/features/ds_vertex_xy[_relpv]_<nocut|eff90|eff70|eff50>.png
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

REL_PV = "--rel-pv" in sys.argv
MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata_f45")
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
SAMPLES = [("Signal_merged", "Signal MC (x weight)", False),
           ("Data_OS_merged", "Data OS (Ds SR blinded)", True),
           ("WtoLNu2Jets_merged", "W+jets MC (x weight x BR rw)", False)]
EFFS = [None, 0.90, 0.70, 0.50]
HW, NB = 2.0, 200   # half-width [cm], bins per axis
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
BR_COLS = ["cand_gen_mode", "cand_gen_common_hf_pdgId"]

THR = {e: common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=e))[0]
       for e in EFFS if e is not None}
print("thresholds:", THR)

# per sample: x, y, score, weight after the base selection
DATA = {}
for name, _label, is_data in SAMPLES:
    cols = ["Ds_vx", "Ds_vy", "primvtx_vx", "primvtx_vy", "score_dsmu", "weight", *common.SELECTION]
    if name.startswith("WtoLNu2Jets"):
        cols += BR_COLS
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(list(dict.fromkeys(cols)), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    w = np.ones(len(m)) if is_data else a["weight"].astype(float)
    if name.startswith("WtoLNu2Jets"):
        w = w * common.br_reweight(a)
    x, y = a["Ds_vx"], a["Ds_vy"]
    if REL_PV:
        x, y = x - a["primvtx_vx"], y - a["primvtx_vy"]
    DATA[name] = (x[m], y[m], a["score_dsmu"][m], w[m])
    print(f"{name}: {len(m):,} rows -> {m.sum():,} after selection_mask & score>-1")

edges = np.linspace(-HW, HW, NB + 1)
for eff in EFFS:
    tag = "nocut" if eff is None else f"eff{round(100 * eff)}"
    fig, axes = plt.subplots(1, len(SAMPLES), figsize=(10 * len(SAMPLES), 9))
    for ax, (name, label, is_data) in zip(axes, SAMPLES):
        x, y, s, w = DATA[name]
        k = np.ones(len(s), bool) if eff is None else s > THR[eff]
        h, _, _ = np.histogram2d(x[k], y[k], bins=[edges, edges], weights=w[k])
        inside = h.sum() / w[k].sum()
        hm = np.where(h > 0, h, np.nan)
        mesh = ax.pcolormesh(edges, edges, hm.T, norm=LogNorm(vmin=np.nanmin(hm), vmax=np.nanmax(hm)),
                             cmap="viridis", rasterized=True)
        fig.colorbar(mesh, ax=ax, label="Candidates / bin" if is_data else "Events / bin")
        ax.set_xlabel("Ds vertex x - PV x [cm]" if REL_PV else "Ds vertex x [cm]")
        ax.set_ylabel("Ds vertex y - PV y [cm]" if REL_PV else "Ds vertex y [cm]")
        ax.set_aspect("equal")
        cut_txt = "no cut" if eff is None else f"{MODEL_KEY} $\\epsilon_{{sig}}$={round(100 * eff)}%"
        ax.text(0.03, 0.97, f"{label}\n{cut_txt}\n{k.sum():,} rows, {100 * inside:.1f}% inside",
                transform=ax.transAxes, va="top", ha="left", fontsize=15,
                bbox=dict(facecolor="white", alpha=0.8, lw=0))
        hep.cms.label(ax=ax, data=is_data, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6,
                      fontsize=14)
    out = common.fig_path("features", f"ds_vertex_xy{'_relpv' if REL_PV else ''}_{tag}.png")
    plt.tight_layout()
    plt.savefig(out, dpi=80, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

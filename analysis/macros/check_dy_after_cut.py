"""Why does DY survive the stage-1 cut 3-7x better than W+jets? (2026-10-01)

Adding DY_merged to the Ds_mass closure removed the low-mass tilt but left MC ~1.5x
above data after the cut. The surviving DY is 91% gen-unmatched combinatorial. Working
hypothesis: Z->mumu with one muon lost -> the lost muon fakes MET -> the event looks
like W->munu. This compares OS data against W+jets / ttbar / DY (stacked, absolute
sigma*L yields) in event-level variables at a stage-1 working point. Ds-mass SR removed
from data AND MC, so both sides are the same event set.

    conda activate ds_dnn && python3 check_dy_after_cut.py [eff_sig]   # default 0.90
"""
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common  # noqa: E402

EFF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.90
PRED = f"{common.out_dir_for('ssdata')}/Prediction"
M_DS, SR_HALF = 1.96835, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
THR = common.stage1_threshold(dict(key="ssdata", tool="DNN", seed=42, eff_sig=EFF))[0]

VARS = [   # (branch or derived name, xlabel, bins)
    ("W_tmass_corr", r"$m_T(\mu, p_T^{miss})$ [GeV]", np.arange(0, 160.1, 8)),
    ("met_pt_corr", r"$p_T^{miss}$ [GeV]", np.arange(0, 120.1, 6)),
    ("mu_pt_corr", r"$p_T(\mu)$ [GeV]", np.arange(26, 126.1, 5)),
    ("n_loose_mu", r"N(loose $\mu$), os+ss", np.arange(-0.5, 5.5, 1)),
    ("nElectron", r"N(electron), nano", np.arange(-0.5, 5.5, 1)),
    ("num_jets_pt30_corr", r"N(jet, $p_T>30$)", np.arange(-0.5, 6.5, 1)),
    ("Dsmu_mass_corr", r"$M(D_s\mu)$ [GeV]", np.arange(0, 125.1, 5)),
    ("Ds_pt", r"$p_T(D_s)$ [GeV]", np.arange(7, 57.1, 2.5)),
]
BR = sorted({v for v, _x, _b in VARS if v != "n_loose_mu"}
            | {"num_os_mu_loose", "num_ss_mu_loose", "Ds_mass", "score_dsmu"}
            | set(common.SELECTION))
MC = [("W+jets", "WtoLNu2Jets_merged", "#5a9e78"), (r"$t\bar{t}$", "TT_merged", "#7f5fbf"),
      ("DY", "DY_merged", "#4a90c2")]


def load(name, weighted):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        BR + (["weight"] if weighted else []), library="np")
    m = (common.selection_mask(a) & (a["score_dsmu"] > THR)
         & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALF))
    out = {k: v[m] for k, v in a.items()}
    out["n_loose_mu"] = out["num_os_mu_loose"] + out["num_ss_mu_loose"]
    out["w"] = out["weight"].astype(float) if weighted else np.ones(m.sum())
    return out


d = load("Data_OS_merged", False)
mc = [(lab, load(f, True), col) for lab, f, col in MC]
print(f"eps_sig={EFF:.2f} thr={THR:.6f}  data rows {len(d['w']):,}  "
      + "  ".join(f"{lab} rows {len(x['w']):,} sumw {x['w'].sum():,.1f}" for lab, x, _c in mc))

fig, axes = plt.subplots(2, 4, figsize=(26, 11))
for ax, (var, xl, bins) in zip(axes.flat, VARS):
    c = 0.5 * (bins[1:] + bins[:-1])
    bottom = np.zeros(len(c))
    tot = np.zeros(len(c))
    for lab, x, col in mc:
        h = np.histogram(np.clip(x[var], bins[0], bins[-1] - 1e-6), bins=bins, weights=x["w"])[0]
        ax.bar(c, h, width=np.diff(bins), bottom=bottom, color=col, alpha=0.9,
               label=f"{lab}  {h.sum():,.0f}", edgecolor="none")
        bottom += h
        tot += h
    hd = np.histogram(np.clip(d[var], bins[0], bins[-1] - 1e-6), bins=bins)[0]
    ax.errorbar(c, hd, yerr=np.sqrt(hd), fmt="o", ms=4, color="black",
                label=f"OS data  {hd.sum():,}")
    ax.set_xlabel(xl)
    ax.set_ylabel("Events (last bin = overflow)")
    ax.legend(fontsize=9)
    means = "  ".join(f"{lab} {np.average(x[var], weights=x['w']):.2f}" for lab, x, _c in mc)
    print(f"  {var:20}  mean: data {d[var].mean():.2f}  {means}   data/MC per bin: "
          + " ".join(f"{a/b:.2f}" if b > 0 else "-" for a, b in zip(hd, tot)))
hep.cms.label(ax=axes[0][0], data=True, text="Preliminary",
              rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)
fig.suptitle(rf"stage-1 $\varepsilon_{{sig}}$={EFF:.0%}, Ds-mass SR removed -- OS data vs W+jets + $t\bar{{t}}$ + DY",
             fontsize=15)
out = common.fig_path("closure", f"dy_after_cut_eff{int(round(EFF * 100))}.png", key="ssdata")
fig.savefig(out, dpi=110, bbox_inches="tight")
print("Saved", out)

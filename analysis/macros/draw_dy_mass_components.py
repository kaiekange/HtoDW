"""DY Ds_mass / phi_mass after the SS-trained (ssdata) vs OS-trained (osdata) cut, per component.

    python3 draw_dy_mass_components.py [KEY_OS]      (default osdata; compared with ssdata)

Same layout as draw_mass_ss_vs_os_dnn.py (ds_mass_ssdata_vs_osdata_DY.png): each panel overlays
the two cuts on the SAME rows, sub-panel KEY_OS / ssdata per bin with the overlap-aware error.
One row per component: GenTag mode groups (cand_gen_mode), then DY jet bins (jet_bin).
Columns eps_sig = 90/70/50%, each model at its own stage1_threshold. DY is in neither model's
training set -> all rows. MC, so the Ds window is not blinded; the grey band marks the SR.
Second set: mode 0 only, one row per jet bin.
Output: os_ss/dy_{ds,phi}_mass_{components,mode0_by_jetbin}_ssdata_vs_<KEY_OS>.png (shared dir).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
import common  # noqa: E402

BASE = "ssdata"
KEY_OS = sys.argv[1] if len(sys.argv) > 1 else "osdata"
EFFS = (0.90, 0.70, 0.50)
M_DS, M_D, M_PHI, SR_HALF = 1.96835, 1.8696, 1.019460, 0.04
PHI_RANGE = (0.984, 1.074)   # clears 2m_K and the M_phi+0.05 post-fit cut
# var: (x label, bin width [GeV] per eps, range)
VARS = {"Ds_mass": (r"$M(K^+K^-\pi^\pm)$ [GeV]", {0.90: 0.020, 0.70: 0.030, 0.50: 0.040},
                    (M_DS - 0.3, M_DS + 0.3)),
        "phi_mass": (r"$M(K^+K^-)$ [GeV]", {0.90: 0.003, 0.70: 0.0045, 0.50: 0.006}, PHI_RANGE)}
COMPONENTS = [   # (label, {branch: values}) -- a row passes if every branch is in its values
    ("comb. (mode 0)", {"cand_gen_mode": [0]}),
    ("partial (mode 8)", {"cand_gen_mode": [8]}),
    ("multi-charm (mode 7)", {"cand_gen_mode": [7]}),
    ("other (mode 6,9,10)", {"cand_gen_mode": [6, 9, 10]}),
    (r"D$\to$K$\pi\pi$ misID (mode 5)", {"cand_gen_mode": [5]}),
    (r"genuine D$\to$KK$\pi$ (mode 3+4)", {"cand_gen_mode": [3, 4]}),
    (r"true D$_s$ (mode 1+2)", {"cand_gen_mode": [1, 2]}),
    ("0J", {"jet_bin": [0]}),
    ("1J", {"jet_bin": [1]}),
    ("2J", {"jet_bin": [2]}),
]
MODE0_BY_JETBIN = [(f"{j}J, comb. (mode 0)", {"jet_bin": [j], "cand_gen_mode": [0]}) for j in (0, 1, 2)]
FIGURES = {"components": COMPONENTS, "mode0_by_jetbin": MODE0_BY_JETBIN}
WANT = list(dict.fromkeys(["Ds_mass", "phi_mass", "weight", "score_dsmu", "cand_gen_mode", "jet_bin",
                           "mu_charge", "pi_charge", *common.SELECTION]))

thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in (BASE, KEY_OS)}
arr = {}
for k in (BASE, KEY_OS):
    with uproot.open(f"{common.out_dir_for(k, 'DNN')}/Prediction/DY_merged.root")[common.TREE] as t:
        arr[k] = t.arrays(WANT, library="np")
a0 = arr[BASE]
assert np.array_equal(arr[KEY_OS]["Ds_mass"], a0["Ds_mass"]), "DY_merged: row order differs"
common.check_charge_config(a0, False, "DY_merged")
sel = common.selection_mask(a0)
w = a0["weight"][sel].astype(float)
pass_b = {e: arr[BASE]["score_dsmu"][sel] > thr[BASE][e] for e in EFFS}
pass_o = {e: arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e] for e in EFFS}
plt.style.use(hep.style.CMS)


def h(x, wt, bins):
    return np.histogram(x, bins=bins, weights=wt)[0], np.histogram(x, bins=bins, weights=wt ** 2)[0]


print(f"DY_merged: {sel.sum():,} rows after selection")
for (ftag, rows), (var, (xlabel, step, (xlo, xhi))) in ((f, v) for f in FIGURES.items() for v in VARS.items()):
    x = a0[var][sel]
    nrow = len(rows)
    fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
    outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25,
                              top=1 - 0.9 / (8.0 * nrow))
    print(f"\n### {ftag}: {var}")
    for row, (clab, cond) in enumerate(rows):
        g = np.ones(len(w), bool)
        for branch, vals in cond.items():
            g &= np.isin(a0[branch][sel], vals)
        for col, e in enumerate(EFFS):
            bins = np.linspace(xlo, xhi, int(round((xhi - xlo) / step[e])) + 1)
            bc = 0.5 * (bins[1:] + bins[:-1])
            inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col],
                                                     height_ratios=[3, 1], hspace=0.05)
            ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            pb, po = g & pass_b[e], g & pass_o[e]
            hb, vb = h(x[pb], w[pb], bins)
            ho, vo = h(x[po], w[po], bins)
            hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)),
                         color="steelblue", marker="o", markersize=5, elinewidth=1.3,
                         label=f"{BASE} (SS DNN) cut  N={hb.sum():,.0f}")
            hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)),
                         color="black", marker="s", markersize=5, elinewidth=1.3,
                         label=f"{KEY_OS} (OS DNN) cut  N={ho.sum():,.0f}")
            if var == "Ds_mass":
                ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
                rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
                ax.axvline(M_D, color="green", ls=":", lw=1.5)
            else:
                ax.axvline(M_PHI, color="green", ls=":", lw=1.5)
            ax.set_xlim(bins[0], bins[-1])
            top = max(hb.max(), ho.max())
            ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * top if top > 0 else 1.0)
            ax.set_ylabel(f"Events / {step[e] * 1000:g} MeV")
            ax.set_title(rf"DY {clab}, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
            ax.legend(fontsize=11, loc="upper right")
            ax.tick_params(labelbottom=False)
            # per-bin ratio, overlap-aware: (C+A)/(C+B) with C both pass, A only OS, B only SS
            c, vc = h(x[pb & po], w[pb & po], bins)
            a, va = h(x[po & ~pb], w[po & ~pb], bins)
            b, vbb = h(x[pb & ~po], w[pb & ~po], bins)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = (c + a) / (c + b)
                er = np.sqrt(vc * ((b - a) / (c + b) ** 2) ** 2 + va / (c + b) ** 2
                             + vbb * ((c + a) / (c + b) ** 2) ** 2)
            rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
            rax.axhline(1.0, color="red", lw=1.0)
            rax.set_ylim(0.0, 1.99)
            rax.set_yticks([0.5, 1.0, 1.5])
            rax.set_xlim(bins[0], bins[-1])
            rax.set_xlabel(xlabel)
            rax.set_ylabel("OS / SS")
            if var == "Ds_mass":
                with np.errstate(divide="ignore", invalid="ignore"):
                    print(f"  {clab.replace('$', ''):34s} eps{100 * e:.0f}: SS-DNN {hb.sum():9,.1f}"
                          f"  OS-DNN {ho.sum():9,.1f}  OS/SS {ho.sum() / hb.sum():.3f}")
    hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
    fig.suptitle(rf"DY$\to\ell\ell$ (0J+1J+2J), OS: {var} by component after SS-trained ({BASE}) "
                 f"vs OS-trained ({KEY_OS}) stage-1 cut", fontsize=18, y=0.995)
    out = common.fig_path("os_ss", f"dy_{var.split('_')[0].lower()}_mass_{ftag}_{BASE}_vs_{KEY_OS}.png")
    plt.savefig(out, dpi=80, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

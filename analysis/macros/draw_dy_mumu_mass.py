"""DY: M(mu, pi as mu) after the SS-trained (ssdata) vs OS-trained (osdata) stage-1 cut, overlaid.

    python3 draw_dy_mumu_mass.py [KEY_OS]      (default osdata; compared with ssdata)

The Ds candidate's pion is given the muon mass and combined with the best muon (corrected pT,
mu_pt_corr) -- if the pion were the Z's other muon this would peak at m_Z. Same layout as
draw_mass_ss_vs_os_dnn.py (ds_mass_ssdata_vs_osdata_DY.png): overlay of the two cuts on the
SAME rows, sub-panel KEY_OS / ssdata per bin with the overlap-aware error. Columns eps_sig =
90/70/50%. DY is in neither model's training set -> all rows.
Output: os_ss/dy_mumu_mass_ssdata_vs_<KEY_OS>.png (shared dir).
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
BIN_GEV = {0.90: 5.0, 0.70: 7.5, 0.50: 10.0}
XMAX, M_MU, M_Z = 150.0, 0.1056583755, 91.1876
WANT = list(dict.fromkeys(["mu_pt_corr", "mu_eta", "mu_phi", "pi_pt", "pi_eta", "pi_phi", "Ds_mass",
                           "weight", "score_dsmu", "mu_charge", "pi_charge", *common.SELECTION]))


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2), px, py, pz


def h(x, w, bins):
    return np.histogram(x, bins=bins, weights=w)[0], np.histogram(x, bins=bins, weights=w ** 2)[0]


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
E1, x1, y1, z1 = p4(a0["mu_pt_corr"][sel], a0["mu_eta"][sel], a0["mu_phi"][sel], M_MU)
E2, x2, y2, z2 = p4(a0["pi_pt"][sel], a0["pi_eta"][sel], a0["pi_phi"][sel], M_MU)
mass = np.sqrt(np.maximum((E1 + E2) ** 2 - (x1 + x2) ** 2 - (y1 + y2) ** 2 - (z1 + z2) ** 2, 0.0))
w = a0["weight"][sel].astype(float)
print(f"DY_merged: {sel.sum():,} rows; M(mu,pi->mu) > {XMAX:.0f} GeV (overflow, not drawn): "
      f"{(mass >= XMAX).mean():.3%}")

plt.style.use(hep.style.CMS)
fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0))
outer = gridspec.GridSpec(1, len(EFFS), figure=fig, hspace=0.32, wspace=0.25, top=0.88)
for col, e in enumerate(EFFS):
    bins = np.linspace(0.0, XMAX, int(round(XMAX / BIN_GEV[e])) + 1)
    bc = 0.5 * (bins[1:] + bins[:-1])
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, col], height_ratios=[3, 1],
                                             hspace=0.05)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    pb = arr[BASE]["score_dsmu"][sel] > thr[BASE][e]
    po = arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e]
    hb, vb = h(mass[pb], w[pb], bins)
    ho, vo = h(mass[po], w[po], bins)
    hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)), color="steelblue",
                 marker="o", markersize=5, elinewidth=1.3, label=f"{BASE} (SS DNN) cut  N={hb.sum():,.0f}")
    hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)), color="black",
                 marker="s", markersize=5, elinewidth=1.3, label=f"{KEY_OS} (OS DNN) cut  N={ho.sum():,.0f}")
    ax.axvline(M_Z, color="green", ls=":", lw=1.5)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * max(hb.max(), ho.max()))
    top = ax.get_ylim()[1]
    ax.text(M_Z + 2, top * 0.75, r"m$_Z$", color="green", ha="left", va="top", fontsize=13)
    ax.set_ylabel(f"Events / {BIN_GEV[e]:g} GeV")
    ax.set_title(rf"DY all modes, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
    ax.legend(fontsize=11, loc="upper right")
    ax.tick_params(labelbottom=False)
    # per-bin ratio, overlap-aware: (C+A)/(C+B) with C both pass, A only OS, B only SS
    c, vc = h(mass[pb & po], w[pb & po], bins)
    a, va = h(mass[po & ~pb], w[po & ~pb], bins)
    b, vbb = h(mass[pb & ~po], w[pb & ~po], bins)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (c + a) / (c + b)
        er = np.sqrt(vc * ((b - a) / (c + b) ** 2) ** 2 + va / (c + b) ** 2
                     + vbb * ((c + a) / (c + b) ** 2) ** 2)
    rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
    rax.axhline(1.0, color="red", lw=1.0)
    rax.set_ylim(0.0, 1.99)
    rax.set_yticks([0.5, 1.0, 1.5])
    rax.set_xlim(bins[0], bins[-1])
    rax.set_xlabel(r"$M(\mu, \pi\to\mu)$ [GeV]")
    rax.set_ylabel("OS / SS")
    zwin = np.abs(bc - M_Z) < 15
    print(f"  eps{100 * e:.0f}: SS-DNN {hb.sum():8,.1f}  OS-DNN {ho.sum():8,.1f}  OS/SS {ho.sum() / hb.sum():.3f}"
          f"  | |M-mZ|<15: SS-DNN {hb[zwin].sum():7,.1f}  OS-DNN {ho[zwin].sum():7,.1f}")
hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
fig.suptitle(rf"DY$\to\ell\ell$ (0J+1J+2J), OS: M($\mu$, $\pi\to\mu$) after SS-trained ({BASE}) "
             f"vs OS-trained ({KEY_OS}) stage-1 cut", fontsize=18, y=0.995)
out = common.fig_path("os_ss", f"dy_mumu_mass_{BASE}_vs_{KEY_OS}.png")
plt.savefig(out, dpi=100, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

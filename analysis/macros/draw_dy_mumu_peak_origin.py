"""Why DY M(mu, pi as mu) peaks near 50 GeV after the stage-1 cut: shape comparison.

    python3 draw_dy_mumu_peak_origin.py [KEY_OS]      (default osdata; with ssdata)

Unit-area shapes of M(mu,pi->mu), M(Ds,mu), Ds pT, pi pT, dphi(mu,pi) and M(mu,pi)/M(Ds,mu) for
DY with no cut, DY after the ssdata and KEY_OS eps_sig=90% cuts, and the signal (no cut).
DY is in neither model's training set -> all rows; signal is not cut (its eps90 shape is ~ the same).
Output: os_ss/dy_mumu_peak_origin_ssdata_<KEY_OS>.png (shared dir).
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

BASE = "ssdata"
KEY_OS = sys.argv[1] if len(sys.argv) > 1 else "osdata"
EFF = 0.90
M_MU = 0.1056583755
WANT = list(dict.fromkeys(["mu_pt_corr", "mu_eta", "mu_phi", "pi_pt", "pi_eta", "pi_phi", "Ds_pt",
                           "Dsmu_mass_corr", "Ds_mass", "weight", "score_dsmu", "mu_charge", "pi_charge",
                           *common.SELECTION]))
VARS = [   # (key, x label, bin edges)
    ("mumu", r"$M(\mu, \pi\to\mu)$ [GeV]", np.linspace(0, 120, 41)),
    ("dsmu", r"$M(D_s, \mu)$ [GeV]", np.linspace(0, 140, 36)),
    ("ds_pt", r"$p_T(D_s)$ [GeV]", np.linspace(0, 80, 41)),
    ("pi_pt", r"$p_T(\pi)$ [GeV]", np.linspace(0, 50, 51)),
    ("dphi", r"$\Delta\phi(\mu, \pi)$", np.linspace(0, np.pi, 33)),
    ("frac", r"$M(\mu, \pi\to\mu)\,/\,M(D_s, \mu)$", np.linspace(0, 1.2, 41)),
]


def p4(pt, eta, phi, m):
    return np.sqrt((pt * np.cosh(eta)) ** 2 + m * m), pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)


def derive(a):
    E1, x1, y1, z1 = p4(a["mu_pt_corr"], a["mu_eta"], a["mu_phi"], M_MU)
    E2, x2, y2, z2 = p4(a["pi_pt"], a["pi_eta"], a["pi_phi"], M_MU)
    mumu = np.sqrt(np.maximum((E1 + E2) ** 2 - (x1 + x2) ** 2 - (y1 + y2) ** 2 - (z1 + z2) ** 2, 0.0))
    dphi = np.abs((a["mu_phi"] - a["pi_phi"] + np.pi) % (2 * np.pi) - np.pi)
    return {"mumu": mumu, "dsmu": a["Dsmu_mass_corr"], "ds_pt": a["Ds_pt"], "pi_pt": a["pi_pt"],
            "dphi": dphi, "frac": mumu / a["Dsmu_mass_corr"]}


def load(key, fname):
    with uproot.open(f"{common.out_dir_for(key, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
        a = t.arrays(WANT, library="np")
    common.check_charge_config(a, False, fname)
    s = common.selection_mask(a)
    return {k: v[s] for k, v in a.items()}


dy = {k: load(k, "DY_merged") for k in (BASE, KEY_OS)}
assert np.array_equal(dy[BASE]["Ds_mass"], dy[KEY_OS]["Ds_mass"]), "DY_merged: row order differs"
sig = load(BASE, "Signal_merged")
thr = {k: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": EFF})[0] for k in (BASE, KEY_OS)}
d_dy, d_sig = derive(dy[BASE]), derive(sig)
w_dy, w_sig = dy[BASE]["weight"].astype(float), sig["weight"].astype(float)
CURVES = [   # (label, values dict, weights, row mask, colour, linestyle)
    ("DY, no cut", d_dy, w_dy, np.ones(len(w_dy), bool), "grey", "-"),
    (rf"DY, {BASE} (SS DNN) $\varepsilon_{{sig}}$={100 * EFF:.0f}%", d_dy, w_dy,
     dy[BASE]["score_dsmu"] > thr[BASE], "steelblue", "-"),
    (rf"DY, {KEY_OS} (OS DNN) $\varepsilon_{{sig}}$={100 * EFF:.0f}%", d_dy, w_dy,
     dy[KEY_OS]["score_dsmu"] > thr[KEY_OS], "black", "-"),
    ("signal (DsW+DsstW), no cut", d_sig, w_sig, np.ones(len(w_sig), bool), "#c0392b", "--"),
]

plt.style.use(hep.style.CMS)
fig, axes = plt.subplots(2, 3, figsize=(25, 15))
for ax, (v, xlabel, bins) in zip(axes.flat, VARS):
    for lab, d, w, m, col, ls in CURVES:
        x, ww = d[v][m], w[m]
        h, _ = np.histogram(x, bins=bins, weights=ww)
        h = h / ww.sum() / np.diff(bins)
        o = np.argsort(x)
        c = np.cumsum(ww[o]) / ww.sum()
        med = x[o][np.searchsorted(c, 0.5)]
        hep.histplot(h, bins, ax=ax, histtype="step", color=col, ls=ls, lw=2.2,
                     label=f"{lab}  (median {med:.2f})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("normalized to unit area")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.45 * ax.get_ylim()[1])
    ax.legend(fontsize=12, loc="upper right" if v != "dphi" else "upper left")
hep.cms.label(ax=axes[0][0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=13)
fig.suptitle(r"DY$\to\ell\ell$ (0J+1J+2J), OS: the stage-1 cut selects signal-like topology", fontsize=20)
plt.tight_layout(rect=(0, 0, 1, 0.97))
out = common.fig_path("os_ss", f"dy_mumu_peak_origin_{BASE}_{KEY_OS}.png")
plt.savefig(out, dpi=90)
plt.close(fig)
print(f"Saved {out}")

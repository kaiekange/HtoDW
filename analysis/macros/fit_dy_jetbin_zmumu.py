"""Z->mumu CR: one normalisation factor per DY generator jet-bin sample (0J/1J/2J), fitted on
the full pT(mumu) spectrum, then every other variable drawn post-fit.

Selection: zmumu CR skim (make_samples_nano.py --region=zmumu), 76 < M(mumu) < 106 and the
analysis phi window |phi_mass - M_phi| < 20 MeV (always; the skim has none); no DNN cut (at eps90 the
window keeps ~400 data events, too few -- user 2026-10-19).

Why pT(mumu): the three DY samples share one M(mumu) shape (the Z peak), so M(mumu) cannot
tell them apart; their pT(mumu) spectra differ (pT<5: 96% 0J, pT>130: ~80% 2J). One fit over
the whole pT(mumu) range, no categories (user 2026-10-19). ttbar and W+jets are fixed MC
(<0.6% of the window). Chi2 with data Poisson + MC-template stat errors, least squares;
the single-factor fit is printed for comparison.

Panels: pT(mumu) pre-fit and post-fit (fitted), then post-fit M(mumu), pT(mu1), pT(mu2),
reco Njet, MET and m_T -- NOT fitted, they show whether the three factors describe them.

    python3 fit_dy_jetbin_zmumu.py [--dir=<CR dir>] [--key=ssdata_phi20]
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import uproot
from scipy.optimize import least_squares

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

ARGS = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
CR_DIR = ARGS.get("dir", f"{common.SAMPLES_DIR}/CR/zmumu")
MODEL_KEY = ARGS.get("key", "ssdata_phi20")
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
PTLL_BINS = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200.])
VARS = [  # (column, x label, bins, per GeV)
    ("ptll", r"$p_T(\mu\mu)$ [GeV]", PTLL_BINS, True),
    ("m_mumu", r"$M(\mu\mu)$ [GeV]", np.arange(ZLO, ZHI + 1e-9, 1.0), False),
    ("mu_pt_corr", r"$p_T(\mu_1)$ [GeV]", np.arange(25.0, 150.0 + 1e-9, 5.0), False),
    ("xmu_pt", r"$p_T(\mu_2)$ [GeV]", np.arange(15.0, 120.0 + 1e-9, 5.0), False),
    ("num_jets_pt30_corr", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 6.5, 1.0), False),
    ("met_pt_corr", r"$p_T^{miss}$ [GeV]", np.arange(0.0, 150.0 + 1e-9, 5.0), False),
    ("W_tmass_corr", r"$m_T$ [GeV]", np.arange(0.0, 150.0 + 1e-9, 5.0), False),
]
COLS = ["m_mumu", "mu_pt_corr", "mu_phi", "xmu_pt", "xmu_phi", "weight", "phi_mass",
        "num_jets_pt30_corr", "met_pt_corr", "W_tmass_corr"]
DY_COL = {0: "#9ecae1", 1: "#4a90c2", 2: "#1f4e79"}
OTHER_COL = "#7f5fbf"


def read(g, extra=()):
    with uproot.open(f"{CR_DIR}/{g}.root")["Events"] as t:
        a = t.arrays(COLS + list(extra), library="np")
    z = ((a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
         & (np.abs(a["phi_mass"] - M_PHI) < 0.020))
    a = {k: v[z] for k, v in a.items()}
    px = a["mu_pt_corr"] * np.cos(a["mu_phi"]) + a["xmu_pt"] * np.cos(a["xmu_phi"])
    py = a["mu_pt_corr"] * np.sin(a["mu_phi"]) + a["xmu_pt"] * np.sin(a["xmu_phi"])
    a["ptll"] = np.hypot(px, py)
    a["weight"] = a["weight"].astype(float)
    return a


data = read("Data_OS_merged")
data["weight"] = np.ones(len(data["m_mumu"]))
dy = read("DY_merged", ["jet_bin"])
other = [read("TT_merged"), read("WtoLNu2Jets_merged")]


def h(a, var, bins, m=None, w2=False):
    m = np.ones(len(a[var]), bool) if m is None else m
    w = a["weight"][m]
    return np.histogram(np.clip(a[var][m], bins[0], bins[-1] - 1e-9), bins,
                        weights=w * w if w2 else w)[0]


# ---- fit on the full pT(mumu) spectrum (last bin holds the overflow) ----
hd = h(data, "ptll", PTLL_BINS)
hn = sum(h(o, "ptll", PTLL_BINS) for o in other)
en = sum(h(o, "ptll", PTLL_BINS, w2=True) for o in other)
T = [h(dy, "ptll", PTLL_BINS, dy["jet_bin"] == j) for j in (0, 1, 2)]
E2 = [h(dy, "ptll", PTLL_BINS, dy["jet_bin"] == j, True) for j in (0, 1, 2)]


def fit(nfac):
    def resid(k):
        k = list(k) * 3 if nfac == 1 else k
        pred = hn + sum(kk * t for kk, t in zip(k, T))
        var = hd + en + sum(kk * kk * e for kk, e in zip(k, E2))
        return (hd - pred) / np.sqrt(np.maximum(var, 1.0))
    r = least_squares(resid, [0.8] * nfac)
    cov = np.linalg.inv(r.jac.T @ r.jac)
    return r.x, np.sqrt(np.diag(cov)), cov, float((r.fun ** 2).sum()), len(hd) - nfac


K, EK, cov, chi3, ndf3 = fit(3)
K1, EK1, _, chi1, ndf1 = fit(1)
print(f"fit on pT(mumu), phi20, no DNN cut, "
      f"data {int(hd.sum()):,}")
print(f"  3-factor: k0J {K[0]:.4f}+-{EK[0]:.4f}  k1J {K[1]:.4f}+-{EK[1]:.4f}  "
      f"k2J {K[2]:.4f}+-{EK[2]:.4f}  chi2/ndf {chi3:.1f}/{ndf3}")
print("  correlation:\n", (cov / np.outer(EK, EK)).round(2))
print(f"  1-factor: k {K1[0]:.4f}+-{EK1[0]:.4f}  chi2/ndf {chi1:.1f}/{ndf1}")

# ---- pictures ----
PANELS = [("pre-fit", VARS[0], (1.0, 1.0, 1.0)), ("post-fit (fitted)", VARS[0], tuple(K))]
PANELS += [("post-fit (not fitted)", v, tuple(K)) for v in VARS[1:]]
ncol = 4
nrow = (len(PANELS) + ncol - 1) // ncol
fig = plt.figure(figsize=(7.4 * ncol, 8.0 * nrow))
outer = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.30, hspace=0.30, top=0.92, bottom=0.07)
for i, (title, (var, xlab, bins, per_gev), k) in enumerate(PANELS):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i],
                                             height_ratios=[3, 1], hspace=0.06)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    norm = wd if per_gev else np.ones_like(wd)
    comps = [(r"$t\bar{t}$ + W+jets", sum(h(o, var, bins) for o in other),
              sum(h(o, var, bins, w2=True) for o in other), OTHER_COL)]
    for j in (0, 1, 2):
        m = dy["jet_bin"] == j
        lab = f"DY {j}J" + (f" x{k[j]:.3f}" if k[j] != 1.0 else "")
        comps.append((lab, k[j] * h(dy, var, bins, m), k[j] ** 2 * h(dy, var, bins, m, True),
                      DY_COL[j]))
    bottom = np.zeros(len(bc))
    tot, tot2 = np.zeros(len(bc)), np.zeros(len(bc))
    for lab, hh, e2, colr in comps:
        ax.bar(bc, hh / norm, width=wd, bottom=bottom, color=colr, edgecolor="none",
               label=f"{lab}  N={hh.sum():,.0f}")
        bottom += hh / norm; tot += hh; tot2 += e2
    terr = np.sqrt(tot2)
    ax.bar(bc, 2 * terr / norm, width=wd, bottom=(tot - terr) / norm, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, label="MC stat")
    dh = h(data, var, bins)
    ax.errorbar(bc, dh / norm, yerr=np.sqrt(dh) / norm, fmt="o", ms=4, color="black",
                label=f"OS data  N={dh.sum():,.0f}")
    ax.set_ylim(0, 1.45 * max((tot / norm).max(), (dh / norm).max()))
    ax.set_ylabel("Events / GeV" if per_gev else "Events / bin", fontsize=16)
    ax.set_xlim(bins[0], bins[-1])
    ax.tick_params(labelbottom=False)
    ax.legend(fontsize=10, loc="upper right")
    ax.set_title(title, fontsize=15, pad=24 if i == 0 else 6)
    with np.errstate(invalid="ignore", divide="ignore"):
        r, re, band = dh / tot, np.sqrt(dh) / tot, terr / tot
    rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////",
            edgecolor="0.35", linewidth=0.0)
    rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
    rax.axhline(1, color="red", lw=1)
    rax.set_ylim(0.5, 1.5)
    rax.set_xlim(bins[0], bins[-1])
    rax.set_xlabel(xlab)
    rax.set_ylabel("data/MC", fontsize=14)
    if i == 0:
        hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}",
                      com=13.6, fontsize=14)

sel_txt = r"$|M_{KK}-M_\phi|$<20 MeV, "
fig.suptitle(rf"Z$\to\mu\mu$ CR, 76<M($\mu\mu$)<106, {sel_txt}no DNN cut: DY 0J/1J/2J fitted on "
             rf"$p_T(\mu\mu)$ -- k$_{{0J}}$={K[0]:.3f}$\pm${EK[0]:.3f}, "
             rf"k$_{{1J}}$={K[1]:.3f}$\pm${EK[1]:.3f}, k$_{{2J}}$={K[2]:.3f}$\pm${EK[2]:.3f}, "
             rf"$\chi^2$/ndf={chi3:.0f}/{ndf3}  (1 factor: {K1[0]:.3f}, "
             rf"$\chi^2$/ndf={chi1:.0f}/{ndf1})", fontsize=16, y=0.975)
out = common.fig_path("sf", "cr_zmumu_dy_jetbin_fit.png", key=MODEL_KEY)
plt.savefig(out, dpi=90, bbox_inches="tight")
print(f"Saved {out}")

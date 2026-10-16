"""tau -> 3-prong background: candidate masses in OS data, signal, DY, TTto2L2Nu and W+jets across stage-1 (ssdata) cuts.

    python3 draw_tau3pi_masses_by_cut.py [--phi-win=0.020]

--phi-win=X adds |phi_mass - M_phi| < X GeV (post-fit branch) on top of the analysis selection
and appends _phi<MeV> to the file names; without it the selection is the producer's
M(KK) < M_phi + 50 MeV only.

Stripped-down copy of draw_tt_m2pi.py: only the SS-data-trained stage-1 DNN (ssdata), no
OS-trained comparison. One figure per variable; rows = samples, columns = no cut / eps_sig
90/70/50%. Absolute yields, y not shared (feedback_plot_style).

Variables (track 4-vectors from pt/eta/phi, mass hypothesis in brackets; Kp and pi are opposite
in charge to the muon, Km same as the muon, so (Km,pi) is the opposite-sign pair):
  Ds_mass   M(KKpi) from the vertex fit (branch)
  Kppi      M(K+[K] pi[pi])            same-sign pair
  Kmpi      M(K-[K] pi[pi])            opposite-sign pair (K*0)
  KpKm      M(K+[K] K-[K])             track level, not the post-fit phi_mass branch
  pipi      M(pi[pi] K-[pi])           opposite-sign pair, K- misID'd as pi (rho0 of tau -> a1 -> rho pi)
  pipipi    M(K+[pi] K-[pi] pi[pi])    all three as pions (a1, <= m_tau)

OS data: the Ds-mass SR (|Ds_mass - M_Ds| < 40 MeV) is removed per event, so it is missing from
EVERY variable, not only Ds_mass -- the data rows of the derived masses are sideband-only while
signal sits almost entirely in the SR. MC (signal, DY, TTto2L2Nu = tt_decay 0, all GenTag modes,
nominal+ext1, OS W+jets 0J+1J+2J) is drawn over the full range. DY/TT/W+jets carry negative weights: errors sqrt(sum w^2),
N_eff and negative-bin counts are printed.
Output: <ssdata model dir>/figures/os_ss/tau3pi_<var>_by_cut.png.
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

KEY = "ssdata"
PHI_WIN = next((float(x.split("=", 1)[1]) for x in sys.argv[1:] if x.startswith("--phi-win=")), None)
TAG = "" if PHI_WIN is None else f"_phi{1e3 * PHI_WIN:.0f}"
CUTS = [(None, "no cut"), (0.90, r"$\varepsilon_{sig}$=90%"), (0.70, r"$\varepsilon_{sig}$=70%"),
        (0.50, r"$\varepsilon_{sig}$=50%")]
M_PI, M_K = 0.13957039, 0.493677
M_DS, M_DP, SR_HALF, M_PHI = 1.96835, 1.8696, 0.040, 1.019460
M_KST, M_RHO, M_TAU = 0.89555, 0.77526, 1.77686
# (tag, xlabel, xrange, bin width [MeV] per cut column, reference lines [(x, label)])
VARS = [
    ("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", (M_DS - 0.3, M_DS + 0.3), (10, 10, 10, 10),
     [(M_DP, r"D$^\pm$"), (M_DS, r"D$_s^\pm$")]),
    ("Kppi", r"$M(K^+\pi)$ (same sign) [GeV]", (0.6, 1.8), (20, 20, 30, 40), []),
    ("Kmpi", r"$M(K^-\pi)$ (opposite sign) [GeV]", (0.6, 1.8), (20, 20, 30, 40), [(M_KST, r"K$^{*0}$")]),
    ("KpKm", r"$M(K^+K^-)$ [GeV]", (0.984, 1.074), (2, 2, 3, 4), [(M_PHI, r"$\phi$")]),
    ("pipi", r"$M(\pi\pi)$ of $\pi K^-$, $K^-\to\pi$ [GeV]", (0.2, 1.64), (20, 20, 30, 40), [(M_RHO, r"$\rho^0$")]),
    ("pipipi", r"$M(\pi\pi\pi)$, $K^\pm\to\pi$ [GeV]", (0.8, 2.3), (20, 20, 30, 40), [(M_TAU, r"m$_\tau$")]),
]
TRK = [f"{p}_{v}" for p in ("Kp", "Km", "pi") for v in ("pt", "eta", "phi")]
WANT = list(dict.fromkeys(["Ds_mass", "phi_mass", "score_dsmu", *TRK, *common.SELECTION]))
# (label, file, extra branches, row mask fn, is data)
SAMPLES = [
    ("OS data (Ds SR excluded)", "Data_OS_merged", [], lambda a: np.abs(a["Ds_mass"] - M_DS) >= SR_HALF, True),
    ("signal (DsW+DsstW)", "Signal_merged", ["weight"], None, False),
    (r"DY$\to\ell\ell$ (0J+1J+2J)", "DY_merged", ["weight"], None, False),
    (r"$t\bar{t}\to 2\ell 2\nu$ (all modes)", "TT_merged", ["weight", "tt_decay"], lambda a: a["tt_decay"] == 0, False),
    (r"W$\to\ell\nu$+jets (0J+1J+2J)", "WtoLNu2Jets_merged", ["weight"], None, False),
]
COLORS = ["black", "#c0392b", "#3a7dc9", "#8e6bb8", "seagreen"]


def mass(a, tracks):
    """Invariant mass from pt/eta/phi; tracks = [(prefix, mass hypothesis)]."""
    E = px = py = pz = 0.0
    for p, m0 in tracks:
        pt, eta, phi = a[f"{p}_pt"], a[f"{p}_eta"], a[f"{p}_phi"]
        E = E + np.sqrt((pt * np.cosh(eta)) ** 2 + m0 ** 2)
        px, py, pz = px + pt * np.cos(phi), py + pt * np.sin(phi), pz + pt * np.sinh(eta)
    return np.sqrt(np.maximum(E ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))


THR = {e: common.stage1_threshold({"key": KEY, "tool": "DNN", "eff_sig": e})[0] for e, _l in CUTS if e}
pred = f"{common.out_dir_for(KEY, 'DNN')}/Prediction"
DATA = []   # per sample: dict(var -> values), w, score
for lab, fname, extra, rowfn, is_data in SAMPLES:   # serial: parallel big reads got OOM-killed before
    with uproot.open(f"{pred}/{fname}.root")[common.TREE] as t:
        a = t.arrays(WANT + extra, library="np")
    common.check_charge_config(a, False, fname)
    k = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if PHI_WIN is not None:
        k &= np.abs(a["phi_mass"] - M_PHI) < PHI_WIN
    if rowfn is not None:
        k &= rowfn(a)
    a = {b: v[k] for b, v in a.items()}
    v = {"Ds_mass": a["Ds_mass"],
         "Kppi": mass(a, [("Kp", M_K), ("pi", M_PI)]),
         "Kmpi": mass(a, [("Km", M_K), ("pi", M_PI)]),
         "KpKm": mass(a, [("Kp", M_K), ("Km", M_K)]),
         "pipi": mass(a, [("pi", M_PI), ("Km", M_PI)]),
         "pipipi": mass(a, [("Kp", M_PI), ("Km", M_PI), ("pi", M_PI)])}
    w = np.ones(k.sum()) if is_data else a["weight"].astype(float)
    d = v["KpKm"] - a["phi_mass"]
    print(f"{fname}: {len(k):,} rows -> {k.sum():,} selected; track-level M(KK) - phi_mass: "
          f"mean {1e3 * d.mean():+.2f} MeV, rms {1e3 * d.std():.2f} MeV")
    for tag, _xl, (lo, hi), *_ in VARS:
        out = (v[tag] < lo) | (v[tag] >= hi)
        print(f"    {tag:8s} outside [{lo:.3f},{hi:.3f}) (not drawn): {w[out].sum() / w.sum():.3%}")
    DATA.append((v, w, a["score_dsmu"]))
    del a

plt.style.use(hep.style.CMS)
for tag, xlab, (lo, hi), bmev, refs in VARS:
    print(f"\n### {tag}")
    nr, nc = len(SAMPLES), len(CUTS)
    fig, axes = plt.subplots(nr, nc, figsize=(8.5 * nc, 7.0 * nr), squeeze=False)
    for r, ((slab, fname, _x, _f, is_data), (v, w, sc), col) in enumerate(zip(SAMPLES, DATA, COLORS)):
        for c, (e, clab) in enumerate(CUTS):
            ax = axes[r, c]
            bins = np.linspace(lo, hi, int(round((hi - lo) * 1000 / bmev[c])) + 1)
            cut = np.ones(len(w), bool) if e is None else sc > THR[e]
            x, ww = v[tag][cut], w[cut]
            hh = np.histogram(x, bins=bins, weights=ww)[0]
            err = np.sqrt(np.histogram(x, bins=bins, weights=ww ** 2)[0])
            if is_data and tag == "Ds_mass":   # blinded bins are empty, not measured zeros
                bc = 0.5 * (bins[1:] + bins[:-1])
                blind = np.abs(bc - M_DS) < SR_HALF
                hh, err = np.where(blind, np.nan, hh), np.where(blind, np.nan, err)
            if is_data:
                hep.histplot(hh, bins, ax=ax, histtype="errorbar", yerr=err, color=col, marker="s",
                             markersize=4, elinewidth=1.1, label=f"N={np.nansum(hh):,.0f}")
            else:
                hep.histplot(hh, bins, ax=ax, histtype="fill", alpha=0.5, color=col,
                             label=f"N={hh.sum():,.3g}")
                bc = 0.5 * (bins[1:] + bins[:-1])
                ax.bar(bc, 2 * err, width=bins[1] - bins[0], bottom=hh - err, color="none", hatch="////",
                       edgecolor="0.35", linewidth=0.0, zorder=3)
            if tag == "Ds_mass":
                ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
            ax.axhline(0, color="black", lw=0.8, alpha=0.5)
            ax.set_xlim(lo, hi)
            top = np.nanmax(hh + err)
            ax.set_ylim(min(0.0, np.nanmin(hh - err) * 1.1), 1.6 * top if top > 0 else 1.0)
            for xr, rl in refs:
                ax.axvline(xr, color="green", ls=":", lw=1.4)
                ax.text(xr, ax.get_ylim()[1] * 0.97, rl, color="green", ha="center", va="top", fontsize=18)
            ax.set_ylabel(f"Events / {bmev[c]} MeV", fontsize=22)
            ax.set_xlabel(xlab, fontsize=20)
            ax.set_title(f"{slab}, {clab}", fontsize=19, pad=30 if r == 0 else 8)   # row 0 clears the CMS label
            ax.legend(fontsize=17, loc="upper left", bbox_to_anchor=(0.0, 0.90))   # below the ref-line labels
            sw, sw2 = ww.sum(), (ww ** 2).sum()
            print(f"  {fname:16s} {'no cut' if e is None else f'eps{100 * e:.0f}':7s}: rows {cut.sum():9,}  "
                  f"sumw {sw:11,.1f}  N_eff {sw ** 2 / sw2 if sw2 else 0:10,.1f}  neg bins {(hh < 0).sum()}")
    hep.cms.label(ax=axes[0, 0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=16)
    fig.suptitle(rf"Stage-1 ({KEY} DNN) cut scan: {xlab.replace(' [GeV]', '')}"
                 + ("" if PHI_WIN is None else rf", $|M(K^+K^-)-m_\phi|<{1e3 * PHI_WIN:.0f}$ MeV"), fontsize=26, y=1.0)
    fig.tight_layout()
    out = common.fig_path("os_ss", f"tau3pi_{tag}_by_cut{TAG}.png", key=KEY)
    plt.savefig(out, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

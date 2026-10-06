"""DY: M(pi pi pi) of the Ds candidate (all three tracks given the pion mass), by best-muon origin.

    python3 draw_dy_m3pi.py [KEY_OS]      (default osdata; compared with ssdata)

A tau -> 3-prong decay feeding all three tracks has M(pipipi) <= m_tau = 1.777 GeV and, by charge
(two tracks opposite the muon, one same), can only form an OS candidate. Same layout as
dy_ds_mass_mu_origin_ssdata_vs_osdata.png: SS-trained vs OS-trained cut overlaid on the SAME rows,
sub-panel KEY_OS / ssdata per bin with the overlap-aware error; one row per muon origin, plus a signal row (DsW+DsstW, all rows) for reference.
Muon origin comes from the cache written by draw_dy_mu_origin.py (run that first).
Output: os_ss/dy_m3pi_mu_origin_ssdata_vs_<KEY_OS>.png (shared dir).
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
BIN_MEV = {0.90: 20, 0.70: 30, 0.50: 40}
XRANGE = (0.8, 2.3)
M_PI, M_TAU = 0.13957039, 1.77686
CACHE = f"{common.MVA_ROOT}/aux/DY_merged_mu_genPartFlav.npz"
ORIGINS = [("muon from tau (flav 15)", [15]), ("prompt muon (flav 1)", [1])]
TRK = [f"{p}_{v}" for p in ("Kp", "Km", "pi") for v in ("pt", "eta", "phi")]
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "luminosityBlock", "event", "mu_charge",
                           "pi_charge", *TRK, *common.SELECTION]))


def h(x, w, bins):
    return np.histogram(x, bins=bins, weights=w)[0], np.histogram(x, bins=bins, weights=w ** 2)[0]


def load(fname):
    """(M(pipipi), weight, pass_ssdata{eps}, pass_KEY_OS{eps}, row-mask a0-aligned arrays) after selection."""
    arr = {}
    for k in (BASE, KEY_OS):
        with uproot.open(f"{common.out_dir_for(k, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
            arr[k] = t.arrays(WANT, library="np")
    a0 = arr[BASE]
    assert np.array_equal(arr[KEY_OS]["Ds_mass"], a0["Ds_mass"]), f"{fname}: row order differs"
    common.check_charge_config(a0, False, fname)
    E = px = py = pz = 0.0
    for p in ("Kp", "Km", "pi"):
        pt, eta, phi = a0[f"{p}_pt"], a0[f"{p}_eta"], a0[f"{p}_phi"]
        E = E + np.sqrt((pt * np.cosh(eta)) ** 2 + M_PI ** 2)
        px, py, pz = px + pt * np.cos(phi), py + pt * np.sin(phi), pz + pt * np.sinh(eta)
    m = np.sqrt(np.maximum(E ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
    sel = common.selection_mask(a0)
    pb = {e: arr[BASE]["score_dsmu"][sel] > thr[BASE][e] for e in EFFS}
    po = {e: arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e] for e in EFFS}
    return a0, sel, m[sel], a0["weight"][sel].astype(float), pb, po


thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in (BASE, KEY_OS)}
dy0, dysel, dy_m, dy_w, dy_pb, dy_po = load("DY_merged")
c = np.load(CACHE)
rowkey = (dy0["luminosityBlock"].astype(np.int64) << 32) | dy0["event"].astype(np.int64)
assert np.array_equal(c["rowkey"], rowkey), f"{CACHE} does not match DY_merged rows"
fl = c["flav"][dysel]
_, _, sg_m, sg_w, sg_pb, sg_po = load("Signal_merged")
# rows: (label, M, w, pass_b, pass_o, row mask, is signal)
ROWS = [(f"DY, {o}", dy_m, dy_w, dy_pb, dy_po, np.isin(fl, v), False) for o, v in ORIGINS]
ROWS.append(("signal (DsW+DsstW)", sg_m, sg_w, sg_pb, sg_po, np.ones(len(sg_w), bool), True))
for lab, m, w, *_ in ROWS[::2]:
    out_rng = (m < XRANGE[0]) | (m >= XRANGE[1])
    print(f"{lab.split(',')[0]}: M(pipipi) outside {XRANGE} (not drawn): {w[out_rng].sum() / w.sum():.3%}")

plt.style.use(hep.style.CMS)
nrow = len(ROWS)
fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25, top=1 - 0.9 / (8.0 * nrow))
for row, (olab, m3pi, w, pass_b, pass_o, g, is_sig) in enumerate(ROWS):
    for col, e in enumerate(EFFS):
        bins = np.linspace(*XRANGE, int(round((XRANGE[1] - XRANGE[0]) * 1000 / BIN_MEV[e])) + 1)
        bc = 0.5 * (bins[1:] + bins[:-1])
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1],
                                                 hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        pb, po = g & pass_b[e], g & pass_o[e]
        hb, vb = h(m3pi[pb], w[pb], bins)
        ho, vo = h(m3pi[po], w[po], bins)
        hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)), color="steelblue",
                     marker="o", markersize=5, elinewidth=1.3, label=f"{BASE} (SS DNN) cut  N={hb.sum():,.0f}")
        hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)), color="black",
                     marker="s", markersize=5, elinewidth=1.3, label=f"{KEY_OS} (OS DNN) cut  N={ho.sum():,.0f}")
        ax.axvline(M_TAU, color="green", ls=":", lw=1.5)
        ax.set_xlim(bins[0], bins[-1])
        top = max(hb.max(), ho.max())
        ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * top if top > 0 else 1.0)
        ax.text(M_TAU + 0.01, ax.get_ylim()[1] * 0.75, r"m$_\tau$", color="green", ha="left", va="top",
                fontsize=13)
        ax.set_ylabel(f"Events / {BIN_MEV[e]} MeV")
        ax.set_title(rf"{olab}, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
        ax.legend(fontsize=11, loc="upper left")
        ax.tick_params(labelbottom=False)
        # per-bin ratio, overlap-aware: (C+A)/(C+B) with C both pass, A only OS, B only SS
        cc, vc = h(m3pi[pb & po], w[pb & po], bins)
        a_, va = h(m3pi[po & ~pb], w[po & ~pb], bins)
        b, vbb = h(m3pi[pb & ~po], w[pb & ~po], bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (cc + a_) / (cc + b)
            er = np.sqrt(vc * ((b - a_) / (cc + b) ** 2) ** 2 + va / (cc + b) ** 2
                         + vbb * ((cc + a_) / (cc + b) ** 2) ** 2)
        rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvline(M_TAU, color="green", ls=":", lw=1.2)
        rax.set_ylim(0.0, 1.99)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(r"$M(\pi\pi\pi)$ [GeV]")
        rax.set_ylabel("OS / SS")
        below = bc < M_TAU
        with np.errstate(divide="ignore", invalid="ignore"):
            print(f"  {olab:26s} eps{100 * e:.0f}: SS-DNN {hb.sum():8,.1f} ({hb[below].sum() / hb.sum():5.1%} below m_tau)"
                  f"  OS-DNN {ho.sum():8,.1f} ({ho[below].sum() / ho.sum():5.1%} below m_tau)")
hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
fig.suptitle(rf"DY$\to\ell\ell$ (0J+1J+2J) by best-muon origin, and signal: M($\pi\pi\pi$) after SS-trained ({BASE}) "
             f"vs OS-trained ({KEY_OS}) stage-1 cut", fontsize=18, y=0.995)
out = common.fig_path("os_ss", f"dy_m3pi_mu_origin_{BASE}_vs_{KEY_OS}.png")
plt.savefig(out, dpi=90, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

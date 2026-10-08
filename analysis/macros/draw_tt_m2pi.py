"""ttbar: M(pi+ pi-) of the two opposite-charge track pairs of the Ds candidate (pion mass), mode 0 by decay.

    python3 draw_tt_m2pi.py [KEY_OS]      (default osdata; compared with ssdata)

Copied from draw_dy_m2pi.py; rows are ttbar OS GenTag mode 0 split by tt_decay (2L2Nu: other
W -> tau -> 3-prong possible; LNu2Q: not, control) instead of the DY muon origin.

Looks for the rho0(770) of tau -> a1 nu -> rho0 pi nu. The opposite-charge pairs are (K+, K-) -- the
phi pair, so its K-mass M(KK) < M_phi+0.05 cut squeezes it to low M(pipi) -- and (pi, K-); (pi, K+)
is same-sign. One figure per pair. Same layout as draw_dy_m3pi.py: SS-trained vs
OS-trained cut overlaid on the SAME rows, sub-panel KEY_OS / ssdata per bin with the overlap-aware
error; rows ttbar 2L2Nu OS mode 0, ttbar LNu2Q OS mode 0, signal (DsW+DsstW, all rows).
Output: os_ss/tt_m2pi_{KpKm,piKm}_mode0_ssdata_vs_<KEY_OS>.png (shared dir).
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
ARGS = [x for x in sys.argv[1:] if not x.startswith("--")]
KEY_OS = ARGS[0] if ARGS else "osdata"
# --decay=<2L2Nu|LNu2Q>: one decay channel, one row per GenTag mode group (all modes first) instead
DECAY = next((x.split("=", 1)[1] for x in sys.argv[1:] if x.startswith("--decay=")), None)
MODE_GROUPS = [("all modes", None), ("mode 0: comb.", [0]), ("mode 7: multi-charm", [7]),
               ("mode 8: partial", [8]), (r"mode 1+2: true D$_s$", [1, 2]),
               (r"mode 3+4: D$\to$KK$\pi$", [3, 4]), (r"mode 5: D$\to$K$\pi\pi$ misID", [5]),
               ("mode 6+9+10: other", [6, 9, 10])]
EFFS = (0.90, 0.70, 0.50)
BIN_MEV = {0.90: 20, 0.70: 30, 0.50: 40}
XRANGE = (0.2, 1.64)
M_PI, M_RHO = 0.13957039, 0.77526
PAIRS = {"KpKm": ("Kp", "Km", r"$K^+K^-$"), "piKm": ("pi", "Km", r"$\pi K^-$")}
DECAYS = [("2L2Nu", 0), ("LNu2Q", 1)]
assert DECAY in (None, *dict(DECAYS)), f"--decay must be one of {[d for d, _ in DECAYS]}"
TRK = [f"{p}_{v}" for p in ("Kp", "Km", "pi") for v in ("pt", "eta", "phi")]
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "mu_charge",
                           "pi_charge", *TRK, *common.SELECTION]))


def h(x, w, bins):
    return np.histogram(x, bins=bins, weights=w)[0], np.histogram(x, bins=bins, weights=w ** 2)[0]


def mass_pi(a, pair):
    E = px = py = pz = 0.0
    for p in pair:
        pt, eta, phi = a[f"{p}_pt"], a[f"{p}_eta"], a[f"{p}_phi"]
        E = E + np.sqrt((pt * np.cosh(eta)) ** 2 + M_PI ** 2)
        px, py, pz = px + pt * np.cos(phi), py + pt * np.sin(phi), pz + pt * np.sinh(eta)
    return np.sqrt(np.maximum(E ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))


def load(fname, extra=()):
    """Selected rows of BASE's file plus pass masks for both models."""
    arr = {}
    for k in (BASE, KEY_OS):
        with uproot.open(f"{common.out_dir_for(k, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
            arr[k] = t.arrays(WANT + list(extra), library="np")
    a0 = arr[BASE]
    assert np.array_equal(arr[KEY_OS]["Ds_mass"], a0["Ds_mass"]), f"{fname}: row order differs"
    common.check_charge_config(a0, False, fname)
    sel = common.selection_mask(a0)
    pb = {e: arr[BASE]["score_dsmu"][sel] > thr[BASE][e] for e in EFFS}
    po = {e: arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e] for e in EFFS}
    return a0, sel, a0["weight"][sel].astype(float), pb, po


thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in (BASE, KEY_OS)}
tt0, ttsel, tt_w, tt_pb, tt_po = load("TT_merged", ("cand_gen_mode", "tt_decay"))
md, dec = tt0["cand_gen_mode"][ttsel], tt0["tt_decay"][ttsel]
sg0, sgsel, sg_w, sg_pb, sg_po = load("Signal_merged")

plt.style.use(hep.style.CMS)
for ptag, (p1, p2, plab) in PAIRS.items():
    tt_m = mass_pi(tt0, (p1, p2))[ttsel]
    sg_m = mass_pi(sg0, (p1, p2))[sgsel]
    # rows: (label, M, w, pass_b, pass_o, row mask)
    if DECAY:
        dmask = dec == dict(DECAYS)[DECAY]
        ROWS = [(rf"$t\bar{{t}}$ {DECAY} OS, {gl}", tt_m, tt_w, tt_pb, tt_po,
                 dmask & (True if mv is None else np.isin(md, mv))) for gl, mv in MODE_GROUPS]
    else:
        ROWS = [(rf"$t\bar{{t}}$ {d} OS, mode 0", tt_m, tt_w, tt_pb, tt_po, (md == 0) & (dec == v)) for d, v in DECAYS]
    ROWS.append(("signal (DsW+DsstW)", sg_m, sg_w, sg_pb, sg_po, np.ones(len(sg_w), bool)))
    print(f"\n### pair {ptag}")
    for lab, m, w, *_, g in ROWS:
        m, w = m[g], w[g]
        out_rng = (m < XRANGE[0]) | (m >= XRANGE[1])
        print(f"  {lab}: M outside {XRANGE} (not drawn): {w[out_rng].sum() / w.sum():.3%}")

    nrow = len(ROWS)
    fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
    outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25,
                              top=1 - 0.9 / (8.0 * nrow))
    for row, (olab, m2, w, pass_b, pass_o, g) in enumerate(ROWS):
        for col, e in enumerate(EFFS):
            bins = np.linspace(*XRANGE, int(round((XRANGE[1] - XRANGE[0]) * 1000 / BIN_MEV[e])) + 1)
            bc = 0.5 * (bins[1:] + bins[:-1])
            inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col],
                                                     height_ratios=[3, 1], hspace=0.05)
            ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            pb, po = g & pass_b[e], g & pass_o[e]
            hb, vb = h(m2[pb], w[pb], bins)
            ho, vo = h(m2[po], w[po], bins)
            hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)), color="steelblue",
                         marker="o", markersize=5, elinewidth=1.3,
                         label=f"{BASE} (SS DNN) cut  N={hb.sum():,.3g}")
            hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)), color="black",
                         marker="s", markersize=5, elinewidth=1.3,
                         label=f"{KEY_OS} (OS DNN) cut  N={ho.sum():,.3g}")
            ax.axvline(M_RHO, color="green", ls=":", lw=1.5)
            ax.set_xlim(bins[0], bins[-1])
            top = max(hb.max(), ho.max())
            ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * top if top > 0 else 1.0)
            ax.text(M_RHO + 0.01, ax.get_ylim()[1] * 0.75, r"m$_\rho$", color="green", ha="left", va="top",
                    fontsize=13)
            ax.set_ylabel(f"Events / {BIN_MEV[e]} MeV")
            ax.set_title(rf"{olab}, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
            ax.legend(fontsize=11, loc="upper right")
            ax.tick_params(labelbottom=False)
            # per-bin ratio, overlap-aware: (C+A)/(C+B) with C both pass, A only OS, B only SS
            cc, vc = h(m2[pb & po], w[pb & po], bins)
            a_, va = h(m2[po & ~pb], w[po & ~pb], bins)
            b, vbb = h(m2[pb & ~po], w[pb & ~po], bins)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = (cc + a_) / (cc + b)
                er = np.sqrt(vc * ((b - a_) / (cc + b) ** 2) ** 2 + va / (cc + b) ** 2
                             + vbb * ((cc + a_) / (cc + b) ** 2) ** 2)
            rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
            rax.axhline(1.0, color="red", lw=1.0)
            rax.axvline(M_RHO, color="green", ls=":", lw=1.2)
            rax.set_ylim(0.0, 1.99)
            rax.set_yticks([0.5, 1.0, 1.5])
            rax.set_xlim(bins[0], bins[-1])
            rax.set_xlabel(rf"$M(\pi^+\pi^-)$ of {plab} [GeV]")
            rax.set_ylabel("OS / SS")
            rho = np.abs(bc - M_RHO) < 0.15
            with np.errstate(divide="ignore", invalid="ignore"):
                print(f"  {olab:34s} eps{100 * e:.0f}: SS-DNN {hb.sum():8,.1f} "
                      f"({hb[rho].sum() / hb.sum():5.1%} in m_rho+-150)  OS-DNN {ho.sum():8,.1f} "
                      f"({ho[rho].sum() / ho.sum():5.1%} in m_rho+-150)")
    hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
    fig.suptitle((rf"$t\bar{{t}}$ {DECAY} OS by GenTag mode" if DECAY else r"$t\bar{t}$ OS mode 0 by decay channel")
                 + rf", and signal: M($\pi^+\pi^-$) of {plab} "
                 f"after SS-trained ({BASE}) vs OS-trained ({KEY_OS}) stage-1 cut", fontsize=18, y=0.995)
    out = common.fig_path("os_ss", f"tt_m2pi_{ptag}_{DECAY + '_by_mode' if DECAY else 'mode0'}_{BASE}_vs_{KEY_OS}.png")
    plt.savefig(out, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

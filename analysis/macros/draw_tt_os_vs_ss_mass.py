"""ttbar OS vs ttbar SS Ds_mass, overlaid, after the ssdata stage-1 cut.

    python3 draw_tt_os_vs_ss_mass.py [KEY] [--modes] [--decay=2L2Nu|LNu2Q]      (default ssdata, both decays)

Columns: no cut, then eps_sig = 90/70/50% (KEY's stage1_threshold). Each panel overlays TT_merged (OS) and
TT_SS_merged (SS), both absolute weighted yields (same norm, Sum w shared by OS/SS), all GenTag
modes (--modes: one row per GenTag mode group instead). Sub-panel OS/SS per bin, errors treated as independent (different files; the two can share
input events, so this is approximate). MC, so the Ds window is not blinded; grey band = SR only.
--decay keeps one tt_decay (0 = TTto2L2Nu, 1 = TTtoLNu2Q, nominal + ext1 alike).
Output: os_ss/ds_mass_tt_os_vs_ss_<KEY>[_<decay>][_by_mode].png (shared dir).
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

ARGS = [x for x in sys.argv[1:] if not x.startswith("--")]
KEY = ARGS[0] if ARGS else "ssdata"
BY_MODE = "--modes" in sys.argv
DECAYS = {"2L2Nu": 0, "LNu2Q": 1}
DECAY = next((x.split("=", 1)[1] for x in sys.argv[1:] if x.startswith("--decay=")), None)
assert DECAY in (None, *DECAYS), f"--decay must be one of {list(DECAYS)}"
# (row label, cand_gen_mode values or None = all)
GROUPS = ([("all modes", None), ("mode 0: comb.", [0]), ("mode 7: multi-charm", [7]), ("mode 8: partial", [8]),
           (r"mode 1+2: true D$_s$", [1, 2]), (r"mode 3+4: D$\to$KK$\pi$", [3, 4]),
           (r"mode 5: D$\to$K$\pi\pi$ misID", [5]), ("mode 6+9+10: other", [6, 9, 10])]
          if BY_MODE else [("all modes", None)])
EFFS = (None, 0.90, 0.70, 0.50)   # None = no cut
BIN_MEV = {None: 10, 0.90: 20, 0.70: 30, 0.50: 40}
M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.04
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "cand_gen_mode", "tt_decay", "mu_charge", "pi_charge",
                           *common.SELECTION]))

thr = {e: common.stage1_threshold({"key": KEY, "tool": "DNN", "eff_sig": e})[0] for e in EFFS if e}
thr[None] = -np.inf
plt.style.use(hep.style.CMS)


def h(x, w, bins):
    return np.histogram(x, bins=bins, weights=w)[0], np.histogram(x, bins=bins, weights=w ** 2)[0]


data = {}
for lab, fname, is_ss in (("OS", "TT_merged", False), ("SS", "TT_SS_merged", True)):
    with uproot.open(f"{common.out_dir_for(KEY, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
        a = t.arrays(WANT, library="np")
    common.check_charge_config(a, is_ss, fname)
    sel = common.selection_mask(a)
    if DECAY:
        sel &= a["tt_decay"] == DECAYS[DECAY]
    data[lab] = (a["Ds_mass"][sel], a["weight"][sel].astype(float), a["score_dsmu"][sel], a["cand_gen_mode"][sel])
    print(f"{lab}: {fname}, {sel.sum():,} rows, no-cut yield {data[lab][1].sum():,.1f}")
print(f"no cut: OS/SS total {data['OS'][1].sum() / data['SS'][1].sum():.3f}")

nrow = len(GROUPS)
fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25, top=1 - 0.9 / (8.0 * nrow))
for row, (glab, modes) in enumerate(GROUPS):
    print(f"--- {glab}")
    for col, e in enumerate(EFFS):
        bins = np.linspace(M_DS - 0.3, M_DS + 0.3, int(round(600 / BIN_MEV[e])) + 1)
        bc = 0.5 * (bins[1:] + bins[:-1])
        bw = 1000 * (bins[1] - bins[0])
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        hist = {}
        for lab, color, marker in (("SS", "steelblue", "o"), ("OS", "black", "s")):
            m, w, s, md = data[lab]
            p = (s > thr[e]) & (True if modes is None else np.isin(md, modes))
            hist[lab] = h(m[p], w[p], bins)
            hep.histplot(hist[lab][0], bins, ax=ax, histtype="errorbar", yerr=np.sqrt(hist[lab][1]),
                         color=color, marker=marker, markersize=5, elinewidth=1.3,
                         label=rf"$t\bar{{t}}$ {lab}  N={hist[lab][0].sum():,.0f}")
        (ho, vo), (hs, vs) = hist["OS"], hist["SS"]
        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0.0, 1.4 * max(ho.max(), hs.max(), 1e-9))
        ax.set_ylabel(f"Events / {bw:.0f} MeV")
        ax.set_title(f"{glab}, " + ("no cut" if e is None else rf"{KEY} cut, $\varepsilon_{{sig}}$={100 * e:.0f}%"), fontsize=15)
        ax.legend(fontsize=11, loc="upper right")
        ax.tick_params(labelbottom=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = ho / hs
            er = np.abs(r) * np.sqrt(vo / ho ** 2 + vs / hs ** 2)   # bins can go negative (genWeight < 0)
        rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        ok = np.isfinite(r) & np.isfinite(er) & (er < r)   # skip empty / pure-noise bins when scaling
        rax.set_ylim(0.0, max(2.0, 1.15 * np.max(r[ok] + er[ok])) if ok.any() else 2.0)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        rax.set_ylabel("OS / SS")
        lo, sr, hi = bc < 1.830, np.abs(bc - M_DS) < SR_HALF, bc > M_DS + SR_HALF
        with np.errstate(divide="ignore", invalid="ignore"):
            print(f"{'no cut' if e is None else f'eps{100 * e:.0f}':6s}: OS {ho.sum():8,.1f}  SS {hs.sum():8,.1f}  OS/SS total {ho.sum() / hs.sum():.3f}"
                  f"  low {ho[lo].sum() / hs[lo].sum():.3f}  SR {ho[sr].sum() / hs[sr].sum():.3f}"
                  f"  high {ho[hi].sum() / hs[hi].sum():.3f}")
hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
fig.suptitle(rf"$t\bar{{t}}$ ({DECAY or '2L2Nu + LNu2Q'}): Ds_mass, OS vs SS after {KEY} stage-1 cut",
             fontsize=18, y=0.995)
out = common.fig_path("os_ss", f"ds_mass_tt_os_vs_ss_{KEY}{'_' + DECAY if DECAY else ''}{'_by_mode' if BY_MODE else ''}.png")
plt.savefig(out, dpi=100, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

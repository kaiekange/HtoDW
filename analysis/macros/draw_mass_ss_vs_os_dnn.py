"""Ds_mass after the SS-trained (ssdata) vs the OS-trained (osdata) stage-1 cut, overlaid.

    python3 draw_mass_ss_vs_os_dnn.py [KEY_OS] [TAG ...] [--modes]   (default osdata vs ssdata, all tags)

Per sample one figure; columns eps_sig = 90/70/50% (each model at its own stage1_threshold).
Each panel overlays the two cuts on the SAME rows; the sub-panel is KEY_OS / ssdata per bin,
error from the both-pass / only-one splitting (the two selections overlap heavily).
DY, ttbar OS/SS, W+jets OS/SS: all GenTag modes together (--modes: plus one row per mode group,
output name gets _by_mode).
Neither model trained on DY, ttbar or W+jets, so all rows are used (no parity cut). MC, so the Ds
window is not blinded; the grey band only marks the SR.
Data OS / data SS: each is one model's training sample (osdata / ssdata, even events), so BOTH
models are compared on the odd half only, x2. Data bins sit on the blind-window edges
(13/10/5 per side = 20/26/52 MeV); the blinded OS Ds window is left empty.
TT_2L2Nu / TT_LNu2Q / TT_SS_2L2Nu / TT_SS_LNu2Q: one ttbar decay channel (tt_decay).
Output: os_ss/ds_mass_ssdata_vs_<KEY_OS>_{DY,TT,TT_SS,WJets,WJets_SS,Data,Data_SS}.png (shared dir).
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
BY_MODE = "--modes" in sys.argv
MODE_GROUPS = [("mode 0: comb.", [0]), ("mode 7: multi-charm", [7]), ("mode 8: partial", [8]),
               (r"mode 1+2: true D$_s$", [1, 2]), (r"mode 3+4: D$\to$KK$\pi$", [3, 4]),
               (r"mode 5: D$\to$K$\pi\pi$ misID", [5]), ("mode 6+9+10: other", [6, 9, 10])]
EFFS = (0.90, 0.70, 0.50)
BIN_MEV = {0.90: 20, 0.70: 30, 0.50: 40}
M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.04
# tag: (file, title, row groups [(label, modes or None=all)])
SAMPLES = {
    "DY": ("DY_merged", r"DY$\to\ell\ell$ (0J+1J+2J), OS", [("all modes", None)]),
    "TT": ("TT_merged", r"$t\bar{t}$ (2L2Nu + LNu2Q), OS", [("all modes", None)]),
    "TT_SS": ("TT_SS_merged", r"$t\bar{t}$ (2L2Nu + LNu2Q), SS", [("all modes", None)]),
    **{f"TT{ss}_{d}": (f"TT{ss}_merged", rf"$t\bar{{t}}$ {d}, {ss[1:] or 'OS'}", [("all modes", None)])
       for ss in ("", "_SS") for d in ("2L2Nu", "LNu2Q")},
    "WJets": ("WtoLNu2Jets_merged", "W+jets (0J+1J+2J), OS", [("all modes", None)]),
    "WJets_SS": ("WtoLNu2Jets_SS_merged", "W+jets (0J+1J+2J), SS", [("all modes", None)]),
    "Data": ("Data_OS_merged", "data OS (odd half x2)", [("data OS", None)]),
    "Data_SS": ("Data_SS_merged", "data SS (odd half x2)", [("data SS", None)]),
}
TT_DECAY = {"2L2Nu": 0, "LNu2Q": 1}   # tt_decay tag; TT*_<decay> keeps one channel
NSIDE_DATA = {0.90: 13, 0.70: 10, 0.50: 5}   # bins per 0.26 GeV side of the blind window
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "event", "mu_charge", "pi_charge",
                           *common.SELECTION]))

thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in (BASE, KEY_OS)}
plt.style.use(hep.style.CMS)


def h(x, w, bins):
    return np.histogram(x, bins=bins, weights=w)[0], np.histogram(x, bins=bins, weights=w ** 2)[0]


TAGS = ARGS[1:] or list(SAMPLES)
for tag in TAGS:
    fname, title, groups = SAMPLES[tag]
    is_data = fname.startswith("Data")
    if BY_MODE and not is_data:
        groups = groups + MODE_GROUPS
    arr = {}
    for k in (BASE, KEY_OS):
        with uproot.open(f"{common.out_dir_for(k, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
            arr[k] = t.arrays(WANT + ([] if is_data else ["cand_gen_mode"])
                              + (["tt_decay"] if fname.startswith("TT") else []), library="np")
    a0 = arr[BASE]
    assert np.array_equal(arr[KEY_OS]["Ds_mass"], a0["Ds_mass"]), f"{fname}: row order differs"
    common.check_charge_config(a0, fname.endswith("_SS_merged"), fname)
    sel = common.selection_mask(a0)
    if is_data:   # odd half, x2: the even half trained one of the two models
        sel &= (a0["event"] % 2) == 1
    if tag.rsplit("_", 1)[-1] in TT_DECAY:
        sel &= a0["tt_decay"] == TT_DECAY[tag.rsplit("_", 1)[-1]]
    mass, w = a0["Ds_mass"][sel], a0["weight"][sel].astype(float) * (2.0 if is_data else 1.0)
    mode = None if is_data else a0["cand_gen_mode"][sel]
    pass_b = {e: arr[BASE]["score_dsmu"][sel] > thr[BASE][e] for e in EFFS}
    pass_o = {e: arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e] for e in EFFS}

    nrow = len(groups)
    fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
    outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25,
                              top=1 - 0.9 / (8.0 * nrow))
    print(f"\n### {tag} ({fname}, {sel.sum():,} rows)")
    for row, (glab, modes) in enumerate(groups):
        g = np.ones(len(mass), bool) if modes is None else np.isin(mode, modes)
        for col, e in enumerate(EFFS):
            if is_data:
                n = NSIDE_DATA[e]
                bins = np.concatenate([np.linspace(M_DS - 0.3, M_DS - SR_HALF, n + 1),
                                       np.linspace(M_DS + SR_HALF, M_DS + 0.3, n + 1)])
            else:
                bins = np.linspace(M_DS - 0.3, M_DS + 0.3, int(round(600 / BIN_MEV[e])) + 1)
            bc = 0.5 * (bins[1:] + bins[:-1])
            blind = (np.abs(bc - M_DS) < SR_HALF) if fname == "Data_OS_merged" else np.zeros(len(bc), bool)
            bw = 1000 * (bins[1] - bins[0])
            inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col],
                                                     height_ratios=[3, 1], hspace=0.05)
            ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            pb, po = g & pass_b[e], g & pass_o[e]
            hb, vb = h(mass[pb], w[pb], bins)
            ho, vo = h(mass[po], w[po], bins)
            for arr_ in (hb, vb, ho, vo):
                arr_[blind] = np.nan
            hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)),
                         color="steelblue", marker="o", markersize=5, elinewidth=1.3,
                         label=f"{BASE} (SS DNN) cut  N={np.nansum(hb):,.0f}")
            hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)),
                         color="black", marker="s", markersize=5, elinewidth=1.3,
                         label=f"{KEY_OS} (OS DNN) cut  N={np.nansum(ho):,.0f}")
            ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
            ax.axvline(M_D, color="green", ls=":", lw=1.5)
            ax.set_xlim(bins[0], bins[-1])
            ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * max(np.nanmax(hb), np.nanmax(ho)))
            ax.set_ylabel(f"Events / {bw:.0f} MeV")
            ax.set_title(rf"{glab}, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
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
            r[blind], er[blind] = np.nan, np.nan
            rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
            rax.axhline(1.0, color="red", lw=1.0)
            rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
            if BY_MODE:   # per-mode ratios can leave [0, 2]; scale on bins whose error < value
                ok = np.isfinite(r) & np.isfinite(er) & (er < np.abs(r))
                rax.set_ylim(0.0, max(1.99, 1.15 * np.max(r[ok] + er[ok])) if ok.any() else 1.99)
            else:
                rax.set_ylim(0.0, 1.99)
                rax.set_yticks([0.5, 1.0, 1.5])
            rax.set_xlim(bins[0], bins[-1])
            rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
            rax.set_ylabel("OS / SS")
            lo, hi = bc < 1.830, bc > M_DS + SR_HALF
            print(f"  {glab:34s} eps{100 * e:.0f}: SS-DNN {np.nansum(hb):9,.1f}  OS-DNN {np.nansum(ho):9,.1f}"
                  f"  OS/SS total {np.nansum(ho) / np.nansum(hb):.3f}  low {ho[lo].sum() / hb[lo].sum():.3f}"
                  f"  high {ho[hi].sum() / hb[hi].sum():.3f}")
    hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=is_data, text="Preliminary", fontsize=11)
    fig.suptitle(f"{title}: Ds_mass after SS-trained ({BASE}) vs OS-trained ({KEY_OS}) stage-1 cut",
                 fontsize=18, y=0.995)
    out = common.fig_path("os_ss", f"ds_mass_{BASE}_vs_{KEY_OS}_{tag}{'_by_mode' if BY_MODE else ''}.png")
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

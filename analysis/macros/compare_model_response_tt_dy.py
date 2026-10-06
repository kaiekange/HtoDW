"""Do the SS-trained and OS-trained stage-1 DNNs respond the same way to OS ttbar->2l2nu and DY?

    python3 compare_model_response_tt_dy.py [key ...]     (default: osdata oswjets; base ssdata)

osdata removes the extra OS preference in data that ssdata keeps (2026-10-06). If that extra
rejection targets a component specific to OS data, a sample without it should be cut the same
by both models. ttbar->2l2nu and DY (Z neutral: no charge correlation between the jets and the
muon) are such samples; W+jets OS and OS data sideband are printed as references.

For each sample and Ds_mass sideband window: eps_key / eps_ssdata at eps_sig = 90/70/50%
(each model's own stage1_threshold). Both models score the SAME rows, so the error splits the
selected weight into both-pass / only-key / only-ssdata pieces (independent) instead of
treating the two efficiencies as independent.
Rows: TT/DY are in no model's training set -> all events. W+jets OS (oswjets trains on its even
half) and OS data (osdata trains on its even half) -> odd events only, for every model.

DY and W+jets OS are also split by GenTag mode group (cand_gen_mode), with each group's share
of the window yield before the cut and after the ssdata cut printed next to the ratio.

Figures (shared dir, multi-model): os_ss/model_response_by_sample.png and
os_ss/model_response_by_mode_<key>.png. A point is drawn only if >= MIN_RAW raw rows pass the
ssdata cut; mode figures clip the y axis to 0-2.5 (negative-weight
noise in small groups); the printed table has every point.
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
KEYS = sys.argv[1:] or ["osdata", "oswjets"]
EFFS = (0.90, 0.70, 0.50)
WINDOWS = {"low <1.830": (0.0, 1.830), "high >2.008": (2.008, 9.0)}
MIN_RAW = 20
# label: (file, extra-row selection on the label branch or None, odd-only, is MC, colour, marker)
SAMPLES = {"TT 2l2nu": ("TT_merged", ("tt_decay", 0), False, True, "#8e44ad", "v"),
           "TT l nu 2q": ("TT_merged", ("tt_decay", 1), False, True, "#c39bd3", "^"),
           "DY": ("DY_merged", None, False, True, "#2e86c1", "o"),
           "W+jets OS": ("WtoLNu2Jets_merged", None, True, True, "seagreen", "D"),
           "data OS": ("Data_OS_merged", None, True, False, "black", "s")}
MODE_SAMPLES = ("DY", "W+jets OS")
MODE_GROUPS = [   # (label, modes, colour) -- colours follow draw_wp_mass.py where they overlap
    ("comb. (0)", [0], "#5a9e78"),
    ("partial (8)", [8], "#5b6fc4"),
    ("multi-charm (7)", [7], "#9aa7e0"),
    (r"true D$_s$ (1+2)", [1, 2], "#c0392b"),
    (r"D$\to$KK$\pi$ (3+4)", [3, 4], "#e07a3c"),
    (r"D$\to$K$\pi\pi$ misID (5)", [5], "#e6c229"),
    ("other (6,9,10)", [6, 9, 10], "grey"),
]
KEYCOLS = ["run", "luminosityBlock", "event", "Ds_mass"]
WANT = list(dict.fromkeys(KEYCOLS + ["weight", "score_dsmu", "mu_charge", "pi_charge",
                                     *common.SELECTION]))


def load(key, fname, extra):
    with uproot.open(f"{common.out_dir_for(key, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
        return t.arrays(list(dict.fromkeys(WANT + extra)), library="np")


def ratio_err(w, pa, pb):
    """eps_a/eps_b on the same rows (same denominator cancels): (C+A)/(C+B)."""
    C, A, B = w[pa & pb], w[pa & ~pb], w[~pa & pb]
    c, a, b = C.sum(), A.sum(), B.sum()
    vc, va, vb = (C ** 2).sum(), (A ** 2).sum(), (B ** 2).sum()
    r = (c + a) / (c + b)
    dr2 = (vc * ((b - a) / (c + b) ** 2) ** 2 + va / (c + b) ** 2
           + vb * ((c + a) / (c + b) ** 2) ** 2)
    return r, np.sqrt(dr2)


thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in [BASE, *KEYS]}
# R[(sample or (sample, group)), key, window] = list over EFFS of (ratio, err, raw rows passing ssdata)
R = {}
for label, (fname, tag, odd, mc, _, _) in SAMPLES.items():
    extra = ([tag[0]] if tag else []) + (["cand_gen_mode"] if label in MODE_SAMPLES else [])
    arr = {k: load(k, fname, extra) for k in [BASE, *KEYS]}
    a0 = arr[BASE]
    for k in KEYS:   # same source file, same row order -> rows must agree exactly
        for c in KEYCOLS:
            assert np.array_equal(arr[k][c], a0[c]), f"{fname}: {k} vs {BASE} differ in {c}"
    common.check_charge_config(a0, False, fname)
    sel = common.selection_mask(a0)
    if tag:
        sel &= a0[tag[0]] == tag[1]
    if odd:
        sel &= (a0["event"] % 2) == 1
    w_all = a0["weight"].astype(float) if mc else np.ones(len(sel))
    subsets = [(label, None)]
    if label in MODE_SAMPLES:
        subsets += [((label, g), modes) for g, modes, _ in MODE_GROUPS]
    print(f"\n### {label}  ({fname}{', odd half' if odd else ''}, {sel.sum():,} rows)")
    for wn, (lo, hi) in WINDOWS.items():
        win = sel & (a0["Ds_mass"] >= lo) & (a0["Ds_mass"] < hi)
        p0_win = {e: arr[BASE]["score_dsmu"] > thr[BASE][e] for e in EFFS}
        tot0 = w_all[win].sum()
        tot_cut = {e: w_all[win & p0_win[e]].sum() for e in EFFS}
        print(f"  --- {wn}")
        for sub, modes in subsets:
            base = win if modes is None else win & np.isin(a0["cand_gen_mode"], modes)
            wb = w_all[base]
            p0 = {e: p0_win[e][base] for e in EFFS}
            name = sub if modes is None else f"  {sub[1]}"
            if modes is None:
                eps0 = "  ".join(f"{wb[p0[e]].sum() / wb.sum():.2e}" for e in EFFS)
                print(f"  {name:28s} eps_{BASE} 90/70/50: {eps0}")
            else:
                share = f"share no cut {wb.sum() / tot0:6.1%} | after {BASE} cut " + " / ".join(
                    f"{wb[p0[e]].sum() / tot_cut[e]:6.1%}" for e in EFFS)
                print(f"  {name:28s} {share}")
            for k in KEYS:
                out, vals = [], []
                for e in EFFS:
                    pk = arr[k]["score_dsmu"][base] > thr[k][e]
                    n_raw = int(p0[e].sum())
                    r, er = ratio_err(wb, pk, p0[e]) if n_raw else (np.nan, np.nan)
                    vals.append((r, er, n_raw))
                    out.append(f"{r:5.2f}±{er:<4.2f}(n={n_raw})")
                R[(sub, k, wn)] = vals
                print(f"  {'':28s}   eps_{k}/eps_{BASE}: " + "  ".join(out))

# ---------------------------------------------------------------- figures
plt.style.use(hep.style.CMS)
XT = np.arange(len(EFFS))
XL = [f"{100 * e:.0f}%" for e in EFFS]


def draw(ax, entries, title, ylabel, ylim=None):
    """entries: (label, [(r, err, n_raw)] over EFFS, colour, marker). ylim fixes the axis."""
    n = len(entries)
    for i, (lab, vals, col, mk) in enumerate(entries):
        x = XT + 0.06 * (i - (n - 1) / 2)
        r = np.array([v[0] if v[2] >= MIN_RAW else np.nan for v in vals])
        e = np.array([v[1] if v[2] >= MIN_RAW else np.nan for v in vals])
        if np.all(np.isnan(r)):
            continue
        ax.errorbar(x, r, yerr=e, color=col, marker=mk, ms=7, lw=1.5, capsize=3, label=lab)
    ax.axhline(1, color="grey", lw=1, ls="--")
    ax.set_xticks(XT)
    ax.set_xticklabels(XL)
    ax.set_xlim(-0.4, len(EFFS) - 0.6)
    ax.set_xlabel(r"stage-1 $\varepsilon_{sig}$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=15)
    lo, hi = ylim if ylim else ax.get_ylim()
    ax.set_ylim(min(lo, 0.0), hi + 0.35 * (hi - lo))
    ax.legend(fontsize=11, loc="upper left", ncol=2)


fig, axes = plt.subplots(len(KEYS), len(WINDOWS), figsize=(9 * len(WINDOWS), 7.5 * len(KEYS)),
                         squeeze=False)
for i, k in enumerate(KEYS):
    for j, wn in enumerate(WINDOWS):
        ent = [(s, R[(s, k, wn)], c, m) for s, (_, _, _, _, c, m) in SAMPLES.items()]
        draw(axes[i][j], ent, f"Ds_mass {wn} GeV",
             rf"$\varepsilon_{{{k}}}\,/\,\varepsilon_{{{BASE}}}$")
hep.cms.label(ax=axes[0][0], rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path("os_ss", "model_response_by_sample.png")
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"\nSaved {out}")

for k in KEYS:
    fig, axes = plt.subplots(len(MODE_SAMPLES), len(WINDOWS),
                             figsize=(9 * len(WINDOWS), 7.5 * len(MODE_SAMPLES)), squeeze=False)
    for i, s in enumerate(MODE_SAMPLES):
        for j, wn in enumerate(WINDOWS):
            ent = [(f"{s}, all", R[(s, k, wn)], "black", "s")]
            ent += [(g, R[((s, g), k, wn)], c, "o") for g, _, c in MODE_GROUPS]
            draw(axes[i][j], ent, f"{s}: Ds_mass {wn} GeV",
                 rf"$\varepsilon_{{{k}}}\,/\,\varepsilon_{{{BASE}}}$", ylim=(0.0, 2.5))
    hep.cms.label(ax=axes[0][0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary",
                  fontsize=11)
    plt.tight_layout()
    out = common.fig_path("os_ss", f"model_response_by_mode_{k}.png")
    plt.savefig(out, dpi=110)
    plt.close(fig)
    print(f"Saved {out}")

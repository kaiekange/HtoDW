"""ttbar control region (b-tagged): what is in it, is there signal, and the Ds / D peaks with
their charm origin.

Input: the ttcr skim (make_samples_nano.py --region=ttcr): pass_mu & ~jvm_veto & n_bjet >= 1 &
phi window, OS and SS files; extra-lepton columns and, for GenTag MC, cand_hf_from_b
(-1 no truth / 0 prompt charm / 1 charm from a b hadron).

Categories (cut here from the skim columns; n_xmu == 0 is exactly pass_muveto):
    1l_b1      n_xmu == 0, n_xele == 0, n_bjet >= 1
    1l_b2      n_xmu == 0, n_xele == 0, n_bjet >= 2              (main: suppresses tW)
    1l_b2j4    n_xmu == 0, n_xele == 0, n_bjet >= 2, Njet(pT>30) >= 4
    mumu       n_xmu == 1 opposite-sign, |M(mumu) - 91| > 15     (MC x w_xmu)
    emu        n_xmu == 0, n_xele >= 1 opposite-sign (leading electron)
    mumu_b2 / emu_b2   the two above with n_bjet >= 2 (opt-in via --cats, not in the default list)
For each category and charge (OS / SS) one figure: n_bjet, Njet, m_T, MET (+ M(ll) for 2l)
as data vs stack {ttbar 2L2Nu, ttbar LNu2Q, W+jets, DY} + signal line; then Ds_mass and
phi_mass with MC stacked by truth: combinatorial (modes 0,6-10), D->Kpipi misID (5),
true D (3+4) prompt / from b, true Ds (1+2) prompt / from b (origin -1 rows of modes 1-4 go
with "prompt" and are counted in the table).

MC weights: skim weight (genWeight*norm*w_pu*w_mu*w_btag_any, b-tag SF method 1a, added
2026-10-20) x common.br_reweight x common.dy_jetbin_sf (DY) x w_xmu (mumu only). No electron SF.
OS data in the Ds SR (|Ds_mass - M_Ds| < 40 MeV) is hidden unless --unblind (user gate:
first show signal yield and S/B in the SR, printed here).

--eff=0.90: stage-1 score cut (score_dsmu > common.stage1_threshold at that signal efficiency of
--key), read from that model's Prediction/ttcr_<sample>.root (predict.py --region ttcr);
figures get an _eff90 suffix.

    python3 draw_cr_ttbar.py [--cats=1l_b2,mumu] [--unblind] [--key=ssdata_phi20] [--eff=0.90]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

ARGS = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
CR_DIR = ARGS.get("dir", f"{common.SAMPLES_DIR}/CR/ttcr")
MODEL_KEY = ARGS.get("key", "ssdata_phi20")
UNBLIND = "--unblind" in sys.argv
EFF = float(ARGS["eff"]) if "eff" in ARGS else None
if EFF is not None:
    THR, _n = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=EFF))
    PRED_DIR = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
    print(f"stage-1 {MODEL_KEY}: eff_sig {EFF} -> score_dsmu > {THR:.6f}")
TAG = "" if EFF is None else f"_eff{round(100 * EFF)}"
M_DS, M_D, SR_HALF, M_Z = 1.96835, 1.86966, 0.040, 91.19
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0

GROUPS = {  # charge -> {role: file}
    "OS": {"data": "Data_OS_merged", "TT": "TT_merged", "W": "WtoLNu2Jets_merged",
           "DY": "DY_merged", "sig": "Signal_merged"},
    "SS": {"data": "Data_SS_merged", "TT": "TT_SS_merged", "W": "WtoLNu2Jets_SS_merged",
           "DY": "DY_SS_merged"},
}
PROC = [  # (label, role, extra mask, colour) -- stack bottom -> top
    ("DY", "DY", None, "#4a90c2"),
    ("W+jets", "W", None, "#5a9e78"),
    (r"$t\bar{t}$ LNu2Q", "TT", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
]
TRUTH = [  # (label, mask on arrays, colour) -- stack bottom -> top
    ("combinatorial (0,6-10)", lambda a: np.isin(a["cand_gen_mode"], [0, 6, 7, 8, 9, 10]), "#9e9e9e"),
    (r"D$\to$K$\pi\pi$ misID (5)", lambda a: a["cand_gen_mode"] == 5, "#e6c229"),
    ("true D (3+4) prompt", lambda a: np.isin(a["cand_gen_mode"], [3, 4]) & (a["cand_hf_from_b"] != 1), "#f4a582"),
    ("true D (3+4) from b", lambda a: np.isin(a["cand_gen_mode"], [3, 4]) & (a["cand_hf_from_b"] == 1), "#e07a3c"),
    (r"true D$_s$ (1+2) prompt", lambda a: np.isin(a["cand_gen_mode"], [1, 2]) & (a["cand_hf_from_b"] != 1), "#f08080"),
    (r"true D$_s$ (1+2) from b", lambda a: np.isin(a["cand_gen_mode"], [1, 2]) & (a["cand_hf_from_b"] == 1), "#c0392b"),
]
OS_LEP = lambda a, q: q == -a["mu_charge"]
CATS = {
    "1l_b1": lambda a: (a["n_xmu"] == 0) & (a["n_xele"] == 0) & (a["n_bjet_corr"] >= 1),
    "1l_b2": lambda a: (a["n_xmu"] == 0) & (a["n_xele"] == 0) & (a["n_bjet_corr"] >= 2),
    "1l_b2j4": lambda a: ((a["n_xmu"] == 0) & (a["n_xele"] == 0) & (a["n_bjet_corr"] >= 2)
                          & (a["num_jets_pt30_corr"] >= 4)),
    "mumu": lambda a: ((a["n_xmu"] == 1) & OS_LEP(a, a["xmu_charge"])
                       & (np.abs(a["m_mumu"] - M_Z) > 15)),
    "emu": lambda a: (a["n_xmu"] == 0) & (a["n_xele"] >= 1) & OS_LEP(a, a["xele_charge"]),
}
OPT_IN = {"mumu_b2", "emu_b2"}
CATS["mumu_b2"] = lambda a: CATS["mumu"](a) & (a["n_bjet_corr"] >= 2)
CATS["emu_b2"] = lambda a: CATS["emu"](a) & (a["n_bjet_corr"] >= 2)
SEL_CATS = [c for c in ARGS.get("cats", ",".join(c for c in CATS if c not in OPT_IN)).split(",") if c]
CTRL = [  # (column, x label, bins)
    ("n_bjet_corr", r"$n_b$ (tight, $p_T>20$)", np.arange(0.5, 5.5, 1.0)),
    ("num_jets_pt30_corr", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 9.5, 1.0)),
    ("W_tmass_corr", r"$m_T$ [GeV]", np.arange(0.0, 200.0 + 1e-9, 10.0)),
    ("met_pt_corr", r"$p_T^{miss}$ [GeV]", np.arange(0.0, 200.0 + 1e-9, 10.0)),
]
LL = {"mumu": ("m_mumu", r"$M(\mu\mu)$ [GeV]"), "emu": ("m_emu", r"$M(e\mu)$ [GeV]")}
LL.update({"mumu_b2": LL["mumu"], "emu_b2": LL["emu"]})
DS_BINS = np.concatenate([np.linspace(M_DS - 0.30, M_DS - SR_HALF, 14),
                          np.linspace(M_DS - SR_HALF, M_DS + SR_HALF, 5)[1:-1],
                          np.linspace(M_DS + SR_HALF, M_DS + 0.30, 14)])   # 20 MeV, SR = 4 bins
PHI_BINS = np.linspace(common.M_PHI - 0.020, common.M_PHI + 0.020, 11)
COLS = ["Ds_mass", "phi_mass", "weight", "mu_charge", "n_xmu", "n_xele", "xmu_charge",
        "xele_charge", "m_mumu", "m_emu", "w_xmu", "n_bjet_corr", "num_jets_pt30_corr",
        "W_tmass_corr", "met_pt_corr"]
TRUTH_COLS = ["cand_gen_mode", "cand_gen_common_hf_pdgId", "cand_hf_from_b"]


def read(fname, role):
    path = f"{CR_DIR}/{fname}.root" if EFF is None else f"{PRED_DIR}/ttcr_{fname}.root"
    if not os.path.exists(path):
        sys.exit(f"missing {path}")
    with uproot.open(path)["Events"] as t:
        keys = set(t.keys())
        want = COLS + [c for c in TRUTH_COLS + ["tt_decay", "jet_bin"] if c in keys]
        if EFF is not None:
            want.append("score_dsmu")
        a = t.arrays(want, library="np")
    if EFF is not None:
        keep = a["score_dsmu"] > THR
        a = {k: v[keep] for k, v in a.items()}
    a["weight"] = a["weight"].astype(float)
    if role in ("TT", "W", "DY"):
        a["weight"] = a["weight"] * common.br_reweight(a)
    if role == "DY":
        a["weight"] = a["weight"] * common.dy_jetbin_sf(a)
    return a


def hist(x, w, bins):
    x = np.clip(x, bins[0], bins[-1] - 1e-9)
    return np.histogram(x, bins, weights=w)[0], np.histogram(x, bins, weights=w * w)[0]


def panel(ax, rax, comps, data, bins, xlab, blind=None, sig=None):
    """comps: [(label, values, weights, colour)]; data: values; blind: (lo, hi) or None."""
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    bottom, tot, tot2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
    for lab, x, w, col in comps:
        hh, e2 = hist(x, w, bins)
        ax.bar(bc, hh, width=wd, bottom=bottom, color=col, edgecolor="none",
               label=f"{lab}  N={hh.sum():,.1f}")
        bottom += hh; tot += hh; tot2 += e2
    terr = np.sqrt(tot2)
    ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////",
           edgecolor="0.35", linewidth=0.0, label="MC stat")
    if sig is not None:
        hs, _ = hist(sig[0], sig[1], bins)
        hep.histplot(hs, bins, ax=ax, color="red", lw=2, label=f"signal  N={hs.sum():.2f}")
    dh = hist(data, np.ones(len(data)), bins)[0]
    vis = np.ones(len(bc), bool) if blind is None else ~((bc > blind[0]) & (bc < blind[1]))
    dshow = np.where(vis, dh, np.nan)
    ax.errorbar(bc, dshow, yerr=np.sqrt(np.where(vis, dh, 0)), fmt="o", ms=4, color="black",
                label=f"data  N={dh[vis].sum():,.0f}" + (" (SR hidden)" if blind else ""))
    if blind is not None:
        for a_ in (ax, rax):
            a_.axvspan(*blind, color="grey", alpha=0.15, zorder=0)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.5 * max(np.nanmax(dshow) if vis.any() else 0, tot.max(), 1e-9))
    ax.tick_params(labelbottom=False)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylabel("Events / bin", fontsize=14)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(vis & (tot > 0), dh / tot, np.nan)
        re = np.where(vis & (tot > 0), np.sqrt(dh) / tot, np.nan)
        band = np.where(tot > 0, terr / tot, np.nan)
    rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////",
            edgecolor="0.35", linewidth=0.0)
    rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
    rax.axhline(1, color="red", lw=1)
    rax.set_ylim(0.0, 2.0)
    rax.set_xlim(bins[0], bins[-1])
    rax.set_xlabel(xlab)
    rax.set_ylabel("data/MC", fontsize=12)


for charge, files in GROUPS.items():
    S = {role: read(f, role) for role, f in files.items()}
    for cat in SEL_CATS:
        cf = CATS[cat]
        m = {role: cf(a) for role, a in S.items()}
        wl = {role: a["weight"] * (a["w_xmu"] if cat.startswith("mumu") and role != "data" else 1.0)
              for role, a in S.items()}
        mc_roles = ("TT", "W", "DY")
        blind = (M_DS - SR_HALF, M_DS + SR_HALF) if (charge == "OS" and not UNBLIND) else None

        # ---- yield table ----
        print(f"\n===== {charge} {cat}: data {int(m['data'].sum()):,}")
        tot = 0.0
        rows = []
        for lab, role, extra, _c in PROC:
            k = m[role] & (extra(S[role]) if extra else True)
            y = wl[role][k].sum()
            rows.append((lab, y)); tot += y
        for lab, y in rows:
            plain = lab.replace(r"$t\bar{t}$", "ttbar")
            print(f"   {plain:14s} {y:11,.1f}  ({100 * y / tot:5.1f}% of MC)")
        print(f"   data/MC = {m['data'].sum() / tot:.3f}")
        sr = lambda a: np.abs(a["Ds_mass"] - M_DS) < SR_HALF
        dwin = lambda a: np.abs(a["Ds_mass"] - M_D) < 0.030
        if "sig" in S:
            ys = wl["sig"][m["sig"] & sr(S["sig"])].sum()
            ybsr = sum(wl[r][m[r] & sr(S[r])].sum() for r in mc_roles)
            print(f"   signal in Ds SR: {ys:.4f}   MC bkg in SR {ybsr:,.1f}   S/B {ys / ybsr:.2e}")
        for name, modes, win in (("true Ds (1+2) in SR", [1, 2], sr), ("true D (3+4) in D+- +-30", [3, 4], dwin)):
            yp = yb = yu = 0.0
            for r in mc_roles:
                a = S[r]
                k = m[r] & win(a) & np.isin(a["cand_gen_mode"], modes)
                o = a["cand_hf_from_b"][k]
                yp += wl[r][k][o == 0].sum(); yb += wl[r][k][o == 1].sum(); yu += wl[r][k][o < 0].sum()
            t_ = yp + yb + yu
            print(f"   {name:26s} MC {t_:9.1f}: prompt {100 * yp / max(t_, 1e-9):5.1f}%  from b "
                  f"{100 * yb / max(t_, 1e-9):5.1f}%  no match {100 * yu / max(t_, 1e-9):4.1f}%")

        # ---- figure: control variables (row 1), Ds_mass and phi_mass by truth (row 2) ----
        ctrl = CTRL + ([(LL[cat][0], LL[cat][1], np.arange(12.0, 252.0 + 1e-9, 8.0))]
                       if cat in LL else [])
        ncol = max(len(ctrl), 2)
        fig = plt.figure(figsize=(6.6 * ncol, 15.5))
        outer = gridspec.GridSpec(2, ncol, figure=fig, wspace=0.32, hspace=0.30,
                                  top=0.93, bottom=0.06)
        for i, (col, xlab, bins) in enumerate(ctrl):
            inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, i],
                                                     height_ratios=[3, 1], hspace=0.06)
            ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            comps = []
            for lab, role, extra, colr in PROC:
                k = m[role] & (extra(S[role]) if extra else True)
                comps.append((lab, S[role][col][k], wl[role][k], colr))
            sig = ((S["sig"][col][m["sig"]], wl["sig"][m["sig"]]) if "sig" in S else None)
            panel(ax, rax, comps, S["data"][col][m["data"]], bins, xlab, sig=sig)
            if i == 0:
                hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}",
                              com=13.6, fontsize=13)
        for i, (col, xlab, bins) in enumerate((("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", DS_BINS),
                                                ("phi_mass", r"$M(K^+K^-)$ [GeV]", PHI_BINS))):
            inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1, i],
                                                     height_ratios=[3, 1], hspace=0.06)
            ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
            comps = []
            for lab, tf, colr in TRUTH:
                xs, ws = [], []
                for r in mc_roles:
                    k = m[r] & tf(S[r])
                    xs.append(S[r][col][k]); ws.append(wl[r][k])
                comps.append((lab, np.concatenate(xs), np.concatenate(ws), colr))
            dmask = m["data"]
            if blind is not None and col == "phi_mass":   # event-level blind: drop SR events
                dmask = dmask & ~sr(S["data"])
            sig = ((S["sig"][col][m["sig"]], wl["sig"][m["sig"]]) if "sig" in S else None)
            panel(ax, rax, comps, S["data"][col][dmask], bins, xlab,
                  blind=blind if col == "Ds_mass" else None, sig=sig)
            if blind is not None and col == "phi_mass":
                ax.set_title("data: Ds SR events removed (MC not)", fontsize=12)
        cut_txt = "" if EFF is None else f" & stage-1 eff {round(100 * EFF)}%"
        fig.suptitle(f"ttbar CR {charge}, {cat}: pass_mu & b-tag & phi window{cut_txt} "
                     f"(MC: BR rw, DY jet-bin SF, b-tag SF)", fontsize=16, y=0.975)
        out = common.fig_path("sf", f"cr_ttbar_{charge}_{cat}{TAG}.png", key=MODEL_KEY)
        plt.savefig(out, dpi=90, bbox_inches="tight")
        plt.close(fig)
        print(f"   Saved {out}")

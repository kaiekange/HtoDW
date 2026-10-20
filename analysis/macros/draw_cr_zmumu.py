"""Z->mumu control region: DY normalisation and what else sits under the Z.

Input: the zmumu CR skim (make_samples_nano.py --region=zmumu; selection documented at
REGIONS there): best muon + exactly one opposite-charge extra muon, ~pass_muveto, b-veto,
no phi window. Data = OS files only (an event can sit in both the OS and the SS file with
different candidates, so using both would double count).

One figure per version, panels = M(mumu), pt(mu1), pt(mu2), pt(mumu), reco Njet, MET, m_T
(all inside 76 < M(mumu) < 106, where k_DY is defined), each
data vs stack {DY, W+jets, ttbar 2L2Nu, ttbar LNu2Q} + signal line, with data/MC below.
mu1 = the best muon (corrected pt; it is the leading one by construction), mu2 = the extra
muon (raw pt):
    phi20       |phi_mass - M_phi| < 20 MeV, the analysis selection (the zmumu skim has no
                phi cut, so it is applied here; the no-phi version was dropped 2026-10-19)
No DNN cut: at eps90 the window keeps ~400 data events, too few (user 2026-10-19).
Printed per version: yield of every process in 76 < M(mumu) < 106 and the single overall
k_DY = (N_data - N_nonDY) / N_DY there. The per-jet-bin DY factors are fitted by
fit_dy_jetbin_zmumu.py.

The second muon carries no ID/iso SF and uses raw pt (first pass); the factor therefore also
contains the data/MC of "a KKpi candidate was found" (custom nano keeps only such events).

Second figure cr_zmumu{TAG}_phi20_dsmass.png: Ds_mass and phi_mass in the Z window (phi20),
MC x common.br_reweight x common.dy_jetbin_sf (the factors fitted here and in Step 1; the
first figure stays on the bare skim weight it was fitted on), stacked by GenTag truth
(no prompt/from-b split: the zmumu skim has no cand_hf_from_b). Not blinded (CR).

--eff=0.90: stage-1 score cut (score_dsmu > common.stage1_threshold at that signal efficiency of
--key), read from that model's Prediction/zmumu_<group>.root; figure names get _eff90.

Third figure cr_zmumu{TAG}_phi20_tau3pi.png: the tau -> 3-prong misID masses of
draw_tau3pi_masses_by_cut.py (track 4-vectors, mass hypothesis in brackets): M(pi[pi] K-[pi])
(rho0), M(K+[pi] K-[pi] pi[pi]) (a1, <= m_tau), M(K-[K] pi) (K*0), M(K+[K] pi); same selection,
weights and truth stack as the Ds_mass figure; 20 MeV bins (80 MeV with --eff).

    python3 draw_cr_zmumu.py [--dir=<CR dir>] [--key=ssdata_phi20] [--eff=0.90]
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
CR_DIR = ARGS.get("dir", f"{common.SAMPLES_DIR}/CR/zmumu")
MODEL_KEY = ARGS.get("key", "ssdata_phi20")
TAG = ARGS.get("tag", "")
EFF = float(ARGS["eff"]) if "eff" in ARGS else None
if EFF is not None:
    THR, _n = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=EFF))
    PRED_DIR = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
    TAG = f"{TAG}_eff{round(100 * EFF)}"
    print(f"stage-1 {MODEL_KEY}: eff_sig {EFF} -> score_dsmu > {THR:.6f}")
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0

BKG = [  # (label, file, extra mask on the file's arrays or None, colour); stack bottom -> top
    (r"$t\bar{t}$ LNu2Q", "TT_merged", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT_merged", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
    ("W+jets", "WtoLNu2Jets_merged", None, "#5a9e78"),
    (r"DY$\to\ell\ell$", "DY_merged", None, "#4a90c2"),
]
VARS = [  # (column, x label, bin edges)
    ("m_mumu", r"$M(\mu\mu)$ [GeV]", np.arange(76.0, 106.0 + 1e-9, 1.0)),
    ("mu_pt_corr", r"$p_T(\mu_1)$ [GeV]", np.arange(20.0, 150.0 + 1e-9, 5.0)),
    ("xmu_pt", r"$p_T(\mu_2)$ [GeV]", np.arange(10.0, 120.0 + 1e-9, 5.0)),
    ("ptll", r"$p_T(\mu\mu)$ [GeV]", np.arange(0.0, 150.0 + 1e-9, 5.0)),
    ("num_jets_pt30_corr", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 6.5, 1.0)),
    ("met_pt_corr", r"$p_T^{miss}$ [GeV]", np.arange(0.0, 150.0 + 1e-9, 5.0)),
    ("W_tmass_corr", r"$m_T$ [GeV]", np.arange(0.0, 150.0 + 1e-9, 5.0)),
]
DERIVED = {"ptll"}   # computed in read() from the two muons
COLS = sorted(({v for v, _l, _b in VARS} - DERIVED)
              | {"phi_mass", "weight", "event", "run", "luminosityBlock",
                 "mu_eta", "mu_phi", "xmu_eta", "xmu_phi", "Ds_mass"}
              | {f"{t}_{v}" for t in ("Kp", "Km", "pi") for v in ("pt", "eta", "phi")})
M_PI_H, M_K_H = 0.13957039, 0.493677
TRUTH_COLS = ["cand_gen_mode", "cand_gen_common_hf_pdgId", "jet_bin"]


def read(group):
    path = f"{CR_DIR}/{group}.root" if EFF is None else f"{PRED_DIR}/zmumu_{group}.root"
    if not os.path.exists(path):
        sys.exit(f"missing {path}")
    with uproot.open(path)["Events"] as t:
        keys = set(t.keys())
        want = [c for c in COLS + ["tt_decay"] + TRUTH_COLS if c in keys]
        if EFF is not None:
            want.append("score_dsmu")
        a = t.arrays(want, library="np")
    if EFF is not None:
        keep = a["score_dsmu"] > THR
        a = {k: v[keep] for k, v in a.items()}
    px = a["mu_pt_corr"] * np.cos(a["mu_phi"]) + a["xmu_pt"] * np.cos(a["xmu_phi"])
    py = a["mu_pt_corr"] * np.sin(a["mu_phi"]) + a["xmu_pt"] * np.sin(a["xmu_phi"])
    a["ptll"] = np.hypot(px, py)
    a["pipi"] = trk_mass(a, [("pi", M_PI_H), ("Km", M_PI_H)])
    a["pipipi"] = trk_mass(a, [("Kp", M_PI_H), ("Km", M_PI_H), ("pi", M_PI_H)])
    a["Kmpi"] = trk_mass(a, [("Km", M_K_H), ("pi", M_PI_H)])
    a["Kppi"] = trk_mass(a, [("Kp", M_K_H), ("pi", M_PI_H)])
    return a


def trk_mass(a, tracks):
    """Invariant mass from pt/eta/phi; tracks = [(prefix, mass hypothesis)] (as draw_tau3pi_masses_by_cut.py)."""
    E = px = py = pz = 0.0
    for p, m0 in tracks:
        pt, eta, phi = a[f"{p}_pt"], a[f"{p}_eta"], a[f"{p}_phi"]
        E = E + np.sqrt((pt * np.cosh(eta)) ** 2 + m0 ** 2)
        px, py, pz = px + pt * np.cos(phi), py + pt * np.sin(phi), pz + pt * np.sinh(eta)
    return np.sqrt(np.maximum(E ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))


def versions(S):
    # phi window always (user 2026-10-19): the skim carries no phi cut, it is applied here
    return [("phi20", {g: np.abs(a["phi_mass"] - M_PHI) < 0.020 for g, a in S.items()})]


def hist(x, w, bins):
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w * w)
    return h, e2


def draw(S, vname, sel):
    d = S["Data_OS_merged"]
    zw = lambda a, m: m & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    print(f"\n== version {vname}: yields in {ZLO:.0f} < M(mumu) < {ZHI:.0f}")
    nd = zw(d, sel["Data_OS_merged"]).sum()
    yields = {}
    for lab, g, extra, _c in BKG:
        a = S[g]
        m = zw(a, sel[g]) & (extra(a) if extra else True)
        yields[lab] = (a["weight"][m].sum(), (a["weight"][m] ** 2).sum())
    s = S["Signal_merged"]
    ysig = s["weight"][zw(s, sel["Signal_merged"])].sum()
    ntot = sum(v for v, _ in yields.values())
    for lab, (v, _e2) in yields.items():
        print(f"   {lab:22s} {v:12,.1f}  ({100 * v / ntot:5.2f}% of MC)")
    print(f"   {'signal':22s} {ysig:12,.3f}  ({100 * ysig / ntot:.4f}% of MC)")
    ndy, e2dy = yields[r"DY$\to\ell\ell$"]
    nnon = sum(v for k, (v, _) in yields.items() if k != r"DY$\to\ell\ell$")
    e2non = sum(e for k, (_, e) in yields.items() if k != r"DY$\to\ell\ell$")
    k = (nd - nnon) / ndy
    ek = np.sqrt(nd + e2non + k * k * e2dy) / ndy
    print(f"   data {nd:,}   MC {ntot:,.1f}   data/MC {nd / ntot:.4f}   "
          f"k_DY = {k:.4f} +- {ek:.4f} (stat)")

    ncol = 4
    nrow = (len(VARS) + ncol - 1) // ncol
    fig = plt.figure(figsize=(7.2 * ncol, 7.6 * nrow))
    outer = gridspec.GridSpec(nrow, ncol, figure=fig, wspace=0.32, hspace=0.30,
                              top=0.92, bottom=0.07)
    for col, (var, xlab, bins) in enumerate(VARS):
        # every panel is drawn inside the Z window, where k_DY is defined (user 2026-10-19):
        # below M ~ 50 the data hold J/psi, Upsilon, low-mass DY and heavy-flavour dimuons that
        # no MC sample here covers (DY is M_ll > 50), and they would distort every shape.
        vs = {g: zw(a_, sel[g]) for g, a_ in S.items()}
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[col],
                                                 height_ratios=[3, 1], hspace=0.06)
        ax, rax = (fig.add_subplot(inner[i]) for i in range(2))
        bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
        bottom, tot, tot_e2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for lab, g, extra, colr in BKG:
            a = S[g]
            m = vs[g] & (extra(a) if extra else True)
            h, e2 = hist(a[var][m], a["weight"][m].astype(float), bins)
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none",
                   label=f"{lab}  N={h.sum():,.0f}")
            bottom += h; tot += h; tot_e2 += e2
        terr = np.sqrt(tot_e2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////",
               edgecolor="0.35", linewidth=0.0, label="MC stat")
        hs, _ = hist(s[var][vs["Signal_merged"]], s["weight"][vs["Signal_merged"]].astype(float), bins)
        hep.histplot(hs, bins, ax=ax, color="red", lw=2, label=f"signal  N={hs.sum():.2f}")
        dh, _ = np.histogram(d[var][vs["Data_OS_merged"]], bins)
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black",
                    label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 1.45 * max(tot.max(), dh.max(), 1))
        ax.legend(fontsize=9, loc="upper right")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel(f"Events / bin", fontsize=16)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////",
                edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.5, 1.5)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}",
                          com=13.6, fontsize=14)
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR ({vname}): k_DY = {k:.3f} $\\pm$ {ek:.3f} (stat)   "
                 f"[2nd $\\mu$: raw $p_T$, no SF]", fontsize=15, y=0.975)
    out = common.fig_path("sf", f"cr_zmumu{TAG}_{vname}.png", key=MODEL_KEY)
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"   Saved {out}")


M_DS, M_D, SR_HALF = 1.96835, 1.86966, 0.040
TRUTH = [  # (label, modes, colour) -- stack bottom -> top, as draw_cr_ttbar.py
    ("combinatorial (0,6-10)", [0, 6, 7, 8, 9, 10], "#9e9e9e"),
    (r"D$\to$K$\pi\pi$ misID (5)", [5], "#e6c229"),
    ("true D (3+4)", [3, 4], "#f4a582"),
    (r"true D$_s$ (1+2)", [1, 2], "#c0392b"),
]
if EFF is None:   # 20 MeV Ds_mass (SR = 4 bins), 2 MeV phi_mass
    _ds_bins = np.concatenate([np.linspace(M_DS - 0.30, M_DS - SR_HALF, 14),
                               np.linspace(M_DS - SR_HALF, M_DS + SR_HALF, 5)[1:-1],
                               np.linspace(M_DS + SR_HALF, M_DS + 0.30, 14)])
    _phi_bins = np.linspace(M_PHI - 0.020, M_PHI + 0.020, 21)
else:             # with a DNN cut: 40 MeV Ds_mass (SR = 2 bins, edges on the SR), 4 MeV phi_mass
    _ds_bins = M_DS - SR_HALF + 0.040 * np.arange(-6, 9)
    _phi_bins = np.linspace(M_PHI - 0.020, M_PHI + 0.020, 11)
MASS = [("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", _ds_bins),
        ("phi_mass", r"$M(K^+K^-)$ [GeV]", _phi_bins)]


_bw = 0.020 if EFF is None else 0.080
TAU = [("pipi", r"$M(\pi\pi)$ of $\pi K^-$, $K^-\to\pi$ [GeV]", np.arange(0.20, 1.64 + 1e-9, _bw)),
       ("pipipi", r"$M(\pi\pi\pi)$, $K^\pm\to\pi$ [GeV]", np.arange(0.80, 2.30 + 1e-9, _bw)),
       ("Kmpi", r"$M(K^-\pi)$ (opposite sign) [GeV]", np.arange(0.60, 1.80 + 1e-9, _bw)),
       ("Kppi", r"$M(K^+\pi)$ (same sign) [GeV]", np.arange(0.60, 1.80 + 1e-9, _bw))]


def draw_mass(S, vname, sel, panels=MASS, suffix="dsmass"):
    zw = lambda a, m: m & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    mc = ("TT_merged", "WtoLNu2Jets_merged", "DY_merged")
    w = {}
    for g in mc:
        a = S[g]
        w[g] = a["weight"].astype(float) * common.br_reweight(a)
        if g == "DY_merged":
            w[g] = w[g] * common.dy_jetbin_sf(a)
    d = S["Data_OS_merged"]
    dm = zw(d, sel["Data_OS_merged"])
    print(f"\n== {vname} {suffix} (MC x BR rw x DY jet-bin SF), Z window")
    for wname, lo, hi in (("Ds SR +-40", M_DS - SR_HALF, M_DS + SR_HALF), ("D+- +-30", M_D - 0.03, M_D + 0.03)):
        win = lambda a: (a["Ds_mass"] > lo) & (a["Ds_mass"] < hi)
        parts = {lab: sum(w[g][zw(S[g], sel[g]) & win(S[g]) & np.isin(S[g]["cand_gen_mode"], modes)].sum()
                          for g in mc) for lab, modes, _c in TRUTH}
        tot = sum(parts.values()); nd = int((dm & win(d)).sum())
        print(f"   {wname}: data {nd:,}  MC {tot:,.1f}  data/MC {nd / tot:.3f}  |  "
              + "  ".join(f"{k.replace('$', '').replace(chr(92), '')} {v:,.1f}" for k, v in parts.items()))
    fig = plt.figure(figsize=(7.2 * len(panels), 7.6))
    outer = gridspec.GridSpec(1, len(panels), figure=fig, wspace=0.32, top=0.90, bottom=0.10)
    for col, (var, xlab, bins) in enumerate(panels):
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[col],
                                                 height_ratios=[3, 1], hspace=0.06)
        ax, rax = (fig.add_subplot(inner[i]) for i in range(2))
        bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
        bottom, tot, tot_e2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for lab, modes, colr in TRUTH:
            xs, ws = [], []
            for g in mc:
                m = zw(S[g], sel[g]) & np.isin(S[g]["cand_gen_mode"], modes)
                xs.append(S[g][var][m]); ws.append(w[g][m])
            h, e2 = hist(np.concatenate(xs), np.concatenate(ws), bins)
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none",
                   label=f"{lab}  N={h.sum():,.1f}")
            bottom += h; tot += h; tot_e2 += e2
        terr = np.sqrt(tot_e2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////",
               edgecolor="0.35", linewidth=0.0, label="MC stat")
        dh, _ = np.histogram(d[var][dm], bins)
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black",
                    label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 1.45 * max(tot.max(), dh.max(), 1))
        ax.legend(fontsize=9, loc="upper right")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Events / bin", fontsize=16)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////",
                edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.5, 1.5)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}",
                          com=13.6, fontsize=14)
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR ({vname}), {ZLO:.0f} < M($\\mu\\mu$) < {ZHI:.0f}: "
                 f"MC x BR rw x DY jet-bin SF, "
                 f"{'no DNN cut' if EFF is None else f'stage-1 eff {round(100 * EFF)}%'}", fontsize=15, y=0.975)
    out = common.fig_path("sf", f"cr_zmumu{TAG}_{vname}_{suffix}.png", key=MODEL_KEY)
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"   Saved {out}")


S = {g: read(g) for g in ("Data_OS_merged", "Signal_merged", "WtoLNu2Jets_merged",
                         "TT_merged", "DY_merged")}
for vname, sel in versions(S):
    draw(S, vname, sel)
    draw_mass(S, vname, sel)
    draw_mass(S, vname, sel, TAU, "tau3pi")

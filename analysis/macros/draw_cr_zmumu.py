"""Z->mumu control region: DY normalisation and what else sits under the Z.

Input: the zmumu CR skim (make_samples_nano.py --region=zmumu; selection documented at
REGIONS there): best muon + exactly one opposite-charge extra muon, ~pass_muveto, b-veto,
no phi window. Data = OS files only (an event can sit in both the OS and the SS file with
different candidates, so using both would double count).

One figure per version, panels = M(mumu), pt(mu1), pt(mu2), pt(mumu), reco Njet, MET, m_T
(all inside 76 < M(mumu) < 106, where k_DY is defined), each
data vs stack {DY, W+jets, ttbar 2L2Nu, ttbar LNu2Q} + signal line, with data/MC below.
mu1 = the best muon (corrected pt; it is the leading one by construction), mu2 = the extra
muon (MuonScaRe-corrected pt, ID x loose-iso SF in the skim weight since the 2026-10-19 skim):
    phi20       |phi_mass - M_phi| < 20 MeV, the analysis selection (the zmumu skim has no
                phi cut, so it is applied here; the no-phi version was dropped 2026-10-19)
No DNN cut: at eps90 the window keeps ~400 data events, too few (user 2026-10-19).
MC weight (2026-10-21): skim weight x common.br_reweight x common.dy_jetbin_sf (DY), as in the
mass figures below, so the printed k_DY = (N_data - N_nonDY) / N_DY is the RESIDUAL factor
after those corrections. --bare: skim weight only (the input the per-jet-bin DY factors were
fitted on by fit_dy_jetbin_zmumu.py); figure name gets _bare.
k_DY also contains the data/MC of "a KKpi candidate was found" (custom nano keeps only such
events).

Second figure cr_zmumu{TAG}_phi20_dsmass.png: Ds_mass and phi_mass in the Z window (phi20),
MC x common.br_reweight x common.dy_jetbin_sf (the factors fitted here and in Step 1), stacked by GenTag truth
(no prompt/from-b split: the zmumu skim has no cand_hf_from_b). Not blinded (CR).

--eff=0.90: stage-1 score cut (score_dsmu > common.stage1_threshold at that signal efficiency of
--key), read from that model's Prediction/zmumu_<group>.root; figure names get _eff90.

Third figure cr_zmumu{TAG}_phi20_tau3pi.png: the tau -> 3-prong misID masses of
draw_tau3pi_masses_by_cut.py (track 4-vectors, mass hypothesis in brackets): M(pi[pi] K-[pi])
(rho0), M(K+[pi] K-[pi] pi[pi]) (a1, <= m_tau), M(K-[K] pi) (K*0), M(K+[K] pi); same selection,
weights and truth stack as the Ds_mass figure; 20 MeV bins (80 MeV with --eff).

    python3 draw_cr_zmumu.py [--dir=<CR dir>] [--key=ssdata_phi20] [--eff=0.90] [--ss]
                             [--key2=resdy --eff2=0.70]   (stage-2 cut on top, see KEY2 below)
--ss (2026-10-21): the same-sign files (Data/WtoLNu2Jets/TT/DY _SS_merged); there is no SS
signal file in this CR, so the signal line is empty. Figure names get _ss.
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
CUT_DESC = "no DNN cut" if EFF is None else f"stage-1 eff {round(100 * EFF)}%"
# --fail (2026-10-26): with --eff, keep the rows the stage-1 cut REJECTS (score_dsmu <= thr, incl.
# the -777 rows), i.e. the complement: pass + fail = no cut. No-cut binning. Not with --key2.
FAIL = "--fail" in sys.argv
if FAIL:
    if EFF is None or "key2" in ARGS:
        sys.exit("--fail needs --eff and no --key2")
    TAG = f"{TAG}_fail"
    CUT_DESC = f"FAILING stage-1 eff {round(100 * EFF)}% (score_dsmu <= {THR:.4f})"
# --key2=<stage-2 BDT key> --eff2=X (2026-10-22): also score_dsmu_res > thr2, read from that BDT's
# Prediction/zmumu_<group>.root (which carries score_dsmu too). thr2 is UNCONDITIONAL, as in
# draw_wp_mass.py: the (1-X) quantile of score_dsmu_res over all truth-matched SR signal
# (selection_mask, both scores > -1, phi20); eps_tot (both cuts, weighted) is printed. The
# zmumu signal is not scored by stage-2 (signal only in the SR), so the signal line is empty.
# Without --eff (user 2026-10-22): the stage-2 cut ALONE, no stage-1 cut; eps_tot = that cut's
# weighted signal efficiency.
KEY2 = ARGS.get("key2")
EFF2 = float(ARGS["eff2"]) if "eff2" in ARGS else None
if (KEY2 is None) != (EFF2 is None):
    sys.exit("--key2 and --eff2 go together")
if KEY2 is not None:
    PRED_DIR = f"{common.out_dir_for(KEY2, tool='BDT')}/Prediction"
    with uproot.open(f"{PRED_DIR}/Signal_merged.root")["Events"] as _t:
        _g = _t.arrays(["score_dsmu", "score_dsmu_res", "phi_mass", "Ds_match", "mu_match", "weight",
                        *common.SELECTION], library="np")
    _m = (common.selection_mask(_g) & (_g["score_dsmu"] > -1) & (_g["score_dsmu_res"] > -1)
          & _g["Ds_match"].astype(bool) & _g["mu_match"].astype(bool)
          & (np.abs(_g["phi_mass"] - 1.019460) < 0.020))
    THR2 = float(np.quantile(_g["score_dsmu_res"][_m], 1 - EFF2))
    _sw = _g["weight"][_m]
    _pass = _g["score_dsmu_res"][_m] > THR2
    if EFF is not None:
        _pass &= _g["score_dsmu"][_m] > THR
    EPS_TOT = float(_sw[_pass].sum() / _sw.sum())
    TAG = f"{TAG}_{KEY2}{round(100 * EFF2)}"
    CUT_DESC = (("" if EFF is None else CUT_DESC + " + ")
                + f"{KEY2} WP2 {round(100 * EFF2)}% (eps_tot {100 * EPS_TOT:.1f}%)"
                + (", no stage-1 cut" if EFF is None else ""))
    print(f"stage-2 {KEY2}: WP2 {EFF2} (unconditional) -> score_dsmu_res > {THR2:.6f}, "
          f"eps_tot {EPS_TOT:.4f}")
FIG_KEY = dict(key=MODEL_KEY) if KEY2 is None else dict(key=KEY2, tool="BDT")
SS = "--ss" in sys.argv
BARE = "--bare" in sys.argv   # first figure only: skim weight without BR rw / DY jet-bin SF
if SS:
    TAG = f"{TAG}_ss"
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0

BKG = [  # (label, file, extra mask on the file's arrays or None, colour); stack bottom -> top
    # DY and ttbar under W+jets (user 2026-10-22)
    (r"DY$\to\ell\ell$", "DY_merged", None, "#4a90c2"),
    (r"$t\bar{t}$ LNu2Q", "TT_merged", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT_merged", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
    ("W+jets", "WtoLNu2Jets_merged", None, "#5a9e78"),
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
    fname = group
    if SS and group != "Signal_merged":
        fname = "Data_SS_merged" if group == "Data_OS_merged" else group.replace("_merged", "_SS_merged")
    scored = EFF is not None or KEY2 is not None
    path = f"{PRED_DIR}/zmumu_{fname}.root" if scored else f"{CR_DIR}/{fname}.root"
    no_s2 = KEY2 is not None and group == "Signal_merged"   # not scored by stage-2: empty line
    if no_s2:   # read the plain skim and keep no row
        path = f"{CR_DIR}/{fname}.root"
    if not os.path.exists(path):
        sys.exit(f"missing {path}")
    with uproot.open(path)["Events"] as t:
        keys = set(t.keys())
        want = [c for c in COLS + ["tt_decay", "mu_charge", "pi_charge"] + TRUTH_COLS if c in keys]
        if EFF is not None and not no_s2:
            want.append("score_dsmu")
        if KEY2 is not None and not no_s2:
            want.append("score_dsmu_res")
        a = t.arrays(want, library="np")
    common.check_charge_config(a, SS and group != "Signal_merged", fname)   # exits on an OS/SS mix-up
    if scored:
        keep = np.ones(len(a["weight"]), dtype=bool)
        if no_s2:
            keep[:] = False
        else:
            if EFF is not None:
                keep &= (a["score_dsmu"] <= THR) if FAIL else (a["score_dsmu"] > THR)
            if KEY2 is not None:
                keep &= a["score_dsmu_res"] > THR2
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


def mc_weight(g, a):
    """skim weight x BR rw x DY jet-bin SF (DY), or the bare skim weight with --bare."""
    w = a["weight"].astype(float)
    if BARE:
        return w
    w = w * common.br_reweight(a)
    return w * common.dy_jetbin_sf(a) if g == "DY_merged" else w


def draw(S, vname, sel):
    d = S["Data_OS_merged"]
    W = {g: mc_weight(g, S[g]) for _l, g, _e, _c in BKG}
    zw = lambda a, m: m & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    print(f"\n== version {vname}: yields in {ZLO:.0f} < M(mumu) < {ZHI:.0f}")
    nd = zw(d, sel["Data_OS_merged"]).sum()
    yields = {}
    for lab, g, extra, _c in BKG:
        a = S[g]
        m = zw(a, sel[g]) & (extra(a) if extra else True)
        yields[lab] = (W[g][m].sum(), (W[g][m] ** 2).sum())
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
            h, e2 = hist(a[var][m], W[g][m], bins)
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
                    label=f"{'SS' if SS else 'OS'} data  N={dh.sum():,}")
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
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR ({vname}): "
                 + ("MC = skim weight only, " if BARE else "MC x BR rw x DY jet-bin SF, residual ")
                 + f"k_DY = {k:.3f} $\\pm$ {ek:.3f} (stat)   "
                 + f"[{CUT_DESC}]",
                 fontsize=15, y=0.975)
    out = common.fig_path("sf", f"cr_zmumu{TAG}_{vname}{'_bare' if BARE else ''}.png", **FIG_KEY)
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
if EFF is None or FAIL:   # 5 MeV Ds_mass (SR = 16 bins), 0.5 MeV phi_mass (user 2026-10-26: stats allow it)
    _ds_bins = np.concatenate([np.linspace(M_DS - 0.30, M_DS - SR_HALF, 53),
                               np.linspace(M_DS - SR_HALF, M_DS + SR_HALF, 17)[1:-1],
                               np.linspace(M_DS + SR_HALF, M_DS + 0.30, 53)])
    _phi_bins = np.linspace(M_PHI - 0.020, M_PHI + 0.020, 81)
else:             # with a DNN cut: 40 MeV Ds_mass (SR = 2 bins, edges on the SR), 4 MeV phi_mass
    _ds_bins = M_DS - SR_HALF + 0.040 * np.arange(-6, 9)
    _phi_bins = np.linspace(M_PHI - 0.020, M_PHI + 0.020, 11)
MASS = [("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", _ds_bins),
        ("phi_mass", r"$M(K^+K^-)$ [GeV]", _phi_bins)]


_bw = 0.020 if EFF is None or FAIL else 0.080
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
                    label=f"{'SS' if SS else 'OS'} data  N={dh.sum():,}")
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
                 f"{CUT_DESC}", fontsize=15, y=0.975)
    out = common.fig_path("sf", f"cr_zmumu{TAG}_{vname}_{suffix}.png", **FIG_KEY)
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"   Saved {out}")


S = {g: read(g) for g in ("Data_OS_merged", "Signal_merged", "WtoLNu2Jets_merged",
                         "TT_merged", "DY_merged")}
if SS:   # no SS signal file: keep the columns, drop every row (empty signal line)
    S["Signal_merged"] = {k: v[:0] for k, v in S["Signal_merged"].items()}
for vname, sel in versions(S):
    draw(S, vname, sel)
    draw_mass(S, vname, sel)
    draw_mass(S, vname, sel, TAU, "tau3pi")

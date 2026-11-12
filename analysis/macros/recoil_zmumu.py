"""Hadronic-recoil calibration of DY MC from the Z->mumu CR (2026-11-12 version: gen boson for the MC, no q_T cut,
JER-in-MET input only; replaces the 11-09 maps).
Input MET = the analysis chain before recoil: Type-1 PuppiMET (JES + JER) + MET xy = skim
met_{pt,phi}_corr_prerecoil; N_jet = num_jets_pt30_corr (JES + JER jets, pt > 30, tight ID).
Recoil u = -(MET + V_vis), projected on q_hat (make_samples_nano.recoil_project):
    data, TT : q = V_vis = the reco dimuon (nano muon pts: the MET carries no muon correction)
    DY       : q = gen boson V, V_vis = V minus the hard-process neutrinos (skim gen_V_*, gen_Vvis_*, from
               make_samples_nano.gen_boson_cols), i.e. the definition used in the SR (W+jets: W, signal: H)
The skim's met_{pt,phi}_corr is already recoil-corrected with recoil_<era>.npz (sha1 in the skim meta), so a
re-derivation must reproduce that map exactly: main() refuses to overwrite a different map (the skim would have
to be re-run) and checks the DY met_pt_corr against the map.
    x_par = u_par + q_T (keeps the q_T spread inside a cell out of the mapping), u_perp.
Method (user, 2026-11-12; JME POG: "map the simulated and measured cumulative distributions, which also accounts
for non-Gaussian shapes"): per era, cells q_T [0,10,20,30,50,80,inf) x N_jet {<=1, 2, >=3}, no q_T cut; in each
cell the DY CDF of x_par and of u_perp is mapped onto the CDF of data - TT (1 GeV histograms in [-300, 300],
cumulative forced monotonic because DY has negative weights). A cell with < MIN_DATA data events uses the
N_jet-inclusive CDFs (data and DY) of its q_T bin. Applied by make_samples_nano.recoil_correct (values outside
[-300, 300] unchanged). Only DY is corrected here. Nominal only: uncertainties and per-JME-variation maps later
(the skim keeps the pre-recoil inputs).
Deployed maps: all DY events -> Samples/2022/CR/zmumu_recoil/recoil_<era>.npz. Closure figures: maps from DY even
events, DY shown = odd events x 2; per-cell closure binned on reco quantities for every sample (like for like).
Selection and weights: zmumu skim, phi20, 76 < M(mumu) < 106, OS data; MC = skim weight x BR rw (DY also x
DY jet-bin SF, common values; only the 0J/1J/2J mixture inside a cell matters).

    python3 recoil_zmumu.py          (conda ds_dnn) maps + figures
Figures MVA/figures/sf/cr_zmumu_recoil_{met,mt,met_qtsplit,uperp}_<era>.png,
cr_zmumu_recoil_cells_{xpar,uperp,met}_<era>.png, cr_zmumu_jme_recoil_final[_<era>].png
"""
import os
import sys

import numpy as np
import uproot

import make_samples_nano as msn

SAMPLES = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples/2022"
CR = f"{SAMPLES}/CR/zmumu"
OUT = f"{SAMPLES}/CR/zmumu_recoil"
GROUPS = ("Data_OS_merged", "TT_merged", "DY_merged")
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
QT = np.array([0, 10, 20, 30, 50, 80, np.inf])
NJ = ((0, 1), (2, 2), (3, 99))
NJ_LAB = ("N$_{jet}\\leq$1", "N$_{jet}$=2", "N$_{jet}\\geq$3")
UB = np.arange(-300, 300 + 1e-9, 1.0)
MIN_DATA = 200   # cells with fewer data events use the N_jet-inclusive CDFs of their q_T bin
ERAS = (("2022", 1), ("2022EE", 2))
COMPS = ("xpar", "uperp")


def njet_class(nj):
    c = np.full(len(nj), -1)
    for k, (lo, hi) in enumerate(NJ):
        c[(nj >= lo) & (nj <= hi)] = k
    return c


def read(group):
    cols = ["run", "luminosityBlock", "event", "era", "weight", "m_mumu", "phi_mass", "mu_pt", "mu_phi",
            "mu_pt_corr", "xmu_pt_raw", "xmu_phi", "mu_charge", "pi_charge", "met_pt_corr_prerecoil", "met_phi_corr_prerecoil",
            "num_jets_pt30_corr", "jet_leadpt_corr", "jet_ht_pt20_corr"]
    with uproot.open(f"{CR}/{group}.root")["Events"] as t:
        ks = set(t.keys())
        a = t.arrays(cols + [c for c in ("cand_gen_mode", "cand_gen_common_hf_pdgId", "jet_bin", "tt_decay")
                             if c in ks], library="np")
    common.check_charge_config(a, False, group)
    a["rqx"] = a["mu_pt"] * np.cos(a["mu_phi"]) + a["xmu_pt_raw"] * np.cos(a["xmu_phi"])
    a["rqy"] = a["mu_pt"] * np.sin(a["mu_phi"]) + a["xmu_pt_raw"] * np.sin(a["xmu_phi"])
    if group == "DY_merged":
        with uproot.open(f"{CR}/{group}.root")["Events"] as t:
            g = t.arrays(list(msn.GEN_V_COLS) + ["met_pt_corr"], library="np")
        a["qx"], a["qy"] = g["gen_V_px"].astype(np.float64), g["gen_V_py"].astype(np.float64)
        a["vx"], a["vy"] = g["gen_Vvis_px"].astype(np.float64), g["gen_Vvis_py"].astype(np.float64)
        a["met_pt_skim"] = g["met_pt_corr"]   # recoil-corrected in the skim, checked in main()
    else:
        a["qx"], a["qy"], a["vx"], a["vy"] = a["rqx"], a["rqy"], a["rqx"], a["rqy"]
    sel = (np.abs(a["phi_mass"] - M_PHI) < 0.020) & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    a = {k: v[sel] for k, v in a.items()}
    if group == "Data_OS_merged":
        a["w"] = np.ones(len(a["event"]))
    else:
        a["w"] = a["weight"].astype(float) * common.br_reweight(a)
        if group == "DY_merged":
            a["w"] = a["w"] * common.dy_jetbin_sf(a)
    if not np.all(np.isfinite(a["w"])):
        sys.exit(f"{group}: non-finite weights")
    a["met_px"] = a["met_pt_corr_prerecoil"].astype(np.float64) * np.cos(a["met_phi_corr_prerecoil"])
    a["met_py"] = a["met_pt_corr_prerecoil"].astype(np.float64) * np.sin(a["met_phi_corr_prerecoil"])
    a["qt"], a["xpar"], a["uperp"] = msn.recoil_project(a["qx"], a["qy"], a["vx"], a["vy"], a["met_px"], a["met_py"])
    a["rqt"] = np.hypot(a["rqx"], a["rqy"])
    a["iq"] = np.clip(np.digitize(a["qt"], QT) - 1, 0, len(QT) - 2)
    a["inj"] = njet_class(a["num_jets_pt30_corr"])
    if (a["inj"] < 0).any():
        sys.exit(f"{group}: N_jet outside the classes")
    return a


def hist(x, w=None):
    return np.histogram(np.clip(x, UB[0], UB[-1] - 1e-9), UB, weights=w)[0].astype(float)


def cdf_of(h):
    """Normalised cumulative on the UB edges, forced monotonic (negative weights / data - TT)."""
    c = np.maximum.accumulate(np.clip(np.cumsum(h), 0.0, None))
    if c[-1] <= 0:
        sys.exit("empty histogram")
    return np.concatenate([[0.0], c / c[-1]])


def med_wid(h):
    """(median, (q84-q16)/2) of a histogram on UB."""
    c = cdf_of(h)
    keep = np.concatenate([[True], np.diff(c) > 0])
    q16, q50, q84 = np.interp([0.16, 0.5, 0.84], c[keep], UB[keep])
    return q50, 0.5 * (q84 - q16)


def derive(D, T, Y):
    """CDF maps for make_samples_nano.recoil_correct (Y = the DY events to derive on)."""
    nq, nk = len(QT) - 1, len(NJ)
    P = {f"{s}_{c}": np.zeros((nk, nq, len(UB))) for s in ("src", "tgt") for c in COMPS}
    P.update(fallback=np.zeros((nk, nq), bool), n_data=np.zeros((nk, nq), int))
    for iq in range(nq):
        for k in range(nk):
            n_d = int(((D["iq"] == iq) & (D["inj"] == k)).sum())
            use = k if n_d >= MIN_DATA else -1
            P["fallback"][k, iq], P["n_data"][k, iq] = use < 0, n_d
            sd, st, sy = ((a["iq"] == iq) & ((a["inj"] == use) if use >= 0 else True) for a in (D, T, Y))
            for c in COMPS:
                P[f"tgt_{c}"][k, iq] = cdf_of(hist(D[c][sd]) - hist(T[c][st], T["w"][st]))
                P[f"src_{c}"][k, iq] = cdf_of(hist(Y[c][sy], Y["w"][sy]))
    P.update(qt_edges=QT, u_edges=UB, nj_lo=np.array([lo for lo, _ in NJ]), nj_hi=np.array([hi for _, hi in NJ]))
    return P


def correct(P, Y):
    px, py = msn.recoil_correct(P, Y["num_jets_pt30_corr"], Y["qx"], Y["qy"], Y["vx"], Y["vy"], Y["met_px"], Y["met_py"])
    return np.hypot(px, py), np.arctan2(py, px)


def wid_w(x, w):
    return med_wid(hist(x, w))[1]


def reco_proj(a, pt, phi):
    """(x_par, u_perp) on the reco dimuon for every sample (display only)."""
    return msn.recoil_project(a["rqx"], a["rqy"], a["rqx"], a["rqy"], pt * np.cos(phi), pt * np.sin(phi))[1:]


def main():
    global plt, hep, common
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
    import common
    if common.SAMPLES_DIR != SAMPLES:
        sys.exit("SAMPLES differs from common.SAMPLES_DIR")
    S = {g: read(g) for g in GROUPS}
    final = {"Data": [], "TT": [], "DY": []}
    os.makedirs(OUT, exist_ok=True)
    for era, eid in ERAS:
        D, T, Y = ({k: v[a["era"] == eid] for k, v in a.items()} for a in (S["Data_OS_merged"], S["TT_merged"], S["DY_merged"]))
        lumi = common.ERAS[eid]["lumi_pb"] / 1000.0
        P = derive(D, T, Y)
        path = f"{OUT}/recoil_{era}.npz"
        if os.path.exists(path):   # the skim was made with this map: a different one must not replace it silently
            old = msn.recoil_load(path)
            same = set(old) == set(P) and all(np.array_equal(old[k], P[k]) for k in P)
            if not same:
                sys.exit(f"{era}: re-derived map differs from {path}, which the skim used; not overwritten "
                         "(re-run the skim if the new map is wanted)")
            print(f"{era}: re-derived map identical to {path}")
        else:
            np.savez(path, **P)
        dev = np.abs(correct(P, Y)[0] - Y["met_pt_skim"])
        if dev.max() > 1e-3:
            sys.exit(f"{era}: DY skim met_pt_corr differs from the map by up to {dev.max():.3g} GeV")
        out_rng = sum(int(((Y[c] <= UB[0]) | (Y[c] >= UB[-1])).sum()) for c in COMPS)
        print(f"{era}: DY skim met_pt_corr == map (max |d| {dev.max():.2g} GeV); DY components outside [{UB[0]:.0f}, {UB[-1]:.0f}] "
              f"(left unchanged): {out_rng}; fallback cells (N_jet x q_T):\n{P['fallback'].astype(int)}\n"
              f"data per cell:\n{P['n_data']}")

        # closure: maps from even DY events, shown on odd ones x 2
        even, odd = Y["event"] % 2 == 0, Y["event"] % 2 == 1
        Ye, Yo = ({k: x[m] for k, x in Y.items()} for m in (even, odd))
        Yo["w"] = 2.0 * Yo["w"]
        new = correct(derive(D, T, Ye), Yo)
        nr = (Yo["met_pt_corr_prerecoil"].astype(float), Yo["met_phi_corr_prerecoil"].astype(float))
        allD, allT, allY = (np.ones(len(a["w"]), bool) for a in (D, T, Yo))
        panels = [("MC without recoil", nr, allD, allT, allY), ("MC with recoil (CDF mapping)", new, allD, allT, allY)]
        title = (f"Z$\\to\\mu\\mu$ CR {era}: recoil calibration (maps from DY even events, DY shown = odd x 2); "
                 "input MET = Type-1 (JES+JER) + xy; MC = DY x BR rw x DY jet-bin SF + TT")
        for q, xlab in (("met", r"$p_T^{miss}$ [GeV]"), ("mt", r"$m_T$ [GeV]")):
            stack_figure(f"cr_zmumu_recoil_{q}_{era}", q, xlab, panels, D, T, Yo, lumi, title)
        split = []
        for lab, cut in (("reco $q_T<20$", lambda a: a["rqt"] < 20), ("reco $q_T\\geq20$", lambda a: a["rqt"] >= 20)):
            split += [(f"{lab}\nMC without recoil", nr, cut(D), cut(T), cut(Yo)),
                      (f"{lab}\nMC with recoil", new, cut(D), cut(T), cut(Yo))]
        stack_figure(f"cr_zmumu_recoil_met_qtsplit_{era}", "met", r"$p_T^{miss}$ [GeV]", split, D, T, Yo, lumi, title)
        uperp_figure(era, lumi, D, T, Yo, new)
        cells_figures(era, lumi, D, T, Yo, new)
        final["Data"].append(final_cols(D, D["met_pt_corr_prerecoil"], D["met_phi_corr_prerecoil"], D["w"]))
        final["TT"].append(final_cols(T, T["met_pt_corr_prerecoil"], T["met_phi_corr_prerecoil"], T["w"], T["tt_decay"]))
        final["DY"].append(final_cols(Yo, *new, Yo["w"]))
        final_figure(f"cr_zmumu_jme_recoil_final_{era}", {k: [x[-1]] for k, x in final.items()}, lumi, era)
    lumi_all = sum(common.ERAS[e]["lumi_pb"] for _, e in ERAS) / 1000.0
    final_figure("cr_zmumu_jme_recoil_final", final, lumi_all, None)


def uperp_figure(era, lumi, D, T, Yo, new):
    """u_perp width data/MC vs reco q_T per N_jet (everything binned on reco quantities, like for like)."""
    qc = [5, 15, 25, 40, 65, 100]
    up_new = reco_proj(Yo, *new)[1]
    up_old = reco_proj(Yo, Yo["met_pt_corr_prerecoil"], Yo["met_phi_corr_prerecoil"])[1]
    up_d = reco_proj(D, D["met_pt_corr_prerecoil"], D["met_phi_corr_prerecoil"])[1]
    up_t = reco_proj(T, T["met_pt_corr_prerecoil"], T["met_phi_corr_prerecoil"])[1]
    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    for k, lab in enumerate(NJ_LAB):
        r0, r1 = [], []
        for iq in range(len(QT) - 1):
            def m(a):
                return (a["rqt"] >= QT[iq]) & (a["rqt"] < QT[iq + 1]) & (a["inj"] == k)
            sd, st, sy = m(D), m(T), m(Yo)
            if sd.sum() < MIN_DATA:
                r0.append(np.nan); r1.append(np.nan)
                continue
            wd = med_wid(hist(up_d[sd]) - hist(up_t[st], T["w"][st]))[1]
            r0.append(wd / wid_w(up_old[sy], Yo["w"][sy]))
            r1.append(wd / wid_w(up_new[sy], Yo["w"][sy]))
        colr = ("#1f77b4", "#ff7f0e", "#2ca02c")[k]
        ax.plot(qc, r0, "o--", mfc="none", ms=6, color=colr, label=f"{lab}, MC without recoil")
        ax.plot(qc, r1, "o-", ms=6, color=colr, label=f"{lab}, MC with recoil")
    ax.axhline(1, color="0.4", ls="--", lw=1)
    ax.set_ylim(0.75, 1.15)
    ax.set_xlabel(r"reco $q_T$ [GeV] (bin centre; last bin $\geq$80)")
    ax.set_ylabel(r"(data$-t\bar{t}$)/DY width $(q_{84}-q_{16})/2$ of $u_\perp$")
    ax.legend(fontsize=10, ncol=2)
    hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{lumi:.1f}", com=13.6, fontsize=13)
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR {era}: $u_\\perp$ width, DY odd events (maps from even events)", fontsize=14)
    out = common.fig_path("sf", f"cr_zmumu_recoil_uperp_{era}.png")
    fig.savefig(out, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print("Saved", out)


def cells_figures(era, lumi, D, T, Yo, new):
    """Per-cell closure (JME POG: u_par, u_perp and |MET| across the parameterisation categories): data/MC with
    and without recoil in every reco q_T x N_jet cell. MC = DY odd x 2 + TT; all quantities on the reco dimuon."""
    vd = dict(zip(COMPS, reco_proj(D, D["met_pt_corr_prerecoil"], D["met_phi_corr_prerecoil"])), met=D["met_pt_corr_prerecoil"])
    vt = dict(zip(COMPS, reco_proj(T, T["met_pt_corr_prerecoil"], T["met_phi_corr_prerecoil"])), met=T["met_pt_corr_prerecoil"])
    vy0 = dict(zip(COMPS, reco_proj(Yo, Yo["met_pt_corr_prerecoil"], Yo["met_phi_corr_prerecoil"])), met=Yo["met_pt_corr_prerecoil"])
    vy1 = dict(zip(COMPS, reco_proj(Yo, *new)), met=new[0])
    for q, xlab, bins in (("xpar", r"$x_\parallel = u_\parallel + q_T$ [GeV]", np.arange(-80, 80 + 1e-9, 10.0)),
                          ("uperp", r"$u_\perp$ [GeV]", np.arange(-80, 80 + 1e-9, 10.0)),
                          ("met", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 10.0))):
        bc = 0.5 * (bins[1:] + bins[:-1])
        fig, axes = plt.subplots(len(NJ), len(QT) - 1, figsize=(5.2 * (len(QT) - 1), 4.2 * len(NJ)),
                                 layout="constrained", sharex=True, sharey=True)
        for k, nlab in enumerate(NJ_LAB):
            for iq in range(len(QT) - 1):
                ax = axes[k, iq]

                def m(a):
                    return (a["rqt"] >= QT[iq]) & (a["rqt"] < QT[iq + 1]) & (a["inj"] == k)
                sd, st, sy = m(D), m(T), m(Yo)
                dh = np.histogram(vd[q][sd], bins)[0]
                ht = np.histogram(vt[q][st], bins, weights=T["w"][st])[0]
                for vy, fmt, lab in ((vy0, "s", "without recoil"), (vy1, "o", "with recoil")):
                    mc = ht + np.histogram(vy[q][sy], bins, weights=Yo["w"][sy])[0]
                    with np.errstate(invalid="ignore", divide="ignore"):
                        r, re = np.where(mc > 0, dh / mc, np.nan), np.where(mc > 0, np.sqrt(dh) / mc, np.nan)
                    ax.errorbar(bc, r, yerr=re, fmt=fmt, ms=4, mfc="none" if fmt == "s" else None,
                                color="#1f77b4" if fmt == "s" else "black", label=f"data/MC {lab}")
                ax.axhline(1, color="red", lw=1)
                ax.set_ylim(0.3, 1.7)
                hi = "\\infty" if np.isinf(QT[iq + 1]) else f"{QT[iq + 1]:.0f}"
                ax.set_title(f"{nlab}, reco $q_T\\in[{QT[iq]:.0f},{hi})$, data N={sd.sum()}", fontsize=11)
                if k == len(NJ) - 1:
                    ax.set_xlabel(xlab, fontsize=13)
                if iq == 0:
                    ax.set_ylabel("data/MC", fontsize=13)
        axes[0, 0].legend(fontsize=10, loc="upper left")
        fig.suptitle(f"Z$\\to\\mu\\mu$ CR {era} ({lumi:.1f} fb$^{{-1}}$): per-cell closure of the recoil calibration "
                     "(maps from DY even events, MC = DY odd x 2 + TT)", fontsize=15)
        out = common.fig_path("sf", f"cr_zmumu_recoil_cells_{q}_{era}.png")
        fig.savefig(out, dpi=80, bbox_inches="tight")
        plt.close(fig)
        print("Saved", out)


def final_cols(a, met_pt, met_phi, w, tt_decay=None):
    """Plotted quantities with every correction (jets JES + JER; MET Type-1 + JER + xy + recoil) and, as
    "<q>_nr", the MET ones without the recoil step. Muon pt corrected."""
    def mt(pt, phi):
        return np.sqrt(2 * a["mu_pt_corr"] * pt * (1 - np.cos(a["mu_phi"] - phi)))
    return {"met": met_pt, "mt": mt(met_pt, met_phi), "metphi": met_phi, "njet": a["num_jets_pt30_corr"],
            "leadpt": a["jet_leadpt_corr"], "ht": a["jet_ht_pt20_corr"],
            "met_nr": a["met_pt_corr_prerecoil"], "mt_nr": mt(a["met_pt_corr_prerecoil"], a["met_phi_corr_prerecoil"]), "metphi_nr": a["met_phi_corr_prerecoil"],
            "w": w, "tt_decay": tt_decay if tt_decay is not None else np.zeros(len(w))}


FINAL_VARS = [("met", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 5.0)),
              ("mt", r"$m_T$ [GeV]", np.arange(0, 150 + 1e-9, 5.0)),
              ("metphi", r"$\phi(p_T^{miss})$", np.linspace(-np.pi, np.pi, 25)),
              ("njet", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 6.5, 1.0)),
              ("leadpt", r"leading jet $p_T$ [GeV]", np.arange(15, 255 + 1e-9, 10.0)),
              ("ht", r"$H_T$ (jets $p_T>20$) [GeV]", np.arange(0, 400 + 1e-9, 20.0))]


def final_figure(fname, parts, lumi, era):
    """Z->mumu CR data/MC, stack style, all corrections applied; DY = odd events x 2 (parameters from even).
    The data/MC panel of the MET variables also shows the same events without the recoil step."""
    P = {k: {c: np.concatenate([p[c] for p in v]) for c in v[0]} for k, v in parts.items()}
    nrow = (len(FINAL_VARS) + 2) // 3
    fig = plt.figure(figsize=(24, 7.5 * nrow), layout="constrained")
    outer = fig.add_gridspec(nrow, 3)
    for i, (q, xlab, bins) in enumerate(FINAL_VARS):
        bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
        inner = outer[i // 3, i % 3].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        comps = [(r"$t\bar{t}$ LNu2Q", P["TT"], P["TT"]["tt_decay"] == 1, "#b39ddb"),
                 (r"$t\bar{t}$ 2L2Nu", P["TT"], P["TT"]["tt_decay"] == 0, "#7f5fbf"),
                 (r"DY$\to\ell\ell$", P["DY"], np.ones(len(P["DY"]["w"]), bool), "#4a90c2")]
        bottom, tot, tot_e2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for clab, a, m, colr in comps:
            h, _ = np.histogram(a[q][m], bins, weights=a["w"][m])
            e2, _ = np.histogram(a[q][m], bins, weights=a["w"][m] ** 2)
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none", label=f"{clab}  N={h.sum():,.0f}")
            bottom += h; tot += h; tot_e2 += e2
        terr = np.sqrt(tot_e2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////", edgecolor="0.35",
               linewidth=0.0, label="MC stat")
        dh, _ = np.histogram(P["Data"][q], bins)
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black", label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 1.45 * max(tot.max(), dh.max(), 1))
        ax.legend(fontsize=10, loc="upper right")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Events / bin", fontsize=16)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        if q + "_nr" in P["DY"]:
            tot0 = sum(np.histogram(a[q + "_nr"][m], bins, weights=a["w"][m])[0] for _, a, m, _ in comps)
            dh0, _ = np.histogram(P["Data"][q + "_nr"], bins)
            with np.errstate(invalid="ignore", divide="ignore"):
                r0, re0 = np.where(tot0 > 0, dh0 / tot0, np.nan), np.where(tot0 > 0, np.sqrt(dh0) / tot0, np.nan)
            rax.errorbar(bc, r0, yerr=re0, fmt="s", ms=5, mfc="none", mew=1.3, color="#1f77b4", label="without recoil")
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black", label="all corrections incl. recoil")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.3, 1.7)
        if i == 0:   # data/MC legend goes into the stack panel, below the stack legend (not over the points)
            leg = ax.get_legend()
            h_, l_ = rax.get_legend_handles_labels()
            ax.legend(h_, [f"data/MC: {x}" for x in l_], fontsize=11, loc="center right")
            ax.add_artist(leg)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if i == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{lumi:.1f}", com=13.6, fontsize=14)
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR{' ' + era if era else ''}. Stack: all corrections (jets JES + JER; MET Type-1 + JER + xy "
                 "+ recoil). MC = DY (odd events x 2) x BR rw x DY jet-bin SF + TT", fontsize=15)
    out = common.fig_path("sf", f"{fname}.png")
    fig.savefig(out, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print("Saved", out)


def stack_figure(fname, q, xlab, panels, D, T, Yo, lumi, title):
    """One column per panel (label, DY (pt, phi), masks for data / TT / DY): stack + data, data/MC below.
    Data and TT always use the MET without recoil (they are not corrected)."""
    bins = np.arange(0, 150 + 1e-9, 5.0)
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)

    def var(a, pt, phi):
        return pt if q == "met" else np.sqrt(2 * a["mu_pt_corr"] * pt * (1 - np.cos(a["mu_phi"] - phi)))
    fig = plt.figure(figsize=(7.2 * len(panels), 7.6), layout="constrained")
    outer = fig.add_gridspec(1, len(panels))
    xT = var(T, T["met_pt_corr_prerecoil"], T["met_phi_corr_prerecoil"])
    xD = var(D, D["met_pt_corr_prerecoil"], D["met_phi_corr_prerecoil"])
    for i, (lab, (ypt, yphi), mD, mT, mY) in enumerate(panels):
        inner = outer[i].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        comps = [(r"$t\bar{t}$ LNu2Q", xT, T["w"], mT & (T["tt_decay"] == 1), "#b39ddb"),
                 (r"$t\bar{t}$ 2L2Nu", xT, T["w"], mT & (T["tt_decay"] == 0), "#7f5fbf"),
                 (r"DY$\to\ell\ell$", var(Yo, ypt, yphi), Yo["w"], mY, "#4a90c2")]
        bottom, tot, tot_e2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for clab, x, w, m, colr in comps:
            h, _ = np.histogram(x[m], bins, weights=w[m])
            e2, _ = np.histogram(x[m], bins, weights=w[m] ** 2)
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none", label=f"{clab}  N={h.sum():,.0f}")
            bottom += h; tot += h; tot_e2 += e2
        terr = np.sqrt(tot_e2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////", edgecolor="0.35",
               linewidth=0.0, label="MC stat")
        dh, _ = np.histogram(xD[mD], bins)
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black", label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 1.45 * max(tot.max(), dh.max(), 1))
        ax.legend(fontsize=9, loc="upper right")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Events / bin", fontsize=16)
        ax.text(0.03, 0.95, lab, transform=ax.transAxes, fontsize=14, va="top", weight="bold")
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.3, 1.7)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if i == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{lumi:.1f}", com=13.6, fontsize=14)
        hi = bins[:-1] >= 60
        print(f"  {fname} | {lab.replace(chr(10), ' / ')}: data/MC {dh.sum() / tot.sum():.3f}, "
              f"60-150 GeV {dh[hi].sum() / tot[hi].sum():.3f}")
    fig.suptitle(title, fontsize=15)
    out = common.fig_path("sf", f"{fname}.png")
    fig.savefig(out, dpi=90, bbox_inches="tight")
    plt.close(fig)
    print("Saved", out)

if __name__ == "__main__":
    a = sys.argv[1:]
    if a == []:
        main()
    else:
        sys.exit(__doc__)

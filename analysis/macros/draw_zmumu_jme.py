"""Z->mumu CR data/MC for the JME steps of jme_zmumu_study.py (2026-11-06).

Same selection as draw_cr_zmumu.py (zmumu skim, phi20, 76 < M(mumu) < 106, OS data), same MC
weight: skim weight x common.br_reweight x common.dy_jetbin_sf (DY). MC = DY + TT (W+jets is
< 1% of the CR in every MET bin and is not drawn). The event set is identical for every step;
only the plotted jet/MET quantity changes (step definitions: jme_zmumu_study.py docstring).

Figures (MVA/figures/sf/cr_zmumu_jme_*.png):
  <q>          stack + data/MC per step, q = met, mt, metphi (MET steps) and njet, leadpt, ht (jet steps)
  ratios       data/MC per step overlaid, all quantities
  met_nano_vs_full   MET stack + data/MC, nano vs full JME chain (S5), two panels for slides
  shift        data and MC separately: each step's MET / m_T histogram divided by step 0's
  recoil       u_par response -<u_par>/q_T and widths sigma(u_par), sigma(u_perp) vs q_T,
               data and MC, per MET step; u_T = -(MET + q_T), q_T from the nano muon pts
               (the MET carries no muon correction)
  pxpy         <p_x^miss>, <p_y^miss> vs N_PV, data and MC, for m0, m0xy, m4, m5 (the MET xy check)
  recoil_njet  per jet multiplicity (j2 N_jet = 0, 1, 2, >=3): robust u_perp width (q84-q16)/2 and
               response vs q_T, data/MC for S0, S1, S4 (is the MC recoil off independent of DY jet modelling?)
  met_qtnj_reweighted  MET data/MC and u_perp width with MC as is vs MC reweighted to data in (q_T, N_jet)
  met_jer_options / uperp_jer_options  MET stack and u_perp width data/MC for: JER not in MET (S1),
               JER scaling only (m4match, JME known-issue option b), JER in MET (S5, current default)
Printed: mean MET / m_T per step for data and MC, and the m4 diagnostics (MC).

    python3 draw_zmumu_jme.py [--era=2022|2022EE]
"""
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
import common  # noqa: E402

ARGS = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
ERA = ARGS.get("era")
ERA_ID = {"2022": 1, "2022EE": 2}.get(ERA) if ERA else None
if ERA and ERA_ID is None:
    sys.exit("--era must be 2022 or 2022EE")
CR = f"{common.SAMPLES_DIR}/CR/zmumu"
JME = f"{common.SAMPLES_DIR}/CR/zmumu_jme_steps"
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
LUMI_FB = ((common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) if ERA_ID is None
           else common.ERAS[ERA_ID]["lumi_pb"]) / 1000.0
TAG = f"_{ERA}" if ERA else ""
WCHAIN = "MC x BR rw x DY jet-bin SF, no DNN cut, same events in every step"

BKG = [  # as draw_cr_zmumu.py, without W+jets
    (r"$t\bar{t}$ LNu2Q", "TT_merged", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT_merged", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
    (r"DY$\to\ell\ell$", "DY_merged", None, "#4a90c2"),
]
GROUPS = ("Data_OS_merged", "TT_merged", "DY_merged")
MC_GROUPS = ("TT_merged", "DY_merged")
MET_STEPS = [("S0 nano", "m0"), ("S1 JES", "m1"), ("S3 +JER AK4", "m3"),
             ("S4 +JER CorrT1", "m4"), ("S5 +xy", "m5")]
JET_STEPS = [("S0 nano", "j0"), ("S1 JES", "j1"), ("S2 +JER", "j2")]
EXTRA_MET = [("S0+xy", "m0xy"), ("S4 matched only", "m4match"), ("S4 thr. on JES pt", "m4nocross"),
             ("S4 jet factor", "m4jf")]
QUANT = [  # (name, x label, bins, steps, column pattern)
    ("met", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 5.0), MET_STEPS, "{}_pt"),
    ("mt", r"$m_T$ [GeV]", np.arange(0, 150 + 1e-9, 5.0), MET_STEPS, "{}_mt"),
    ("metphi", r"$\phi(p_T^{miss})$", np.linspace(-np.pi, np.pi, 25), MET_STEPS, "{}_phi"),
    ("njet", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 6.5, 1.0), JET_STEPS, "{}_num_jets_pt30"),
    ("leadpt", r"leading jet $p_T$ [GeV]", np.arange(15, 255 + 1e-9, 10.0), JET_STEPS, "{}_jet_leadpt"),
    ("ht", r"$H_T$ (jets $p_T>20$) [GeV]", np.arange(0, 400 + 1e-9, 20.0), JET_STEPS, "{}_jet_ht_pt20"),
]
QT_BINS = np.array([0, 10, 20, 30, 40, 60, 80, 120, 200.0])
NPV_BINS = np.arange(0, 60 + 1e-9, 5.0)
MIN_EVENTS = 100   # pxpy: bins with fewer (unweighted) events are not drawn


def read(group):
    base = ["run", "luminosityBlock", "event", "weight", "phi_mass", "m_mumu", "era",
            "mu_pt", "mu_phi", "mu_pt_corr", "xmu_pt_raw", "xmu_phi", "mu_charge", "pi_charge"]
    with uproot.open(f"{CR}/{group}.root")["Events"] as t:
        keys = set(t.keys())
        want = base + [c for c in ("tt_decay", "jet_bin", "cand_gen_mode", "cand_gen_common_hf_pdgId")
                       if c in keys]
        a = t.arrays(want, library="np")
    common.check_charge_config(a, False, group)
    with uproot.open(f"{JME}/{group}.root")["Events"] as t:
        j = t.arrays(library="np")
    for k in ("run", "luminosityBlock", "event"):   # row alignment, fail loud
        if not np.array_equal(a[k], j[k]):
            sys.exit(f"{group}: JME steps friend not row-aligned on {k}")
    a.update({k: v.astype(np.float64) if v.dtype == np.float32 else v
              for k, v in j.items() if k not in ("run", "luminosityBlock", "event")})
    sel = (np.abs(a["phi_mass"] - M_PHI) < 0.020) & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    if ERA_ID is not None:
        sel &= a["era"] == ERA_ID
    a = {k: v[sel] for k, v in a.items()}
    if group == "Data_OS_merged":
        a["w"] = np.ones(len(a["event"]))
    else:
        w = a["weight"].astype(float) * common.br_reweight(a)
        a["w"] = w * common.dy_jetbin_sf(a) if group == "DY_merged" else w
    if not np.all(np.isfinite(a["w"])):
        sys.exit(f"{group}: non-finite weights in the Z window")
    # derived per MET step: m_T (corrected best muon, as W_tmass_corr), p_x, p_y, recoil
    qx = a["mu_pt"] * np.cos(a["mu_phi"]) + a["xmu_pt_raw"] * np.cos(a["xmu_phi"])
    qy = a["mu_pt"] * np.sin(a["mu_phi"]) + a["xmu_pt_raw"] * np.sin(a["xmu_phi"])
    a["qt"] = np.hypot(qx, qy)
    ux_hat, uy_hat = qx / a["qt"], qy / a["qt"]
    for _l, s in MET_STEPS + EXTRA_MET:
        pt, ph = a[f"{s}_pt"], a[f"{s}_phi"]
        a[f"{s}_mt"] = np.sqrt(2.0 * a["mu_pt_corr"] * pt * (1.0 - np.cos(a["mu_phi"] - ph)))
        a[f"{s}_px"], a[f"{s}_py"] = pt * np.cos(ph), pt * np.sin(ph)
        ux, uy = -(a[f"{s}_px"] + qx), -(a[f"{s}_py"] + qy)
        a[f"{s}_upar"] = ux * ux_hat + uy * uy_hat
        a[f"{s}_uperp"] = -ux * uy_hat + uy * ux_hat
    return a


def hist(x, w, bins):
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w * w)
    return h, e2


def label_title(ax):
    hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6, fontsize=14)


def save(fig, name):
    out = common.fig_path("sf", f"cr_zmumu_jme_{name}{TAG}.png")
    fig.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


S = {g: read(g) for g in GROUPS}
d = S["Data_OS_merged"]
MC = {k: np.concatenate([S[g][k] for g in MC_GROUPS]) for k in S["DY_merged"] if k in S["TT_merged"]}
print(f"selected (phi20, Z window{', ' + ERA if ERA else ''}): data {len(d['event'])}, "
      + ", ".join(f"{g} {S[g]['w'].sum():.1f}" for g in MC_GROUPS))

# ---- stack + data/MC per step; collect ratios ----
def stack_figure(name, xlab, bins, steps, pat, title):
    """One column per step: stack + data, data/MC below. Returns [(label, ratio, ratio error)]."""
    ncol = len(steps)
    fig = plt.figure(figsize=(7.2 * ncol, 7.6))
    outer = gridspec.GridSpec(1, ncol, figure=fig, wspace=0.32, top=0.86, bottom=0.12)
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    out = []
    for i, (vlab, s) in enumerate(steps):
        col = pat.format(s)
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i], height_ratios=[3, 1], hspace=0.06)
        ax, rax = (fig.add_subplot(inner[k]) for k in range(2))
        bottom, tot, tot_e2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for lab, g, extra, colr in BKG:
            a = S[g]
            m = extra(a) if extra else np.ones(len(a["w"]), bool)
            h, e2 = hist(a[col][m], a["w"][m], bins)
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none", label=f"{lab}  N={h.sum():,.0f}")
            bottom += h; tot += h; tot_e2 += e2
        terr = np.sqrt(tot_e2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////",
               edgecolor="0.35", linewidth=0.0, label="MC stat")
        dh, _ = np.histogram(d[col], bins)
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black", label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 1.45 * max(tot.max(), dh.max(), 1))
        ax.legend(fontsize=9, loc="upper right")
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Events / bin", fontsize=16)
        ax.text(0.03, 0.95, vlab, transform=ax.transAxes, fontsize=16, va="top", weight="bold")
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
            rerr = np.sqrt(dh / tot ** 2 + dh ** 2 * tot_e2 / tot ** 4)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.3, 1.7)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if i == 0:
            label_title(ax)
        out.append((vlab, r, rerr))
    fig.suptitle(title, fontsize=15, y=0.95)
    save(fig, name)
    return out


ERA_TXT = f" ({ERA})" if ERA else ""
summary = {name: stack_figure(name, xlab, bins, steps, pat, f"Z$\\to\\mu\\mu$ CR, JME steps{ERA_TXT}: {WCHAIN}")
           for name, xlab, bins, steps, pat in QUANT}
# for slides: MET before any correction vs after the full chain (the analysis default)
stack_figure("met_nano_vs_full", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 5.0),
             [("before: nano PuppiMET", "m0"), ("after: JES + JER + xy", "m5")], "{}_pt",
             f"Z$\\to\\mu\\mu$ CR{ERA_TXT}: MC x BR rw x DY jet-bin SF, same events in both")


# ---- data/MC per step overlaid ----
fig, axes = plt.subplots(2, 3, figsize=(24, 13))
colors = ["black", "#d62728", "#1f77b4", "#ff7f0e", "#2ca02c"]
for ax, (name, xlab, bins, _s, _p) in zip(axes.flat, QUANT):
    bc = 0.5 * (bins[1:] + bins[:-1])
    off = np.diff(bins)[0] * 0.08
    for k, (vlab, r, re) in enumerate(summary[name]):
        ax.errorbar(bc + (k - 2) * off, r, yerr=re, fmt="o-", ms=3.5, lw=1, color=colors[k], label=vlab)
    ax.axhline(1, color="0.5", lw=1, ls="--")
    ax.set_ylim(0.3, 1.7)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(xlab)
    ax.set_ylabel("data/MC")
    ax.legend(fontsize=12)
label_title(axes.flat[0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR data/MC per JME step{' (' + ERA + ')' if ERA else ''} ({WCHAIN})", fontsize=17)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "ratios")

# ---- shift: each step / step 0, data and MC separately ----
fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex="col")
for c, (q, xlab) in enumerate((("pt", r"$p_T^{miss}$ [GeV]"), ("mt", r"$m_T$ [GeV]"))):
    bins = np.arange(0, 150 + 1e-9, 5.0)
    bc = 0.5 * (bins[1:] + bins[:-1])
    for r, (lab, a) in enumerate((("OS data", d), ("MC (DY+TT)", MC))):
        ax = axes[r, c]
        h0, _ = np.histogram(a[f"m0_{q}"], bins, weights=a["w"])
        for k, (vlab, s) in enumerate(MET_STEPS[1:] + EXTRA_MET[:1]):
            h, _ = np.histogram(a[f"{s}_{q}"], bins, weights=a["w"])
            with np.errstate(invalid="ignore", divide="ignore"):
                ax.plot(bc, np.where(h0 > 0, h / h0, np.nan), "o-", ms=3, lw=1,
                        color=(colors[1:] + ["0.5"])[k], label=vlab)
        ax.axhline(1, color="0.5", lw=1, ls="--")
        ax.set_ylim(0.3, 2.5)
        ax.set_ylabel(f"{lab}: step / S0")
        ax.legend(fontsize=11)
        if r == 1:
            ax.set_xlabel(xlab)
label_title(axes[0, 0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR: change of each JME step w.r.t. nano, data and MC separately"
             f"{' (' + ERA + ')' if ERA else ''}", fontsize=16)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "shift")


# ---- recoil: response and resolution vs q_T ----
def wstats(x, w, sel):
    ww, xx = w[sel], x[sel]
    if ww.sum() <= 0:
        return np.nan, np.nan
    mu = np.average(xx, weights=ww)
    return mu, np.sqrt(np.average((xx - mu) ** 2, weights=ww))


qc = 0.5 * (QT_BINS[1:] + QT_BINS[:-1])
fig, axes = plt.subplots(2, 3, figsize=(24, 12), sharex=True)
for k, (vlab, s) in enumerate(MET_STEPS):
    res = {}
    for lab, a in (("data", d), ("MC", MC)):
        rows = []
        for lo, hi in zip(QT_BINS[:-1], QT_BINS[1:]):
            sel = (a["qt"] >= lo) & (a["qt"] < hi)
            mpar, spar = wstats(a[f"{s}_upar"], a["w"], sel)
            _m, sperp = wstats(a[f"{s}_uperp"], a["w"], sel)
            qm = np.average(a["qt"][sel], weights=a["w"][sel]) if a["w"][sel].sum() > 0 else np.nan
            rows.append((-mpar / qm, spar, sperp))
        res[lab] = np.array(rows)
    for c, (ttl, ylab) in enumerate(((r"response $-\langle u_\parallel\rangle/q_T$", "response"),
                                     (r"$\sigma(u_\parallel)$ [GeV]", r"$\sigma(u_\parallel)$"),
                                     (r"$\sigma(u_\perp)$ [GeV]", r"$\sigma(u_\perp)$"))):
        axes[0, c].plot(qc, res["MC"][:, c], "o-", ms=3.5, lw=1, color=colors[k], label=f"MC {vlab}")
        if k in (0, 1):   # data only changes at S1 (JES) and S5 (xy, ~no effect on widths)
            axes[0, c].plot(qc, res["data"][:, c], "s--", ms=5, lw=1, mfc="none", color=colors[k],
                            label=f"data {vlab}")
        axes[1, c].plot(qc, res["data"][:, c] / res["MC"][:, c], "o-", ms=3.5, lw=1, color=colors[k], label=vlab)
        axes[0, c].set_title(ttl, fontsize=16)
        axes[0, c].set_ylabel(ylab)
        axes[1, c].set_ylabel("data/MC")
        axes[1, c].axhline(1, color="0.5", lw=1, ls="--")
        axes[1, c].set_xlabel(r"$q_T$ [GeV]")
for ax in axes.flat:
    ax.legend(fontsize=10)
label_title(axes[0, 0])
fig.suptitle(r"Z$\to\mu\mu$ CR hadronic recoil $u_T=-(p_T^{miss}+q_T)$, per JME step "
             f"{'(' + ERA + ') ' if ERA else ''}(weighted mean / RMS per $q_T$ bin; {WCHAIN})", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.95))
save(fig, "recoil")

# ---- MET xy check: <p_x>, <p_y> vs N_PV ----
nc = 0.5 * (NPV_BINS[1:] + NPV_BINS[:-1])
fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex=True)
for r, (lab, a) in enumerate((("OS data", d), ("MC (DY+TT)", MC))):
    for c, comp in enumerate(("px", "py")):
        ax = axes[r, c]
        for k, (vlab, s) in enumerate((("S0 nano", "m0"), ("S0+xy", "m0xy"), ("S4", "m4"), ("S5 = S4+xy", "m5"))):
            mean = []
            for lo, hi in zip(NPV_BINS[:-1], NPV_BINS[1:]):
                sel = (a["npv"] >= lo) & (a["npv"] < hi)
                mean.append(np.average(a[f"{s}_{comp}"][sel], weights=a["w"][sel])
                            if sel.sum() >= MIN_EVENTS else np.nan)
            ax.plot(nc, mean, "o-" if k % 2 == 0 else "s--", ms=4, lw=1, color=["black", "0.5", "#ff7f0e", "#2ca02c"][k],
                    label=vlab)
        ax.axhline(0, color="0.5", lw=1, ls=":")
        ax.set_ylabel(rf"{lab}: $\langle p_{{{comp[1]}}}^{{miss}}\rangle$ [GeV]")
        ax.legend(fontsize=11)
        if r == 1:
            ax.set_xlabel(r"$N_{PV}$ (good)")
label_title(axes[0, 0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR MET xy check{' (' + ERA + ')' if ERA else ''} (bins with < {MIN_EVENTS} events not drawn)",
             fontsize=16)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "pxpy")


# ---- is it the DY modelling? recoil per N_jet, and MC reweighted to data in (q_T, N_jet) ----
def wquant(x, w, qs):
    """Weighted quantiles (cumulative weight; negative MC weights included as they are)."""
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return np.interp(np.asarray(qs) * c[-1], c, x[o])


def half_width(x, w):
    lo, hi = wquant(x, w, (0.16, 0.84))
    return 0.5 * (hi - lo)


NJ_COL = "j2_num_jets_pt30"   # JES + JER jets (mandatory corrections, = the analysis default)
NJ_CATS = [("N$_{jet}$=0", 0, 0), ("N$_{jet}$=1", 1, 1), ("N$_{jet}$=2", 2, 2), (r"N$_{jet}\geq$3", 3, 99)]
RECOIL_STEPS = [("S0 nano", "m0"), ("S1 JES", "m1"), ("S4 +JER", "m4")]
QT_COARSE = np.array([0, 10, 20, 30, 50, 80, 200.0])
qcc = 0.5 * (QT_COARSE[1:] + QT_COARSE[:-1])
fig, axes = plt.subplots(3, len(NJ_CATS), figsize=(7 * len(NJ_CATS), 15), sharex=True)
print(f"\nrecoil per N_jet ({NJ_COL}); robust width (q84-q16)/2 of u_perp, data/MC:")
for c, (njlab, lo_nj, hi_nj) in enumerate(NJ_CATS):
    cat = {lab: (a[NJ_COL] >= lo_nj) & (a[NJ_COL] <= hi_nj) for lab, a in (("data", d), ("MC", MC))}
    print(f"  {njlab:14s} data {int(cat['data'].sum()):6d}  MC {MC['w'][cat['MC']].sum():9.1f}")
    for k, (vlab, s) in enumerate(RECOIL_STEPS):
        wid = {"data": [], "MC": []}
        resp = {"data": [], "MC": []}
        for lab, a in (("data", d), ("MC", MC)):
            for lo, hi in zip(QT_COARSE[:-1], QT_COARSE[1:]):
                sel = cat[lab] & (a["qt"] >= lo) & (a["qt"] < hi)
                if sel.sum() < MIN_EVENTS:
                    wid[lab].append(np.nan); resp[lab].append(np.nan)
                    continue
                wid[lab].append(half_width(a[f"{s}_uperp"][sel], a["w"][sel]))
                resp[lab].append(-np.average(a[f"{s}_upar"][sel], weights=a["w"][sel])
                                 / np.average(a["qt"][sel], weights=a["w"][sel]))
        wd_, wm = np.array(wid["data"]), np.array(wid["MC"])
        axes[0, c].plot(qcc, wm, "o-", ms=4, lw=1, color=colors[k], label=f"MC {vlab}")
        if k < 2:
            axes[0, c].plot(qcc, wd_, "s--", ms=5, lw=1, mfc="none", color=colors[k], label=f"data {vlab}")
        axes[1, c].plot(qcc, wd_ / wm, "o-", ms=4, lw=1, color=colors[k], label=vlab)
        axes[2, c].plot(qcc, np.array(resp["data"]) / np.array(resp["MC"]), "o-", ms=4, lw=1,
                        color=colors[k], label=vlab)
        print(f"    {vlab:9s} " + " ".join(f"{x:.3f}" for x in wd_ / wm) + f"   (q_T bins {QT_COARSE.tolist()})")
    axes[0, c].set_title(njlab, fontsize=18)
    axes[0, c].set_ylabel(r"$(q_{84}-q_{16})/2$ of $u_\perp$ [GeV]")
    axes[1, c].set_ylabel(r"data/MC width of $u_\perp$")
    axes[2, c].set_ylabel(r"data/MC response $-\langle u_\parallel\rangle/q_T$")
    axes[2, c].set_xlabel(r"$q_T$ [GeV]")
    for r in (1, 2):
        axes[r, c].axhline(1, color="0.5", lw=1, ls="--")
        axes[r, c].set_ylim(0.75, 1.15)
for ax in axes.flat:
    ax.legend(fontsize=10)
label_title(axes[0, 0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR recoil per jet multiplicity ({NJ_COL}){ERA_TXT}; {WCHAIN}; "
             f"q_T bins with < {MIN_EVENTS} events not drawn (response at q_T<10 GeV is resolution-dominated)",
             fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "recoil_njet")

# MC reweighted to data in (q_T, N_jet): removes the DY q_T / jet-multiplicity mismodelling
QT_RW = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 1e9])
NJ_RW = np.array([-0.5, 0.5, 1.5, 2.5, 3.5, 99])
hd2, _, _ = np.histogram2d(d["qt"], d[NJ_COL], (QT_RW, NJ_RW))
hm2, _, _ = np.histogram2d(MC["qt"], MC[NJ_COL], (QT_RW, NJ_RW), weights=MC["w"])
empty = hm2 <= 0
if (empty & (hd2 > 0)).any():
    sys.exit(f"(q_T, N_jet) reweighting: {int((empty & (hd2 > 0)).sum())} cells with data but no MC")
f2 = np.where(empty, 1.0, hd2 / np.where(empty, 1.0, hm2))   # cells empty in both: factor 1
iq = np.clip(np.digitize(MC["qt"], QT_RW) - 1, 0, len(QT_RW) - 2)
inj = np.clip(np.digitize(MC[NJ_COL], NJ_RW) - 1, 0, len(NJ_RW) - 2)
w_rw = MC["w"] * f2[iq, inj]
hm2_rw, _, _ = np.histogram2d(MC["qt"], MC[NJ_COL], (QT_RW, NJ_RW), weights=w_rw)
if not np.allclose(hm2_rw[~empty], hd2[~empty], rtol=1e-9):
    sys.exit("(q_T, N_jet) reweighting does not close")
print(f"\n(q_T, N_jet) reweighting factors: min {f2.min():.3f} max {f2.max():.3f}; MC sum "
      f"{MC['w'].sum():.1f} -> {w_rw.sum():.1f} (data {len(d['w'])})")
bins = np.arange(0, 150 + 1e-9, 5.0)
bc = 0.5 * (bins[1:] + bins[:-1])
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
for k, (vlab, s) in enumerate([("S0 nano", "m0"), ("S1 JES", "m1"), ("S5 = default", "m5")]):
    hd, _ = np.histogram(d[f"{s}_pt"], bins)
    for ax, w_, ttl in ((axes[0], MC["w"], "MC as is"), (axes[1], w_rw, r"MC reweighted to data in ($q_T$, N$_{jet}$)")):
        hm, _ = np.histogram(MC[f"{s}_pt"], bins, weights=w_)
        with np.errstate(invalid="ignore", divide="ignore"):
            ax.errorbar(bc, hd / hm, yerr=np.sqrt(hd) / hm, fmt="o-", ms=3.5, lw=1, color=colors[[0, 1, 4][k]], label=vlab)
        ax.set_title(ttl, fontsize=16)
    wid_d, wid_m, wid_rw = [], [], []
    for lo, hi in zip(QT_COARSE[:-1], QT_COARSE[1:]):
        sd, sm = (d["qt"] >= lo) & (d["qt"] < hi), (MC["qt"] >= lo) & (MC["qt"] < hi)
        wid_d.append(half_width(d[f"{s}_uperp"][sd], d["w"][sd]))
        wid_m.append(half_width(MC[f"{s}_uperp"][sm], MC["w"][sm]))
        wid_rw.append(half_width(MC[f"{s}_uperp"][sm], w_rw[sm]))
    axes[2].plot(qcc, np.array(wid_d) / np.array(wid_m), "o--", ms=4, lw=1, mfc="none", color=colors[[0, 1, 4][k]],
                 label=f"{vlab}, MC as is")
    axes[2].plot(qcc, np.array(wid_d) / np.array(wid_rw), "o-", ms=4, lw=1, color=colors[[0, 1, 4][k]],
                 label=f"{vlab}, MC reweighted")
    print(f"  {vlab:13s} u_perp width data/MC  as is " + " ".join(f"{x:.3f}" for x in np.array(wid_d) / np.array(wid_m))
          + "  reweighted " + " ".join(f"{x:.3f}" for x in np.array(wid_d) / np.array(wid_rw)))
for ax in axes[:2]:
    ax.set_xlabel(r"$p_T^{miss}$ [GeV]")
    ax.set_ylabel("data/MC")
    ax.set_ylim(0.3, 1.7)
for ax in axes:
    ax.axhline(1, color="0.5", lw=1, ls="--")
    ax.legend(fontsize=11)
axes[2].set_xlabel(r"$q_T$ [GeV]")
axes[2].set_ylabel(r"data/MC width of $u_\perp$")
axes[2].set_ylim(0.75, 1.15)
label_title(axes[0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR: MET data/MC before / after reweighting MC to data in ($q_T$, {NJ_COL}){ERA_TXT}",
             fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.94))
save(fig, "met_qtnj_reweighted")

# ---- MET options for JER (2026-11-09): not propagated / scaling only (JME known-issue option b) / default ----
JER_OPTIONS = [("Type-1 MET, JES only\n(JER not propagated)", "m1"),
               ("JER propagated only for\ngen-matched jets", "m4match"),
               ("JER propagated for\nall jets", "m4")]   # all three without MET xy: only the JER treatment differs
stack_figure("met_jer_options", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 5.0), JER_OPTIONS, "{}_pt",
             f"Z$\\to\\mu\\mu$ CR{ERA_TXT}: MET with three treatments of the MC JER smearing (Type-1, no MET xy in any); "
             "MC x BR rw x DY jet-bin SF, same events in all")
fig, ax = plt.subplots(figsize=(10, 7))
for k, (vlab, s) in enumerate(JER_OPTIONS):
    r_ = []
    for lo, hi in zip(QT_COARSE[:-1], QT_COARSE[1:]):
        sd, sm = (d["qt"] >= lo) & (d["qt"] < hi), (MC["qt"] >= lo) & (MC["qt"] < hi)
        r_.append(half_width(d[f"{s}_uperp"][sd], d["w"][sd]) / half_width(MC[f"{s}_uperp"][sm], MC["w"][sm]))
    ax.plot(qcc, r_, "o-", ms=5, lw=1.5, color=["#1f77b4", "#ff7f0e", "#2ca02c"][k], label=vlab.replace("\n", " "))
ax.axhline(1, color="0.5", lw=1, ls="--")
ax.set_ylim(0.75, 1.05)
ax.set_xlabel(r"$q_T$ [GeV]")
ax.set_ylabel(r"data/MC width $(q_{84}-q_{16})/2$ of $u_\perp$", fontsize=18)
ax.legend(fontsize=13, loc="lower right")
label_title(ax)
ax.set_title(f"Z$\\to\\mu\\mu$ CR{ERA_TXT}: MET resolution data/MC per JER treatment", fontsize=14, pad=28)
save(fig, "uperp_jer_options")

# ---- printed summary ----
print("\nmean MET / m_T [GeV] per step (Z window, weighted):")
print(f"  {'step':20s} {'data MET':>9s} {'MC MET':>9s} {'data mT':>9s} {'MC mT':>9s}")
for vlab, s in MET_STEPS + EXTRA_MET:
    print(f"  {vlab:20s} {np.average(d[s + '_pt']):9.3f} {np.average(MC[s + '_pt'], weights=MC['w']):9.3f}"
          f" {np.average(d[s + '_mt']):9.3f} {np.average(MC[s + '_mt'], weights=MC['w']):9.3f}")
dy = S["DY_merged"]
print("\nDY: MET jets crossing 15 GeV per event (weighted mean): "
      + ", ".join(f"{k} {np.average(dy[k], weights=dy['w']):.3f}"
                  for k in ("ncross_up_ct1", "ncross_up_ak4", "ncross_dn_ct1", "ncross_dn_ak4")))

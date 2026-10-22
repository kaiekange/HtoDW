"""Z->mumu CR data/MC for the JME variants of jme_zmumu_study.py (2026-10-22).

Same selection as draw_cr_zmumu.py (zmumu skim, phi20, 76 < M(mumu) < 106, OS data), same MC
weight: skim weight x common.br_reweight x common.dy_jetbin_sf (DY). Signal not drawn (signal
is only used in the SR). The event set is identical in every column; only the plotted
jet/MET quantity changes:
    nano      the skim's *_corr columns (= nano jets/MET since 2026-10-21)
    JER       JME JES+JER jets, Type-1 PuppiMET with JER propagated (MC)
    noJER     same jets, Type-1 PuppiMET with JES only
    +xy       the above + MET xy correction
    idfix     JME jets with the NanoAODv12 jet-ID bugfix (jets only)
One figure per quantity (columns = variants, stack + data/MC) and a summary figure overlaying
the data/MC ratios. Figures in MVA/figures/sf/cr_zmumu_jme_*.png.

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
JME = f"{common.SAMPLES_DIR}/CR/zmumu_jme"
M_PHI, ZLO, ZHI = 1.019460, 76.0, 106.0
LUMI_FB = (sum(e["lumi_pb"] for e in common.ERAS.values()) if ERA_ID is None
           else common.ERAS[ERA_ID]["lumi_pb"]) / 1000.0
TAG = f"_{ERA}" if ERA else ""

BKG = [  # as draw_cr_zmumu.py
    (r"$t\bar{t}$ LNu2Q", "TT_merged", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT_merged", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
    ("W+jets", "WtoLNu2Jets_merged", None, "#5a9e78"),
    (r"DY$\to\ell\ell$", "DY_merged", None, "#4a90c2"),
]
GROUPS = ("Data_OS_merged", "TT_merged", "WtoLNu2Jets_merged", "DY_merged")
MET_V = [("nano", "met_pt_corr", "met_phi_corr", "W_tmass_corr")] + [
    (lab, f"jme_met_pt_{v}", f"jme_met_phi_{v}", f"jme_W_tmass_{v}")
    for lab, v in (("JER", "jer"), ("noJER", "nojer"), ("JER+xy", "jer_xy"), ("noJER+xy", "nojer_xy"))]
JET_V = [("nano", "{}_corr"), ("JME", "jme_{}"), ("JME idfix", "jme_{}_idfix")]
QUANT = [  # (name, x label, bins, list of (variant label, column))
    ("met", r"$p_T^{miss}$ [GeV]", np.arange(0, 150 + 1e-9, 5.0), [(l, p) for l, p, _f, _m in MET_V]),
    ("mt", r"$m_T$ [GeV]", np.arange(0, 150 + 1e-9, 5.0), [(l, m) for l, _p, _f, m in MET_V]),
    ("metphi", r"$\phi(p_T^{miss})$", np.linspace(-np.pi, np.pi, 25), [(l, f) for l, _p, f, _m in MET_V]),
    ("njet", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 6.5, 1.0),
     [(l, c.format("num_jets_pt30")) for l, c in JET_V]),
    ("leadpt", r"leading jet $p_T$ [GeV]", np.arange(15, 255 + 1e-9, 10.0),
     [(l, c.format("jet_leadpt")) for l, c in JET_V]),
    ("ht", r"$H_T$ (jets $p_T>20$) [GeV]", np.arange(0, 400 + 1e-9, 20.0),
     [(l, c.format("jet_ht_pt20")) for l, c in JET_V]),
]


def read(group):
    base = ["run", "luminosityBlock", "event", "weight", "phi_mass", "m_mumu", "era",
            "met_pt_corr", "met_phi_corr", "W_tmass_corr", "num_jets_pt30_corr",
            "jet_leadpt_corr", "jet_ht_pt20_corr", "mu_charge", "pi_charge"]
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
            sys.exit(f"{group}: JME friend not row-aligned on {k}")
    a.update({k: v for k, v in j.items() if k.startswith("jme_")})
    sel = (np.abs(a["phi_mass"] - M_PHI) < 0.020) & (a["m_mumu"] > ZLO) & (a["m_mumu"] < ZHI)
    if ERA_ID is not None:
        sel &= a["era"] == ERA_ID
    a = {k: v[sel] for k, v in a.items()}
    if group == "Data_OS_merged":
        a["w"] = np.ones(len(a["event"]))
    else:
        w = a["weight"].astype(float) * common.br_reweight(a)
        a["w"] = w * common.dy_jetbin_sf(a) if group == "DY_merged" else w
    return a


def hist(x, w, bins):
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w * w)
    return h, e2


def label_title(ax):
    hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6, fontsize=14)


S = {g: read(g) for g in GROUPS}
d = S["Data_OS_merged"]
print(f"selected (phi20, Z window{', ' + ERA if ERA else ''}): data {len(d['event'])}, "
      + ", ".join(f"{g} {S[g]['w'].sum():.1f}" for g in GROUPS[1:]))
summary = {}
for name, xlab, bins, cols in QUANT:
    ncol = len(cols)
    fig = plt.figure(figsize=(7.2 * ncol, 7.6))
    outer = gridspec.GridSpec(1, ncol, figure=fig, wspace=0.32, top=0.86, bottom=0.12)
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    summary[name] = []
    for i, (vlab, col) in enumerate(cols):
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
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.5, 1.5)
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(xlab)
        rax.set_ylabel("data/MC", fontsize=14)
        if i == 0:
            label_title(ax)
        summary[name].append((vlab, r, np.sqrt(dh / np.maximum(tot, 1e-12) ** 2 + dh ** 2 * tot_e2 / np.maximum(tot, 1e-12) ** 4)))
        print(f"  {name:7s} {vlab:10s} data {dh.sum():8,}  MC {tot.sum():11,.1f}  data/MC {dh.sum() / tot.sum():.4f}"
              f"  mean data {np.average(d[col][(d[col] >= bins[0]) & (d[col] < bins[-1])]):.2f}")
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR, JME variants{' (' + ERA + ')' if ERA else ''}: "
                 "MC x BR rw x DY jet-bin SF, no DNN cut, same events in every column", fontsize=15, y=0.95)
    out = common.fig_path("sf", f"cr_zmumu_jme_{name}{TAG}.png")
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

# summary: data/MC ratio per variant overlaid (stat of data and MC combined)
fig, axes = plt.subplots(2, 3, figsize=(24, 13))
colors = ["black", "#d62728", "#1f77b4", "#ff7f0e", "#2ca02c"]
for ax, (name, xlab, bins, _cols) in zip(axes.flat, QUANT):
    bc = 0.5 * (bins[1:] + bins[:-1])
    off = np.diff(bins)[0] * 0.08
    for k, (vlab, r, re) in enumerate(summary[name]):
        ax.errorbar(bc + (k - 2) * off, r, yerr=re, fmt="o-", ms=3.5, lw=1, color=colors[k], label=vlab)
    ax.axhline(1, color="0.5", lw=1, ls="--")
    ax.set_ylim(0.5, 1.5)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(xlab)
    ax.set_ylabel("data/MC")
    ax.legend(fontsize=12)
label_title(axes.flat[0])
fig.suptitle(f"Z$\\to\\mu\\mu$ CR data/MC per JME variant{' (' + ERA + ')' if ERA else ''} "
             "(MC x BR rw x DY jet-bin SF, no DNN cut)", fontsize=17)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = common.fig_path("sf", f"cr_zmumu_jme_ratios{TAG}.png")
plt.savefig(out, dpi=100, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

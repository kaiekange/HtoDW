"""Z->mumu CR: is the MC MET too hard because of the JER propagation into MET?

Also reco Njet (pT>30): producer `num_jets_pt30` (GT JEC, no JER) vs `num_jets_pt30_corr`.

Compares data/MC for the producer's MET (`met_pt`: slimmedMETsPuppi, Type-1 with the GT JEC,
no JER smearing) and the analysis MET (`met_pt_corr`: jerc-application-tutorial
get_corrected_met, Type-1 with JES + JER smearing in MC), same for m_T (`W_tmass` vs
`W_tmass_corr`). If the raw MET agrees much better, the culprit is JER in MET; if both are
equally off, it is the unclustered / soft activity that Type-1 does not touch.

Selection as fit_dy_jetbin_zmumu.py: zmumu CR skim, 76 < M(mumu) < 106, phi window 20 MeV,
no DNN cut; DY weighted by common.dy_jetbin_sf. Both MC stacks are normalised to the data
count, so only the SHAPES are compared.

    python3 draw_zmumu_met_raw_vs_corr.py [--key=ssdata_phi20]
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
import common

MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata_phi20")
CR_DIR = f"{common.SAMPLES_DIR}/CR/zmumu"
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
MASS_BINS = np.arange(0.0, 150.0 + 1e-9, 5.0)
PAIRS = [  # (raw column, corrected column, x label, bins)
    ("met_pt", "met_pt_corr", r"$p_T^{miss}$ [GeV]", MASS_BINS),
    ("W_tmass", "W_tmass_corr", r"$m_T$ [GeV]", MASS_BINS),
    ("num_jets_pt30", "num_jets_pt30_corr", r"$N_{jet}$ ($p_T>30$ GeV)", np.arange(-0.5, 7.5, 1.0)),
]


def read(g, extra=()):
    with uproot.open(f"{CR_DIR}/{g}.root")["Events"] as t:
        a = t.arrays(["m_mumu", "phi_mass", "weight", "met_pt", "met_pt_corr", "W_tmass",
                      "W_tmass_corr", "num_jets_pt30", "num_jets_pt30_corr", *extra],
                     library="np")
    m = (a["m_mumu"] > 76) & (a["m_mumu"] < 106) & (np.abs(a["phi_mass"] - common.M_PHI) < 0.020)
    return {k: v[m] for k, v in a.items()}


data = read("Data_OS_merged")
mc = [read("DY_merged", ["jet_bin"]), read("TT_merged"), read("WtoLNu2Jets_merged")]
mc[0]["weight"] = mc[0]["weight"] * common.dy_jetbin_sf(mc[0])


def h(a, col, bins, w=None):
    x = np.clip(a[col], bins[0], bins[-1] - 1e-9)
    return np.histogram(x, bins, weights=w)[0]


fig = plt.figure(figsize=(9 * len(PAIRS), 8.5))
outer = gridspec.GridSpec(1, len(PAIRS), figure=fig, wspace=0.25)
for i, (raw, cor, xlab, BINS) in enumerate(PAIRS):
    bc = 0.5 * (BINS[1:] + BINS[:-1])
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i], height_ratios=[3, 1.3],
                                             hspace=0.06)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    print(f"== {raw} vs {cor}: data/MC per bin (MC normalised to data)")
    for col, colr, lab in ((raw, "#e07a3c", "raw (GT JEC, no JER)"),
                           (cor, "#1f4e79", "corrected (JES+JER)")):
        dh = h(data, col, BINS).astype(float)
        mh = sum(h(a, col, BINS, a["weight"].astype(float)) for a in mc)
        mh = mh * dh.sum() / mh.sum()
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color=colr, label=f"data, {lab}")
        hep.histplot(mh, BINS, ax=ax, color=colr, lw=2, label=f"MC, {lab}")
        with np.errstate(invalid="ignore", divide="ignore"):
            r, re = dh / mh, np.sqrt(dh) / mh
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=4, color=colr, label=lab)
        print(f"   {lab:22s}: " + " ".join(f"{x:.2f}" for x in r[:20]))
    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_ylim(0, None)
    ax.set_ylabel("Events / bin", fontsize=16)
    ax.tick_params(labelbottom=False)
    ax.legend(fontsize=12)
    rax.axhline(1, color="red", lw=1)
    rax.set_ylim(0.3, 1.7)
    rax.set_xlim(BINS[0], BINS[-1])
    rax.set_xlabel(xlab)
    rax.set_ylabel("data/MC (shape)", fontsize=14)
    if i == 0:
        hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6,
                      fontsize=14)
fig.suptitle(r"Z$\to\mu\mu$ CR, 76<M($\mu\mu$)<106, $|M_{KK}-M_\phi|$<20 MeV: raw vs corrected MET, m_T, N_jet "
             r"(MC normalised to data, DY x jet-bin SF)", fontsize=16)
out = common.fig_path("sf", "cr_zmumu_met_raw_vs_corr.png", key=MODEL_KEY)
plt.savefig(out, dpi=100, bbox_inches="tight")
print(f"Saved {out}")

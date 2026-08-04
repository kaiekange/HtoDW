"""Signal plots for H->DsW: OS/SS lepton multiplicity and the phi / Ds mass peaks.

Usage (needs the ds_dnn conda env for mplhep, same as the MVA scripts):
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_lepton_counts.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

INFILE = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_phipimu/HtoDsW.root"
TREE = "HphipimuDsW/Events"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"

M_PHI, M_DS, M_K = 1.019460, 1.96835, 0.493677

# (panel title, OS branch, SS branch)
PANELS = [
    ("muon - loose",         "num_os_mu_loose",     "num_ss_mu_loose"),
    ("muon - medium",        "num_os_mu_medium",    "num_ss_mu_medium"),
    ("muon - tight",         "num_os_mu_tight",     "num_ss_mu_tight"),
    ("electron - cut loose",  "num_os_e_cutloose",  "num_ss_e_cutloose"),
    ("electron - cut medium", "num_os_e_cutmedium", "num_ss_e_cutmedium"),
    ("electron - cut tight",  "num_os_e_cuttight",  "num_ss_e_cuttight"),
    ("electron - MVA WP80",   "num_os_e_mvawp80",   "num_ss_e_mvawp80"),
    ("electron - MVA WP90",   "num_os_e_mvawp90",   "num_ss_e_mvawp90"),
]

lep_branches = sorted({b for _, o, s in PANELS for b in (o, s)})
a = uproot.open(INFILE)[TREE].arrays(
    lep_branches + ["phi_mass", "Ds_mass", "run_id", "lumi_id", "event_id",
                    "cand_id", "event_match"],
    library="np")

# event_match (== Ds_match && mu_match) picks the one truth-matched Ds candidate per
# event, but the gamma loop writes one row per (Ds, gamma) pair -- so that candidate can
# span several rows. Deduplicate on (event, cand_id) to get one row per candidate.
m = a["event_match"].astype(bool)
key = np.stack([a["run_id"][m], a["lumi_id"][m], a["event_id"][m], a["cand_id"][m]])
_, first = np.unique(key, axis=1, return_index=True)
sel = np.flatnonzero(m)[np.sort(first)]
n = sel.size
print(f"rows: {len(m)} | event_match rows: {m.sum()} | truth-matched candidates: {n}")

# ---------- OS vs SS lepton multiplicity ----------
maxval = max(int(a[b][sel].max()) for b in lep_branches)
bins = np.arange(-0.5, maxval + 1.5, 1.0)

fig, axes = plt.subplots(3, 3, figsize=(24, 18))
axes = axes.flatten()
for i, (title, ob, sb) in enumerate(PANELS):
    ax = axes[i]
    os_v, ss_v = a[ob][sel], a[sb][sel]
    ax.hist(os_v, bins=bins, alpha=0.6, color="steelblue", label=f"OS  (mean {os_v.mean():.3f})")
    ax.hist(ss_v, bins=bins, alpha=0.6, color="firebrick", label=f"SS  (mean {ss_v.mean():.3f})")
    ax.set_xlabel("N leptons / event")
    ax.set_ylabel("# events")
    ax.set_title(title, fontsize=20)
    ax.set_xticks(np.arange(0, maxval + 1))
    ax.legend(fontsize=15)
    if i == 0:
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False, fontsize=16)
for j in range(len(PANELS), len(axes)):
    axes[j].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/lepton_counts.png", dpi=150)
plt.close(fig)
print(f"Saved {OUTDIR}/lepton_counts.png")

# ---------- phi mass ----------
# Plot exactly the range the tuplizer keeps: the selection applies an upper cut at
# M_phi + 50 MeV, and the lower side is bounded by the KK threshold (2*m_K).
v = a["phi_mass"][sel]
fig, ax = plt.subplots(figsize=(10, 8))
ax.hist(v, bins=100, range=(2 * M_K, M_PHI + 0.05), color="steelblue")
ax.set_xlabel(r"$M(K^{+}K^{-})$ [GeV]")
ax.set_ylabel("# events")
ax.set_ylim(0, 1.3 * ax.get_ylim()[1])
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/phi_mass.png", dpi=150)
plt.close(fig)
print(f"Saved {OUTDIR}/phi_mass.png  (median {np.median(v):.5f})")

# ---------- Ds mass ----------
# Plot exactly the range the tuplizer keeps: |M - M_Ds| < 0.3.
v = a["Ds_mass"][sel]
fig, ax = plt.subplots(figsize=(10, 8))
ax.hist(v, bins=120, range=(M_DS - 0.3, M_DS + 0.3), color="steelblue")
ax.set_xlabel(r"$M(K^{+}K^{-}\pi^{\pm})$ [GeV]")
ax.set_ylabel("# events")
ax.set_ylim(0, 1.3 * ax.get_ylim()[1])
hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=False)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/Ds_mass.png", dpi=150)
plt.close(fig)
print(f"Saved {OUTDIR}/Ds_mass.png  (median {np.median(v):.5f})")

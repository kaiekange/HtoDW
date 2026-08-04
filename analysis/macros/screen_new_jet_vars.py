"""Screen the 18 jet/b-tag variables added to the tuplizer on 2026-07-23 before
committing them to the MVA feature set.

Three things decide whether a variable earns its slot:
  1. separation      -- does signal differ from background?
  2. data/MC modelling -- does the background MC reproduce data? A variable the MC
     gets wrong is a liability in a classifier no matter how discriminating it looks.
  3. mass sculpting  -- correlation with Ds_mass / phi_mass. The existing 37-feature
     set was chosen to keep this near zero so the MVA does not carve a fake peak into
     the fit variable; the new 18 have never been screened this way.

Data is the OS SR-blinded sample, so it has no |M-M_Ds|<40 MeV events. To compare
like with like, the same sideband requirement is applied to the background MC in
the data/MC panels (its size is reported so the effect is visible, not assumed).

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 screen_new_jet_vars.py
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
OUTDIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/HtoDsW/figures"
M_DS = 1.96835
SR_HALFWIDTH = 0.04

NEW_VARS = [
    "Ds_nearestjet_bscore", "Ds_nearestjet_dR",
    "jet_ht", "jet_ht_pt20", "jet_leadpt", "jet_maxbscore",
    "n_b_jets_loose_dR03", "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR03", "n_b_jets_medium_dR04",
    "n_b_jets_tight_dR03", "n_b_jets_tight_dR04",
    "n_b_jets_loose_pt20_dR03", "n_b_jets_loose_pt20_dR04",
    "n_b_jets_medium_pt20_dR03", "n_b_jets_medium_pt20_dR04",
    "n_b_jets_tight_pt20_dR03", "n_b_jets_tight_pt20_dR04",
]

WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}

SEL = ["num_os_mu_tight", "num_ss_mu_tight"]
EXTRA = ["Ds_mass", "phi_mass"]


def one_tight_mu(a):
    return (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)


def load(name, truth_filter=False):
    br = NEW_VARS + SEL + EXTRA + (["Ds_match", "mu_match"] if truth_filter else [])
    a = uproot.open(f"{SAMPLES_DIR}/{name}.root")["Events"].arrays(br, library="np")
    m = one_tight_mu(a)
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    return {k: v[m] for k, v in a.items()}


print("Loading...")
sig = load("HtoDsW", truth_filter=True)
data = load("Data_2022")
bkg_parts, bkg_w_parts = [], []
for name, wgt in WJETS_WEIGHTS.items():
    d = load(name)
    bkg_parts.append(d)
    bkg_w_parts.append(np.full(len(d["Ds_mass"]), wgt))
bkg = {k: np.concatenate([p[k] for p in bkg_parts]) for k in bkg_parts[0]}
bkg_w = np.concatenate(bkg_w_parts)

print(f"signal {len(sig['Ds_mass']):,} | background {len(bkg['Ds_mass']):,} | data {len(data['Ds_mass']):,}")

# background restricted to the same sidebands as the blinded data
bkg_sb = np.abs(bkg["Ds_mass"] - M_DS) >= SR_HALFWIDTH
print(f"background in sidebands (matching blinded data): {bkg_sb.sum():,} "
      f"({100 * bkg_sb.mean():.1f}%)")

# ---------- sentinel audit ----------
print("\n=== sentinel fractions (-777 = not computable, ~-10 = tagger unavailable) ===")
for v in NEW_VARS:
    s777 = 100 * np.mean(bkg[v] == -777)
    sneg = 100 * np.mean((bkg[v] < -1) & (bkg[v] != -777))
    if s777 > 0.01 or sneg > 0.01:
        print(f"  {v:32s} bkg: -777 {s777:5.2f}%   other-negative {sneg:5.2f}%")

# ---------- sculpting check ----------
# Pearson r against the two fit variables, on background only (that is where a
# sculpted shape would fake a peak). |r| well under ~0.1 is the working target.
print("\n=== sculpting risk: corr(feature, fit variable) on background ===")
print(f"{'feature':32s} {'r(Ds_mass)':>11s} {'r(phi_mass)':>12s}")
sculpt = {}
for v in NEW_VARS:
    ok = bkg[v] > -1          # drop sentinels, they would fake a correlation
    if ok.sum() < 100:
        continue
    r_ds = np.corrcoef(bkg[v][ok], bkg["Ds_mass"][ok])[0, 1]
    r_phi = np.corrcoef(bkg[v][ok], bkg["phi_mass"][ok])[0, 1]
    sculpt[v] = (r_ds, r_phi)
    flag = "  <-- check" if max(abs(r_ds), abs(r_phi)) > 0.05 else ""
    print(f"{v:32s} {r_ds:11.4f} {r_phi:12.4f}{flag}")

# ---------- comparison plot ----------
print("\nDrawing...")
N_COLS, N_ROWS = 5, 4
fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(26, 19))
axes_flat = axes.flatten()

for i, v in enumerate(NEW_VARS):
    ax = axes_flat[i]
    s, b, d = sig[v], bkg[v], data[v]
    bw = bkg_w
    # sentinels are excluded from the plotted range but their rate is in the title
    ok_all = np.concatenate([s[s > -1], b[b > -1], d[d > -1]])
    if v.startswith("n_b_jets"):
        bins = np.arange(-0.5, 5.5, 1.0)
    else:
        lo, hi = np.percentile(ok_all, [0.5, 99.5])
        bins = np.linspace(lo, hi, 50) if hi > lo else np.linspace(lo - 1, hi + 1, 50)

    ax.hist(b[bkg_sb], bins=bins, weights=bw[bkg_sb], density=True, histtype="step",
            lw=2, color="steelblue", label="W+jets (sideband)")
    ax.hist(d, bins=bins, density=True, histtype="step", lw=2, ls="--",
            color="black", label="Data (SR blinded)")
    ax.hist(s, bins=bins, density=True, alpha=0.45, color="firebrick", label="HtoDsW")
    r = sculpt.get(v, (np.nan, np.nan))[0]
    ax.set_title(f"{v}\n" + r"$r(M_{D_s})$=" + f"{r:+.3f}", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("a.u.")
    ax.legend(fontsize=11)
    if i == 0:
        hep.cms.label(ax=ax, rlabel="2022 (13.6 TeV)", data=True, text="Preliminary", fontsize=12)

for j in range(len(NEW_VARS), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.tight_layout()
out = f"{OUTDIR}/new_jet_vars_screening.png"
plt.savefig(out, dpi=130)
print(f"Saved {out}")

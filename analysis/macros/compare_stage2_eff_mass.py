"""At matched (WP1, stage-2 eff_sig), compare eff_bkg and Ds_mass shape between the
two stage-2 training populations -- "wjetsres" (trained behind stage-1 eff_sig=0.90)
vs "wjetsres_nostage1" (trained on all of GenTag modes 1-5, no stage-1 pre-filter).

2026-09-18: compare_stage2_wp_fom.py compared the two via a Punzi-style FOM
(eps/sqrt(B) and eps/(1.5+sqrt(B))), which folds signal efficiency and background
size into one number. The user asked for the more direct, more interpretable
comparison instead: hold WP1 and stage-2's OWN (conditional) signal efficiency
fixed and look at eff_bkg alone (not a significance figure of merit), plus the
Ds_mass shape of what survives -- a FOM can be equal while the shape underneath
differs, and a FOM ratio does not show that.

eff_bkg is reported two ways, same background definitions as compare_stage2_wp_fom.py:
    eff_bkg(OS data sb) = N_pass / N_total on the OS mass sideband
                          (|M(KKpi)-M_Ds| >= 40 MeV, OS is SR-blind anyway)
    eff_bkg(GenTag MC)  = weighted sum(modes 1-5, SR window) after cut
                          / weighted sum(modes 1-5, SR window) before any stage-1/2 cut
                          (denominator is OS+b-veto selected, no score cut)

wjetsres_nostage1 has stage1_cut=None, so its own Prediction/ files carry
score_dsmu_res but not score_dsmu (stage-1's score never entered that path). Both
this variant's and stage-1's own Prediction/<sample>.root are scored from the same
SAMPLES_DIR/<sample>.root with no row filtering before scoring, so they share row
count and on-disk order -- asserted here via (run, luminosityBlock, event) rather
than assumed (see compare_stage2_wp_fom.py's load_joined for the same pattern).

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_stage2_eff_mass.py [WP1] [eff_sig2]   # defaults 0.90 0.70, table always prints all of WP1_GRID x EFF2_GRID
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

STAGE1_KEY = "ssdata"
STAGE2_KEYS = ("wjetsres", "wjetsres_nostage1")
STAGE2_TOOL = "BDT"
STAGE1_PRED = f"{common.out_dir_for(STAGE1_KEY, tool='DNN')}/Prediction"

MC_SAMPLE = "WtoLNu2Jets_merged"
DATA_SAMPLE = "Data_OS_merged"

M_DS, SR_HALF = 1.96835, 0.040
NONRES = [0, 6, 7, 8, 9, 10]
RES = [1, 2, 3, 4, 5]

WP1_GRID = [0.90, 0.70, 0.50]
EFF2_GRID = [0.90, 0.70, 0.50]

IDS = ["run", "luminosityBlock", "event"]
SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO]


def sel_mask(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
            & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1) & (a["score_dsmu_res"] > -1))


def load_joined(sample, stage2_key, extra_cols):
    with uproot.open(f"{STAGE1_PRED}/{sample}.root")["Events"] as t:
        a1 = t.arrays(["score_dsmu"] + IDS + SEL + extra_cols, library="np")
    s2_dir = f"{common.out_dir_for(stage2_key, tool=STAGE2_TOOL)}/Prediction"
    with uproot.open(f"{s2_dir}/{sample}.root")["Events"] as t:
        a2 = t.arrays(["score_dsmu_res"] + IDS, library="np")
    if len(a1["event"]) != len(a2["event"]):
        sys.exit(f"{sample}: row count mismatch, stage-1 {len(a1['event']):,} vs "
                 f"{stage2_key} {len(a2['event']):,}")
    for c in IDS:
        if not np.array_equal(a1[c], a2[c]):
            sys.exit(f"{sample}: column {c} disagrees between stage-1 and {stage2_key} "
                     f"Prediction files -- not row-aligned.")
    out = dict(a1)
    out["score_dsmu_res"] = a2["score_dsmu_res"]
    return out


def load_variant(stage2_key):
    """Everything one variant needs: truth-matched signal (s1,s2,sw), W+jets MC
    (b1,b2,mode,mass,weight), OS data (d1,d2,mass)."""
    s1, s2, sw = [], [], []
    for s in common.SIGNALS:
        a = load_joined(s, stage2_key, ["Ds_match", "mu_match", "weight"])
        m = sel_mask(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
        s1.append(a["score_dsmu"][m]); s2.append(a["score_dsmu_res"][m])
        sw.append(a["weight"][m])
    s1, s2, sw = np.concatenate(s1), np.concatenate(s2), np.concatenate(sw)

    b = load_joined(MC_SAMPLE, stage2_key, ["cand_gen_mode", "Ds_mass", "weight"])
    kb = sel_mask(b)
    b1, b2 = b["score_dsmu"][kb], b["score_dsmu_res"][kb]
    b_mode, b_mass, b_w = b["cand_gen_mode"][kb], b["Ds_mass"][kb], b["weight"][kb]

    d = load_joined(DATA_SAMPLE, stage2_key, ["Ds_mass"])
    kd = sel_mask(d)
    d1, d2, d_mass = d["score_dsmu"][kd], d["score_dsmu_res"][kd], d["Ds_mass"][kd]

    return dict(s1=s1, s2=s2, sw=sw, b1=b1, b2=b2, b_mode=b_mode, b_mass=b_mass, b_w=b_w,
               d1=d1, d2=d2, d_mass=d_mass)


print("Loading both variants (this reads stage-1 + each stage-2's Prediction/ for "
      "Signal_merged, WtoLNu2Jets_merged, Data_OS_merged -- a few minutes)...")
V = {k: load_variant(k) for k in STAGE2_KEYS}

# ---------- eff_bkg denominators: OS+b-veto selected, no score cut at all ----------
D_TOT = {k: len(V[k]["d1"]) for k in STAGE2_KEYS}
b_res_mask = {k: np.isin(V[k]["b_mode"], RES) for k in STAGE2_KEYS}
b_sr_mask = {k: np.abs(V[k]["b_mass"] - M_DS) < SR_HALF for k in STAGE2_KEYS}
MC_RES_SR_TOT = {k: float(V[k]["b_w"][b_res_mask[k] & b_sr_mask[k]].sum()) for k in STAGE2_KEYS}
for k in STAGE2_KEYS:
    print(f"{k}: D_TOT(OS, no score cut)={D_TOT[k]:,}  "
          f"MC_RES_SR_TOT(modes 1-5, SR, no score cut)={MC_RES_SR_TOT[k]:.6g}")

# ---------- table: eff_bkg at matched (WP1, eff_sig2), both variants ----------
print(f"\n{'WP1':>5} {'eff2':>5} "
      f"{'eff_bkg(OS sb) a':>17} {'eff_bkg(OS sb) b':>17} "
      f"{'eff_bkg(MC res) a':>18} {'eff_bkg(MC res) b':>18}")
results = {}
for wp1 in WP1_GRID:
    for eff2 in EFF2_GRID:
        row = {}
        for k in STAGE2_KEYS:
            d = V[k]
            thr1 = float(np.quantile(d["s1"], 1 - wp1))
            s2_cond = d["s2"][d["s1"] > thr1]
            thr2 = float(np.quantile(s2_cond, 1 - eff2))

            pass_d = (d["d1"] > thr1) & (d["d2"] > thr2)
            sb = np.abs(d["d_mass"] - M_DS) >= SR_HALF
            n_sb_pass = int(np.sum(pass_d & sb))
            n_sb_tot = int(np.sum(sb))
            eff_data = n_sb_pass / n_sb_tot if n_sb_tot > 0 else np.nan

            pass_b = (d["b1"] > thr1) & (d["b2"] > thr2)
            w_pass = d["b_w"][pass_b & b_res_mask[k] & b_sr_mask[k]]
            eff_mc = float(w_pass.sum()) / MC_RES_SR_TOT[k] if MC_RES_SR_TOT[k] > 0 else np.nan

            row[k] = dict(thr1=thr1, thr2=thr2, eff_data=eff_data, eff_mc=eff_mc,
                         n_sb_pass=n_sb_pass, n_sb_tot=n_sb_tot)
        results[(wp1, eff2)] = row
        a, b = row["wjetsres"], row["wjetsres_nostage1"]
        print(f"{wp1:>5.2f} {eff2:>5.2f} "
              f"{a['eff_data']:>17.5f} {b['eff_data']:>17.5f} "
              f"{a['eff_mc']:>18.5f} {b['eff_mc']:>18.5f}")

# ---------- Ds_mass shape at one chosen operating point ----------
WP1_PLOT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.90
EFF2_PLOT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.70
if (WP1_PLOT, EFF2_PLOT) not in results:
    sys.exit(f"(WP1={WP1_PLOT}, eff2={EFF2_PLOT}) not in the scanned grid "
             f"{WP1_GRID} x {EFF2_GRID}")
row = results[(WP1_PLOT, EFF2_PLOT)]

BINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.020)
bc = 0.5 * (BINS[1:] + BINS[:-1])
width = BINS[1] - BINS[0]
blind = np.abs(bc - M_DS) < SR_HALF
vis = ~blind

COLORS = {"wjetsres": "firebrick", "wjetsres_nostage1": "steelblue"}
LABELS = {"wjetsres": "wjetsres (stage-1 90% pre-filter)",
          "wjetsres_nostage1": "wjetsres_nostage1 (no pre-filter)"}

fig, axes = plt.subplots(1, 2, figsize=(19, 8))

# left: OS data sideband (SR blind, absolute counts, both variants overlaid)
ax = axes[0]
for k in STAGE2_KEYS:
    d = V[k]
    thr1, thr2 = row[k]["thr1"], row[k]["thr2"]
    pass_d = (d["d1"] > thr1) & (d["d2"] > thr2)
    vals = d["d_mass"][pass_d]
    h, _ = np.histogram(vals, bins=BINS)
    h_m = np.where(vis, h, np.nan)
    ax.errorbar(bc, h_m, yerr=np.sqrt(np.where(vis, h, 0)), fmt="o", ms=4,
                color=COLORS[k], elinewidth=1.2,
                label=f"{LABELS[k]}  N={h[vis].sum():,}")
ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
ax.set_xlim(BINS[0], BINS[-1])
ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
ax.set_title("OS data (SR blind)")
ax.legend(fontsize=10)
hep.cms.label(ax=ax, data=True, text="Preliminary",
             rlabel=rf"{(common.ERAS[1]['lumi_pb']+common.ERAS[2]['lumi_pb'])/1000:.1f} fb$^{{-1}}$ (13.6 TeV)",
             fontsize=10)

# right: W+jets GenTag modes 1-5 (weighted, SR not blinded -- it's MC)
ax = axes[1]
for k in STAGE2_KEYS:
    d = V[k]
    thr1, thr2 = row[k]["thr1"], row[k]["thr2"]
    pass_b = (d["b1"] > thr1) & (d["b2"] > thr2) & b_res_mask[k]
    vals, w = d["b_mass"][pass_b], d["b_w"][pass_b]
    h, _ = np.histogram(vals, bins=BINS, weights=w)
    e2, _ = np.histogram(vals, bins=BINS, weights=w ** 2)
    ax.errorbar(bc, h, yerr=np.sqrt(e2), fmt="o", ms=4, color=COLORS[k], elinewidth=1.2,
                label=f"{LABELS[k]}  N(weighted)={h.sum():,.1f}")
ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
ax.set_xlim(BINS[0], BINS[-1])
ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
ax.set_title("W+jets GenTag modes 1-5 (true D$_s$/D, weighted)")
ax.legend(fontsize=10)
hep.cms.label(ax=ax, data=False, text="Simulation", fontsize=10)

fig.suptitle(rf"WP1 $\varepsilon_1$={WP1_PLOT:.0%}, stage-2 $\varepsilon_2$={EFF2_PLOT:.0%} "
             r"(conditional) -- both stage-2 variants at matched (WP1,$\varepsilon_2$)",
             fontsize=13)
plt.tight_layout()
out = common.fig_path("comparisons", f"stage2_variant_mass_wp1_{WP1_PLOT:.2f}_eff2_{EFF2_PLOT:.2f}.png",
                      key="wjetsres", tool="BDT")
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

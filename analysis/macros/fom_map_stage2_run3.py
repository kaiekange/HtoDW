"""2D (WP1, WP2) Punzi FOM maps for the two-stage MVA, scaled to Run-3 (326 /fb).

    FOM = eps_tot / (1.5 + sqrt(B))        Punzi a=3

Punzi is BR-INDEPENDENT, which is why it is the right choice: the signal
normalisation carries whatever BR(H->DsW) the sample was generated with, so any
S/sqrt(B) map would inherit that assumption.

TWO MAPS, ONE PER BACKGROUND DEFINITION (2026-09-22 request):
  B_data    OS data in the CONTROL REGION, |M(KKpi) - M_Ds| >= 40 MeV, counted as it
            stands. The SR stays blind and is NOT extrapolated into. Model-free, and
            it contains tt/DY/QCD, at the cost of not being the SR.
  B_wjets   WtoLNu2Jets_merged ENTIRE, weight-weighted: ALL cand_gen_mode (resonant
            AND non-resonant) AND the full M(KKpi) range, no mass window of any kind.
            tt/DY/QCD are absent entirely (the pre-cut data/MC gap is 1.44 and rises
            steeply with jet multiplicity), so this UNDERSTATES B.

THE TWO MAPS ARE NOT COMPARABLE IN ABSOLUTE VALUE. B_wjets spans the full 600 MeV
mass range while B_data has the 80 MeV blind window punched out of it, so B_wjets is
the larger of the two on window alone. What IS comparable across the two is WHERE IN
(WP1, WP2) the optimum sits. Data cannot be given the same full-range treatment --
its SR is blind, and the control region is all there is.

WHY NOT compare_stage2_wp_fom.py: that script exists to compare the two stage-2
TRAINING POPULATIONS (wjetsres vs wjetsres_nostage1) and is the record of the
2026-09-22 finding that they are statistically indistinguishable. That question is
settled -- only the stage-1-based `wjetsres` is used now -- so this script scans the
single surviving model instead, at Run-3 scale, with the backgrounds the 0922
request asked for. fom_scan_stage2d.py is the pre-migration (2026-09-16) historical
record and is likewise left alone.

WP2 IS AN UNCONDITIONAL EFFICIENCY (2026-09-22): its threshold is a quantile of ALL
truth-matched signal, not of the subset that already passed WP1. So WP2=0.70 is the
same stage-2 cut in every row of the map, which is what makes a column readable.
The price is that eps_tot is NOT eps_1*eps_2 -- the two cuts are correlated -- so
eps_tot is counted directly from the AND and printed in every cell. Guard 3 checks
the WP2 quantile against its nominal efficiency instead.

WP1 IS CAPPED AT 0.90: `wjetsres` was trained behind the stage-1 eff_sig=0.90
pre-filter, so above that it is being asked about a population it never saw.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 fom_map_stage2_run3.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common
import eras

plt.style.use(hep.style.CMS)

STAGE1_KEY = "ssdata"
STAGE2_KEY = "wjetsres"
STAGE2_TOOL = "BDT"
PRED = f"{common.out_dir_for(STAGE2_KEY, tool=STAGE2_TOOL)}/Prediction"

MC_SAMPLE = "WtoLNu2Jets_merged"
DATA_SAMPLE = "Data_OS_merged"

M_DS, SR_HALF = 1.96835, 0.040
PUNZI_A = 3.0

LUMI_NOW_PB = sum(e["lumi_pb"] for e in eras.ERAS.values())     # 2022 preEE+postEE
LUMI_RUN3_PB = 326000.0
SCALE = LUMI_RUN3_PB / LUMI_NOW_PB

WP1_GRID = [round(0.90 - 0.05 * i, 2) for i in range(13)]   # 0.90 .. 0.30
WP2_GRID = [round(0.95 - 0.05 * i, 2) for i in range(14)]   # 0.95 .. 0.30
SEL = [*common.SELECTION]


def sel_mask(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a) & (a["score_dsmu"] > -1) & (a["score_dsmu_res"] > -1)


def neff(w):
    """Kish effective sample size; the relative error on sum(w) is 1/sqrt(neff)."""
    s2 = np.sum(w ** 2)
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


def punzi(eps, B):
    return eps / (PUNZI_A / 2 + np.sqrt(B)) if B > 0 else np.nan


def load(sample, extra):
    """Both scores live in stage-2's own Prediction file: `wjetsres` is scored from
    stage-1's Prediction dir, so score_dsmu rides along. No positional join needed
    (unlike wjetsres_nostage1, which is scored straight from SAMPLES_DIR)."""
    cols = ["score_dsmu", "score_dsmu_res"] + SEL + extra
    return uproot.open(f"{PRED}/{sample}.root")["Events"].arrays(cols, library="np")


# ---------- guard 1: the stage-1 cut stage-2 trained behind is still the same ----------
cfg = common.resolve(STAGE2_KEY)
thr_now, _ = common.stage1_threshold(cfg["stage1_cut"])
with open(f"{common.out_dir_for(STAGE2_KEY, tool=STAGE2_TOOL)}/stage1_filter.txt") as f:
    stamped = dict(line.split() for line in f)
assert abs(thr_now - float(stamped["thr"])) < 1e-6, (
    f"stage-1 threshold drift: stamped {float(stamped['thr']):.6f} vs recomputed "
    f"{thr_now:.6f} -- stage-1 was re-scored after stage-2 trained behind another cut.")
print(f"guard 1 OK: stage1_filter.txt thr={float(stamped['thr']):.6f} matches recomputed")
print(f"Luminosity: {LUMI_NOW_PB/1000:.2f} /fb -> {LUMI_RUN3_PB/1000:.0f} /fb "
      f"(scale {SCALE:.4f})")

# ---------- signal ----------
s1, s2, sw = [], [], []
for s in common.SIGNALS:
    a = load(s, ["Ds_match", "mu_match", "weight"])
    m = sel_mask(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s1.append(a["score_dsmu"][m])
    s2.append(a["score_dsmu_res"][m])
    sw.append(a["weight"][m])
s1, s2, sw = np.concatenate(s1), np.concatenate(s2), np.concatenate(sw)
S_TOT = sw.sum()
print(f"truth-matched signal: {len(s1):,} rows, weighted {S_TOT:.4f} events "
      f"-> {S_TOT*SCALE:.3f} at Run-3 (context only; Punzi is BR-independent)")
THR1 = {wp1: float(np.quantile(s1, 1 - wp1)) for wp1 in WP1_GRID}

# ---------- background A: the ENTIRE W+jets sample (all modes, full mass range) ----------
b = load(MC_SAMPLE, ["cand_gen_mode", "Ds_mass", "weight"])
kb = sel_mask(b)
b1, b2 = b["score_dsmu"][kb], b["score_dsmu_res"][kb]
b_w, b_mode = b["weight"][kb], b["cand_gen_mode"][kb]
uniq = np.unique(b_mode)
assert np.isin(uniq, np.arange(0, 11)).all(), f"cand_gen_mode outside 0..10: {uniq}"
print(f"guard 2 OK: cand_gen_mode present = {uniq.tolist()}; ALL of them enter B_wjets")
print(f"{MC_SAMPLE}: {kb.sum():,} rows after selection, ALL of them enter B_wjets "
      f"(no mass window), weighted {b_w.sum()*SCALE:,.1f} at Run-3")

# ---------- background B: OS data control region (SR left blind) ----------
d = load(DATA_SAMPLE, ["Ds_mass"])
kd = sel_mask(d)
d1, d2 = d["score_dsmu"][kd], d["score_dsmu_res"][kd]
d_cr = np.abs(d["Ds_mass"][kd] - M_DS) >= SR_HALF
print(f"{DATA_SAMPLE}: {kd.sum():,} rows after selection, "
      f"{int(d_cr.sum()):,} in the control region\n")

# ---------- scan ----------
# THR2 on ALL truth-matched signal, so a WP2 column means one fixed stage-2 cut.
THR2 = {wp2: float(np.quantile(s2, 1 - wp2)) for wp2 in WP2_GRID}
rows = []
for wp1 in WP1_GRID:
    thr1 = THR1[wp1]
    p1s = s1 > thr1
    p1b, p1d = b1 > thr1, d1 > thr1
    for wp2 in WP2_GRID:
        thr2 = THR2[wp2]
        eps_and = float(sw[p1s & (s2 > thr2)].sum() / S_TOT)
        eps2 = float(sw[s2 > thr2].sum() / S_TOT)

        w_sr = b_w[p1b & (b2 > thr2)]
        n_cr = int(np.sum(p1d & (d2 > thr2) & d_cr))
        rows.append(dict(
            wp1=wp1, wp2=wp2, eps_and=eps_and, eps2=eps2,
            B_w=float(w_sr.sum()) * SCALE, raw_w=int(w_sr.size), neff_w=neff(w_sr),
            n_cr=n_cr, B_d=n_cr * SCALE))
for r in rows:
    r["fom_w"] = punzi(r["eps_and"], r["B_w"])
    r["fom_d"] = punzi(r["eps_and"], r["B_d"])

worst = max(rows, key=lambda r: abs(r["eps2"] - r["wp2"]))
diff = abs(worst["eps2"] - worst["wp2"])
assert diff < 5e-3, (f"WP2 quantile off by {diff:.4f} at WP2={worst['wp2']}: the "
                     f"threshold is not delivering its nominal unconditional efficiency")
print(f"guard 3 OK: max |eps_2 - WP2| = {diff:.5f} across the grid "
      f"(WP2 is unconditional)\n")

hdr = (f"{'WP1':>5} {'WP2':>5} {'eps_tot':>8} | {'B_wjets':>9} {'rawB':>6} {'NeffB':>6} "
       f"{'FOM_w':>8} | {'n_CR':>6} {'B_data':>10} {'FOM_d':>8}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    print(f"{r['wp1']:>5.2f} {r['wp2']:>5.2f} {r['eps_and']:>8.4f} | "
          f"{r['B_w']:>9.1f} {r['raw_w']:>6d} {r['neff_w']:>6.0f} {r['fom_w']:>8.5f} | "
          f"{r['n_cr']:>6d} {r['B_d']:>10.1f} {r['fom_d']:>8.5f}")

for key, name in (("fom_w", "B = W+jets entire (all modes, full mass range)"),
                  ("fom_d", "B = OS data (control region)")):
    bst = max(rows, key=lambda r: (r[key] if not np.isnan(r[key]) else -np.inf))
    print(f"\n{name}: max FOM {bst[key]:.5f} at WP1={bst['wp1']:.2f}, WP2={bst['wp2']:.2f} "
          f"(eps_tot={bst['eps_and']:.4f}, B={bst[key.replace('fom', 'B')]:.1f}, "
          f"rawB/n_CR={bst['raw_w'] if key == 'fom_w' else bst['n_cr']})")


# ---------- figures: one per background ----------
def grid_of(key):
    return np.array([[next(r[key] for r in rows if r["wp1"] == w1 and r["wp2"] == w2)
                      for w2 in WP2_GRID] for w1 in WP1_GRID])


eps_tot = grid_of("eps_and")
RLABEL = rf"{LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$ (13.6 TeV)"
PANELS = (
    ("wjets", "fom_w", False,
     r"B = W+jets entire: all GenTag modes, full $M(KK\pi)$ range, no mass window"),
    ("data", "fom_d", True,
     r"B = OS data, control region ($|M-M_{D_s}|\geq 40$ MeV, SR blind)"),
)
for tag, key, is_data, subtitle in PANELS:
    Z = grid_of(key)
    fig, ax = plt.subplots(figsize=(16, 11))
    im = ax.imshow(Z, cmap="viridis", aspect="auto", origin="upper")
    for i in range(len(WP1_GRID)):
        for j in range(len(WP2_GRID)):
            color = "white" if Z[i, j] < 0.6 * np.nanmax(Z) else "black"
            ax.text(j, i - 0.15, f"{Z[i, j]:.4f}", ha="center", va="center",
                    fontsize=13, color=color, fontweight="bold")
            ax.text(j, i + 0.20, r"$\varepsilon_{tot}$=" + f"{eps_tot[i, j]:.2f}",
                    ha="center", va="center", fontsize=8.5, color=color)
    ax.set_xticks(range(len(WP2_GRID)))
    ax.set_xticklabels([f"{w:.2f}" for w in WP2_GRID], fontsize=9, rotation=90)
    ax.set_yticks(range(len(WP1_GRID)))
    ax.set_yticklabels([f"{w:.2f}" for w in WP1_GRID], fontsize=9)
    ax.set_xlabel(r"WP2 $\varepsilon_2$")
    ax.set_ylabel(r"WP1 $\varepsilon_1$")
    ax.set_title(r"Punzi FOM $=\varepsilon_{tot}/(1.5+\sqrt{B})$, "
                 rf"scaled to {LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$" "\n" + subtitle,
                 fontsize=13)
    fig.colorbar(im, ax=ax, label="Punzi FOM")
    # ALWAYS text="Preliminary"; `data=` alone decides the "Simulation" prefix, giving
    # "CMS Preliminary" / "CMS Simulation Preliminary". Never text="Simulation" (that
    # doubles it) and never drop text= (that loses Preliminary). See feedback_plot_style.
    hep.cms.label(ax=ax, rlabel=RLABEL, data=is_data, text="Preliminary",
                  fontsize=11)
    plt.tight_layout()
    out = common.fig_path("comparisons", f"fom_map_run3_{tag}.png",
                          key=STAGE2_KEY, tool=STAGE2_TOOL)
    plt.savefig(out, dpi=120)
    plt.close(fig)
    print(f"Saved {out}")

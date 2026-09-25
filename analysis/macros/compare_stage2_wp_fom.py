"""2D (WP1, WP2) working-point FOM scan, run for BOTH stage-2 training populations
and compared side by side -- the NanoMVA-migration analogue of fom_scan_stage2d.py
(2026-09-16, pre-migration), generalized to compare "wjetsres" (trained behind
stage-1 eff_sig=0.90) against "wjetsres_nostage1" (trained on all of GenTag modes
1-5, no stage-1 pre-filter) rather than scanning just one of them.

WHY A SEPARATE SCRIPT, NOT AN EDIT OF fom_scan_stage2d.py: that script is a valid
historical record of the pre-migration WP1=50%/WP2=75% decision and still reads the
old sample names (WtoLNu2Jets_GenTag, Data, Lumi_reweight) -- editing it in place
would silently invalidate anyone reading it as "what we did on 2026-09-16".

SAMPLE NAMES (NanoMVA convention, 2026-09-17 migration):
    WtoLNu2Jets_merged (was WtoLNu2Jets_GenTag -- cand_gen_mode lives in this file
    directly now), Data_OS_merged (was Data), weight (was Lumi_reweight).

THE JOIN PROBLEM AND WHY IT DOES NOT NEED ONE: "wjetsres" reads its rows from
stage-1's own Prediction/ dir (so score_dsmu and score_dsmu_res already coexist in
the same file/row). "wjetsres_nostage1" has stage1_cut=None, so predict.py scored
it directly from SAMPLES_DIR -- its Prediction/ files carry score_dsmu_res but NOT
score_dsmu. Both this variant's Prediction/<name>.root and stage-1's own
Prediction/<name>.root are scored from the SAME underlying SAMPLES_DIR/<name>.root
with no row filtering before scoring (predict.py's stated design: "every row gets a
score"), so they have the same row count in the same on-disk order -- a positional
zip is enough, PROVIDED the (run, luminosityBlock, event) columns actually agree,
which this script asserts rather than assumes (this is exactly the kind of
"looks safe, verify anyway" case flagged in feedback_enumerate_all_paths_before_concluding).

Two background columns per cell (same as fom_scan_stage2d.py, same caveats):
    B_data_sb_OS  OS-sideband data survival (|M(KKpi)-M_Ds| >= 40 MeV) -- a
                  SIDEBAND PROXY for the SR background, not the SR background
                  itself; score cuts are known to reshape the mass spectrum.
    B_res         WtoLNu2Jets_merged modes {1,2,3,4,5}, weight-weighted, SR window
                  (|M(KKpi)-M_Ds| < 40 MeV) -- MC-only, but isolates what stage-2
                  is actually for. B_nonres (modes {0,6,7,8,9,10}) is also read for
                  the guard-2 exhaustiveness check but not folded into the FOM.

WP1 grid is capped at 0.90 for BOTH variants, even though wjetsres_nostage1 has no
training-population reason to stop there: the point of this script is a fair
side-by-side, and wjetsres genuinely cannot be evaluated above WP1=0.90 (out of its
training distribution), so extending wjetsres_nostage1 further would not give a
comparable cell.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 compare_stage2_wp_fom.py
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
SIGNAL_SAMPLE = common.SIGNALS  # ["Signal_merged"]

M_DS, SR_HALF = 1.96835, 0.040
NONRES = [0, 6, 7, 8, 9, 10]
RES = [1, 2, 3, 4, 5]

WP1_GRID = [round(0.90 - 0.05 * i, 2) for i in range(13)]   # 0.90 .. 0.30, step 0.05
WP2_GRID = [round(0.95 - 0.05 * i, 2) for i in range(10)]   # 0.95 .. 0.50, step 0.05
NEFF_RES_FLOOR = 50.0
NEFF_MC_FLOOR = 50.0
N_DATA_FLOOR = 10

IDS = ["run", "luminosityBlock", "event"]
SEL = [*common.SELECTION]


def sel_mask(a, has_score2=True):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if has_score2:
        m &= (a["score_dsmu_res"] > -1)
    return m


def neff(w):
    s2 = np.sum(w ** 2)
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


def load_joined(sample, stage2_key, extra_cols):
    """(score_dsmu, score_dsmu_res, extra columns), row-aligned, from stage-1's
    Prediction/<sample>.root and stage2_key's own Prediction/<sample>.root.

    Both files are scored from the same SAMPLES_DIR/<sample>.root with no
    pre-filtering, so they share row count and on-disk order -- asserted here via
    the (run, luminosityBlock, event) columns rather than trusted blindly.
    """
    with uproot.open(f"{STAGE1_PRED}/{sample}.root")["Events"] as t:
        a1 = t.arrays(["score_dsmu"] + IDS + SEL + extra_cols, library="np")

    s2_dir = f"{common.out_dir_for(stage2_key, tool=STAGE2_TOOL)}/Prediction"
    with uproot.open(f"{s2_dir}/{sample}.root")["Events"] as t:
        a2 = t.arrays(["score_dsmu_res"] + IDS, library="np")

    if len(a1["event"]) != len(a2["event"]):
        sys.exit(f"{sample}: stage-1 Prediction has {len(a1['event']):,} rows, "
                 f"{stage2_key} Prediction has {len(a2['event']):,} -- not row-aligned, "
                 f"a positional zip would silently pair the wrong rows.")
    for c in IDS:
        if not np.array_equal(a1[c], a2[c]):
            sys.exit(f"{sample}: column {c} disagrees between stage-1 and {stage2_key} "
                     f"Prediction files at some row -- not row-aligned.")

    out = dict(a1)
    out["score_dsmu_res"] = a2["score_dsmu_res"]
    return out


def run_scan(stage2_key):
    print(f"\n{'='*70}\n{stage2_key}\n{'='*70}")
    cfg = common.resolve(stage2_key)

    # ---------- guard 1 (only meaningful when stage-2 trained behind a stage-1 cut) ----------
    if cfg.get("stage1_cut"):
        thr_now, _ = common.stage1_threshold(cfg["stage1_cut"])
        s2_dir = common.out_dir_for(stage2_key, tool=STAGE2_TOOL)
        with open(f"{s2_dir}/stage1_filter.txt") as f:
            stamped = dict(line.split() for line in f)
        thr_stamped = float(stamped["thr"])
        assert abs(thr_now - thr_stamped) < 1e-6, (
            f"{stage2_key}: stage-1 threshold drift: stamped {thr_stamped:.6f} vs "
            f"recomputed {thr_now:.6f} -- stage-1 was re-scored after stage-2 trained "
            f"behind a different cut.")
        print(f"guard 1 OK: stage1_filter.txt thr={thr_stamped:.6f} matches recomputed")
    else:
        print("guard 1 skipped: stage1_cut=None, no training-time stage-1 threshold to check")

    # ---------- signal: truth-matched score_dsmu / score_dsmu_res, both channels ----------
    s1, s2, sw = [], [], []
    for s in SIGNAL_SAMPLE:
        a = load_joined(s, stage2_key, ["Ds_match", "mu_match", "weight"])
        m = sel_mask(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
        s1.append(a["score_dsmu"][m])
        s2.append(a["score_dsmu_res"][m])
        sw.append(a["weight"][m])
    s1 = np.concatenate(s1)
    s2 = np.concatenate(s2)
    sw = np.concatenate(sw)
    S_TOT = sw.sum()
    print(f"truth-matched signal: {len(s1):,} rows, weighted {S_TOT:.6g}")

    THR1 = {wp1: float(np.quantile(s1, 1 - wp1)) for wp1 in WP1_GRID}

    # ---------- background: WtoLNu2Jets_merged (resonant + non-resonant, cand_gen_mode) ----------
    b = load_joined(MC_SAMPLE, stage2_key, ["cand_gen_mode", "Ds_mass", "weight"])
    kb = sel_mask(b)
    b1, b2, b_mode = b["score_dsmu"][kb], b["score_dsmu_res"][kb], b["cand_gen_mode"][kb]
    b_mass, b_w = b["Ds_mass"][kb], b["weight"][kb]
    b_sr = np.abs(b_mass - M_DS) < SR_HALF
    b_nonres = np.isin(b_mode, NONRES)
    b_res = np.isin(b_mode, RES)
    print(f"WtoLNu2Jets_merged scored rows after selection: {kb.sum():,}")

    uniq, cnt = np.unique(b_mode, return_counts=True)
    stray = uniq[(~np.isin(uniq, NONRES + RES)) & (cnt > 0)]
    assert stray.size == 0, (
        f"cand_gen_mode values outside {{0..10}} at nonzero count: {stray.tolist()}")
    print(f"guard 2 OK: cand_gen_mode values present = {uniq.tolist()}, all in NONRES+RES")

    # ---------- OS sideband data ----------
    d = load_joined(DATA_SAMPLE, stage2_key, ["Ds_mass"])
    kd = sel_mask(d)
    d1, d2, d_mass = d["score_dsmu"][kd], d["score_dsmu_res"][kd], d["Ds_mass"][kd]
    d_sb = np.abs(d_mass - M_DS) >= SR_HALF
    print(f"OS data scored rows after selection: {kd.sum():,}\n")

    # ---------- scan ----------
    rows = []
    for wp1 in WP1_GRID:
        thr1 = THR1[wp1]
        s2_cond = s2[s1 > thr1]
        sw_cond = sw[s1 > thr1]
        for wp2 in WP2_GRID:
            thr2 = float(np.quantile(s2_cond, 1 - wp2))

            pass_sig_and = (s1 > thr1) & (s2 > thr2)
            eps_and = float(sw[pass_sig_and].sum() / S_TOT)
            eps1 = float(sw[s1 > thr1].sum() / S_TOT)
            eps2_cond = float(sw_cond[s2_cond > thr2].sum() / sw_cond.sum())
            eps_product = eps1 * eps2_cond

            pass_b = (b1 > thr1) & (b2 > thr2)
            w_res = b_w[pass_b & b_sr & b_res]
            w_nonres = b_w[pass_b & b_sr & b_nonres]
            B_res, B_nonres = float(w_res.sum()), float(w_nonres.sum())
            Neff_res = neff(w_res)
            Neff_mc = neff(np.concatenate([w_res, w_nonres]))
            B_mc_total = B_res + B_nonres

            pass_d = (d1 > thr1) & (d2 > thr2)
            B_data_sb = int(np.sum(pass_d & d_sb))

            vB_sb = np.sqrt(B_data_sb) if B_data_sb > 0 else np.nan
            fom_punzi3_sb = eps_and / (1.5 + vB_sb) if not np.isnan(vB_sb) else np.nan
            fom_naive_sb = eps_and / vB_sb if vB_sb > 0 else np.nan
            vB_mc = np.sqrt(B_mc_total) if B_mc_total > 0 else np.nan
            fom_punzi3_mc = eps_and / (1.5 + vB_mc) if not np.isnan(vB_mc) else np.nan
            fom_naive_mc = eps_and / vB_mc if vB_mc > 0 else np.nan

            rows.append(dict(wp1=wp1, wp2=wp2, thr1=thr1, thr2=thr2,
                             eps_and=eps_and, eps_product=eps_product,
                             B_res=B_res, Neff_res=Neff_res, B_nonres=B_nonres,
                             B_mc_total=B_mc_total, Neff_mc=Neff_mc,
                             B_data_sb=B_data_sb,
                             fom_naive_sb=fom_naive_sb, fom_punzi3_sb=fom_punzi3_sb,
                             fom_naive_mc=fom_naive_mc, fom_punzi3_mc=fom_punzi3_mc))

    worst = max(rows, key=lambda r: abs(r["eps_and"] - r["eps_product"]))
    diff = abs(worst["eps_and"] - worst["eps_product"])
    assert diff < 5e-3, (
        f"eps_total(AND) vs eps_total(product) disagree by {diff:.4f} at "
        f"WP1={worst['wp1']}, WP2={worst['wp2']}")
    print(f"guard 3 OK: max |eps_and - eps_product| = {diff:.5f} across the grid\n")

    hdr = (f"{'WP1':>5} {'WP2':>5} {'eps_tot':>8} "
           f"{'B_sb':>8} {'FOMsb_p3':>9} {'B_mc':>8} {'FOMmc_p3':>9} "
           f"{'B_res':>8} {'Neff_res':>9} {'Neff_mc':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        flag = " <-- Neff_res low" if r["Neff_res"] < NEFF_RES_FLOOR else ""
        flag += " <-- data stat low" if r["B_data_sb"] < N_DATA_FLOOR else ""
        print(f"{r['wp1']:>5.2f} {r['wp2']:>5.2f} {r['eps_and']:>8.4f} "
              f"{r['B_data_sb']:>8d} {r['fom_punzi3_sb']:>9.5f} "
              f"{r['B_mc_total']:>8.1f} {r['fom_punzi3_mc']:>9.5f} "
              f"{r['B_res']:>8.1f} {r['Neff_res']:>9.0f} {r['Neff_mc']:>8.0f}{flag}")

    return rows


def grid_of(rows, key):
    Z = np.full((len(WP1_GRID), len(WP2_GRID)), np.nan)
    for r in rows:
        Z[WP1_GRID.index(r["wp1"]), WP2_GRID.index(r["wp2"])] = r[key]
    return Z


def draw_panel(ax_, fom, low_stat, eps_tot, title):
    disp = np.where(low_stat, np.nan, fom)
    im = ax_.imshow(disp, origin="lower", aspect="auto", cmap="viridis")
    for i in range(len(WP1_GRID)):
        for j in range(len(WP2_GRID)):
            if low_stat[i, j]:
                ax_.text(j, i, "low\nstat", ha="center", va="center", fontsize=8,
                         color="grey", style="italic")
                continue
            color = "white" if fom[i, j] < 0.6 * np.nanmax(disp) else "black"
            ax_.text(j, i - 0.15, f"{fom[i, j]:.4f}", ha="center", va="center",
                     fontsize=13, color=color, fontweight="bold")
            ax_.text(j, i + 0.20, r"$\varepsilon_{tot}$=" + f"{eps_tot[i, j]:.2f}",
                     ha="center", va="center", fontsize=8.5, color=color)
    ax_.set_xticks(range(len(WP2_GRID)))
    ax_.set_xticklabels([f"{w:.2f}" for w in WP2_GRID], fontsize=9, rotation=90)
    ax_.set_yticks(range(len(WP1_GRID)))
    ax_.set_yticklabels([f"{w:.2f}" for w in WP1_GRID], fontsize=9)
    ax_.set_xlabel(r"WP2 $\varepsilon_2$ (conditional on WP1)")
    ax_.set_ylabel(r"WP1 $\varepsilon_1$")
    ax_.set_title(title, fontsize=13)
    return im


ALL_ROWS = {}
for key in STAGE2_KEYS:
    ALL_ROWS[key] = run_scan(key)

# ---------- figures: one file per (stage2 variant, FOM definition) ----------
FOM_PANELS = (
    ("epsB_sb", "fom_naive_sb", "B_data_sb", N_DATA_FLOOR,
     r"$\varepsilon_{sig}/\sqrt{B}$, B = B$_{\rm data\ sb}^{\rm OS}$ (sideband proxy)"),
    ("epsB_mc", "fom_naive_mc", "Neff_mc", NEFF_MC_FLOOR,
     r"$\varepsilon_{sig}/\sqrt{B}$, B = B$_{\rm res}$+B$_{\rm nonres}$ (GenTag MC)"),
    ("punzi3_sb", "fom_punzi3_sb", "B_data_sb", N_DATA_FLOOR,
     r"$\varepsilon_{sig}/(1.5+\sqrt{B})$, B = B$_{\rm data\ sb}^{\rm OS}$ (sideband proxy)"),
    ("punzi3_mc", "fom_punzi3_mc", "Neff_mc", NEFF_MC_FLOOR,
     r"$\varepsilon_{sig}/(1.5+\sqrt{B})$, B = B$_{\rm res}$+B$_{\rm nonres}$ (GenTag MC)"),
)

for key in STAGE2_KEYS:
    rows = ALL_ROWS[key]
    eps_tot = grid_of(rows, "eps_and")
    for tag, fom_key, stat_key, floor, title in FOM_PANELS:
        Z = grid_of(rows, fom_key)
        low_stat = grid_of(rows, stat_key) < floor
        fig, ax_ = plt.subplots(figsize=(13, 11))
        im = draw_panel(ax_, Z, low_stat, eps_tot, f"{title}\n({key})")
        fig.colorbar(im, ax=ax_, label="FOM")
        hep.cms.label(ax=ax_, data=False, text="Preliminary", fontsize=11)
        plt.tight_layout()
        out = common.fig_path("comparisons", f"fom_scan_stage2d_{tag}.png",
                              key=key, tool=STAGE2_TOOL)
        plt.savefig(out, dpi=120)
        plt.close(fig)
        print(f"Saved {out}")

# ---------- direct comparison: same-cell FOM ratio, wjetsres_nostage1 / wjetsres ----------
# Only meaningful where BOTH variants have usable statistics in that cell.
print(f"\n{'='*70}\nwjetsres_nostage1 / wjetsres, same (WP1,WP2) cell\n{'='*70}")
a_rows, b_rows = ALL_ROWS["wjetsres"], ALL_ROWS["wjetsres_nostage1"]
hdr = f"{'WP1':>5} {'WP2':>5} {'FOMsb_p3 a':>11} {'FOMsb_p3 b':>11} {'ratio b/a':>10} {'FOMmc_p3 a':>11} {'FOMmc_p3 b':>11} {'ratio b/a':>10}"
print(hdr)
print("-" * len(hdr))
for ra, rb in zip(a_rows, b_rows):
    assert ra["wp1"] == rb["wp1"] and ra["wp2"] == rb["wp2"]
    fa_sb, fb_sb = ra["fom_punzi3_sb"], rb["fom_punzi3_sb"]
    fa_mc, fb_mc = ra["fom_punzi3_mc"], rb["fom_punzi3_mc"]
    r_sb = fb_sb / fa_sb if (fa_sb and not np.isnan(fa_sb) and not np.isnan(fb_sb)) else np.nan
    r_mc = fb_mc / fa_mc if (fa_mc and not np.isnan(fa_mc) and not np.isnan(fb_mc)) else np.nan
    print(f"{ra['wp1']:>5.2f} {ra['wp2']:>5.2f} {fa_sb:>11.5f} {fb_sb:>11.5f} {r_sb:>10.3f} "
          f"{fa_mc:>11.5f} {fb_mc:>11.5f} {r_mc:>10.3f}")


# ---------- background survival at MATCHED signal efficiency, split res / non-res ----------
# The FOM table above folds efficiency and background into one number. This one answers
# the more basic question directly: at the same (WP1, WP2) cell -- hence the same nominal
# signal efficiency -- how much background is LEFT, and is the difference in the resonant
# part (real D/Ds, what stage-2 exists for) or the non-resonant part (combinatorial,
# which is stage-1's job)?
#
# CAVEAT on "matched": thr2 is an UNWEIGHTED quantile of the signal score while eps_tot is
# WEIGHTED, so the realized eps_tot of the two arms differs slightly in the same cell. The
# eps_tot columns are printed side by side so the residual mismatch is visible rather than
# assumed; measured median 1.4e-3 relative, which propagates to ~0.4% on B (local
# dlnB/dlnE ~ 2). Ratios far from 1 survive that; ratios within ~1% of 1 do not.
print(f"\n{'='*70}\nBackground left at the same cell, b=wjetsres_nostage1 vs a=wjetsres"
      f"\n(b/a < 1 means b leaves LESS background = rejects more)\n{'='*70}")
hdr2 = (f"{'WP1':>5} {'WP2':>5} {'eps_a':>7} {'eps_b':>7} | {'Bres a':>8} {'Bres b':>8} {'b/a':>6} | "
        f"{'Bnres a':>8} {'Bnres b':>8} {'b/a':>6} | {'Bsb a':>7} {'Bsb b':>7} {'b/a':>6} | "
        f"{'res% a':>6} {'res% b':>6} {'Neff':>5}")
print(hdr2)
print("-" * len(hdr2))
for ra, rb in zip(a_rows, b_rows):
    anr, bnr = ra["B_nonres"], rb["B_nonres"]
    q = lambda x, y: (y / x) if x else np.nan
    low = " *" if min(ra["Neff_res"], rb["Neff_res"]) < NEFF_RES_FLOOR else ""
    print(f"{ra['wp1']:>5.2f} {ra['wp2']:>5.2f} {ra['eps_and']:>7.4f} {rb['eps_and']:>7.4f} | "
          f"{ra['B_res']:>8.1f} {rb['B_res']:>8.1f} {q(ra['B_res'], rb['B_res']):>6.3f} | "
          f"{anr:>8.1f} {bnr:>8.1f} {q(anr, bnr):>6.3f} | "
          f"{ra['B_data_sb']:>7d} {rb['B_data_sb']:>7d} {q(ra['B_data_sb'], rb['B_data_sb']):>6.3f} | "
          f"{100 * ra['B_res'] / ra['B_mc_total']:>5.1f}% {100 * rb['B_res'] / rb['B_mc_total']:>5.1f}% "
          f"{min(ra['Neff_res'], rb['Neff_res']):>5.0f}{low}")

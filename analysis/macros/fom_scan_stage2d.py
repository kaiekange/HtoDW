"""2D (WP1, WP2) working-point scan for the two-stage MVA (follow-up to fom_scan_stage1.py).

Stage-1 (score_dsmu, DNN 'ssdata') rejects NON-resonant background. Stage-2
(score_dsmu_res, BDT 'wjetsres') rejects RESONANT background (real Ds/D in
GenTag modes 1-5) and was trained on rows that already pass stage-1 at
eff_sig=0.90 -- it has never seen rows below that cut. The grid is therefore
restricted to WP1 <= 0.90 (a subset of the training population, safe); WP1
looser than 0.90 would score stage-2 out of its training distribution and is
not attempted here.

WP2 is a CONDITIONAL quantile: thr2(WP1, WP2) is computed from the
truth-matched signal that already passes WP1's cut, so eps_total = eps1 * eps2
by construction and the AND of both cuts reproduces the analysis selection.
This means "WP2 = 0.70" is a DIFFERENT score_dsmu_res threshold at every WP1 --
read the printed thr2, not the WP2 label, as "the cut".

Two background columns, both first-class (2026-09-16: the user's stated goal
is high eps_sig / low OS data background, but B_res is kept alongside because
the OS sideband barely moves under WP2 -- real Ds/D from W+jets sits under the
peak, mostly inside the blind SR, so B_res is where WP2's payoff actually shows
up):
  B_data_sb_OS   OS-sideband data survival count (|M(KKpi)-M_Ds| >= 40 MeV,
                 OS is SR-blind anyway). This is a SIDEBAND PROXY for the SR
                 background, not the SR background itself -- prior sessions
                 established the score cuts reshape the mass spectrum (low/high
                 mass double ratio 1.02 -> 1.29 @ eps90 -> 1.78 @ eps70), so a
                 cell that minimizes this count is not guaranteed to minimize
                 the actual (blind) SR background. Extrapolation to SR is NOT
                 established here -- do not read this column as "the" answer.
  B_res          WtoLNu2Jets_GenTag modes {1,2,3,4,5}, Lumi_reweight-weighted,
                 SR mass window (|M(KKpi)-M_Ds| < 40 MeV). MC-only (no truth
                 tag exists on data), but it is the one column that isolates
                 what stage-2 is actually for.
  B_nonres       GenTag modes {0,6,7,8,9,10}, same weighting/window. Printed
                 for context (compare to B_data_sb_OS as a normalization
                 cross-check, same convention as fom_scan_stage1.py) -- NOT
                 summed with B_data_sb_OS (that would double-count: modes 1-5
                 are physically inside the data sideband already).

No single FOM is computed -- the user has not decided whether Punzi is the
right quantification yet. This prints the full grid and draws eps_total
contours over each background column so the trade-off can be read directly.

Three guards before trusting the table (see project CLAUDE.md Rule 8):
  1. stage1_filter.txt in the stage-2 dir is asserted to match
     stage1_threshold() recomputed now -- catches stage-1 having been
     re-scored after stage-2 trained behind a different cut.
  2. cand_gen_mode values are asserted to be a subset of {0..10} with nothing
     else present at nonzero count -- catches a non-exhaustive mode partition
     silently dropping a background component from B_nonres + (nothing).
  3. eps_total is computed two ways (product of marginals vs direct AND) and
     the two are asserted to agree -- catches a sentinel/mask bug in either
     score branch (score == -777 rows must be excluded from BOTH, not just
     one).

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 fom_scan_stage2d.py
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

STAGE1_KEY, STAGE2_KEY = "ssdata", "wjetsres"
STAGE2_TOOL = "BDT"
PRED = f"{common.out_dir_for(STAGE2_KEY, tool=STAGE2_TOOL)}/Prediction"
M_DS, SR_HALF = 1.96835, 0.040
NONRES = [0, 6, 7, 8, 9, 10]
RES = [1, 2, 3, 4, 5]

WP1_GRID = [round(0.90 - 0.05 * i, 2) for i in range(13)]   # 0.90 .. 0.30, step 0.05
WP2_GRID = [round(0.95 - 0.05 * i, 2) for i in range(10)]   # 0.95 .. 0.50, step 0.05
NEFF_RES_FLOOR = 50.0     # readability floor for THIS grid, NOT the 1500 training
                           # floor in fom_scan_stage1.py -- that one sizes the BDT's
                           # training set; this one is "is sqrt(B_res) still a number
                           # or is it noise" for a single scan cell.
NEFF_MC_FLOOR = 50.0      # same idea, for the combined B_res+B_nonres column
N_DATA_FLOOR = 10         # raw (unweighted) OS-sideband data count floor per cell

SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "score_dsmu", "score_dsmu_res"]


def sel(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
            & (a[common.B_VETO] == 0)
            & (a["score_dsmu"] > -1) & (a["score_dsmu_res"] > -1))


def neff(w):
    s2 = np.sum(w ** 2)
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


# ---------- guard 1: stage-1 threshold behind stage-2's training hasn't moved ----------
stage1_cut = common.resolve(STAGE2_KEY)["stage1_cut"]
thr_now, _ = common.stage1_threshold(stage1_cut)
with open(f"{common.out_dir_for(STAGE2_KEY, tool=STAGE2_TOOL)}/stage1_filter.txt") as f:
    stamped = dict(line.split() for line in f)
thr_stamped = float(stamped["thr"])
assert abs(thr_now - thr_stamped) < 1e-6, (
    f"stage-1 threshold drift: stamped {thr_stamped:.6f} vs recomputed {thr_now:.6f} -- "
    f"stage-1 was re-scored after stage-2 trained behind a different cut; the grid's "
    f"WP1=0.90 cell would not match stage-2's training population.")
print(f"guard 1 OK: stage1_filter.txt thr={thr_stamped:.6f} matches recomputed thr={thr_now:.6f}")

# ---------- signal: truth-matched score_dsmu / score_dsmu_res, both eras ----------
s1, s2, sw = [], [], []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "score_dsmu_res", "Ds_match", "mu_match", "Lumi_reweight"] + SEL,
        library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s1.append(a["score_dsmu"][m])
    s2.append(a["score_dsmu_res"][m])
    sw.append(a["Lumi_reweight"][m])
s1 = np.concatenate(s1)
s2 = np.concatenate(s2)
sw = np.concatenate(sw)
S_TOT = sw.sum()
print(f"truth-matched signal: {len(s1):,} rows, weighted {S_TOT:,.1f}")

THR1 = {wp1: float(np.quantile(s1, 1 - wp1)) for wp1 in WP1_GRID}

# ---------- background: GenTag (resonant + non-resonant, same tree) ----------
b = uproot.open(f"{PRED}/WtoLNu2Jets_GenTag.root")["Events"].arrays(
    ["score_dsmu", "score_dsmu_res", "cand_gen_mode", "Ds_mass", "Lumi_reweight"] + SEL,
    library="np")
kb = sel(b)
b1, b2, b_mode, b_mass, b_w = (b["score_dsmu"][kb], b["score_dsmu_res"][kb],
                               b["cand_gen_mode"][kb], b["Ds_mass"][kb],
                               b["Lumi_reweight"][kb])
b_sr = np.abs(b_mass - M_DS) < SR_HALF
b_nonres = np.isin(b_mode, NONRES)
b_res = np.isin(b_mode, RES)
print(f"GenTag scored rows after selection: {kb.sum():,}")

# ---------- guard 2: mode partition is exhaustive ----------
uniq, cnt = np.unique(b_mode, return_counts=True)
stray = uniq[(~np.isin(uniq, NONRES + RES)) & (cnt > 0)]
assert stray.size == 0, (
    f"cand_gen_mode values outside {{0..10}} at nonzero count: {stray.tolist()} -- "
    f"B_nonres + B_res would silently drop this component.")
print(f"guard 2 OK: cand_gen_mode values present = {uniq.tolist()}, all in NONRES+RES")

# ---------- cross-check: OS sideband data ----------
d = uproot.open(f"{PRED}/Data.root")["Events"].arrays(
    ["score_dsmu", "score_dsmu_res", "Ds_mass"] + SEL, library="np")
kd = sel(d)
d1, d2, d_mass = d["score_dsmu"][kd], d["score_dsmu_res"][kd], d["Ds_mass"][kd]
d_sb = np.abs(d_mass - M_DS) >= SR_HALF
print(f"OS data scored rows after selection: {kd.sum():,}\n")

# ---------- scan ----------
rows = []
for wp1 in WP1_GRID:
    thr1 = THR1[wp1]
    # WP2 thresholds are conditional on WP1: quantiles of score_dsmu_res among
    # truth-matched signal that already passes thr1.
    s2_cond = s2[s1 > thr1]
    sw_cond = sw[s1 > thr1]
    for wp2 in WP2_GRID:
        thr2 = float(np.quantile(s2_cond, 1 - wp2))

        # eps_total, two independent ways
        pass_sig_and = (s1 > thr1) & (s2 > thr2)
        eps_and = float(sw[pass_sig_and].sum() / S_TOT)
        eps1 = float(sw[s1 > thr1].sum() / S_TOT)
        eps2_cond = float(sw_cond[s2_cond > thr2].sum() / sw_cond.sum())
        eps_product = eps1 * eps2_cond

        pass_b = (b1 > thr1) & (b2 > thr2)
        w_res = b_w[pass_b & b_sr & b_res]
        w_nonres = b_w[pass_b & b_sr & b_nonres]
        B_res = float(w_res.sum())
        B_nonres = float(w_nonres.sum())
        Neff_res = neff(w_res)
        Neff_mc = neff(np.concatenate([w_res, w_nonres]))
        B_mc_total = B_res + B_nonres   # same tree, disjoint modes -- see guard 2

        pass_d = (d1 > thr1) & (d2 > thr2)
        B_data_sb = int(np.sum(pass_d & d_sb))

        # Two FOM definitions, same eps_and numerator, different background:
        #   *_sb  : B_data_sb_OS, the sideband proxy the user asked for.
        #   *_mc  : B_res + B_nonres, the full GenTag MC background in the SR
        #           window -- includes the peaking component the sideband can't see.
        vB_sb = np.sqrt(B_data_sb) if B_data_sb > 0 else np.nan
        fom_naive_sb = eps_and / vB_sb if vB_sb > 0 else np.nan
        fom_punzi3_sb = eps_and / (1.5 + vB_sb) if not np.isnan(vB_sb) else np.nan
        vB_mc = np.sqrt(B_mc_total) if B_mc_total > 0 else np.nan
        fom_naive_mc = eps_and / vB_mc if vB_mc > 0 else np.nan
        fom_punzi3_mc = eps_and / (1.5 + vB_mc) if not np.isnan(vB_mc) else np.nan

        rows.append(dict(wp1=wp1, wp2=wp2, thr1=thr1, thr2=thr2,
                         eps_and=eps_and, eps_product=eps_product,
                         B_res=B_res, Neff_res=Neff_res, B_nonres=B_nonres,
                         B_mc_total=B_mc_total, Neff_mc=Neff_mc,
                         B_data_sb=B_data_sb,
                         fom_naive_sb=fom_naive_sb, fom_punzi3_sb=fom_punzi3_sb,
                         fom_naive_mc=fom_naive_mc, fom_punzi3_mc=fom_punzi3_mc))

# ---------- guard 3: the two eps_total computations agree ----------
worst = max(rows, key=lambda r: abs(r["eps_and"] - r["eps_product"]))
diff = abs(worst["eps_and"] - worst["eps_product"])
assert diff < 5e-3, (
    f"eps_total(AND) vs eps_total(product) disagree by {diff:.4f} at "
    f"WP1={worst['wp1']}, WP2={worst['wp2']} ({worst['eps_and']:.4f} vs "
    f"{worst['eps_product']:.4f}) -- likely a sentinel/mask mismatch between "
    f"score_dsmu and score_dsmu_res.")
print(f"guard 3 OK: max |eps_and - eps_product| = {diff:.5f} across the grid\n")

# ---------- print table ----------
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

# ---------- figure: two FOM definitions side by side, same eps_tot numerator ----------
# Left: B_data_sb_OS (sideband proxy, what the user asked for). Right: B_res+B_nonres
# (full GenTag MC background in the SR window, including the peaking component the
# sideband can't see). Put side by side deliberately -- see the 2026-09-16 exchange:
# the sideband-based FOM is monotone AGAINST tightening WP2 (it can't see what WP2
# kills), so if the two panels disagree, that disagreement -- not either number
# alone -- is the finding, and it maps onto the still-open question of whether the
# final background comes from a data-sideband fit or an MC shape.
wp1s = np.array(WP1_GRID)
wp2s = np.array(WP2_GRID)


def grid_of(key):
    Z = np.full((len(wp1s), len(wp2s)), np.nan)
    for r in rows:
        Z[WP1_GRID.index(r["wp1"]), WP2_GRID.index(r["wp2"])] = r[key]
    return Z


FOM_SB_P3, FOM_SB_NV = grid_of("fom_punzi3_sb"), grid_of("fom_naive_sb")
FOM_MC_P3, FOM_MC_NV = grid_of("fom_punzi3_mc"), grid_of("fom_naive_mc")
BDATA, NEFF_MC = grid_of("B_data_sb"), grid_of("Neff_mc")
EPS_TOT = grid_of("eps_and")   # combined (WP1 AND WP2) signal efficiency, same in every panel

low_stat_sb = BDATA < N_DATA_FLOOR
low_stat_mc = NEFF_MC < NEFF_MC_FLOOR


def draw_panel(ax_, fom, low_stat, title):
    disp = np.where(low_stat, np.nan, fom)
    im = ax_.imshow(disp, origin="lower", aspect="auto", cmap="viridis")
    for i in range(len(wp1s)):
        for j in range(len(wp2s)):
            if low_stat[i, j]:
                ax_.text(j, i, "low\nstat", ha="center", va="center", fontsize=8,
                         color="grey", style="italic")
                continue
            color = "white" if fom[i, j] < 0.6 * np.nanmax(disp) else "black"
            ax_.text(j, i - 0.15, f"{fom[i, j]:.4f}", ha="center", va="center",
                     fontsize=13, color=color, fontweight="bold")
            ax_.text(j, i + 0.20, r"$\varepsilon_{tot}$=" + f"{EPS_TOT[i, j]:.2f}",
                     ha="center", va="center", fontsize=8.5, color=color)
    ax_.set_xticks(range(len(wp2s)))
    ax_.set_xticklabels([f"{w:.2f}" for w in wp2s], fontsize=9, rotation=90)
    ax_.set_yticks(range(len(wp1s)))
    ax_.set_yticklabels([f"{w:.2f}" for w in wp1s], fontsize=9)
    ax_.set_xlabel(r"WP2 $\varepsilon_2$ (conditional on WP1)")
    ax_.set_ylabel(r"WP1 $\varepsilon_1$")
    ax_.set_title(title, fontsize=13)
    return im


# Four separate figures (2026-09-16: one file per FOM definition, not four subplots
# in one file) -- title states the formula actually used, not a name for it.
FOM_PANELS = (
    ("epsB_sb", FOM_SB_NV, low_stat_sb,
     r"$\varepsilon_{sig}/\sqrt{B}$, B = B$_{\rm data\ sb}^{\rm OS}$ (sideband proxy)"),
    ("epsB_mc", FOM_MC_NV, low_stat_mc,
     r"$\varepsilon_{sig}/\sqrt{B}$, B = B$_{\rm res}$+B$_{\rm nonres}$ (GenTag MC)"),
    ("punzi3_sb", FOM_SB_P3, low_stat_sb,
     r"$\varepsilon_{sig}/(1.5+\sqrt{B})$, B = B$_{\rm data\ sb}^{\rm OS}$ (sideband proxy)"),
    ("punzi3_mc", FOM_MC_P3, low_stat_mc,
     r"$\varepsilon_{sig}/(1.5+\sqrt{B})$, B = B$_{\rm res}$+B$_{\rm nonres}$ (GenTag MC)"),
)
for tag, Z, low_stat, title in FOM_PANELS:
    fig, ax_ = plt.subplots(figsize=(13, 11))
    im = draw_panel(ax_, Z, low_stat, title)
    fig.colorbar(im, ax=ax_, label="FOM")
    hep.cms.label(ax=ax_, data=False, text="Simulation", fontsize=11)
    plt.tight_layout()
    out = common.fig_path("comparisons", f"fom_scan_stage2d_{tag}.png",
                          key=STAGE2_KEY, tool=STAGE2_TOOL)
    plt.savefig(out, dpi=120)
    plt.close(fig)
    print(f"Saved {out}")

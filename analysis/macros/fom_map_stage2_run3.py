"""2D (WP1, WP2) Punzi FOM maps for the two-stage MVA, scaled to Run-3 (326 /fb).

    FOM = eps_tot / (1.5 + sqrt(B))        Punzi a=3

Punzi is BR-INDEPENDENT, which is why it is the right choice: the signal
normalisation carries whatever BR(H->DsW) the sample was generated with, so any
S/sqrt(B) map would inherit that assumption.

2026-10-26: B_wjets -> B_mc = W+jets + DY + tt (each x br_reweight, DY also x
dy_jetbin_sf), all modes, counted ONLY in the Ds SR |M(KKpi) - M_Ds| < 40 MeV (user,
1026): the full range was 70% outside the SR at the optimum (D+- peak, misID, sidebands are
separate fit components), so it rewarded stage-2 for removing things that do not sit under
the Ds peak. In the SR, MC NeffB is small (~10-20 near the optimum) -- read NeffB per cell. Takes one or more
stage-2 keys (default resdy); each key gets its own maps (never merged into one figure),
and with two keys the best cell of each is printed side by side. Each key's eps_tot is its own (same WP2 = same unconditional
stage-2 signal efficiency, not exactly the same eps_tot), printed in every cell.

    python3 fom_map_stage2_run3.py resdy wjetsres_v2

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

# --lo=X: lowest WP1/WP2 in the grid (default 0.30, step 0.05); file names get _loX
LO = next((float(a.split("=", 1)[1]) for a in sys.argv[1:] if a.startswith("--lo=")), 0.30)
KEYS = [a for a in sys.argv[1:] if not a.startswith("--")] or ["resdy"]
assert len(KEYS) in (1, 2), "one stage-2 key, or two to print their best cells side by side"
STAGE2_TOOL = "BDT"

MC_SAMPLES = ["WtoLNu2Jets_merged", "DY_merged", "TT_merged"]
DATA_SAMPLE = "Data_OS_merged"

M_DS, SR_HALF = 1.96835, 0.040
PUNZI_A = 3.0

LUMI_NOW_PB = sum(e["lumi_pb"] for e in eras.ERAS.values())     # 2022 preEE+postEE
LUMI_RUN3_PB = 326000.0
SCALE = LUMI_RUN3_PB / LUMI_NOW_PB

WP1_GRID = [round(0.90 - 0.05 * i, 2) for i in range(int(round((0.90 - LO) / 0.05)) + 1)]   # 0.90 .. LO
WP2_GRID = [round(0.95 - 0.05 * i, 2) for i in range(int(round((0.95 - LO) / 0.05)) + 1)]   # 0.95 .. LO
LO_TAG = "" if abs(LO - 0.30) < 1e-9 else f"_lo{LO:.2f}"
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


def load(pred, sample, extra):
    """Both scores live in stage-2's own Prediction file: `wjetsres` is scored from
    stage-1's Prediction dir, so score_dsmu rides along. No positional join needed
    (unlike wjetsres_nostage1, which is scored straight from SAMPLES_DIR)."""
    cols = ["score_dsmu", "score_dsmu_res"] + SEL + extra
    return uproot.open(f"{pred}/{sample}.root")["Events"].arrays(cols, library="np")


def scan(key):
    """All rows of the (WP1, WP2) grid for one stage-2 key."""
    pred = f"{common.out_dir_for(key, tool=STAGE2_TOOL)}/Prediction"
    # ---------- guard 1: the stage-1 cut stage-2 trained behind is still the same ----------
    cfg = common.resolve(key)
    thr_now, _ = common.stage1_threshold(cfg["stage1_cut"])
    with open(f"{common.out_dir_for(key, tool=STAGE2_TOOL)}/stage1_filter.txt") as f:
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
        a = load(pred, s, ["Ds_match", "mu_match", "weight"])
        m = sel_mask(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
        s1.append(a["score_dsmu"][m])
        s2.append(a["score_dsmu_res"][m])
        sw.append(a["weight"][m])
    s1, s2, sw = np.concatenate(s1), np.concatenate(s2), np.concatenate(sw)
    S_TOT = sw.sum()
    print(f"truth-matched signal: {len(s1):,} rows, weighted {S_TOT:.4f} events "
          f"-> {S_TOT*SCALE:.3f} at Run-3 (context only; Punzi is BR-independent)")
    THR1 = {wp1: float(np.quantile(s1, 1 - wp1)) for wp1 in WP1_GRID}

    # ---------- background A: W+jets + DY + tt in the Ds SR (all modes) ----------
    b1, b2, b_w = [], [], []
    for smp in MC_SAMPLES:
        extra = ["cand_gen_mode", "cand_gen_common_hf_pdgId", "weight", "Ds_mass"] + (["jet_bin"] if smp.startswith("DY") else [])
        b = load(pred, smp, extra)
        kb = sel_mask(b) & (np.abs(b["Ds_mass"] - M_DS) < SR_HALF)
        b = {k: v[kb] for k, v in b.items()}
        w = b["weight"].astype(np.float64) * common.br_reweight(b)
        if smp.startswith("DY"):
            w = w * common.dy_jetbin_sf(b)
        uniq = np.unique(b["cand_gen_mode"])
        assert np.isin(uniq, np.arange(0, 11)).all(), f"{smp}: cand_gen_mode outside 0..10: {uniq}"
        print(f"{smp}: {kb.sum():,} rows after selection, in the Ds SR, ALL modes enter B_mc, "
              f"weighted {w.sum()*SCALE:,.1f} at Run-3 (x BR rw" + (" x DY SF)" if smp.startswith("DY") else ")"))
        b1.append(b["score_dsmu"]); b2.append(b["score_dsmu_res"]); b_w.append(w)
    b_src = np.concatenate([np.full(len(x), i) for i, x in enumerate(b1)])
    b1, b2, b_w = np.concatenate(b1), np.concatenate(b2), np.concatenate(b_w)

    # ---------- background B: OS data control region (SR left blind) ----------
    d = load(pred, DATA_SAMPLE, ["Ds_mass"])
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

            pb = p1b & (b2 > thr2)
            w_sr = b_w[pb]
            n_cr = int(np.sum(p1d & (d2 > thr2) & d_cr))
            rows.append(dict(
                wp1=wp1, wp2=wp2, eps_and=eps_and, eps2=eps2,
                B_w=float(w_sr.sum()) * SCALE, raw_w=int(w_sr.size), neff_w=neff(w_sr),
                n_cr=n_cr, B_d=n_cr * SCALE,
                B_parts=[float(b_w[pb & (b_src == i)].sum()) * SCALE for i in range(len(MC_SAMPLES))]))
    for r in rows:
        r["fom_w"] = punzi(r["eps_and"], r["B_w"])
        r["fom_d"] = punzi(r["eps_and"], r["B_d"])

    worst = max(rows, key=lambda r: abs(r["eps2"] - r["wp2"]))
    diff = abs(worst["eps2"] - worst["wp2"])
    assert diff < 5e-3, (f"WP2 quantile off by {diff:.4f} at WP2={worst['wp2']}: the "
                         f"threshold is not delivering its nominal unconditional efficiency")
    print(f"guard 3 OK: max |eps_2 - WP2| = {diff:.5f} across the grid "
          f"(WP2 is unconditional)\n")
    return rows


RESULTS = {}
for KEY in KEYS:
    print(f"\n========== stage-2 {KEY} ==========")
    rows = scan(KEY)
    RESULTS[KEY] = rows
    hdr = (f"{'WP1':>5} {'WP2':>5} {'eps_tot':>8} | {'B_mc':>9} {'rawB':>6} {'NeffB':>6} "
           f"{'FOM_w':>8} | {'n_CR':>6} {'B_data':>10} {'FOM_d':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['wp1']:>5.2f} {r['wp2']:>5.2f} {r['eps_and']:>8.4f} | "
              f"{r['B_w']:>9.1f} {r['raw_w']:>6d} {r['neff_w']:>6.0f} {r['fom_w']:>8.5f} | "
              f"{r['n_cr']:>6d} {r['B_d']:>10.1f} {r['fom_d']:>8.5f}")

    for key, name in (("fom_w", "B = W+jets+DY+tt in the Ds SR (all modes)"),
                      ("fom_d", "B = OS data (control region)")):
        bst = max(rows, key=lambda r: (r[key] if not np.isnan(r[key]) else -np.inf))
        print(f"\n{name}: max FOM {bst[key]:.5f} at WP1={bst['wp1']:.2f}, WP2={bst['wp2']:.2f} "
              f"(eps_tot={bst['eps_and']:.4f}, B={bst[key.replace('fom', 'B')]:.1f}, "
              f"rawB/n_CR={bst['raw_w'] if key == 'fom_w' else bst['n_cr']}"
              + (f", B W+jets/DY/tt = {'/'.join(f'{x:.1f}' for x in bst['B_parts'])}" if key == 'fom_w' else "") + ")")


# ---------- figures: one per background ----------
def grid_of(rows, key):
    return np.array([[next(r[key] for r in rows if r["wp1"] == w1 and r["wp2"] == w2)
                      for w2 in WP2_GRID] for w1 in WP1_GRID])


RLABEL = rf"{LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$ (13.6 TeV)"


for KEY, rows in RESULTS.items():
    eps_tot = grid_of(rows, "eps_and")
    PANELS = (
        ("mcsr", "fom_w", False,
         r"B = W+jets + DY + $t\bar{t}$ in the D$_s$ SR ($|M-M_{D_s}|<40$ MeV; x BR rw, DY x jet-bin SF), all modes"),
        ("data", "fom_d", True,
         r"B = OS data, control region ($|M-M_{D_s}|\geq 40$ MeV, SR blind)"),
    )
    for tag, key, is_data, subtitle in PANELS:
        Z = grid_of(rows, key)
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
        ax.set_title(rf"stage-2 {KEY}: " r"Punzi FOM $=\varepsilon_{tot}/(1.5+\sqrt{B})$, "
                     rf"scaled to {LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$" "\n" + subtitle,
                     fontsize=13)
        fig.colorbar(im, ax=ax, label="Punzi FOM")
        # ALWAYS text="Preliminary"; `data=` alone decides the "Simulation" prefix, giving
        # "CMS Preliminary" / "CMS Simulation Preliminary". Never text="Simulation" (that
        # doubles it) and never drop text= (that loses Preliminary). See feedback_plot_style.
        hep.cms.label(ax=ax, rlabel=RLABEL, data=is_data, text="Preliminary",
                      fontsize=11)
        plt.tight_layout()
        out = common.fig_path("comparisons", f"fom_map_run3_{tag}{LO_TAG}.png",
                              key=KEY, tool=STAGE2_TOOL)
        plt.savefig(out, dpi=120)
        plt.close(fig)
        print(f"Saved {out}")


# ---------- two keys: print each model's best cell (no combined figure -- user, 1026) ----------
if len(KEYS) == 2:
    KA, KB = KEYS
    for tag, key in (("mcsr", "fom_w"), ("data", "fom_d")):
        bA, bB = (max(RESULTS[k], key=lambda r: (r[key] if not np.isnan(r[key]) else -np.inf))
                  for k in (KA, KB))
        print(f"\n{tag}: max FOM {KA} {bA[key]:.5f} (WP1 {bA['wp1']:.2f}, WP2 {bA['wp2']:.2f}, "
              f"eps_tot {bA['eps_and']:.3f}) vs {KB} {bB[key]:.5f} (WP1 {bB['wp1']:.2f}, "
              f"WP2 {bB['wp2']:.2f}, eps_tot {bB['eps_and']:.3f}) -> max/max {bA[key] / bB[key]:.3f}")

"""Stage-1 standalone working-point scan, Punzi a=3, scaled to Run-3 (326 /fb).

Stage-1's job is to reject NON-resonant background (combinatorial + partial-match +
misID). Its working point is therefore chosen against that background only -- the
resonant piece (GenTag modes 1-5) is stage-2's problem and is NOT counted in B.

MIGRATED 2026-09-22 to the NanoMVA samples (it read Lumi_reweight /
WtoLNu2Jets_GenTag / Data, none of which exist any more, so it could not run at
all) and reworked to the 2026-09-22 request: Punzi a=3 only, the two backgrounds
scanned SEPARATELY rather than one folded FOM, everything scaled to 326 /fb.

    FOM = eps_sig / (1.5 + sqrt(B))        Punzi a=3

Punzi is BR-INDEPENDENT, which is exactly why it is the right choice here: the
signal normalisation carries whatever BR(H->DsW) the sample was generated with, so
any S/sqrt(B) number would inherit that assumption. The scaled signal yield is
printed once as context and deliberately kept out of the comparison.

TWO BACKGROUNDS, SCANNED SEPARATELY
  B_wjets   WtoLNu2Jets_merged, cand_gen_mode in {0,6,7,8,9,10}, weight-weighted,
            SR window |M(KKpi) - M_Ds| < 40 MeV. Pure non-resonant W+jets -- but
            tt/DY/QCD are absent entirely (the pre-cut data/MC gap is 1.44 and
            rises steeply with jet multiplicity), so this UNDERSTATES B.
  B_data    OS data in the CONTROL REGION, |M(KKpi) - M_Ds| >= 40 MeV, counted as
            it stands. The SR stays blind and is NOT extrapolated into. Model-free
            and it does contain tt/DY/QCD, at the cost of not being the SR.

THE TWO ARMS ARE NOT COMPARABLE IN ABSOLUTE VALUE. B_wjets is an 80 MeV SR yield,
B_data a 520 MeV sideband yield, so B_data is larger by roughly the width ratio and
its Punzi value is correspondingly smaller. What IS comparable, and the only thing
read across the two arms, is WHERE IN WP1 the optimum sits.

WHY THE 326 /fb SCALING MOVES THE ANSWER: Punzi's a/2 term is what creates an
interior maximum. At 34.67 /fb it is a ~12% regularizer at WP1=0.55; at 326 /fb
the same cell has sqrt(B) three times larger and 1.5 is ~4%, so Punzi collapses
toward eps/sqrt(B), which has NO interior maximum. Scaling multiplies weights but
adds no MC rows and no data events, so what bounds the scan is the statistics
floor below, NOT the FOM. Any maximum sitting on the floor must be reported as
"still rising where the background estimate runs out of events".

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 fom_scan_stage1.py
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

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL, tool='DNN')}/Prediction"
MC_SAMPLE = "WtoLNu2Jets_merged"
DATA_SAMPLE = "Data_OS_merged"

M_DS, SR_HALF = 1.96835, 0.040
NONRES = [0, 6, 7, 8, 9, 10]
RES = [1, 2, 3, 4, 5]
PUNZI_A = 3.0

LUMI_NOW_PB = sum(e["lumi_pb"] for e in eras.ERAS.values())     # 2022 preEE+postEE
LUMI_RUN3_PB = 326000.0
SCALE = LUMI_RUN3_PB / LUMI_NOW_PB

WP1_GRID = [0.99, 0.97, 0.95, 0.92, 0.90, 0.87, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60,
            0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]
# Statistics floors. Lumi scaling adds no events, so these are unchanged by it and
# they -- not the FOM -- are what bounds the usable scan.
#
# NEFF_B_FLOOR is a ROUND NUMBER CHOSEN FOR ~10% Punzi error (Kish N_eff = 25 -> 20%
# on B -> 10% on sqrt(B)), not derived from anything. It is NOT what determines the
# answer, and the floor-sensitivity block below prints the proof: the argmax stays at
# WP1=0.55 for every floor from 5 to 40, and only moves once cells with N_eff ~ 3 are
# admitted -- where the paired bootstrap puts them at +0.4 sigma, i.e. nowhere.
NEFF_B_FLOOR = 25.0                 # effective MC rows behind B_wjets
N_CR_FLOOR = 50                     # raw OS data events in the control region
NEFF_RES_FLOOR = 1500.0             # stage-2 training background (guard rail only)
SEL = [*common.SELECTION]


def sel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a) & (a["score_dsmu"] > -1)


def neff(w):
    """Kish effective sample size: the relative statistical error on sum(w) is
    1/sqrt(neff). Rows are candidates, not events, so this would understate the
    error if one event contributed several rows -- measured 2026-09-22 on the
    non-resonant SR subset: 490,965 rows for 490,965 unique (run, lumi, event),
    exactly 1.000 rows/event, so rows ARE independent events here and both this
    and the row-wise Poisson bootstrap below are valid."""
    s2 = np.sum(w ** 2)
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


def punzi(eps, B):
    return eps / (PUNZI_A / 2 + np.sqrt(B)) if B > 0 else np.nan


# ---------- signal: thresholds and weighted efficiency ----------
s_score, s_w = [], []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match", "weight"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s_score.append(a["score_dsmu"][m])
    s_w.append(a["weight"][m])
s_score = np.concatenate(s_score)
s_w = np.concatenate(s_w)
S_TOT = s_w.sum()
print(f"Luminosity: {LUMI_NOW_PB/1000:.2f} /fb -> {LUMI_RUN3_PB/1000:.0f} /fb  "
      f"(scale {SCALE:.4f})")
print(f"truth-matched signal: {len(s_score):,} rows, weighted {S_TOT:.4f} events "
      f"-> {S_TOT*SCALE:.3f} at Run-3.\n  CONTEXT ONLY: that yield carries the "
      f"generated BR(H->DsW); Punzi is BR-independent and is the number to compare.")

# Thresholds are defined on the FULL-range truth-matched signal, same as every other
# stage-1 WP in this analysis (common.stage1_threshold), so WP1=0.90 here is the
# same cut as the stage-2 pre-filter.
THR = {wp: float(np.quantile(s_score, 1 - wp)) for wp in WP1_GRID}

# ---------- background A: W+jets non-resonant, SR ----------
b = uproot.open(f"{PRED}/{MC_SAMPLE}.root")["Events"].arrays(
    ["score_dsmu", "cand_gen_mode", "Ds_mass", "weight"] + SEL, library="np")
k = sel(b)
b_sc, b_mode, b_mass, b_w = (b["score_dsmu"][k], b["cand_gen_mode"][k],
                             b["Ds_mass"][k], b["weight"][k])
b_nr_sr = np.isin(b_mode, NONRES) & (np.abs(b_mass - M_DS) < SR_HALF)
b_res = np.isin(b_mode, RES)
print(f"{MC_SAMPLE}: {k.sum():,} rows after selection, "
      f"{int(b_nr_sr.sum()):,} non-resonant in SR")

# ---------- background B: OS data control region (SR left blind) ----------
d = uproot.open(f"{PRED}/{DATA_SAMPLE}.root")["Events"].arrays(
    ["score_dsmu", "Ds_mass"] + SEL, library="np")
kd = sel(d)
d_sc, d_mass = d["score_dsmu"][kd], d["Ds_mass"][kd]
d_cr = np.abs(d_mass - M_DS) >= SR_HALF
print(f"{DATA_SAMPLE}: {kd.sum():,} rows after selection, "
      f"{int(d_cr.sum()):,} in the control region\n")

# ---------- scan ----------
rows = []
for wp in WP1_GRID:
    thr = THR[wp]
    w_nr = b_w[(b_sc > thr) & b_nr_sr]
    n_cr = int(np.sum((d_sc > thr) & d_cr))
    rows.append(dict(
        wp=wp, thr=thr,
        eps=float(s_w[s_score > thr].sum() / S_TOT),
        B_w=float(w_nr.sum()) * SCALE, raw_w=int(w_nr.size), neff_w=neff(w_nr),
        neff_res=neff(b_w[(b_sc > thr) & b_res]),
        n_cr=n_cr, B_d=n_cr * SCALE,
    ))

hdr = (f"{'WP1':>5} {'thr':>8} {'eps':>7} | {'B_wjets':>9} {'rawB':>6} {'NeffB':>6} "
       f"{'Punzi_w':>8} | {'n_CR':>7} {'B_data':>10} {'Punzi_d':>8} | {'Neff(1-5)':>9}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    ok_w, ok_d = r["neff_w"] >= NEFF_B_FLOOR, r["n_cr"] >= N_CR_FLOOR
    flag = "" if (ok_w and ok_d) else ("   below floor: " + ",".join(
        x for x, o in (("MC", ok_w), ("data", ok_d)) if not o))
    print(f"{r['wp']:>5.2f} {r['thr']:>8.4f} {r['eps']:>7.4f} | "
          f"{r['B_w']:>9.1f} {r['raw_w']:>6d} {r['neff_w']:>6.0f} "
          f"{punzi(r['eps'], r['B_w']):>8.5f} | "
          f"{r['n_cr']:>7d} {r['B_d']:>10.1f} {punzi(r['eps'], r['B_d']):>8.5f} | "
          f"{r['neff_res']:>9.0f}{flag}")


def best(key, name, floor_ok):
    usable = [r for r in rows if floor_ok(r)]
    vals = [punzi(r["eps"], r[key]) for r in usable]
    i = int(np.nanargmax(vals))
    r = usable[i]
    edge = ("AT THE STATISTICS FLOOR -- the FOM is still rising where the background "
            "estimate runs out of events" if i == len(usable) - 1 else "interior")
    print(f"{name:>46}: max {vals[i]:.5f} at WP1 = {r['wp']:.2f} "
          f"(thr {r['thr']:.4f}, B = {r[key]:.1f})  [{edge}]")


print(f"\n---------- best WP1, Punzi a={PUNZI_A:.0f}, {LUMI_RUN3_PB/1000:.0f} /fb "
      f"(inside the statistics floor) ----------")
best("B_w", f"B = W+jets non-resonant SR (NeffB>={NEFF_B_FLOOR:.0f})",
     lambda r: r["neff_w"] >= NEFF_B_FLOOR)
best("B_d", f"B = OS data control region (n_CR>={N_CR_FLOOR})",
     lambda r: r["n_cr"] >= N_CR_FLOOR)

wr = [r for r in rows if r["neff_res"] >= NEFF_RES_FLOOR]
if wr:
    print(f"\nGuard rail (stage-2 training, not part of this FOM): tightest WP1 with "
          f"Neff(modes 1-5) >= {NEFF_RES_FLOOR:.0f} is {min(r['wp'] for r in wr):.2f}")

# ---------- is the floor value doing the work? ----------
# The point of a floor is to stop the scan where the background estimate runs out of
# events. If the reported argmax moved around as the floor moved, the floor would BE
# the answer -- so print the dependence rather than assert independence.
print("\nFloor sensitivity, B = W+jets non-resonant: argmax of Punzi among WP1 with "
      "NeffB >= floor")
for fl in (50, 40, 30, 25, 20, 15, 10, 5, 3, 0):
    ok = [r for r in rows if r["neff_w"] >= fl]
    if not ok:
        continue
    bst = max(ok, key=lambda r: punzi(r["eps"], r["B_w"]))
    print(f"  floor {fl:>3}: usable down to WP1={min(r['wp'] for r in ok):.2f}, "
          f"argmax = {bst['wp']:.2f} (Punzi {punzi(bst['eps'], bst['B_w']):.5f})")

# ---------- paired Poisson bootstrap: is the maximum distinguishable? ----------
print(f"\n---------- paired Poisson bootstrap ({400} replicas, one resample shared "
      f"across all WP1) ----------")
rng = np.random.default_rng(42)
NBOOT = 400
BG = list(WP1_GRID)
EPS = {r["wp"]: r["eps"] for r in rows}

# Each arm resamples ONE pool of rows and reuses that single resample across every
# WP1, so the cut-to-cut differences carry the correct (correlated) error -- nested
# cuts share events, and treating the WP1 columns as independent would inflate the
# difference errors by roughly sqrt(2).
ARMS = (
    ("B = W+jets non-resonant SR", b_w[b_nr_sr],
     {wp: (b_sc > THR[wp])[b_nr_sr] for wp in BG}),
    ("B = OS data control region", np.full(int(d_cr.sum()), 1.0),
     {wp: (d_sc > THR[wp])[d_cr] for wp in BG}),
)
for label, pool_w, pool_pass in ARMS:
    fom = {wp: np.empty(NBOOT) for wp in BG}
    for i in range(NBOOT):
        c = rng.poisson(1.0, len(pool_w))
        for wp in BG:
            fom[wp][i] = punzi(EPS[wp], float((pool_w * c)[pool_pass[wp]].sum()) * SCALE)
    means = {wp: float(np.nanmean(fom[wp])) for wp in BG}
    wstar = max(BG, key=lambda wp: means[wp])
    print(f"\n{label}  (bootstrap-mean peak at WP1 = {wstar:.2f})")
    print(f"{'WP1':>5} {'FOM':>9} {'+-':>8} {'rel':>6} | {'vs peak':>9} "
          f"{'+-(paired)':>11} {'sigma':>6}")
    for wp in BG:
        dd = fom[wp] - fom[wstar]
        sd = float(np.nanstd(dd))
        sig = float(np.nanstd(fom[wp]))
        print(f"{wp:>5.2f} {means[wp]:>9.5f} {sig:>8.5f} {sig/means[wp]:>6.1%} | "
              f"{float(np.nanmean(dd)):>9.5f} {sd:>11.5f} "
              f"{(float(np.nanmean(dd))/sd if sd > 0 else 0):>6.2f}")

# ---------- figure ----------
eps = np.array([r["eps"] for r in rows])
RLABEL = rf"{LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$ (13.6 TeV)"
fig, ax = plt.subplots(1, 2, figsize=(16, 7))
for i, (key, c, lab) in enumerate((
        ("B_w", "seagreen", r"B = W+jets non-resonant (SR, 80 MeV)"),
        ("B_d", "steelblue", r"B = OS data (control region, SR blind)"))):
    y = np.array([punzi(r["eps"], r[key]) for r in rows])
    ax[i].plot(eps, y, "-", marker="o", color=c, label=lab)
    ax[i].set_xlabel(r"stage-1 $\varepsilon_{sig}$ (WP1)")
    ax[i].set_ylabel(r"Punzi FOM $=\varepsilon_{sig}/(1.5+\sqrt{B})$")
    ax[i].legend(fontsize=11, loc="upper left")
    ax[i].invert_xaxis()
    # text=, NOT label= -- mplhep deprecated the `label` kwarg of cms.label(), so it
    # emits a FutureWarning and SILENTLY DROPS the word. Every other macro in this
    # directory uses text=; match them.
    hep.cms.label(ax=ax[i], rlabel=RLABEL, data=(key == "B_d"),
                  text="Preliminary", fontsize=11)
fig.suptitle(rf"stage-1 alone, Punzi a=3, scaled to {LUMI_RUN3_PB/1000:.0f} fb$^{{-1}}$"
             "   (the two panels share no normalisation -- compare positions, not values)",
             fontsize=14)
plt.tight_layout()
out = common.fig_path("comparisons", "fom_scan_stage1.png", key=MODEL)
plt.savefig(out, dpi=120)
plt.close(fig)
print(f"\nSaved {out}")

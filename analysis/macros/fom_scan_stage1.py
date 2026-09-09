"""Stage-1 working-point scan for the two-stage MVA (see the stage-2 plan, section 4a).

Stage-1's job is to reject NON-resonant background (combinatorial + partial-match +
misID). Its working point is therefore chosen against that background only -- the
resonant piece (GenTag modes 1-5) is stage-2's problem and is NOT counted in B here.

For a grid of stage-1 signal efficiencies WP1:
  eps_sig(WP1)  Lumi_reweight-weighted survival of truth-matched HtoDsW+HtoDsstW
                (OS + tight b-veto). S has no absolute scale (signal_era weights are
                relative), so S proportional to eps_sig and the FOM numerator is eps_sig.
  B_nonres      lumi-weighted SR yield (|M(KKpi) - M_Ds| < 40 MeV, both eras) of
                WtoLNu2Jets_GenTag with cand_gen_mode in {0,6,7,8,9,10}.
  B_data_sb     OS-sideband data count surviving the cut -- a second, model-free
                probe of "how fast does non-resonant background fall" (cross-check
                only; different normalisation, compare where the curves flatten).

  FOM (Punzi)   eps_sig / (a/2 + sqrt(B)),  tabulated for a in {2, 3, 5}
  FOM (naive)   eps_sig / sqrt(B)

Guard-rail columns, because the Punzi a/2 term is negligible here (B is in the
hundreds+, sqrt(B) >> a/2) so the FOM has no interior maximum -- what actually
bounds WP1 is:
  N_eff(1-5, full)   effective size of the stage-2 TRAINING background after this
                     cut = (sum w)^2 / sum w^2 for modes 1-5, full mass range.
                     Pick the TIGHTEST WP1 that keeps this above ~1500.
  N_eff(B_nonres)     if this is small the FOM value itself is noise.

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

plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL)}/Prediction"
M_DS, SR_HALF = 1.96835, 0.040
NONRES = [0, 6, 7, 8, 9, 10]
RES = [1, 2, 3, 4, 5]
WP1_GRID = [0.99, 0.97, 0.95, 0.92, 0.90, 0.87, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50]
NEFF_FLOOR = 1500.0
SEL = ["num_os_mu_tight", "num_ss_mu_tight", common.B_VETO, "score_dsmu"]


def sel(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
            & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1))


def neff(w):
    s2 = np.sum(w ** 2)
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


# ---------- signal: thresholds and weighted efficiency ----------
s_score, s_w = [], []
for s in common.SIGNALS:
    a = uproot.open(f"{PRED}/{s}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match", "Lumi_reweight"] + SEL, library="np")
    m = sel(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    s_score.append(a["score_dsmu"][m])
    s_w.append(a["Lumi_reweight"][m])
s_score = np.concatenate(s_score)
s_w = np.concatenate(s_w)
S_TOT = s_w.sum()
print(f"truth-matched signal: {len(s_score):,} rows, weighted {S_TOT:,.1f}")

THR = {wp: float(np.quantile(s_score, 1 - wp)) for wp in WP1_GRID}

# ---------- background: GenTag ----------
b = uproot.open(f"{PRED}/WtoLNu2Jets_GenTag.root")["Events"].arrays(
    ["score_dsmu", "cand_gen_mode", "Ds_mass", "Lumi_reweight"] + SEL, library="np")
k = sel(b)
b_sc, b_mode, b_mass, b_w = (b["score_dsmu"][k], b["cand_gen_mode"][k],
                             b["Ds_mass"][k], b["Lumi_reweight"][k])
b_sr = np.abs(b_mass - M_DS) < SR_HALF
b_nonres = np.isin(b_mode, NONRES)
b_res = np.isin(b_mode, RES)
print(f"GenTag scored rows after selection: {k.sum():,}")

# ---------- cross-check background: OS sideband data ----------
d = uproot.open(f"{PRED}/Data.root")["Events"].arrays(
    ["score_dsmu", "Ds_mass"] + SEL, library="np")
kd = sel(d)
d_sc, d_mass = d["score_dsmu"][kd], d["Ds_mass"][kd]
d_sb = np.abs(d_mass - M_DS) >= SR_HALF     # OS data is SR-blind anyway
print(f"OS data scored rows after selection: {kd.sum():,}\n")

# ---------- scan ----------
rows = []
for wp in WP1_GRID:
    thr = THR[wp]
    pas = b_sc > thr
    w_sr_nr = b_w[pas & b_sr & b_nonres]
    B = float(w_sr_nr.sum())
    raw_nr = int((pas & b_sr & b_nonres).sum())
    w_res_full = b_w[pas & b_res]
    neff_res = neff(w_res_full)
    Bd = int(np.sum((d_sc > thr) & d_sb))
    eps = float(s_w[s_score > thr].sum() / S_TOT)
    rows.append(dict(wp=wp, thr=thr, eps=eps, B=B, raw=raw_nr, neff_B=neff(w_sr_nr),
                     neff_res=neff_res, Bd=Bd))

hdr = (f"{'WP1':>5} {'thr':>8} {'eps_sig':>8} {'B_nonres':>10} {'rawB':>7} "
       f"{'NeffB':>7} {'Neff(1-5)':>10} {'Bdata_sb':>9} "
       f"{'FOM a2':>8} {'FOM a3':>8} {'FOM a5':>8} {'eps/vB':>8}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    vB = np.sqrt(r["B"]) if r["B"] > 0 else np.nan
    fa = {a: r["eps"] / (a / 2 + vB) for a in (2, 3, 5)}
    flag = "  <-- Neff(1-5) below floor" if r["neff_res"] < NEFF_FLOOR else ""
    print(f"{r['wp']:>5.2f} {r['thr']:>8.4f} {r['eps']:>8.4f} {r['B']:>10.1f} "
          f"{r['raw']:>7d} {r['neff_B']:>7.0f} {r['neff_res']:>10.0f} {r['Bd']:>9d} "
          f"{fa[2]:>8.5f} {fa[3]:>8.5f} {fa[5]:>8.5f} "
          f"{(r['eps']/vB if vB>0 else np.nan):>8.5f}{flag}")

ok = [r for r in rows if r["neff_res"] >= NEFF_FLOOR]
if ok:
    wp_star = min(ok, key=lambda r: r["wp"])["wp"]
    print(f"\nTightest WP1 with Neff(modes 1-5, full) >= {NEFF_FLOOR:.0f}: "
          f"eps_sig = {wp_star:.2f}  (thr = {THR[wp_star]:.4f})")
    print("This is the WP1* candidate -- fill it into TRAININGS['wjetsres'] "
          "stage1_cut.eff_sig and out_name '_wp<NN>'.")
else:
    print(f"\nNo WP1 keeps Neff(modes 1-5) above {NEFF_FLOOR:.0f} -- loosen the floor "
          "or reconsider stage-2 (BDT on fewer features, or a cut set).")

# ---------- figure ----------
wps = np.array([r["wp"] for r in rows])
eps = np.array([r["eps"] for r in rows])
B = np.array([r["B"] for r in rows])
Bd = np.array([float(r["Bd"]) for r in rows])
neff_res = np.array([r["neff_res"] for r in rows])
vB = np.sqrt(np.where(B > 0, B, np.nan))
vBd = np.sqrt(np.where(Bd > 0, Bd, np.nan))

fig, ax = plt.subplots(1, 2, figsize=(17, 7))
for a_, ls in ((2, "-"), (5, "--")):
    ax[0].plot(eps, eps / (a_ / 2 + vB), ls, marker="o", color="seagreen",
               label=rf"B$_{{\rm nonres}}$ (Punzi a={a_})")
    ax[0].plot(eps, eps / (a_ / 2 + vBd), ls, marker="s", color="steelblue",
               label=rf"B$_{{\rm data\ sb}}$ (Punzi a={a_})")
ax[0].set_xlabel(r"stage-1 $\varepsilon_{sig}$")
ax[0].set_ylabel(r"FOM $= \varepsilon_{sig} / (a/2 + \sqrt{B})$")
ax[0].legend(fontsize=11)
ax[0].invert_xaxis()
hep.cms.label(ax=ax[0], data=False, text="Simulation", fontsize=11)

ax[1].plot(eps, neff_res, "o-", color="#D1495B")
ax[1].axhline(NEFF_FLOOR, color="black", ls=":", label=f"floor {NEFF_FLOOR:.0f}")
ax[1].set_xlabel(r"stage-1 $\varepsilon_{sig}$")
ax[1].set_ylabel(r"N$_{\rm eff}$(modes 1-5, full range) = stage-2 training background")
ax[1].legend(fontsize=11)
ax[1].invert_xaxis()

plt.tight_layout()
out = common.fig_path("comparisons", "fom_scan_stage1.png", key=MODEL)
plt.savefig(out, dpi=120)
plt.close(fig)
print(f"\nSaved {out}")

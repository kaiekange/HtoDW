"""Corrected Higgs mass m(D_s mu nu) at a fixed 2D working point: signal vs W+jets.

The "correction" is the neutrino: its px/py come from the corrected Puppi MET and its
pz from the W-mass constraint, so the Ds+mu system is completed into a full H
candidate instead of the visible-only m(Ds mu). Exactly the construction used by
check_dsmunu_mass.py / compare_dsmunu_signal_gentag.py, whose functions are imported
rather than re-derived (complex discriminant -> real part, two real roots -> smaller
|pz|). m(Ds mu)_corr is drawn alongside as the baseline it is meant to improve on.

W+jets is the absolute lumi-weighted yield, stacked by the same four GenTag truth
groups as draw_wp_mass.py. Signal is a LINE scaled by a single printed factor -- its
absolute yield carries whatever BR(H->DsW) the sample was generated with, so only its
SHAPE is meaningful here.

No mass window and no blinding: this is a signal-vs-background shape comparison, not
a data closure test, and OS data is not drawn at all.

WP2 is unconditional (same convention as fom_map_stage2_run3.py / draw_wp_mass.py);
WP2 >= 1 means stage-1 only.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_higgs_mass_wp.py [WP1] [WP2]      # defaults 0.40 0.55
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
from check_dsmunu_mass import M_MU, four, inv_mass, nu_pz  # noqa: E402
import common  # noqa: E402
import eras  # noqa: E402

plt.style.use(hep.style.CMS)

WP1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.40
WP2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.55

PRED = f"{common.out_dir_for('wjetsres', tool='BDT')}/Prediction"
LUMI = sum(e["lumi_pb"] for e in eras.ERAS.values())
SEL = [*common.SELECTION]
KIN = ["mu_px", "mu_py", "mu_pz", "mu_pt", "mu_pt_corr", "Dsmu_px", "Dsmu_py",
       "Dsmu_pz", "Dsmu_mass", "Dsmu_mass_corr", "met_pt_corr", "met_phi_corr"]

GROUPS = [   # stack bottom -> top, identical to draw_wp_mass.py
    (r"other (mode 0,6-10)",         [0, 6, 7, 8, 9, 10], "#5a9e78"),
    (r"D$\to$K$\pi\pi$ misID (5)",   [5],                 "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4],              "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2],              "#c0392b"),
]
VARS = [  # (key, xlabel, bins)
    ("m", r"corrected $M(D_s^\pm\mu^\mp\nu)$ [GeV]", np.arange(40.0, 260.0 + 1e-9, 5.0)),
    ("base", r"$M(D_s^\pm\mu^\mp)$ [GeV]  (baseline)", np.arange(0.0, 130.0 + 1e-9, 2.5)),
]


def ok(a):
    return (common.selection_mask(a) & (a["score_dsmu"] > -1)
            & (a["score_dsmu_res"] > -1))


def load(name, extra):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        ["score_dsmu", "score_dsmu_res", "weight"] + SEL + KIN + extra, library="np")
    f = {k: a[k].astype(np.float64) for k in KIN}
    # The muon is rescaled to its corrected pt, but Dsmu is the UNCORRECTED four-vector,
    # so the Ds is recovered by subtracting the UNCORRECTED muon -- mixing the two would
    # leave a spurious pt shift inside the Ds.
    mu_raw = four(f["mu_px"], f["mu_py"], f["mu_pz"], M_MU)
    k = f["mu_pt_corr"] / f["mu_pt"]
    mu = four(k * f["mu_px"], k * f["mu_py"], k * f["mu_pz"], M_MU)
    ds = four(f["Dsmu_px"], f["Dsmu_py"], f["Dsmu_pz"], f["Dsmu_mass"]) - mu_raw
    nux = f["met_pt_corr"] * np.cos(f["met_phi_corr"])
    nuy = f["met_pt_corr"] * np.sin(f["met_phi_corr"])
    pz1, pz2, cplx = nu_pz(mu, nux, nuy)
    pz = np.where(np.abs(pz1) < np.abs(pz2), pz1, pz2)
    out = dict(m=inv_mass(ds + mu + four(nux, nuy, pz, 0.0)),
               base=f["Dsmu_mass_corr"], w=a["weight"].astype(np.float64),
               cplx=cplx, keep=ok(a), s1=a["score_dsmu"], s2=a["score_dsmu_res"])
    out.update({c: a[c] for c in extra})
    return out


sig = load("Signal_merged", ["Ds_match", "mu_match"])
ms = sig["keep"] & sig["Ds_match"].astype(bool) & sig["mu_match"].astype(bool)
thr1 = float(np.quantile(sig["s1"][ms], 1 - WP1))
STAGE1_ONLY = WP2 >= 1.0
thr2 = -np.inf if STAGE1_ONLY else float(np.quantile(sig["s2"][ms], 1 - WP2))
eps_tot = float(sig["w"][ms & (sig["s1"] > thr1) & (sig["s2"] > thr2)].sum()
                / sig["w"][ms].sum())
print(f"WP1={WP1:.3f} (thr1={thr1:.4f})  "
      + ("stage-1 ONLY" if STAGE1_ONLY else f"WP2={WP2:.2f} (thr2={thr2:.4f})")
      + f"  eps_tot={eps_tot:.4f}")

sk = ms & (sig["s1"] > thr1) & (sig["s2"] > thr2)
bkg = load("WtoLNu2Jets_merged", ["cand_gen_mode"])
bk = bkg["keep"] & (bkg["s1"] > thr1) & (bkg["s2"] > thr2)
print(f"signal rows {int(sk.sum()):,} (weighted {sig['w'][sk].sum():.4f}),  "
      f"W+jets rows {int(bk.sum()):,} (weighted {bkg['w'][bk].sum():.1f})")
for nm, d, kk in (("signal", sig, sk), ("W+jets", bkg, bk)):
    print(f"  {nm}: complex-pz (no real W-constraint solution) "
          f"{100 * np.average(d['cplx'][kk], weights=np.abs(d['w'][kk])):.1f}% of weight")

fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.0))
for ax, (key, xlabel, bins) in zip(axes, VARS):
    bc = 0.5 * (bins[1:] + bins[:-1])
    width = bins[1] - bins[0]
    sv, swv = sig[key][sk], sig["w"][sk]
    bv, bwv, bm = bkg[key][bk], bkg["w"][bk], bkg["cand_gen_mode"][bk]

    hs, errs = [], []
    for _, modes, _ in GROUPS:
        mm = np.isin(bm, modes)
        h, _ = np.histogram(bv[mm], bins=bins, weights=bwv[mm])
        e2, _ = np.histogram(bv[mm], bins=bins, weights=bwv[mm] ** 2)
        hs.append(h); errs.append(np.sqrt(e2))
    tot = np.sum(hs, axis=0)
    tot_err = np.sqrt(np.sum(np.square(errs), axis=0))
    bottom = np.zeros(len(bc))
    for (lab, _, cl), h in zip(GROUPS, hs):
        ax.bar(bc, h, width=width, bottom=bottom, color=cl, alpha=0.9,
               label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
        bottom = bottom + h
    ax.bar(bc, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
           label="W+jets MC stat")

    hsig, _ = np.histogram(sv, bins=bins, weights=swv)
    # One scale factor for the whole panel: the signal's absolute normalisation carries
    # the generated BR, so only its shape is being compared here.
    sf = (0.55 * tot.max() / hsig.max()) if hsig.max() > 0 else 1.0
    ax.step(bins, np.r_[hsig[0], hsig] * sf, where="pre", color="black", lw=2.0,
            label=rf"signal $\times${sf:,.0f}  ($\varepsilon_{{tot}}$={eps_tot:.1%})")

    lo, hi = (np.percentile(sv, [16, 84]) if key == "m" else (np.nan, np.nan))
    if key == "m":
        ax.axvspan(lo, hi, color="black", alpha=0.07, zorder=0)
        frac = float(bwv[(bv >= lo) & (bv <= hi)].sum() / bwv.sum())
        print(f"  signal central 68% of corrected M(Dsmunu): [{lo:.1f}, {hi:.1f}] GeV; "
              f"W+jets fraction inside = {frac:.3f}")
        # in-axes text, not set_title(): the title row is occupied by hep.cms.label()
        ax.text(0.03, 0.95, rf"signal 68% = [{lo:.0f}, {hi:.0f}] GeV" "\n"
                rf"W+jets inside = {frac:.1%}", transform=ax.transAxes,
                ha="left", va="top", fontsize=11)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, max(bottom.max(), (hsig * sf).max()) * 1.45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"Events / {width:g} GeV")
    ax.legend(fontsize=8.5, loc="upper right")
    hep.cms.label(ax=ax, data=False, text="Preliminary",
                  rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

fig.suptitle((rf"stage-1 only, WP1 $\varepsilon_1$={WP1:.1%}" if STAGE1_ONLY else
              rf"WP1 $\varepsilon_1$={WP1:.0%}  +  WP2 $\varepsilon_2$={WP2:.0%}")
             + rf"  $\to$  $\varepsilon_{{tot}}$={eps_tot:.1%}", fontsize=14)
plt.tight_layout()
name = (f"higgs_mass_wp1_{WP1:.2f}_stage1only.png" if STAGE1_ONLY
        else f"higgs_mass_wp1_{WP1:.2f}_wp2_{WP2:.2f}.png")
out = common.fig_path("comparisons", name, key="wjetsres", tool="BDT")
plt.savefig(out, dpi=120)
plt.close(fig)
print(f"\nSaved {out}")

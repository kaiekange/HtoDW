"""ttbar CR: TTto2L2Nu vs TTtoLNu2Q, in ==1l and ==2l, Ds_mass / phi_mass by truth.

Input: the ttcr skim (pass_mu & ~jvm_veto & n_bjet >= 1 & phi window), OS and SS ttbar MC only.
No n_b split beyond the skim's n_bjet >= 1.
    1l   n_xmu == 0 & n_xele == 0                         (same as 1l_b1 in draw_cr_ttbar.py)
    2l   n_xmu + n_xele == 1, extra lepton opposite-sign to the analysis muon,
         |M(mumu) - 91| > 15 for the mumu case (as draw_cr_ttbar.py's mumu)
MC weight: skim weight x common.br_reweight x w_xmu (mumu rows only), as draw_cr_ttbar.py.

Per charge and category one figure: Ds_mass and phi_mass stacked by truth for 2L2Nu (col 1) and
LNu2Q (col 2), and the two shapes normalised to unit area (col 3). Yield table printed.

    python3 draw_cr_ttbar_decay.py
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

CR_DIR = f"{common.SAMPLES_DIR}/CR/ttcr"
MODEL_KEY = "ssdata_f45"   # only picks the figure directory
M_DS, M_D, SR_HALF, M_Z = 1.96835, 1.86966, 0.040, 91.19
FILES = {"OS": "TT_merged", "SS": "TT_SS_merged"}
DECAYS = [("2L2Nu", 0), ("LNu2Q", 1)]
COLS = ["Ds_mass", "phi_mass", "weight", "mu_charge", "n_xmu", "n_xele", "xmu_charge",
        "xele_charge", "m_mumu", "w_xmu", "tt_decay", "cand_gen_mode",
        "cand_gen_common_hf_pdgId", "cand_hf_from_b"]
TRUTH = [  # same split and colours as draw_cr_ttbar.py, stack bottom -> top
    ("combinatorial (0,6-10)", lambda a: np.isin(a["cand_gen_mode"], [0, 6, 7, 8, 9, 10]), "#9e9e9e"),
    (r"D$\to$K$\pi\pi$ misID (5)", lambda a: a["cand_gen_mode"] == 5, "#e6c229"),
    ("true D (3+4) prompt", lambda a: np.isin(a["cand_gen_mode"], [3, 4]) & (a["cand_hf_from_b"] != 1), "#f4a582"),
    ("true D (3+4) from b", lambda a: np.isin(a["cand_gen_mode"], [3, 4]) & (a["cand_hf_from_b"] == 1), "#e07a3c"),
    (r"true D$_s$ (1+2) prompt", lambda a: np.isin(a["cand_gen_mode"], [1, 2]) & (a["cand_hf_from_b"] != 1), "#f08080"),
    (r"true D$_s$ (1+2) from b", lambda a: np.isin(a["cand_gen_mode"], [1, 2]) & (a["cand_hf_from_b"] == 1), "#c0392b"),
]
CATS = {
    "1l": lambda a: (a["n_xmu"] == 0) & (a["n_xele"] == 0),
    "2l": lambda a: ((a["n_xmu"] + a["n_xele"]) == 1)
                    & (np.where(a["n_xmu"] == 1, a["xmu_charge"], a["xele_charge"]) == -a["mu_charge"])
                    & ((a["n_xmu"] == 0) | (np.abs(a["m_mumu"] - M_Z) > 15)),
}
DS_BINS = np.arange(M_DS - 0.30, M_DS + 0.30 + 1e-9, 0.020)
PHI_BINS = np.linspace(common.M_PHI - 0.020, common.M_PHI + 0.020, 11)
VARS = [("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]", DS_BINS), ("phi_mass", r"$M(K^+K^-)$ [GeV]", PHI_BINS)]


def read(fname):
    with uproot.open(f"{CR_DIR}/{fname}.root")["Events"] as t:
        a = t.arrays(COLS, library="np")
    w = a["weight"].astype(float) * common.br_reweight(a)
    a["w"] = np.where(a["n_xmu"] == 1, w * a["w_xmu"], w)
    return a


def hist(x, w, bins):
    x = np.clip(x, bins[0], bins[-1] - 1e-9)
    return np.histogram(x, bins, weights=w)[0], np.histogram(x, bins, weights=w * w)[0]


for charge, fname in FILES.items():
    A = read(fname)
    for cat, cf in CATS.items():
        m = cf(A)
        sr = np.abs(A["Ds_mass"] - M_DS) < SR_HALF
        dwin = np.abs(A["Ds_mass"] - M_D) < 0.030
        ytot = A["w"][m].sum()
        print(f"\n===== {charge} {cat}: ttbar MC {ytot:,.1f}  (raw rows {m.sum():,})")
        for dname, code in DECAYS:
            k = m & (A["tt_decay"] == code)
            y = A["w"][k].sum()
            print(f"  {dname}: {y:11,.1f} ({100 * y / ytot:5.1f}% of ttbar)  raw rows {k.sum():,}")
            for lab, tf, _c in TRUTH:
                kk = k & tf(A)
                plain = lab.replace("$", "").replace(r"\to", "->").replace(r"\pi", "pi")
                print(f"      {plain:24s} all {A['w'][kk].sum():10,.1f} ({100 * A['w'][kk].sum() / y:5.1f}%)"
                      f"   Ds SR {A['w'][kk & sr].sum():8,.1f}   D+- +-30 {A['w'][kk & dwin].sum():8,.1f}")

        fig, axes = plt.subplots(2, 3, figsize=(21, 13))
        for row, (col, xlab, bins) in enumerate(VARS):
            bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
            shapes = {}
            for ic, (dname, code) in enumerate(DECAYS):
                ax = axes[row, ic]
                k = m & (A["tt_decay"] == code)
                bottom, t2 = np.zeros(len(bc)), np.zeros(len(bc))
                for lab, tf, colr in TRUTH:
                    kk = k & tf(A)
                    hh, e2 = hist(A[col][kk], A["w"][kk], bins)
                    ax.bar(bc, hh, width=wd, bottom=bottom, color=colr, edgecolor="none",
                           label=f"{lab}  N={hh.sum():,.0f}")
                    bottom += hh; t2 += e2
                err = np.sqrt(t2)
                ax.bar(bc, 2 * err, width=wd, bottom=bottom - err, color="none", hatch="////",
                       edgecolor="0.35", linewidth=0.0, label="MC stat")
                shapes[dname] = (bottom, err)
                ax.set_xlim(bins[0], bins[-1]); ax.set_ylim(0, 1.45 * bottom.max())
                ax.set_xlabel(xlab); ax.set_ylabel("Events / bin")
                ax.text(0.04, 0.88, rf"$t\bar{{t}}$ {dname}", transform=ax.transAxes, fontsize=18)
                ax.legend(fontsize=9, loc="upper right")
                if row == 0 and ic == 0:
                    hep.cms.label(ax=ax, data=False, text="Preliminary", com=13.6, fontsize=13,
                                  lumi=f"{(common.ERAS[1]['lumi_pb'] + common.ERAS[2]['lumi_pb']) / 1000:.1f}")
            ax = axes[row, 2]
            for (dname, _code), colr in zip(DECAYS, ("#7f5fbf", "#b39ddb")):
                h, e = shapes[dname]
                s = h.sum()
                hep.histplot(h / s, bins, yerr=e / s, ax=ax, color=colr, lw=2, label=f"{dname} (unit area)")
            ax.set_xlim(bins[0], bins[-1]); ax.set_ylim(0, None)
            ax.set_xlabel(xlab); ax.set_ylabel("Fraction / bin")
            ax.set_title("shape comparison", fontsize=16)
            ax.legend(fontsize=11)
            if col == "Ds_mass":
                for a_ in axes[row]:
                    for x0 in (M_DS, M_D):
                        a_.axvline(x0, color="k", ls=":", lw=1)
        fig.suptitle(f"ttbar CR {charge}, =={cat}: pass_mu & n_b>=1 & phi window "
                     f"(MC: skim weight x BR rw" + (" x 2nd-mu SF" if cat == "2l" else "") + ")",
                     fontsize=17)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        out = common.fig_path("sf", f"cr_ttbar_decay_{charge}_{cat}.png", key=MODEL_KEY)
        plt.savefig(out, dpi=90, bbox_inches="tight")
        plt.close(fig)
        print(f"   Saved {out}")

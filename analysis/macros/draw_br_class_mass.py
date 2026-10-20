"""Where the Pythia->PDG BR correction (common.br_reweight) lands, per sample.

For W+jets, ttbar and DY (OS, analysis selection + the model's phi window, Ds_mass full
range): one figure per sample, one panel per true-charm GenTag mode 1-4, Ds_mass split by
the common heavy-flavour ancestor (Ds 431 / D+ 411 / D0 421), dashed = Pythia BR, solid =
after br_reweight. A yield table per (mode, hf) is printed at no cut and at the eps_sig
working points, with the after/before ratio, which must equal the BR_FACTORS entry (1 for
classes not in the table) -- that is the closure check of the reweighting itself.

    python3 draw_br_class_mass.py [--key=ssdata_phi20]
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

MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata_phi20")
PHI_WIN = common.resolve(MODEL_KEY).get("phi_win")
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS = 1.96835
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
SAMPLES = [("WtoLNu2Jets_merged", "W+jets"), ("TT_merged", r"$t\bar{t}$"), ("DY_merged", "DY")]
MODES = [(1, r"mode 1: $\phi\pi$"), (2, r"mode 2: KK$\pi$ non-res."),
         (3, r"mode 3: $\phi\pi$"), (4, r"mode 4: KK$\pi$ non-res.")]
HF = [(431, r"$D_s$", "#c0392b"), (411, r"$D^+$", "#e07a3c"), (421, r"$D^0$", "#4a90c2")]
EFFS = (None, 0.90, 0.70, 0.50)
BINS = np.arange(M_DS - 0.30, M_DS + 0.30 + 1e-9, 0.010)

THR = {None: -1.0}
for e in EFFS[1:]:
    THR[e], _ = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=e))

for fname, label in SAMPLES:
    a = uproot.open(f"{PRED}/{fname}.root")["Events"].arrays(
        list(dict.fromkeys(["Ds_mass", "phi_mass", "score_dsmu", "weight", "cand_gen_mode",
                            "cand_gen_common_hf_pdgId", *common.SELECTION])), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if PHI_WIN is not None:
        m &= np.abs(a["phi_mass"] - common.M_PHI) < PHI_WIN
    mass, score = a["Ds_mass"][m], a["score_dsmu"][m]
    mode = a["cand_gen_mode"][m]
    hf = np.abs(a["cand_gen_common_hf_pdgId"][m])
    w0 = a["weight"][m].astype(float)
    w1 = w0 * common.br_reweight({"cand_gen_mode": mode, "cand_gen_common_hf_pdgId": hf})

    print(f"\n== {label} ({fname}): {m.sum():,} rows after selection"
          + (f" + |phi_mass-M_phi|<{PHI_WIN}" if PHI_WIN else ""))
    tot0 = {e: w0[score > THR[e]].sum() for e in EFFS}
    for md, _ in MODES:
        for h, hl, _c in HF:
            k = (mode == md) & (hf == h)
            if not k.any():
                continue
            fac = common.BR_FACTORS.get((md, h), 1.0)
            cells = []
            for e in EFFS:
                ke = k & (score > THR[e])
                b, af = w0[ke].sum(), w1[ke].sum()
                assert abs(af - fac * b) <= 1e-9 * max(1.0, abs(b)), (fname, md, h, e, af, b)
                cells.append(f"{b:10.1f} ({100*b/tot0[e]:5.2f}%)")
            print(f"  mode {md} hf {h}: factor {fac:.4f} | yield (share of all MC rows) "
                  f"no cut / e90 / e70 / e50: " + " | ".join(cells))
    rest = ~np.isin(mode, [1, 2, 3, 4])
    print(f"  modes outside 1-4: after/before = {w1[rest].sum()/w0[rest].sum():.6f} (must be 1)")

    fig, axes = plt.subplots(1, 4, figsize=(26, 6.5))
    for ax, (md, mlab) in zip(axes, MODES):
        for h, hl, col in HF:
            k = (mode == md) & (hf == h)
            if not k.any():
                continue
            h0, _ = np.histogram(mass[k], BINS, weights=w0[k])
            h1, _ = np.histogram(mass[k], BINS, weights=w1[k])
            fac = common.BR_FACTORS.get((md, h), 1.0)
            hep.histplot(h0, BINS, ax=ax, color=col, ls="--", lw=1.5, label=f"{hl} Pythia BR")
            hep.histplot(h1, BINS, ax=ax, color=col, lw=2, label=f"{hl} x{fac:.3f}")
        ax.axvline(M_DS, color="grey", ls=":", lw=1)
        ax.axvline(1.86966, color="grey", ls=":", lw=1)
        ax.set_title(mlab, fontsize=15, pad=24)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        ax.legend(fontsize=11)
    axes[0].set_ylabel("Events / 10 MeV")
    hep.cms.label(ax=axes[0], data=False, text="Preliminary", lumi=f"{LUMI_FB:.1f}",
                  com=13.6, fontsize=14)
    fig.suptitle(f"{label}: true-charm modes, Pythia BR vs PDG-reweighted (OS, no DNN cut)",
                 fontsize=15)
    out = common.fig_path("sf", f"br_class_mass_{fname.replace('_merged', '')}.png",
                          key=MODEL_KEY)
    plt.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")

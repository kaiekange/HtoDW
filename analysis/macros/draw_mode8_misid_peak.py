"""GenTag mode 8 单独看: KpToPi 互换质量里有没有 D 峰。

**为什么问这个**: mode 8 在 HDsPhiPiMuGenTagProducer.cc:925 是 `!one_origin`，读起来像
"真重味 + 无关径迹" 的部分重味组合。但实测它的 `cand_gen_nmatched` **只取 1 和 2、
从不取 3**，而且 `n_hf` 恒等于 `nmatched` —— 也就是说 mode 8 里**每一条配上 gen 的
径迹都有重味祖先**，那些"没有重味祖先"的径迹是**根本没配上 gen 粒子**的。
所以 mode 8 不是一个物理类别，是 **partial truth-matching 的失败桶**：真的 D/Ds 衰变
只要有 1-2 条女儿没配上，就会掉进这里而不是 mode 1-5。

如果是这样，mode 8 里应当能看到和 mode 5 同源的 D 峰，那部分就是**从 misID 反射里
漏出去的**，直接影响"MC 把反射高估多少"这个因子的分母。

按 nmatched 拆: nmatched=2 (只丢 1 条径迹) 最可能是漏掉的真单强子衰变;
nmatched=1 (丢 2 条) 就难说了。这两条曲线形状不同本身就是判据。

    conda activate ds_dnn && python3 draw_mode8_misid_peak.py
"""
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common  # noqa: E402

MODEL_KEY = "ssdata"
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
M_DS, SR_HW, M_D = 1.96835, 0.040, 1.8696
M_K, M_PI = 0.493677, 0.13957039
BINS = np.arange(1.30, 2.301, 0.010)
CENT = 0.5 * (BINS[1:] + BINS[:-1])
WIDTH = BINS[1] - BINS[0]
SEL = [*common.SELECTION]
MOM = [f"{p}_p{x}" for p in ("Kp", "Km", "pi") for x in ("x", "y", "z")]

a = uproot.open(f"{PRED}/WtoLNu2Jets_merged.root")["Events"].arrays(
    SEL + MOM + ["Ds_mass", "score_dsmu", "weight", "cand_gen_mode", "cand_gen_nmatched"],
    library="np")
m = common.selection_mask(a) & (a["score_dsmu"] > -1) & (np.abs(a["Ds_mass"] - M_DS) >= SR_HW)
a = {k: v[m] for k, v in a.items()}
e = sum(np.sqrt(a[f"{p}_px"]**2 + a[f"{p}_py"]**2 + a[f"{p}_pz"]**2 + mm**2)
        for p, mm in (("Kp", M_PI), ("Km", M_K), ("pi", M_PI)))
px, py, pz = (sum(a[f"{p}_p{x}"] for p in ("Kp", "Km", "pi")) for x in ("x", "y", "z"))
a["mis"] = np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))

COLS = [("no cut", None), (r"$\varepsilon_{sig}$=90%", 0.90),
        (r"$\varepsilon_{sig}$=70%", 0.70), (r"$\varepsilon_{sig}$=50%", 0.50)]
SERIES = [  # (label, mask fn, colour)
    ("mode 8, nmatched=2 (1 trk unmatched)", lambda x: (x["cand_gen_mode"] == 8) & (x["cand_gen_nmatched"] == 2), "#5b6fc4"),
    ("mode 8, nmatched=1 (2 trk unmatched)", lambda x: (x["cand_gen_mode"] == 8) & (x["cand_gen_nmatched"] == 1), "#9aa7e0"),
    ("mode 5 (D$\\to$K$\\pi\\pi$ misID, reference)", lambda x: x["cand_gen_mode"] == 5, "#e6c229"),
]

fig, axes = plt.subplots(1, len(COLS), figsize=(6.2 * len(COLS), 4.6))
pk = (CENT > M_D - 0.02) & (CENT < M_D + 0.02)
sb = ((CENT > M_D - 0.10) & (CENT < M_D - 0.04)) | ((CENT > M_D + 0.04) & (CENT < M_D + 0.10))
for j, (lbl, eff) in enumerate(COLS):
    thr = None if eff is None else common.stage1_threshold(
        dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))[0]
    s = np.ones(len(a["mis"]), bool) if thr is None else a["score_dsmu"] > thr
    ax = axes[j]
    print(f"\n=== {lbl} ===")
    for name, fn, col in SERIES:
        k = s & fn(a)
        h, _ = np.histogram(a["mis"][k], bins=BINS, weights=a["weight"][k])
        h2, _ = np.histogram(a["mis"][k], bins=BINS, weights=a["weight"][k]**2)
        ax.errorbar(CENT, h, yerr=np.sqrt(h2), fmt="o-", ms=2.5, lw=1.2, color=col,
                    elinewidth=0.8, label=f"{name}  N={h.sum():,.0f}")
        # peak-over-sideband: local linear-free estimate, sideband = +-40..100 MeV
        bkg = h[sb].mean() * pk.sum()
        exc = h[pk].sum() - bkg
        err = np.sqrt(h2[pk].sum() + h2[sb].sum() * (pk.sum() / max(sb.sum(), 1))**2)
        print(f"  {name:42} peak={h[pk].sum():>9,.0f}  sb-pred={bkg:>9,.0f}  "
              f"excess={exc:>9,.0f} +- {err:>7,.0f}  ({exc/max(err,1e-9):>5.1f} sigma)")
    ax.axvline(M_D, color="navy", ls=":", lw=1.4)
    ax.set_xlim(BINS[0], BINS[-1])
    ax.set_ylim(bottom=0)
    ax.set_title(lbl, fontsize=13, pad=6)
    ax.set_xlabel(r"M(K$^+$K$^-\pi^\pm$), K$^+\to\pi$  [GeV]")
    ax.legend(fontsize=8, loc="upper right")
    if j == 0:
        ax.set_ylabel(f"W+jets events / {WIDTH*1000:.0f} MeV")
        # data=False ALREADY renders "Simulation"; passing text="Simulation" on top of
        # it prints "CMS Simulation Simulation" (mplhep signature is label(text, label)).
        hep.cms.label(ax=ax, data=False,
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

fig.suptitle("GenTag mode 8 = partial truth-matching: does it hide the D misID reflection?",
             fontsize=15)
out = common.fig_path("misid", "mode8_misid_peak.png", key=MODEL_KEY)
fig.savefig(out, dpi=140, bbox_inches="tight")
print("\nSaved", out)

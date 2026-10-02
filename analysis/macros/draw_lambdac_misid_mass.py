"""KpToP 互换后的 M(KK pi): OS data vs W+jets GenTag，找 Lambda_c+ -> p K- pi+ 反射。

动机：上侧 sideband 已由 D->Kpipi（pi 被当成 K）解释（draw_misid_mass_data_vs_wjets.py），
下侧 sideband 的 data 超出没有解释。能把质量往**下**拉的 misID 需要一个比 K 重的真粒子
被套上 K 质量，最自然的就是 p->K。Lambda_c+ -> p K- pi+（BR ~6%）的电荷正好对上 OS：
W- -> mu-、c -> Lambda_c+，p+ 与 pi+ 都与 muon 反号、K- 与 muon 同号，所以 p 落在 Kp slot。
把 Kp 换回 proton 质量就该复原出 Lambda_c 峰（2.28646）。

对照：同样的互换施加到 Km slot（与 muon 同号），Lambda_c 不会把 p 放在那里，应当无峰。

互换配方沿用 draw_misid_mass_data_vs_wjets.py：用原始动量分量，
    E = sum sqrt(px^2+py^2+pz^2+m^2)
SR 在原始 Ds_mass 上挖，两个样本同样施加（理由见那个脚本的 docstring）。
绝对产额，不归一化。

MC 堆叠额外单列「charm-baryon 起源」：三条径迹任一条的 nearest-HF 祖先是 4xxx 重子
（GenTag 的 mode 只认 431/411/421，所以全匹配的 Lambda_c 落 mode 6、部分匹配落 mode 8/0，
这里把它们从各自 mode 里拿出来）。一次运行同时出 Kp 互换和 Km 对照两张图。

    conda activate ds_dnn && python3 draw_lambdac_misid_mass.py
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
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, SR_HW = 1.96835, 0.040
M_LC = 2.28646
M_XIC = 2.46771   # Xi_c+ -> p K- pi+ lands on the same swap
M_K, M_PI, M_P = 0.493677, 0.13957039, 0.93827209
# 40 MeV: 10 and 20 MeV both left W+jets MC with too few events per bin past eps 70%.
# Edges are aligned so M_LC sits at a bin CENTRE -> the +-20 MeV peak window is exactly one bin.
BINS = np.arange(M_LC - 0.020 - 0.040 * 9, 2.901, 0.040)
CENT = 0.5 * (BINS[1:] + BINS[:-1])

# The Km swap is the control: Lambda_c+ never puts its proton in that slot.
# One MC histogram, no mode split (user, 0928): per-mode stacks were unreadable once
# W+jets runs out of events past eps 70%. The charm-baryon share is still printed.
GROUPS = [(r"W+jets MC", "all", "#5a9e78")]
HF = [f"{p}_gen_hf_pdgId" for p in ("Kp", "Km", "pi")]
SEL = [*common.SELECTION]
MOM = [f"{p}_p{x}" for p in ("Kp", "Km", "pi") for x in ("x", "y", "z")]


def swap_mass(a, swap):
    """M(KK pi) with the `swap` track re-assigned the proton mass."""
    mass = {"Kp": M_K, "Km": M_K, "pi": M_PI}
    mass[swap] = M_P
    e = sum(np.sqrt(a[f"{p}_px"]**2 + a[f"{p}_py"]**2 + a[f"{p}_pz"]**2 + m**2)
            for p, m in mass.items())
    px, py, pz = (sum(a[f"{p}_p{x}"] for p in mass) for x in ("x", "y", "z"))
    return np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))


def load(name, extra=()):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        SEL + MOM + ["Ds_mass", "score_dsmu"] + list(extra), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HW
    out = {k: v[m].astype(np.float64) if v.dtype.kind == "f" else v[m] for k, v in a.items()}
    for sw in ("Kp", "Km"):
        out[f"mis_{sw}"] = swap_mass(out, sw)
    if HF[0] in out:
        hf = np.stack([np.abs(out[b]) for b in HF])
        out["baryon"] = ((hf > 4000) & (hf < 5000)).any(axis=0)
    print(f"  {name:22} {m.sum():>9,} rows")
    return out


def peak_excess(h, var):
    """Signal window +-20 MeV around M_LC, background from 40-100 MeV sidebands
    on both sides (40-120 MeV), scaled by width (3 bins per side-bin pair are equal width)."""
    sig = np.abs(CENT - M_LC) < 0.020
    sb = (np.abs(CENT - M_LC) > 0.040) & (np.abs(CENT - M_LC) < 0.110)
    s, b = h[sig].sum(), h[sb].sum() * sig.sum() / sb.sum()
    vs, vb = var[sig].sum(), var[sb].sum() * (sig.sum() / sb.sum()) ** 2
    return s - b, np.sqrt(vs + vb), b


print("Loading (selection + score>-1 + Ds-SR removed on ORIGINAL mass)...")
d = load("Data_OS_merged")
w = load("WtoLNu2Jets_merged", ["weight", "cand_gen_mode"] + HF)

COLS = [("no cut", None), (r"$\varepsilon_{sig}$=90%", 0.90),
        (r"$\varepsilon_{sig}$=70%", 0.70), (r"$\varepsilon_{sig}$=50%", 0.50)]
THR = {lbl: None if eff is None else common.stage1_threshold(
    dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))[0] for lbl, eff in COLS}

# Where does MC's charm-baryon component sit in the ORIGINAL mass? That is what matters
# for the lower sideband, independent of whether a peak shows up in the swapped mass.
WIN = [("low <1.830", w["Ds_mass"] < 1.830),
       ("D+- 1.83-1.91", (w["Ds_mass"] >= 1.830) & (w["Ds_mass"] < 1.910)),
       ("high >2.008", w["Ds_mass"] > M_DS + SR_HW)]
for lbl, _eff in COLS:
    sw = np.ones(len(w["weight"]), bool) if THR[lbl] is None else w["score_dsmu"] > THR[lbl]
    row = []
    for wn, wm in WIN:
        tot_w = w["weight"][sw & wm].sum()
        bar_w = w["weight"][sw & wm & w["baryon"]].sum()
        row.append(f"{wn}: {bar_w:8.1f}/{tot_w:10.1f} ({100*bar_w/max(tot_w,1e-9):.2f}%)")
    print(f"MC charm-baryon origin by ORIGINAL mass, {lbl:>22} | " + " | ".join(row))

for SWAP in ("Kp", "Km"):
  CONTROL = SWAP == "Km"
  d["mis"], w["mis"] = d[f"mis_{SWAP}"], w[f"mis_{SWAP}"]
  print(f"=== swap {SWAP} -> proton {'(CONTROL)' if CONTROL else ''} ===")
  fig, axes = plt.subplots(2, len(COLS), figsize=(6.2 * len(COLS), 7.4),
                           gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06, wspace=0.22))
  width = BINS[1] - BINS[0]
  for j, (lbl, _eff) in enumerate(COLS):
      thr = THR[lbl]
      ax, rax = axes[0][j], axes[1][j]
      sd = np.ones(len(d["mis"]), bool) if thr is None else d["score_dsmu"] > thr
      sw = np.ones(len(w["mis"]), bool) if thr is None else w["score_dsmu"] > thr

      hs, h2s = [], []
      for _n, modes, _c in GROUPS:
          k = sw
          hs.append(np.histogram(w["mis"][k], bins=BINS, weights=w["weight"][k])[0])
          h2s.append(np.histogram(w["mis"][k], bins=BINS, weights=w["weight"][k] ** 2)[0])
      tot, tot2 = np.sum(hs, axis=0), np.sum(h2s, axis=0)
      tot_err = np.sqrt(tot2)

      bottom = np.zeros(len(CENT))
      for (lab, _m, cl), h in zip(GROUPS, hs):
          ax.bar(CENT, h, width=width, bottom=bottom, color=cl, alpha=0.9,
                 label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
          bottom = bottom + h
      ax.bar(CENT, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
             hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3, label="W+jets MC stat")

      hd = np.histogram(d["mis"][sd], bins=BINS)[0]
      n_out = sd.sum() - hd.sum()
      ax.errorbar(CENT, hd, yerr=np.sqrt(hd), fmt="o", ms=4, color="black", elinewidth=1.2,
                  label=f"OS data (Ds-SR excl.)  N={hd.sum():,.0f}")
      for mref in (M_LC, M_XIC):
          ax.axvline(mref, color="navy", ls=":", lw=1.4)
          rax.axvline(mref, color="navy", ls=":", lw=1.4)

      ax.set_xlim(BINS[0], BINS[-1])
      peak = max(bottom.max(), hd.max())
      ax.set_ylim(0, peak * 1.5 if peak > 0 else 1.0)
      ax.tick_params(labelbottom=False)
      ax.legend(fontsize=8, loc="upper right")
      ax.set_title(lbl, fontsize=13, pad=6)
      if j == 0:
          ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
          hep.cms.label(ax=ax, data=True, text="Preliminary",
                        rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

      with np.errstate(invalid="ignore", divide="ignore"):
          r = np.where(tot > 0, hd / tot, np.nan)
          re = np.where(tot > 0, np.sqrt(hd) / tot, np.nan)
          band = np.where(tot > 0, tot_err / tot, np.nan)
      rax.bar(CENT, 2 * band, width=width, bottom=1 - band, color="none",
              hatch="////", edgecolor="0.35", linewidth=0.0)
      rax.errorbar(CENT, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
      rax.axhline(1.0, color="red", lw=1.0)
      rax.set_ylim(0.0, 1.99)
      rax.set_yticks([0.5, 1.0, 1.5])
      rax.set_xlim(BINS[0], BINS[-1])
      if j == 0:
          rax.set_ylabel("data / MC")

      ed, sed, bd = peak_excess(hd.astype(float), hd.astype(float))
      em, sem, bm = peak_excess(tot, tot2)
      kb = sw & w["baryon"]
      per6 = peak_excess(np.histogram(w["mis"][kb], bins=BINS, weights=w["weight"][kb])[0],
                         np.histogram(w["mis"][kb], bins=BINS, weights=w["weight"][kb] ** 2)[0])[0]
      print(f"{lbl:>22}: out-of-range data {n_out:,}/{sd.sum():,} | Lc +-20MeV excess: "
            f"data {ed:8.1f}+-{sed:6.1f} ({ed/max(sed,1e-9):5.1f}s, S/B={ed/max(bd,1e-9):.3f})  "
            f"MC {em:8.1f}+-{sem:6.1f} ({em/max(sem,1e-9):5.1f}s, S/B={em/max(bm,1e-9):.3f}, "
            f"baryon part {per6:.1f})  data/MC={ed/em if em else float('nan'):.3f}")

  fig.supxlabel(rf"M(K$^+$K$^-\pi^\pm$) with the {'K$^-$' if CONTROL else 'K$^+$'} track "
                r"re-assigned the proton mass  [GeV]", fontsize=15, y=0.02)
  fig.suptitle(r"misID reflection: K$\to$p swap "
               + ("(CONTROL: same-sign-as-$\\mu$ slot)" if CONTROL else r"-- $\Lambda_c^+\to$pK$^-\pi^+$")
               + " -- OS data vs W+jets", fontsize=15)
  out = common.fig_path("misid", f"misid_{SWAP}ToP_mass_data_vs_wjets.png", key=MODEL_KEY)
  fig.savefig(out, dpi=140, bbox_inches="tight")
  print("Saved", out)

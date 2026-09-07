"""挖掉 misID D+- 反射峰之后, 重看 Ds_mass 上的 data vs W+jets。

    python3 check_misid_veto.py            # bins, veto 半宽 40 MeV
    python3 check_misid_veto.py incl 0.02

veto 定义: 丢掉 |Ds_mass_KpToPi - 1.8696| < VETO_HW 的候选 —— 即那些「把 K+ 换成 pi+
之后落在 D+- 质量上」的候选, 它们是真实的 D+- -> K-+ pi+- pi+- 被误指认进来的。

要一起回答的两件事:
  1. 挖掉之后 data/MC 闭合有没有变好? MC 高估这个反射约 2 倍(双比值 0.43-0.54),
     所以挖掉它应该把高质量端的过预言拉回来一些。
  2. **代价是多少信号?** 信号的 misID 质量中位数 1.764、95 分位 1.861, 与 D+- 窗口
     重叠, 所以这个 veto 一定会吃掉信号。不报这个数就没法判断 veto 能不能用。

Ds_mass 是原始质量。OS data 的 SR 被 blind, W+jets 不 blind ——
按画图标准, 谱图上只有 OS data 挖 SR, W+jets 画满并把 blind 区涂灰;
产额和比值一律只在 sideband 上引用, 两边同样处理。
"""
import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt

import mplhep as hep
plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

MODEL = "ssdata"
PRED = f"{common.out_dir_for(MODEL)}/Prediction"

M_DS, SR_HW, M_D = 1.96835, 0.04, 1.8696
LOW_MAX, HIGH_MIN = 1.830, 2.008
EFFS = [1.0, 0.90, 0.70, 0.50]
MIN_NEFF = 30      # 低于此值加权直方图会出假结构(负权重抵消), 表里标记出来

argv = sys.argv[1:]
SRC = "bins"
VETO_HW = 0.04
BW = 0.01
SIGN = "OS"
PHI_HW = None
for a in list(argv):
    if a in ("bins", "incl"):
        SRC = a
    elif a in ("OS", "SS"):
        SIGN = a
    elif a.startswith("phi="):
        PHI_HW = float(a.split("=")[1])
    elif a.startswith("bw="):
        BW = float(a.split("=")[1])
    elif a.startswith("--effs"):
        EFFS = [float(x) for x in a.split("=", 1)[1].split(",")]
    else:
        VETO_HW = float(a)

PHITAG = "" if PHI_HW is None else f"_phi{PHI_HW*1000:.0f}"
TAG = f"_{SIGN}{PHITAG}" + ("" if EFFS == [1.0, 0.90, 0.70, 0.50] else "_loweff") + ("" if BW == 0.01 else f"_bw{BW*1000:.0f}")

ERAS = {"2022 (preEE)": ("2022", "7.99 fb$^{-1}$"),
        "2022EE (postEE)": ("2022EE", "26.68 fb$^{-1}$")}

M_PHI = 1.019460
D = {}
for era, (tag, _) in ERAS.items():
    D[era] = {
        "wj": uproot.open(f"{PRED}/WtoLNu2Jets_{SRC}_{tag}{'' if SIGN=='OS' else '_SS'}_misID.root")["Events"].arrays(library="np"),
        "data": uproot.open(f"{PRED}/Data_{tag}{'' if SIGN=='OS' else '_SS'}_misID.root")["Events"].arrays(library="np"),
    }
    if PHI_HW is not None:
        for k in ("wj", "data"):
            keep = np.abs(D[era][k]["phi_mass"] - M_PHI) < PHI_HW
            D[era][k] = {kk: vv[keep] for kk, vv in D[era][k].items()}

# 信号: 走 Prediction 主文件(misID 文件只为 W+jets/data 建), 就地算 KpToPi。
M_K, M_PI = 0.493677, 0.13957039
MOM = [f"{p}_{c}" for p in ("Kp", "Km", "pi") for c in ("px", "py", "pz")]
sp = []
for sn in common.SIGNALS:
    a = uproot.open(f"{PRED}/{sn}.root")["Events"].arrays(
        MOM + ["score_dsmu", "Ds_match", "mu_match", "num_os_mu_tight",
               "num_ss_mu_tight", common.B_VETO], library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1)
         & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool))
    sp.append({k: v[m] for k, v in a.items()})
S = {k: np.concatenate([p[k] for p in sp]) for k in sp[0]}
g = {k: S[k].astype(np.float64) for k in MOM}
e = sum(np.sqrt(g[f"{p}_px"]**2 + g[f"{p}_py"]**2 + g[f"{p}_pz"]**2 + m**2)
        for p, m in [("Kp", M_PI), ("Km", M_K), ("pi", M_PI)])
px = sum(g[f"{p}_px"] for p in ("Kp", "Km", "pi"))
py = sum(g[f"{p}_py"] for p in ("Kp", "Km", "pi"))
pz = sum(g[f"{p}_pz"] for p in ("Kp", "Km", "pi"))
S["mis"] = np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))
THR = {x: (-1.0 if x >= 1.0 else float(np.quantile(S["score_dsmu"], 1 - x))) for x in EFFS}

sig_keep = np.abs(S["mis"] - M_D) >= VETO_HW
print(f"veto: |M(KKpi with K+->pi+) - {M_D}| < {VETO_HW*1000:.0f} MeV   source={SRC}")
print(f"\n信号代价 (truth-matched, {len(S['mis']):,} 行):")
for x in EFFS:
    k = S["score_dsmu"] > THR[x]
    lbl = "no cut" if x >= 1.0 else f"{x:.0%}"
    print(f"  eps_sig={lbl:>6}: veto 保留 {sig_keep[k].mean()*100:5.2f}% 的信号  "
          f"(丢 {100-sig_keep[k].mean()*100:.2f}%)")


def cell(era, x, veto):
    thr = THR[x]
    out = {}
    for kind in ("data", "wj"):
        s = D[era][kind]
        m = (s["score_dsmu"] > thr) & (np.abs(s["Ds_mass"] - M_DS) >= SR_HW)
        if veto:
            m &= np.abs(s["Ds_mass_KpToPi"] - M_D) >= VETO_HW
        w = s["Lumi_reweight"][m] if kind == "wj" else None
        mass = s["Ds_mass"][m]
        n = w.sum() if w is not None else float(len(mass))
        v = (w**2).sum() if w is not None else float(len(mass))
        lo = (w[mass < LOW_MAX].sum() if w is not None else (mass < LOW_MAX).sum())
        hi = (w[mass > HIGH_MIN].sum() if w is not None else (mass > HIGH_MIN).sum())
        lo2 = ((w[mass < LOW_MAX]**2).sum() if w is not None else lo)
        hi2 = ((w[mass > HIGH_MIN]**2).sum() if w is not None else hi)
        out[kind] = (n, v, lo, hi, lo2, hi2)
    (nd, vd, dl, dh, dl2, dh2) = out["data"]
    (nw, vw, wl, wh, wl2, wh2) = out["wj"]
    r = nd / nw
    re_ = r * np.sqrt(vd / nd**2 + vw / nw**2)
    dr = (dl / wl) / (dh / wh)
    dre = dr * np.sqrt(dl2 / dl**2 + dh2 / dh**2 + wl2 / wl**2 + wh2 / wh**2)
    return r, re_, dr, dre, nd, nw, nw**2 / vw


print("\n" + "=" * 100)
print("OS sideband, 绝对归一 (W+jets 带 sigma*L, 未归到 data)")
print("=" * 100)
print(f"{'era':<17}{'eps':>7} | {'N_eff(MC,veto)':>14} | {'data/MC 无veto':>18} "
      f"{'data/MC 有veto':>18} | {'低/高 无veto':>16} {'低/高 有veto':>16}")
for era in ERAS:
    for x in EFFS:
        a0 = cell(era, x, False)
        a1 = cell(era, x, True)
        lbl = "no cut" if x >= 1.0 else f"{x:.0%}"
        flag = " <-N_eff低" if a1[6] < MIN_NEFF else ""
        print(f"{era:<17}{lbl:>7} | {a1[6]:14,.0f} | {a0[0]:8.3f}+-{a0[1]:<7.3f} "
              f"{a1[0]:8.3f}+-{a1[1]:<7.3f} | {a0[2]:7.3f}+-{a0[3]:<7.3f} "
              f"{a1[2]:7.3f}+-{a1[3]:<7.3f}{flag}")
    print()

# ---------- 图 ----------
mbins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, BW)   # blind 窗仍落在 bin 边界
mcent_mass = 0.5 * (mbins[1:] + mbins[:-1])
fig, axes = plt.subplots(len(ERAS), len(EFFS),
                         figsize=(8.5 * len(EFFS), 7.5 * len(ERAS)), squeeze=False)
for row, (era, (_, lumi)) in enumerate(ERAS.items()):
    for col, x in enumerate(EFFS):
        ax = axes[row][col]
        thr = THR[x]
        w = D[era]["wj"]
        kw = (w["score_dsmu"] > thr) & (np.abs(w["Ds_mass_KpToPi"] - M_D) >= VETO_HW)
        hw, _ = np.histogram(w["Ds_mass"][kw], mbins, weights=w["Lumi_reweight"][kw])
        hw2, _ = np.histogram(w["Ds_mass"][kw], mbins, weights=w["Lumi_reweight"][kw] ** 2)
        d = D[era]["data"]
        kd = (d["score_dsmu"] > thr) & (np.abs(d["Ds_mass_KpToPi"] - M_D) >= VETO_HW)
        if SIGN == "OS":
            kd = kd & (np.abs(d["Ds_mass"] - M_DS) >= SR_HW)   # OS data 才 blind
        hd, _ = np.histogram(d["Ds_mass"][kd], mbins)
        hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     edgecolor="darkgreen", linewidth=0.8,
                     label=f"W+jets {SRC} {SIGN} MC  N={hw.sum():,.0f}")
        # mplhep's histtype="fill" silently ignores yerr (only "step"/"errorbar"/
        # "band" draw it -- checked in the installed mplhep source), so the W+jets
        # uncertainty has to be a separate manual errorbar on top of the fill.
        ax.errorbar(mcent_mass, hw, yerr=np.sqrt(hw2), fmt="none", ecolor="darkgreen",
                    elinewidth=1.0, capsize=0, alpha=0.8)
        hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=np.sqrt(hd),
                     color="black", marker="s", markersize=4.5, elinewidth=1.4,
                     label=f"{SIGN} data{' (SR blind)' if SIGN=='OS' else ' (odd parity)'}  N={hd.sum():,.0f}")
        ax.axvspan(M_DS - SR_HW, M_DS + SR_HW, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        top = ax.get_ylim()[1]
        ax.text(M_D, top * 0.97, r"D$^\pm$", color="green", ha="center", va="top", fontsize=13)
        ax.text(M_DS, top * 0.97, r"D$_s^\pm$", color="dimgrey", ha="center", va="top", fontsize=13)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel(f"Events / {BW*1000:.0f} MeV")
        lbl = "no cut" if x >= 1.0 else rf"$\varepsilon_{{sig}}$={100*x:.0f}%"
        ax.set_title(f"{era} -- {lbl}\nmisID D$^\\pm$ veto "
                     rf"($\pm${VETO_HW*1000:.0f} MeV)"
                     + ("" if PHI_HW is None else rf",  $|M_{{KK}}-M_\phi|<${PHI_HW*1000:.0f} MeV"),
                     fontsize=12)
        ax.legend(fontsize=11, loc="upper right")
        if col == 0:
            hep.cms.label(ax=ax, rlabel=f"{lumi} (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path("misid", f"misid_veto_{SRC}{TAG}_mass.png", key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

# ---------- 逐 bin data/MC, 绝对归一, 只画 veto 之后 ----------
# 不重归一化: W+jets 带 sigma*L, 高度本身就是被检验的量。
# veto 前后叠画过一版, 太挤; 前后对比看上面那张数值表(低/高双比值那两列)。
# SR 里 OS data 是空的(blind), 比值无定义, 自动落在 ok 掩码外。
RATIO_YMAX = 2.5
mcent = 0.5 * (mbins[1:] + mbins[:-1])
fig, axes = plt.subplots(len(ERAS), len(EFFS),
                         figsize=(8.5 * len(EFFS), 6.0 * len(ERAS)), squeeze=False)
for row, era in enumerate(ERAS):
    for col, x in enumerate(EFFS):
        ax = axes[row][col]
        thr = THR[x]
        for veto, col_, mk, lab in [(True, "black", "s",
                                     f"veto $\\pm${VETO_HW*1000:.0f} MeV")]:
            w, d = D[era]["wj"], D[era]["data"]
            kw = w["score_dsmu"] > thr
            kd = (d["score_dsmu"] > thr) & (np.abs(d["Ds_mass"] - M_DS) >= SR_HW)
            if veto:
                kw = kw & (np.abs(w["Ds_mass_KpToPi"] - M_D) >= VETO_HW)
                kd = kd & (np.abs(d["Ds_mass_KpToPi"] - M_D) >= VETO_HW)
            hw, _ = np.histogram(w["Ds_mass"][kw], mbins, weights=w["Lumi_reweight"][kw])
            hw2, _ = np.histogram(w["Ds_mass"][kw], mbins, weights=w["Lumi_reweight"][kw] ** 2)
            hd, _ = np.histogram(d["Ds_mass"][kd], mbins)
            ok = (hd > 0) & (hw > 0)
            r = np.where(ok, hd / np.where(hw > 0, hw, 1), np.nan)
            re_ = np.where(ok, r * np.sqrt(1 / np.maximum(hd, 1)
                                           + hw2 / np.maximum(hw, 1e-9) ** 2), np.nan)
            ins = ok & (r <= RATIO_YMAX)
            off = 0.0
            ax.errorbar(mcent[ins] + off, r[ins], yerr=re_[ins], fmt=mk, color=col_,
                        markersize=4.5, elinewidth=1.3, label=lab)
            ab = ok & (r > RATIO_YMAX)
            if ab.any():
                ax.plot(mcent[ab] + off, np.full(ab.sum(), RATIO_YMAX * 0.97), "^",
                        color=col_, markersize=7)
        ax.set_ylim(0, RATIO_YMAX)
        ax.axhline(1.0, color="firebrick", lw=1.4)
        ax.axvspan(M_DS - SR_HW, M_DS + SR_HW, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        if col == 0:
            ax.set_ylabel(f"data / W+jets {SRC}")
        lbl = "no cut" if x >= 1.0 else rf"$\varepsilon_{{sig}}$={100*x:.0f}%"
        ax.set_title(f"{era} -- {lbl}", fontsize=15)
        ax.legend(fontsize=11, loc="upper left")
plt.tight_layout()
out = common.fig_path("misid", f"misid_veto_{SRC}{TAG}_ratio.png", key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"Saved {out}")

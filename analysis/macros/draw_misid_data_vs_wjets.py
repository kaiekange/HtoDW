"""misID 后的质量谱: OS data vs W+jets, 按 era 和 score cut 铺开。

    python3 draw_misid_data_vs_wjets.py            # bins + cr (默认)
    python3 draw_misid_data_vs_wjets.py incl high

横轴是 KpToPi 互换后的 M(KKpi): 真实的 D+- -> K-+ pi+- pi+- 里有一个 pi 被误指认成 K
才进了我们的 KKpi 候选池, 把它换回 pi 质量就复原出 D+- 峰。所以这张图上 1.8696 处的
峰**就是**那批被错指认的真 D+-。

两件必须一起做的事:

* **SS 与 OS 的两处结构性差异**: (1) SS data **没有 blind**(SS 结构上排除信号),
  所以那里的 SR 不是空的; (2) SS data 有一半(偶数 event_id)参与了 ssdata 模型训练,
  画图只用**奇数半**, 否则是自证自洁 —— 这一步在 make_misid.py 建文件时就做掉了。
  W+jets SS 从未参与训练, 不切 parity。region 仍照 OS 施加, 只为让两者可比。

* **原始质量上必须挖掉 SR, 三个样本同样处理。** OS data 的 SR 被 blind 而 W+jets
  没有; blind 在原始质量上是干净窗口, 映射到 misID 质量上却是弥散的, 所以不挖就是
  data 在 misID 谱一大片区域里系统性缺事例而 MC 不缺。默认 region=cr 只挖 SR、两侧
  sideband 都留 —— 那是 blind 唯一真正要求的东西。region=high 另留了一个只看高侧的
  变体(数值抽取用, 反射几乎全在高侧), 但对画图是过度限制。

* **绝对产额, 不归一化。** W+jets 带 sigma*L, 它的高度本身就是被检验的量; 归到 data
  等于把要测的东西变成一条贴着 1 的平线。

**只画 KpToPi。** KmToPi (wrong-sign) 是数值上的对照, 但**不能画在同一张图上**:
两种指认的连续本底形状本身就不同, 叠上去眼睛会把形状差异读成峰 —— 实测 KmToPi 的
峰超出是 -0.5/1.4 sigma(无峰), 与图上的印象相反。对照结论只从数值表读。
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
M_DS, SR_HW = 1.96835, 0.04
M_D = 1.8696
PK = (1.8496, 1.8896)
EFFS = [1.0, 0.90, 0.70, 0.50]

argv = sys.argv[1:]
SRC, REGION, VAR, BW, SIGN, PHI_HW = "bins", "cr", "kp", 0.01, "OS", None
for a in list(argv):
    if a in ("bins", "incl"):
        SRC = a
        argv.remove(a)
    elif a in ("cr", "high", "low"):
        REGION = a
        argv.remove(a)
    elif a in ("kp", "both"):
        VAR = a
        argv.remove(a)
    elif a in ("OS", "SS"):
        SIGN = a
        argv.remove(a)
    elif a.startswith("phi="):
        PHI_HW = float(a.split("=")[1])
        argv.remove(a)
    elif a.startswith("bw="):
        BW = float(a.split("=")[1])
        argv.remove(a)

# 原始 Ds_mass 上的区域。**必须**加一个, 因为 OS data 的 SR 被 blind 而 W+jets 没有,
# 而 blind 在原始质量上是干净窗口、映射到 misID 质量上却是弥散的 -- 不挖就是 data 在
# misID 谱一大片区域里系统性缺事例而 MC 不缺。三个样本一律同样处理。
#
#   cr    完整 control region: 只挖 SR, 两侧 sideband 都留。**默认, 也是 blind 唯一
#         真正要求的东西。**
#   high  只留高侧 (Ds_mass > 2.00835)。这是数值抽取时的选择 -- 反射的原始质量中位数
#         2.051、5/25 分位 1.974/2.009, 只留高侧能装下绝大部分反射, 峰最纯。对画图
#         是过度限制: 它把整个低侧 sideband 的统计量白扔掉, 而 blind 并不要求这么做。
#   low   只留低侧 (Ds_mass < 1.92835)。high 的补集, 两者拼起来就是 cr。这一档是
#         **对照**: 反射在原始质量上的中位数是 2.051、5 分位 1.974, 所以低侧几乎不含
#         反射, D+- 峰应当明显更弱。若它和 high 一样强, "反射落在高侧" 这个说法就错了。
REGION_DEF = {
    "cr":   (lambda m: np.abs(m - M_DS) >= SR_HW,
             rf"$|M-M_{{D_s}}|\geq${SR_HW*1000:.0f} MeV (SR masked)"),
    "high": (lambda m: m > M_DS + SR_HW,
             rf"$M>${M_DS + SR_HW:.5f} GeV (high sideband)"),
    "low":  (lambda m: m < M_DS - SR_HW,
             rf"$M<${M_DS - SR_HW:.5f} GeV (low sideband)"),
}
REG_FN, REG_LABEL = REGION_DEF[REGION]

ERAS = {"2022 (preEE)": ("2022", "7.99 fb$^{-1}$"),
        "2022EE (postEE)": ("2022EE", "26.68 fb$^{-1}$")}

# bin 宽 10 MeV(2026-08-25 实测后从 20 改过来): D+- 反射峰在 20 MeV 下只跨 **1 个 bin**,
# 那不叫分辨出峰; 10 MeV 下跨 3 个(FWHM), 而峰顶 bin 相对误差在主工作点仍 <10%
# (eps_sig=90%: 6.1-9.9%)。峰高/连续本底比也随之从 5.90 升到 6.32 —— 20 MeV 在把峰
# 抹进邻近连续本底, 是在损失真实显著性而不只是难看。5 MeV 进一步到 5 个 bin, 但
# eps_sig=50% 时 preEE 峰顶误差到 32%, 不划算。
#
# 注意这只针对 misID 谱(峰窄)。原始 Ds_mass 的图仍用 20 MeV 标准 —— 那里的限制是
# 统计量而非分辨率, 且 20 MeV 是 compare_models_mass.png 定下的画图标准。
#
# 哪一种误指认。x 轴范围随之变: 每换一个 kaon 质量, 谱就整体左移约 0.2 GeV
# (truth-matched 信号 3body 中位数 1.9684 -> 1.7637 -> 1.5424)。
#   kp    K+ -> pi+ 单互换。对应 D+- -> K-+ pi+- pi+-, 净电荷等于候选电荷。**有峰。**
#   both  两个 kaon 都当 pion。对应真三 pi 末态。结构上可达(6.5% 的行越过 1.8696),
#         但 phi 窗映射到 pipi 假设下只有 M(pipi) ~ 0.30-0.47, 而 D+- -> pi pi pi
#         主要走 rho(770)pi -> 落在 rho 窗内 0.00%。实测无峰(六格全 <1 sigma)。
VARS = {
    "kp":   ("Ds_mass_KpToPi",   r"K$^+\to\pi^+$",              (1.50, 2.10)),
    "both": ("Ds_mass_BothToPi", r"K$^+$K$^-\to\pi^+\pi^-$",   (1.25, 1.95)),
}
MVAR, VAR_LABEL, (XLO, XHI) = VARS[VAR]

mbins = np.arange(XLO, XHI + 1e-9, BW)
mcent = 0.5 * (mbins[1:] + mbins[:-1])


def load(name):
    a = uproot.open(f"{PRED}/{name}_misID.root")["Events"].arrays(library="np")
    m = REG_FN(a["Ds_mass"])
    # phi 窗收紧: tuplizer 自带的是 M(KK) < Mphi+0.05 (83 MeV 宽)。真 phi 是窄峰几乎
    # 全保留; 而 D->Kpipi 反射的 "M(KK)" 是给 pi 套 K 质量算出来的假质量, 在窗内近似平坦,
    # 所以收紧会按窗宽比例压低它 —— 这本身就是一个把反射与真 phi 分开的手段。
    if PHI_HW is not None:
        m &= np.abs(a["phi_mass"] - 1.019460) < PHI_HW
    return {k: v[m] for k, v in a.items()}


D = {}
for era, (tag, _) in ERAS.items():
    sfx = "" if SIGN == "OS" else "_SS"
    D[era] = dict(wj=load(f"WtoLNu2Jets_{SRC}_{tag}{sfx}"), data=load(f"Data_{tag}{sfx}"))
    print(f"{era:18s} W+jets {len(D[era]['wj']['Ds_mass']):>9,} | "
          f"OS data {len(D[era]['data']['Ds_mass']):>9,} rows  [{REGION}]")

sig = []
for sn in common.SIGNALS:
    d = uproot.open(f"{PRED}/{sn}.root")["Events"].arrays(
        ["score_dsmu", "Ds_match", "mu_match", *common.SELECTION], library="np")
    m = (common.selection_mask(d) & (d["score_dsmu"] > -1)
         & d["Ds_match"].astype(bool) & d["mu_match"].astype(bool))
    sig.append(d["score_dsmu"][m])
sig = np.concatenate(sig)
THR = {e: (-1.0 if e >= 1.0 else float(np.quantile(sig, 1 - e))) for e in EFFS}


def hist(s, thr, var, weighted):
    k = s["score_dsmu"] > thr
    w = s["Lumi_reweight"][k] if weighted else None
    h, _ = np.histogram(s[var][k], bins=mbins, weights=w)
    v, _ = np.histogram(s[var][k], bins=mbins,
                        weights=None if w is None else w ** 2)
    return h, np.sqrt(v if w is not None else h)


fig, axes = plt.subplots(len(ERAS), len(EFFS),
                         figsize=(8.5 * len(EFFS), 7.5 * len(ERAS)), squeeze=False)
for row, (era, (_, lumi)) in enumerate(ERAS.items()):
    for col, eff in enumerate(EFFS):
        ax = axes[row][col]
        thr = THR[eff]
        hw, _ = hist(D[era]["wj"], thr, MVAR, True)
        hd, ed = hist(D[era]["data"], thr, MVAR, False)
        hep.histplot(hw, mbins, ax=ax, histtype="fill", alpha=0.5, color="seagreen",
                     label=f"W+jets {SRC} {SIGN} MC  N={hw.sum():,.0f}")
        hep.histplot(hd, mbins, ax=ax, histtype="errorbar", yerr=ed, color="black",
                     marker="s", markersize=4.5, elinewidth=1.4,
                     label=f"{SIGN} data{' (odd parity)' if SIGN=='SS' else ''}  N={hd.sum():,.0f}")
        ax.axvspan(*PK, color="grey", alpha=0.15)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        ax.text(M_D, ax.get_ylim()[1] * 0.97, r"D$^\pm$", color="green",
                ha="center", va="top", fontsize=13)
        ax.set_xlabel(rf"$M(K^+K^-\pi^\pm)$ with {VAR_LABEL}  [GeV]")
        if col == 0:
            ax.set_ylabel(f"Events / {BW*1000:.0f} MeV")
        lbl = "no cut" if eff >= 1.0 else rf"$\varepsilon_{{sig}}$={100*eff:.0f}%"
        ax.set_title(f"{era} -- {lbl}\n{REG_LABEL}" + ("" if PHI_HW is None else rf",  $|M_{{KK}}-M_\phi|<${PHI_HW*1000:.0f} MeV"), fontsize=12)
        ax.legend(fontsize=10, loc="upper left")
        if col == 0:
            hep.cms.label(ax=ax, rlabel=f"{lumi} (13.6 TeV)", data=True,
                          text="Preliminary", fontsize=11)
plt.tight_layout()
out = common.fig_path(
    "misid",
    f"misid_{VAR}_{SRC}_{SIGN}{'' if PHI_HW is None else f'_phi{PHI_HW*1000:.0f}'}_{REGION}_mass{'' if BW==0.01 else f'_bw{BW*1000:.0f}'}.png",
    key=MODEL)
plt.savefig(out, dpi=110)
plt.close(fig)
print(f"\nSaved {out}")

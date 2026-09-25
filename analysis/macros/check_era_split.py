"""postEE 的 data/MC 闭合比 preEE 差 —— 把它拆成 data 侧和 MC 侧, 再问是不是单个
pileup 敏感特征就能复现。

    python3 check_era_split.py            # DNN:ssdata seed42
    python3 check_era_split.py bins       # W+jets 用三个 jet bin (默认) 或 incl

曾经有两个互相竞争的假说, 符号相同, 所以 era split 本身分不开它们:

  (A) 训练集 era 不平衡泄漏。**这条已经死了, 别再引用。** 它的前提是「信号 100% preEE」,
      而 postEE 信号 MC 在 2026-08-25..09-04 已产出并合进 Samples/。2026-09-08 实测
      两个训练类的 era 构成逐位相同 —— 信号(按 Lumi_reweight 加权, 这才是训练用的)
      23.2/76.8, 背景 Data_SS parity=0 也是 23.2/76.8。era 不携带类别信息, 学它没有收益。
  (B) MC 把 pileup/conditions 建模错了。这个假说**不需要**任何训练不平衡。

  (C) 真正还在的偏差, 同样不需要 era 不平衡: **signal 是 MC 而 background 是真 data**,
      于是任何 data/MC 失配都是免费判别量。实测 <Ds_IsoR04_sumPUPt> = 9.65 GeV (signal MC)
      vs 11.99 (preEE data) / 14.03 (postEE data), 而**全仓库没有任何 pileup reweighting**。
      MC 逐位 era-flat 而 data 涨 17%, 所以 data/MC 的 gap 在 preEE 是 4%、postEE 是 22%
      —— 这就是 DNN cut 多砍掉 postEE data 的原因。

「分别给 preEE / postEE 的 SS data 打分, 看到分裂就说明模型学到了 era」**这个读法是错的**,
这一条仍然有效: postEE 本来 pileup 就更高, 而 pileup 会推动 num_jets / jet_ht / isolation
—— 这些是合法判别量。看到分裂是必要条件, 不是充分条件。

能区分的是这里的单特征对照: 只用一个 pileup 敏感特征切, 在**匹配的 MC 效率**上,
重算 era 双比值。2026-09-08 结果: Ds_IsoR04_sumPUPt 单独切给 0.533/0.531/0.570,
**比 DNN 同档的 0.752/0.707/0.686 还极端**, 而对照特征 jet_ht 给 1.310 (反方向)
⇒ pileup 建模误差单独就足够解释全部效应, 且有余。
"""
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HW = 1.96835, 0.04
EFFS = [0.90, 0.70, 0.50]
SOURCE = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ("bins", "incl") else "bins"
WJ = {"bins": "WtoLNu2Jets", "incl": "WtoLNu2Jets_incl"}[SOURCE]

# pileup 敏感的单特征。方向: 全部是「越大越像背景」, 所以 cut 是保留 <= thr。
#
# 两组, 是对照不是挑选。分组依据是 2026-09-08 实测的逐特征 era 漂移
# (preEE->postEE, 以 preEE data sigma 为单位, d(data) vs d(MC)):
#
#   MISMODELLED -- data 动了而 MC 没跟上, 残差最大的三个:
#     Ds_IsoR04_sumPUPt             +0.276 vs -0.001   (残差 +0.276)
#     Ds_IsoR04_sumNeutralHadronEt  +0.272 vs +0.092   (残差 +0.180)
#     n_b_jets_loose_dR04           +0.153 vs +0.012   (残差 +0.141)
#
#   MODELLED (原来那三个, 留作对照) -- data 和 MC 一起动, 残差近零:
#     jet_ht                        -0.035 vs -0.028   (残差 -0.007)
#     num_jets_pt30                 +0.034 vs +0.011   (残差 +0.023)
#     Ds_IsoR04_sumChargedHadronPt  -0.053 vs -0.008   (残差 -0.045)
#
# 预期: 第一组单独切就该复现 DNN 的 ~0.75, 第二组不该。两组一起跑才是受控结论,
# 只跑第一组是「挑到了一个能复现的特征」。
FEATS_MISMODELLED = ["Ds_IsoR04_sumPUPt", "Ds_IsoR04_sumNeutralHadronEt",
                     "n_b_jets_loose_dR04"]
FEATS_MODELLED = ["jet_ht", "num_jets_pt30", "Ds_IsoR04_sumChargedHadronPt"]
FEATS = FEATS_MISMODELLED + FEATS_MODELLED
BR = (["score_dsmu", "Ds_mass", "Lumi_reweight", "era", "num_os_mu_tight",
       "num_ss_mu_tight", common.B_VETO] + FEATS)
BR_SIG = [b for b in BR if b not in ("Lumi_reweight", "era")] + ["Ds_match", "mu_match"]
PRED = f"{common.out_dir_for('ssdata')}/Prediction"


def load(name, *, sig=False, era=None):
    p = f"{PRED}/{name}.root"
    if not os.path.exists(p):
        sys.exit(f"missing {p}")
    a = uproot.open(p)["Events"].arrays(BR_SIG if sig else BR, library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    else:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HW      # OS sideband, 两侧同样处理
    if era is not None:
        m &= a["era"] == era
    return {k: v[m] for k, v in a.items()}


sig_score = np.concatenate([load(n, sig=True)["score_dsmu"] for n in common.SIGNALS])
THR = {e: float(np.quantile(sig_score, 1 - e)) for e in EFFS}

ERAS = [1, 2]                                     # 1 = 2022 preEE, 2 = 2022EE postEE
D = {e: dict(data=load("Data", era=e), wj=load(WJ, era=e)) for e in ERAS}
print(f"W+jets source: {SOURCE}   signal (truth-matched): {len(sig_score):,} rows")
for e in ERAS:
    print(f"  era {e}: data {len(D[e]['data']['Ds_mass']):>9,}   "
          f"W+jets {len(D[e]['wj']['Ds_mass']):>9,} rows")


def eff(d, keep, w=None):
    """(efficiency, its relative error) for a boolean mask over a (weighted) sample."""
    if w is None:
        n, tot = keep.sum(), len(keep)
        return n / tot, (1 / max(n, 1)) ** 0.5
    num, den = w[keep].sum(), w.sum()
    rel = np.sqrt((w[keep] ** 2).sum() / max(num, 1e-9) ** 2
                  + (w ** 2).sum() / max(den, 1e-9) ** 2)
    return num / den, rel


def era_split(mask_fn):
    """mask_fn(sample_dict) -> boolean keep-mask. Returns the era decomposition."""
    out = {}
    for kind, wkey in [("data", None), ("MC", "Lumi_reweight")]:
        e_, r_ = {}, {}
        for era in ERAS:
            s = D[era]["data" if kind == "data" else "wj"]
            w = None if wkey is None else s[wkey]
            e_[era], r_[era] = eff(s, mask_fn(s), w)
        out[kind] = (e_[2] / e_[1],
                     (e_[2] / e_[1]) * np.hypot(r_[2], r_[1]),
                     e_[1], e_[2])
    dr = out["data"][0] / out["MC"][0]
    dre = dr * np.hypot(out["data"][1] / out["data"][0], out["MC"][1] / out["MC"][0])
    return out, dr, dre


print("\n" + "=" * 96)
print("era 分解: postEE/preEE 的效率比, data 侧与 MC 侧分开 (OS sideband)")
print("=" * 96)
print(f"{'cut':>26} {'eps_MC(preEE)':>14} {'data post/pre':>16} {'MC post/pre':>16} "
      f"{'双比值':>16}")
rows = []
for e in EFFS:
    o, dr, dre = era_split(lambda s, t=THR[e]: s["score_dsmu"] > t)
    print(f"{'DNN eps_sig=' + f'{e:.0%}':>26} {o['MC'][2]:14.5f} "
          f"{o['data'][0]:8.3f}+-{o['data'][1]:<6.3f} {o['MC'][0]:8.3f}+-{o['MC'][1]:<6.3f} "
          f"{dr:8.3f}+-{dre:.3f}")
    rows.append((e, o["MC"][2], dr, dre))

print("\n" + "=" * 96)
print("单特征对照: 同一个特征单独切, 阈值调到与上表同一个 preEE MC 效率")
print("(报告的是**实际达到的** eps, 不是目标 eps -- 有单值尖峰的特征切不动)")
print("=" * 96)
for f in FEATS:
    tag = "MISMODELLED" if f in FEATS_MISMODELLED else "MODELLED (对照)"
    print(f"\n--- {f}  (保留 <= thr)   [{tag}]")
    pre = D[1]["wj"]
    v, w = pre[f], pre["Lumi_reweight"]
    order = np.argsort(v)
    cw = np.cumsum(w[order]) / w.sum()
    for e, target, dr_dnn, _ in rows:
        i = np.searchsorted(cw, target)
        if i >= len(order):
            print(f"    target eps_MC={target:.5f}: 特征无法达到 (最大 {cw[-1]:.5f})")
            continue
        thr = v[order[i]]
        o, dr, dre = era_split(lambda s, t=thr, ff=f: s[ff] <= t)
        print(f"    thr={thr:10.4f}  achieved eps_MC={o['MC'][2]:.5f} (target {target:.5f})"
              f"   data {o['data'][0]:.3f}  MC {o['MC'][0]:.3f}"
              f"   双比值 {dr:.3f}+-{dre:.3f}   [DNN 同档 {dr_dnn:.3f}]")

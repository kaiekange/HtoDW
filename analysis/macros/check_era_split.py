"""postEE 的 data/MC 闭合比 preEE 差 —— 把它拆成 data 侧和 MC 侧, 再问是不是单个
pileup 敏感特征就能复现。

    python3 check_era_split.py            # DNN:ssdata seed42
    python3 check_era_split.py bins       # W+jets 用三个 jet bin (默认) 或 incl

两个互相竞争的假说, 符号相同, 所以 era split 本身分不开它们:

  (A) 训练集 era 不平衡泄漏。背景 ~76% 是 postEE SS data 而信号 100% preEE
      (common.py 的 KNOWN BIAS), 网络可以学到「像 postEE ⇒ 像背景」。
  (B) W+jets postEE MC 把 pileup/conditions 建模错了。这个假说**不需要**任何训练
      不平衡, era 完全平衡的训练集也会给出同样的符号。

CLAUDE.md 里写的那个检验(「分别给 preEE / postEE 的 SS data 打分, 看到分裂就说明模型
学到了 era」)**读法是错的**: postEE 本来 pileup 就更高, 而 pileup 会推动
num_jets / jet_ht / isolation —— 这些是合法判别量。看到分裂是必要条件, 不是充分条件。

能区分的是这里的单特征对照: 只用一个 pileup 敏感特征切, 在**匹配的 MC 效率**上,
重算 era 双比值。若单个特征就能复现 MVA 的 ~0.77, 那 MVA 只是 pileup 建模误差的代理,
(A) 就不必要。
真正决定性的检验是另训一个只用 preEE 背景的模型(samples={"Data_2022_SS": 1.0}),
本脚本不做。
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
WJ = {"bins": ["WtoLNu2Jets_0J_{e}", "WtoLNu2Jets_1J_{e}", "WtoLNu2Jets_2J_{e}"],
      "incl": ["WtoLNu2Jets_incl_{e}"]}[SOURCE]

# pileup 敏感的单特征。方向: 这三个都是「越大越像背景」, 所以 cut 是保留 <= thr。
FEATS = ["num_jets_pt30", "jet_ht", "Ds_IsoR04_sumChargedHadronPt"]
BR = (["score_dsmu", "Ds_mass", "Lumi_reweight", "num_os_mu_tight",
       "num_ss_mu_tight", common.B_VETO] + FEATS)
BR_SIG = [b for b in BR if b != "Lumi_reweight"] + ["Ds_match", "mu_match"]
PRED = f"{common.out_dir_for('ssdata')}/Prediction"


def load(names, *, sig=False):
    if isinstance(names, str):
        names = [names]
    ps = []
    for n in names:
        p = f"{PRED}/{n}.root"
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        ps.append(uproot.open(p)["Events"].arrays(BR_SIG if sig else BR, library="np"))
    a = ps[0] if len(ps) == 1 else {k: np.concatenate([q[k] for q in ps]) for k in ps[0]}
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a[common.B_VETO] == 0) & (a["score_dsmu"] > -1))
    if sig:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    else:
        m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HW      # OS sideband, 两侧同样处理
    return {k: v[m] for k, v in a.items()}


sig_score = np.concatenate([load(n, sig=True)["score_dsmu"] for n in common.SIGNALS])
THR = {e: float(np.quantile(sig_score, 1 - e)) for e in EFFS}

D = {}
for era in ["2022", "2022EE"]:
    D[era] = dict(data=load(f"Data_{era}"), wj=load([n.format(e=era) for n in WJ]))
print(f"W+jets source: {SOURCE}   signal (truth-matched): {len(sig_score):,} rows")


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
        for era in ["2022", "2022EE"]:
            s = D[era]["data" if kind == "data" else "wj"]
            w = None if wkey is None else s[wkey]
            e_[era], r_[era] = eff(s, mask_fn(s), w)
        out[kind] = (e_["2022EE"] / e_["2022"],
                     (e_["2022EE"] / e_["2022"]) * np.hypot(r_["2022EE"], r_["2022"]),
                     e_["2022"], e_["2022EE"])
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
    print(f"\n--- {f}  (保留 <= thr)")
    pre = D["2022"]["wj"]
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

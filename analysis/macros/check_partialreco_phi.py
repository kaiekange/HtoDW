"""partial-reco charm 在 lower sideband 里能不能解释 data 超出？(2026-09-28)

partial-reco = 三条径迹来自同一个 charm 强子但丢了粒子(pi0/gamma/多余的 pi)，
M(KKpi) 被压低到 M_parent - m_missing 以下(Ds->phi rho+ 端点 1.833，D+->phi pi pi0 端点 1.735)。
主导道都带一个**真 phi**，所以两个检验：

A. MC truth(GenTag)：按物理重分类，而不是按 cand_gen_mode(mode 6 只接全匹配的，π 没配上的
   会漏进 mode 8)。KK 同一母亲且是 phi 为「真 phi」；phi 的 nearest-HF 祖先存在为「HF phi」，
   再按 pi 分：pi 与 phi 同一个 HF 且是 HF 直接女儿 = 全重建；同一 HF 但不是直接女儿 = partial-reco
   (pi 来自 rho 等)；pi 不来自同一个 HF(含没配上) = HF phi + 别的 pi。
   **更正**：pi 是 HF 直接女儿不代表没丢粒子(D+->phi pi+ pi0 的 pi+ 也是直接女儿)，
   所以只分「phi 与 pi 同一个 HF」和「不同/没配上」，完整与否交给原始质量窗口判断。
   pruning 下 prompt(非 HF) phi 的 K 配不上，所以 truth 的「真 phi」只含 HF phi，
   比拟合出的 phi 峰少；拟合 vs truth 的差就是 truth 看不见的 phi。

B. data 侧：每个原始质量窗口 x score cut，对 data 和 MC 分别拟合 M(KK)，拆成 phi 峰 + 非共振，
   分别算 data/MC。lower sideband 的缺口若集中在 phi 峰部分 -> 缺的是带真 phi 的东西。
   MC 的 phi 峰拟合产额与 truth 真 phi 产额对照，作为拟合的 closure。

C. 原始 Ds_mass 上 data - MC 残差，和 MC partial-reco 成分的形状放在一起看有没有端点。

    conda activate ds_dnn && python3 check_partialreco_phi.py
"""
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import least_squares
from scipy.special import voigt_profile

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common  # noqa: E402

MODEL_KEY = "ssdata"
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, SR_HW = 1.96835, 0.040
M_PHI, G_PHI, M_K = 1.019461, 0.004249, 0.493677
KK_BINS = np.arange(2 * M_K, M_PHI + 0.0500001, 0.002)
KK_C = 0.5 * (KK_BINS[1:] + KK_BINS[:-1])
DS_BINS = np.arange(1.668, 2.2701, 0.020)
DS_C = 0.5 * (DS_BINS[1:] + DS_BINS[:-1])

SEL = [*common.SELECTION]
GEN = ["cand_gen_KK_mom_pdgId", "Kp_gen_hf_key", "Km_gen_hf_key", "pi_gen_hf_key",
       "Kp_gen_hf_pdgId", "pi_gen_hf_isMom"]


def load(name, extra=()):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        SEL + ["Ds_mass", "phi_mass", "score_dsmu"] + list(extra), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HW   # OS data blind; applied to MC too
    print(f"  {name:22} {m.sum():>9,} rows")
    return {k: v[m] for k, v in a.items()}


print("Loading...")
d = load("Data_OS_merged")
w = load("WtoLNu2Jets_merged", ["weight"] + GEN)
d["weight"] = np.ones(len(d["Ds_mass"]))

# ---- A. truth categories (MC) ----
realphi = w["cand_gen_KK_mom_pdgId"] == 333
hfphi = realphi & (w["Kp_gen_hf_key"] >= 0) & (w["Kp_gen_hf_key"] == w["Km_gen_hf_key"])
pi_same = hfphi & (w["pi_gen_hf_key"] == w["Kp_gen_hf_key"])
CATS = {
    "HF phi, pi same HF":    pi_same,
    "HF phi + other pi":     hfphi & ~pi_same,
    "no truth phi":          ~hfphi,
}
hfid = np.abs(w["Kp_gen_hf_pdgId"])

WIN = [("low <1.830", lambda a: a["Ds_mass"] < 1.830),
       ("D+- 1.83-1.91", lambda a: (a["Ds_mass"] >= 1.830) & (a["Ds_mass"] < 1.910)),
       ("high >2.008", lambda a: a["Ds_mass"] > M_DS + SR_HW)]
COLS = [("no cut", None), ("eps90", 0.90), ("eps70", 0.70), ("eps50", 0.50)]
THR = {lbl: None if e is None else common.stage1_threshold(
    dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=e))[0] for lbl, e in COLS}


def cutmask(a, lbl):
    return np.ones(len(a["Ds_mass"]), bool) if THR[lbl] is None else a["score_dsmu"] > THR[lbl]


print("\n=== A. MC truth composition (weighted %), per original-mass window ===")
for wn, wf in WIN:
    for lbl, _ in COLS:
        k = cutmask(w, lbl) & wf(w)
        tot = w["weight"][k].sum()
        parts = " | ".join(f"{c}: {100*w['weight'][k & m].sum()/tot:5.2f}" for c, m in CATS.items())
        print(f"{wn:14} {lbl:7} MCtot={tot:10.1f} | {parts}")
    print()
for wn, wf in WIN[:1]:
    for lbl, _ in COLS:
        k = cutmask(w, lbl) & wf(w) & CATS["HF phi, pi same HF"]
        u, idx = np.unique(hfid[k], return_inverse=True)
        s = np.bincount(idx, weights=w["weight"][k])
        o = np.argsort(-s)[:4]
        print(f"  phi+pi same-HF parent |pdgId| in {wn}, {lbl}: "
              + ", ".join(f"{u[i]}:{100*s[i]/s.sum():.0f}%" for i in o))


# ---- B. M(KK) fit: phi peak + threshold background ----
def model(p, x):
    ns, nb, sig, a, c = p
    q = np.clip(x - 2 * M_K, 0, None)
    sh_s = voigt_profile(x - M_PHI, sig, G_PHI / 2)
    sh_b = q ** a * (1 + c * q)
    return ns * sh_s / sh_s.sum() + nb * sh_b / sh_b.sum()


def fit_kk(x, wt):
    h = np.histogram(x, bins=KK_BINS, weights=wt)[0]
    h2 = np.histogram(x, bins=KK_BINS, weights=wt ** 2)[0]
    err = np.sqrt(np.maximum(h2, np.max(wt) ** 2 if len(wt) else 1.0))
    n = h.sum()
    p0 = [0.1 * n, 0.9 * n, 0.0015, 0.5, 0.0]
    r = least_squares(lambda p: (model(p, KK_C) - h) / err, p0,
                      bounds=([0, 0, 0.0003, 0.0, -30], [2 * n + 10, 2 * n + 10, 0.006, 3.0, 300]))
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, 2 * r.cost / max(len(h) - 5, 1))
        e = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        e = np.full(5, np.nan)
    return r.x, e, h, err


print("\n=== B. M(KK) fit: phi-peak vs non-resonant, data vs MC ===")
print(f"{'window':14} {'cut':6} | {'data phi':>14} {'MC phi':>14} {'MCtruth phi':>11} {'r_phi':>12} | "
      f"{'data non':>14} {'MC non':>14} {'r_non':>12} | gap(data-MC) phi / non")
FITS = {}
for wn, wf in WIN:
    for lbl, _ in COLS:
        kd, kw = cutmask(d, lbl) & wf(d), cutmask(w, lbl) & wf(w)
        pd, ed, hd, errd = fit_kk(d["phi_mass"][kd], d["weight"][kd])
        pw, ew, hw, errw = fit_kk(w["phi_mass"][kw], w["weight"][kw])
        FITS[(wn, lbl)] = (pd, hd, errd, pw, hw, errw)
        tru = w["weight"][kw & realphi].sum()
        rp = pd[0] / pw[0]
        rpe = rp * np.hypot(ed[0] / pd[0], ew[0] / pw[0])
        rn = pd[1] / pw[1]
        rne = rn * np.hypot(ed[1] / pd[1], ew[1] / pw[1])
        print(f"{wn:14} {lbl:6} | {pd[0]:8.0f}+-{ed[0]:5.0f} {pw[0]:8.0f}+-{ew[0]:5.0f} {tru:11.0f} "
              f"{rp:6.3f}+-{rpe:5.3f} | {pd[1]:8.0f}+-{ed[1]:5.0f} {pw[1]:8.0f}+-{ew[1]:5.0f} "
              f"{rn:6.3f}+-{rne:5.3f} | {pd[0]-pw[0]:8.0f} / {pd[1]-pw[1]:8.0f}")
    print()

# ---- figures ----
fig, axes = plt.subplots(2, len(COLS), figsize=(6.2 * len(COLS), 9), sharex=True)
for j, (lbl, _) in enumerate(COLS):
    pd, hd, errd, pw, hw, errw = FITS[("low <1.830", lbl)]
    for i, (p, h, e, name, col) in enumerate([(pd, hd, errd, "OS data", "black"),
                                               (pw, hw, errw, "W+jets MC", "#5a9e78")]):
        ax = axes[i][j]
        ax.errorbar(KK_C, h, yerr=e, fmt="o", ms=3, color=col, label=f"{name}, M(KK$\\pi$)<1.830")
        ax.plot(KK_C, model(p, KK_C), color="red", lw=1.5,
                label=f"fit: N$_\\phi$={p[0]:.0f}, N$_{{non}}$={p[1]:.0f}")
        ax.plot(KK_C, model([0, *p[1:]], KK_C), color="red", ls="--", lw=1.0)
        ax.set_ylim(0, 1.45 * h.max() if h.max() > 0 else 1)
        ax.legend(fontsize=9, loc="upper right")
        if i == 0:
            ax.set_title(lbl, fontsize=13)
        if j == 0:
            ax.set_ylabel("Events / 2 MeV")
    hep.cms.label(ax=axes[0][0], data=True, text="Preliminary",
                  rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)
fig.supxlabel(r"M(K$^+$K$^-$) [GeV]", fontsize=15, y=0.02)
out = common.fig_path("closure", "partialreco_phi_fit_lowsideband.png", key=MODEL_KEY)
fig.savefig(out, dpi=140, bbox_inches="tight")
print("Saved", out)

# C. data - MC residual in original mass, with the MC partial-reco component shape
fig, axes = plt.subplots(1, len(COLS), figsize=(6.2 * len(COLS), 5))
for j, (lbl, _) in enumerate(COLS):
    ax = axes[j]
    kd, kw = cutmask(d, lbl), cutmask(w, lbl)
    hd = np.histogram(d["Ds_mass"][kd], bins=DS_BINS)[0].astype(float)
    hw = np.histogram(w["Ds_mass"][kw], bins=DS_BINS, weights=w["weight"][kw])[0]
    hw2 = np.histogram(w["Ds_mass"][kw], bins=DS_BINS, weights=w["weight"][kw] ** 2)[0]
    blind = np.abs(DS_C - M_DS) < SR_HW
    res = np.where(blind, np.nan, hd - hw)
    ax.errorbar(DS_C, res, yerr=np.sqrt(hd + hw2), fmt="o", ms=4, color="black",
                label="OS data $-$ W+jets MC")
    for c, col in (("HF phi, pi same HF", "#8e44ad"), ("HF phi + other pi", "#5b6fc4")):
        k = kw & CATS[c]
        hp = np.histogram(w["Ds_mass"][k], bins=DS_BINS, weights=w["weight"][k])[0]
        ax.step(DS_BINS[:-1], hp, where="post", color=col, lw=1.6, label=f"MC {c} (as is)")
    ax.axhline(0, color="red", lw=1)
    ax.axvspan(M_DS - SR_HW, M_DS + SR_HW, color="0.85", zorder=0)
    for e in (1.833, 1.735):
        ax.axvline(e, color="navy", ls=":", lw=1.3)
    ax.set_title(lbl, fontsize=13)
    ax.legend(fontsize=8, loc="upper right")
    if j == 0:
        ax.set_ylabel("Events / 20 MeV")
        hep.cms.label(ax=ax, data=True, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)
fig.supxlabel(r"M(K$^+$K$^-\pi^\pm$) [GeV]", fontsize=15, y=0.0)
out = common.fig_path("closure", "partialreco_residual_dsmass.png", key=MODEL_KEY)
fig.savefig(out, dpi=140, bbox_inches="tight")
print("Saved", out)

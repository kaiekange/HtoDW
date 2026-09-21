"""KpToPi 互换后的 M(KK pi): OS data vs W+jets GenTag，按 stage-1 score cut 铺开。

上侧 sideband 里 W+jets 预测得比 data 多（2026-09-21 实测 data/MC = 0.86/0.89/0.81
at eps_sig = 90/70/50%），而下侧反过来。这张图针对上侧：那里的本底主要是真
D+- -> K-+ pi+- pi+- 里有一个 pi 被误指认成 K 才进了 KKpi 候选池，把它换回 pi 质量
就该复原出 D+- 峰（1.8696），且该峰应当几乎全部来自 GenTag mode 5。

互换配方与 check_misid_veto.py:92-97 逐字相同 —— Kp 取 PION 质量，Km 保持 KAON，
pi 保持 PION，用原始动量分量（不是 vertex-fit refit 的）：
    E = sum sqrt(px^2+py^2+pz^2+m^2),  m = (M_PI, M_K, M_PI) for (Kp, Km, pi)
为什么是 Kp：producer 里 K+ 与 pi 都与 muon 反号、K- 与 muon 同号，所以在真实的
D+ -> K- pi+ pi+ 中，Km 是真 K-，而 Kp 与 pi 是那两个真 pi+ —— 被误指认的是 Kp。

必须做对的一点（旧 draw_misid_data_vs_wjets.py 的 docstring 记过）：
**SR 要在原始 Ds_mass 上挖，并且对两个样本同样施加。** OS data 的 SR 是 blind 的而
W+jets 不是；blind 在原始质量上是一个干净窗口，映射到互换质量上却被抹开成弥散的一片，
所以不挖就等于 data 在互换谱的一大片区域里系统性缺事例而 MC 不缺 —— 那会被读成一个
假的形状差异。

绝对产额，不归一化：W+jets 带 sigma*L，它的高度本身就是被检验的量。

    conda activate ds_dnn && python3 draw_misid_mass_data_vs_wjets.py [--split8]
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
# Luminosity from eras.py, never hardcoded -- a literal here goes stale silently the
# moment an era is added, and nothing would flag it.
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
PRED = f"{common.out_dir_for(MODEL_KEY)}/Prediction"
M_DS, SR_HW = 1.96835, 0.040
M_D = 1.8696
M_K, M_PI = 0.493677, 0.13957039
# [1.30, 2.30] covers 97.6% of the swapped-mass spectrum. An earlier [1.60, 2.20]
# covered only 55.9% -- the median is 1.631, so it cut away the continuum BELOW the D
# peak and made the peak look far more prominent than it is.
BINS = np.arange(1.30, 2.301, 0.010)
CENT = 0.5 * (BINS[1:] + BINS[:-1])

# --split8 pulls GenTag mode 8 out of "other" into its own stack entry.
#
# What mode 8 IS (HDsPhiPiMuGenTagProducer.cc:923-936, read rather than assumed):
# the cascade reaches mode 8 as `!one_origin` AFTER mode 0 (n_hf == 0) and mode 7
# (n_distinct >= 2) have been taken, and one_origin = (n_hf == 3 && n_distinct == 1).
# So mode 8 is exactly: 1 or 2 of the three tracks come from a heavy-flavour hadron
# and they all share the SAME one, while the remaining track(s) do not come from
# heavy flavour at all. I.e. a REAL charm hadron plus one or two unrelated tracks --
# partially-reconstructed charm, not pure combinatorial (that is mode 0) and not a
# two-hadron mixture (that is mode 7). It is the largest single component hiding
# inside "other" after the selection.
SPLIT8 = "--split8" in sys.argv

_OTHER = [0, 6, 7, 9, 10] if SPLIT8 else [0, 6, 7, 8, 9, 10]
GROUPS = [   # stack bottom -> top, same order/colours as draw_stage1_mass_by_score_cut.py
    (rf"other (mode {','.join(map(str, _OTHER))})", _OTHER, "#5a9e78"),
]
if SPLIT8:
    # Blue-violet: has to be unmistakable against the green/yellow/orange/red already
    # in the stack, since the whole point of this variant is reading mode 8's shape.
    GROUPS.append((r"partial HF, 1-2 trk (mode 8)", [8], "#5b6fc4"))
GROUPS += [
    (r"D$\to$K$\pi\pi$ misID (5)",   [5],    "#e6c229"),
    (r"genuine D$\to$KK$\pi$ (3+4)", [3, 4], "#e07a3c"),
    (r"true D$_s$ (mode 1+2)",       [1, 2], "#c0392b"),
]
SEL = [*common.SELECTION]
MOM = [f"{p}_p{x}" for p in ("Kp", "Km", "pi") for x in ("x", "y", "z")]


def misid_mass(a):
    """M(KK pi) with the Kp track re-assigned the pion mass."""
    e = sum(np.sqrt(a[f"{p}_px"]**2 + a[f"{p}_py"]**2 + a[f"{p}_pz"]**2 + m**2)
            for p, m in (("Kp", M_PI), ("Km", M_K), ("pi", M_PI)))
    px, py, pz = (sum(a[f"{p}_p{x}"] for p in ("Kp", "Km", "pi")) for x in ("x", "y", "z"))
    return np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))


def load(name, extra=()):
    a = uproot.open(f"{PRED}/{name}.root")["Events"].arrays(
        SEL + MOM + ["Ds_mass", "score_dsmu"] + list(extra), library="np")
    m = common.selection_mask(a) & (a["score_dsmu"] > -1)
    # SR removed on the ORIGINAL mass, for BOTH samples -- see module docstring.
    m &= np.abs(a["Ds_mass"] - M_DS) >= SR_HW
    out = {k: v[m] for k, v in a.items()}
    out["mis"] = misid_mass(out)
    print(f"  {name:22} {m.sum():>9,} rows")
    return out


print("Loading (current selection + score>-1 + Ds-SR removed on ORIGINAL mass)...")
d = load("Data_OS_merged")
w = load("WtoLNu2Jets_merged", ["weight", "cand_gen_mode"])

COLS = [("no cut", None), (r"$\varepsilon_{sig}$=90%", 0.90),
        (r"$\varepsilon_{sig}$=70%", 0.70), (r"$\varepsilon_{sig}$=50%", 0.50)]
THR = {}
for lbl, eff in COLS:
    THR[lbl] = None if eff is None else common.stage1_threshold(
        dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=eff))[0]

fig, axes = plt.subplots(2, len(COLS), figsize=(6.2 * len(COLS), 7.4),
                         gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06, wspace=0.22))
for j, (lbl, _eff) in enumerate(COLS):
    thr = THR[lbl]
    ax, rax = axes[0][j], axes[1][j]
    sd = np.ones(len(d["mis"]), bool) if thr is None else d["score_dsmu"] > thr
    sw = np.ones(len(w["mis"]), bool) if thr is None else w["score_dsmu"] > thr
    width = BINS[1] - BINS[0]

    hs, errs = [], []
    for _n, modes, _c in GROUPS:
        k = sw & np.isin(w["cand_gen_mode"], modes)
        h, _ = np.histogram(w["mis"][k], bins=BINS, weights=w["weight"][k])
        h2, _ = np.histogram(w["mis"][k], bins=BINS, weights=w["weight"][k] ** 2)
        hs.append(h)
        errs.append(np.sqrt(h2))
    tot = np.sum(hs, axis=0)
    tot_err = np.sqrt(np.sum(np.square(errs), axis=0))

    bottom = np.zeros(len(CENT))
    for (lab, _m, cl), h in zip(GROUPS, hs):
        ax.bar(CENT, h, width=width, bottom=bottom, color=cl, alpha=0.9,
               label=f"{lab}  N={h.sum():,.0f}", edgecolor="none")
        bottom = bottom + h
    ax.bar(CENT, 2 * tot_err, width=width, bottom=tot - tot_err, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3, label="W+jets MC stat")

    hd, _ = np.histogram(d["mis"][sd], bins=BINS)
    n_drop = sd.sum() - hd.sum()
    print(f"{lbl:>22}: out-of-range OS data {n_drop:,} of {sd.sum():,} "
          f"({100*n_drop/max(sd.sum(),1):.1f}%)")
    ax.errorbar(CENT, hd, yerr=np.sqrt(hd), fmt="o", ms=4, color="black", elinewidth=1.2,
                label=f"OS data (Ds-SR excl.)  N={hd.sum():,.0f}")
    ax.axvline(M_D, color="navy", ls=":", lw=1.4)

    ax.set_xlim(BINS[0], BINS[-1])
    peak = max(bottom.max(), hd.max())
    ax.set_ylim(0, peak * 1.5 if peak > 0 else 1.0)
    ax.tick_params(labelbottom=False)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(lbl, fontsize=13, pad=6)
    if j == 0:
        ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
        # text=, NOT label=: mplhep deprecated `label` and only emits a FutureWarning,
        # so `label="Preliminary"` renders a bare "CMS" and the sub-label is silently
        # dropped -- same class of bug as histplot(histtype="fill") swallowing yerr.
        hep.cms.label(ax=ax, data=True, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

    # Ratio: data points carry the DATA stat only; MC stat is the band around 1, so a
    # deviation can be read as data- or MC-limited. Reference: draw_stage1_mass_by_score_cut.py.
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(tot > 0, hd / tot, np.nan)
        re = np.where(tot > 0, np.sqrt(hd) / tot, np.nan)
        mc_band = np.where(tot > 0, tot_err / tot, np.nan)
    rax.bar(CENT, 2 * mc_band, width=width, bottom=1 - mc_band, color="none",
            hatch="////", edgecolor="0.35", linewidth=0.0)
    rax.errorbar(CENT, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
    rax.axhline(1.0, color="red", lw=1.0)
    rax.axvline(M_D, color="navy", ls=":", lw=1.4)
    rax.set_ylim(0.0, 1.99)
    rax.set_yticks([0.5, 1.0, 1.5])
    rax.set_xlim(BINS[0], BINS[-1])
    if j == 0:
        rax.set_ylabel("data / MC")

    pk = (CENT > M_D - 0.02) & (CENT < M_D + 0.02)
    per = [f"{n.split('(')[0].strip()}={h[pk].sum():.0f}" for (n, _m, _c), h in zip(GROUPS, hs)]
    print(f"{'':>22}  D-peak +-20MeV: data={hd[pk].sum():>7,.0f} MC={tot[pk].sum():>9,.0f} "
          f"ratio={hd[pk].sum()/max(tot[pk].sum(),1e-9):.3f} | " + ", ".join(per))

fig.supxlabel(r"M(K$^+$K$^-\pi^\pm$) with the K$^+$ track re-assigned the $\pi$ mass  [GeV]",
              fontsize=15, y=0.02)
fig.suptitle(r"misID reflection: K$^+\to\pi^+$ swap -- OS data vs W+jets GenTag", fontsize=15)
out = common.fig_path("misid",
                      "misid_KpToPi_mass_data_vs_wjets"
                      + ("_split8" if SPLIT8 else "") + ".png", key=MODEL_KEY)
fig.savefig(out, dpi=140, bbox_inches="tight")
print("Saved", out)

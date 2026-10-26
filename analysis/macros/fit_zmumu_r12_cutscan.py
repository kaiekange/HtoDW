# Quick (scipy chi2, not RooFit) r12 vs stage-1 cut, 2026-10-26. Run with the ds_dnn python from analysis/macros/.
# Postfit: zmumu_r12_fit_cut_<cut>_pol1.png; summary: zmumu_r12_vs_cut.png (under ssdata_f45/figures/sf/).
# Z->mumu CR (76<M(mumu)<106, phi20): r12 vs stage-1 ssdata_f45 cut (score_dsmu > thr(eps_sig)).
# Same OS+SS joint chi2 fit as fit_zmumu_r12_quick.py: per charge n_c * MC_cont * (1 + a1 x) + r12 * Ds + r34 * D,
# range 1.75-2.20, peak shapes fixed to MC. Reads Prediction/zmumu_*.root (rows = skim rows + score).
import sys
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import numpy as np, uproot, common
from scipy.optimize import least_squares
KEY = "ssdata_f45"
P = f"{common.out_dir_for(KEY)}/Prediction"
M_DS, M_PHI = 1.96835, 1.019460
LO, HI = 1.75, 2.20
EFFS = [None, 0.99, 0.98, 0.97, 0.95, 0.93, 0.90]   # tighter than eps90 not considered (user 1026)
thr = {e: common.stage1_threshold(dict(key=KEY, tool="DNN", seed=42, eff_sig=e))[0] for e in EFFS if e}
cols = ["Ds_mass", "phi_mass", "m_mumu", "weight", "score_dsmu", "cand_gen_mode", "cand_gen_common_hf_pdgId"]
def get(f, mc, dy=False):
    t = uproot.open(f"{P}/zmumu_{f}.root")["Events"]
    a = t.arrays([c for c in cols + ["jet_bin"] if c in t.keys()], library="np")
    m = ((a["m_mumu"] > 76) & (a["m_mumu"] < 106) & (np.abs(a["phi_mass"] - M_PHI) < 0.020)
         & (a["Ds_mass"] > LO) & (a["Ds_mass"] < HI) & (a["score_dsmu"] > -1))
    a = {k: v[m] for k, v in a.items()}
    w = np.ones(m.sum())
    if mc:
        w = a["weight"].astype(float) * common.br_reweight(a)
        if dy: w = w * common.dy_jetbin_sf(a)
    return a, w
RAW = {}
for ch, sfx in (("OS", ""), ("SS", "_SS")):
    RAW[ch] = {"data": get("Data_OS_merged" if ch == "OS" else "Data_SS_merged", False),
               "mc": [get(f + sfx + "_merged", True, dy) for f, dy in (("DY", True), ("TT", False), ("WtoLNu2Jets", False))]}

def hists(eff, bw):
    bins = np.arange(LO, HI + 1e-9, bw); bc = 0.5 * (bins[1:] + bins[:-1])
    H = {}
    for ch in ("OS", "SS"):
        a, _ = RAW[ch]["data"]
        k = np.ones(len(a["Ds_mass"]), bool) if eff is None else a["score_dsmu"] > thr[eff]
        d = np.histogram(a["Ds_mass"][k], bins)[0].astype(float)
        T = {x: np.zeros(len(bc)) for x in ("cont", "Ds", "D", "var")}
        for a, w in RAW[ch]["mc"]:
            k = np.ones(len(w), bool) if eff is None else a["score_dsmu"] > thr[eff]
            md = a["cand_gen_mode"]
            for x, sel in (("cont", ~np.isin(md, [1, 2, 3, 4])), ("Ds", np.isin(md, [1, 2])), ("D", np.isin(md, [3, 4]))):
                T[x] += np.histogram(a["Ds_mass"][k & sel], bins, weights=w[k & sel])[0]
            T["var"] += np.histogram(a["Ds_mass"][k], bins, weights=w[k] ** 2)[0]
        H[ch] = (d, T)
    return bc, H

def fit(bc, H, deg):
    x = bc - M_DS
    def res(p):
        rr = []
        for i, ch in enumerate(("OS", "SS")):
            d, T = H[ch]
            corr = p[2 + i] * (1 + (p[4 + i] * x if deg == 1 else 0))
            mu = corr * T["cont"] + p[0] * T["Ds"] + p[1] * T["D"]
            rr.append((d - mu) / np.sqrt(np.maximum(d + T["var"], 1.0)))
        return np.concatenate(rr)
    p0 = np.r_[1, 1, 1, 1, np.zeros(2 * deg)]
    s = least_squares(res, p0); cov = np.linalg.inv(s.jac.T @ s.jac)
    fit.last = s.x
    return s.x[0], np.sqrt(cov[0, 0]), np.sum(s.fun ** 2), len(s.fun) - len(p0)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
LUMI = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
SUMMARY = []


def draw(lab, bc, bw, H, p, r12, e12, chi2, ndf):
    """Postfit (pol1): stack + data/postfit + data minus fitted continuum, OS and SS."""
    edges = np.r_[bc - bw / 2, bc[-1] + bw / 2]
    x = bc - M_DS
    fig, axes = plt.subplots(3, 2, figsize=(18, 16), gridspec_kw=dict(height_ratios=[3, 1, 2], hspace=0.08), sharex=True)
    for j, ch in enumerate(("OS", "SS")):
        d, T = H[ch]
        cont = p[2 + j] * (1 + p[4 + j] * x) * T["cont"]
        Dp, Dsp = p[1] * T["D"], p[0] * T["Ds"]
        tot = cont + Dp + Dsp
        sig = np.sqrt(np.maximum(d + T["var"], 1.0))
        ax = axes[0, j]
        hep.histplot([cont, Dp, Dsp], edges, stack=True, histtype="fill", ax=ax,
                     color=["#9e9e9e", "#e07a3c", "#c0392b"],
                     label=[f"continuum (MC shape x pol1, fitted)  {cont.sum():,.0f}",
                            f"true D$^\\pm$ x r34={p[1]:.2f}  {Dp.sum():.0f}",
                            f"true D$_s$ x r12={p[0]:.2f}  {Dsp.sum():.0f}"])
        hep.histplot(T["cont"] + T["D"] + T["Ds"], edges, ax=ax, color="blue", lw=1.5, ls="--", label="prefit MC (r=1)")
        ax.errorbar(bc, d, yerr=np.sqrt(d), fmt="o", ms=4, color="k", label=f"data  {d.sum():,.0f}")
        ax.set_ylim(0, 1.5 * max(d.max(), tot.max())); ax.set_ylabel(f"Events / {1000*bw:.0f} MeV")
        ax.legend(fontsize=11, loc="upper right")
        ax.text(0.04, 0.90, f"Z$\\to\\mu\\mu$ CR, {ch}", transform=ax.transAxes, fontsize=18)
        hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI:.1f}", com=13.6, fontsize=14)
        rax = axes[1, j]
        with np.errstate(invalid="ignore", divide="ignore"):
            ok = tot > 0
            rax.errorbar(bc[ok], d[ok] / tot[ok], yerr=sig[ok] / tot[ok], fmt="o", ms=3.5, color="k")
        rax.axhline(1, color="r", lw=1); rax.set_ylim(0.5, 1.5); rax.set_ylabel("data/postfit", fontsize=13)
        sax = axes[2, j]
        sax.errorbar(bc, d - cont, yerr=sig, fmt="o", ms=4, color="k", label="data $-$ fitted continuum")
        hep.histplot(Dp + Dsp, edges, ax=sax, color="#c0392b", lw=2, label=f"fitted peaks (r12={p[0]:.2f}, r34={p[1]:.2f})")
        hep.histplot(T["D"] + T["Ds"], edges, ax=sax, color="blue", lw=1.5, ls="--", label="MC peaks (r=1)")
        sax.axhline(0, color="0.5", lw=0.8)
        sax.set_ylabel("Events / bin"); sax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        sax.legend(fontsize=11, loc="upper left"); sax.set_xlim(LO, HI)
    cut_txt = "no DNN cut" if lab == "no cut" else f"stage-1 {KEY} {lab}"
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR OS+SS joint fit, 76<M($\\mu\\mu$)<106, $\\phi$20, {cut_txt}, range {LO}-{HI}, pol1 per charge:  "
                 f"r12 = {r12:.2f} $\\pm$ {e12:.2f},  $\\chi^2$/ndf = {chi2:.0f}/{ndf}\n"
                 "(MC: skim weight x BR rw x DY jet-bin SF; peak shapes fixed to MC; errors data + MC stat)", fontsize=15)
    tag = "nocut" if lab == "no cut" else lab
    out = common.fig_path("sf", f"zmumu_r12_fit_cut_{tag}_pol1.png", key=KEY)
    plt.savefig(out, dpi=90, bbox_inches="tight"); plt.close(fig)
    print("   Saved", out)


print(f"{'cut':>7s} {'thr':>8s} {'data':>7s} {'MC trueDs':>9s} {'Ds share':>8s} {'bin':>4s} | "
      f"{'r12 pol0':>15s} {'chi2/ndf':>9s} | {'r12 pol1':>15s} {'chi2/ndf':>9s}")
for eff in EFFS:
    bc0, H0 = hists(eff, 0.010)
    ndata = sum(H0[c][0].sum() for c in H0)
    bw = 0.010 if ndata > 20000 else 0.020
    bc, H = hists(eff, bw)
    win = np.abs(bc - M_DS) < 0.040
    tds = sum(H[c][1]["Ds"][win].sum() for c in H)
    tot = sum((H[c][1]["cont"] + H[c][1]["Ds"] + H[c][1]["D"])[win].sum() for c in H)
    r0, e0, c0, n0 = fit(bc, H, 0)
    r1, e1, c1, n1 = fit(bc, H, 1)
    lab = "no cut" if eff is None else f"eps{round(100*eff)}"
    print(f"{lab:>7s} {thr.get(eff, -1):8.4f} {ndata:7.0f} {sum(H[c][1]['Ds'].sum() for c in H):9.1f} {tds/tot:8.3f} {1000*bw:3.0f}M | "
          f"{r0:6.3f} +- {e0:5.3f} {c0:5.1f}/{n0:<3d} | {r1:6.3f} +- {e1:5.3f} {c1:5.1f}/{n1:<3d}")
    fit(bc, H, 1)
    SUMMARY.append((lab, eff, r0, e0, r1, e1, tds / tot))
    draw(lab, bc, bw, H, fit.last, r1, e1, c1, n1)

# summary: r12 vs cut (points are NESTED subsets -> errors strongly correlated between points)
fig, ax = plt.subplots(figsize=(11, 8))
xs = np.arange(len(SUMMARY))
ax.errorbar(xs - 0.08, [r[2] for r in SUMMARY], yerr=[r[3] for r in SUMMARY], fmt="o", color="#4a90c2", label="continuum: MC shape (pol0)")
ax.errorbar(xs + 0.08, [r[4] for r in SUMMARY], yerr=[r[5] for r in SUMMARY], fmt="s", color="#c0392b", label="continuum: MC shape x pol1")
for i, r in enumerate(SUMMARY):
    ax.text(xs[i], 0.42, f"{100*r[6]:.0f}%", ha="center", fontsize=11, color="0.3")
ax.text(0.99, 0.97, "grey numbers: true D/Ds share of MC\nin the Ds window (+-40 MeV)", transform=ax.transAxes, fontsize=11, color="0.3", ha="right", va="top")
ax.axhline(1, color="k", lw=1, ls=":")
ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in SUMMARY], fontsize=14)
ax.set_ylim(0.35, 1.7); ax.set_ylabel("r12 (data / MC true Ds)"); ax.set_xlabel(f"stage-1 {KEY} cut")
ax.legend(fontsize=13, loc="upper left")
fig.suptitle("Z$\\to\\mu\\mu$ CR, OS+SS joint fit, 1.75-2.20 GeV (nested cuts: errors correlated)", fontsize=14, y=1.0)
hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI:.1f}", com=13.6, fontsize=14)
out = common.fig_path("sf", "zmumu_r12_vs_cut.png", key=KEY)
plt.savefig(out, dpi=90, bbox_inches="tight"); plt.close(fig)
print("Saved", out)

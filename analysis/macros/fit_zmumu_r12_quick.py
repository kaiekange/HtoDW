# Quick (scipy chi2, not RooFit) r12 fit, 2026-10-26. Run with the ds_dnn python from analysis/macros/.
# Postfit figures: <key ssdata_f45>/figures/sf/zmumu_r12_fit_<range>_pol<deg>.png
# Z->mumu CR (76<M(mumu)<106, phi20, no DNN cut): OS+SS joint binned chi2 fit of Ds_mass.
# Per charge c: mu_c = n_c * cont_c * (1 + a1 x + a2 x^2) + r12 * Ds_c + r34 * D_c,  x = m - M_Ds
#   cont = MC modes 0,5,6-10 (DY+TT+W x BR rw x DY SF), Ds = modes 1+2, D = modes 3+4.
# r12, r34 shared by OS and SS. sigma^2 = data + MC variance of the template sum (fixed, at r=1).
import sys
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import numpy as np, uproot, common
from scipy.optimize import least_squares
D_ = f"{common.SAMPLES_DIR}/CR/zmumu"
M_DS, M_PHI = 1.96835, 1.019460
BW = 0.010
bins = np.arange(1.67, 2.27 + 1e-9, BW); bc = 0.5 * (bins[1:] + bins[:-1])
cols = ["Ds_mass", "phi_mass", "m_mumu", "weight", "cand_gen_mode", "cand_gen_common_hf_pdgId"]
def get(f, mc, dy=False):
    t = uproot.open(f"{D_}/{f}.root")["Events"]
    a = t.arrays([c for c in cols + ["jet_bin"] if c in t.keys()], library="np")
    m = (a["m_mumu"] > 76) & (a["m_mumu"] < 106) & (np.abs(a["phi_mass"] - M_PHI) < 0.020)
    a = {k: v[m] for k, v in a.items()}
    w = a["weight"].astype(float)
    if mc:
        w = w * common.br_reweight(a)
        if dy: w = w * common.dy_jetbin_sf(a)
    return a, w
H = {}
for ch, sfx in (("OS", ""), ("SS", "_SS")):
    a, _ = get("Data_OS_merged" if ch == "OS" else "Data_SS_merged", False)
    d = np.histogram(a["Ds_mass"], bins)[0].astype(float)
    T = {k: np.zeros(len(bc)) for k in ("cont", "Ds", "D", "var")}
    for f, dy in (("DY", True), ("TT", False), ("WtoLNu2Jets", False)):
        a, w = get(f + sfx + "_merged", True, dy)
        md = a["cand_gen_mode"]
        for k, sel in (("cont", ~np.isin(md, [1, 2, 3, 4])), ("Ds", np.isin(md, [1, 2])), ("D", np.isin(md, [3, 4]))):
            T[k] += np.histogram(a["Ds_mass"][sel], bins, weights=w[sel])[0]
        T["var"] += np.histogram(a["Ds_mass"], bins, weights=w * w)[0]
    H[ch] = (d, T)
    print(f"{ch}: data {d.sum():.0f}  MC cont {T['cont'].sum():.0f}  Ds {T['Ds'].sum():.1f}  D {T['D'].sum():.1f}")

def fit(deg, shared_shape, lo, hi):
    r = (bc > lo) & (bc < hi); x = bc[r] - M_DS
    chs = ["OS", "SS"]
    npoly = deg if shared_shape else 2 * deg
    def unpack(p):
        r12, r34, nOS, nSS = p[:4]; poly = p[4:]
        out = {}
        for i, ch in enumerate(chs):
            co = poly[:deg] if shared_shape else poly[i * deg:(i + 1) * deg]
            corr = 1 + sum(co[j] * x ** (j + 1) for j in range(deg))
            out[ch] = (nOS if ch == "OS" else nSS) * corr
        return r12, r34, out
    def res(p):
        r12, r34, corr = unpack(p)
        rr = []
        for ch in chs:
            d, T = H[ch]
            mu = corr[ch] * T["cont"][r] + r12 * T["Ds"][r] + r34 * T["D"][r]
            rr.append((d[r] - mu) / np.sqrt(d[r] + T["var"][r]))
        return np.concatenate(rr)
    p0 = np.r_[1, 1, 1, 1, np.zeros(npoly)]
    s = least_squares(res, p0)
    J = s.jac; cov = np.linalg.inv(J.T @ J)
    chi2 = np.sum(s.fun ** 2); ndf = len(s.fun) - len(p0)
    fit.last = (s.x, unpack, r)
    return s.x[0], np.sqrt(cov[0, 0]), s.x[1], np.sqrt(cov[1, 1]), chi2, ndf, cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])

print(f"\nbin {1000*BW:.0f} MeV")
for lo, hi in ((1.67, 2.27), (1.75, 2.20)):
    for shared in (False, True):
        for deg in (0, 1, 2):
            r12, e12, r34, e34, chi2, ndf, rho = fit(deg, shared, lo, hi)
            print(f"range {lo}-{hi}  pol{deg} {'shared ' if shared else 'per-chg'}  "
                  f"r12 = {r12:.3f} +- {e12:.3f}   r34 = {r34:.3f} +- {e34:.3f}  rho {rho:+.2f}  chi2/ndf {chi2:.1f}/{ndf}")
# OS-only and SS-only (pol1 per charge, full range) to test that OS and SS agree
for only in ("OS", "SS"):
    keep = H.copy(); H = {only: keep[only]}
    r = (bc > 1.67) & (bc < 2.27); x = bc[r] - M_DS
    d, T = H[only]
    def res1(p):
        mu = p[2] * (1 + p[3] * x) * T["cont"][r] + p[0] * T["Ds"][r] + p[1] * T["D"][r]
        return (d[r] - mu) / np.sqrt(d[r] + T["var"][r])
    s = least_squares(res1, [1, 1, 1, 0]); cov = np.linalg.inv(s.jac.T @ s.jac)
    print(f"{only} only, pol1: r12 = {s.x[0]:.3f} +- {np.sqrt(cov[0,0]):.3f}   r34 = {s.x[1]:.3f} +- {np.sqrt(cov[1,1]):.3f}  chi2/ndf {np.sum(s.fun**2):.1f}/{len(s.fun)-4}")
    H = keep


# ---------------- postfit figures ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
LUMI = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
for lo, hi, deg in ((1.67, 2.27, 0), (1.67, 2.27, 1), (1.67, 2.27, 2), (1.75, 2.20, 0), (1.75, 2.20, 1), (1.75, 2.20, 2)):
    r12, e12, r34, e34, chi2, ndf, _ = fit(deg, False, lo, hi)
    px, unpack, rm = fit.last
    _r12, _r34, corr = unpack(px)
    edges = bins[np.r_[rm, False] | np.r_[False, rm]]
    xc = bc[rm]
    fig, axes = plt.subplots(3, 2, figsize=(18, 16), gridspec_kw=dict(height_ratios=[3, 1, 2], hspace=0.08), sharex=True)
    for j, ch in enumerate(("OS", "SS")):
        d, T = H[ch]
        dd, cont = d[rm], corr[ch] * T["cont"][rm]
        Dp, Dsp = r34 * T["D"][rm], r12 * T["Ds"][rm]
        tot = cont + Dp + Dsp
        ax = axes[0, j]
        hep.histplot([cont, Dp, Dsp], edges, stack=True, histtype="fill", ax=ax,
                     color=["#9e9e9e", "#e07a3c", "#c0392b"],
                     label=[f"continuum (MC shape" + ("" if deg == 0 else f" x pol{deg}") + f", fitted)  {cont.sum():,.0f}",
                            f"true D$^\pm$ x r34={r34:.2f}  {Dp.sum():.0f}",
                            f"true D$_s$ x r12={r12:.2f}  {Dsp.sum():.0f}"])
        prefit = T["cont"][rm] + T["D"][rm] + T["Ds"][rm]
        hep.histplot(prefit, edges, ax=ax, color="blue", lw=1.5, ls="--", label="prefit MC (r=1)")
        ax.errorbar(xc, dd, yerr=np.sqrt(dd), fmt="o", ms=4, color="k", label=f"data  {dd.sum():,.0f}")
        ax.set_ylim(0, 1.45 * dd.max()); ax.set_ylabel(f"Events / {1000*BW:.0f} MeV")
        ax.legend(fontsize=11, loc="upper right")
        ax.text(0.04, 0.90, f"Z$\\to\\mu\\mu$ CR, {ch}", transform=ax.transAxes, fontsize=18)
        hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI:.1f}", com=13.6, fontsize=14)
        rax = axes[1, j]
        sig = np.sqrt(dd + T["var"][rm])
        rax.errorbar(xc, dd / tot, yerr=sig / tot, fmt="o", ms=3.5, color="k")
        rax.axhline(1, color="r", lw=1); rax.set_ylim(0.9, 1.1); rax.set_ylabel("data/postfit", fontsize=13)
        sax = axes[2, j]
        sax.errorbar(xc, dd - cont, yerr=sig, fmt="o", ms=4, color="k", label="data $-$ fitted continuum")
        hep.histplot(Dp + Dsp, edges, ax=sax, color="#c0392b", lw=2, label=f"fitted peaks (r12={r12:.2f}, r34={r34:.2f})")
        hep.histplot(T["D"][rm] + T["Ds"][rm], edges, ax=sax, color="blue", lw=1.5, ls="--", label="MC peaks (r=1)")
        sax.axhline(0, color="0.5", lw=0.8)
        sax.set_ylabel("Events / bin"); sax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        sax.legend(fontsize=11, loc="upper left")
        sax.set_xlim(edges[0], edges[-1])
    fig.suptitle(f"Z$\\to\\mu\\mu$ CR OS+SS joint fit, 76<M($\\mu\\mu$)<106, $\\phi$20, no DNN cut, range {lo}-{hi}, pol{deg} per charge:  "
                 f"r12 = {r12:.2f} $\\pm$ {e12:.2f},  r34 = {r34:.2f} $\\pm$ {e34:.2f},  $\\chi^2$/ndf = {chi2:.0f}/{ndf}\n"
                 "(MC: skim weight x BR rw x DY jet-bin SF; peak shapes fixed to MC; errors data + MC stat)", fontsize=15)
    out = common.fig_path("sf", f"zmumu_r12_fit_{lo:.2f}-{hi:.2f}_pol{deg}.png", key="ssdata_f45")
    plt.savefig(out, dpi=90, bbox_inches="tight"); plt.close(fig)
    print("Saved", out)

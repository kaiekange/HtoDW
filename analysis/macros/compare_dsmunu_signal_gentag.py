#!/usr/bin/env python3
"""m(Ds mu nu) shape: signal vs W+jets GenTag modes 1-5, from the merged NanoAOD samples.

m(Ds mu nu) is built exactly as in check_dsmunu_mass.py (functions imported from there):
corrected muon (stored direction x mu_pt_corr), Ds = stored (Ds+mu) - uncorrected mu,
nu px/py from the corrected MET, nu pz from the W-mass constraint (complex -> real part,
two real -> smaller |pz|). Baseline for comparison: the current feature Dsmu_mass_corr.

Inputs: NanoMVA/Samples/{Signal,WtoLNu2Jets}_merged.root (selection already applied),
weights = `weight` (absolute for both). GenTag groups (cand_gen_mode): 1+2 real Ds -> KKpi,
3+4 D+ -> KKpi, 5 D+ -> Kpipi reflection. Shapes are normalised; negative W+jets weights
are kept signed.

Plots (left m(Ds mu), right m(Ds mu nu), both from corrected inputs) follow draw_wp_mass.py: absolute lumi-weighted W+jets stack by GenTag mode with MC
stat hatch; signal drawn as lines scaled by one common factor (stated in the legend); the
signal central 68% interval is shaded. One figure for the full Ds-mass window, one for SR.

Metric (m is not monotonic, so no AUC): the central 68% interval of the signal
distribution (weighted 16th-84th percentile) and the weighted fraction of each background
group that falls inside it -- for m(Ds mu nu) and for the baseline.

Run in conda ds_dnn (check_dsmunu_mass imports common, which needs torch).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
from check_dsmunu_mass import M_MU, four, inv_mass, nu_pz, common  # noqa: E402
from eras import M_DS, SR_HALF_WIDTH  # noqa: E402

plt.style.use(hep.style.CMS)

SAMPLES = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples"
COLS = ["mu_px", "mu_py", "mu_pz", "mu_pt", "mu_pt_corr", "Dsmu_px", "Dsmu_py", "Dsmu_pz",
        "Dsmu_mass", "Dsmu_mass_corr", "met_pt_corr", "met_phi_corr", "Ds_mass", "weight"]


def load(name, extra):
    t = uproot.open(f"{SAMPLES}/{name}.root")["Events"]
    a = t.arrays(COLS + extra, library="np")
    f = {k: a[k].astype(np.float64) for k in COLS}
    mu_raw = four(f["mu_px"], f["mu_py"], f["mu_pz"], M_MU)
    k = f["mu_pt_corr"] / f["mu_pt"]
    mu = four(k * f["mu_px"], k * f["mu_py"], k * f["mu_pz"], M_MU)
    ds = four(f["Dsmu_px"], f["Dsmu_py"], f["Dsmu_pz"], f["Dsmu_mass"]) - mu_raw
    nux = f["met_pt_corr"] * np.cos(f["met_phi_corr"])
    nuy = f["met_pt_corr"] * np.sin(f["met_phi_corr"])
    pz1, pz2, cplx = nu_pz(mu, nux, nuy)
    pz = np.where(np.abs(pz1) < np.abs(pz2), pz1, pz2)
    out = dict(m=inv_mass(ds + mu + four(nux, nuy, pz, 0.0)), base=f["Dsmu_mass_corr"],
               w=f["weight"], cplx=cplx, sr=np.abs(f["Ds_mass"] - M_DS) < SR_HALF_WIDTH)
    out.update({c: a[c] for c in extra})
    return out


def wquantile(x, w, qs):
    """Weighted quantiles. Signal has a handful of negative genWeight rows (7 of 869k), so
    the running sum is forced monotonic; the effect is at the 1e-5 level."""
    o = np.argsort(x)
    c = np.maximum.accumulate(np.cumsum(w[o]))
    return np.interp(np.asarray(qs) * c[-1], c, x[o])


def wfrac(x, w, lo, hi):
    return w[(x >= lo) & (x <= hi)].sum() / w.sum()


def main():
    sig = load("Signal_merged", ["channel"])
    bkg = load("WtoLNu2Jets_merged", ["cand_gen_mode"])
    print(f"signal rows with negative weight: {(sig['w'] < 0).sum()} of {len(sig['w'])}")

    mode = bkg["cand_gen_mode"]
    groups = [
        ("Signal DsW", sig, sig["channel"] == 0),
        ("Signal DsstW", sig, sig["channel"] == 1),
        ("W+jets mode 1+2 (Ds->KKpi)", bkg, (mode == 1) | (mode == 2)),
        ("W+jets mode 3+4 (D->KKpi)", bkg, (mode == 3) | (mode == 4)),
        ("W+jets mode 5 (D->Kpipi)", bkg, mode == 5),
        ("W+jets mode 1-5", bkg, (mode >= 1) & (mode <= 5)),
    ]
    sig_all = np.ones(len(sig["w"]), dtype=bool)

    for region, rmask_of in (("full Ds-mass window", lambda d: np.ones(len(d["w"]), bool)),
                             (f"SR |M(KKpi)-M_Ds|<{SR_HALF_WIDTH}", lambda d: d["sr"])):
        print(f"\n===== {region} =====")
        for var, label in (("m", "m(Ds mu nu)"), ("base", "Dsmu_mass_corr")):
            s = sig_all & rmask_of(sig)
            lo, med, hi = wquantile(sig[var][s], sig["w"][s], [0.16, 0.5, 0.84])
            print(f" {label}: signal (DsW+DsstW) median {med:6.1f}, central 68% [{lo:6.1f}, {hi:6.1f}]")
            for name, d, m in groups:
                mm = m & rmask_of(d)
                print(f"   {name:28s} rows {mm.sum():8d}  sum w {d['w'][mm].sum():12.5g}"
                      f"  frac in signal 68% interval {wfrac(d[var][mm], d['w'][mm], lo, hi):6.3f}"
                      + (f"  complex-pz {100 * np.average(d['cplx'][mm], weights=np.abs(d['w'][mm])):4.1f}%"
                         if var == "m" else ""))

    lumi_fb = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000
    stack = [   # bottom -> top, colours as draw_wp_mass.py
        (r"D$\to$K$\pi\pi$ misID (5)", mode == 5, "#e6c229"),
        (r"genuine D$\to$KK$\pi$ (3+4)", (mode == 3) | (mode == 4), "#e07a3c"),
        (r"true D$_s$ (mode 1+2)", (mode == 1) | (mode == 2), "#c0392b"),
    ]
    lines = [(r"H$\to$D$_s$W", sig["channel"] == 0, "#1f4e9c", "-"),
             (r"H$\to$D$_s^*$W", sig["channel"] == 1, "#17becf", "--")]
    variables = [("base", r"$m(D_s\mu)$ [GeV]", np.arange(0, 130 + 1e-9, 5.0), None),
                 ("m", r"$m(D_s\mu\nu)$ [GeV]", np.arange(75, 300 + 1e-9, 5.0), 125.0)]
    for tag, title, rmask_of in (("", "full Ds-mass window", lambda d: np.ones(len(d["w"]), bool)),
                                 ("_SR", rf"SR $|M(KK\pi)-M_{{D_s}}|<{SR_HALF_WIDTH*1000:.0f}$ MeV",
                                  lambda d: d["sr"])):
        fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.5), gridspec_kw=dict(wspace=0.26))
        rs, rb = rmask_of(sig), rmask_of(bkg)
        for ax, (var, xlabel, bins, refline) in zip(axes, variables):
            bc, width = 0.5 * (bins[1:] + bins[:-1]), bins[1] - bins[0]

            def hist(x, w):
                x = np.clip(x, bins[0], bins[-1] - 1e-9)   # last bin includes overflow
                h, _ = np.histogram(x, bins=bins, weights=w)
                e2, _ = np.histogram(x, bins=bins, weights=w ** 2)
                assert np.isclose(h.sum(), w.sum()), "events lost outside the bin range"
                return h, np.sqrt(e2)

            bottom, err2 = np.zeros(len(bc)), np.zeros(len(bc))
            for lab, m, cl in stack:
                mm = m & rb
                h, e = hist(bkg[var][mm], bkg["w"][mm])
                ax.bar(bc, h, width=width, bottom=bottom, color=cl, alpha=0.9, edgecolor="none",
                       label=f"W+jets {lab}  N={h.sum():,.0f}")
                bottom, err2 = bottom + h, err2 + e ** 2
            err = np.sqrt(err2)
            ax.bar(bc, 2 * err, width=width, bottom=bottom - err, color="none", hatch="////",
                   edgecolor="0.35", linewidth=0.0, zorder=3, label="W+jets MC stat")

            s_all = rs
            scale = bottom.sum() / sig["w"][s_all].sum()
            scale = float(f"{scale:.1g}")
            for lab, m, cl, ls in lines:
                mm = m & rs
                h, _ = hist(sig[var][mm], sig["w"][mm])
                ax.step(bins, np.append(h, h[-1]) * scale, where="post", color=cl, ls=ls, lw=2.2,
                        zorder=4, label=f"{lab} ($\\times${scale:.0e})  N={h.sum():.3g}")

            lo, hi = wquantile(sig[var][s_all], sig["w"][s_all], [0.16, 0.84])
            b15 = (mode >= 1) & (mode <= 5) & rb
            fb = wfrac(bkg[var][b15], bkg["w"][b15], lo, hi)
            ax.axvspan(lo, hi, color="grey", alpha=0.12, zorder=0)
            ax.text(0.03, 0.72, f"signal central 68%: [{lo:.0f}, {hi:.0f}] GeV\n"
                    f"W+jets mode 1-5 inside: {fb:.3f}", transform=ax.transAxes,
                    fontsize=11, va="top", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
            if refline:
                ax.axvline(refline, color="green", ls=":", lw=1.4)
            ax.set_xlim(bins[0], bins[-1])
            ax.set_ylim(0, bottom.max() * 1.9)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(f"Events / {width:.0f} GeV")
            ax.legend(fontsize=9, loc="upper right")
            hep.cms.label(ax=ax, data=False, text="Preliminary",
                          rlabel=rf"{lumi_fb:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=11)
        fig.suptitle(f"Signal vs W+jets GenTag mode 1-5, {title}", fontsize=14)
        out = common.fig_path("features", f"dsmunu_mass_signal_vs_gentag{tag}.png")
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()

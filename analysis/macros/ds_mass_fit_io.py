"""Ds_mass template fit of MC scale factors, I/O half (conda ds_dnn): build inputs, draw postfit.

    conda activate ds_dnn
    python3 ds_mass_fit_io.py inputs --key=ssdata_phi20 --cuts=0.70,0.60,0.50 --joint=slices --float=r34,r5,k
    (cmsenv) python3 fit_ds_mass_roofit.py <fit dir printed above>
    python3 ds_mass_fit_io.py draw   --key=ssdata_phi20 --cuts=0.70,0.60,0.50 --joint=slices --float=r34,r5,k
or all three at once: run_ds_mass_fit.sh <same arguments>

The fit itself (RooFit/HistFactory) is fit_ds_mass_roofit.py, run under cmsenv: ds_dnn has no
ROOT and cmsenv has no torch (common.py imports it), so the two halves talk through
<model dir>/fit/ds_mass/<tag>/{inputs.npz, meta.json, result.json}.

Model (same as fit_sideband_scales_stage1.py, 2026-10-05, extended by the user 2026-10-16):
    OS data = W+jets(mode 0,6-10) + TT                                   [fixed at MC]
              + r12 W+jets(1+2) + r34 W+jets(3+4) + r5 W+jets(5) + k DY   [floated if in --float]
A component whose factor is not in --float stays at MC (factor 1). r12 is NOT floated by default
(user, 2026-10-16): mode 1+2 sits ~90% inside the blinded SR, so the sidebands cannot constrain it --
floated, it was pulled to ~0 by the peak tails (MINOS [0, ~0.6-0.9]) and dragged k via rho = -0.4. Factors are SHARED by every
fitted region (feedback_shared_scale_factors); only this joint fit is done, no per-region fits (user,
2026-10-16). Only the Ds-mass sidebands are fitted (|M-M_Ds|>=40 MeV,
the OS blind), 20 MeV bins by default (--bin=0.010 for 10 MeV; the SR stays whole bins).

--cuts=e1,e2,... (eps_sig, descending) and --joint picks the fitted regions:
  slices  DISJOINT score slices [thr(e1), thr(e2)), ..., [thr(en), inf) -- each event enters once;
  nested  the nested cuts eps>=e1, ..., eps>=en together. They share events (eps50 is a subset of
          eps70), so an event above the tightest cut enters the likelihood once per region containing
          it: errors are underestimated and the tight end is over-weighted. Kept because the user asked
          for it as an option; the figure title says so.
The draw step always shows the nested cuts with the shared factors; with --joint=slices it also
shows the slices that were actually fitted.

Selection = common.selection_mask + score_dsmu > -1, plus the key's phi_win if its training
config has one (same population the stage-1 model was trained on); thresholds from
common.stage1_threshold (which applies the same window).
"""
import json
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common  # noqa: E402

ARG = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[2:] if "=" in a}
MODE = sys.argv[1] if len(sys.argv) > 1 else ""
assert MODE in ("inputs", "draw"), "usage: ds_mass_fit_io.py inputs|draw --key=K --cuts=e1,e2 --joint=slices|nested --float=p1,p2"
KEY = ARG.get("--key", "ssdata_phi20")
EFFS = [float(x) for x in ARG.get("--cuts", "0.70,0.60,0.50").split(",")]
assert EFFS == sorted(EFFS, reverse=True), "--cuts must be descending eps_sig"
JOINT = ARG.get("--joint", "slices")
assert JOINT in ("slices", "nested"), "--joint must be slices or nested"
FLOAT = [x for x in ARG.get("--float", "r34,r5,k").split(",") if x]
TAG = ARG.get("--tag", f"{JOINT}{'-'.join(f'{100 * e:.0f}' for e in EFFS)}_f{'-'.join(FLOAT)}"
              + ("" if ARG.get("--bin", "0.020") == "0.020" else f"_b{1e3 * float(ARG['--bin']):.0f}"))
FITDIR = f"{common.out_dir_for(KEY)}/fit/ds_mass/{TAG}"
PRED = f"{common.out_dir_for(KEY)}/Prediction"
PHI_WIN = common.resolve(KEY).get("phi_win")

M_DS, M_DP, SR_HALF = 1.96835, 1.8696, 0.040
BIN_W = float(ARG.get("--bin", "0.020"))   # GeV; must divide the 40 MeV SR half-width so the blind is whole bins
assert abs(round(SR_HALF / BIN_W) * BIN_W - SR_HALF) < 1e-9, "--bin must divide 0.040"
BINS = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, BIN_W)
BC = 0.5 * (BINS[1:] + BINS[:-1])
VIS = np.abs(BC - M_DS) >= SR_HALF
# (name, legend, colour, factor); stack bottom -> top, colours as fit_sideband_scales_stage1.py
COMPS = [
    ("comb", r"W+jets other (mode 0,6-10)", "#5a9e78", None),
    ("TT",   r"$t\bar{t}$",                 "#7f5fbf", None),
    ("DY",   r"DY$\to\ell\ell$",            "#4a90c2", "k"),
    ("5",    r"W+jets D$\to$K$\pi\pi$ misID (5)",   "#e6c229", "r5"),
    ("34",   r"W+jets genuine D$\to$KK$\pi$ (3+4)", "#e07a3c", "r34"),
    ("12",   r"W+jets true D$_s$ (1+2)",            "#c0392b", "r12"),
]
COMPS = [(g, l, c, p if p in FLOAT else None) for g, l, c, p in COMPS]
bad = set(FLOAT) - {p for *_x, p in COMPS if p}
assert not bad, f"--float: unknown factor(s) {bad}; choices r12,r34,r5,k"
WMODES = {"comb": [0, 6, 7, 8, 9, 10], "12": [1, 2], "34": [3, 4], "5": [5]}
PAR_TEX = {"r12": "r12 (W+jets 1+2)", "r34": "r34 (W+jets 3+4)", "r5": "r5 (W+jets 5)", "k": "k (DY)"}


def make_inputs():
    def rd(name, extra):
        with uproot.open(f"{PRED}/{name}.root")[common.TREE] as t:
            a = t.arrays(["Ds_mass", "phi_mass", "score_dsmu"] + extra + list(common.SELECTION), library="np")
        m = common.selection_mask(a) & (a["score_dsmu"] > -1)
        if PHI_WIN is not None:
            m &= np.abs(a["phi_mass"] - common.M_PHI) < PHI_WIN
        out = {x: a[x][m] for x in ["Ds_mass", "score_dsmu"] + extra}
        n_out = (~((out["Ds_mass"] >= BINS[0]) & (out["Ds_mass"] < BINS[-1]))).sum()
        print(f"  {name:20s} {len(m):>9,} -> {m.sum():>9,} rows ({n_out} outside the mass range)")
        return out

    print(f"key {KEY}  phi_win {PHI_WIN}  -> {FITDIR}")
    W = rd("WtoLNu2Jets_merged", ["weight", "cand_gen_mode"])
    assert np.isin(W["cand_gen_mode"], sum(WMODES.values(), [])).all(), "cand_gen_mode outside the groups"
    S = {"TT": rd("TT_merged", ["weight"]), "DY": rd("DY_merged", ["weight"])}
    X = rd("Data_OS_merged", [])
    X = {k: v[np.abs(X["Ds_mass"] - M_DS) >= SR_HALF] for k, v in X.items()}   # blind, event level

    thr = {e: common.stage1_threshold(dict(key=KEY, tool="DNN", seed=42, eff_sig=e))[0] for e in EFFS}
    edges = [thr[e] for e in EFFS] + [np.inf]
    regions = {f"slice{i}": (edges[i], edges[i + 1]) for i in range(len(EFFS))}
    regions.update({f"nested{i}": (thr[e], np.inf) for i, e in enumerate(EFFS)})

    def h(D, m):
        w = D["weight"][m].astype(float) if "weight" in D else np.ones(m.sum())
        return (np.histogram(D["Ds_mass"][m], BINS, weights=w)[0],
                np.histogram(D["Ds_mass"][m], BINS, weights=w ** 2)[0])

    arr = {"bins": BINS, "vis": VIS}
    for r, (lo, hi) in regions.items():
        sel = lambda D: (D["score_dsmu"] > lo) & (D["score_dsmu"] <= hi)
        arr[f"{r}_data"] = np.histogram(X["Ds_mass"][sel(X)], BINS)[0].astype(float)
        for g, ms in WMODES.items():
            arr[f"{r}_{g}_w"], arr[f"{r}_{g}_w2"] = h(W, sel(W) & np.isin(W["cand_gen_mode"], ms))
        for g in ("TT", "DY"):
            arr[f"{r}_{g}_w"], arr[f"{r}_{g}_w2"] = h(S[g], sel(S[g]))
        tot = sum(arr[f"{r}_{g}_w"] for g, *_ in COMPS)
        neg = {g: int((arr[f"{r}_{g}_w"][VIS] < 0).sum()) for g, *_ in COMPS}
        f12 = arr[f"{r}_12_w"][VIS].sum() / arr[f"{r}_12_w"].sum()
        print(f"  {r:8s} score in ({lo:.6f}, {hi:.6f}]: data(SB) {arr[f'{r}_data'][VIS].sum():8,.0f}  "
              f"MC(SB) {tot[VIS].sum():10,.1f}  mode1+2 outside SR {100 * f12:.1f}%  "
              f"negative SB bins {', '.join(f'{g}:{n}' for g, n in neg.items() if n) or 'none'}")
        assert (tot[VIS] > 0).all(), f"{r}: total MC <= 0 in a sideband bin"

    os.makedirs(FITDIR, exist_ok=True)
    np.savez(f"{FITDIR}/inputs.npz", **arr)
    meta = dict(key=KEY, phi_win=PHI_WIN, effs=EFFS, thresholds=[thr[e] for e in EFFS], float=FLOAT,
                comps=[[g, p] for g, _l, _c, p in COMPS], slices=[f"slice{i}" for i in range(len(EFFS))],
                nested=[f"nested{i}" for i in range(len(EFFS))], joint=JOINT,
                fit_regions=[f"{'slice' if JOINT == 'slices' else 'nested'}{i}" for i in range(len(EFFS))],
                bin_width=BIN_W, sr_half=SR_HALF)
    with open(f"{FITDIR}/meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"Wrote {FITDIR}/inputs.npz, meta.json")


def draw():
    """The shared-factor joint fit, drawn:
    ds_mass_fit_<tag>.png         nested cuts eps>=e1..en with the shared factors;
    ds_mass_fit_<tag>_slices.png  (--joint=slices only) the disjoint slices that were fitted."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.style.CMS)

    A = np.load(f"{FITDIR}/inputs.npz")
    with open(f"{FITDIR}/result.json") as f:
        R = json.load(f)
    lumi = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
    width = BINS[1] - BINS[0]

    def panel(fig, spec, r, P, title, cms):
        """Postfit stack of region r with factors P (result.json params dict); returns (chi2, chi2_pre)."""
        S = {p: P[p]["val"] for p in FLOAT}
        dh = A[f"{r}_data"]
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=spec, height_ratios=[3, 1], hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        bottom, post, post_e2, pre = (np.zeros(len(BC)) for _ in range(4))
        for g, lab, cl, p in COMPS:
            s = 1.0 if p is None else S[p]
            hh, e2 = A[f"{r}_{g}_w"], A[f"{r}_{g}_w2"]
            ax.bar(BC, s * hh, width=width, bottom=bottom, color=cl, alpha=0.9, edgecolor="none",
                   label=lab + (f" ×{s:.2f}" if p else "") + f"  N={s * hh.sum():,.0f}")
            bottom += s * hh
            post += s * hh; post_e2 += s * s * e2; pre += hh
        pe = np.sqrt(post_e2)
        ax.bar(BC, 2 * pe, width=width, bottom=post - pe, color="none", hatch="////", edgecolor="0.35",
               linewidth=0.0, zorder=3, label="MC stat (postfit)")
        dm = np.where(VIS, dh, np.nan)
        ax.errorbar(BC, dm, yerr=np.sqrt(np.where(VIS, dh, 0)), fmt="o", ms=4, color="black",
                    label=f"OS data (Ds SR excl.)  N={np.nansum(dm):,.0f}")
        sig = np.sqrt(np.maximum(dh + post_e2, 1.0))[VIS]
        chi2 = (((dh - post)[VIS] / sig) ** 2).sum()
        chi2_pre = (((dh - pre)[VIS] / np.sqrt(np.maximum(dh + sum(A[f"{r}_{g}_w2"] for g, *_ in COMPS), 1.0))[VIS]) ** 2).sum()
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)
        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
        ax.axvline(M_DP, color="green", ls=":", lw=1.4)
        ax.set_xlim(BINS[0], BINS[-1])
        ax.set_ylim(0, max(bottom.max(), np.nanmax(dm)) * 1.9)
        ax.tick_params(labelbottom=False)
        ax.legend(fontsize=10, loc="upper right")
        ax.text(0.03, 0.97, f"$\\chi^2$/nbins = {chi2:.1f}/{VIS.sum()}\n(prefit {chi2_pre:.1f}/{VIS.sum()})",
                transform=ax.transAxes, va="top", ha="left", fontsize=13)
        ax.set_title(title, fontsize=16, pad=24)
        ax.set_ylabel(f"Events / {width * 1000:.0f} MeV", fontsize=18)
        if cms:
            hep.cms.label(ax=ax, data=True, text="Preliminary",
                          rlabel=rf"{lumi / 1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=13)
        with np.errstate(invalid="ignore", divide="ignore"):
            ok = VIS & (post > 0)
            r_post, r_pre = np.where(ok, dh / post, np.nan), np.where(ok, dh / pre, np.nan)
            re_post, re_pre = np.where(ok, np.sqrt(dh) / post, np.nan), np.where(ok, np.sqrt(dh) / pre, np.nan)
            band = np.where(ok, pe / post, np.nan)
        rax.bar(BC, 2 * band, width=width, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(BC - 0.002, r_pre, yerr=re_pre, fmt="o", ms=3.5, mfc="white", color="0.45", elinewidth=1.0, label="prefit")
        rax.errorbar(BC + 0.002, r_post, yerr=re_post, fmt="o", ms=3.5, color="black", elinewidth=1.0, label="postfit")
        n_off = int(np.sum(VIS & ((np.nan_to_num(r_pre, nan=1) > 1.99) | (np.nan_to_num(r_post, nan=1) > 1.99))))
        if n_off:
            rax.text(0.02, 0.92, f"{n_off} point(s) off scale", transform=rax.transAxes, fontsize=11, va="top")
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        rax.set_ylim(0.0, 1.99)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.set_xlim(BINS[0], BINS[-1])
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]", fontsize=18)
        rax.set_ylabel("data / MC", fontsize=16)
        rax.legend(fontsize=11, loc="lower right", ncol=2)
        return chi2, chi2_pre

    def fmt(P, p):
        q = P[p]
        lo, hi = q.get("minos_lo"), q.get("minos_hi")
        e = f"$^{{{hi:+.2f}}}_{{{lo:+.2f}}}$" if lo is not None else f" ± {q['err']:.2f}"
        return f"{PAR_TEX[p]} = {q['val']:.2f}{e}" + (" [bound]" if q.get("at_bound") else "")

    PS = R["shared"]["params"]
    sl_names = [f"{100 * a:.0f}-{100 * b:.0f}%" for a, b in zip(EFFS, EFFS[1:])] + [f">{100 * EFFS[-1]:.0f}%"]
    phi_txt = "" if PHI_WIN is None else rf", $|M(KK)-m_\phi|<{1e3 * PHI_WIN:.0f}$ MeV"
    how = (f"joint fit on disjoint slices {', '.join(sl_names)}" if JOINT == "slices" else
           f"joint fit on nested cuts {', '.join(f'>={100 * e:.0f}%' for e in EFFS)} "
           "(regions share events: errors underestimated)")
    head = (f"Ds-mass sideband fit ({KEY}{phi_txt}), {how}\n"
            "shared: " + ",  ".join(fmt(PS, p) for p in FLOAT)
            + f"   ($\\chi^2$/ndf = {R['shared']['chi2']:.1f}/{R['shared']['ndf']}, MC stat: Barlow-Beeston-lite)")

    figs = [("", [(f"nested{i}", rf"$\varepsilon_{{sig}}$$\geq${100 * e:.0f}% — shared-factor postfit")
                  for i, e in enumerate(EFFS)])]
    if JOINT == "slices":
        figs.append(("_slices", [(f"slice{i}", f"slice {n} — shared-factor postfit") for i, n in enumerate(sl_names)]))
    for suffix, cols in figs:
        fig = plt.figure(figsize=(8.0 * len(cols), 8.0))
        outer = gridspec.GridSpec(1, len(cols), figure=fig, wspace=0.30)
        for col, (r, ttl) in enumerate(cols):
            c2, c2p = panel(fig, outer[col], r, PS, ttl, col == 0)
            print(f"{r}: chi2/nbins post {c2:.1f}/{VIS.sum()}  pre {c2p:.1f}/{VIS.sum()}")
        fig.suptitle(head, fontsize=16, y=1.08)
        out = common.fig_path("closure", f"ds_mass_fit_{TAG}{suffix}.png", key=KEY)
        plt.savefig(out, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {out}")


make_inputs() if MODE == "inputs" else draw()

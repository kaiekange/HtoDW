"""Stage-1 DNN features in the Z->mumu CR, data vs MC, no cut and stage-1 eps90, with the
SR signal shape as a reference.

Z window 76 < M(mumu) < 106, |phi_mass - M_phi| < 20 MeV, OS data, MC = DY + ttbar + W+jets
x common.br_reweight x common.dy_jetbin_sf (as draw_cr_zmumu.py's mass figures), read from
the ssdata_phi20 Prediction/zmumu_<group>.root files (score_dsmu needed for the eps90 panel).
Signal = the SR Signal_merged (Prediction/Signal_merged.root), selected as in
common.stage1_threshold (selection_mask, Ds/mu truth match, score > -1, phi window), drawn as
a line normalised to the data in each panel -- a shape reference only, it is not in this CR.
Sentinel / non-finite values (< -5, e.g. bscore -10, -777) are left out of the histograms and
their fractions printed. Integer-valued features get unit bins; the rest 40 bins over the
1-99% quantiles of data and signal together (edge bins hold under/overflow).
Written 2026-10-20 after the feature-by-feature reweighting scan ranked
mu_IsoR04_sumNeutralHadronEt first.

    python3 draw_zmumu_feature.py [--var=mu_IsoR04_sumNeutralHadronEt | --all] [--key=ssdata_phi20] [--nocut]
Output: <model>/figures/sf/zmumu_features/<var>.png
--nocut: only the no-cut panel, written to <model>/figures/sf/zmumu_features_nocut/<var>.png
--raw-vs-corr: every corrected feature (the *_corr ones and the three dphi_met_*), no cut,
  left = before correction (raw nano branch; dphi from raw met_phi / Dsmu_phi), right = after,
  same binning; written to <model>/figures/sf/zmumu_features_raw_vs_corr/<var>.png
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

ARGS = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
MODEL_KEY = ARGS.get("key", "ssdata_phi20")
MODEL_DIR = common.out_dir_for(MODEL_KEY)
PRED = f"{MODEL_DIR}/Prediction"
FEATS = [l.strip() for l in open(f"{MODEL_DIR}/feature_names.txt") if l.strip()]
VARS = FEATS if "--all" in sys.argv else [ARGS.get("var", "mu_IsoR04_sumNeutralHadronEt")]
RAWCORR = "--raw-vs-corr" in sys.argv
DER = dict(common.resolve(MODEL_KEY)["features"].derived)   # this model's derived features
PLOT = VARS
if RAWCORR:
    DER.update({
        "dphi_met_mu_raw":   (lambda a: common._dphi(a["met_phi"], a["mu_phi"]), ("met_phi", "mu_phi")),
        "dphi_met_Ds_raw":   (lambda a: common._dphi(a["met_phi"], a["Ds_phi"]), ("met_phi", "Ds_phi")),
        "dphi_met_Dsmu_raw": (lambda a: common._dphi(a["met_phi"], a["Dsmu_phi"]), ("met_phi", "Dsmu_phi")),
    })
    RAW_OF = {f: f[:-len("_corr")] for f in FEATS if f.endswith("_corr")}
    RAW_OF.update({f: f + "_raw" for f in common.DERIVED})
    PLOT = [ARGS["var"]] if "var" in ARGS else [f for f in FEATS if f in RAW_OF]
    VARS = PLOT + [RAW_OF[f] for f in PLOT]
THR, _n = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=0.90))
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
SENT = -5.0
PROC = [  # (label, group, extra mask, colour) -- stack bottom -> top, as draw_cr_zmumu.py
    (r"$t\bar{t}$ LNu2Q", "TT_merged", lambda a: a["tt_decay"] == 1, "#b39ddb"),
    (r"$t\bar{t}$ 2L2Nu", "TT_merged", lambda a: a["tt_decay"] == 0, "#7f5fbf"),
    ("W+jets", "WtoLNu2Jets_merged", None, "#5a9e78"),
    (r"DY$\to\ell\ell$", "DY_merged", None, "#4a90c2"),
]
BRANCHES = set()
for f in VARS:
    BRANCHES |= set(DER[f][1]) if f in DER else {f}


def feats(a):
    return {f: (np.asarray(DER[f][0](a), float) if f in DER
                else a[f].astype(float)) for f in VARS}


def read_cr(group, mc):
    with uproot.open(f"{PRED}/zmumu_{group}.root")["Events"] as t:
        keys = set(t.keys())
        want = BRANCHES | {"phi_mass", "m_mumu", "weight", "score_dsmu"}
        if mc:
            want |= {"tt_decay", "cand_gen_mode", "cand_gen_common_hf_pdgId", "jet_bin"}
        a = t.arrays([c for c in want if c in keys], library="np")
    m = (np.abs(a["m_mumu"] - 91.19) < 15) & (np.abs(a["phi_mass"] - common.M_PHI) < 0.020)
    a = {k: v[m] for k, v in a.items()}
    w = np.ones(m.sum()) if not mc else a["weight"].astype(float) * common.br_reweight(a)
    if mc and group == "DY_merged":
        w = w * common.dy_jetbin_sf(a)
    out = feats(a)
    out.update(w=w, score=a["score_dsmu"], tt_decay=a.get("tt_decay"))
    return out


def read_signal():
    with uproot.open(f"{PRED}/Signal_merged.root")[common.TREE] as t:
        want = BRANCHES | {"phi_mass", "weight", "score_dsmu", "Ds_match", "mu_match", *common.SELECTION}
        a = t.arrays(list(want), library="np")
    m = (common.selection_mask(a) & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
         & (a["score_dsmu"] > -1) & (np.abs(a["phi_mass"] - common.M_PHI) < 0.020))
    a = {k: v[m] for k, v in a.items()}
    out = feats(a)
    out.update(w=a["weight"].astype(float), score=a["score_dsmu"])
    return out


D = read_cr("Data_OS_merged", False)
S = {g: read_cr(g, True) for g in ("TT_merged", "WtoLNu2Jets_merged", "DY_merged")}
SIG = read_signal()
NOCUT = "--nocut" in sys.argv
OUT_DIR = (os.path.dirname(common.fig_path("sf", "x", key=MODEL_KEY))
           + ("/zmumu_features_raw_vs_corr" if RAWCORR else "/zmumu_features_nocut" if NOCUT
              else "/zmumu_features"))
os.makedirs(OUT_DIR, exist_ok=True)
CUTS = (("no cut", None),) if NOCUT else (("no cut", None), (r"stage-1 $\varepsilon_{sig}$ 90%", THR))


def panels(pvar):  # (legend title, variable drawn, score cut)
    if RAWCORR:
        return (("no cut, before correction", RAW_OF[pvar], None), ("no cut, after correction", pvar, None))
    return tuple((cname, pvar, cut) for cname, cut in CUTS)
print(f"{'feature':32s} | {'cut':7s} | mean data / MC / sig | frac sentinel data / MC / sig")


def make_bins(vs):
    xs = [x[np.isfinite(x) & (x >= SENT)] for v in vs for x in (D[v], SIG[v])]
    allx = np.concatenate(xs)
    if np.all(np.equal(np.mod(allx[:20000], 1), 0)):   # integer-valued
        lo, hi = np.floor(np.quantile(allx, 0.001)), np.ceil(np.quantile(allx, 0.999))
        return np.arange(lo - 0.5, hi + 1.5, 1.0)
    lo, hi = np.quantile(allx, 0.01), np.quantile(allx, 0.99)
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, 41)


for pvar in PLOT:
    PAN = panels(pvar)
    bins = make_bins([p[1] for p in PAN])
    lo, hi = bins[0], bins[-1]
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    clip = lambda x: np.clip(x, lo, hi - 1e-9 * (hi - lo))
    good = lambda x: np.isfinite(x) & (x >= SENT)
    fig = plt.figure(figsize=(7.2 * len(PAN), 7.6))
    outer = gridspec.GridSpec(1, len(PAN), figure=fig, wspace=0.30, top=0.90, bottom=0.10)
    for col, (cname, var, cut) in enumerate(PAN):
        sel = lambda a: np.ones(len(a["w"]), bool) if cut is None else a["score"] > cut
        ds, ss = sel(D), sel(SIG)
        mc_x = np.concatenate([S[g][var][sel(S[g])] for g in S])
        mc_w = np.concatenate([S[g]["w"][sel(S[g])] for g in S])
        gd, gm, gs = good(D[var][ds]), good(mc_x), good(SIG[var][ss])
        print(f"{var:32s} | {('noCut' if cut is None else 'eps90'):7s} | "
              f"{D[var][ds][gd].mean():9.3g} {np.average(mc_x[gm], weights=mc_w[gm]):9.3g} "
              f"{np.average(SIG[var][ss][gs], weights=SIG['w'][ss][gs]):9.3g} | "
              f"{1 - gd.mean():.3f} {mc_w[~gm].sum() / mc_w.sum():.3f} "
              f"{SIG['w'][ss][~gs].sum() / SIG['w'][ss].sum():.3f}")
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[col], height_ratios=[3, 1], hspace=0.06)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        bottom, tot, tot2 = np.zeros(len(bc)), np.zeros(len(bc)), np.zeros(len(bc))
        for lab, g, extra, colr in PROC:
            a = S[g]
            k = sel(a) & good(a[var]) & (extra(a) if extra else True)
            h = np.histogram(clip(a[var][k]), bins, weights=a["w"][k])[0]
            e2 = np.histogram(clip(a[var][k]), bins, weights=a["w"][k] ** 2)[0]
            ax.bar(bc, h, width=wd, bottom=bottom, color=colr, edgecolor="none", label=f"{lab}  N={h.sum():,.1f}")
            bottom += h; tot += h; tot2 += e2
        terr = np.sqrt(tot2)
        ax.bar(bc, 2 * terr, width=wd, bottom=tot - terr, color="none", hatch="////", edgecolor="0.35",
               linewidth=0.0, label="MC stat")
        dh = np.histogram(clip(D[var][ds][gd]), bins)[0]
        hs = np.histogram(clip(SIG[var][ss][gs]), bins, weights=SIG["w"][ss][gs])[0]
        if hs.sum() > 0:
            hs = hs * dh.sum() / hs.sum()
        hep.histplot(hs, bins, ax=ax, color="red", lw=2, label="SR signal (shape, norm. to data)")
        ax.errorbar(bc, dh, yerr=np.sqrt(dh), fmt="o", ms=4, color="black", label=f"OS data  N={dh.sum():,}")
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, 1.35 * max(tot.max(), dh.max(), hs.max(), 1))
        ax.legend(fontsize=9, loc="upper right", title=cname, title_fontsize=11)
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("Events / bin", fontsize=16)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(tot > 0, dh / tot, np.nan)
            re = np.where(tot > 0, np.sqrt(dh) / tot, np.nan)
            band = np.where(tot > 0, terr / tot, np.nan)
        rax.bar(bc, 2 * band, width=wd, bottom=1 - band, color="none", hatch="////", edgecolor="0.35", linewidth=0.0)
        rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black")
        rax.axhline(1, color="red", lw=1)
        rax.set_ylim(0.0, 2.0)
        rax.set_xlim(lo, hi)
        rax.set_xlabel(var, fontsize=14)
        rax.set_ylabel("data/MC", fontsize=14)
        if col == 0:
            hep.cms.label(ax=ax, data=True, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6, fontsize=14)
    fig.suptitle(r"Z$\to\mu\mu$ CR, 76 < M($\mu\mu$) < 106, $\phi$20: MC x BR rw x DY jet-bin SF;"
                 r" red = SR signal shape", fontsize=15, y=0.975)
    plt.savefig(f"{OUT_DIR}/{pvar}.png", dpi=90, bbox_inches="tight")
    plt.close(fig)
print(f"\nSaved {len(PLOT)} figures in {OUT_DIR}")

"""DY MC only: Ds_mass split by where the best muon comes from (Muon_genPartFlav: 15 = from
tau, 1 = prompt, 0/3/4/5 = other), at no cut and eps_sig 90/70/50 of a stage-1 model.

The per-row muon origin comes from the cache draw_dy_mu_origin.py built on the pre-phi20
DY_merged (MVA_ROOT/aux/DY_merged_mu_genPartFlav.npz, key (lumi<<32)|event, run == 1 for MC).
The current DY_merged is a subset of those events (phi window + the 65 dropped JOBIDs), so
the rows are looked up by that key; every row must be found and keys must be unique.

Selection: common.selection_mask, score > -1, the model's phi window; the Ds SR is NOT
removed (MC only). DY weights include common.dy_jetbin_sf (Z->mumu CR factors).
Printed: per window (low <1.830 / D+- 1.83-1.91 / SR / high >2.008) the yield of each origin
and its share of DY.

    python3 draw_dy_tau_fraction.py [--key=ssdata_phi20]
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

MODEL_KEY = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--key=")), "ssdata_phi20")
PHI_WIN = common.resolve(MODEL_KEY).get("phi_win")
CACHE = f"{common.MVA_ROOT}/aux/DY_merged_mu_genPartFlav.npz"
M_DS, SR_HALF = 1.96835, 0.040
LUMI_FB = (common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]) / 1000.0
CUTS = [(None, "no cut"), (0.90, r"$\varepsilon_{sig}$=90%"), (0.70, r"$\varepsilon_{sig}$=70%"),
        (0.50, r"$\varepsilon_{sig}$=50%")]
DS_NSIDE = {None: 26, 0.90: 13, 0.70: 10, 0.50: 5}   # same widths as draw_stage1_mass_by_score_cut
ORIGINS = [("other (flav 0,3,4,5)", [0, 3, 4, 5], "#bdbdbd"),
           ("prompt $\\mu$ (flav 1)", [1], "#4a90c2"),
           ("$\\mu$ from $\\tau$ (flav 15)", [15], "#e07a3c")]
WINDOWS = [("low <1.830", lambda x: x < 1.830),
           ("D+- 1.83-1.91", lambda x: (x >= 1.830) & (x < 1.910)),
           ("SR", lambda x: np.abs(x - M_DS) < SR_HALF),
           ("high >2.008", lambda x: x > M_DS + SR_HALF)]


def ds_edges(eff):
    n = DS_NSIDE[eff]
    w = (0.3 - SR_HALF) / n
    n_sr = int(round(2 * SR_HALF / w))
    n_sr = n_sr if abs(n_sr * w - 2 * SR_HALF) < 1e-9 else 1
    return np.concatenate([np.linspace(M_DS - 0.3, M_DS - SR_HALF, n + 1),
                           np.linspace(M_DS - SR_HALF, M_DS + SR_HALF, n_sr + 1)[1:-1],
                           np.linspace(M_DS + SR_HALF, M_DS + 0.3, n + 1)])


with uproot.open(f"{common.out_dir_for(MODEL_KEY)}/Prediction/DY_merged.root")[common.TREE] as t:
    a = t.arrays(list(dict.fromkeys(["Ds_mass", "phi_mass", "score_dsmu", "weight", "jet_bin",
                                     "run", "luminosityBlock", "event", *common.SELECTION])),
                 library="np")
assert np.all(a["run"] == 1), "MC run number is not 1 -- the cache key would be ambiguous"
rowkey = (a["luminosityBlock"].astype(np.int64) << 32) | a["event"].astype(np.int64)
c = np.load(CACHE)
order = np.argsort(c["rowkey"])
ck = c["rowkey"][order]
pos = np.searchsorted(ck, rowkey)
ok = (pos < len(ck)) & (ck[np.minimum(pos, len(ck) - 1)] == rowkey)
assert ok.all() and len(np.unique(ck)) == len(ck), \
    f"{(~ok).sum()} DY rows not in {CACHE} (or duplicate keys) -- rebuild it"
flav = c["flav"][order][pos]

m = common.selection_mask(a) & (a["score_dsmu"] > -1)
if PHI_WIN is not None:
    m &= np.abs(a["phi_mass"] - common.M_PHI) < PHI_WIN
mass, score, fl = a["Ds_mass"][m], a["score_dsmu"][m], flav[m]
w = a["weight"][m].astype(float) * common.dy_jetbin_sf({"jet_bin": a["jet_bin"][m]})
print(f"DY rows after selection: {m.sum():,}  (all found in the muon-origin cache)")

thr = {None: -1.0}
for e, _ in CUTS[1:]:
    thr[e], _ = common.stage1_threshold(dict(key=MODEL_KEY, tool="DNN", seed=42, eff_sig=e))

fig, axes = plt.subplots(1, len(CUTS), figsize=(7.2 * len(CUTS), 6.8))
for ax, (eff, lbl) in zip(axes, CUTS):
    k = score > thr[eff]
    print(f"\n== {lbl.replace('$', '').replace(chr(92), '')}: DY N={w[k].sum():,.1f}")
    for wn, wf in WINDOWS:
        kw = k & wf(mass)
        tot = w[kw].sum()
        parts = "  ".join(f"{n.split(' (')[0].replace('$', '').replace(chr(92), '')} "
                          f"{w[kw & np.isin(fl, f)].sum():9,.1f} ({100 * w[kw & np.isin(fl, f)].sum() / tot:5.1f}%)"
                          for n, f, _c in ORIGINS[::-1])
        print(f"   {wn:14s} DY {tot:10,.1f}  |  {parts}")
    bins = ds_edges(eff)
    bc, wd = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    scale = wd[0] / wd
    bottom = np.zeros(len(bc))
    for n, f, col in ORIGINS:
        sel = k & np.isin(fl, f)
        hh, _ = np.histogram(mass[sel], bins, weights=w[sel])
        ax.bar(bc, hh * scale, width=wd, bottom=bottom, color=col, edgecolor="none",
               label=f"{n}  N={hh.sum():,.0f}")
        bottom += hh * scale
    ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.45 * bottom.max())
    ax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
    ax.set_ylabel(f"Events / {wd[0] * 1000:.0f} MeV", fontsize=16)
    ax.set_title(lbl, fontsize=15, pad=24 if eff is None else 6)
    ax.legend(fontsize=11, loc="upper right")
hep.cms.label(ax=axes[0], data=False, text="Preliminary", lumi=f"{LUMI_FB:.1f}", com=13.6,
              fontsize=14)
fig.suptitle(r"DY$\to\ell\ell$ MC (x jet-bin SF from Z$\to\mu\mu$ CR): $D_s$ mass by best-$\mu$ origin"
             + (f", |M(KK)-M$_\\phi$|<{PHI_WIN * 1000:.0f} MeV" if PHI_WIN else ""), fontsize=15)
out = common.fig_path("sf", "dy_ds_mass_mu_origin.png", key=MODEL_KEY)
plt.savefig(out, dpi=100, bbox_inches="tight")
print(f"\nSaved {out}")

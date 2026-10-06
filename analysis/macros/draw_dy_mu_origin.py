"""DY split by where the best muon comes from (Muon_genPartFlav): SS-trained vs OS-trained cut.

    python3 draw_dy_mu_origin.py [KEY_OS]      (default osdata; compared with ssdata)

DY_merged has no muon truth, but the raw custom nano does (standard NanoAOD Muon_genPartFlav:
1 prompt, 15 from tau, 5 b, 4 c, 3 light, 0 unmatched). Step 1 reads run/lumi/event and the
Muon_* collection from every raw DY nano file (local process pool, only events present in
DY_merged), joins on (era, jet_bin, lumi, event) and picks the best muon by matching the
candidate's uncorrected mu_pt/eta/phi (dR < 0.005, |dpt/pt| < 1e-3); every row must match.
The per-row flag is cached in MVA_ROOT/aux/ so step 2 can be re-run alone.
Step 2: same layout as ds_mass_ssdata_vs_osdata_DY.png, one row per muon origin.
Output: os_ss/dy_ds_mass_mu_origin_ssdata_vs_<KEY_OS>.png (shared dir).
"""
import glob
import os
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
import common  # noqa: E402

BASE = "ssdata"
KEY_OS = sys.argv[1] if len(sys.argv) > 1 else "osdata"
EFFS = (0.90, 0.70, 0.50)
BIN_MEV = {0.90: 20, 0.70: 30, 0.50: 40}
M_DS, M_D, SR_HALF = 1.96835, 1.8696, 0.04
NANO = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
ERA_DIR = {1: "2022", 2: "2022EE"}   # DY_merged `era` branch -> nano directory
CACHE = f"{common.MVA_ROOT}/aux/DY_merged_mu_genPartFlav.npz"
ORIGINS = [("muon from tau (flav 15)", [15]), ("prompt muon (flav 1)", [1]),
           ("other (flav 0,3,4,5)", [0, 3, 4, 5])]
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "era", "jet_bin", "run", "luminosityBlock",
                           "event", "mu_pt", "mu_eta", "mu_phi", "mu_charge", "pi_charge",
                           *common.SELECTION]))


def key(lumi, evt):
    return (lumi.astype(np.int64) << 32) | evt.astype(np.int64)


def read_file(args):
    path, wanted = args
    with uproot.open(path)["Events"] as t:
        a = t.arrays(["run", "luminosityBlock", "event", "Muon_pt", "Muon_eta", "Muon_phi",
                      "Muon_genPartFlav"], library="np")
    k = key(a["luminosityBlock"], a["event"])
    keep = np.isin(k, wanted)
    return (a["run"][keep], k[keep], a["Muon_pt"][keep], a["Muon_eta"][keep], a["Muon_phi"][keep],
            a["Muon_genPartFlav"][keep])


arr = {}
for k in (BASE, KEY_OS):
    with uproot.open(f"{common.out_dir_for(k, 'DNN')}/Prediction/DY_merged.root")[common.TREE] as t:
        arr[k] = t.arrays(WANT, library="np")
a0 = arr[BASE]
assert np.array_equal(arr[KEY_OS]["Ds_mass"], a0["Ds_mass"]), "DY_merged: row order differs"
common.check_charge_config(a0, False, "DY_merged")
assert np.all(a0["run"] == 1), "MC run number is not 1 -- the join key would be ambiguous"
rowkey = key(a0["luminosityBlock"], a0["event"])

# ---------------- step 1: muon origin per row (cached) ----------------
if os.path.exists(CACHE):
    c = np.load(CACHE)
    assert np.array_equal(c["rowkey"], rowkey), f"{CACHE} is stale (row keys differ); delete it"
    flav = c["flav"]
    print(f"loaded {CACHE}")
else:
    flav = np.full(len(rowkey), -1, dtype=np.int16)
    for era, edir in ERA_DIR.items():
        for jb in (0, 1, 2):
            rows = np.where((a0["era"] == era) & (a0["jet_bin"] == jb))[0]
            files = sorted(glob.glob(f"{NANO}/{edir}/DYto2L2Jets_{jb}J_H_phipimu_GenTag/output_*.root"))
            wanted = np.unique(rowkey[rows])
            with Pool(8) as pool:
                parts = pool.map(read_file, [(f, wanted) for f in files], chunksize=4)
            ev_key = np.concatenate([p[1] for p in parts])
            assert np.all(np.concatenate([p[0] for p in parts]) == 1)
            assert len(np.unique(ev_key)) == len(ev_key), f"{edir} {jb}J: duplicate event keys in nano"
            mpt, meta, mphi, mfl = (np.concatenate([p[i] for p in parts]) for i in (2, 3, 4, 5))
            idx = {kk: i for i, kk in enumerate(ev_key)}
            missing = 0
            for r in rows:
                i = idx.get(rowkey[r])
                if i is None:
                    missing += 1
                    continue
                dphi = np.abs((mphi[i] - a0["mu_phi"][r] + np.pi) % (2 * np.pi) - np.pi)
                dr = np.hypot(meta[i] - a0["mu_eta"][r], dphi)
                j = int(np.argmin(dr)) if len(dr) else -1
                if j >= 0 and dr[j] < 0.005 and abs(mpt[i][j] / a0["mu_pt"][r] - 1) < 1e-3:
                    flav[r] = mfl[i][j]
            print(f"  {edir} {jb}J: {len(files)} files, {len(rows):,} rows, events missing in nano "
                  f"{missing}, rows without a matched muon {(flav[rows] < 0).sum()}")
    assert np.all(flav >= 0), f"{(flav < 0).sum()} rows without a muon match -- refusing to cache"
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez(CACHE, rowkey=rowkey, flav=flav)
    print(f"cached {CACHE}")

# ---------------- step 2: SS vs OS DNN per muon origin ----------------
thr = {k: {e: common.stage1_threshold({"key": k, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
       for k in (BASE, KEY_OS)}
sel = common.selection_mask(a0)
mass, w, fl = a0["Ds_mass"][sel], a0["weight"][sel].astype(float), flav[sel]
pass_b = {e: arr[BASE]["score_dsmu"][sel] > thr[BASE][e] for e in EFFS}
pass_o = {e: arr[KEY_OS]["score_dsmu"][sel] > thr[KEY_OS][e] for e in EFFS}
print("\nmuon origin share of DY (weighted):")
for lab, keep in [("no cut", np.ones(len(w), bool))] + [
        (f"{k} eps{100 * e:.0f}", (pass_b if k == BASE else pass_o)[e]) for k in (BASE, KEY_OS) for e in EFFS]:
    print(f"  {lab:14s} " + "  ".join(f"{o.split(' (')[0]} {w[keep & np.isin(fl, v)].sum() / w[keep].sum():6.1%}"
                                       for o, v in ORIGINS))


def h(x, ww, bins):
    return np.histogram(x, bins=bins, weights=ww)[0], np.histogram(x, bins=bins, weights=ww ** 2)[0]


plt.style.use(hep.style.CMS)
nrow = len(ORIGINS)
fig = plt.figure(figsize=(8.5 * len(EFFS), 8.0 * nrow))
outer = gridspec.GridSpec(nrow, len(EFFS), figure=fig, hspace=0.32, wspace=0.25, top=1 - 0.9 / (8.0 * nrow))
print()
for row, (olab, vals) in enumerate(ORIGINS):
    g = np.isin(fl, vals)
    for col, e in enumerate(EFFS):
        bins = np.linspace(M_DS - 0.3, M_DS + 0.3, int(round(600 / BIN_MEV[e])) + 1)
        bc = 0.5 * (bins[1:] + bins[:-1])
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=[3, 1],
                                                 hspace=0.05)
        ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
        pb, po = g & pass_b[e], g & pass_o[e]
        hb, vb = h(mass[pb], w[pb], bins)
        ho, vo = h(mass[po], w[po], bins)
        hep.histplot(hb, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vb)), color="steelblue",
                     marker="o", markersize=5, elinewidth=1.3, label=f"{BASE} (SS DNN) cut  N={hb.sum():,.0f}")
        hep.histplot(ho, bins, ax=ax, histtype="errorbar", yerr=np.sqrt(np.abs(vo)), color="black",
                     marker="s", markersize=5, elinewidth=1.3, label=f"{KEY_OS} (OS DNN) cut  N={ho.sum():,.0f}")
        ax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15, zorder=0)
        ax.axvline(M_D, color="green", ls=":", lw=1.5)
        ax.set_xlim(bins[0], bins[-1])
        top = max(hb.max(), ho.max())
        ax.set_ylim(min(0.0, ax.get_ylim()[0]), 1.4 * top if top > 0 else 1.0)
        ax.set_ylabel(f"Events / {BIN_MEV[e]} MeV")
        ax.set_title(rf"DY, {olab}, $\varepsilon_{{sig}}$={100 * e:.0f}%", fontsize=15)
        ax.legend(fontsize=11, loc="upper right")
        ax.tick_params(labelbottom=False)
        # per-bin ratio, overlap-aware: (C+A)/(C+B) with C both pass, A only OS, B only SS
        c, vc = h(mass[pb & po], w[pb & po], bins)
        a_, va = h(mass[po & ~pb], w[po & ~pb], bins)
        b, vbb = h(mass[pb & ~po], w[pb & ~po], bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (c + a_) / (c + b)
            er = np.sqrt(vc * ((b - a_) / (c + b) ** 2) ** 2 + va / (c + b) ** 2
                         + vbb * ((c + a_) / (c + b) ** 2) ** 2)
        rax.errorbar(bc, r, yerr=er, fmt="o", ms=4, color="black", elinewidth=1.1)
        rax.axhline(1.0, color="red", lw=1.0)
        rax.axvspan(M_DS - SR_HALF, M_DS + SR_HALF, color="grey", alpha=0.15)
        rax.set_ylim(0.0, 1.99)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.set_xlim(bins[0], bins[-1])
        rax.set_xlabel(r"$M(K^+K^-\pi^\pm)$ [GeV]")
        rax.set_ylabel("OS / SS")
        with np.errstate(divide="ignore", invalid="ignore"):
            print(f"  {olab:26s} eps{100 * e:.0f}: SS-DNN {hb.sum():8,.1f}  OS-DNN {ho.sum():8,.1f}"
                  f"  OS/SS {ho.sum() / hb.sum():.3f}")
hep.cms.label(ax=fig.axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=11)
fig.suptitle(rf"DY$\to\ell\ell$ (0J+1J+2J), OS: Ds_mass by best-muon origin after SS-trained ({BASE}) "
             f"vs OS-trained ({KEY_OS}) stage-1 cut", fontsize=18, y=0.995)
out = common.fig_path("os_ss", f"dy_ds_mass_mu_origin_{BASE}_vs_{KEY_OS}.png")
plt.savefig(out, dpi=90, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

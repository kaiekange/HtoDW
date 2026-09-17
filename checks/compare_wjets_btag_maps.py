#!/usr/bin/env python3
"""Is the pooled W+jets b-tag efficiency map justified? Compare jet bins cell by cell.

Usage: compare_wjets_btag_maps.py [era]      (default 2022)
Reads  btag_efficiency/WtoLNu2Jets_{0J,1J,2J}_<era>_counts.json (written by
build_btag_efficiency_map.py) and, for each jet bin j and flavour, compares its efficiency
in every (|eta|, pt) cell with the efficiency of the OTHER two jet bins combined -- two
independent samples, unlike j versus the pooled map, which contains j. Binomial errors;
cells with fewer than MIN_JETS jets on either side are skipped.
Also checks that the three per-bin counts add up exactly to the pooled map's counts.
"""
import json
import sys

import numpy as np

OUT_DIR = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/btag_efficiency"
MIN_JETS = 10
FLAV = {"0": "light", "4": "c", "5": "b"}


def load(name):
    with open(f"{OUT_DIR}/{name}_counts.json") as f:
        c = json.load(f)
    return {k: np.array(v) for k, v in c["total"].items()}, {k: np.array(v) for k, v in c["passed"].items()}, c


def eff_err(p, t):
    e = np.divide(p, t, out=np.zeros_like(p, dtype=float), where=t > 0)
    err = np.sqrt(np.divide(e * (1 - e), t, out=np.zeros_like(e), where=t > 0))
    return e, err


def main():
    era = sys.argv[1] if len(sys.argv) > 1 else "2022"
    bins = {jb: load(f"WtoLNu2Jets_{jb}_{era}") for jb in ("0J", "1J", "2J")}
    pooled_t, pooled_p, meta = load(f"WtoLNu2Jets_{era}")

    for fl in FLAV:
        st = sum(bins[jb][0][fl] for jb in bins)
        sp = sum(bins[jb][1][fl] for jb in bins)
        same = np.array_equal(st, pooled_t[fl]) and np.array_equal(sp, pooled_p[fl])
        print(f"[closure] {FLAV[fl]:5s}: sum of per-bin counts == pooled counts: {same}")

    print(f"\nabseta edges {meta['abseta_edges']}, pt edges {meta['pt_edges']}; cells need >= {MIN_JETS} jets on both sides")
    for fl in FLAV:
        print(f"\n=== {FLAV[fl]} (flav {fl}) ===")
        pe, _ = eff_err(pooled_p[fl], pooled_t[fl])
        print(f"  pooled: jets {int(pooled_t[fl].sum())}, tagged {int(pooled_p[fl].sum())},"
              f" overall eff {pooled_p[fl].sum() / max(pooled_t[fl].sum(), 1):.4f}")
        for jb in bins:
            t, p = bins[jb][0][fl], bins[jb][1][fl]
            rt = sum(bins[o][0][fl] for o in bins if o != jb)
            rp = sum(bins[o][1][fl] for o in bins if o != jb)
            e, s = eff_err(p, t)
            re, rs = eff_err(rp, rt)
            ok = (t >= MIN_JETS) & (rt >= MIN_JETS) & ((s > 0) | (rs > 0))
            pull = np.where(ok, (e - re) / np.sqrt(np.where(ok, s ** 2 + rs ** 2, 1.0)), 0.0)
            chi2, ndf = float(np.sum(pull[ok] ** 2)), int(ok.sum())
            worst = np.unravel_index(np.argmax(np.abs(np.where(ok, pull, 0))), pull.shape) if ndf else None
            ratio = (p.sum() / max(t.sum(), 1)) / max(rp.sum() / max(rt.sum(), 1), 1e-12)
            line = (f"  {jb}: jets {int(t.sum()):8d}  overall eff {p.sum() / max(t.sum(), 1):.4f}"
                    f" (others {rp.sum() / max(rt.sum(), 1):.4f}, ratio {ratio:.3f})"
                    f"  cells compared {ndf:2d}  chi2/ndf {chi2:7.2f}/{ndf}")
            if ndf:
                i, k = worst
                line += (f"  largest |pull| {abs(pull[i, k]):.2f} at abseta bin {i} pt bin {k}"
                         f" ({e[i, k]:.4f}+-{s[i, k]:.4f} vs {re[i, k]:.4f}+-{rs[i, k]:.4f})")
            print(line)
        # cells actually used by the corrections: per-cell pooled efficiency for reference
        print("  pooled eff per cell [abseta rows x pt columns]:")
        for i in range(pe.shape[0]):
            print("    " + " ".join(f"{x:.4f}" for x in pe[i]))


if __name__ == "__main__":
    main()

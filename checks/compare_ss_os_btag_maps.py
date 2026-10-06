#!/usr/bin/env python3
"""Does the same-sign-pi (SS) W+jets sample need its own b-tag efficiency map? Compare it with OS.

Usage: compare_ss_os_btag_maps.py <era> <jet bin> [<jet bin> ...]    e.g. 2022 0J 2J
Reads  btag_efficiency/WtoLNu2Jets_<jb>_<era>_counts.json (OS) and
       btag_efficiency/WtoLNu2Jets_<jb>_SS_<era>_counts.json (SS)
and, per flavour, compares the efficiency in every (|eta|, pt) cell. Binomial errors, treated
as independent -- they are not quite: an event can carry both an OS and an SS candidate, so
the two maps share some jets. Independent errors overestimate sigma, so the pulls are
conservative (too small); a small chi2 alone does not prove the maps agree.
Cells with fewer than MIN_JETS jets on either side are skipped.
"""
import sys

import numpy as np

from compare_wjets_btag_maps import FLAV, MIN_JETS, eff_err, load


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    era = sys.argv[1]
    for jb in sys.argv[2:]:
        os_t, os_p, meta = load(f"WtoLNu2Jets_{jb}_{era}")
        ss_t, ss_p, _ = load(f"WtoLNu2Jets_{jb}_SS_{era}")
        print(f"\n##### {jb} {era}: abseta edges {meta['abseta_edges']}, pt edges {meta['pt_edges']};"
              f" cells need >= {MIN_JETS} jets on both sides")
        for fl in FLAV:
            e_os, s_os = eff_err(os_p[fl], os_t[fl])
            e_ss, s_ss = eff_err(ss_p[fl], ss_t[fl])
            ok = (os_t[fl] >= MIN_JETS) & (ss_t[fl] >= MIN_JETS) & ((s_os > 0) | (s_ss > 0))
            pull = np.where(ok, (e_ss - e_os) / np.sqrt(np.where(ok, s_os ** 2 + s_ss ** 2, 1.0)), 0.0)
            chi2, ndf = float(np.sum(pull[ok] ** 2)), int(ok.sum())
            o_all = os_p[fl].sum() / max(os_t[fl].sum(), 1)
            s_all = ss_p[fl].sum() / max(ss_t[fl].sum(), 1)
            line = (f"  {FLAV[fl]:5s}: jets OS {int(os_t[fl].sum()):8d} SS {int(ss_t[fl].sum()):8d}"
                    f"  overall eff OS {o_all:.4f} SS {s_all:.4f} (SS/OS {s_all / max(o_all, 1e-12):.3f})"
                    f"  chi2/ndf {chi2:7.2f}/{ndf}")
            if ndf:
                i, k = np.unravel_index(np.argmax(np.abs(np.where(ok, pull, 0))), pull.shape)
                line += (f"  largest |pull| {abs(pull[i, k]):.2f} at abseta bin {i} pt bin {k}"
                         f" (SS {e_ss[i, k]:.4f}+-{s_ss[i, k]:.4f} vs OS {e_os[i, k]:.4f}+-{s_os[i, k]:.4f})")
            print(line)
            print("    pull per cell [abseta rows x pt columns], '.' = skipped:")
            for i in range(pull.shape[0]):
                print("      " + " ".join(f"{pull[i, k]:+6.2f}" if ok[i, k] else "     ." for k in range(pull.shape[1])))


if __name__ == "__main__":
    main()

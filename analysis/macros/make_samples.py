"""Build every file in MVA/Samples/ from one table: cut, dedup, tag with era/jet_bin, merge.

    python3 make_samples.py                # everything missing; existing outputs skipped
    python3 make_samples.py --only Data    # just one output
    python3 make_samples.py --dry-run      # print the plan, touch nothing
    python3 make_samples.py --force        # rebuild even if the output exists

WHY THIS EXISTS. Until 2026-09-04 the 28 files in Samples/ were produced by hand, one
highest_pt_dedup.py command each -- `grep -rn highest_pt_dedup` found eight hits, all of
them comments, and zero call sites. The recipe for each sample survived only as the
content of the resulting file. This table IS the recipe now.

WHAT CHANGED IN THE OUTPUT. era and jet bin used to be encoded in the FILENAME
(Data_2022 vs Data_2022EE, WtoLNu2Jets_0J_2022 vs _1J_...). They are now BRANCHES, so
28 files collapse to 9 and a new era is one more `piece` entry rather than a new file
plus edits in common.py, both predict.py's, and ~18 macros.

THE THREE AXES THAT ARE **NOT** BRANCHES, and why:

  * OS vs SS -- these are two different candidate reconstructions, not a cut. The
    tuplizer flips the required pi charge (H_phipimu.cc:744-747), so the same event can
    appear in both files with a DIFFERENT Ds; measured overlap is 12.9% of events. A
    merged file would break the one-row-per-event invariant this script exists to
    establish. Worse, num_os_mu_tight/num_ss_mu_tight are defined relative to the pi
    (H_phipimu.cc:749-759), so the same muon counts as "os" in an OS file and "ss" in an
    SS one -- merged, those two branches become meaningless.
  * incl vs the 0J/1J/2J bins -- ALTERNATIVE descriptions of one physical process, each
    with its own sumGenWeights_ and cross section (63990 = 51846+8860+3284 exactly).
    Never components. A jet_bin value of "incl" would invite exactly the double-count.
  * GenTag -- the SAME events as 0J+1J+2J, re-tuplized with truth branches. Summing it
    with them double-counts. It also carries a different branch set (+30 truth, -14
    match), so it cannot share a file with them anyway.

Lumi_reweight MEANS THREE DIFFERENT THINGS depending on weight_mode, which is why the
Meta tree records the mode per file rather than leaving it to a comment:
    data        1.0                        -- the observation carries no weight
    xsec        gen_weight * sigma*L/sum(w) -- absolute expected yield; summable
    signal_era  SIGNAL_ERA_WEIGHT[era]      -- RELATIVE only; summing gives no yield
"""

import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eras import ERAS, SIGNAL_ERA_WEIGHT

from highest_pt_dedup import (apply_row_cuts, data_weight, dedup_best_indices,
                              sum_gen_weights_of, xsec_weight)

BASE = "/eos/home-k/kakang/H2DW/Analysis"
OUTDIR = f"{BASE}/MVA/Samples"
DEFAULT_TREE = "Hphipimu/Events"

# W+jets cross sections in pb, reweighted 2026-08-25. The bins partition the inclusive
# sample exactly under these values: 63990 = 51846 + 8860 + 3284.
#
# DO NOT APPLY AN EXTRA W+jets SCALE FACTOR ON TOP OF THESE. The 0.929803403 factor
# that has been passed around is the rescaling that PRODUCED this set from the earlier
# one -- it is already inside these numbers, and multiplying by it again double-counts.
# (Separately, a provisional SF = 20.79/21.03 was carried on 2026-08-25 and removed the
# same day, never traced to a source; see compare_incl_by_era.py's header.)
#
# So Lumi_reweight has exactly three inputs, all recorded in the Meta tree:
#     Lumi_reweight_i = gen_weight_i * (sigma * L / sum(gen_weight))
XSEC_JETBIN = {0: 51846.0, 1: 8860.0, 2: 3284.0}
XSEC_INCL = 63990.0


def _wj(era, jetbin, ss):
    tag = ERAS[era]["tag"]
    sfx = "_SS" if ss else ""
    return (f"{BASE}/MC/{tag}/WtoLNu2Jets_{jetbin}J_H_phipimu{sfx}"
            f"/WtoLNu2Jets_{jetbin}J_{tag}{sfx}.root")


def _wj_incl(era, ss):
    tag = ERAS[era]["tag"]
    sfx = "_SS" if ss else ""
    return (f"{BASE}/MC/{tag}/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"
            f"/WtoLNu2Jets_{tag}_H_phipimu{sfx}/WtoLNu2Jets_{tag}{sfx}.root")


def _wj_bin_pieces(ss):
    return [dict(era=e, jet_bin=j, xsec=XSEC_JETBIN[j], path=_wj(e, j, ss))
            for e in (1, 2) for j in (0, 1, 2)]


# Every output file, and the input pieces it is built from. `tree` defaults to
# Hphipimu/Events; the two preEE signals predate the 2026-08-18 module-label
# unification and still live under HphipimuDsW / HphipimuDsstW.
SAMPLES = {
    "HtoDsW": dict(
        weight_mode="signal_era", cuts=("veto_pass", "os", "b_veto"),
        pieces=[
            dict(era=1, tree="HphipimuDsW/Events",
                 path=f"{BASE}/MC/2022/Signal_H_phipimu/HtoDsW_2022.root"),
            dict(era=2, path=f"{BASE}/MC/2022EE/Signal_H_phipimu/HtoDsW_2022EE.root"),
        ]),
    "HtoDsstW": dict(
        weight_mode="signal_era", cuts=("veto_pass", "os", "b_veto"),
        pieces=[
            dict(era=1, tree="HphipimuDsstW/Events",
                 path=f"{BASE}/MC/2022/Signal_DsstW_H_phipimu/HtoDsstW_2022.root"),
            dict(era=2, path=f"{BASE}/MC/2022EE/Signal_DsstW_H_phipimu/HtoDsstW_2022EE.root"),
        ]),
    # OS data is blinded in the Ds signal region; SS is not, because a same-sign pi
    # excludes real signal by charge conservation, so there is nothing to hide.
    "Data": dict(
        weight_mode="data", cuts=("veto_pass", "os", "b_veto", "blind"),
        pieces=[dict(era=e, path=f"{BASE}/Data/{ERAS[e]['tag']}/Data_{ERAS[e]['tag']}.root")
                for e in (1, 2)]),
    "Data_SS": dict(
        weight_mode="data", cuts=("veto_pass", "ss", "b_veto"),
        pieces=[dict(era=e, path=f"{BASE}/Data/{ERAS[e]['tag']}/Data_{ERAS[e]['tag']}_SS.root")
                for e in (1, 2)]),
    "WtoLNu2Jets": dict(
        weight_mode="xsec", cuts=("veto_pass", "os", "b_veto"),
        pieces=_wj_bin_pieces(ss=False)),
    "WtoLNu2Jets_SS": dict(
        weight_mode="xsec", cuts=("veto_pass", "ss", "b_veto"),
        pieces=_wj_bin_pieces(ss=True)),
    "WtoLNu2Jets_incl": dict(
        weight_mode="xsec", cuts=("veto_pass", "os", "b_veto"),
        pieces=[dict(era=e, xsec=XSEC_INCL, path=_wj_incl(e, ss=False)) for e in (1, 2)]),
    "WtoLNu2Jets_incl_SS": dict(
        weight_mode="xsec", cuts=("veto_pass", "ss", "b_veto"),
        pieces=[dict(era=e, xsec=XSEC_INCL, path=_wj_incl(e, ss=True)) for e in (1, 2)]),
    # preEE only: the postEE GenTag jet bins are produced on disk but not yet hadd'ed
    # (0J and 1J have output_*.root but no merged file as of 2026-09-04).
    "WtoLNu2Jets_GenTag": dict(
        weight_mode="xsec", cuts=("veto_pass", "os", "b_veto"),
        pieces=[dict(era=1, jet_bin=j, xsec=XSEC_JETBIN[j],
                     path=f"{BASE}/MC/2022/WtoLNu2Jets_{j}J_H_phipimu_GenTag"
                          f"/WtoLNu2Jets_{j}J_GenTag_2022.root")
                for j in (0, 1, 2)]),
}


def _assert_matches_common():
    """The producer's table and the consumer's list must name the same files.

    Imported lazily and tolerantly: common.py pulls in torch, which this script does
    not need, so a missing torch must not stop a rebuild -- it only skips the check.
    """
    try:
        sys.path.insert(0, os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
        import common
    except ImportError as e:
        print(f"(skipping the common.ALL_SAMPLES cross-check: {e})")
        return
    if sorted(common.ALL_SAMPLES) != sorted(SAMPLES):
        sys.exit("SAMPLES here and common.ALL_SAMPLES disagree:\n"
                 f"  only here          : {sorted(set(SAMPLES) - set(common.ALL_SAMPLES))}\n"
                 f"  only in common.py  : {sorted(set(common.ALL_SAMPLES) - set(SAMPLES))}")


def build(name, spec, dry_run=False):
    out_path = f"{OUTDIR}/{name}.root"
    has_jet_bin = any("jet_bin" in p for p in spec["pieces"])
    cuts = dict(veto_pass="veto_pass" in spec["cuts"],
                blind="blind" in spec["cuts"],
                b_veto="b_veto" in spec["cuts"],
                lepton="ss" if "ss" in spec["cuts"] else ("os" if "os" in spec["cuts"] else None))

    print(f"\n{'='*78}\n{name}  ->  {out_path}")
    print(f"  weight_mode={spec['weight_mode']}  cuts={'+'.join(spec['cuts'])}"
          f"  jet_bin branch={'yes' if has_jet_bin else 'no'}")
    for p in spec["pieces"]:
        jb = f" jet_bin={p['jet_bin']}" if "jet_bin" in p else ""
        xs = f" xsec={p['xsec']:g}" if "xsec" in p else ""
        missing = "" if os.path.exists(p["path"]) else "   *** MISSING ***"
        print(f"    era={p['era']}{jb}{xs}  {p['path']}{missing}")
    if dry_run:
        return

    missing = [p["path"] for p in spec["pieces"] if not os.path.exists(p["path"])]
    if missing:
        sys.exit(f"{name}: input(s) missing:\n  " + "\n  ".join(missing))

    tmp_path = out_path + ".tmp"
    meta = {k: [] for k in ("era", "jet_bin", "sum_gen_weights", "n_rows_in",
                            "n_rows_out", "xsec", "lumi_pb")}
    ref_types = None
    total_out = 0

    # One piece at a time: read, cut, weight, dedup, tag, extend. Never two pieces in
    # memory at once -- Data_2022EE alone is 3.8 GB on disk and many times that as numpy,
    # so holding both eras would not fit in this machine's RAM.
    with uproot.recreate(tmp_path) as f_out:
        for piece in spec["pieces"]:
            era, path = piece["era"], piece["path"]
            tree = piece.get("tree", DEFAULT_TREE)
            print(f"\n  -- era {era}"
                  + (f" jet_bin {piece['jet_bin']}" if "jet_bin" in piece else "")
                  + f": {os.path.basename(path)} [{tree}]")
            with uproot.open(path)[tree] as t:
                data = t.arrays(library="np")
            n_in = len(data["event_id"])
            print(f"     rows in: {n_in:,}")

            apply_row_cuts(data, **cuts)
            n_cut = len(data["event_id"])

            lumi = ERAS[era]["lumi_pb"]
            sgw = float("nan")
            if spec["weight_mode"] == "data":
                w = data_weight(data, n_cut)
            elif spec["weight_mode"] == "xsec":
                sgw = sum_gen_weights_of(path, tree)
                w = xsec_weight(data, piece["xsec"], lumi, sgw)
            elif spec["weight_mode"] == "signal_era":
                w = np.full(n_cut, SIGNAL_ERA_WEIGHT[era], dtype=np.float64)
                print(f"     signal era weight (RELATIVE): {SIGNAL_ERA_WEIGHT[era]:.6f}")
            else:
                sys.exit(f"unknown weight_mode {spec['weight_mode']!r}")
            data["Lumi_reweight"] = w

            best = dedup_best_indices(data)
            print(f"     dedup: {n_cut:,} -> {best.size:,}")
            out = {k: v[best] for k, v in data.items()}
            del data

            out["era"] = np.full(best.size, era, dtype=np.int32)
            if has_jet_bin:
                if "jet_bin" not in piece:
                    sys.exit(f"{name}: piece era={era} has no jet_bin but others do")
                out["jet_bin"] = np.full(best.size, piece["jet_bin"], dtype=np.int32)

            types = {k: v.dtype for k, v in out.items()}
            if ref_types is None:
                ref_types = types
                f_out.mktree("Events", ref_types)
            elif types != ref_types:
                only_a = set(types) - set(ref_types)
                only_b = set(ref_types) - set(types)
                sys.exit(f"{name}: branch set differs between pieces.\n"
                         f"  only in this piece: {sorted(only_a)}\n"
                         f"  only in earlier   : {sorted(only_b)}")
            f_out["Events"].extend(out)

            meta["era"].append(era)
            meta["jet_bin"].append(piece.get("jet_bin", -1))
            meta["sum_gen_weights"].append(sgw)
            meta["n_rows_in"].append(n_in)
            meta["n_rows_out"].append(best.size)
            meta["xsec"].append(piece.get("xsec", float("nan")))
            meta["lumi_pb"].append(lumi)
            total_out += best.size
            del out

        # Meta restores what dedup drops: highest_pt_dedup writes only Events, so
        # sumGenWeights_ -- the normalisation denominator -- is otherwise unrecoverable
        # from a Samples/ file and Lumi_reweight can never be re-derived or rescaled.
        f_out.mktree("Meta", {"era": np.int32, "jet_bin": np.int32,
                              "sum_gen_weights": np.float64, "n_rows_in": np.int64,
                              "n_rows_out": np.int64, "xsec": np.float64,
                              "lumi_pb": np.float64})
        f_out["Meta"].extend({k: np.asarray(v) for k, v in meta.items()})

    os.replace(tmp_path, out_path)
    print(f"\n  wrote {total_out:,} rows -> {out_path} "
          f"({os.path.getsize(out_path)/1e9:.2f} GB)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only", help="build just this output (a key of SAMPLES)")
    p.add_argument("--dry-run", action="store_true", help="print the plan, touch nothing")
    p.add_argument("--force", action="store_true", help="rebuild even if the output exists")
    args = p.parse_args()

    _assert_matches_common()

    todo = list(SAMPLES)
    if args.only:
        if args.only not in SAMPLES:
            sys.exit(f"unknown sample {args.only!r}; known: {', '.join(SAMPLES)}")
        todo = [args.only]

    os.makedirs(OUTDIR, exist_ok=True)
    for name in todo:
        out_path = f"{OUTDIR}/{name}.root"
        if os.path.exists(out_path) and not args.force and not args.dry_run:
            print(f"{name}: exists, skipping (use --force to rebuild)")
            continue
        build(name, SAMPLES[name], dry_run=args.dry_run)
    print("\nDone.")


if __name__ == "__main__":
    main()

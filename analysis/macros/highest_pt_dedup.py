"""Dedup an H_phipimu (OS or sameSignPi) flat ntuple to one row per event: keep the
candidate with the highest Ds_pt, breaking ties on gamma_pt (rows are per-(Ds,gamma)
candidate: Ds_pt is identical across every gamma row of the same Ds, so when the
winning Ds has several gamma rows they all tie on Ds_pt -- gamma_pt picks the best one).

Works on any H_phipimu(_SS) output, signal or background, OS or SS -- just point
it at the right TDirectory/tree path. Every cfg uses the module label Hphipimu since
2026-08-18, so that path is "Hphipimu/Events" everywhere; tuples produced before then have
the signal trees under HphipimuDsW / HphipimuDsstW instead.

Also writes a Lumi_reweight branch. For MC it is the per-event luminosity weight

    Lumi_reweight_i = gen_weight_i * sigma * L / sum(gen_weight over all generated)

The denominator is sumGenWeights_ from the Info tree, summed over its entries (hadd
concatenates Info rather than adding, so one row per merged input file). The per-event
gen_weight factor is what carries the SIGN: amcatnloFXFX samples are ~15.6% negative
weight, and dropping the factor would count those events as positive -- the exact defect
the gen_weight branch was added to fix. Summing the output weight over selected events
gives sigma * L * (sum_sel / sum_all), i.e. the expected yield.

For data --data writes Lumi_reweight = 1.0: data is the observation and carries no
luminosity weight.

Mode is explicit, never inferred: forgetting --xsec on an MC sample would otherwise
silently produce an unnormalised sample, which is invisible on a plot.

This file is BOTH a CLI and a library. make_samples.py imports apply_row_cuts(),
dedup_best_indices() and xsec_weight() to build the merged per-era Samples/ files;
the CLI below stays for one-off and debugging runs on a single input.

Usage:
    python3 highest_pt_dedup.py <in.root> <tree_path> <out.root> --data [--veto-pass] [--blind]
    python3 highest_pt_dedup.py <in.root> <tree_path> <out.root> --xsec PB --lumi IPB [--veto-pass]
    python3 highest_pt_dedup.py <in.root> <tree_path> <out.root> --no-lumi-weight [--veto-pass]

Examples:
    python3 highest_pt_dedup.py \\
        /eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_phipimu/HtoDsW.root \\
        Hphipimu/Events best_HtoDsW.root --xsec 0.0123 --lumi 7.9804

    python3 highest_pt_dedup.py Data_2022EE_raw.root Hphipimu/Events Data_2022EE.root --data
"""

import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eras import B_VETO_BRANCH, LEPTON_BRANCHES, M_DS, SR_HALF_WIDTH


# ---------------------------------------------------------------------------
# Library half
# ---------------------------------------------------------------------------
def apply_row_cuts(all_data, *, veto_pass=False, blind=False, lepton=None,
                   b_veto=False, verbose=True):
    """Filter all_data in place. `lepton` is None, "os" or "ss".

    Row cuts that make the output equivalent to the pre-2026-08-17 tuples. Applied
    before the dedup so the "highest-pt Ds" is chosen among rows that survive them.
    Returns the number of surviving rows.
    """
    n_total = len(all_data["event_id"])
    keep = np.ones(n_total, dtype=bool)

    if veto_pass:
        if "jet_vetomap_pass" not in all_data:
            sys.exit("--veto-pass given but there is no jet_vetomap_pass branch; this tuple "
                     "predates 2026-08-17 and already had the veto applied in the tuplizer.")
        keep &= all_data["jet_vetomap_pass"].astype(bool)
        if verbose:
            print(f"  jet_vetomap_pass : {n_total:,} -> {keep.sum():,}")
    if blind:
        keep &= np.abs(all_data["Ds_mass"] - M_DS) >= SR_HALF_WIDTH
        if verbose:
            print(f"  blind |M-Ds|<{SR_HALF_WIDTH}: -> {keep.sum():,}")
    if lepton is not None:
        if lepton not in ("os", "ss"):
            sys.exit(f"lepton must be None, 'os' or 'ss'; got {lepton!r}")
        for b in LEPTON_BRANCHES:
            if b not in all_data:
                sys.exit(f"--{lepton} given but there is no {b} branch.")
        n_os, n_ss = all_data[LEPTON_BRANCHES[0]], all_data[LEPTON_BRANCHES[1]]
        keep &= (n_ss == 1) & (n_os == 0) if lepton == "ss" else (n_os == 1) & (n_ss == 0)
        if verbose:
            print(f"  lepton {lepton.upper()}       : -> {keep.sum():,}")
    if b_veto:
        if B_VETO_BRANCH not in all_data:
            sys.exit(f"--b-veto given but there is no {B_VETO_BRANCH} branch.")
        keep &= all_data[B_VETO_BRANCH] == 0
        if verbose:
            print(f"  {B_VETO_BRANCH}==0: -> {keep.sum():,}")

    if not keep.all():
        # One branch at a time, rebinding in place: a dict comprehension would hold the
        # whole old and new sets simultaneously, which doubles peak RSS on the multi-GB
        # data files.
        for k in list(all_data):
            all_data[k] = all_data[k][keep]
        n_total = len(all_data["event_id"])
        if verbose:
            print(f"Rows after cuts: {n_total}")
    return n_total


def dedup_best_indices(all_data):
    """Sorted row indices of the highest-Ds_pt candidate of each event.

    Group by (run_id, lumi_id, event_id) -- event_id alone is not safe across
    runs/lumis. Ties on Ds_pt are broken on gamma_pt.
    """
    run_id, lumi_id, event_id = all_data["run_id"], all_data["lumi_id"], all_data["event_id"]
    order = np.lexsort((event_id, lumi_id, run_id))
    r, l, e = run_id[order], lumi_id[order], event_id[order]
    group_start = np.ones(len(order), dtype=bool)
    group_start[1:] = (r[1:] != r[:-1]) | (l[1:] != l[:-1]) | (e[1:] != e[:-1])
    starts = np.flatnonzero(group_start)
    ends = np.append(starts[1:], len(order))

    ds_pt = all_data["Ds_pt"]
    gamma_pt = all_data["gamma_pt"]
    best = np.empty(len(starts), dtype=np.int64)
    for i, (s, en) in enumerate(zip(starts, ends)):
        idx = order[s:en]
        pt = ds_pt[idx]
        top = idx[pt == pt.max()]
        best[i] = top[0] if top.size == 1 else top[np.argmax(gamma_pt[top])]
    best.sort()
    return best


def sum_gen_weights_of(input_file, tree_path):
    """sum(Info.sumGenWeights_) for one input file. Exits loudly if unavailable.

    hadd concatenates Info (one row per merged input file), so the rows are summed.
    This number is the normalisation denominator and is NOT recoverable from a
    dedup'd Samples/ file -- highest_pt_dedup writes only the Events tree -- which
    is why make_samples.py records it in a Meta tree.
    """
    info_path = tree_path.rsplit("/", 1)[0] + "/Info" if "/" in tree_path else "Info"
    try:
        with uproot.open(input_file)[info_path] as info:
            if "sumGenWeights_" not in info.keys():
                sys.exit(f"{info_path} has no sumGenWeights_ branch. This tuple predates the "
                         f"2026-08-17 tuplizer and cannot be luminosity-weighted correctly; "
                         f"reproduce it before using --xsec.")
            return float(np.sum(info["sumGenWeights_"].array(library="np")))
    except uproot.KeyInFileError:
        sys.exit(f"No {info_path} tree in {input_file}; cannot get sum(gen_weight).")


def xsec_weight(all_data, xsec, lumi, sum_gen_weights, verbose=True):
    """Per-event Lumi_reweight = gen_weight * sigma * L / sum(gen_weight)."""
    if sum_gen_weights == 0.0:
        sys.exit("sumGenWeights_ sums to zero -- refusing to divide by it.")
    if "gen_weight" not in all_data:
        sys.exit("Events tree has no gen_weight branch. This tuple predates the 2026-08-17 "
                 "tuplizer; reproduce it before using --xsec.")
    norm = xsec * lumi / sum_gen_weights
    if verbose:
        print(f"MC: sigma={xsec:g}  L={lumi:g}  sum(gen_weight)={sum_gen_weights:.6g}")
        print(f"    sigma*L/sum(gen_weight) = {norm:.6g}   "
              f"(per-event weight also carries gen_weight's sign)")
    return all_data["gen_weight"].astype(np.float64) * norm


def data_weight(all_data, n_total, verbose=True):
    """Lumi_reweight = 1.0, with the guard that this really is data."""
    if "gen_weight" in all_data:
        uniq = np.unique(all_data["gen_weight"])
        if not np.array_equal(uniq, np.array([1.0])):
            sys.exit(f"--data given but gen_weight is not all 1.0 (values {uniq[:5]}); "
                     f"this looks like MC, refusing to write an unweighted sample.")
    if verbose:
        print("Data: Lumi_reweight = 1.0")
    return np.ones(n_total, dtype=np.float64)


# ---------------------------------------------------------------------------
# CLI half
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_file")
    p.add_argument("tree_path", help='"Hphipimu/Events" for anything produced since 2026-08-18')
    p.add_argument("out_file")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--data", action="store_true",
                   help="real data: Lumi_reweight = 1.0, no cross section involved")
    g.add_argument("--xsec", type=float,
                   help="MC cross section in pb (same units as --lumi's inverse)")
    g.add_argument("--no-lumi-weight", action="store_true",
                   help="MC, but write no Lumi_reweight branch at all. For samples used at "
                        "weight 1.0 per event, which is how the pre-2026-08 Samples/ files "
                        "were built. Explicit rather than a default so an unnormalised "
                        "background cannot slip through.")
    p.add_argument("--veto-pass", action="store_true",
                   help="keep only jet_vetomap_pass==True rows. Needed for anything produced "
                        "since 2026-08-17, where JvmAppliedEventFilter stopped dropping those "
                        "events; it makes the output equivalent to the older tuples that had "
                        "the veto applied inside the tuplizer. Errors out if the branch is "
                        "absent, so 'I meant to apply it' can never silently mean 'I did not'.")
    p.add_argument("--blind", action="store_true",
                   help="drop rows with |Ds_mass - 1.96835| < 0.04 (the signal region). Applied "
                        "to OS data only; SS structurally excludes signal so it stays open.")
    g2 = p.add_mutually_exclusive_group()
    g2.add_argument("--os", action="store_true",
                    help="keep only num_os_mu_tight==1 & num_ss_mu_tight==0 (the opposite-sign "
                         "lepton configuration: signal, OS data, OS background)")
    g2.add_argument("--ss", action="store_true",
                    help="the same-sign mirror: num_ss_mu_tight==1 & num_os_mu_tight==0")
    p.add_argument("--b-veto", action="store_true",
                   help=f"keep only {B_VETO_BRANCH}==0. Note this branch is per-CANDIDATE (it "
                        "is recomputed from each row's own Ds), so applying it here -- before "
                        "the dedup -- is not identical to applying it afterwards: an event "
                        "whose highest-pt candidate fails can still contribute a lower-pt one. "
                        "Measured on preEE DsW the difference is +23 events (+0.02%%).")
    p.add_argument("--lumi", type=float,
                   help="integrated luminosity in /pb, of the era this MC represents "
                        "(2022 preEE and postEE are different numbers)")
    args = p.parse_args()

    if args.xsec is not None and args.lumi is None:
        p.error("--xsec also needs --lumi")

    out_tree = args.tree_path.split("/")[-1]

    print(f"Reading {args.input_file} [{args.tree_path}]...")
    with uproot.open(args.input_file)[args.tree_path] as t:
        all_data = t.arrays(library="np")

    n_total = len(all_data["event_id"])
    print(f"Total rows: {n_total}")

    n_total = apply_row_cuts(all_data, veto_pass=args.veto_pass, blind=args.blind,
                             lepton="ss" if args.ss else ("os" if args.os else None),
                             b_veto=args.b_veto)

    # Luminosity weight. Resolved before the dedup so a misconfigured run fails
    # immediately rather than after the expensive grouping pass.
    if args.no_lumi_weight:
        per_event_weight = None
        print("MC: no Lumi_reweight branch written (--no-lumi-weight)")
    elif args.data:
        per_event_weight = data_weight(all_data, n_total)
    else:
        sgw = sum_gen_weights_of(args.input_file, args.tree_path)
        per_event_weight = xsec_weight(all_data, args.xsec, args.lumi, sgw)

    if per_event_weight is not None:
        all_data["Lumi_reweight"] = per_event_weight

    best = dedup_best_indices(all_data)
    print(f"Unique events kept: {best.size}  (dropped {n_total - best.size} lower-Ds_pt rows)")

    out_branches = {k: v[best] for k, v in all_data.items()}
    branch_types = {k: v.dtype for k, v in out_branches.items()}

    # Dedup keeps one row per event and gen_weight is per event, so the surviving weights
    # are a subset of the generated ones -- summing them gives the expected yield directly.
    if not args.data and not args.no_lumi_weight:
        print(f"Expected yield after selection: {out_branches['Lumi_reweight'].sum():.4g} events")

    print(f"Writing {args.out_file} [{out_tree}]...")
    with uproot.recreate(args.out_file) as f_out:
        f_out.mktree(out_tree, branch_types)
        f_out[out_tree].extend(out_branches)

    print("Done.")


if __name__ == "__main__":
    main()

"""Does the NanoAOD-based Samples production reproduce the old private-tuplizer
Samples, feature by feature, on the same physical events?

    python3 check_nano_migration_parity.py

Step 1 of the NanoMVA MVA-pipeline migration (see the plan / CLAUDE.md). The old
`analysis/DNN/common.py` FEATURE_NAMES were screened and the models trained on
`/eos/.../Analysis/MVA/Samples/*.root`, produced by the private CMSSW_13_0_13/15_0_5
tuplizer (H_phipimu). The new `/eos/.../Analysis/NanoMVA/Samples/*.root` come from a
different reconstruction path (custom NanoAOD + a separate corrections friend tree,
CMSSW_13_0_13 branch project_nanoaod_migration) that only claims 1e-6-level agreement
with the old tuplizer for the vertex fit (S5 validation) -- NOT the 1e-9 level that
was appropriate for same-code/different-run comparisons (e.g. the era-E condor
regression). Anything worse than ~1e-4 relative here is a real discrepancy to chase,
not floating-point noise.

This does a ROW-LEVEL join on (era[, channel][, jet_bin], run, event) -- the old
tuplizer writes run_id/event_id (no luminosityBlock), the new one writes standard
NanoAOD run/luminosityBlock/event, so luminosityBlock is dropped from the key on
both sides. This was checked to still be unique on both sides for every one of the
five sample pairs before writing this script (all four MVA/Samples + NanoMVA/Samples
productions are already 1 row/event; no dedup ambiguity). Row-level join means the
population-selection mismatch trap (old vs new files carry different upstream
skims) does not apply here -- we only ever compare the SAME physical row on both
sides, never two differently-selected populations.

Usage:
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_nano_migration_parity.py
"""

import sys

import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

OLD_DIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/Samples"
NEW_DIR = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples"

# FEATURE_NAMES as currently declared (pre-migration: 43 features, includes jet_ht
# and num_jets which Step 2 is going to drop -- checking them here anyway confirms
# there really is no usable _corr replacement, closing that question empirically
# rather than by branch-name inspection alone).
FEATURES = common.FEATURE_NAMES

REL_TOL_PASS = 1e-6   # at or below this: identical up to float noise
REL_TOL_FLAG = 1e-4   # above this: real discrepancy, print every time


def pack_key(event, *small_fields):
    """One int64 join key. `small_fields` are each assumed < 100 (era, channel,
    jet_bin, run -- run is always 1 for private MC, small for real data eras)."""
    acc = np.zeros_like(event, dtype=np.int64)
    for f in small_fields:
        acc = acc * 100 + f.astype(np.int64)
    return acc * np.int64(10_000_000_000) + event.astype(np.int64)


def overlap_indices(key_old, key_new):
    """Indices into (old, new) arrays for keys present in both, aligned so
    idx_old[i] and idx_new[i] refer to the same physical event."""
    order_old = np.argsort(key_old, kind="mergesort")
    order_new = np.argsort(key_new, kind="mergesort")
    ko, kn = key_old[order_old], key_new[order_new]
    common_keys, io, in_ = np.intersect1d(ko, kn, assume_unique=True, return_indices=True)
    return order_old[io], order_new[in_], common_keys


def print_feature_table(a_old, a_new, io, in_):
    print(f"  {'feature':32s} {'median|rel diff|':>18s} {'p95':>12s} {'p99':>12s} {'n_finite':>10s}")
    for feat in FEATURES:
        vo = a_old[feat][io].astype(np.float64)
        vn = a_new[feat][in_].astype(np.float64)
        finite = np.isfinite(vo) & np.isfinite(vn)
        denom = np.maximum(np.abs(vo[finite]), 1e-12)
        rel = np.abs(vo[finite] - vn[finite]) / denom
        med = np.median(rel) if len(rel) else np.nan
        p95 = np.percentile(rel, 95) if len(rel) else np.nan
        p99 = np.percentile(rel, 99) if len(rel) else np.nan
        flag = "  <-- FLAG" if p95 > REL_TOL_FLAG else ("  ok(corr?)" if med > REL_TOL_PASS else "")
        print(f"  {feat:32s} {med:18.3e} {p95:12.3e} {p99:12.3e} {finite.sum():10,d}{flag}")


def main():
    # ---- signal: HtoDsW.root / HtoDsstW.root  <->  Signal_merged.root (channel 0/1) ----
    for old_name, channel in (("HtoDsW", 0), ("HtoDsstW", 1)):
        old_path = f"{OLD_DIR}/{old_name}.root"
        new_path = f"{NEW_DIR}/Signal_merged.root"
        a_old = uproot.open(old_path)["Events"].arrays(
            list(dict.fromkeys(FEATURES + ["event_id", "era", "run_id"])), library="np")
        a_new_full = uproot.open(new_path)["Events"].arrays(
            list(dict.fromkeys(FEATURES + ["event", "era", "run", "channel"])), library="np")
        chan_mask = a_new_full["channel"] == channel
        a_new = {k: v[chan_mask] for k, v in a_new_full.items()}

        key_old = pack_key(a_old["event_id"].astype(np.int64), a_old["era"], a_old["run_id"])
        key_new = pack_key(a_new["event"].astype(np.int64), a_new["era"], a_new["run"])
        assert len(np.unique(key_old)) == len(key_old)
        assert len(np.unique(key_new)) == len(key_new)
        io, in_, keys = overlap_indices(key_old, key_new)
        print(f"\n=== signal: {old_name}.root <-> Signal_merged.root (channel={channel}) ===")
        print(f"  n_old={len(key_old):>10,}  n_new(channel filtered)={len(key_new):>10,}  "
              f"n_overlap={len(keys):>10,} ({100*len(keys)/min(len(key_old), len(key_new)):.1f}%)")
        if len(keys) == 0:
            print("  NO OVERLAP")
            continue
        print_feature_table(a_old, a_new, io, in_)

    # ---- Data / Data_SS: (run_id,event_id) is already unique on the old side ----
    for old_name, new_name in (("Data", "Data_OS_merged"), ("Data_SS", "Data_SS_merged")):
        old_path, new_path = f"{OLD_DIR}/{old_name}.root", f"{NEW_DIR}/{new_name}.root"
        a_old = uproot.open(old_path)["Events"].arrays(
            list(dict.fromkeys(FEATURES + ["event_id", "run_id"])), library="np")
        a_new = uproot.open(new_path)["Events"].arrays(
            list(dict.fromkeys(FEATURES + ["event", "run"])), library="np")
        key_old = pack_key(a_old["event_id"].astype(np.int64), a_old["run_id"])
        key_new = pack_key(a_new["event"].astype(np.int64), a_new["run"])
        assert len(np.unique(key_old)) == len(key_old)
        assert len(np.unique(key_new)) == len(key_new)
        io, in_, keys = overlap_indices(key_old, key_new)
        print(f"\n=== {old_name}.root <-> {new_name}.root ===")
        print(f"  n_old={len(key_old):>10,}  n_new={len(key_new):>10,}  "
              f"n_overlap={len(keys):>10,} ({100*len(keys)/min(len(key_old), len(key_new)):.1f}%)")
        if len(keys) == 0:
            print("  NO OVERLAP")
            continue
        print_feature_table(a_old, a_new, io, in_)

    # ---- WtoLNu2Jets: old (run_id,event_id) has 2642 collisions -> key needs era+jet_bin ----
    old_path = f"{OLD_DIR}/WtoLNu2Jets.root"
    new_path = f"{NEW_DIR}/WtoLNu2Jets_merged.root"
    a_old = uproot.open(old_path)["Events"].arrays(
        list(dict.fromkeys(FEATURES + ["event_id", "run_id", "era", "jet_bin"])), library="np")
    a_new = uproot.open(new_path)["Events"].arrays(
        list(dict.fromkeys(FEATURES + ["event", "run", "era", "jet_bin"])), library="np")
    key_old = pack_key(a_old["event_id"].astype(np.int64), a_old["era"], a_old["jet_bin"], a_old["run_id"])
    key_new = pack_key(a_new["event"].astype(np.int64), a_new["era"], a_new["jet_bin"], a_new["run"])
    assert len(np.unique(key_old)) == len(key_old)
    assert len(np.unique(key_new)) == len(key_new)
    io, in_, keys = overlap_indices(key_old, key_new)
    print(f"\n=== WtoLNu2Jets.root <-> WtoLNu2Jets_merged.root ===")
    print(f"  n_old={len(key_old):>10,}  n_new={len(key_new):>10,}  "
          f"n_overlap={len(keys):>10,} ({100*len(keys)/min(len(key_old), len(key_new)):.1f}%)")
    if len(keys) > 0:
        print_feature_table(a_old, a_new, io, in_)


if __name__ == "__main__":
    main()

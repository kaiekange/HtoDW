#!/usr/bin/env python3
"""One-off check: new GenTag MC nano (phi20 inside the producer, 1026-1029 condor)
vs the old nano (<dir>_pre_phi20nano), one sample directory per call.

  check_mc_phi20nano_vs_old.py <year> <dirname> <outdir> [--exclude=J1,J2,...] [--max=N]

The rerun used the same JOBIDs and the same input slicing, so files are paired by
JOBID (output_<J>.root) and events joined by (luminosityBlock, event) inside each pair.
--exclude = JOBIDs dropped from the new sample (bad jobs); they must be missing from
new or are skipped, and must exist in old.  --max is for testing only.

Hard checks (exit 1 if any fails, report is written first):
  JOBID set new == old - excluded; Runs genEventCount/genEventSumw new == old per JOBID;
  every new row inside phi20; new events subset of old; no duplicate keys; run == 1;
  every lost event on the phi20 edge; same events bitwise equal outside EXPECTED_DIFF.
Writes <outdir>/<year>_<dirname>.txt and .json (per-JOBID Runs sums, for OS==SS later).
"""
import json
import sys
import time
from collections import Counter

import awkward as ak
import numpy as np
import uproot
from XRootD import client

SERVER = "root://eosuser.cern.ch"
BASE = "/eos/user/k/kakang/H2DW/Analysis/NanoMVA/MC"
MPHI = 1.019460
PHI_WIN = 0.020
EDGE = 0.01999  # a lost event must have |dM| above this (float32 rounding at the edge)
SKIP_PREFIX = ("HLT_", "L1_", "Flag_", "TrigObj", "SoftActivityJet", "L1Reco", "HLTriggerF")
JAGGED_SUMS = ["Jet_pt", "Muon_pt", "Electron_pt"]
# cand_id = running index over passing candidates, shifts when phi20 removes earlier ones
EXPECTED_DIFF = ("nCandidates", "cand_id", "gamma_uc_dM", "gamma_uc_dR_Ds", "gamma_uc_pt")

out_lines = []
hard_fail = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out_lines.append(s)


def fail(msg):
    hard_fail.append(msg)
    log("FAIL:", msg)


def list_jobids(d):
    st, ls = client.FileSystem(SERVER).dirlist(d)
    assert st.ok, (d, st.message)
    ids = {}
    for e in ls:
        n = e.name
        if n.startswith("output_") and n.endswith(".root"):
            ids[int(n[len("output_"):-len(".root")])] = f"{d}/{n}"
    assert ids, d
    return ids


def read(path, branches):
    for attempt in range(4):
        try:
            with uproot.open(SERVER + "/" + path) as f:
                r = f["Runs"].arrays(["genEventCount", "genEventSumw"], library="np")
                t = f["Events"]
                a = t.arrays(branches, library="np")
                j = t.arrays(JAGGED_SUMS, library="ak")
            break
        except OSError as e:
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1} {path}: {e}", flush=True)
            time.sleep(30)
    for b in JAGGED_SUMS:
        a["sum_" + b] = ak.to_numpy(ak.sum(j[b], axis=1)).astype(np.float64)
    return int(r["genEventCount"].sum()), float(r["genEventSumw"].sum()), a


def flat_branches(path):
    with uproot.open(SERVER + "/" + path) as f:
        t = f["Events"]
        return [k for k, b in t.items() if not k.startswith(SKIP_PREFIX)
                and isinstance(b.interpretation, uproot.interpretation.numerical.AsDtype)]


def key(a):
    lumi, evt = a["luminosityBlock"].astype(np.int64), a["event"].astype(np.int64)
    assert lumi.max(initial=0) < 2**29 and evt.max(initial=0) < 2**34, (lumi.max(), evt.max())
    return (lumi << 34) | evt


def main():
    year, name, outdir = sys.argv[1:4]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[4:])
    excl = {int(x) for x in opts.get("exclude", "").split(",") if x}
    dn, do = f"{BASE}/{year}/{name}", f"{BASE}/{year}/{name}_pre_phi20nano"
    fn, fo = list_jobids(dn), list_jobids(do)
    log(f"=== {year} {name}: new {len(fn)} files, old {len(fo)} files, excluded {sorted(excl)}")

    # --- JOBID sets
    if not excl <= set(fo):
        fail(f"excluded JOBIDs not in old: {sorted(excl - set(fo))}")
    skipped = sorted(excl & set(fn))
    if skipped:
        log(f"[jobid] excluded but present in new (skipped): {skipped}")
    want = set(fo) - excl
    have = set(fn) - excl
    if have != want:
        fail(f"JOBID sets: only new {sorted(have - want)}, only old {sorted(want - have)}")
    jobs = sorted(have & want)
    if "max" in opts:  # testing only
        jobs = jobs[:int(opts["max"])]
    log(f"[jobid] compared {len(jobs)} JOBIDs")

    bn, bo = flat_branches(fn[jobs[0]]), flat_branches(fo[jobs[0]])
    log(f"flat branches: new {len(bn)}, old {len(bo)}; only new {sorted(set(bn) - set(bo))}; "
        f"only old {sorted(set(bo) - set(bn))}")
    br = [b for b in bn if b in set(bo)]
    cmp_br = br + ["sum_" + b for b in JAGGED_SUMS]
    assert "cand_gen_mode" in br and "phi_mass" in br

    runs = {}
    n = Counter()
    mode = {k: Counter() for k in ("old_all", "old_pass", "new", "same", "recovered", "lost")}
    nd = Counter()
    maxd = Counter()
    nanmis = Counter()
    max_dphi_new = 0.0
    lost_list, notold_list = [], []
    for i, J in enumerate(jobs):
        cn, sn, N = read(fn[J], br)
        co, so, O = read(fo[J], br)
        runs[J] = {"new": [cn, sn], "old": [co, so]}
        if cn != co or sn != so:
            fail(f"JOBID {J} Runs new {cn} {sn!r} vs old {co} {so!r}")
        for nm, a in (("new", N), ("old", O)):
            if len(a["run"]) and not (a["run"] == 1).all():
                fail(f"JOBID {J} {nm}: run != 1")
        kn, ko = key(N), key(O)
        if len(np.unique(kn)) != len(kn) or len(np.unique(ko)) != len(ko):
            fail(f"JOBID {J}: duplicate keys new {len(kn) - len(np.unique(kn))} old {len(ko) - len(np.unique(ko))}")
        dmn = np.abs(N["phi_mass"].astype(np.float64) - MPHI)
        dmo = np.abs(O["phi_mass"].astype(np.float64) - MPHI)
        if len(dmn):
            max_dphi_new = max(max_dphi_new, dmn.max())
        n["new_rows"] += len(kn)
        n["old_rows"] += len(ko)
        n["new_rows_ge20"] += int((dmn >= PHI_WIN).sum())
        opass = dmo < PHI_WIN
        n["old_pass"] += int(opass.sum())
        mode["old_all"].update(O["cand_gen_mode"].tolist())
        mode["old_pass"].update(O["cand_gen_mode"][opass].tolist())
        mode["new"].update(N["cand_gen_mode"].tolist())

        both, iN, iO = np.intersect1d(kn, ko, assume_unique=True, return_indices=True)
        notold = np.setdiff1d(kn, ko)
        if len(notold):
            notold_list += [(J, int(k >> 34), int(k & (2**34 - 1))) for k in notold[:5]]
        n["new_not_in_old"] += len(notold)
        s = opass[iO]
        sN, sO, rN = iN[s], iO[s], iN[~s]
        n["same"] += len(sN)
        n["recovered"] += len(rN)
        mode["same"].update(N["cand_gen_mode"][sN].tolist())
        mode["recovered"].update(N["cand_gen_mode"][rN].tolist())
        if len(rN):
            dpt = N["Ds_pt"][rN].astype(np.float64) - O["Ds_pt"][iO[~s]].astype(np.float64)
            n["recovered_pt_not_lower"] += int((dpt >= 0).sum())
        lost = opass & ~np.isin(ko, kn)
        n["lost"] += int(lost.sum())
        mode["lost"].update(O["cand_gen_mode"][lost].tolist())
        for k, d in zip(ko[lost], dmo[lost]):
            lost_list.append((J, int(k >> 34), int(k & (2**34 - 1)), float(d)))

        for b in cmp_br:
            a, c = N[b][sN], O[b][sO]
            if a.dtype.kind == "f":
                eq = (a == c) | (np.isnan(a) & np.isnan(c))
                nanmis[b] += int((np.isnan(a) != np.isnan(c)).sum())
                fin = np.isfinite(a) & np.isfinite(c)
                d = np.abs(a[fin].astype(np.float64) - c[fin].astype(np.float64))
            else:
                eq = a == c
                d = np.abs(a.astype(np.float64) - c.astype(np.float64))
            nd[b] += int((~eq).sum())
            if d.size:
                maxd[b] = max(maxd[b], float(d.max()))
        if i % 50 == 0:
            print(f"  {i}/{len(jobs)} JOBID {J}", flush=True)

    # --- report
    log(f"\n[runs] genEventCount total new {sum(v['new'][0] for v in runs.values())}, "
        f"old {sum(v['old'][0] for v in runs.values())}")
    log(f"[rows] new {n['new_rows']}, old {n['old_rows']}, old passing phi20 {n['old_pass']} "
        f"({n['old_pass'] / max(n['old_rows'], 1):.4f}); new/old-phi20 {n['new_rows'] / max(n['old_pass'], 1):.5f}")
    log(f"[phi20] new rows max |M(KK)-Mphi| = {max_dphi_new * 1e3:.5f} MeV; rows >= 20 MeV: {n['new_rows_ge20']}")
    if max_dphi_new > PHI_WIN + 1e-6:
        fail(f"new row outside phi20: {max_dphi_new * 1e3:.5f} MeV")
    log(f"[subset] new events not in old: {n['new_not_in_old']} (expect 0) {notold_list[:10]}")
    if n["new_not_in_old"]:
        fail("new events not in old")
    log(f"[split] same {n['same']}, recovered {n['recovered']}, lost {n['lost']}; "
        f"anchor same+recovered {n['same'] + n['recovered']} =? new {n['new_rows']}")
    log(f"[split] recovered / old-phi20 = {n['recovered'] / max(n['old_pass'], 1):.5f}; "
        f"recovered with new Ds_pt >= old: {n['recovered_pt_not_lower']} (expect 0)")
    for J, l, e, d in lost_list[:30]:
        log(f"[lost]   JOBID {J} lumi {l} evt {e} |dM| {d * 1e3:.6f} MeV")
    bad_lost = [x for x in lost_list if x[3] < EDGE]
    if bad_lost:
        fail(f"{len(bad_lost)} lost events away from the phi20 edge")

    modes = sorted(set().union(*[set(c) for c in mode.values()]))
    log("\n[mode] rows by cand_gen_mode (old_pass/same/lost: old mode; new/recovered: new mode)")
    log(f"[mode] {'mode':>4} " + " ".join(f"{k:>10}" for k in mode) + f" {'rec/old_pass':>13}")
    for m in modes:
        log(f"[mode] {m:>4} " + " ".join(f"{mode[k][m]:>10}" for k in mode)
            + f" {mode['recovered'][m] / max(mode['old_pass'][m], 1):>13.5f}")

    log(f"\n[same] per-branch comparison on {n['same']} same events (listed only if not all equal)")
    n_eq = 0
    for b in cmp_br:
        if nd[b] == 0:
            n_eq += 1
            continue
        tag = "  (expected)" if b in EXPECTED_DIFF else ""
        log(f"[same]   {b:40s} n_diff {nd[b]:8d} ({nd[b] / max(n['same'], 1):.6f})  max|d| {maxd[b]:.3g}"
            f"  NaN mismatch {nanmis[b]}{tag}")
        if not tag:
            fail(f"branch {b} differs on same events")
    log(f"[same] branches fully bitwise equal: {n_eq} / {len(cmp_br)}")

    log("\nRESULT: " + ("FAIL (" + "; ".join(hard_fail[:5]) + ")" if hard_fail else "PASS"))
    stem = f"{outdir}/{year}_{name}"
    with open(stem + ".txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    with open(stem + ".json", "w") as f:
        json.dump({"runs": runs, "counts": n, "modes": {k: dict(v) for k, v in mode.items()},
                   "excluded": sorted(excl), "hard_fail": hard_fail}, f, indent=1)
    print("DONE", flush=True)
    sys.exit(1 if hard_fail else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""One-off check: new 2022 preEE data OS nano (phi20 inside the producer, 1026 CRAB)
vs the old nano (*_pre_phi20nano, 0915 CRAB).

  check_data_phi20nano_vs_old.py <C|D> <outdir> [max files, testing only]

Events are joined by (run, luminosityBlock, event), never file by file
(unitsPerJob changed, so the file splitting differs).  Writes <outdir>/Run2022<X>.txt
and a small .npz (keys + a few columns) for plotting.  Prints no Ds SR counts
(OS data is blind).
"""
import glob
import sys
import time

import awkward as ak
import numpy as np
import uproot

BASE = "/eos/user/k/kakang/H2DW/Analysis/NanoMVA/Data/2022/Muon"
MPHI = 1.019460
MDS = 1.96835
PHI_WIN = 0.020
SKIP_PREFIX = ("HLT_", "L1_", "Flag_", "PPS", "Proton", "TrigObj", "SoftActivityJet",
               "L1Reco", "HLTriggerF")
JAGGED_SUMS = ["Jet_pt", "Muon_pt", "Electron_pt"]
# cand_id = running index over passing candidates, shifts when phi20 removes earlier ones
EXPECTED_DIFF = ("nCandidates", "cand_id", "gamma_uc_dM", "gamma_uc_dR_Ds", "gamma_uc_pt")

out_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out_lines.append(s)


def files(era, tag):
    d = f"{BASE}/Muon_Run2022{era}_NanoH_phipimu{tag}"
    stamps = glob.glob(f"{d}/*/")
    assert len(stamps) == 1, (d, stamps)
    fs = sorted(glob.glob(f"{d}/*/*/nano_*.root"))
    assert fs, d
    if len(sys.argv) > 3:  # testing only: first N files
        fs = fs[:int(sys.argv[3])]
    return stamps[0], fs


def make_key(run, lumi, evt):
    run = run.astype(np.int64) - 355000
    lumi = lumi.astype(np.int64)
    evt = evt.astype(np.int64)
    assert run.min() >= 0 and run.max() < 2**13, (run.min(), run.max())
    assert lumi.max() < 2**14, lumi.max()
    assert evt.max() < 2**36, evt.max()
    return (run << 50) | (lumi << 36) | evt


def flat_branches(fname):
    t = uproot.open("root://eosuser.cern.ch/" + fname)["Events"]
    out = []
    for k, b in t.items():
        if k.startswith(SKIP_PREFIX):
            continue
        if isinstance(b.interpretation, uproot.interpretation.numerical.AsDtype):
            out.append(k)
    return out


def read(fs, branches):
    lumis, cols = [], {b: [] for b in branches}
    for b in JAGGED_SUMS:
        cols["sum_" + b] = []
    for i, f in enumerate(fs):
        # read via xrootd; EOS gave transient I/O errors (1028), so retry, then fail loud
        for attempt in range(4):
            try:
                with uproot.open("root://eosuser.cern.ch/" + f) as rf:
                    lb = rf["LuminosityBlocks"].arrays(["run", "luminosityBlock"], library="np")
                    t = rf["Events"]
                    a = t.arrays(branches, library="np") if t.num_entries else None
                    j = t.arrays(JAGGED_SUMS, library="ak") if t.num_entries else None
                break
            except OSError as e:
                if attempt == 3:
                    raise
                print(f"  retry {attempt + 1} {f}: {e}", flush=True)
                time.sleep(30)
        lumis.append(make_key(lb["run"], lb["luminosityBlock"], np.zeros(len(lb["run"]))))
        if a is not None:
            for b in branches:
                cols[b].append(a[b])
            for b in JAGGED_SUMS:
                cols["sum_" + b].append(ak.to_numpy(ak.sum(j[b], axis=1)).astype(np.float64))
        if i % 50 == 0:
            print(f"  read {i}/{len(fs)}", flush=True)
    cols = {k: np.concatenate(v) for k, v in cols.items()}
    return np.concatenate(lumis), cols


def main():
    era, outdir = sys.argv[1], sys.argv[2]
    sn, fn = files(era, "")
    so, fo = files(era, "_pre_phi20nano")
    log(f"=== Run2022{era} OS: new {sn} ({len(fn)} files) vs old {so} ({len(fo)} files)")

    bn, bo = flat_branches(fn[0]), flat_branches(fo[0])
    only_n, only_o = sorted(set(bn) - set(bo)), sorted(set(bo) - set(bn))
    log(f"flat branches: new {len(bn)}, old {len(bo)}; only new {only_n}; only old {only_o}")
    br = [b for b in bn if b in set(bo)]

    print("reading new", flush=True)
    ln, N = read(fn, br)
    print("reading old", flush=True)
    lo, O = read(fo, br)

    # --- lumi completeness
    log("\n[lumi] LuminosityBlocks entries new/old:", len(ln), len(lo))
    un, uo = np.unique(ln), np.unique(lo)
    log(f"[lumi] duplicated lumis new {len(ln) - len(un)}, old {len(lo) - len(uo)}")
    common_l = np.intersect1d(un, uo)
    log(f"[lumi] unique new {len(un)}, old {len(uo)}, common {len(common_l)}, "
        f"only new {len(un) - len(common_l)}, only old {len(uo) - len(common_l)}")
    for nm, s in (("only new", np.setdiff1d(un, uo)), ("only old", np.setdiff1d(uo, un))):
        if len(s):
            runs = (s >> 50) + 355000
            r, c = np.unique(runs, return_counts=True)
            log(f"[lumi]   {nm} by run: " + ", ".join(f"{a}:{b}" for a, b in zip(r, c)))

    kn = make_key(N["run"], N["luminosityBlock"], N["event"])
    ko = make_key(O["run"], O["luminosityBlock"], O["event"])
    lk = lambda k: k & ~np.int64(2**36 - 1)
    inl_n, inl_o = np.isin(lk(kn), common_l), np.isin(lk(ko), common_l)
    log(f"\n[events] rows new {len(kn)}, old {len(ko)}; in common lumis new {inl_n.sum()}, old {inl_o.sum()}")
    log(f"[events] duplicate keys new {len(kn) - len(np.unique(kn))}, old {len(ko) - len(np.unique(ko))}")

    # --- phi20 on new rows
    dmn = np.abs(N["phi_mass"].astype(np.float64) - MPHI)
    dmo = np.abs(O["phi_mass"].astype(np.float64) - MPHI)
    log(f"\n[phi20] new rows max |M(KK)-Mphi| = {dmn.max()*1e3:.5f} MeV; rows >=20 MeV: {(dmn >= PHI_WIN).sum()}")
    log(f"[phi20] old rows passing phi20 on the best: {(dmo < PHI_WIN).sum()} / {len(dmo)} = {(dmo < PHI_WIN).mean():.4f}")

    # restrict to common lumis
    selN, selO = np.where(inl_n)[0], np.where(inl_o)[0]
    kN, kO = kn[selN], ko[selO]
    on, oo = np.argsort(kN), np.argsort(kO)
    kN, kO, selN, selO = kN[on], kO[oo], selN[on], selO[oo]
    both, iN, iO = np.intersect1d(kN, kO, assume_unique=True, return_indices=True)
    rN, rO = selN[iN], selO[iO]  # row indices of matched events
    log(f"\n[subset] new events not in old: {len(kN) - len(both)}  (expect 0)")
    if len(kN) - len(both):
        miss = np.setdiff1d(kN, kO)[:20]
        log("[subset]   first ones (run,lumi,event):",
            [((m >> 50) + 355000, (m >> 36) & (2**14 - 1), m & (2**36 - 1)) for m in miss])

    old_pass = dmo[rO] < PHI_WIN
    same_rN, same_rO = rN[old_pass], rO[old_pass]
    rec_rN, rec_rO = rN[~old_pass], rO[~old_pass]
    oldpass_all = selO[dmo[selO] < PHI_WIN]
    lost = np.setdiff1d(ko[oldpass_all], kn[selN])
    log(f"[split] old best passes phi20 & in new (same): {len(same_rN)}")
    log(f"[split] old best passes phi20 but NOT in new (lost): {len(lost)}")
    if len(lost):
        lm = np.isin(ko, lost)
        for k, d in list(zip(ko[lm], dmo[lm]))[:30]:
            log(f"[split]   lost run {(k >> 50) + 355000} lumi {(k >> 36) & (2**14-1)} evt {k & (2**36-1)} |dM| {d*1e3:.6f} MeV")
    log(f"[split] old best fails phi20, new has a phi20 cand (recovered): {len(rec_rN)}")
    n_oldpass = len(oldpass_all)
    log(f"[split] anchor: new rows in common lumis {len(kN)} =? same {len(same_rN)} + recovered {len(rec_rN)}"
        f" = {len(same_rN) + len(rec_rN)};  recovered / old-phi20 = {len(rec_rN) / max(n_oldpass, 1):.4f}")

    # recovered: new best pt must be <= old best pt (old candidates are a superset)
    if len(rec_rN):
        dpt = N["Ds_pt"][rec_rN].astype(np.float64) - O["Ds_pt"][rec_rO].astype(np.float64)
        log(f"[recovered] new Ds_pt - old Ds_pt: max {dpt.max():.4g} (expect <0), n>=0: {(dpt >= 0).sum()}")
        dnc = N["nCandidates"][rec_rN].astype(int) - O["nCandidates"][rec_rO].astype(int)
        log(f"[recovered] nCandidates new-old: max {dnc.max()}")
    if len(same_rN):
        dnc = N["nCandidates"][same_rN].astype(int) - O["nCandidates"][same_rO].astype(int)
        log(f"[same] nCandidates new-old: max {dnc.max()}, frac <0 {(dnc < 0).mean():.4f}")

    # --- branch-by-branch on the same events
    log(f"\n[same] per-branch comparison on {len(same_rN)} same events "
        f"(bitwise equal / NaN-aware; listed only if not all equal)")
    n_all_eq = 0
    for b in br + ["sum_" + b for b in JAGGED_SUMS]:
        a, c = N[b][same_rN], O[b][same_rO]
        if a.dtype.kind == "f":
            eq = (a == c) | (np.isnan(a) & np.isnan(c))
            af, cf = a.astype(np.float64), c.astype(np.float64)
            fin = np.isfinite(af) & np.isfinite(cf)
            d = np.abs(af - cf)[fin]
            rel = (d / np.maximum(np.abs(cf[fin]), 1e-30))
            maxd = d.max() if d.size else 0.0
            maxr = rel.max() if rel.size else 0.0
            nanmis = (np.isnan(a) != np.isnan(c)).sum()
        else:
            eq = a == c
            maxd = np.abs(a.astype(np.float64) - c.astype(np.float64)).max() if a.size else 0
            maxr, nanmis = float("nan"), 0
        if eq.all():
            n_all_eq += 1
            continue
        tag = "  (expected)" if b in EXPECTED_DIFF else ""
        log(f"[same]   {b:40s} equal {eq.mean():.6f}  n_diff {(~eq).sum():7d}  max|d| {maxd:.3g}"
            f"  max rel {maxr:.3g}  NaN mismatch {nanmis}{tag}")
    log(f"[same] branches fully bitwise equal: {n_all_eq} / {len(br) + len(JAGGED_SUMS)}")

    # --- per-row: how many same events are bitwise identical in all non-expected branches
    rows_eq = np.ones(len(same_rN), bool)
    for b in br:
        if b in EXPECTED_DIFF:
            continue
        a, c = N[b][same_rN], O[b][same_rO]
        e = a == c
        if a.dtype.kind == "f":
            e |= np.isnan(a) & np.isnan(c)
        rows_eq &= e
    log(f"[same] events identical in every branch except {EXPECTED_DIFF}: {rows_eq.sum()} / {len(rows_eq)}")

    cat = np.full(len(kn), -1, np.int8)
    cat[same_rN] = 0
    cat[rec_rN] = 1
    np.savez_compressed(
        f"{outdir}/Run2022{era}.npz",
        new_key=kn, new_cat=cat, new_Ds_mass=N["Ds_mass"], new_phi_mass=N["phi_mass"], new_Ds_pt=N["Ds_pt"],
        old_key=ko, old_Ds_mass=O["Ds_mass"], old_phi_mass=O["phi_mass"], old_Ds_pt=O["Ds_pt"],
    )
    with open(f"{outdir}/Run2022{era}.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()

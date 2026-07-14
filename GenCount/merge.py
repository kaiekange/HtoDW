"""Sum the per-shard counts written by gencount.py and report the gen-level totals.

Usage: python3 merge.py
"""
import glob
import os

BASE = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/GenCount_phipimu"

NJOBS = 50  # must match condor_submit.py (1000 files / GroupSize 20)

# nEventsGen_ from the tuplizer Info tree: the same gen filter (munu only), but
# applied *after* jvmAppliedEventFilter. The pre-veto munu counts must exceed these.
POSTVETO_MUGEN = {"DsW": 225720, "DsstW": 232565}

for channel in ("DsW", "DsstW"):
    outdir = os.path.join(BASE, channel)
    files = glob.glob(os.path.join(outdir, "count_*.txt"))

    total = mu = mu_gen = e = e_gen = 0
    for path in files:
        with open(path) as f:
            t, m, mg, el, eg = (int(x) for x in f.read().split())
        total += t
        mu += m
        mu_gen += mg
        e += el
        e_gen += eg

    print(f"=== {channel} ===")
    if len(files) != NJOBS:
        print(f"  !! INCOMPLETE: {len(files)}/{NJOBS} shards -- totals below are NOT final")
    if total == 0:
        print("  no events counted")
        continue

    print(f"  MiniAOD events (pre veto)     {total:>10d}")
    print(f"  W->munu                       {mu:>10d}   ({mu/total:.4f} of total)")
    print(f"    & Ds->phipi->KKpi           {mu_gen:>10d}   ({mu_gen/mu:.4f} of W->munu)" if mu else "")
    print(f"  W->enu                        {e:>10d}   ({e/total:.4f} of total)")
    print(f"    & Ds->phipi->KKpi           {e_gen:>10d}   ({e_gen/e:.4f} of W->enu)" if e else "")
    print(f"  munu+enu                      {mu+e:>10d}   ({(mu+e)/total:.4f} of total)")

    postveto = POSTVETO_MUGEN[channel]
    if mu_gen:
        print(f"  tuplizer nEventsGen_ (post veto, munu) {postveto:>7d}   "
              f"(veto kept {postveto/mu_gen:.4f})")
    if mu_gen < postveto:
        print("  !! pre-veto munu count is BELOW the post-veto count -- replication is wrong")

"""Submit gencount.py, one job per GROUPSIZE MiniAOD files, one cluster per sample.

    python3 condor_submit_gencount.py <DsW|DsstW>_<2022|2022EE> [...] [--njobs N]

Run from this directory on lxplus8 (el8, CMSSW_13_0_13). The MiniAOD file lists are the
ones the nano production used (CMSSW_13_0_13 .../condor/filelists/<channel>_<era>.txt),
so the counted files are exactly the files behind Runs.genEventSumw.
"""
import argparse
import math
import os
import subprocess

GROUPSIZE = 10
FILELISTS = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/condor/filelists"

ap = argparse.ArgumentParser()
ap.add_argument("targets", nargs="+")
ap.add_argument("--njobs", type=int, default=0, help="override job count (test job)")
args = ap.parse_args()

os.makedirs("condor_logs", exist_ok=True)
for target in args.targets:
    channel, era = target.split("_")
    assert channel in ("DsW", "DsstW") and era in ("2022", "2022EE"), target
    filelist = f"{FILELISTS}/{channel}_{era}.txt"
    nfiles = sum(1 for l in open(filelist) if l.strip())
    njobs = args.njobs or math.ceil(nfiles / GROUPSIZE)
    sched = f"condor_logs/HTCondor_gencount_{target}.log"
    if os.path.exists(sched):
        os.remove(sched)
    desc = f"""
executable = run_gencount.sh
arguments = {channel} {era} {filelist} $(Process) {GROUPSIZE}
log = {sched}
JobBatchName = gencount_{target}
request_cpus = 1
request_memory = 2G
request_disk = 100M
+JobFlavour = "longlunch"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue {njobs}
"""
    print(f"Submitting {njobs} gencount jobs for {target} ({nfiles} files)")
    subprocess.run(["condor_submit"], input=desc.encode(), check=True)

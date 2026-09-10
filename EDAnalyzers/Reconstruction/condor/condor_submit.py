"""Submit custom-NanoAOD production jobs.

    python3 condor_submit.py DsW 2022                 # signal, GROUPSIZE PAT files per job
    python3 condor_submit.py DsstW 2022EE
    python3 condor_submit.py WJ1 2022 --gentag        # W+jets 1-jet bin, GenTag tuplizer
    python3 condor_submit.py WJ1 2022 --ss            # same-sign-pi control region
    python3 condor_submit.py DsW 2022 --njobs 1       # smoke test: only the first job

The input file list is resolved HERE (on lxplus, where xrootd tools work) and written to
filelists/<sample>_<year>.txt as root:// URLs, numerically ordered. Each job slices it by
JOBID/GROUPSIZE. Re-run this to pick up files that finished since (private signal
productions may still be growing) -- same model as the static W+jets dumps.
"""
import argparse
import math
import os
import re
import subprocess

GROUPSIZE = 10  # PAT files per job; 2-thread nano ~1 min/file + ~2 min startup

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = {
    "2022":   dict(ver="20260717", patyear="2022"),
    "2022EE": dict(ver="20260825", patyear="2022EE"),
}
SIGNAL_DAS = {"DsW": "ggHtoDsW_phipi", "DsstW": "ggHtoDsstW_phipi"}
WJ_BIN = {"WJ0": "0J", "WJ1": "1J", "WJ2": "2J"}
WJ_FILELIST_DIR = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WtoLNu2Jets"
EOS_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA"


def build_filelist(sample, year):
    """Write filelists/<sample>_<year>.txt (root:// URLs, numeric order) and return its path."""
    era = ERA[year]
    os.makedirs(f"{HERE}/filelists", exist_ok=True)
    out = f"{HERE}/filelists/{sample}_{year}.txt"

    if sample in SIGNAL_DAS:
        patdir = f"/store/user/kakang/Analysis/{SIGNAL_DAS[sample]}/{era['ver']}/Simulation/PAT/{era['patyear']}"
        ls_out = subprocess.check_output(["xrdfs", "maite.iihe.ac.be", "ls", patdir]).decode()
        idx = sorted(int(m) for m in re.findall(r"output_(\d+)\.root", ls_out))
        urls = [f"root://maite.iihe.ac.be/{patdir}/output_{i}.root" for i in idx]
    elif sample in WJ_BIN:
        src = f"{WJ_FILELIST_DIR}/filelist_{WJ_BIN[sample]}_{year}.txt"
        with open(src) as fh:
            urls = ["root://xrootd-cms.infn.it/" + ln.strip() for ln in fh if ln.strip()]
    else:
        raise SystemExit(f"unknown sample {sample!r}")

    with open(out, "w") as fh:
        fh.write("\n".join(urls) + "\n")
    print(f"{out}: {len(urls)} files")
    return out, len(urls)


def out_subdir(sample, ss, gentag):
    if sample in SIGNAL_DAS:
        base = f"Signal_{sample}_H_phipimu"
    else:
        base = f"WtoLNu2Jets_{WJ_BIN[sample]}_H_phipimu"
    if gentag:
        base += "_GenTag"
    if ss:
        base += "_SS"
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample", choices=list(SIGNAL_DAS) + list(WJ_BIN))
    ap.add_argument("year", choices=list(ERA))
    ap.add_argument("--ss", action="store_true", help="same-sign-pi control region")
    ap.add_argument("--gentag", action="store_true", help="HDsPhiPiMuGenTagProducer")
    ap.add_argument("--njobs", type=int, default=0, help="override job count (smoke test)")
    args = ap.parse_args()

    filelist, nfiles = build_filelist(args.sample, args.year)
    njobs = args.njobs or math.ceil(nfiles / GROUPSIZE)

    outdir = f"{EOS_BASE}/MC/{args.year}/{out_subdir(args.sample, args.ss, args.gentag)}"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir + "/log", exist_ok=True)
    os.makedirs(f"{HERE}/condor_logs", exist_ok=True)

    ss = "1" if args.ss else "0"
    tup = "gentag" if args.gentag else "nominal"
    batch = out_subdir(args.sample, args.ss, args.gentag) + f"_{args.year}"
    sched_log = f"{HERE}/condor_logs/HTCondor_{batch}.log"
    if os.path.exists(sched_log):
        os.remove(sched_log)

    submit = f"""
executable = {HERE}/run.sh
arguments = {outdir} {GROUPSIZE} $(Process) {args.sample} {args.year} {ss} {tup} {filelist}
log = {sched_log}
JobBatchName = {batch}
request_cpus = 2
request_memory = 4G
request_disk = 2G
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue {njobs}
"""
    print(submit)
    subprocess.run(["condor_submit"], input=submit.encode(), check=True)


if __name__ == "__main__":
    main()

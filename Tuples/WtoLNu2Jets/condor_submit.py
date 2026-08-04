"""Submit the WtoLNu2Jets tuplizer jobs for one jet bin, OS or SS.

Replaces the eight copies of this script there used to be (inclusive/0J/1J/2J x
OS/SS), which differed only in the output directory, the file list, the batch name
and the scheduler log.

    python3 condor_submit.py <incl|0J|1J|2J> [SS]

Names are reproduced exactly as the per-bin scripts built them, so the already
produced output directories stay the ones being written to.
"""
import sys
import subprocess
import os
import math

# Per jet bin: file-list suffix and files-per-job. The inclusive sample has far
# fewer, larger files than the jet-binned ones, hence the smaller group size.
JETBINS = {
    "incl": ("",    5),
    "0J":   ("_0J", 20),
    "1J":   ("_1J", 20),
    "2J":   ("_2J", 20),
}


def submit_jobs(jetbin, same_sign):
    suffix, GroupSize = JETBINS[jetbin]
    ss_tag = "_SS" if same_sign else ""

    outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/2022/WtoLNu2Jets{suffix}_H_phipimu{ss_tag}"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir + "/log", exist_ok=True)

    scheduler_log = f"./HTCondor{suffix}{ss_tag}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    NFiles = sum(1 for _ in open(f"filelist{suffix}.txt"))
    njobs = math.ceil(NFiles / GroupSize)

    submit_description = f"""
executable = run.sh
arguments = {outdir} {GroupSize} $(Process) {jetbin} {1 if same_sign else 0}
log = {scheduler_log}
JobBatchName = WtoLNu2Jets{suffix}_H_phipimu{ss_tag}
request_cpus = 1
request_memory = 8G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue {njobs}
"""
    print(f"{jetbin}{ss_tag}: {NFiles} files / {GroupSize} per job -> {njobs} jobs -> {outdir}")
    subprocess.run(["condor_submit"], input=submit_description.encode())


if __name__ == "__main__":
    args = sys.argv[1:]
    # Validate before makedirs: an unrecognised token would otherwise silently create
    # a directory nothing reads and submit a few hundred jobs writing into it.
    if not args or args[0] not in JETBINS:
        sys.exit(f"usage: python3 condor_submit.py <{'|'.join(JETBINS)}> [SS]")
    sys.exit(submit_jobs(args[0], len(args) > 1 and args[1].upper() == "SS"))

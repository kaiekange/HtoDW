"""Submit one phipimu DNN training+scoring job to HTCondor.

    python3 condor_submit.py ssdata
    python3 condor_submit.py wjets

One job per background configuration. Each gets its own scheduler log so two
concurrent submissions do not interleave into the same file.
"""

import os
import subprocess
import sys

import common

if len(sys.argv) != 2:
    sys.exit(f"usage: python3 {os.path.basename(__file__)} <{'|'.join(common.TRAININGS)}>")
BKG = sys.argv[1]
cfg = common.resolve(BKG)          # validate the token before creating anything

OUTDIR = cfg["out_dir"]
os.makedirs(OUTDIR, exist_ok=True)


def submit_jobs():
    submit_description = f"""
executable = run.sh
arguments = {BKG}
log = ./HTCondor_{BKG}.log
JobBatchName = phipimu_DNN_{BKG}
request_cpus = 4
request_memory = 8G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1
"""
    print(f"Submitting {BKG} ({cfg['label']}) -> {OUTDIR}")
    subprocess.run(["condor_submit"], input=submit_description.encode())


if __name__ == "__main__":
    sys.exit(submit_jobs())

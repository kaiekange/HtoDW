import sys
import subprocess
import os

OUTDIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets"

os.makedirs(OUTDIR, exist_ok=True)

def submit_jobs():

    submit_description = f"""
executable = run.sh
log = ./HTCondor.log
JobBatchName = phipimu_DNN
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
    subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    sys.exit(submit_jobs())

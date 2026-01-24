import sys
import subprocess
import os;

def submit_jobs():

    year = 2022
    outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/{year}/match_Info"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    scheduler_log = f"./HTCondor.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {outdir} $(Process)
log = {scheduler_log}
JobBatchName = MC_H2DW_{year}_match_Info
request_cpus = 1
request_memory = 8G
request_disk = 10M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1000
"""
    subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    sys.exit(submit_jobs())
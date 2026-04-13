import sys
import subprocess
import os;

groupsize=20

def submit_jobs():

    outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    scheduler_log = f"./HTCondor.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {outdir} {groupsize} $(Process)
log = {scheduler_log}
request_cpus = 1
request_memory = 8G
request_disk = 10M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 50
"""
    subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    sys.exit(submit_jobs())

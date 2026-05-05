import sys
import subprocess
import os

GroupSize = 1

def submit_jobs():

    outdir="/eos/home-k/kakang/H2DW/Analysis/MC/2022/RecoAll_DsW_KKpimu"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    scheduler_log = f"./HTCondor.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    njobs = 1000/GroupSize

    submit_description = f"""
executable = run.sh
arguments = {outdir} {GroupSize} $(Process)
log = {scheduler_log}
JobBatchName = RecoAll_DsW_KKpimu
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
    subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    sys.exit(submit_jobs())

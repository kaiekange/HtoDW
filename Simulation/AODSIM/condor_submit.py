import sys
import subprocess
import os;

def submit_jobs():

    ver = 20260119
    year = 2022
    indir = f"/eos/home-k/kakang/H2DW/Datasets/{ver}/Simulation/LHEGS/{year}/260119_142940/0000"
    outdir = f"/eos/home-k/kakang/H2DW/Datasets/{ver}/Simulation/AODSIM/{year}"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    scheduler_log = f"./HTCondor_1.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {indir} {outdir} $(Process)
log = {scheduler_log}
JobBatchName = MC_ggHtoDsW_AODSIM_{year}_{ver}
request_cpus = 1
request_memory = 8G
request_disk = 10M
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
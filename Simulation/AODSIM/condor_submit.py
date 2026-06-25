import sys
import subprocess
import os

def submit_jobs():

    ver = 20260624
    year = 2022
    # mode = "ggHtoDsW_KKpi"
    # mode = "ggHtoDsstW_KKpi"
    # mode = "ggHtoDsW_phipi"
    mode = "ggHtoDsstW_phipi"

    indir = f"/store/user/kakang/Analysis/{mode}/{ver}/Simulation/LHEGS/{year}"
    outdir = f"/store/user/kakang/Analysis/{mode}/{ver}/Simulation/AODSIM/{year}"

    os.makedirs("./tmp", exist_ok=True)

    subprocess.run(["xrdfs", "root://maite.iihe.ac.be", "mkdir", "-p", f"{outdir}/log"])
    subprocess.run(["./mkconfig.sh", mode])

    scheduler_log = f"condor_logs/HTCondor_{mode}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {mode} {indir} {outdir} $(Process)
log = {scheduler_log}
JobBatchName = MC_{mode}_AODSIM_{year}_{ver}
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

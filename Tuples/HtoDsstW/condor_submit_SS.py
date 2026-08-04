import re
import sys
import subprocess
import os

GroupSize = 10

# 20260717 production caps out short of 1000 (2 AODSIM jobs still running as of
# 2026-08-04, 6 AODSIM outputs are broken and will never produce PAT) -- checked live
# so njobs tracks the sample actually on storage instead of assuming all 1000 exist.
PATDIR = "/store/user/kakang/Analysis/ggHtoDsstW_phipi/20260717/Simulation/PAT/2022"

def submit_jobs():

    outdir="/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_DsstW_H_phipimu_SS"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    scheduler_log = "./HTCondor_SS.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    ls_out = subprocess.check_output(["xrdfs", "maite.iihe.ac.be", "ls", PATDIR]).decode()
    max_index = max(int(m) for m in re.findall(r"output_(\d+)\.root", ls_out))
    njobs = (max_index + GroupSize - 1) // GroupSize

    submit_description = f"""
executable = run_SS.sh
arguments = {outdir} {GroupSize} $(Process)
log = {scheduler_log}
JobBatchName = Signal_DsstW_H_phipimu_SS
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

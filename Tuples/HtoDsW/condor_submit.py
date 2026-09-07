import re
import sys
import subprocess
import os

GroupSize = 10

# Single YEAR knob drives PATDIR + ver + outdir together (same principle as applyJmeEra /
# the grid-production mkconfig.sh scripts): querying njobs live from PATDIR still tracks
# the actual sample size instead of an assumed file count.
ERA = {
    "2022":   dict(ver="20260717", patyear="2022"),
    "2022EE": dict(ver="20260825", patyear="2022EE"),
}

def submit_jobs(year):
    era = ERA[year]
    patdir = f"/store/user/kakang/Analysis/ggHtoDsW_phipi/{era['ver']}/Simulation/PAT/{era['patyear']}"
    outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/{year}/Signal_H_phipimu"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    # Scheduler logs go in condor_logs/ rather than next to the script: one per
    # (era, bin, sign) piles up fast in the working directory.
    os.makedirs("condor_logs", exist_ok=True)
    scheduler_log = f"condor_logs/HTCondor_{year}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    ls_out = subprocess.check_output(["xrdfs", "maite.iihe.ac.be", "ls", patdir]).decode()
    max_index = max(int(m) for m in re.findall(r"output_(\d+)\.root", ls_out))
    njobs = (max_index + GroupSize - 1) // GroupSize

    submit_description = f"""
executable = run.sh
arguments = {outdir} {GroupSize} $(Process) {year}
log = {scheduler_log}
JobBatchName = Signal_H_phipimu_{year}
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
    if len(sys.argv) != 2 or sys.argv[1] not in ERA:
        sys.exit(f"usage: python3 condor_submit.py <{'|'.join(ERA)}>")
    sys.exit(submit_jobs(sys.argv[1]))

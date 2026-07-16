import sys
import subprocess
import os

model_name = sys.argv[1] if len(sys.argv) > 1 else "v3"

os.makedirs(f"/eos/home-k/kakang/H2DW/Analysis/MVA/BDT/{model_name}", exist_ok=True)

def submit_jobs():
    submit_description = f"""
executable = run.sh
arguments = {model_name}
log = ./HTCondor.log
JobBatchName = phipimu_BDT
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

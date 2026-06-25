import sys
import subprocess
import os

model_name = "v1"

os.makedirs(f"/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/{model_name}", exist_ok=True)

def submit_jobs():
    submit_description = f"""
executable = run.sh
arguments = {model_name}
log = ./HTCondor.log
JobBatchName = Gamma_DNN
request_cpus = 1
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

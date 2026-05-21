import sys
import subprocess
import os;

model_name = "Model_200ep_0p3do_4ly"

os.makedirs(f"/eos/user/k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/DNN/{model_name}", exist_ok=True)

def submit_jobs():

    submit_description = f"""
executable = run.sh
arguments = {model_name}
log = ./HTCondor.log
JobBatchName = Ds_DNN
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

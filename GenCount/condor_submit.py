import sys
import subprocess
import os

GroupSize = 20

# channel -> highest MiniAOD file index (DsstW is missing output_188, see gencount.py)
CHANNELS = {
    "DsW":   1000,
    "DsstW": 1000,
}

def submit_jobs():

    for channel, nfiles in CHANNELS.items():

        outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/2022/GenCount_phipimu/{channel}"

        os.makedirs(outdir, exist_ok=True)
        os.makedirs(outdir+"/log", exist_ok=True)

        scheduler_log = f"./HTCondor_{channel}.log"
        if os.path.exists(scheduler_log):
            os.remove(scheduler_log)

        njobs = -(-nfiles // GroupSize)  # ceil, so the last short shard is kept

        submit_description = f"""
executable = run.sh
arguments = {outdir} {GroupSize} $(Process) {channel}
log = {scheduler_log}
JobBatchName = GenCount_{channel}_phipimu
request_cpus = 1
request_memory = 2G
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

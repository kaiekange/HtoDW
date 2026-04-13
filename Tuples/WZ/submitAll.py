import sys
import math
import subprocess
import os
from pathlib import Path

groupsize = 20

def save_das_files(dataset, outfile):
    cmd = f'dasgoclient --query="file dataset={dataset}"'

    out = subprocess.check_output(cmd, shell=True, text=True)

    lines = [l.strip() for l in out.splitlines() if l.strip()]
    Path(outfile).write_text("\n".join(lines))

    return len(lines)

def submitAll():
    dataset = "/WZtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM"
    outdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/WZ"
    filelist = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WZ/filelist.txt"

    njobs = save_das_files(dataset, filelist)
    ngroups = math.ceil(njobs / groupsize)

    os.makedirs(f"{outdir}/log", exist_ok=True)

    scheduler_log = "HTCondor.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = runAll.sh
arguments = {filelist} {outdir} {groupsize} $(Process)
log = {scheduler_log}
request_cpus = 1
request_memory = 8G
request_disk = 10M
+JobFlavour = "tomorrow"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue {ngroups}
"""
    subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    submitAll()
    sys.exit(0)

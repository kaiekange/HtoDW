"""Submit HTCondor jobs running apply_corrections.py, one job per nano output file.

    python3 condor_submit_apply_corrections.py <era>/<sample dir> [<era>/<sample dir> ...]
    e.g. python3 condor_submit_apply_corrections.py 2022/Signal_DsW_H_phipimu 2022EE/Signal_DsW_H_phipimu

One cluster per sample directory, one job per output_*.root in it. Friend trees go to
.../NanoMVA/MC/<era>/Corrections/<sample dir>/ (overwritten). Scheduler log per sample
under condor_logs/, per-job output under condor_logs/apply_corrections/.
"""
import glob
import os
import subprocess
import sys

NANO = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"

if len(sys.argv) < 2:
    sys.exit(__doc__)

os.makedirs("condor_logs/apply_corrections", exist_ok=True)

for target in sys.argv[1:]:
    era, sample = target.strip("/").split("/")
    files = sorted(glob.glob(f"{NANO}/{era}/{sample}/output_*.root"))
    if not files:
        sys.exit(f"no output_*.root under {NANO}/{era}/{sample}")
    tag = f"{era}_{sample}"
    file_list = f"condor_logs/apply_corrections/{tag}.txt"
    with open(file_list, "w") as f:
        f.write("\n".join(files) + "\n")
    scheduler_log = f"condor_logs/HTCondor_apply_corrections_{tag}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run_apply_corrections.sh
arguments = $(INFILE)
log = {scheduler_log}
JobBatchName = apply_corrections_{tag}
request_cpus = 1
request_memory = 4G
request_disk = 200M
+JobFlavour = "longlunch"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue INFILE from {file_list}
"""
    print(f"Submitting {len(files)} apply_corrections jobs for {tag}")
    subprocess.run(["condor_submit"], input=submit_description.encode())

"""Submit HTCondor jobs running build_btag_efficiency_map.py, one job per sample group.

    python3 condor_submit_btag_eff.py [group ...]

With no arguments, submits all groups. Run from this directory (the executable and
condor_logs/ are relative paths); each job gets its own scheduler log under condor_logs/.
"""
import os
import subprocess
import sys

ALL_GROUPS = [f"{s}_{e}" for e in ("2022", "2022EE") for s in ("DsW", "DsstW", "WtoLNu2Jets")]
GROUPS = sys.argv[1:] if len(sys.argv) > 1 else ALL_GROUPS
for g in GROUPS:
    if g not in ALL_GROUPS:
        sys.exit(f"unknown group {g!r}; choices are {ALL_GROUPS}")

os.makedirs("condor_logs", exist_ok=True)

for group in GROUPS:
    scheduler_log = f"condor_logs/HTCondor_btag_eff_{group}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run_btag_eff.sh
arguments = {group}
log = {scheduler_log}
JobBatchName = btag_eff_{group}
request_cpus = 1
request_memory = 4G
request_disk = 200M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1
"""
    print(f"Submitting b-tag efficiency map job for {group}")
    subprocess.run(["condor_submit"], input=submit_description.encode())

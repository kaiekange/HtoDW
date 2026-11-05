"""Submit HTCondor jobs running build_btag_efficiency_map.py, one job per sample group.

    python3 condor_submit_btag_eff.py [--jme] [group ...]

--jme: maps on the JME-corrected jet pt (btag_efficiency_jme/); with no group given, only the
groups apply_corrections uses (not the pooled WtoLNu2Jets_<era> check map).

With no arguments, submits all groups. Run from this directory (the executable and
condor_logs/ are relative paths); each job gets its own scheduler log under condor_logs/.
"""
import os
import subprocess
import sys

ALL_GROUPS = [f"{s}_{e}" for e in ("2022", "2022EE")
              for s in ("DsW", "DsstW", "WtoLNu2Jets", "WtoLNu2Jets_0J", "WtoLNu2Jets_1J", "WtoLNu2Jets_2J",
                        "WtoLNu2Jets_0J_SS", "WtoLNu2Jets_1J_SS", "WtoLNu2Jets_2J_SS",
                        "TTto2L2Nu", "TTtoLNu2Q", "TTto2L2Nu_SS", "TTtoLNu2Q_SS", "DYto2L2Jets_0J", "DYto2L2Jets_1J", "DYto2L2Jets_2J",
                        "DYto2L2Jets_0J_SS", "DYto2L2Jets_1J_SS", "DYto2L2Jets_2J_SS")]
JME = "--jme" in sys.argv[1:]
ARGV = [x for x in sys.argv[1:] if x != "--jme"]
DEFAULT = [g for g in ALL_GROUPS if not g.startswith("WtoLNu2Jets_2022")] if JME else ALL_GROUPS
GROUPS = ARGV if ARGV else DEFAULT
for g in GROUPS:
    if g not in ALL_GROUPS:
        sys.exit(f"unknown group {g!r}; choices are {ALL_GROUPS}")

os.makedirs("condor_logs", exist_ok=True)

for group in GROUPS:
    tag = "_jme" if JME else ""
    scheduler_log = f"condor_logs/HTCondor_btag_eff_{group}{tag}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run_btag_eff.sh
arguments = {group}{" --jme" if JME else ""}
log = {scheduler_log}
JobBatchName = btag_eff{tag}_{group}
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

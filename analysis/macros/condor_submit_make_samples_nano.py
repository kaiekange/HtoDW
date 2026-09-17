"""Submit HTCondor jobs running make_samples_nano.py worker, FILES_PER_JOB nano files per job.

    python3 condor_submit_make_samples_nano.py <group> [--pieces=<era>/<sample>,...] [--correct] [--dry-run]
    python3 condor_submit_make_samples_nano.py <group> [--pieces=...] [--correct] --only-missing
    <group>: Signal_merged | WtoLNu2Jets_merged | Data_OS_merged | Data_SS_merged

Default: for every piece of the group, rewrite the chunk lists (numeric order of
output_N, FILES_PER_JOB per chunk), delete that piece's old skim chunks, submit every chunk.
--only-missing: keep the existing lists, submit only chunks whose .json is absent
(the worker writes the .json last, so absent = not finished successfully).
--correct: each job first runs apply_corrections.py on its files (writing the friends),
then skims -- no separate corrections submission; friends need not exist beforehand.
--pieces: restrict to these pieces of the group (e.g. a data era whose CRAB task is done);
merge still requires every piece of the group.

Then, once all jobs are done:  make_samples_nano.py merge <group>  (inside cmssw-el9).
One cluster per group. Chunk lists under condor_logs/make_samples_nano/lists/, scheduler
log under condor_logs/, per-job output under condor_logs/make_samples_nano/.
"""
import glob
import os
import subprocess
import sys

import make_samples_nano as msn

args = [a for a in sys.argv[1:] if not a.startswith("--")]
flags = {a.split("=", 1)[0]: a.split("=", 1)[1] if "=" in a else None
         for a in sys.argv[1:] if a.startswith("--")}
if len(args) != 1 or args[0] not in msn.GROUPS or set(flags) - {"--dry-run", "--only-missing",
                                                                 "--pieces", "--correct"}:
    sys.exit(__doc__)
group = args[0]
dry_run = "--dry-run" in flags
only_missing = "--only-missing" in flags
correct = "--correct" in flags
pieces = msn.GROUPS[group]
if "--pieces" in flags:
    wanted = [tuple(x.split("/")) for x in (flags["--pieces"] or "").split(",") if x]
    unknown = [w for w in wanted if w not in pieces]
    if not wanted or unknown:
        sys.exit(f"--pieces must name pieces of {group}; unknown: {unknown}")
    pieces = [p for p in pieces if p in wanted]

os.makedirs(msn.LIST_DIR, exist_ok=True)
jobs = []
for era, sample in pieces:
    if only_missing:
        lists = msn.piece_chunks(era, sample)
        if not lists:
            sys.exit(f"--only-missing but no chunk lists for {era}/{sample}")
        for lp in lists:
            k = int(lp.rsplit("_chunk", 1)[1][:-len(".txt")])
            out = msn.skim_path(era, sample, k)
            if not os.path.exists(out[:-len(".root")] + ".json"):
                jobs.append((lp, out))
        continue

    files = msn.piece_nano_files(era, sample)
    if not files:
        sys.exit(f"no nano files for {era}/{sample}")
    is_data = msn.parse_nano_path(files[0])[3]
    norm = f"{msn.MC_BASE}/{era}/Corrections/{sample}/normalization.json"
    has_norm = is_data or os.path.exists(norm)
    missing = [] if correct else [f for f in files
                                  if not os.path.exists(msn.parse_nano_path(f)[2])]
    if not has_norm or missing:
        sys.exit(f"{era}/{sample}: normalization.json present={has_norm}, "
                 f"{len(missing)} friend files missing (first: {missing[:1]})")
    chunks = [files[i:i + msn.FILES_PER_JOB] for i in range(0, len(files), msn.FILES_PER_JOB)]
    print(f"{era}/{sample}: {len(files)} files -> {len(chunks)} chunks")
    if dry_run:
        continue
    for old in glob.glob(f"{msn.LIST_DIR}/{era}_{sample}_chunk*.txt"):
        os.remove(old)
    for old in glob.glob(f"{msn.SKIM_DIR}/{era}/{sample}/chunk_*"):
        os.remove(old)
    for k, chunk in enumerate(chunks):
        lp = f"{msn.LIST_DIR}/{msn.chunk_name(era, sample, k)}.txt"
        with open(lp, "w") as f:
            f.write("\n".join(chunk) + "\n")
        jobs.append((lp, msn.skim_path(era, sample, k)))

if dry_run:
    sys.exit(0)
if not jobs:
    sys.exit("nothing to submit")

os.makedirs("condor_logs/make_samples_nano", exist_ok=True)
tag = group if "--pieces" not in flags else f"{group}_" + "_".join(s for _, s in pieces)
job_table = f"condor_logs/make_samples_nano/{tag}_jobs.txt"
with open(job_table, "w") as f:
    f.write("\n".join(f"{lp} {out}" for lp, out in jobs) + "\n")
scheduler_log = f"condor_logs/HTCondor_make_samples_nano_{tag}.log"
if os.path.exists(scheduler_log):
    os.remove(scheduler_log)

submit_description = f"""
executable = run_make_samples_nano.sh
arguments = $(LIST) $(OUT){" --correct" if correct else ""}
log = {scheduler_log}
JobBatchName = make_samples_nano_{tag}
request_cpus = 1
request_memory = 4G
request_disk = 2G
+JobFlavour = "longlunch"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue LIST,OUT from {job_table}
"""
print(f"Submitting {len(jobs)} make_samples_nano jobs for {tag}")
subprocess.run(["condor_submit"], input=submit_description.encode())

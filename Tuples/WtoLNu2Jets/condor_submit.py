"""Submit the WtoLNu2Jets tuplizer jobs for one jet bin, OS or SS.

Replaces the copies of this script there used to be, which differed only in the
output directory, the file list, the batch name and the scheduler log.

    python3 condor_submit.py <2022|2022EE> <0J|1J|2J> [SS] [--gentag] [--jobs 27,190,429]

--jobs resubmits only those JOBIDs (targeted refill after transient xrootd
FileReadErrors), leaving every other output file alone. run.sh rm -f's its own
output first, so a refilled job cleanly replaces the truncated file it left.

--gentag runs H_phipimu_GenTag instead of H_phipimu: same reco content and the same
cuts, plus the per-candidate truth tagging of what the three KKpi tracks really are.
It writes to a parallel *_GenTag output directory, so the two productions never mix.

Names are reproduced exactly as the per-bin scripts built them, so the already
produced output directories stay the ones being written to.
"""
import sys
import subprocess
import os
import math

# Per jet bin: file-list suffix and files-per-job.
JETBINS = {
    "0J":   ("_0J", 20),
    "1J":   ("_1J", 20),
    "2J":   ("_2J", 20),
}

# 2022 = Run3Summer22 (preEE), 2022EE = Run3Summer22EE (postEE). The era is the top-level
# output directory, never a suffix on the sample name: preEE output keeps the paths it has
# always had, postEE goes to a parallel MC/2022EE tree. Mixing the two is the one thing
# that must not happen silently, so it is visible in the path rather than in a comment.
ERAS = ("2022", "2022EE")


def submit_jobs(era, jetbin, same_sign, jobids=None, gentag=False):
    suffix, GroupSize = JETBINS[jetbin]
    ss_tag = "_SS" if same_sign else ""
    # Separate output tree, never a suffix inside the existing one: the GenTag tuple carries
    # 28 branches the plain one does not, and hadd'ing the two together would fail or, worse,
    # silently drop them.
    gt_tag = "_GenTag" if gentag else ""
    tuplizer = "H_phipimu_GenTag" if gentag else "H_phipimu"

    outdir = f"/eos/home-k/kakang/H2DW/Analysis/MC/{era}/WtoLNu2Jets{suffix}_H_phipimu{ss_tag}{gt_tag}"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir + "/log", exist_ok=True)

    # Scheduler logs go in condor_logs/ rather than next to the script: one per
    # (era, bin, sign) piles up fast in the working directory.
    os.makedirs("condor_logs", exist_ok=True)
    scheduler_log = f"condor_logs/HTCondor_{era}{suffix}{ss_tag}{gt_tag}{'_refill' if jobids else ''}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    NFiles = sum(1 for _ in open(f"filelist_{jetbin}_{era}.txt"))
    njobs = math.ceil(NFiles / GroupSize)

    if jobids:
        out_of_range = [j for j in jobids if not 0 <= j < njobs]
        if out_of_range:
            sys.exit(f"JOBID(s) {out_of_range} outside 0..{njobs-1} for {era} {jetbin}")
        # $(jid) rather than $(Process): Process would renumber the refilled jobs 0,1,2
        # and overwrite the wrong output files.
        queue_stmt = "queue jid in (" + ", ".join(str(j) for j in jobids) + ")"
        jobid_macro = "$(jid)"
    else:
        queue_stmt = f"queue {njobs}"
        jobid_macro = "$(Process)"

    submit_description = f"""
executable = run.sh
arguments = {outdir} {GroupSize} {jobid_macro} {jetbin} {1 if same_sign else 0} {era} {tuplizer}
log = {scheduler_log}
JobBatchName = WtoLNu2Jets{suffix}_H_phipimu{ss_tag}{gt_tag}_{era}
request_cpus = 1
request_memory = 8G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
{queue_stmt}
"""
    what = f"refill {len(jobids)} job(s) {jobids}" if jobids else f"{njobs} jobs"
    print(f"{era} {jetbin}{ss_tag}{gt_tag} [{tuplizer}]: {NFiles} files / {GroupSize} per job "
          f"-> {what} -> {outdir}")
    subprocess.run(["condor_submit"], input=submit_description.encode())


if __name__ == "__main__":
    args = sys.argv[1:]
    # Validate before makedirs: an unrecognised token would otherwise silently create
    # a directory nothing reads and submit a few hundred jobs writing into it.
    jobids = None
    if "--jobs" in args:
        i = args.index("--jobs")
        if i + 1 >= len(args):
            sys.exit("--jobs needs a comma-separated JOBID list, e.g. --jobs 27,190,429")
        jobids = sorted({int(x) for x in args[i + 1].split(",") if x.strip()})
        del args[i:i + 2]
    gentag = "--gentag" in args
    if gentag:
        args.remove("--gentag")
    if len(args) < 2 or args[0] not in ERAS or args[1] not in JETBINS:
        sys.exit(f"usage: python3 condor_submit.py <{'|'.join(ERAS)}> <{'|'.join(JETBINS)}> "
                 f"[SS] [--gentag] [--jobs 27,190,429]")
    sys.exit(submit_jobs(args[0], args[1], len(args) > 2 and args[2].upper() == "SS",
                         jobids, gentag))

#!/usr/bin/env python3
"""Resubmit the background GenTag nano (OS+SS) in slots of <=200 jobs, 10 min apart.

2026-10-26: phi20 moved into the producer, so every background nano is redone. Only the
JOBIDs that exist in the old (renamed *_pre_phi20nano) directories are resubmitted -- same
file slices, and the inputs dropped earlier as bad stay dropped. The JOBID lists come from
a JSON {"<sample>|<year>": [jobids]} (identical for OS and SS, checked before writing it).

    python3 resubmit_bkg_batches.py <jobids.json> <state.txt> [--dryrun]

Resumable: every successfully submitted (sample, year, ss, jobid) is appended to the state
file and skipped on a rerun. Stops at the first failed condor_submit. After the first slot it
waits until a job of that cluster is running, then submits one slot every 10 minutes.
"""
import filecmp
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SLOT, PAUSE = 200, 600
ORDER = ["WJ0", "WJ1", "WJ2", "DY0J", "DY1J", "DY2J",
         "TTto2L2Nu", "TTto2L2Nu_ext1", "TTtoLNu2Q", "TTtoLNu2Q_ext1"]
GROUPS = [["WJ0", "WJ1", "WJ2"], ["DY0J", "DY1J", "DY2J"],
          ["TTto2L2Nu", "TTto2L2Nu_ext1", "TTtoLNu2Q", "TTtoLNu2Q_ext1"]]


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def main():
    jobids_json, state, dry = sys.argv[1], sys.argv[2], "--dryrun" in sys.argv
    ids = json.load(open(jobids_json))
    filelist_ref = sys.argv[3] if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else None
    opt = dict(a[2:].split("=", 1) for a in sys.argv if a.startswith("--") and "=" in a)
    only_samples = opt["samples"].split(",") if "samples" in opt else ORDER
    only_years = opt["years"].split(",") if "years" in opt else ["2022", "2022EE"]
    only_ss = [int(x) for x in opt["ss"].split(",")] if "ss" in opt else [0, 1]
    slot_size = int(opt.get("slot", SLOT))  # --slot=N overrides the 200-job slot (user's call)

    queue = []  # (sample, year, ss, jobid): W+jets, then DY, then ttbar; preEE before postEE
    for group in GROUPS:
        for year in ("2022", "2022EE"):
            if year not in only_years:
                continue
            for sample in [x for x in group if x in only_samples]:
                for ss in only_ss:
                    queue += [(sample, year, ss, j) for j in ids["%s|%s" % (sample, year)]]
    done = set()
    if os.path.exists(state):
        done = {ln.strip() for ln in open(state) if ln.strip()}
    todo = [q for q in queue if "%s|%s|%d|%d" % q not in done]
    log("queue %d, already submitted %d, to do %d" % (len(queue), len(queue) - len(todo), len(todo)))

    first_slot = len(done) == 0
    for s in range(0, len(todo), slot_size):
        slot = todo[s:s + slot_size]
        groups = []  # consecutive runs of the same (sample, year, ss)
        for q in slot:
            if groups and groups[-1][0] == q[:3]:
                groups[-1][1].append(q[3])
            else:
                groups.append((q[:3], [q[3]]))
        clusters = []
        for (sample, year, ss), jids in groups:
            cmd = ["python3", "%s/condor_submit.py" % HERE, sample, year, "--gentag",
                   "--jobids", ",".join(map(str, jids))] + (["--ss"] if ss else [])
            log("submit %s %s ss=%d: %d jobs (JOBID %d..%d)" % (sample, year, ss, len(jids), jids[0], jids[-1]))
            if dry:
                continue
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            print(p.stdout, flush=True)
            m = re.search(r"(\d+) job\(s\) submitted to cluster (\d+)", p.stdout)
            if p.returncode != 0 or not m or int(m.group(1)) != len(jids):
                log("SUBMIT FAILED (rc=%d) -- stopping" % p.returncode)
                sys.exit(1)
            if filelist_ref:
                fl = "%s_%s.txt" % (sample, year)
                if not filecmp.cmp("%s/filelists/%s" % (HERE, fl), "%s/%s" % (filelist_ref, fl), shallow=False):
                    log("FILELIST CHANGED: %s -- stopping" % fl)
                    sys.exit(1)
            clusters.append(m.group(2))
            with open(state, "a") as fh:
                for j in jids:
                    fh.write("%s|%s|%d|%d\n" % (sample, year, ss, j))
        if dry:
            continue
        log("slot done, clusters %s; submitted so far %d/%d" % (clusters, len(queue) - len(todo) + s + len(slot), len(queue)))
        if s + slot_size >= len(todo):
            break
        if first_slot:  # the first slot must be running before the 10-minute cadence starts
            first_slot = False
            while True:
                time.sleep(120)
                q = subprocess.run(["condor_q", clusters[0], "-af", "JobStatus"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, universal_newlines=True)
                if q.returncode == 0 and "2" in q.stdout.split():
                    log("first slot is running")
                    break
                log("first slot not running yet (condor_q rc=%d)" % q.returncode)
        time.sleep(PAUSE)
    log("all submitted" if not dry else "dryrun done")


if __name__ == "__main__":
    main()

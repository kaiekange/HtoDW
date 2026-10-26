"""Submit one phipimu BDT training+scoring job to HTCondor.

    python3 condor_submit.py ssdata
    python3 condor_submit.py wjets
    python3 condor_submit.py ssdata 101 --samples HtoDsW,HtoDsstW,Data_2022
    python3 condor_submit.py resdy --train-only
    python3 condor_submit.py resdy --predict --samples Signal_merged,Data_OS_merged

--train-only runs run_train.sh (BDT.py only, no scoring) instead of run.sh; scheduler log
condor_logs/HTCondor_train_BDT_<tag>.log. --predict scores an already-trained model only
(../DNN/run_predict.sh BDT ..., never retrains): ONE job per sample in --samples (required),
scheduler log condor_logs/HTCondor_predict_BDT_<tag>_<sample>.log. --region R (with --predict)
scores the control region into Prediction/<R>_<name>.root. Same as DNN/condor_submit.py.

One job per (background, seed) combination. Each gets its own scheduler log so
concurrent submissions do not interleave into the same file.

Must be submitted with the ds_dnn conda env active -- this imports common, which
imports torch (needed by the DNN half of that module, not by the BDT).
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNN"))
import common

USAGE = (f"usage: python3 {os.path.basename(__file__)} "
         f"<{'|'.join(common.TRAININGS)}> [seed] [--samples a,b,c] [--train-only | --predict]")

argv = sys.argv[1:]
TRAIN_ONLY = "--train-only" in argv
PREDICT = "--predict" in argv
argv = [x for x in argv if x not in ("--train-only", "--predict")]
REGION = None
if "--region" in argv:   # --predict only: score the control region (predict.py --region)
    i = argv.index("--region")
    REGION = argv[i + 1]
    argv = argv[:i] + argv[i + 2:]
    if not PREDICT:
        sys.exit("--region only makes sense with --predict")
if TRAIN_ONLY and PREDICT:
    sys.exit("--train-only and --predict are exclusive")
SAMPLES = ""
if "--samples" in argv:
    i = argv.index("--samples")
    if i + 1 >= len(argv):
        sys.exit(USAGE)
    SAMPLES = argv[i + 1]
    argv = argv[:i] + argv[i + 2:]

if len(argv) not in (1, 2):
    sys.exit(USAGE)
BKG = argv[0]
cfg = common.resolve(BKG)          # validate the token before creating anything
SEED = int(argv[1]) if len(argv) == 2 else common.SEED_DEFAULT

OUTDIR = common.out_dir_for(BKG, tool="BDT", seed=SEED)
os.makedirs(OUTDIR, exist_ok=True)

TAG = BKG if SEED == common.SEED_DEFAULT else f"{BKG}_seed{SEED}"
# The --samples pair is appended only when non-empty rather than passed as an empty
# positional: HTCondor's new-style argument syntax reads '' inside the quoted
# arguments string as an escaped literal quote, not as an empty argument.
ARGS = f"{BKG} {SEED}" + (f" --samples {SAMPLES}" if SAMPLES else "")
if TRAIN_ONLY and SAMPLES:
    sys.exit("--samples selects what to score; it has no meaning with --train-only")
EXE = "run_train.sh" if TRAIN_ONLY else "run.sh"
LOG = f"condor_logs/HTCondor_train_BDT_{TAG}.log" if TRAIN_ONLY else f"./HTCondor_{TAG}.log"


def submit_jobs():
    os.makedirs("condor_logs", exist_ok=True)
    submit_description = f"""
executable = {EXE}
arguments = "{ARGS}"
log = {LOG}
JobBatchName = phipimu_BDT_{"train_" if TRAIN_ONLY else ""}{TAG}
request_cpus = 4
request_memory = 24G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1
"""
    print(f"Submitting {BKG} ({cfg['label']}), seed={SEED}, {EXE}"
          f"{', samples=' + SAMPLES if SAMPLES else ''} -> {OUTDIR}")
    subprocess.run(["condor_submit"], input=submit_description.encode())


def submit_predict():
    samples = [x for x in SAMPLES.split(",") if x]
    if not samples:
        sys.exit("--predict needs --samples a,b,c")
    bad = [x for x in samples if x not in common.ALL_SAMPLES]
    if bad:
        sys.exit(f"--samples: unknown {bad}; choices {common.ALL_SAMPLES}")
    if not os.path.exists(f"{OUTDIR}/bdt_model.pkl"):
        sys.exit(f"no trained model in {OUTDIR}")
    os.makedirs("condor_logs", exist_ok=True)
    rarg = f" --region {REGION}" if REGION else ""
    rtag = f"{REGION}_" if REGION else ""
    for smp in samples:
        submit_description = f"""
executable = ../DNN/run_predict.sh
arguments = "BDT {BKG} {SEED} --samples {smp}{rarg}"
log = condor_logs/HTCondor_predict_BDT_{TAG}_{rtag}{smp}.log
JobBatchName = phipimu_BDT_predict_{TAG}_{rtag}{smp}
request_cpus = 8
request_memory = 24G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1
"""
        print(f"Submitting predict {BKG} seed={SEED} on {smp} -> {OUTDIR}/Prediction")
        subprocess.run(["condor_submit"], input=submit_description.encode())


if __name__ == "__main__":
    sys.exit(submit_predict() if PREDICT else submit_jobs())

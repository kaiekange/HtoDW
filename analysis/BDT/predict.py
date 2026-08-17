"""Score every current Samples/ file with a trained BDT and write a copy with
score_dsmu added.

    python3 predict.py ssdata                 # score with the SS-data-trained BDT
    python3 predict.py wjets                  # score with the W+jets-trained BDT
    python3 predict.py ssdata 101             # score the seed-101 BDT
    python3 predict.py ssdata --samples HtoDsW,HtoDsstW,Data_2022

The model directory is resolved from the same common.TRAININGS table the training
used, so the model and its Prediction/ output can never be crossed.

Both OS and SS variants are scored, and data as well as MC: the diagnostics need the
score distribution in the SS control region and in data sidebands, not just in the
samples the MVA was trained on.

The event selection (1 tight muon + b-jet veto) is deliberately NOT applied here --
every row gets a score and the selection is applied downstream, so a cut can be
revisited without re-scoring. Rows with non-finite features get score_dsmu = -777.
All input branches are copied through, so event_id survives into Prediction/ and
downstream plots can honour the training parity split.

The feature list and the preprocessing chain come from ../DNN/common.py, the same
module BDT.py trains with. They used to be a second hand-typed copy here, carrying a
comment claiming they matched DNN/DNN.py -- which was already false.
"""

import os
import sys

import numpy as np
import uproot
import joblib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNN"))
import common
from common import FEATURE_NAMES, SAMPLES_DIR, TREE

USAGE = (f"usage: python3 {os.path.basename(__file__)} "
         f"<{'|'.join(common.TRAININGS)}> [seed] [--samples a,b,c]")

# Current Samples/ contents. The old inclusive WtoLNu2Jets{,_SS} are deliberately
# left out -- superseded by the 0J/1J/2J jet-binned samples. Identical to the list in
# DNN/predict.py: the two must score the same files for the head-to-head to work.
ALL_SAMPLES = [
    "HtoDsW", "HtoDsW_SS",
    "HtoDsstW", "HtoDsstW_SS",
    "WtoLNu2Jets_0J", "WtoLNu2Jets_1J", "WtoLNu2Jets_2J",
    "WtoLNu2Jets_0J_SS", "WtoLNu2Jets_1J_SS", "WtoLNu2Jets_2J_SS",
    "Data_2022", "Data_2022_SS",
]

argv = sys.argv[1:]
SAMPLES = ALL_SAMPLES
if "--samples" in argv:
    i = argv.index("--samples")
    if i + 1 >= len(argv):
        sys.exit(USAGE)
    SAMPLES = [s for s in argv[i + 1].split(",") if s]
    unknown = [s for s in SAMPLES if s not in ALL_SAMPLES]
    if unknown:
        sys.exit(f"unknown sample(s): {', '.join(unknown)}\n"
                 f"known: {', '.join(ALL_SAMPLES)}")
    argv = argv[:i] + argv[i + 2:]

if len(argv) not in (1, 2):
    sys.exit(USAGE)
BKG = argv[0]
cfg = common.resolve(BKG)
SEED = int(argv[1]) if len(argv) == 2 else common.SEED_DEFAULT

MODEL_DIR = common.out_dir_for(BKG, tool="BDT", seed=SEED)
OUT_DIR = f"{MODEL_DIR}/Prediction"
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Loading model from {MODEL_DIR} (background={BKG}, seed={SEED})...")
if SAMPLES is not ALL_SAMPLES:
    print(f"  scoring subset only: {', '.join(SAMPLES)}")

# Refuses to proceed if this dir's feature stamp disagrees with common.py -- the
# alternative is silently wrong scores, since feature order is only implicitly bound
# to the tree splits.
common.check_feature_stamp(MODEL_DIR)
model = joblib.load(f"{MODEL_DIR}/bdt_model.pkl")
common.check_n_features(MODEL_DIR, model.n_features_in_, "bdt_model.pkl")

for name in SAMPLES:
    path = f"{SAMPLES_DIR}/{name}.root"
    if not os.path.exists(path):
        print(f"  {name:24s} MISSING -- skipped")
        continue

    with uproot.open(path)[TREE] as t:
        data = t.arrays(library="np")

    X = np.column_stack([data[b].astype(np.float64) for b in FEATURE_NAMES])
    n = len(X)

    X, mask = common.preprocess(X)
    scores = np.full(n, -777.0, dtype=np.float32)
    scores[mask] = model.predict_proba(X[mask])[:, 1].astype(np.float32)

    data["score_dsmu"] = scores
    out_file = f"{OUT_DIR}/{name}.root"
    with uproot.recreate(out_file) as f:
        f.mktree(TREE, {k: v.dtype for k, v in data.items()})
        f[TREE].extend(data)

    print(f"  {name:24s} {n:>9,} events, {int((~mask).sum()):>5} non-finite -> {out_file}")

print("Done.")

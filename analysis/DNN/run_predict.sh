#!/bin/bash
#   run_predict.sh <DNN|BDT> <training key> [seed] --samples a,b,c
#
# Scoring ONLY, with an already-trained model -- for adding a new sample (e.g.
# TT_merged) to an existing model's Prediction/ dir. run.sh is not usable for that:
# it retrains first and would overwrite the model.
#
# A stage-2 key reads its stage-1 model's Prediction/ file, so stage-1 must have
# scored the sample first (predict.py exits otherwise).

set -euo pipefail

TOOL=${1:?usage: run_predict.sh <DNN|BDT> <key> [seed] --samples a,b,c}
KEY=${2:?usage: run_predict.sh <DNN|BDT> <key> [seed] --samples a,b,c}
shift 2

source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn

ANALYSIS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis
OUTDIR=$(python3 -c "import sys; sys.path.insert(0,'${ANALYSIS}/DNN'); import common; print(common.out_dir_for('${KEY}', tool='${TOOL}'))")

python3 ${ANALYSIS}/${TOOL}/predict.py ${KEY} "$@" &> ${OUTDIR}/predict_$(date +%Y%m%d_%H%M%S).log

#!/bin/bash
#   run_train.sh <training key> [seed]
#
# Training ONLY -- the counterpart of run_predict.sh. run.sh trains and then scores
# every sample in the same job; this stops after DNN.py so scoring can be submitted
# separately, one condor job per sample (faster than one job scoring them in series).

set -euo pipefail

BKG=${1:?usage: run_train.sh <training key> [seed]}
SEED=${2:-42}

source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn

SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN

OUTDIR=$(python3 -c "import sys; sys.path.insert(0,'${SCRIPTS}'); import common; print(common.out_dir_for('${BKG}', tool='DNN', seed=${SEED}))")
mkdir -p "${OUTDIR}"

python3 ${SCRIPTS}/DNN.py ${BKG} ${SEED} &> ${OUTDIR}/training.log

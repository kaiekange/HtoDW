#!/bin/bash
#   score_only.sh <bkg> [seed] [--samples a,b,c]
#
# predict.py ONLY -- no retraining. run.sh always retrains first, which would
# overwrite an existing canonical model; this is for re-scoring / back-filling
# Prediction/ with a model that must not change (e.g. stage-1 for the stage-2 build).
set -euo pipefail

BKG=${1:?usage: score_only.sh <bkg> [seed] [--samples a,b,c]}
SEED=${2:-42}
shift $(( $# > 2 ? 2 : $# ))
PREDICT_ARGS=("$@")

source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn

SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN
OUTDIR=$(python3 -c "import sys; sys.path.insert(0,'${SCRIPTS}'); import common; print(common.out_dir_for('${BKG}', tool='DNN', seed=${SEED}))")
mkdir -p "${OUTDIR}"

python3 ${SCRIPTS}/predict.py ${BKG} ${SEED} ${PREDICT_ARGS[@]+"${PREDICT_ARGS[@]}"} \
    &> ${OUTDIR}/score_only.log

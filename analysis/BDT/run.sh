#!/bin/bash
#   run.sh <ssdata|wjets> [seed] [--samples a,b,c]
#
# Trains the phipimu BDT against the named background and then scores every Samples/
# file with it. Both steps resolve the output directory from the same
# common.TRAININGS table, so they cannot end up pointing at different models.
#
#   ssdata  signal(OS) vs Data_2022_SS (even event_id half; odd half held out)
#   wjets   signal(OS) vs WtoLNu2Jets 0J+1J+2J (lumi-reweighted)
#
# Mirrors DNN/run.sh exactly, on purpose: as of 2026-08-14 BDT and DNN share the
# feature set, preprocessing, selection, signal mixture, class weighting, split and
# seed via ../DNN/common.py, so the only difference between the two is the algorithm.

set -euo pipefail

BKG=${1:?usage: run.sh <ssdata|wjets> [seed] [--samples a,b,c]}
SEED=${2:-42}
shift $(( $# > 2 ? 2 : $# ))
PREDICT_ARGS=("$@")

source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn

SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/BDT
COMMON=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN

# Resolve (and validate) the output dir before redirecting the log into it.
OUTDIR=$(python3 -c "import sys; sys.path.insert(0,'${COMMON}'); import common; print(common.out_dir_for('${BKG}', tool='BDT', seed=${SEED}))")
mkdir -p "${OUTDIR}"

{
    python3 ${SCRIPTS}/BDT.py     ${BKG} ${SEED}
    # No best_cand step: the Samples tuples are already one row per event.
    # ${a[@]+"${a[@]}"} rather than "${a[@]}": under `set -u`, expanding an empty
    # array the plain way is an unbound-variable error on bash < 4.4.
    python3 ${SCRIPTS}/predict.py ${BKG} ${SEED} ${PREDICT_ARGS[@]+"${PREDICT_ARGS[@]}"}
} &> ${OUTDIR}/training.log

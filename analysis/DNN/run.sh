#!/bin/bash
#   run.sh <ssdata|wjets>
#
# Trains the phipimu DNN against the named background and then scores every
# Samples/ file with it. Both steps resolve the output directory from the same
# common.TRAININGS table, so they cannot end up pointing at different models.
#
#   ssdata  signal(OS) vs Data_2022_SS (even event_id half; odd half held out)
#   wjets   signal(OS) vs WtoLNu2Jets 0J+1J+2J (lumi-reweighted)
#
# The two runs differ ONLY in the background -- same 43 features, same seed, same
# network -- which is what makes them comparable. See DNN.py's docstring.

set -euo pipefail

BKG=${1:?usage: run.sh <ssdata|wjets>}

source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn

SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN

# Resolve (and validate) the output dir before redirecting the log into it.
OUTDIR=$(python3 -c "import sys; sys.path.insert(0,'${SCRIPTS}'); import common; print(common.resolve('${BKG}')['out_dir'])")
mkdir -p "${OUTDIR}"

{
    python3 ${SCRIPTS}/DNN.py     ${BKG}
    # No best_cand step: the Samples tuples are already one row per event.
    python3 ${SCRIPTS}/predict.py ${BKG}
} &> ${OUTDIR}/training.log

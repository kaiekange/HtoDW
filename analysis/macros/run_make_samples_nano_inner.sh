#!/bin/bash
# Runs inside the cmssw-el9 container (invoked by run_make_samples_nano.sh).
set -euo pipefail

LIST=${1:?usage: run_make_samples_nano_inner.sh <list.txt> <out.root>}
OUT=${2:?usage: run_make_samples_nano_inner.sh <list.txt> <out.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${REPO}"
eval `scramv1 runtime -sh`
python3 analysis/macros/make_samples_nano.py worker "${LIST}" "${OUT}" ${3:-}

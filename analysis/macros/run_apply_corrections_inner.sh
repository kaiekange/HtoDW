#!/bin/bash
# Runs inside the cmssw-el9 container (invoked by run_apply_corrections.sh).
set -euo pipefail

INFILE=${1:?usage: run_apply_corrections_inner.sh <nano output_N.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${REPO}"
eval `scramv1 runtime -sh`
python3 analysis/macros/apply_corrections.py "${INFILE}"

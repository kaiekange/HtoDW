#!/bin/bash
# Runs inside the cmssw-el9 container (invoked by run_apply_corrections.sh). CMSSW_15_0_5 is
# used only as a Python stack (correctionlib 2.6.4, uproot, awkward, scipy, ROOT); no CMSSW
# code is involved.
set -euo pipefail

INFILE=${1:?usage: run_apply_corrections_inner.sh <nano output_N.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/Corrections

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
eval `scramv1 runtime -sh`
python3 "${REPO}/corrections/apply_corrections.py" "${INFILE}"

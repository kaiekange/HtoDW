#!/bin/bash
# Runs inside the cmssw-el9 container (invoked by run_btag_eff.sh).
set -euo pipefail

GROUP=${1:?usage: run_btag_eff_inner.sh <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${REPO}"
eval `scramv1 runtime -sh`
python3 analysis/macros/build_btag_efficiency_map.py "${GROUP}"

#!/bin/bash
# HTCondor executable for build_btag_efficiency_map.py -- one job per sample group.
#   run_btag_eff.sh <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>
# Same structure as run_corrections.sh: the work runs inside cmssw-el9 (correctionlib
# 2.6.4) via a separate inner script file, see run_corrections_inner.sh for why.
set -euo pipefail

GROUP=${1:?usage: run_btag_eff.sh <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
mkdir -p "${REPO}/analysis/macros/condor_logs"

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/analysis/macros/run_btag_eff_inner.sh" "${GROUP}" \
    &> "${REPO}/analysis/macros/condor_logs/build_btag_efficiency_map_${GROUP}.log"

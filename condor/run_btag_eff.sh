#!/bin/bash
# HTCondor executable for build_btag_efficiency_map.py -- one job per sample group.
#   run_btag_eff.sh <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>
set -euo pipefail

GROUP=${1:?usage: run_btag_eff.sh <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/Corrections
mkdir -p "${REPO}/condor/condor_logs"

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/condor/run_btag_eff_inner.sh" "${GROUP}" \
    &> "${REPO}/condor/condor_logs/build_btag_efficiency_map_${GROUP}.log"

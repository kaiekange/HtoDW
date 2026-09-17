#!/bin/bash
# HTCondor executable for apply_corrections.py -- one job per nano output file.
#   run_apply_corrections.sh <.../NanoMVA/MC/<era>/<sample dir>/output_N.root>
# The work runs inside cmssw-el9 via a separate inner script file (cmssw-el9 re-quotes an
# inline `bash -c` string through two shell layers, which mangles multi-line commands).
set -euo pipefail

INFILE=${1:?usage: run_apply_corrections.sh <nano output_N.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/Corrections
if [[ "${INFILE}" == */NanoMVA/Data/* ]]; then
    # CRAB layout .../Data/<era>/Muon/<task>/<timestamp>/<NNNN>/nano_N.root
    SAMPLE=$(basename "$(dirname "$(dirname "$(dirname "${INFILE}")")")")
    ERA=data_$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$(dirname "${INFILE}")")")")")")
else
    ERA=$(basename "$(dirname "$(dirname "${INFILE}")")")
    SAMPLE=$(basename "$(dirname "${INFILE}")")
fi
STEM=$(basename "${INFILE}" .root)
LOGDIR="${REPO}/condor/condor_logs/apply_corrections"
mkdir -p "${LOGDIR}"

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/condor/run_apply_corrections_inner.sh" "${INFILE}" \
    &> "${LOGDIR}/${ERA}_${SAMPLE}_${STEM}.log"

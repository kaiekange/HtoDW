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
LOG="${LOGDIR}/${ERA}_${SAMPLE}_${STEM}.log"
mkdir -p "${LOGDIR}"

# The log is kept only for a failed job (2026-11-04: ~2e4 per-job logs filled the AFS
# directory, new files failed with "File too large"). It is written on the node's scratch
# disk and copied to LOGDIR on failure; on success any older log of this file is removed.
LOCAL_LOG="${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/${ERA}_${SAMPLE}_${STEM}.log"
rc=0
/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/condor/run_apply_corrections_inner.sh" "${INFILE}" \
    &> "${LOCAL_LOG}" || rc=$?
if [[ ${rc} -ne 0 ]]; then
    cp "${LOCAL_LOG}" "${LOG}"
    exit "${rc}"
fi
rm -f "${LOG}" "${LOCAL_LOG}"

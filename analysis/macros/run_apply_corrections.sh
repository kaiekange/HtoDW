#!/bin/bash
# HTCondor executable for apply_corrections.py -- one job per nano output file.
#   run_apply_corrections.sh <.../NanoMVA/MC/<era>/<sample dir>/output_N.root>
# Same structure as run_btag_eff.sh: the work runs inside cmssw-el9 (correctionlib 2.6.4)
# via a separate inner script file, see run_corrections_inner.sh for why.
set -euo pipefail

INFILE=${1:?usage: run_apply_corrections.sh <nano output_N.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
ERA=$(basename "$(dirname "$(dirname "${INFILE}")")")
SAMPLE=$(basename "$(dirname "${INFILE}")")
STEM=$(basename "${INFILE}" .root)
LOGDIR="${REPO}/analysis/macros/condor_logs/apply_corrections"
mkdir -p "${LOGDIR}"

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/analysis/macros/run_apply_corrections_inner.sh" "${INFILE}" \
    &> "${LOGDIR}/${ERA}_${SAMPLE}_${STEM}.log"

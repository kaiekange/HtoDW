#!/bin/bash
# HTCondor executable for make_samples_nano.py worker -- one job per chunk of nano files.
#   run_make_samples_nano.sh <chunk list.txt> <skim chunk_K.root> [--correct] [--region=R]
# Same structure as run_apply_corrections.sh: the work runs inside cmssw-el9.
set -euo pipefail

LIST=${1:?usage: run_make_samples_nano.sh <list.txt> <out.root>}
OUT=${2:?usage: run_make_samples_nano.sh <list.txt> <out.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
LOGDIR="${REPO}/analysis/macros/condor_logs/make_samples_nano"
mkdir -p "${LOGDIR}"
# CR jobs reuse the SR chunk lists: suffix the log with the region so they do not overwrite.
REGION=$(printf "%s\n" "$@" | sed -n "s/^--region=//p")

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/analysis/macros/run_make_samples_nano_inner.sh" "${LIST}" "${OUT}" "${@:3}" \
    &> "${LOGDIR}/$(basename "${LIST}" .txt)${REGION:+_${REGION}}.log"

#!/bin/bash
# HTCondor executable for jme_zmumu_study.py worker -- one job per zmumu skim chunk.
#   run_jme_zmumu_study.sh <chunk list.txt> <out chunk_K.root>
# Same structure as run_make_samples_nano.sh: the work runs inside cmssw-el9.
set -euo pipefail

LIST=${1:?usage: run_jme_zmumu_study.sh <list.txt> <out.root>}
OUT=${2:?usage: run_jme_zmumu_study.sh <list.txt> <out.root>}
REPO=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
LOGDIR="${REPO}/analysis/macros/condor_logs/jme_zmumu_study"
mkdir -p "${LOGDIR}"

/cvmfs/cms.cern.ch/common/cmssw-el9 -- "${REPO}/analysis/macros/run_jme_zmumu_study_inner.sh" "${LIST}" "${OUT}" \
    &> "${LOGDIR}/$(basename "${LIST}" .txt).log"

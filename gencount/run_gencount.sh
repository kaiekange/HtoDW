#!/bin/bash
# HTCondor executable for gencount.py (lxplus8 / el8, CMSSW_13_0_13 FWLite).
#   run_gencount.sh <DsW|DsstW> <era> <filelist> <JOBID> <GROUPSIZE>
set -uo pipefail
CHANNEL=$1; ERA=$2; FILELIST=$3; JOBID=$4; GROUPSIZE=$5
export X509_USER_PROXY=/afs/cern.ch/user/k/kakang/tmp/x509up_u158898
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
OUTDIR=/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/${ERA}/Corrections/Signal_${CHANNEL}_H_phipimu/gencount
mkdir -p "${OUTDIR}/log"
OUT="${OUTDIR}/count_${JOBID}.json"
LOG="${OUTDIR}/log/run_${JOBID}.log"
rm -f "${OUT}" "${LOG}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src
eval `scramv1 runtime -sh`
python3 /afs/cern.ch/work/k/kakang/H2DW/Corrections/gencount/gencount.py \
    "${CHANNEL}" "${FILELIST}" "${JOBID}" "${GROUPSIZE}" "${OUT}" &> "${LOG}"
STATUS=$?
exit ${STATUS}

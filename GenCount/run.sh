#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3
CHANNEL=$4

OUTFILE="${OUTDIR}/count_${JOBID}.txt"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f $OUTFILE $LOGFILE

SCRIPT="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/GenCount/gencount.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
python3 ${SCRIPT} ${CHANNEL} ${JOBID} ${GROUPSIZE} "${OUTFILE}" &> ${LOGFILE}

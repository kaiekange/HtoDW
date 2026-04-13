#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

FILELIST=$1
OUTDIR=$2
GROUPSIZE=$3
PROCESS=$4

STARTLINE=$(( PROCESS * GROUPSIZE + 1 ))
ENDLINE=$(( STARTLINE + GROUPSIZE - 1 ))

INFILE="$(sed -n "${STARTLINE}, ${ENDLINE}p" "${FILELIST}")"

echo $INFILE

OUTFILE="${OUTDIR}/output_$(( PROCESS + 1 )).root"
LOGFILE="${OUTDIR}/log/output_$(( PROCESS + 1 )).log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/FullReco.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} INFILE="${INFILE}" OUTFILE="${OUTFILE}" &> ${LOGFILE}

#!/bin/bash

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3

OUTFILE="${OUTDIR}/output_${JOBID}.root"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/SignalReco.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
cmsRun ${CFGFILE} JOBID=${JOBID} GROUPSIZE=${GROUPSIZE} OUTFILE="${OUTFILE}" &> ${LOGFILE}

#!/bin/bash

OUTDIR=$1
# JOBID=$(( $2 + 1 ))
JOBID=$2

outfile="file:${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
rm -f $outfile $logfile

cfgfile="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/RecoKineFit_test.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
cmsRun ${cfgfile} ${JOBID} ${outfile} &> ${logfile}

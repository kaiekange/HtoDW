#!/bin/bash

INDIR=$1
OUTDIR=$2
JOBID=$(( $3 + 1000 ))

infile="file:${INDIR}/output_${JOBID}.root"
outfile="file:${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
rm -f $outfile $logfile

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/Simulation/PAT/PAT_cfg.py" ${infile} ${outfile} &> ${logfile}

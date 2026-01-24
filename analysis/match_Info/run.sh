#!/bin/bash

OUTDIR=$1
JOBID=$(( $2 + 1 ))

infile="file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_${JOBID}.root"
outfile="file:${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
rm -f $outfile $logfile

cfgfile="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/TrigAnalyzer/test/match_Info.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
cmsRun ${cfgfile} ${infile} ${outfile} &> ${logfile}

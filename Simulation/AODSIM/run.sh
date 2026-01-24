#!/bin/bash

INDIR=$1
OUTDIR=$2
JOBID=$(( $3 + 185 ))

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

infile="file:${INDIR}/output_${JOBID}.root"
outfile="file:${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
rm -f $outfile $logfile

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/Simulation/AODSIM/AODSIM_cfg.py" ${infile} ${outfile} &> ${logfile}

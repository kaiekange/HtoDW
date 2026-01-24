#!/bin/bash

OUTDIR=$1
GROUPID=$2
PROCESS=$3
JOBID=$(( $GROUPID + $PROCESS ))

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

TIMESTAMP=$(date +%s)
SEED=$(( ( $JOBID * 100000 + $TIMESTAMP % 100000) % 900000000 + 1 ))

outfile="file:${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
rm -f $outfile $logfile

workdir="${OUTDIR}/workdir/WD_${JOBID}"
rm -rf $workdir
mkdir $workdir

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cd $workdir
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/Simulation/ggHtoDsW_cfg.py" ${outfile} ${SEED} &> ${logfile}

rm -rf $workdir

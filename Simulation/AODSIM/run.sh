#!/bin/bash

MODE=$1
INDIR=$2
OUTDIR=$3
JOBID=$4

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

infile="${INDIR}/output_${JOBID}.root"
outfile="${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
xrdfs maite.iihe.ac.be rm -f $outfile
xrdfs maite.iihe.ac.be rm -f $logfile

tmplog="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/Simulation/AODSIM/tmp/run_${MODE}_${JOBID}.log"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/Simulation/AODSIM/AODSIM_${MODE}_cfg.py" "root://maite.iihe.ac.be/${infile}" "root://maite.iihe.ac.be/${outfile}" &> ${tmplog}

xrdcp -f ${tmplog} root://maite.iihe.ac.be/${logfile}
rm -f ${tmplog}

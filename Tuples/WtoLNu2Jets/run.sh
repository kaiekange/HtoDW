#!/bin/bash
# HTCondor executable for every WtoLNu2Jets tuplizer job: jet bin and OS/SS are
# arguments now, replacing the eight copies of this script there used to be.
#   run.sh <outdir> <groupsize> <jobid> <jetbin: incl|0J|1J|2J> <ss: 0|1>

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3
JETBIN=$4
SS=$5

OUTFILE="${OUTDIR}/output_${JOBID}.root"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/WtoLNu2Jets_H_phipimu.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} JOBID=${JOBID} GROUPSIZE=${GROUPSIZE} JETBIN=${JETBIN} SS=${SS} OUTFILE="${OUTFILE}" &> ${LOGFILE}

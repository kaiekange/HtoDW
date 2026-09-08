#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"
# The worker's own /etc/grid-security/certificates carries 883 CAs, the CVMFS bundle 1050.
# Without this, xrootd can fail to verify the SERVER's certificate and reports it as
# "Auth failed: No protocols left to try", which reads like a proxy problem and is not one.
# Same line already in WtoLNu2Jets/run.sh and the PAT and AODSIM run.sh.
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3

OUTFILE="${OUTDIR}/output_${JOBID}.root"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/SIG_DsstW_H_phipimu_SS.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} JOBID=${JOBID} GROUPSIZE=${GROUPSIZE} OUTFILE="${OUTFILE}" &> ${LOGFILE}

#!/bin/bash
# HTCondor executable for every WtoLNu2Jets tuplizer job: jet bin and OS/SS are
# arguments now, replacing the eight copies of this script there used to be.
#   run.sh <outdir> <groupsize> <jobid> <jetbin: incl|0J|1J|2J> <ss: 0|1> <era: 2022|2022EE>
#          [tuplizer: H_phipimu|H_phipimu_GenTag]
#
# The 7th argument is optional and defaults to H_phipimu, so submissions that pass six
# arguments behave exactly as before.

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"
# The worker's own /etc/grid-security/certificates carries 883 CAs, the CVMFS bundle 1050.
# Without this, xrootd can fail to verify the SERVER's certificate and reports it as
# "Auth failed: No protocols left to try", which reads like a proxy problem and is not one.
# Same line already in the PAT and AODSIM run.sh.
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3
JETBIN=$4
SS=$5
ERA=$6
TUPLIZER=${7:-H_phipimu}

OUTFILE="${OUTDIR}/output_${JOBID}.root"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/WtoLNu2Jets_H_phipimu.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} JOBID=${JOBID} GROUPSIZE=${GROUPSIZE} JETBIN=${JETBIN} SS=${SS} ERA=${ERA} TUPLIZER=${TUPLIZER} OUTFILE="${OUTFILE}" &> ${LOGFILE}

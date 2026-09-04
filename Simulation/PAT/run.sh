#!/bin/bash

MODE=$1
INDIR=$2
OUTDIR=$3
JOBID=$4

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

# The worker nodes' own CA bundle (/etc/grid-security/certificates, 883 certs) cannot
# verify the IIHE server certificate, so GSI fails with
#   Secgsi: ErrParseBuffer: unknown CA: cannot verify server certificate
# and xrootd reports "Auth failed: No protocols left to try". That message names the
# client credential, but the proxy is fine -- it is the server cert that cannot be
# checked. The CVMFS bundle (1050 certs) has the missing CA.
# Verified 2026-09-02 on one worker node: the same xrdfs command failed with the
# default CA dir and listed the directory with this one; a full PAT job then ran to
# completion (exit 0, 166 MB output, log uploaded).
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates

infile="${INDIR}/output_${JOBID}.root"
outfile="${OUTDIR}/output_${JOBID}.root"
logfile="${OUTDIR}/log/run_${JOBID}.log"
xrdfs maite.iihe.ac.be rm -f $outfile
xrdfs maite.iihe.ac.be rm -f $logfile

tmplog="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/Simulation/PAT/tmp/run_${MODE}_${JOBID}.log"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/Simulation/PAT/PAT_${MODE}_cfg.py" "root://maite.iihe.ac.be/${infile}" "root://maite.iihe.ac.be/${outfile}" &> ${tmplog}
CMSRUN_STATUS=$?

xrdcp -f ${tmplog} root://maite.iihe.ac.be/${logfile}
rm -f ${tmplog}

exit ${CMSRUN_STATUS}

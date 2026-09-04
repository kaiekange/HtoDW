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
# Verified 2026-09-02 on the PAT side (same run.sh structure); see PAT/run.sh.
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates

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
CMSRUN_STATUS=$?

# Delete the local copy only once the upload has succeeded. Deleting it
# unconditionally destroys the evidence in exactly the case where it is needed: a job
# that cannot reach the storage writes no output AND uploads no log, so the only
# record of why it died was this file. That is how 147 AODSIM failures on 2026-09-02
# ended up undiagnosable. A leftover file here means "this job failed, look at me".
xrdcp -f ${tmplog} root://maite.iihe.ac.be/${logfile} && rm -f ${tmplog}

exit ${CMSRUN_STATUS}

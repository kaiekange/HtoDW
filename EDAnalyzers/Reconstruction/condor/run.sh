#!/bin/bash
# Custom NanoAOD production worker. One output file per job.
#   run.sh <OUTDIR> <GROUPSIZE> <JOBID> <SAMPLE> <YEAR> <SS> <TUPLIZER> <FILELIST>

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"
# Worker /etc/grid-security has fewer CAs than CVMFS; without this xrootd can fail to verify
# the SERVER cert and report "Auth failed: No protocols left to try" (not a proxy problem).
export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
export SCRAM_ARCH=el8_amd64_gcc11

OUTDIR=$1
GROUPSIZE=$2
JOBID=$3
SAMPLE=$4
YEAR=$5
SS=$6
TUPLIZER=$7
FILELIST=$8

OUTFILE="${OUTDIR}/output_${JOBID}.root"
LOGFILE="${OUTDIR}/log/run_${JOBID}.log"
rm -f "$OUTFILE" "$LOGFILE"

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/NANO_H_phipimu.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
eval `scramv1 runtime -sh`
voms-proxy-info -all -file "$X509_USER_PROXY"

cmsRun "${CFGFILE}" SAMPLE="${SAMPLE}" YEAR="${YEAR}" JOBID="${JOBID}" GROUPSIZE="${GROUPSIZE}" \
    SS="${SS}" TUPLIZER="${TUPLIZER}" NTHREADS=2 FILELIST="${FILELIST}" \
    OUTFILE="file:${OUTFILE}" &> "${LOGFILE}"
CMSRUN_STATUS=$?

# Keep the (possibly partial) output on failure as evidence -- do not rm.
exit $CMSRUN_STATUS

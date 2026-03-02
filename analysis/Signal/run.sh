#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/Signal"

INFILE="file_list = ['file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_{}.root'.format(i) for i in range(1, 11) if i != 185]"
OUTFILE="file:${OUTDIR}/test.root"

LOGFILE="${OUTDIR}/test.log"
rm -f $OUTFILE $LOGFILE

cfgfile="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/FullReco_test.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${cfgfile} ${INFILE} ${OUTFILE} &> ${LOGFILE}

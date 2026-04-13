#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/ZZto2L2Q"

INFILE="/store/mc/Run3Summer22MiniAODv4/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/0ebb2c98-1462-4c76-888d-60bcc2476f4e.root"
OUTFILE="${OUTDIR}/test.root"
LOGFILE="${OUTDIR}/test.log"

rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/FullReco.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} INFILE="${INFILE}" OUTFILE="${OUTFILE}" &> ${LOGFILE}
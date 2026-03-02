#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/TW"

INFILE="/store/mc/Run3Summer22MiniAODv4/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/239b7947-60df-430e-a5a8-94a067482b0a.root"
OUTFILE="file:${OUTDIR}/test.root"

LOGFILE="${OUTDIR}/test.log"
rm -f $OUTFILE $LOGFILE

cfgfile="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/FullReco_test.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${cfgfile} ${INFILE} ${OUTFILE} &> ${LOGFILE}

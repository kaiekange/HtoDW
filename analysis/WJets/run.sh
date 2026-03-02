#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/WJets"

INFILE="/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/03d57826-8e8d-4b11-8683-c75d0310b16f.root"
OUTFILE="file:${OUTDIR}/test.root"

LOGFILE="${OUTDIR}/test.log"
rm -f $OUTFILE $LOGFILE

cfgfile="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/FullReco_test.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${cfgfile} ${INFILE} ${OUTFILE} &> ${LOGFILE}

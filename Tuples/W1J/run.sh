#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/W1J"

INFILE="/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2540000/06118768-c07c-43c4-9aec-790384565132.root"
OUTFILE="${OUTDIR}/test.root"
LOGFILE="${OUTDIR}/test.log"
rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/FullReco.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} INFILE="${INFILE}" OUTFILE="${OUTFILE}" &> ${LOGFILE}

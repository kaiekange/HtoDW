#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

OUTDIR="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WWto2L2Nu_SS_EW"

INFILE="/store/mc/Run3Summer22MiniAODv4/WWto2L2Nu-2Jets_SS_noTop_EW_TuneCP5_13p6TeV_madgraph-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v3/2530000/a596b84b-253e-4440-9086-3ab5fe36d60a.root"
OUTFILE="${OUTDIR}/test.root"
LOGFILE="${OUTDIR}/test.log"

rm -f $OUTFILE $LOGFILE

CFGFILE="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/FullReco.py"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CFGFILE} INFILE="${INFILE}" OUTFILE="${OUTFILE}" &> ${LOGFILE}

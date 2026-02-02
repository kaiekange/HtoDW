#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/AllRecoTest.py" &> "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/RecoCompare/test.log"

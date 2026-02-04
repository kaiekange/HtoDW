#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/
cmsenv
cmsRun "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/TrigAnalyzer/test/KineFit_test.py" &> "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/KineFit/test.log"

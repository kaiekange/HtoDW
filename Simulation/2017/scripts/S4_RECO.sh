#!/bin/bash

export HOME=/user/kakang

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

date=$1
process=$2
events=$3
nThreads=$4

cd ../$date/samples/process_$process

stepname=RECO
cfgfile=cfg/S4_${stepname}_cfg.py
infile=output/S3_HLT.root
outfile=output/S4_${stepname}.root
logfile=log/S4_${stepname}.log

{
    cmsDriver.py \
        --python_filename $cfgfile \
        --eventcontent AODSIM \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier AODSIM \
        --filein file:$infile \
        --fileout file:$outfile \
        --conditions 106X_mc2017_realistic_v6 \
        --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
        --beamspot Realistic25ns13TeVEarly2017Collision \
        --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI \
        --geometry DB:Extended \
        --era Run2_2017 \
        --no_exec \
        --mc \
        --nThreads=$nThreads \
        -n $events || exit $? ;
    
    cmsRun $cfgfile
} &> $logfile

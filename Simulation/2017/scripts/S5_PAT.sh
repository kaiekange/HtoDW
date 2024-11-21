#!/bin/bash

export HOME=/user/kakang

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

date=$1
process=$2
events=$3
nThreads=$4

cd ../$date/samples/process_$process

stepname=PAT
cfgfile=cfg/S5_${stepname}_cfg.py
infile=output/S4_RECO.root
outfile=output/S5_${stepname}.root
logfile=log/S5_${stepname}.log

{
    cmsDriver.py \
        --python_filename $cfgfile \
        --eventcontent MINIAODSIM \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier MINIAODSIM \
        --filein file:$infile \
        --fileout file:$outfile \
        --conditions 106X_mc2017_realistic_v6 \
        --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
        --beamspot Realistic25ns13TeVEarly2017Collision \
        --step PAT \
        --procModifiers run2_miniAOD_UL \
        --geometry DB:Extended \
        --era Run2_2017 \
        --runUnscheduled \
        --no_exec \
        --mc \
        --nThreads=$nThreads \
        -n $events || exit $? ;
    
    cmsRun $cfgfile
} &> $logfile

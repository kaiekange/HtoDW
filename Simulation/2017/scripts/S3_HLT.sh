#!/bin/bash

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

date=$1
process=$2
events=$3
nThreads=$4

cd ../$date/samples/process_$process

stepname=HLT
cfgfile=cfg/S3_${stepname}_cfg.py
infile=output/S2_DIGI.root
outfile=output/S3_${stepname}.root
logfile=log/S3_${stepname}.log

{
    cmsDriver.py \
        --python_filename $cfgfile \
        --eventcontent RAWSIM \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier GEN-SIM-RAW \
        --filein file:$infile \
        --fileout file:$outfile \
        --conditions 106X_mc2017_realistic_v6 \
        --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
        --beamspot Realistic25ns13TeVEarly2017Collision \
        --step HLT:@relval2017 \
        --geometry DB:Extended \
        --era Run2_2017 \
        --no_exec \
        --mc \
        --nThreads=$nThreads \
        -n $events || exit $? ;
    
    cmsRun $cfgfile
} &> $logfile

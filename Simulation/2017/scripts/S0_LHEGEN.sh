#!/bin/bash

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

date=$1
process=$2
events=$3
nThreads=$4

mkdir -p ../$date/samples/process_$process
cd ../$date/samples/process_$process
mkdir -p log cfg output

stepname=LHEGEN
cfgfile=cfg/S0_${stepname}_cfg.py
outfile=output/S0_${stepname}.root
logfile=log/S0_${stepname}.log

{
    cmsDriver.py Configuration/GenProduction/python/gridpacktogen_ggHccDW_fragment_cff.py \
        --python_filename $cfgfile \
        --eventcontent RAWSIM \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier GEN-LHE \
        --fileout file:$outfile \
        --conditions 106X_mc2017_realistic_v6 \
        --beamspot Realistic25ns13TeVEarly2017Collision \
        --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" \
        --step LHE,GEN \
        --geometry DB:Extended \
        --era Run2_2017 \
        --no_exec \
        --mc \
        --nThreads=$nThreads \
        -n $events || exit $? ;
    
    cmsRun $cfgfile
} &> $logfile 

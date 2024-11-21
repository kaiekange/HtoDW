#!/bin/bash

export HOME=/user/kakang

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

export X509_USER_PROXY=$1
date=$2
process=$3
events=$4
nThreads=$5

cd ../$date/samples/process_$process

stepname=DIGI
cfgfile=cfg/S2_${stepname}_cfg.py
infile=output/S1_SIM.root
outfile=output/S2_${stepname}.root
logfile=log/S2_${stepname}.log

{
    voms-proxy-info -all
    voms-proxy-info -all -file $X509_USER_PROXY

    cmsDriver.py \
        --python_filename $cfgfile \
        --eventcontent PREMIXRAW \
        --customise Configuration/DataProcessing/Utils.addMonitoring \
        --datatier GEN-SIM-DIGI \
        --filein file:$infile \
        --fileout file:$outfile \
        --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3/PREMIX" \
        --conditions 106X_mc2017_realistic_v6 \
        --beamspot Realistic25ns13TeVEarly2017Collision \
        --step DIGI,DATAMIX,L1,DIGI2RAW \
        --procModifiers premix_stage2 \
        --datamix PreMix \
        --geometry DB:Extended \
        --era Run2_2017 \
        --no_exec \
        --mc \
        --nThreads=$nThreads \
        -n $events || exit $? ;
    
    cmsRun $cfgfile
} &> $logfile 

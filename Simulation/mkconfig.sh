#!/bin/bash

set -e

ver=20260114
year=2022
outdir="/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="LHEGS"

cfgfile=${mode}_cfg.py

cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_ggHtoDsW_fragment_cff.py" \
    --eventcontent PREMIXRAW \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM-RAW \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step LHE,GEN,SIM,DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v12 \
    --procModifiers premix_stage2,siPixelQualityRawToDigi \
    --datamix PreMix \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX" \
    --geometry DB:Extended \
    --era Run3 \
    --fileout file:test.root \
    --python_filename "${cfgfile}" \
    --number 10 \
    --no_exec \
    --mc || exit $? ;

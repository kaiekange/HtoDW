#!/bin/bash

set -e

ver=20260114
year=2022
outdir="/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="NANO"

cfgfile=${mode}_cfg.py

# cmsDriver.py \
#     --eventcontent MINIAODSIM \
#     --customise Configuration/DataProcessing/Utils.addMonitoring \
#     --datatier MINIAODSIM \
#     --conditions 124X_mcRun3_2022_realistic_v12 \
#     --step PAT \
#     --geometry DB:Extended \
#     --era Run3 \
#     --filein file:test.root \
#     --fileout file:MiniAOD.root \
#     --python_filename "${cfgfile}" \
#     --number 10 \
#     --no_exec \
#     --mc || exit $? ;

cmsDriver.py \
    --eventcontent NANOEDMAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier NANOEDMAODSIM \
    --conditions 130X_mcRun3_2022_realistic_v5 \
    --step NANO \
    --era Run3,run3_miniAOD_12X \
    --scenario pp \
    --filein file:MiniAOD.root \
    --fileout file:NanoAOD.root \
    --python_filename "${cfgfile}" \
    --number 10 \
    --no_exec \
    --mc || exit $? ;

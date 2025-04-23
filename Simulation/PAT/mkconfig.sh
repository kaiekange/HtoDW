#!/bin/bash

set -e

export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="PAT"

mkdir -p "${pnfspath}/FullGEN/${mode}/log"

cmsDriver.py \
    --python_filename config_temp.py \
    --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier  MINIAODSIM \
    --filein "file:infile" \
    --fileout "file:outfile" \
    --conditions 106X_mc2017_realistic_v9 \
    --beamspot Realistic25ns13TeVEarly2017Collision \
    --step PAT \
    --procModifiers run2_miniAOD_UL \
    --geometry DB:Extended \
    --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' \
    --era Run2_2017 \
    --runUnscheduled \
    --no_exec \
    --mc \
    --nThreads=4 \
    -n -1 

sed -e "s%'file:infile'%sys.argv[2]%g" \
    -e "s%'file:outfile'%sys.argv[3]%g" \
    config_temp.py > config.py

sed -i "6iimport sys" config.py

rm config_temp.py

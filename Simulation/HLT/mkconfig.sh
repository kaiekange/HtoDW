#!/bin/bash

set -e

export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="HLT"

mkdir -p "${pnfspath}/FullGEN/${mode}/log"

cmsDriver.py \
    --python_filename config_temp.py \
    --eventcontent RAWSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM-RAW \
    --filein "file:infile" \
    --fileout "file:outfile" \
    --conditions 106X_mc2017_realistic_v9 \
    --beamspot Realistic25ns13TeVEarly2017Collision \
    --step HLT:Fake2 \
    --geometry DB:Extended \
    --era Run2_2017 \
    --no_exec \
    --mc \
    --nThreads=4 \
    -n -1 

sed -e "s%'file:infile'%sys.argv[2]%g" \
    -e "s%'file:outfile'%sys.argv[3]%g" \
    config_temp.py > config.py

sed -i "6iimport sys" config.py

rm config_temp.py

#!/bin/bash

set -e

export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="LHEGEN"

cfgfile=${mode}_cfg.py

cmsDriver.py Configuration/GenProduction/python/gridpacktogen_ggHtoDsW_fragment_cff.py \
    --python_filename ${cfgfile} \
    --eventcontent RAWSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-LHE \
    --fileout "file:output.root" \
    --conditions 106X_mc2017_realistic_v9 \
    --beamspot Realistic25ns13TeVEarly2017Collision \
    --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" \
    --step LHE,GEN \
    --geometry DB:Extended \
    --era Run2_2017 \
    --no_exec \
    --mc \
    -n 10

REQUESTNAME="MC_H2DW_${year}_${ver}"
UNITSPERJOB="200"
TOTALUNITS="100000"
OUTPUTDATASETTAG="H2DW"
OUTLFNDIRBASE="${pnfspath}"

sed -e "s%REQUESTNAME%${REQUESTNAME}%g" \
    -e "s%UNITSPERJOB%${UNITSPERJOB}%g" \
    -e "s%TOTALUNITS%${TOTALUNITS}%g" \
    -e "s%OUTPUTDATASETTAG%${OUTPUTDATASETTAG}%g" \
    -e "s%OUTLFNDIRBASE%${OUTLFNDIRBASE}%g" \
    "crabConfigTemplate.py" > "crabConfig.py"

source /cvmfs/cms.cern.ch/crab3/crab.sh

crab submit -c "crabConfig.py"

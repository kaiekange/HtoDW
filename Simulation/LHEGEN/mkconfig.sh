#!/bin/bash

set -e

# ver=20260114
# year=2022
# pnfspath="/store/user/kakang/Analysis/Simulation/${ver}/${year}"

SEED=$(( 10000000 + $1 * 1000 ))
echo $SEED

mode="SIM"

cfgfile=${mode}_cfg.py

cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_ggHtoDsW_fragment_cff.py" \
    --eventcontent RAWSIM,LHE \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM,LHE \
    --conditions 130X_mcRun3_2022_realistic_v5 \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step LHE,GEN,SIM \
    --geometry DB:Extended \
    --era Run3 \
    --fileout file:test.root \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(${SEED})" \
    --python_filename "${cfgfile}" \
    --number 10 \
    --no_exec \
    --mc || exit $? ;
    
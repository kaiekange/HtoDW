#!/bin/bash

mode="LHEGS"

ver=20260119
year=2022
outdir="/store/user/kakang/H2DW/Datasets/${ver}/Simulation"


cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_ggHtoDsW_fragment_cff.py" \
    --python_filename "${mode}_cfg.py" \
    --eventcontent RAWSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM \
    --fileout "output.root" \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step LHE,GEN,SIM \
    --geometry DB:Extended \
    --era Run3 \
    --number 10 \
    --no_exec \
    --mc || exit $? ;

REQUESTNAME="MC_ggHtoDsW_${year}_${ver}"
UNITSPERJOB="500"
TOTALUNITS="500000"
OUTPUTDATASETTAG="${year}"
OUTLFNDIRBASE=${outdir}

sed -e "s%REQUESTNAME%${REQUESTNAME}%g" \
    -e "s%UNITSPERJOB%${UNITSPERJOB}%g" \
    -e "s%TOTALUNITS%${TOTALUNITS}%g" \
    -e "s%OUTPUTDATASETTAG%${OUTPUTDATASETTAG}%g" \
    -e "s%OUTLFNDIRBASE%${OUTLFNDIRBASE}%g" \
    "crabConfigTemplate.py" > "crabConfig.py"

source /cvmfs/cms.cern.ch/crab3/crab.sh

# crab submit -c "crabConfig.py"

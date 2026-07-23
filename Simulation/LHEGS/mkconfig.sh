#!/bin/bash

ver=20260717

STEP="LHEGS"

YEAR="2022"
CONDITIONS="124X_mcRun3_2022_realistic_v12"
BEAMSPOT="Realistic25ns13p6TeVEarly2022Collision"

# YEAR="2022EE"
# CONDITIONS="124X_mcRun3_2022_realistic_postEE_v1"
# BEAMSPOT="Realistic25ns13p6TeVEarly2022Collision"

# MODE="ggHtoDsW_phipi"
# UNITSPERJOB="3200"
# TOTALUNITS="1600000"

# outdir="/store/user/kakang/Analysis/${MODE}/${ver}/Simulation"

# cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_${MODE}_fragment_cff.py" \
#     --python_filename "${STEP}_${MODE}_cfg.py" \
#     --eventcontent RAWSIM \
#     --customise Configuration/DataProcessing/Utils.addMonitoring \
#     --datatier GEN-SIM \
#     --fileout "output.root" \
#     --conditions ${CONDITIONS} \
#     --beamspot  ${BEAMSPOT}\
#     --step LHE,GEN,SIM \
#     --geometry DB:Extended \
#     --era Run3 \
#     --number 3200 \
#     --no_exec \
#     --mc || exit $? ;

MODE="ggHtoDsstW_phipi"
UNITSPERJOB="3000"
TOTALUNITS="3000000"

outdir="/store/user/kakang/Analysis/${MODE}/${ver}/Simulation"

cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_${MODE}_fragment_cff.py" \
    --python_filename "${STEP}_${MODE}_cfg.py" \
    --eventcontent RAWSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM \
    --fileout "output.root" \
    --conditions ${CONDITIONS} \
    --beamspot  ${BEAMSPOT}\
    --step LHE,GEN,SIM \
    --geometry DB:Extended \
    --era Run3 \
    --number 3000 \
    --no_exec \
    --mc || exit $? ;

REQUESTNAME="MC_${MODE}_${YEAR}_${ver}"
OUTPUTDATASETTAG="${YEAR}"
OUTLFNDIRBASE=${outdir}

sed -e "s%CONFIGFILE%${STEP}_${MODE}_cfg.py%g" \
    -e "s%REQUESTNAME%${REQUESTNAME}%g" \
    -e "s%UNITSPERJOB%${UNITSPERJOB}%g" \
    -e "s%TOTALUNITS%${TOTALUNITS}%g" \
    -e "s%OUTPUTDATASETTAG%${OUTPUTDATASETTAG}%g" \
    -e "s%OUTLFNDIRBASE%${OUTLFNDIRBASE}%g" \
    "crabConfigTemplate.py" > "crabConfig.py"

# source /cvmfs/cms.cern.ch/crab3/crab.sh

# crab submit -c "crabConfig.py"

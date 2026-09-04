#!/bin/bash

STEP="LHEGS"

MODE=$1
YEAR=$2

usage() {
    echo "usage: $0 <ggHtoDsW_phipi|ggHtoDsstW_phipi> <2022|2022EE>" >&2
    exit 1
}

# The GlobalTag is era-specific; the beamspot is the same for both 2022 eras.
BEAMSPOT="Realistic25ns13p6TeVEarly2022Collision"

# ver belongs in this table too: preEE was produced under 20260717, postEE under
# 20260825, and both appear in the output path. Kept outside, one era's tag would
# silently pair with the other era's GlobalTag.
case "${YEAR}" in
    2022)   CONDITIONS="124X_mcRun3_2022_realistic_v12"       ; ver=20260717 ;;
    2022EE) CONDITIONS="124X_mcRun3_2022_realistic_postEE_v1" ; ver=20260825 ;;
    *)      usage ;;
esac

# UNITSPERJOB doubles as cmsDriver's --number: the pset has to generate exactly
# the number of events CRAB accounts to one job, so the two cannot drift apart.
case "${MODE}" in
    ggHtoDsW_phipi)   UNITSPERJOB="3200" ; TOTALUNITS="1600000" ;;
    ggHtoDsstW_phipi) UNITSPERJOB="3000" ; TOTALUNITS="3000000" ;;
    *)                usage ;;
esac

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
    --number ${UNITSPERJOB} \
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

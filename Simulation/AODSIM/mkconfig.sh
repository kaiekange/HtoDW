#!/bin/bash

STEP="AODSIM"
MODE=$1
YEAR=$2

usage() {
    echo "usage: $0 <ggHtoDsW_phipi|ggHtoDsstW_phipi> <2022|2022EE>" >&2
    exit 1
}

[ -n "${MODE}" ] || usage

# The GlobalTag and the HLT menu are both era-specific and have to move together:
# a postEE GlobalTag paired with the preEE menu yields a sample labelled 2022EE
# whose trigger conditions are preEE, and nothing throws. Deriving both from YEAR
# makes that combination unrepresentable. The beamspot and the premix pileup
# dataset are the same for both 2022 eras, so they stay literal below.
case "${YEAR}" in
    2022)   CONDITIONS="124X_mcRun3_2022_realistic_v12"       ; HLTMENU="2022v12" ;;
    2022EE) CONDITIONS="124X_mcRun3_2022_realistic_postEE_v1" ; HLTMENU="2022v14" ;;
    *)      usage ;;
esac

cmsDriver.py \
    --python_filename "${STEP}_${MODE}_cfg.py" \
    --eventcontent AODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier AODSIM \
    --filein INPUTFILE \
    --fileout OUTPUTFILE \
    --conditions ${CONDITIONS} \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:${HLTMENU},RAW2DIGI,L1Reco,RECO,RECOSIM \
    --procModifiers premix_stage2,siPixelQualityRawToDigi \
    --datamix PreMix \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX" \
    --geometry DB:Extended \
    --era Run3 \
    --number -1 \
    --no_exec \
    --mc || exit $? ;

sed -i "6iimport sys" "${STEP}_${MODE}_cfg.py"
sed -i "s|'INPUTFILE'|sys.argv[2]|g" "${STEP}_${MODE}_cfg.py"
sed -i "s|'OUTPUTFILE'|sys.argv[3]|g" "${STEP}_${MODE}_cfg.py"

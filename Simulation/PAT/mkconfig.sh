#!/bin/bash

step="PAT"
mode=$1
YEAR=$2
cfgfile=${step}_${mode}_cfg.py

usage() {
    echo "usage: $0 <ggHtoDsW_phipi|ggHtoDsstW_phipi> <2022|2022EE>" >&2
    exit 1
}

[ -n "${mode}" ] || usage

# MiniAODv4 GlobalTag, era-specific. No default: an unset era must fail loudly
# rather than silently stamp preEE conditions onto a postEE sample.
case "${YEAR}" in
    2022)   CONDITIONS="130X_mcRun3_2022_realistic_v5" ;;
    2022EE) CONDITIONS="130X_mcRun3_2022_realistic_postEE_v6" ;;
    *)      usage ;;
esac

cmsDriver.py \
    --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier MINIAODSIM \
    --conditions ${CONDITIONS} \
    --step PAT \
    --geometry DB:Extended \
    --era Run3,run3_miniAOD_12X \
    --filein INPUTFILE \
    --fileout OUTPUTFILE \
    --python_filename "${cfgfile}" \
    --number -1 \
    --no_exec \
    --mc || exit $? ;

sed -i "6iimport sys" "${cfgfile}"
sed -i "s|'INPUTFILE'|sys.argv[2]|g" "${cfgfile}"
sed -i "s|'OUTPUTFILE'|sys.argv[3]|g" "${cfgfile}"

#!/bin/bash

step="PAT"
mode=$1
cfgfile=${step}_${mode}_cfg.py

cmsDriver.py \
    --eventcontent MINIAODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier MINIAODSIM \
    --conditions 130X_mcRun3_2022_realistic_v5 \
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

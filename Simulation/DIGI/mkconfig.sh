#!/bin/bash

set -e

export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="DIGI"

mkdir -p "${pnfspath}/FullGEN/${mode}/log"

voms-proxy-init --voms cms --hours 168

cmsDriver.py \
    --python_filename config.py \
    --eventcontent PREMIXRAW \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --datatier GEN-SIM-DIGI \
    --filein "file:infile" \
    --fileout "file:outfile" \
    --pileup_input "dbs:/Neutrino_E-10_gun/RunIIFall17FSPrePremix-PUFSUL17CP5_106X_mc2017_realistic_v9-v2/PREMIX" \
    --conditions 106X_mc2017_realistic_v9 \
    --beamspot Realistic25ns13TeVEarly2017Collision \
    --step DIGI,DATAMIX,L1,DIGI2RAW \
    --procModifiers premix_stage2 \
    --datamix PreMix \
    --geometry DB:Extended \
    --era Run2_2017 \
    --no_exec \
    --mc \
    --nThreads=4 \
    -n -1

linenumber=$(grep -nF "process.mixData.input.fileNames = cms.untracked.vstring([" config.py | cut -d: -f1)
head -n $(( $linenumber - 1 )) config.py > config_temp.py
pileup_files=$(echo -n "process.mixData.input.fileNames = cms.untracked.vstring([" >> config_temp.py; dasgoclient --query="file dataset=/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3/PREMIX site=T2_CH_CERN" | while IFS= read -r line; do echo -n "'$line',"; done)
pileup_files="${pileup_files%?}"
pileup_files="${pileup_files}])"
echo $pileup_files >> config_temp.py
tail -n +$(( $linenumber + 1 )) config.py >> config_temp.py

sed -e "s%'file:infile'%sys.argv[2]%g" \
    -e "s%'file:outfile'%sys.argv[3]%g" \
    config_temp.py > config.py

sed -i "6iimport sys" config.py

cp /tmp/x509up_u23273 /user/kakang/tmp/

rm config_temp.py

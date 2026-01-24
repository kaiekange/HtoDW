#!/bin/bash

cmsDriver.py "Configuration/GenProduction/python/gridpacktogen_ggHtoDsW_fragment_cff.py" \
    --python_filename "ggHtoDsW_cfg.py" \
    --eventcontent AODSIM \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=31415926 \
    --datatier AODSIM \
    --fileout OUTPUTFILE \
    --conditions 124X_mcRun3_2022_realistic_v12 \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --step LHE,GEN,SIM,DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v12,RAW2DIGI,L1Reco,RECO,RECOSIM \
    --procModifiers premix_stage2,siPixelQualityRawToDigi \
    --datamix PreMix \
    --pileup_input "dbs:/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX" \
    --geometry DB:Extended \
    --era Run3 \
    --number 500 \
    --no_exec \
    --mc || exit $? ;

sed -i "6iimport sys" ggHtoDsW_cfg.py
sed -i "s|'OUTPUTFILE'|sys.argv[2]|g" ggHtoDsW_cfg.py
sed -i "s|initialSeed=31415926|initialSeed=int(sys.argv[3])|g" ggHtoDsW_cfg.py

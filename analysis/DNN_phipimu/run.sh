#!/bin/bash

model_name=$1
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn

    SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN_phipimu
    BASE=/eos/home-k/kakang/H2DW/Analysis/MC/2022
    PRED=/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/${model_name}/Prediction

    mkdir -p ${PRED}

    # --- training ---
    python3 ${SCRIPTS}/DNN_phipimu.py ${model_name}

    # --- predict ---
    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/Signal_DsstW_H_phipimu/Signal_DsstW_H_phipimu_bestgamma.root \
        HphipimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_phipimu_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/Signal_H_phipimu/Signal_H_phipimu_bestgamma.root \
        HphipimuDsW/Events \
        ${PRED}/Signal_H_phipimu_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/WJets_H_phipimu/WJets_H_phipimu_bestgamma.root \
        Hphipimu/Events \
        ${PRED}/WJets_H_phipimu_scored.root

    # --- best candidate ---
    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_DsstW_H_phipimu_scored.root \
        HphipimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_phipimu_best.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_H_phipimu_scored.root \
        HphipimuDsW/Events \
        ${PRED}/Signal_H_phipimu_best.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/WJets_H_phipimu_scored.root \
        Hphipimu/Events \
        ${PRED}/WJets_H_phipimu_best.root

} &> /eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/${model_name}/training.log

#!/bin/bash

model_name=$1
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn

    SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN_Dsmu
    BASE=/eos/home-k/kakang/H2DW/Analysis/MC/2022
    PRED=/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/${model_name}/Prediction

    mkdir -p ${PRED}

    # --- training ---
    python3 ${SCRIPTS}/DNN_Dsmu.py ${model_name}

    # --- predict ---
    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/Signal_DsstW_H_KKpimu/Signal_DsstW_H_KKpimu_bestgamma.root \
        HKKpimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_KKpimu_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/Signal_H_KKpimu/Signal_H_KKpimu_bestgamma.root \
        HKKpimuDsW/Events \
        ${PRED}/Signal_H_KKpimu_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${BASE}/WJets_H_KKpimu/WJets_H_KKpimu_bestgamma.root \
        HKKpimu/Events \
        ${PRED}/WJets_H_KKpimu_scored.root

    # --- best candidate ---
    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_DsstW_H_KKpimu_scored.root \
        HKKpimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_KKpimu_best.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_H_KKpimu_scored.root \
        HKKpimuDsW/Events \
        ${PRED}/Signal_H_KKpimu_best.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/WJets_H_KKpimu_scored.root \
        HKKpimu/Events \
        ${PRED}/WJets_H_KKpimu_best.root

} &> /eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/${model_name}/training.log

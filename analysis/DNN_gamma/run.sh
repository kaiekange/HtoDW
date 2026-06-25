#!/bin/bash

model_name=$1
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn

    SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN_gamma
    PRED_DSMU=/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/v1/Prediction
    PRED=/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/${model_name}/Prediction

    mkdir -p ${PRED}

    # --- training ---
    python3 ${SCRIPTS}/DNN_gamma.py ${model_name}

    # --- predict (input: DNN_Dsmu _best files) ---
    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${PRED_DSMU}/Signal_DsstW_H_KKpimu_best.root \
        HKKpimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_KKpimu_best_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${PRED_DSMU}/Signal_H_KKpimu_best.root \
        HKKpimuDsW/Events \
        ${PRED}/Signal_H_KKpimu_best_scored.root

    python3 ${SCRIPTS}/predict.py ${model_name} \
        ${PRED_DSMU}/WJets_H_KKpimu_best.root \
        HKKpimu/Events \
        ${PRED}/WJets_H_KKpimu_best_scored.root

    # --- best gamma per event ---
    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_DsstW_H_KKpimu_best_scored.root \
        HKKpimuDsstW/Events \
        ${PRED}/Signal_DsstW_H_KKpimu_best_gamma.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/Signal_H_KKpimu_best_scored.root \
        HKKpimuDsW/Events \
        ${PRED}/Signal_H_KKpimu_best_gamma.root

    python3 ${SCRIPTS}/best_cand.py \
        ${PRED}/WJets_H_KKpimu_best_scored.root \
        HKKpimu/Events \
        ${PRED}/WJets_H_KKpimu_best_gamma.root

} &> /eos/home-k/kakang/H2DW/Analysis/MVA/gamma/${model_name}/training.log

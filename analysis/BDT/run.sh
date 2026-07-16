#!/bin/bash

model_name=$1
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn

    SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/BDT

    # --- training: signal vs W+jets, 26 features ---
    python3 ${SCRIPTS}/BDT_phipimu.py ${model_name}

    # --- score all 36 Samples (2 signals + 34 backgrounds) ---
    # No best_cand step: the Samples tuples are already one row per event.
    python3 ${SCRIPTS}/predict.py ${model_name}

} &> /eos/home-k/kakang/H2DW/Analysis/MVA/BDT/${model_name}/training.log

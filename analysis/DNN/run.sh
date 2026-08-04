#!/bin/bash

OUTDIR=/eos/home-k/kakang/H2DW/Analysis/MVA/DNN/dsw_dsstw_wjets
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn

    SCRIPTS=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN

    # --- training: HtoDsW+HtoDsstW vs lumi-reweighted WtoLNu2Jets (0J+1J+2J), 40 features ---
    # NOT comparable to BDT/BDT.py any more: that one is still on the 42-feature set
    # with a single signal channel.
    # Takes no arguments: the old v1/v2/v3 background-mixing configs are gone along
    # with the four-background Samples/ files they referred to.
    python3 ${SCRIPTS}/DNN.py

    # --- score the current Samples/ files (signal + 3 jet bins + data, OS and SS) ---
    # No best_cand step: the Samples tuples are already one row per event.
    python3 ${SCRIPTS}/predict.py

} &> ${OUTDIR}/training.log

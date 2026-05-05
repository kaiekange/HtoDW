#!/bin/bash

model_name=$1
{
    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/DNN_ver2.py ${model_name}
} &> /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsW_KKpimu/DNN/${model_name}/training.log

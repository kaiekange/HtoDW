#!/bin/bash
# One-off condor wrapper for check_data_phi20nano_vs_old.py (delete after use).
set -euo pipefail
source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
conda activate ds_dnn
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/macros
python3 -u check_data_phi20nano_vs_old.py "$1" "$2"

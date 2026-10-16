#!/bin/bash
#   run_ds_mass_fit.sh --key=ssdata_phi20 --cuts=0.70,0.60,0.50 --joint=slices|nested --float=r34,r5,k [--tag=T]
#
# Ds-mass sideband fit, all three steps, so no fit is ever left without its figures:
#   1. ds_mass_fit_io.py inputs   (conda ds_dnn: needs common.py -> torch)
#   2. fit_ds_mass_roofit.py      (cmsenv: needs ROOT; ds_dnn has none)
#   3. ds_mass_fit_io.py draw     (conda ds_dnn: mplhep 1.1.2, the version the label rules were checked on)
# Arguments are passed unchanged to steps 1 and 3; step 2 gets the fit dir printed by step 1.

set -euo pipefail
HERE=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/macros
PY=/afs/cern.ch/work/k/kakang/miniconda3/envs/ds_dnn/bin/python3
cd "${HERE}"   # never /tmp: a stray /tmp/inspect.py shadows the stdlib there

LOG=$(mktemp)
"${PY}" -u ds_mass_fit_io.py inputs "$@" 2>&1 | tee "${LOG}"
FITDIR=$(sed -n 's|^Wrote \(.*\)/inputs.npz, meta.json$|\1|p' "${LOG}")
rm -f "${LOG}"
[ -n "${FITDIR}" ] || { echo "run_ds_mass_fit.sh: inputs step did not report a fit dir" >&2; exit 1; }

(
  set +u   # cmsset_default.sh / scramv1 reference unset variables
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
  eval "$(scramv1 runtime -sh)"
  cd "${HERE}"
  python3 -u fit_ds_mass_roofit.py "${FITDIR}"
)

"${PY}" -u ds_mass_fit_io.py draw "$@"

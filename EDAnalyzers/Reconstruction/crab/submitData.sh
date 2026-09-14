#!/bin/bash
#
# Submit the 2022 Muon data listed in ../../../../CMSSW_15_0_5/src/Tuples/Data/DatasetList/data_2022.txt
# (shared with the old H_phipimu-tree workflow -- it's a pure dataset/GlobalTag/Year table,
# not release-specific) through the new custom-NanoAOD pset (../test/NANO_DATA_H_phipimu.py).
#
#   ./submitData.sh OS                    # opposite-sign (SS=0), every active line
#   ./submitData.sh SS                    # same-sign-pi control region (SS=1)
#   ./submitData.sh OS --dryrun           # crab submit --dryrun
#   ./submitData.sh OS --print            # write generated/<requestName>.py, submit nothing
#   ./submitData.sh OS submit Run2022C    # 4th arg: only datasets whose era contains this
#   ./submitData.sh OS --dryrun Run2022C  # e.g. a single-era smoke test before the full table
#
# Comment a line out in the txt (leading '#') to skip that dataset -- that is how you
# choose what to submit, rather than editing this script. See that file's header for the
# GlobalTag/Year column semantics (two independent axes, not one preEE/postEE flag).
#
# Adapted from Tuples/Data/submitData.sh -- same mechanics (sed-substituted template,
# .requestcache success check, mandatory cmsenv before crab.sh), pointed at the nano pset
# and CMSSW_13_0_13.

set -euo pipefail

CHANNEL="${1:-}"
MODE="${2:-submit}"
ONLY="${3:-}"   # optional: substring match against the era (e.g. "Run2022C")

case "${CHANNEL}" in
    OS) SS_VALUE="0"; SFX="" ;;
    SS) SS_VALUE="1"; SFX="_SS" ;;
    *)  echo "usage: $0 <OS|SS> [--dryrun|--print]" >&2; exit 1 ;;
esac

case "${MODE}" in
    submit|--dryrun|--print) ;;
    *) echo "usage: $0 <OS|SS> [--dryrun|--print]" >&2; exit 1 ;;
esac

cd "$(dirname "$0")"

datasetlist="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Data/DatasetList/data_2022.txt"
configtemplate="crabConfigTemplateData.py"
PSETNAME="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/NANO_DATA_H_phipimu.py"

[[ -f "${datasetlist}"     ]] || { echo "missing ${datasetlist}" >&2;     exit 1; }
[[ -f "${configtemplate}"  ]] || { echo "missing ${configtemplate}" >&2;  exit 1; }
[[ -f "${PSETNAME}"        ]] || { echo "missing pset ${PSETNAME}" >&2;   exit 1; }

CMSSW_SRC="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src"

if [[ "${MODE}" != "--print" ]]; then
    # cmsenv BEFORE crab.sh, and unconditionally rather than relying on the caller's shell:
    # `crab submit` aborts with "CMSSW is missing. You must do cmsenv first" otherwise, and
    # it does so while still returning 0, so set -e does not catch it.
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    pushd "${CMSSW_SRC}" >/dev/null
    eval "$(scramv1 runtime -sh)"
    popd >/dev/null
    [[ -n "${CMSSW_BASE:-}" ]] || { echo "cmsenv failed: CMSSW_BASE unset" >&2; exit 1; }

    # 'prod' must be passed explicitly: `source` hands the *current* script's positional
    # parameters to the sourced one, so a bare `source crab.sh` makes crab.sh read this
    # script's $1 ("OS"/"SS") as its own type argument and fail with
    #   ERROR: Invalid CRAB type 'OS' provided, valid values are prod, pre and dev.
    source /cvmfs/cms.cern.ch/crab3/crab.sh prod
    command -v crab >/dev/null 2>&1 || {
        echo "crab is not on PATH after sourcing crab.sh -- refusing to continue" >&2
        exit 1
    }
fi

# Read the whole file before submitting anything, so a malformed line fails the run
# before half the datasets are already queued.
datasets=(); globaltags=(); years=()
lineno=0
while IFS= read -r line || [[ -n "${line}" ]]; do
    lineno=$((lineno + 1))
    line="${line%%#*}"                      # strip comments, including trailing ones
    [[ -z "${line// }" ]] && continue       # skip blank / comment-only lines
    read -r ds gt yr extra <<< "${line}"
    if [[ -z "${ds}" || -z "${gt}" || -z "${yr}" || -n "${extra:-}" ]]; then
        echo "${datasetlist}:${lineno}: expected exactly 3 columns <dataset> <GlobalTag> <Year>, got: ${line}" >&2
        exit 1
    fi
    if [[ -n "${ONLY}" && "${ds}" != *"${ONLY}"* ]]; then
        continue
    fi
    datasets+=("${ds}"); globaltags+=("${gt}"); years+=("${yr}")
done < "${datasetlist}"

if [[ ${#datasets[@]} -eq 0 ]]; then
    if [[ -n "${ONLY}" ]]; then
        echo "no active line in ${datasetlist} matches '${ONLY}'" >&2
    else
        echo "no active lines in ${datasetlist}" >&2
    fi
    exit 1
fi
echo "[submitData] ${CHANNEL}, ${#datasets[@]} dataset(s)${ONLY:+ matching '${ONLY}'}, mode=${MODE}"

for i in "${!datasets[@]}"; do

    dataset="${datasets[i]}"
    GLOBALTAG="${globaltags[i]}"
    YEAR="${years[i]}"

    # /Muon/Run2022C-22Sep2023-v1/MINIAOD -> primary=Muon, era=Run2022C
    split=(${dataset//\// })
    primary="${split[0]}"
    era="${split[1]%%-*}"

    # Output directory era, from the Year column. Same 2022 / 2022EE names the MC tree uses.
    case "${YEAR}" in
        2022Pre)  ERADIR="2022"   ;;
        2022Post) ERADIR="2022EE" ;;
        *) echo "no output directory mapped for Year='${YEAR}'" >&2; exit 1 ;;
    esac

    # NanoH prefix (not H_phipimu) to keep crab_*/ dirs and DBS output-dataset tags from
    # colliding with the old TFileService-tree production's crab_<Muon_Run2022C_H_phipimu>/.
    REQUESTNAME="${primary}_${era}_NanoH_phipimu${SFX}"
    OUTPUTDATASETTAG="${REQUESTNAME}"

    rm -f crabConfig.py*

    # '|' as the sed delimiter throughout: dataset names and pset paths contain '/'.
    sed -e "s|__REQUESTNAME__|${REQUESTNAME}|g" \
        -e "s|__PSETNAME__|${PSETNAME}|g" \
        -e "s|__GLOBALTAG_VALUE__|${GLOBALTAG}|g" \
        -e "s|__SS_VALUE__|${SS_VALUE}|g" \
        -e "s|__INPUTDATASET__|${dataset}|g" \
        -e "s|__ERADIR__|${ERADIR}|g" \
        -e "s|__OUTPUTDATASETTAG__|${OUTPUTDATASETTAG}|g" \
        "${configtemplate}" > crabConfig.py

    if grep -q '__[A-Z_]*__' crabConfig.py; then
        echo "unsubstituted placeholder left in crabConfig.py for ${dataset}:" >&2
        grep -n '__[A-Z_]*__' crabConfig.py >&2
        exit 1
    fi

    echo "[submitData] ${REQUESTNAME}  GT=${GLOBALTAG}  SS=${SS_VALUE}  -> NanoMVA/Data/${ERADIR}"

    case "${MODE}" in
        submit)    crab submit -c crabConfig.py ;;
        --dryrun)  crab submit -c crabConfig.py --dryrun ;;
        --print)   mkdir -p generated
                   cp crabConfig.py "generated/${REQUESTNAME}.py"
                   echo "           -> generated/${REQUESTNAME}.py" ;;
    esac

    # `crab` is a shell function from crab.sh and returns 0 even when it aborts -- the
    # "CMSSW is missing. You must do cmsenv first" abort exits 0 and set -e sails past it.
    # Check an artefact instead of the exit status.
    #
    # The artefact is crab_<name>/.requestcache, NOT the directory: crab creates
    # crab_<name>/ with a crab.log before it attempts proxy delegation, so a submission
    # that dies at "Problems delegating My-proxy" still leaves the directory behind.
    # .requestcache is written only once the task is actually registered.
    if [[ "${MODE}" == "submit" && ! -f "crab_${REQUESTNAME}/.requestcache" ]]; then
        echo "no crab_${REQUESTNAME}/.requestcache -- submission FAILED, stopping." >&2
        echo "see crab_${REQUESTNAME}/crab.log" >&2
        exit 1
    fi

done

#!/bin/bash
#
# Submit the 2022 Muon data listed in DatasetList/data_2022.txt.
#
#   ./submitData.sh OS              # opposite-sign, DATA_H_phipimu.py
#   ./submitData.sh SS              # same-sign control region, DATA_H_phipimu_SS.py
#   ./submitData.sh OS --dryrun     # crab submit --dryrun
#   ./submitData.sh OS --print      # write generated/<requestName>.py, submit nothing
#
# Comment a line out in the txt (leading '#') to skip that dataset -- that is how you
# choose what to submit, rather than editing this script.
#
# DatasetList/*.txt is NOT under version control (.gitignore takes *.txt), so its format
# is recorded here in case it has to be rebuilt. Whitespace-separated, '#' comments:
#
#   <dataset>  <GlobalTag>  <Year>
#
#   GlobalTag: read off DAS for the dataset itself, never inferred --
#     dasgoclient -query="config dataset=<ds> system=dbs3 | grep config.global_tag"
#     For 2022 Muon 22Sep2023: C/D/E -> 130X_dataRun3_v2, F/G ->
#     130X_dataRun3_PromptAnalysis_v1 (F and G were not re-reconstructed in that
#     campaign). It reaches the kinematic vertex fit via TransientTrackBuilder.
#   Year: 2022Pre for C/D, 2022Post for E/F/G. Selects both the JEC payload set and the
#     jet veto map, on three separate modules inside the pset.
#
#   So the two columns are independent: Run2022E is postEE in Year but shares C/D's
#   GlobalTag, and F/G differ from E in GlobalTag but not in Year.
#
# Deviation from crab/submitData.sh in the IPres package, on purpose: OS/SS is a required
# argument instead of an uncomment-one-of-several block, because there are 5 datasets per
# sign and an uncommenting mistake submits all 10 at once. Note also that the GlobalTag
# lives per-line in the txt rather than in a per-YEAR block up here -- it has to, since
# Run2022E and Run2022F share a Year but not a GlobalTag.

set -euo pipefail

CHANNEL="${1:-}"
MODE="${2:-submit}"

case "${CHANNEL}" in
    OS) PSETNAME="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/DATA_H_phipimu.py";    SFX="" ;;
    SS) PSETNAME="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/DATA_H_phipimu_SS.py"; SFX="_SS" ;;
    *)  echo "usage: $0 <OS|SS> [--dryrun|--print]" >&2; exit 1 ;;
esac

case "${MODE}" in
    submit|--dryrun|--print) ;;
    *) echo "usage: $0 <OS|SS> [--dryrun|--print]" >&2; exit 1 ;;
esac

cd "$(dirname "$0")"

datasetlist="DatasetList/data_2022.txt"
configtemplate="crabConfigTemplateData.py"

[[ -f "${datasetlist}"     ]] || { echo "missing ${datasetlist}" >&2;     exit 1; }
[[ -f "${configtemplate}"  ]] || { echo "missing ${configtemplate}" >&2;  exit 1; }
[[ -f "${PSETNAME}"        ]] || { echo "missing pset ${PSETNAME}" >&2;   exit 1; }

CMSSW_SRC="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src"

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
    datasets+=("${ds}"); globaltags+=("${gt}"); years+=("${yr}")
done < "${datasetlist}"

[[ ${#datasets[@]} -gt 0 ]] && echo "[submitData] ${CHANNEL}, ${#datasets[@]} dataset(s), mode=${MODE}" \
                            || { echo "no active lines in ${datasetlist}" >&2; exit 1; }

for i in "${!datasets[@]}"; do

    dataset="${datasets[i]}"
    GLOBALTAG="${globaltags[i]}"
    YEAR="${years[i]}"

    # /Muon/Run2022C-22Sep2023-v1/MINIAOD -> primary=Muon, era=Run2022C
    split=(${dataset//\// })
    primary="${split[0]}"
    era="${split[1]%%-*}"

    # Output directory era, from the Year column. Same 2022 / 2022EE names the MC tree uses,
    # so Data/ and MC/ stay searchable with one convention.
    case "${YEAR}" in
        2022Pre)  ERADIR="2022"   ;;
        2022Post) ERADIR="2022EE" ;;
        *) echo "no output directory mapped for Year='${YEAR}'" >&2; exit 1 ;;
    esac

    # Reproduces the names already used for C/D (Muon_Run2022C_H_phipimu). Do not
    # "simplify" this to include the processing version: outLFNDirBase/<primary>/<tag>/
    # is where the 2026-08-03 output lives, and a changed tag scatters new output into
    # new directories and breaks the merge step without any error.
    REQUESTNAME="${primary}_${era}_H_phipimu${SFX}"
    OUTPUTDATASETTAG="${REQUESTNAME}"

    rm -f crabConfig.py*

    # '|' as the sed delimiter throughout: dataset names and pset paths contain '/'.
    sed -e "s|__REQUESTNAME__|${REQUESTNAME}|g" \
        -e "s|__PSETNAME__|${PSETNAME}|g" \
        -e "s|__GLOBALTAG_VALUE__|${GLOBALTAG}|g" \
        -e "s|__YEAR_VALUE__|${YEAR}|g" \
        -e "s|__INPUTDATASET__|${dataset}|g" \
        -e "s|__ERADIR__|${ERADIR}|g" \
        -e "s|__OUTPUTDATASETTAG__|${OUTPUTDATASETTAG}|g" \
        "${configtemplate}" > crabConfig.py

    if grep -q '__[A-Z_]*__' crabConfig.py; then
        echo "unsubstituted placeholder left in crabConfig.py for ${dataset}:" >&2
        grep -n '__[A-Z_]*__' crabConfig.py >&2
        exit 1
    fi

    echo "[submitData] ${REQUESTNAME}  GT=${GLOBALTAG}  Year=${YEAR}  -> Data/${ERADIR}"

    case "${MODE}" in
        submit)    crab submit -c crabConfig.py ;;
        --dryrun)  crab submit -c crabConfig.py --dryrun ;;
        --print)   mkdir -p generated
                   cp crabConfig.py "generated/${REQUESTNAME}.py"
                   echo "           -> generated/${REQUESTNAME}.py" ;;
    esac

    # `crab` is a shell function from crab.sh and returns 0 even when it aborts -- the
    # "CMSSW is missing. You must do cmsenv first" abort exits 0 and set -e sails past it,
    # so three datasets were once "submitted" with nothing created. Check an artefact
    # instead of the exit status, and stop rather than repeating the failure per dataset.
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

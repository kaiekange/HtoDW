#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="PVStudy_noBS"

mkdir -p "${pnfspath}/tuples/${mode}/log"

infile="${pnfspath}/FullGEN/PAT/output_${process}.root"
outfile="${pnfspath}/tuples/${mode}/output_${process}.root"
logfile="${pnfspath}/tuples/${mode}/log/output_${process}.log"
cfgfile="EDAnalyzers/RecoAnalyzer/test/${mode}_cfg.py"

if [ ! -f "${logfile}" ] || [ ! -f "${outfile}" ] || ! grep -q "Begin processing the 200th record" "${logfile}" || ! grep -q "dropped waiting message count 0" "${logfile}"; then

    echo "Tuple ${process} not produced, now producing tuple" >> job_noBS.out
    rm -f "${logfile}" "${outfile}"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cmssw-el7 -- "cd ${HOME}Analysis/CMSSW_10_6_30_patch1/src/; cmsenv; cmsRun ${cfgfile} file:${infile} file:${outfile}" &> ${logfile}

    if grep -q "Begin processing the 200th record" "${logfile}" && grep -q "dropped waiting message count 0" ${logfile}; then
        echo "Tuple ${process} successfully produced" >> job_noBS.out
    else
        echo "ERROR: Tuple ${process} production failed" >> job_noBS.out
    fi

else
    echo "Tuple ${process} already successfully produced" >> job_noBS.out
fi

#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode=SelectionStudy

outfile="${pnfspath}/tuples/${mode}/checkpt_${process}.root"
logfile="${pnfspath}/tuples/${mode}/log/checkpt_${process}.log"

if [ ! -f "${logfile}" ] || [ ! -f "${outfile}" ] || ! grep -q "(int) 0" "${logfile}"; then

    echo "Tuple ${process} not checkpt, now producing tuple" >> job_checkpt.out
    rm -f "${logfile}" "${outfile}"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cmssw-el7 -- "cd ${HOME}Analysis/CMSSW_10_6_30_patch1/src/; cmsenv; cd scripts/nonmatch/; root -l -b -q checkpt.cc\(${process}\)" &> ${logfile}

    if grep -q "(int) 0" "${logfile}"; then
        echo "Tuple ${process} successfully checkpt" >> job_checkpt.out
    else
        echo "ERROR: Tuple ${process} checkpt failed" >> job_checkpt.out
    fi

else
    echo "Tuple ${process} already successfully checkpt" >> job_checkpt.out
fi

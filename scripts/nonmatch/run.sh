#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode=SelectionStudy

outfile="${pnfspath}/tuples/${mode}/skimmed_${process}.root"
logfile="${pnfspath}/tuples/${mode}/log/skimmed_${process}.log"

if [ ! -f "${logfile}" ] || [ ! -f "${outfile}" ] || ! grep -q "(int) 0" "${logfile}"; then

    echo "Tuple ${process} not skimmed, now producing tuple" >> job_skimmed.out
    rm -f "${logfile}" "${outfile}"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cmssw-el7 -- "cd ${HOME}Analysis/CMSSW_10_6_30_patch1/src/; cmsenv; cd scripts/nonmatch/; root -l -b -q skimmer.cc\(${process}\)" &> ${logfile}

    if grep -q "(int) 0" "${logfile}"; then
        echo "Tuple ${process} successfully skimmed" >> job_skimmed.out
    else
        echo "ERROR: Tuple ${process} skimming failed" >> job_skimmed.out
    fi

else
    echo "Tuple ${process} already successfully skimmed" >> job_skimmed.out
fi

#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/${ver}/${year}"

mode="RECO"
premode="HLT"

infile="${pnfspath}/FullGEN/${premode}/output_${process}.root"
outfile="${pnfspath}/FullGEN/${mode}/output_${process}.root"
logfile="${pnfspath}/FullGEN/${mode}/log/output_${process}.log"

if [ ! -f $logfile ] || [ ! -f $outfile ] || ! grep -q "Begin processing the 200th record" $logfile || ! grep -q "dropped waiting message count" $logfile; then

    echo "Tuple $process not produced, now producing tuple" >> job.out
    rm -f "$logfile" "$outfile"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cmssw-el7 -- "cd ${HOME}Analysis/CMSSW_10_6_30_patch1/src/; cmsenv; cd Simulation/${mode}; cmsRun config.py file:${infile} file:${outfile}" &> ${logfile}

    if grep -q "Begin processing the 200th record" $logfile && grep -q "dropped waiting message count" $logfile; then
        echo "Tuple $process successfully produced" >> job.out
    else
        echo "ERROR: Tuple $process production failed" >> job.out
    fi

else
    echo "Tuple $process already successfully produced" >> job.out
fi

#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/

ver=20250202
year=2017
pnfspath="/pnfs/iihe/cms/store/user/kakang/Higgs_charm/Simulation/${ver}/${year}"

mode=GenMatch

mkdir -p "${pnfspath}/tuples/${mode}/log"

infile="${pnfspath}/FullGEN/PAT/output_${process}.root"
outfile="${pnfspath}/tuples/${mode}/output_${process}.root"
logfile="${pnfspath}/tuples/${mode}/log/output_${process}.log"
cfgfile="EDAnalyzers/GenParticleAnalyzer/python/${mode}_cfg.py"
# cfgfile="${HOME}Higgs_charm/CMSSW_10_6_30_patch1/src/EDAnalyzers/GenParticleAnalyzer/python/${mode}_cfg.py"

if [ ! -f "${logfile}" ] || [ ! -f "${outfile}" ] || ! grep -q "dropped waiting message count 0" "${logfile}"; then

    echo "Tuple ${process} not produced, now producing tuple" >> GenMatch.out
    rm -f "${logfile}" "${outfile}"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    
    
    cmssw-el7 -- "cd ${HOME}Higgs_charm/CMSSW_10_6_30_patch1/src/; cmsenv; cmsRun ${cfgfile} file:${infile} file:${outfile}" &> ${logfile}
    # cd H2DW/scripts

    if grep -q "dropped waiting message count 0" ${logfile}; then
        echo "Tuple ${process} successfully produced" >> GenMatch.out
    else
        echo "ERROR: Tuple ${process} production failed" >> GenMatch.out
    fi

else
    echo "Tuple ${process} already successfully produced" >> GenMatch.out
fi

#!/bin/bash

process=$(( $1 + 1 ))
export HOME=/user/kakang/
export X509_USER_PROXY=/user/kakang/tmp/x509up_u23273

ver=20250417
year=2017UL
pnfspath="/pnfs/iihe/cms/store/user/kakang/Analysis/Reconstruction/${ver}/${year}/WJetsBKG"

mkdir -p "${pnfspath}/tuples/${mode}/log"

# infile="${pnfspath}/FullGEN/PAT/output_${process}.root"
infile=$(sed -n "${process}p" filelist.txt)
outfile="${pnfspath}/tuples/${mode}/output_${process}.root"
logfile="${pnfspath}/tuples/${mode}/log/output_${process}.log"
cfgfile="EDAnalyzers/RecoAnalyzer/test/RecoBestAnalyzer_cfg2.py"

if [ ! -f "${logfile}" ] || [ ! -f "${outfile}" ] || ! grep -q "Begin processing the 2000th record" "${logfile}" || ! grep -q "dropped waiting message count 0" "${logfile}"; then

    echo "Tuple ${process} not produced, now producing tuple" >> job.out
    rm -f "${logfile}" "${outfile}"

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    voms-proxy-info -all
    voms-proxy-info -all -file $X509_USER_PROXY

    cmssw-el7 -- "cd ${HOME}Analysis/CMSSW_10_6_30_patch1/src/; cmsenv; cmsRun ${cfgfile} ${infile} file:${outfile}" &>> ${logfile}

    if grep -q "Begin processing the 2000th record" "${logfile}" && grep -q "dropped waiting message count 0" ${logfile}; then
        echo "Tuple ${process} successfully produced" >> job.out
    else
        echo "ERROR: Tuple ${process} production failed" >> job.out
    fi

else
    echo "Tuple ${process} already successfully produced" >> job.out
fi

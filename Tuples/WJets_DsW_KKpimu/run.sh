#!/bin/bash

GroupSize=$8
Process=$9

StartLine=$(( Process * GroupSize + 1 ))
EndLine=$(( StartLine + GroupSize - 1 ))

CfgFile="/afs/cern.ch/work/k/kakang/IPres/CMSSW_15_0_16/src/TrackingAnalysis/EDAnalyzers/test/residuals_condor_group.py"
LogFile="${OutputDir}/log/output_$(( Process + 1 )).log"
InputFile="$(sed -n "${StartLine}, ${EndLine}p" "${InputList}")"
OutputFile="${OutputDir}/output_$(( Process + 1 )).root"

if [ -e "${OutputFile}" ]; then
    echo "Removing existing ${OutputFile}"
    rm -f "${OutputFile}"
fi

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/IPres/CMSSW_15_0_16/src/
cmsenv
voms-proxy-info -all -file $X509_USER_PROXY
cmsRun ${CfgFile} RunOnData=${RunOnData} GlobalTag=${GlobalTag} EventScale=${EventScale} EventModulo=${EventModulo} SampleType="${SampleType}" InputFile="${InputFile}" OutputFile="${OutputFile}" GoldenJSON="${GoldenJSON}" &> ${LogFile}
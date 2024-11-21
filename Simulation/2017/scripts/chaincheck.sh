#!/bin/bash

export HOME=/user/kakang

source $VO_CMS_SW_DIR/cmsset_default.sh
cmsenv

process=$1

cfgfile=/user/kakang/H2cc/CMSSW_10_6_30_patch1/src/EDAnalyzer/GenParticleAnalyzer/python/GenPart_cfg.py
logfile=/user/kakang/H2cc/CMSSW_10_6_30_patch1/src/Simulation/2017/20241113/analysis/chaincheck/log/tuple_${process}.log
{
    cmsRun $cfgfile $process
} &> $logfile 

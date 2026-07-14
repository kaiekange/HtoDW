#!/bin/bash

export X509_USER_PROXY="/afs/cern.ch/user/k/kakang/tmp/x509up_u158898"

SAMPLE=$1
WEIGHT=$2
INFILE=$3
SUBDIR=$4

OUTDIR="/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
MACRO="/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/Samples/best_cand_tuple.cc"

OUTFILE="${OUTDIR}/${SAMPLE}.root"
LOGFILE="${OUTDIR}/log/${SAMPLE}.log"
rm -f $OUTFILE $LOGFILE

root -l -b -q "${MACRO}(\"${INFILE}\", \"${OUTFILE}\", \"${SUBDIR}\", ${WEIGHT}, \"Hphipimu\")" &> ${LOGFILE}

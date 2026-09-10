"""WtoLNu2Jets tuplizer cfg, parametrized over jet bin and OS/SS.

Replaces the near-identical cfgs this used to be (one per jet bin x OS/SS), which
differed only in the file list, the xrootd redirector and the process name.

    cmsRun WtoLNu2Jets_H_phipimu.py ERA=2022 JETBIN=0J SS=1 JOBID=0 GROUPSIZE=20 OUTFILE=...

ERA is 2022 (Run3Summer22, preEE) or 2022EE (Run3Summer22EE, postEE); JETBIN is one of
0J/1J/2J; SS=1 selects the same-sign control region. Process name and module label are
both `Hphipimu` in every combination, so the tree path is always Hphipimu/Events.

The filelist*.txt dumps are gitignored (*.txt) and must exist locally; they are
static dasgoclient dumps, the datasets being published and fixed (unlike the
still-running private signal productions, which need a live existence check).
"""
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra

options = VarParsing('analysis')
options.register("JOBID", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Job ID")
options.register("GROUPSIZE", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Group size")
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.register('JETBIN', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, '0J | 1J | 2J')
# int, not bool: VarParsing's bool type reads SS=0 as truthy, which would silently
# write OS data into an SS output directory.
options.register('SS', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'same-sign control region (0/1)')
# No default: the two eras need different GlobalTags AND different JEC/veto-map eras, and
# getting that wrong produces plausible numbers rather than an error. Required, exactly
# like GLOBALTAG/YEAR in the data cfgs.
options.register('ERA', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 '2022 (Run3Summer22, preEE) | 2022EE (Run3Summer22EE, postEE)')
# Which analyzer to run. Parametrized rather than given its own cfg on purpose: the era ->
# (GlobalTag, JEC era) table below, the redirector and the file-list slicing must stay a
# single source of truth. A second copy of this cfg that drifted on the era table would
# produce plausible numbers with the wrong JEC rather than an error.
options.register('TUPLIZER', 'H_phipimu', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'H_phipimu (default) | H_phipimu_GenTag (background KKpi truth tagging)')
options.parseArguments()

TUPLIZERS = ("H_phipimu", "H_phipimu_GenTag")
if options.TUPLIZER not in TUPLIZERS:
    sys.exit(f"TUPLIZER='{options.TUPLIZER}' is not one of {sorted(TUPLIZERS)}")

# Era -> (GlobalTag, JME era key). Two things that must move together:
#   2022    Run3Summer22MiniAODv4    130X_mcRun3_2022_realistic_v5          JEC 2022Pre
#   2022EE  Run3Summer22EEMiniAODv4  130X_mcRun3_2022_realistic_postEE_v6   JEC 2022Post
# The JME half is guarded: the cfi Year default is empty and all three JME modules throw on
# it, so a 2022EE run that forgot to switch cannot silently use preEE JEC or the preEE jet
# veto map. The GlobalTag half has no such guard -- hence no default here either.
ERAS = {
    "2022":   ("130X_mcRun3_2022_realistic_v5",        "2022Pre"),
    "2022EE": ("130X_mcRun3_2022_realistic_postEE_v6", "2022Post"),
}
if options.ERA not in ERAS:
    sys.exit(f"ERA='{options.ERA}' is not one of {sorted(ERAS)}")
GLOBALTAG, JEC_YEAR = ERAS[options.ERA]

# Per (era, jet bin): the file-list suffix and the redirector that can serve it. The
# jet-binned datasets are not replicated at T2_BE_IIHE, so they go through the global AAA
# redirector, which finds whichever site actually holds each file. File lists are named
# filelist_<bin>_<era>.txt throughout. The inclusive sample was dropped on 2026-08-18:
# nothing read it (common.py merges 0J+1J+2J), and it duplicated the same physical process.
AAA = "root://xrootd-cms.infn.it/"
DATASETS = {
    ("2022",   "0J"):   ("_0J_2022",   AAA),
    ("2022",   "1J"):   ("_1J_2022",   AAA),
    ("2022",   "2J"):   ("_2J_2022",   AAA),
    ("2022EE", "0J"):   ("_0J_2022EE", AAA),
    ("2022EE", "1J"):   ("_1J_2022EE", AAA),
    ("2022EE", "2J"):   ("_2J_2022EE", AAA),
}
key = (options.ERA, options.JETBIN)
if key not in DATASETS:
    sys.exit(f"(ERA, JETBIN) = {key} is not one of {sorted(DATASETS)}")
suffix, redirector = DATASETS[key]
sameSign = bool(options.SS)

start = options.JOBID * options.GROUPSIZE
end = start + options.GROUPSIZE

FILELIST = f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WtoLNu2Jets/filelist{suffix}.txt"
with open(FILELIST) as fh:
    ALL_FILES = [line.strip() for line in fh if line.strip()]
INFILE = [redirector + f for f in ALL_FILES[start:end]]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

# Nothing downstream reads the process name (the tree path follows the module label), and a
# composed one must stay alphanumeric. Fixed value; bin/era/sign live in the output path.
process = cms.Process("Hphipimu")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, GLOBALTAG, "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")

# Sets Year on all three JME modules from the era chosen above. The cfi default is empty
# and every JME plugin throws without it, so a 2022EE run cannot silently inherit the preEE
# JEC payloads or the preEE jet veto map.
applyJmeEra(process, year=JEC_YEAR, isData=False)
# Module label is Hphipimu either way, so the tree path stays Hphipimu/Events and the
# downstream tooling (highest_pt_dedup.py, predict.py) needs no change. The two tuples are
# told apart by their output directory, not by the tree path.
if options.TUPLIZER == "H_phipimu_GenTag":
    # No doTruth/channel: this analyzer has no signal gen filter at all. It tags what the
    # three reconstructed tracks actually are, which is orthogonal to the H->Ds(*)W topology.
    process.Hphipimu = cms.EDAnalyzer("H_phipimu_GenTag",
        sameSignPi = cms.bool(sameSign)
    )
else:
    process.Hphipimu = cms.EDAnalyzer("H_phipimu",
        sameSignPi = cms.bool(sameSign),
        doTruth = cms.bool(False),
        channel = cms.string("DsW")
    )
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

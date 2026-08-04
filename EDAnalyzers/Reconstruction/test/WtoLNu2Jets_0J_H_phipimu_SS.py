import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing('analysis')
options.register("JOBID", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Job ID")
options.register("GROUPSIZE", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Group size")
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.parseArguments()

start = options.JOBID * options.GROUPSIZE
end = start + options.GROUPSIZE

# WtoLNu-2Jets_0J is an official grid dataset (4333 files, UUID names), NOT
# replicated at T2_BE_IIHE (checked via dasgoclient site query) -- unlike the
# inclusive WtoLNu-2Jets sample. Read via the global AAA redirector instead of the
# IIHE-local one so xrootd can find whichever site actually holds each file. The
# file list is a static dasgoclient dump (dataset is published/fixed).
FILELIST = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WtoLNu2Jets/filelist_0J.txt"
with open(FILELIST) as fh:
    ALL_FILES = [line.strip() for line in fh if line.strip()]
INFILE = ["root://xrootd-cms.infn.it/" + f for f in ALL_FILES[start:end]]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("HphipimuWJets0JSS")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")
process.Hphipimu = cms.EDAnalyzer("H_phipimu_SS",
    doTruth = cms.bool(False),
    channel = cms.string("DsW")
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

import FWCore.ParameterSet.Config as cms

# The process object
process = cms.Process("MyNanoAOD")
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.noTimeStamps = True

# Import other attributes and functions
# from jerc-application-tutorial.ApplyOnMiniAOD.JetMETExtra_cff import *

# Input root files and number of events
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring("file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_1.root"),
    # fileNames = cms.untracked.vstring("/store/mc/Run3Summer22EEMiniAODv4/G-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/70000/257eacb5-2985-43fb-8db7-66256a019f0e.root")
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Output file
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName = cms.string("outputTree.root")

# Global tags
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("130X_mcRun3_2022_realistic_v5")

# Inputs for MyMiniAOD EDAnalyser
process.load("jerc-application-tutorial.ApplyOnMiniAOD.JecAppliedJetProducer_cfi")
process.load("jerc-application-tutorial.ApplyOnMiniAOD.JecAppliedMetProducer_cfi")
process.load("jerc-application-tutorial.ApplyOnMiniAOD.JvmAppliedEventFilter_cfi")
process.load("jerc-application-tutorial.ApplyOnMiniAOD.MyMiniAodAnalyzer_cfi")

process.myMiniAodAnalyzer.jetAK4Src = cms.InputTag("jecAppliedJetProducer:CorrectedAK4")
process.myMiniAodAnalyzer.metSrc = cms.InputTag("jecAppliedMetProducer:CorrectedMet")


# Add ED Filters, Producers, Analysers in the cms Path
process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.myMiniAodAnalyzer)
process.schedule = cms.Schedule(process.p)

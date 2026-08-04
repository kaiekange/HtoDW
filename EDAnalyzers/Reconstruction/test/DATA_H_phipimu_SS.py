import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

# source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("/store/data/Run2022C/Muon/MINIAOD/22Sep2023-v1/..."))
# OUTFILE = "file:output.root"

source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())
OUTFILE = "output.root"

process = cms.Process("HphipimuSS")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_dataRun3_v2", "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")

# Data JEC differs from MC: no JER smearing/JES-syst (both gated off by IsData in the
# plugin), plus the L2L3Residual correction, read from the "Era2022PreAll" payload in
# JecConfigAK4.json (the only data-era key under "2022Pre" -- covers Run2022C+D
# together, matching the "Run3-22CDSep23" JERC campaign the JSON points at).
# tightJetId/JvmAppliedEventFilter have no IsData knob -- same config for data and MC.
process.jecAppliedJetProducer.Jets.IsData = cms.bool(True)
process.jecAppliedJetProducer.Jets.Era    = cms.string("Era2022PreAll")
process.jecAppliedMetProducer.Jets.IsData = cms.bool(True)
process.jecAppliedMetProducer.Jets.Era    = cms.string("Era2022PreAll")
process.Hphipimu = cms.EDAnalyzer("H_phipimu_SS",
    doTruth = cms.bool(False),
    channel = cms.string("DsW"),
    # No electron ID overrides: this 22Sep2023 data reprocessing's pat::Electron
    # objects carry RunIIIWinter22-V1 (cut-based) and RunIIIWinter22-iso-V1 (MVA)
    # natively, matching the plugin's MC defaults exactly (verified via FWLite on
    # a real file -- 27Jun2023 lacks these IDs entirely, which is why data was
    # switched to 22Sep2023).
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra

# source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/03d57826-8e8d-4b11-8683-c75d0310b16f.root"))
# OUTFILE = "file:output.root"

source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())
OUTFILE = "output.root"

process = cms.Process("Hphipimu")
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

# All MC in this analysis is Run3Summer22MiniAODv4, i.e. preEE. This call is what makes
# that explicit -- the cfi default is empty and every JME plugin throws without it, so a
# future Summer22EE sample cannot silently inherit preEE JEC and veto map.
applyJmeEra(process, year="2022Pre", isData=False)
process.Hphipimu = cms.EDAnalyzer("H_phipimu_GenTag",
    sameSignPi = cms.bool(True)   # same-sign control region
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

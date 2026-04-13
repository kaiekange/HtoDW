import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_KKpimu/20260408/Simulation/PAT/2022/output_1.root"))

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer", maxEventsToPrint=cms.untracked.int32(-1), printVertex=cms.untracked.bool(False), src=cms.InputTag("prunedGenParticles"))

process.p = cms.Path(process.printTree)

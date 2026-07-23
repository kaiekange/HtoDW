import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:/afs/cern.ch/work/k/kakang/H2DW/CMSSW_12_4_14_patch3/src/Simulation/AODSIM/output.root"))

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleListDrawer", maxEventsToPrint=cms.untracked.int32(10), printVertex=cms.untracked.bool(False), src=cms.InputTag("genParticles"))

process.p = cms.Path(process.printTree)

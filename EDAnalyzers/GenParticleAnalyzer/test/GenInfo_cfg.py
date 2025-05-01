import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

process = cms.Process("GenPartAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

myinfile = sys.argv[2]
myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.prunedGenParticles = cms.EDAnalyzer('GenInfo',
    pruned = cms.InputTag('prunedGenParticles')
)

process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile))

process.p = cms.Path(process.prunedGenParticles)

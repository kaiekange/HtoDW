import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

process = cms.Process("GenPartAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

myinfile = "file:/user/kakang/H2cc/CMSSW_10_6_30_patch1/src/Simulation/2017/20241113/samples/process_" + sys.argv[2] + "/output/S5_PAT.root"
myoutfile = "file:/user/kakang/H2cc/CMSSW_10_6_30_patch1/src/Simulation/2017/20241113/analysis/chaincheck/tuple_" + sys.argv[2] + ".root"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.prunedGenParticles = cms.EDAnalyzer('GenParticleAnalyzer',
    pruned = cms.InputTag('prunedGenParticles')
)

print(myoutfile)
process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile))

process.p = cms.Path(process.prunedGenParticles)

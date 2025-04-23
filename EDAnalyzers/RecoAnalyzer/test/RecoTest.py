import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

process = cms.Process("GenPartAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

myinfile = "file:/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/FullGEN/PAT/output_1.root"
myoutfile = "file:test.root"
# myinfile = sys.argv[2]
# myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.reconstruction = cms.EDAnalyzer('RecoAnalyzer',
    packed = cms.InputTag('packedPFCandidates')
)

process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile))

process.p = cms.Path(process.reconstruction)

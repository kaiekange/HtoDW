import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from Configuration.AlCa.GlobalTag import GlobalTag
import sys

process = cms.Process("GenPartAnalysis")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag =  GlobalTag(process.GlobalTag, "106X_mc2017_realistic_v9")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

# myinfile = "file:/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/FullGEN/PAT/output_1.root"
# myoutfile = "file:test.root"
myinfile = sys.argv[2]
myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.PVStudy = cms.EDAnalyzer('PVStudy',
    prunedGenParticles = cms.InputTag('prunedGenParticles'),
    packedPFCandidates = cms.InputTag('packedPFCandidates')
)

process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile))

process.p = cms.Path(process.PVStudy)

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from Configuration.AlCa.GlobalTag import GlobalTag
import sys

process = cms.Process("MuonAnalysis")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag =  GlobalTag(process.GlobalTag, "106X_mc2017_realistic_v9")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

# myinfile = "file:/user/kakang/Analysis/CMSSW_10_6_30_patch1/src/tuples/PAT/output_1.root"
# myoutfile = "file:test.root"
myinfile = sys.argv[2]
myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.MuonStudy = cms.EDAnalyzer('MuonStudy',
    genPart = cms.InputTag('prunedGenParticles'),
    pfCands = cms.InputTag('packedPFCandidates'),
    muons = cms.InputTag('slimmedMuons'),
    primvtx = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamspot = cms.InputTag("offlineBeamSpot"),
)

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string(myoutfile), closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.MuonStudy)

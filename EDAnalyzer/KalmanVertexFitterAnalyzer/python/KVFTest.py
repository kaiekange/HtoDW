import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

process = cms.Process("KVFTest")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag =  GlobalTag(process.GlobalTag, "106X_mc2017_realistic_v6")

process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.Services_cff")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) ) 

myinfile = sys.argv[2]
#myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(myinfile)
)

process.KVFTest = cms.EDAnalyzer('KalmanVertexFitterAnalyzer',
    packed = cms.InputTag('packedPFCandidates'),
    KVFParameters = cms.PSet(
        maxDistance = cms.double(0.01),
        maxNbrOfIterations = cms.int32(10)
    )
)

#process.TFileService = cms.Service(
#    "TFileService", fileName=cms.string(myoutfile)
#)

process.p = cms.Path(process.KVFTest)

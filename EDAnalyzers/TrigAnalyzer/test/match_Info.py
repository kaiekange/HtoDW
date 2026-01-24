import sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

myinfile = sys.argv[2]
myoutfile = sys.argv[3]

process = cms.Process("TriggerTest")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(myinfile))

process.TriggerTest = cms.EDAnalyzer("MatchAnalyzer")

# process.residuals.stageL1Trigger = cms.uint32(2)

process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile), closeFileFast=cms.untracked.bool(True))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.p = cms.Path(process.TriggerTest)

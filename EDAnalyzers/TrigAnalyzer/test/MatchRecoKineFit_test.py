import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

file_list = ["file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_{}.root".format(i) for i in range(1, 1001) if i != 185]
# file_list = ["file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_{}.root".format(i) for i in range(1, 21)]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*file_list))

# source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_1.root"))

process = cms.Process("MatchRecoKineFit")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.MatchRecoKineFit = cms.EDAnalyzer("MatchRecoKineFit")

# process.residuals.stageL1Trigger = cms.uint32(2)

process.TFileService = cms.Service("TFileService", fileName=cms.string("file:/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/MatchRecoKineFit/MatchRecoKineFit.root"), closeFileFast=cms.untracked.bool(True))
# process.TFileService = cms.Service("TFileService", fileName=cms.string("file:./Kinetest.root"), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.p = cms.Path(process.MatchRecoKineFit)

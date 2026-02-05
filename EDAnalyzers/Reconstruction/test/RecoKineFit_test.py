import sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

JOBID = int(sys.argv[2])
start = JOBID * 20 + 1
end = start + 20

myoutfile = sys.argv[3]
# myoutfile = "file:./Kinetest.root"
# myoutfile = "file:/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/RecoKine/RecoKine.root"

file_list = ["file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_{}.root".format(i) for i in range(start, end) if i != 185]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*file_list))

process = cms.Process("RecoKineFit")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.RecoKineFit = cms.EDAnalyzer("RecoKineFit")

# process.residuals.stageL1Trigger = cms.uint32(2)

process.TFileService = cms.Service("TFileService", fileName=cms.string(myoutfile), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.p = cms.Path(process.RecoKineFit)

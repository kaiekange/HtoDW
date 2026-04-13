import sys
import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

# INFILE = ["file:/eos/home-k/kakang/H2DW/Datasets/20260119/Simulation/PAT/2022/output_1.root"]
# OUTFILE = "signal_jettest.root"
INFILE = ["/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v3/2540000/047c5864-607c-4943-a200-d5d26f8a56a0.root"]
OUTFILE = "W0J_jettest.root"
# INFILE = ["/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2540000/06118768-c07c-43c4-9aec-790384565132.root"]
# OUTFILE = "W1J_jettest.root"
# INFILE = ["/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2540000/14e868f4-f92a-431f-8b10-e9932878c3e4.root"]
# OUTFILE = "W2J_jettest.root"
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("CheckJet")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.Reconstruction.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.Reconstruction.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.Reconstruction.JvmAppliedEventFilter_cfi")
process.CheckJet = cms.EDAnalyzer("CheckJet")

process.TFileService = cms.Service("TFileService", fileName=cms.string(OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Run the JEC/JVM producers and filter before the CheckJet analyzer
process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.CheckJet)
process.schedule = cms.Schedule(process.p)

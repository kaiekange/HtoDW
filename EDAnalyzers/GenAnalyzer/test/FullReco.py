import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing('analysis')
options.register('INFILE', 'file:input.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Input file')
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.parseArguments()

input_files = [f.strip() for f in options.INFILE.split() if f.strip()]

source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(input_files))

process = cms.Process("FullReco")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.GenAnalyzer.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.GenAnalyzer.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.GenAnalyzer.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.GenAnalyzer.tightJetId_cfi")
process.FullReco = cms.EDAnalyzer("FullReco")

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Run the JEC/JVM producers and filter before the FullReco analyzer
process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.FullReco)
process.schedule = cms.Schedule(process.p)

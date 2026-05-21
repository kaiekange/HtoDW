import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing('analysis')
options.register("JOBID", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Job ID")
options.register("GROUPSIZE", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Group size")
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.parseArguments()

start = options.JOBID * options.GROUPSIZE + 1
end = start + options.GROUPSIZE

INFILE = ["root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_KKpimu/20260408/Simulation/PAT/2022/output_{}.root".format(i) for i in range(start, end)]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("DsstWKKpimuonDsW")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_v5", "")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")
process.DsstWKKpimuonDsW = cms.EDAnalyzer("DsstW_KKpimu_onDsW")
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Run the JEC/JVM producers and filter before the DsstWKKpimuonDsW analyzer
process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.DsstWKKpimuonDsW)
process.schedule = cms.Schedule(process.p)

import re
import subprocess
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

# old sample (20260624): W->mu/e only, Ds*->gamma only generation fragment. Superseded by the
# 20260717 inclusive W->e/mu/tau + Ds*->gamma/pi0 sample below; kept here for reference.
# INFILE = ["root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_phipi/20260624/Simulation/PAT/2022/output_{}.root".format(i) for i in range(start, end)]

# 20260717 production is still running (not all 500 output files exist yet), so check what's
# actually on storage and only queue the ones present -- checked live so a rerun later picks up
# files that finished in the meantime.
PATDIR = "/store/user/kakang/Analysis/ggHtoDsW_phipi/20260717/Simulation/PAT/2022"
ls_out = subprocess.check_output(["xrdfs", "maite.iihe.ac.be", "ls", PATDIR]).decode()
existing = set(int(m) for m in re.findall(r"output_(\d+)\.root", ls_out))
INFILE = ["root://maite.iihe.ac.be/{}/output_{}.root".format(PATDIR, i) for i in range(start, end) if i in existing]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("HphipimuDsWSS")
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
process.HphipimuDsW = cms.EDAnalyzer("H_phipimu_SS",
    doTruth = cms.bool(True),
    channel = cms.string("DsW")
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.HphipimuDsW)
process.schedule = cms.Schedule(process.p)

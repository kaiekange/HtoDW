import re
import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra

options = VarParsing('analysis')
options.register("JOBID", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Job ID")
options.register("GROUPSIZE", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Group size")
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.register('YEAR', '2022', VarParsing.multiplicity.singleton, VarParsing.varType.string, '2022 (preEE) or 2022EE (postEE)')
options.parseArguments()

start = options.JOBID * options.GROUPSIZE + 1
end = start + options.GROUPSIZE

# old sample (20260624): W->mu/e only, Ds*->gamma only generation fragment. Superseded by the
# 20260717 inclusive W->e/mu/tau + Ds*->gamma/pi0 sample below; kept here for reference.
# INFILE = ["root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsstW_phipi/20260624/Simulation/PAT/2022/output_{}.root".format(i) for i in range(start, end) if i != 188]

# Single YEAR knob drives ver/GT/JME era together so they can't desync (same principle as
# applyJmeEra and the grid-production mkconfig.sh scripts). 2022EE (postEE) added 2026-09-01.
ERA = {
    "2022":   dict(ver="20260717", gt="130X_mcRun3_2022_realistic_v5",         jme="2022Pre"),
    "2022EE": dict(ver="20260825", gt="130X_mcRun3_2022_realistic_postEE_v6", jme="2022Post"),
}[options.YEAR]

# Neither era reaches all 1000 files (6 AODSIM outputs are broken per era and will never
# produce PAT) -- check what's actually on storage and only queue the ones present, so this
# is self-healing on reruns as late files trickle in (though those 6 indices never will).
PATDIR = "/store/user/kakang/Analysis/ggHtoDsstW_phipi/{}/Simulation/PAT/{}".format(ERA["ver"], options.YEAR)
ls_out = subprocess.check_output(["xrdfs", "maite.iihe.ac.be", "ls", PATDIR]).decode()
existing = set(int(m) for m in re.findall(r"output_(\d+)\.root", ls_out))
INFILE = ["root://maite.iihe.ac.be/{}/output_{}.root".format(PATDIR, i) for i in range(start, end) if i in existing]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("Hphipimu")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, ERA["gt"], "")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")

# cfi default is empty and every JME plugin throws without it, so a mismatched YEAR
# cannot silently inherit the wrong JEC/veto map.
applyJmeEra(process, year=ERA["jme"], isData=False)
process.Hphipimu = cms.EDAnalyzer("H_phipimu",
    doTruth = cms.bool(True),
    channel = cms.string("DsstW")
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

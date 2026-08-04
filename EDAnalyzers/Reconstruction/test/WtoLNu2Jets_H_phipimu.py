"""WtoLNu2Jets tuplizer cfg, parametrized over jet bin and OS/SS.

Replaces the eight near-identical cfgs this used to be (inclusive/0J/1J/2J x OS/SS),
which differed only in the file list, the xrootd redirector and the process name.

    cmsRun WtoLNu2Jets_H_phipimu.py JETBIN=0J SS=1 JOBID=0 GROUPSIZE=20 OUTFILE=...

JETBIN is one of incl/0J/1J/2J; SS=1 selects the same-sign control region. Process
names are unchanged from the per-bin cfgs, and the module label stays `Hphipimu`
in every combination, so the tree path is always Hphipimu/Events -- already
produced files and every downstream script keep working.

The filelist*.txt dumps are gitignored (*.txt) and must exist locally; they are
static dasgoclient dumps, the datasets being published and fixed (unlike the
still-running private signal productions, which need a live existence check).
"""
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag

options = VarParsing('analysis')
options.register("JOBID", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Job ID")
options.register("GROUPSIZE", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Group size")
options.register('OUTFILE', 'file:output.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output file')
options.register('JETBIN', 'incl', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'incl | 0J | 1J | 2J')
# int, not bool: VarParsing's bool type reads SS=0 as truthy, which would silently
# write OS data into an SS output directory.
options.register('SS', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'same-sign control region (0/1)')
options.parseArguments()

# Per dataset: the file-list suffix, and which redirector can serve it. Only the
# inclusive sample is replicated at T2_BE_IIHE (checked via dasgoclient site query);
# the three jet-binned ones are not, so they go through the global AAA redirector,
# which finds whichever site actually holds each file.
DATASETS = {
    "incl": ("",    "root://maite.iihe.ac.be/"),
    "0J":   ("_0J", "root://xrootd-cms.infn.it/"),
    "1J":   ("_1J", "root://xrootd-cms.infn.it/"),
    "2J":   ("_2J", "root://xrootd-cms.infn.it/"),
}
if options.JETBIN not in DATASETS:
    sys.exit(f"JETBIN='{options.JETBIN}' is not one of {sorted(DATASETS)}")
suffix, redirector = DATASETS[options.JETBIN]
sameSign = bool(options.SS)

start = options.JOBID * options.GROUPSIZE
end = start + options.GROUPSIZE

FILELIST = f"/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WtoLNu2Jets/filelist{suffix}.txt"
with open(FILELIST) as fh:
    ALL_FILES = [line.strip() for line in fh if line.strip()]
INFILE = [redirector + f for f in ALL_FILES[start:end]]
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INFILE))

process = cms.Process("HphipimuWJets" + suffix.lstrip("_") + ("SS" if sameSign else ""))
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

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")
process.Hphipimu = cms.EDAnalyzer("H_phipimu",
    sameSignPi = cms.bool(sameSign),
    doTruth = cms.bool(False),
    channel = cms.string("DsW")
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

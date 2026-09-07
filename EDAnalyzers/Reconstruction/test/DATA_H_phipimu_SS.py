import sys

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing
from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra

# Local test:
#   cmsRun DATA_H_phipimu_SS.py GLOBALTAG=130X_dataRun3_v2 YEAR=2022Pre \
#          inputFiles=root://cms-xrd-global.cern.ch//store/data/... OUTFILE=file:output.root

source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())

options = VarParsing("analysis")
options.register("GLOBALTAG", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Conditions GlobalTag, e.g. 130X_dataRun3_v2")
options.register("YEAR", "", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JEC / jet-veto-map era key: 2022Pre | 2022Post")
options.register("OUTFILE", "output.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output file")
options.parseArguments()
OUTFILE = options.OUTFILE

# Local testing only. Left empty under CRAB, which injects the per-job file list itself;
# the guard keeps the expanded config identical to what it was before this was added.
if options.inputFiles:
    source.fileNames = cms.untracked.vstring(options.inputFiles)

# GLOBALTAG and YEAR are per-dataset and come from Tuples/Data/DatasetList/data_2022.txt
# (one line per dataset) via submitData.sh -> crabConfigTemplateData.py -> pyCfgParams.
# That table is deliberately NOT duplicated here: a copy in each of the two data psets
# is the same drift hazard that the single FEATURE_NAMES source in analysis/DNN/common.py
# exists to prevent. It is also why OS and SS can no longer drift apart on era -- both
# read the same line of the same file.
#
# The two are INDEPENDENT axes, which is why the txt carries both as separate columns
# rather than one pre/postEE flag: Run2022E is postEE for JEC and the jet veto map but
# shares Run2022C/D's GlobalTag, while F and G report 130X_dataRun3_PromptAnalysis_v1
# (they were not re-reconstructed in the 22Sep2023 campaign the way C/D/E were). Both
# columns were read off DAS, not inferred:
#   dasgoclient -query="config dataset=<ds> system=dbs3 | grep config.global_tag"
#
# YEAR -> JEC data-era key. Each Year has exactly one data-era key in JecConfigAK4.json,
# so this map is a property of that file rather than of any dataset -- it stays here
# instead of becoming a fourth txt column that could disagree with YEAR.
JEC_ERA_KEY = {
    "2022Pre":  "Era2022PreAll",   # Run2022C+D,   JERC campaign Run3-22CDSep23-Summer22
    "2022Post": "Era2022PostAll",  # Run2022E+F+G, campaign Run3-22EFGSep23-Summer22EE
}
# Both required, no defaults: a silently inherited era is the exact failure this
# parametrisation exists to prevent, and a bad GlobalTag otherwise surfaces only when a
# job reaches Frontier, potentially hours in.
if options.YEAR not in JEC_ERA_KEY:
    raise ValueError("YEAR must be one of %s, got '%s'" % (sorted(JEC_ERA_KEY), options.YEAR))
if not options.GLOBALTAG:
    raise ValueError("GLOBALTAG is required (see Tuples/Data/DatasetList/data_2022.txt)")
GLOBALTAG, JEC_YEAR, JEC_ERA = options.GLOBALTAG, options.YEAR, JEC_ERA_KEY[options.YEAR]
# stderr, not stdout: edmConfigDump writes the expanded config to stdout, and a stray
# print there makes that dump non-reexecutable.
sys.stderr.write("[DATA pset] GlobalTag=%s  Year=%s  JecEra=%s\n"
                 % (GLOBALTAG, JEC_YEAR, JEC_ERA))

process = cms.Process("Hphipimu")
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag = GlobalTag(process.GlobalTag, GLOBALTAG, "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = source

process.load("EDAnalyzers.JMEcorr.JecAppliedJetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JecAppliedMetProducer_cfi")
process.load("EDAnalyzers.JMEcorr.JvmAppliedEventFilter_cfi")
process.load("EDAnalyzers.JMEcorr.tightJetId_cfi")

# One call sets Year on all three JME modules (both JEC producers and the jet veto map
# filter) plus the data-only IsData/Era. Doing it by hand is what this function exists to
# prevent: the veto map has no run-number dimension, so a Year set on two modules out of
# three is silently wrong rather than an error.
#
# Data JEC differs from MC: no JER smearing/JES-syst (both gated off by IsData in the
# plugin), plus the L2L3Residual correction. Each Year has exactly one data-era key in
# JecConfigAK4.json -- "Era2022PreAll" under "2022Pre" (covers Run2022C+D, JERC campaign
# Run3-22CDSep23-Summer22) and "Era2022PostAll" under "2022Post" (covers Run2022E+F+G,
# campaign Run3-22EFGSep23-Summer22EE).
# tightJetId has neither an IsData nor a Year knob -- same config for every era and MC.
applyJmeEra(process, year=JEC_YEAR, isData=True, era=JEC_ERA)
process.Hphipimu = cms.EDAnalyzer("H_phipimu",
    sameSignPi = cms.bool(True),   # same-sign control region (merged into the OS plugin)
    doTruth = cms.bool(False),
    channel = cms.string("DsW"),
    # No electron ID overrides: this 22Sep2023 data reprocessing's pat::Electron
    # objects carry RunIIIWinter22-V1 (cut-based) and RunIIIWinter22-iso-V1 (MVA)
    # natively, matching the plugin's MC defaults exactly (verified via FWLite on
    # a real file -- 27Jun2023 lacks these IDs entirely, which is why data was
    # switched to 22Sep2023).
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(OUTFILE), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

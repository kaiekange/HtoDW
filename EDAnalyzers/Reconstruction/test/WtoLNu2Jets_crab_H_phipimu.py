"""INCLUSIVE WtoLNu-2Jets tuplizer cfg for CRAB, parametrized over era and OS/SS.

Sibling of WtoLNu2Jets_H_phipimu.py. The two differ on ONE axis -- where the input
files come from -- and that is why they are separate files rather than one cfg with
another switch:

    WtoLNu2Jets_H_phipimu.py       HTCondor. Reads filelist_<bin>_<era>.txt locally and
                                   slices it by JOBID/GROUPSIZE. Jet-binned (0J/1J/2J).
    WtoLNu2Jets_crab_H_phipimu.py  CRAB. Source is left empty; CRAB fills fileNames per
                                   job from Data.inputDataset. Inclusive (no jet bin).

    crab submit -c Tuples/WtoLNu2Jets/crabConfig_incl_2022.py

ERA is 2022 (Run3Summer22, preEE) or 2022EE (Run3Summer22EE, postEE); SS=1 selects the
same-sign control region. Both are passed through JobType.pyCfgParams. Process name and
module label are both `Hphipimu`, so the tree path is always Hphipimu/Events -- which is
what highest_pt_dedup.py and every downstream consumer expect.

The inclusive datasets are the union of the three jet bins, i.e. the SAME physical
process the jet-binned samples cover. They were dropped on 2026-08-18 in favour of the
bins; this cfg exists to produce them again through CRAB.
"""
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra

options = VarParsing('analysis')
# int, not bool: VarParsing's bool type reads SS=0 as truthy, which would silently
# write OS data into an SS output directory.
options.register('SS', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'same-sign control region (0/1)')
# No default: the two eras need different GlobalTags AND different JEC/veto-map eras, and
# getting that wrong produces plausible numbers rather than an error. Required, exactly
# like GLOBALTAG/YEAR in the data cfgs.
options.register('ERA', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 '2022 (Run3Summer22, preEE) | 2022EE (Run3Summer22EE, postEE)')
options.parseArguments()

# Era -> (GlobalTag, JME era key). Two things that must move together:
#   2022    Run3Summer22MiniAODv4    130X_mcRun3_2022_realistic_v5          JEC 2022Pre
#   2022EE  Run3Summer22EEMiniAODv4  130X_mcRun3_2022_realistic_postEE_v6   JEC 2022Post
# The JME half is guarded: the cfi Year default is empty and all three JME modules throw on
# it, so a 2022EE run that forgot to switch cannot silently use preEE JEC or the preEE jet
# veto map. The GlobalTag half has no such guard -- hence no default here either.
ERAS = {
    "2022":   ("130X_mcRun3_2022_realistic_v5",        "2022Pre"),
    "2022EE": ("130X_mcRun3_2022_realistic_postEE_v6", "2022Post"),
}
if options.ERA not in ERAS:
    sys.exit(f"ERA='{options.ERA}' is not one of {sorted(ERAS)}")
GLOBALTAG, JEC_YEAR = ERAS[options.ERA]

sameSign = bool(options.SS)

# Left EMPTY on purpose: CRAB overwrites process.source.fileNames per job from
# Data.inputDataset. No file list, no redirector, no JOBID/GROUPSIZE slicing -- that is
# the whole difference from the HTCondor sibling. Running this cfg by hand with cmsRun
# therefore reads nothing; use WtoLNu2Jets_H_phipimu.py for a local test.
source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())

# Nothing downstream reads the process name (the tree path follows the module label), and a
# composed one must stay alphanumeric. Fixed value; era and sign live in the output path.
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

# Sets Year on all three JME modules from the era chosen above. The cfi default is empty
# and every JME plugin throws without it, so a 2022EE run cannot silently inherit the preEE
# JEC payloads or the preEE jet veto map.
applyJmeEra(process, year=JEC_YEAR, isData=False)
process.Hphipimu = cms.EDAnalyzer("H_phipimu",
    sameSignPi = cms.bool(sameSign),
    doTruth = cms.bool(False),
    channel = cms.string("DsW")
)
# Plain literal, NOT a VarParsing option: CRAB does not pass OUTFILE, it collects whatever
# the job wrote by this name. A "file:" prefix here would make the written name and the
# name CRAB looks for disagree, which shows up as a task that finishes with no output.
process.TFileService = cms.Service("TFileService", fileName=cms.string("output.root"), closeFileFast=cms.untracked.bool(True))

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.p = cms.Path(process.jecAppliedJetProducer * process.jecAppliedMetProducer * process.jvmAppliedEventFilter * process.tightJetId * process.Hphipimu)
process.schedule = cms.Schedule(process.p)

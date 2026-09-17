"""Custom NanoAOD production cfg for 2022 Muon data: standard bare NANO + the H -> Ds(*)W
phi-pi candidate table. Submitted via CRAB (../crab/submitData.sh), which injects the
per-job file list and applies the Golden-JSON lumi mask itself (config.Data.lumiMask) --
this config never lists input files or a lumi mask on its own.

Local smoke test:
    cmsRun NANO_DATA_H_phipimu.py GLOBALTAG=130X_dataRun3_v2 YEAR=2022Pre OUTFILE=file:nano_data.root \
        inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2022C/Muon/MINIAOD/22Sep2023-v1/...

Knobs (VarParsing 'analysis'):
  GLOBALTAG  required, no default -- per-dataset, from Tuples/Data/DatasetList/data_2022.txt
             (C/D/E: 130X_dataRun3_v2; F/G: 130X_dataRun3_PromptAnalysis_v1 -- GT and the
             preEE/postEE output split are independent axes, see that file's header comment)
  OUTFILE    output file (file:...)
  SS         0 | 1   -- same-sign-pi control region (int, not bool: bool reads SS=0 as truthy)
  NTHREADS   cmsRun threads
  YEAR       2022Pre | 2022Post, required -- Year column of data_2022.txt; picks the ParticleNet WPs

doTruth is forced False and only the nominal (non-GenTag) tuplizer runs: GenTag reads
prunedGenParticles, which doesn't exist on data. No JOBID/GROUPSIZE/FILELIST here -- CRAB
does its own LumiBased splitting and overrides process.source.fileNames per job.
"""
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3

options = VarParsing('analysis')
options.register('GLOBALTAG', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'Conditions GlobalTag, e.g. 130X_dataRun3_v2 (required, no default)')
options.register('OUTFILE', 'file:nano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'output file')
options.register('SS', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'same-sign-pi control region (0/1)')
options.register('NTHREADS', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'cmsRun threads')
options.register('YEAR', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 '2022Pre | 2022Post (required) -- the Year column of data_2022.txt; selects the ParticleNet WPs')
options.parseArguments()

if not options.GLOBALTAG:
    sys.exit("GLOBALTAG is required (see Tuples/Data/DatasetList/data_2022.txt)")
_YEAR_TO_ERA = {"2022Pre": "2022", "2022Post": "2022EE"}
if options.YEAR not in _YEAR_TO_ERA:
    sys.exit("YEAR=%r must be one of %s (Year column of data_2022.txt)" % (options.YEAR, sorted(_YEAR_TO_ERA)))

process = cms.Process('NANO', Run3)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),  # empty under CRAB; it injects this
    secondaryFileNames=cms.untracked.vstring(),
)
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(options.NTHREADS),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary=cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm=cms.untracked.string('LZMA'),
    compressionLevel=cms.untracked.int32(9),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string('NANOAOD'),
        filterName=cms.untracked.string(''),
    ),
    fileName=cms.untracked.string(options.OUTFILE),
    outputCommands=process.NANOAODEventContent.outputCommands,
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GLOBALTAG, '')

process.nanoAOD_step = cms.Path(process.nanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon
process = nanoAOD_customizeCommon(process)

# --- custom H -> Ds(*)W phi-pi candidate table (nominal tuplizer only, doTruth forced False) ---
process.load('EDAnalyzers.Reconstruction.nano_phipimu_cff')
from EDAnalyzers.Reconstruction.nano_phipimu_cff import customise
process.phipimuCand.sameSignPi = cms.bool(bool(options.SS))
process.phipimuCand.doTruth = cms.bool(False)
process = customise(process, tuplizer="nominal", is_data=True)
from EDAnalyzers.Reconstruction.nano_phipimu_cff import set_btag_wp
process = set_btag_wp(process, _YEAR_TO_ERA[options.YEAR])  # era ParticleNet WPs, required

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

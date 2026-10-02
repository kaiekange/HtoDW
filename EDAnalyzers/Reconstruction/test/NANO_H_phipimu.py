"""Custom NanoAOD production cfg: standard NANO + the H -> Ds(*)W phi-pi candidate table.

    cmsRun NANO_H_phipimu.py SAMPLE=DsW YEAR=2022 JOBID=0 GROUPSIZE=1 OUTFILE=file:nano.root

Knobs (VarParsing 'analysis'):
  SAMPLE    DsW | DsstW | WJ0-2 | TTto2L2Nu | TTtoLNu2Q | DY0J-2J   (or leave empty and pass inputFiles=...)
  YEAR      2022 (Run3Summer22, preEE) | 2022EE (Run3Summer22EE, postEE)
  JOBID     job index; GROUPSIZE files per job -> [JOBID*GROUPSIZE, +GROUPSIZE)
  OUTFILE   output file (file:...)
  SS        0 | 1   -- same-sign-pi control region (int, not bool: bool reads SS=0 as truthy)
  TUPLIZER  nominal | gentag   -- HDsPhiPiMuProducer vs HDsPhiPiMuGenTagProducer

Input files come from FILELIST (one root:// URL per line, sliced by JOBID/GROUPSIZE) or
explicit inputFiles=. condor_submit.py builds the FILELIST at submit time (xrootd ls for
signal, the static dasgoclient dumps for W+jets/ttbar/DY). Data is not handled yet (per-era-subset GT).
"""
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3

options = VarParsing('analysis')
options.register('JOBID', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'job index')
options.register('GROUPSIZE', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'files per job')
options.register('OUTFILE', 'file:nano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'output file')
options.register('SAMPLE', 'DsW', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'DsW|DsstW|WJ0|WJ1|WJ2 (empty -> use inputFiles=)')
options.register('YEAR', '2022', VarParsing.multiplicity.singleton, VarParsing.varType.string, '2022 (preEE) | 2022EE (postEE)')
options.register('SS', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'same-sign-pi control region (0/1)')
options.register('TUPLIZER', 'nominal', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'nominal | gentag')
options.register('NTHREADS', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'cmsRun threads (match condor request_cpus)')
options.register('FILELIST', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'file with one root:// URL per line, sliced by JOBID/GROUPSIZE (0-based)')
options.parseArguments()

# Era -> (signal production version, MC GlobalTag). The version only tracks preEE/postEE
# (both channels: preEE = 20260717/.../2022, postEE = 20260825/.../2022EE).
ERAS = {
    "2022":   ("20260717", "130X_mcRun3_2022_realistic_v5"),
    "2022EE": ("20260825", "130X_mcRun3_2022_realistic_postEE_v6"),
}
if options.YEAR not in ERAS:
    sys.exit("YEAR=%r not in %s" % (options.YEAR, sorted(ERAS)))
_ver, GLOBALTAG = ERAS[options.YEAR]  # _ver: informational, condor_submit.py resolves files

if options.TUPLIZER not in ("nominal", "gentag"):
    sys.exit("TUPLIZER=%r must be 'nominal' or 'gentag'" % options.TUPLIZER)

# Input files. File-list RESOLUTION (xrootd ls / DAS) is done at submit time by
# condor_submit.py, on lxplus where the tools work -- NOT here: `xrdfs` on an el8 worker
# fails to load libreadline.so.7, and a subprocess in a cmsRun config is fragile anyway.
# This config just slices an already-resolved FILELIST (or takes explicit inputFiles=).
if options.inputFiles:
    infiles = list(options.inputFiles)
elif options.FILELIST:
    with open(options.FILELIST) as fh:
        allf = [ln.strip() for ln in fh if ln.strip()]
    s = options.JOBID * options.GROUPSIZE
    infiles = allf[s:s + options.GROUPSIZE]
else:
    sys.exit("give FILELIST=<path> or inputFiles=...")

if not infiles:
    sys.exit("no input files for JOBID=%d GROUPSIZE=%d (list has %d)"
             % (options.JOBID, options.GROUPSIZE, len(open(options.FILELIST).readlines()) if options.FILELIST else 0))

process = cms.Process('NANO', Run3)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(infiles),
    secondaryFileNames=cms.untracked.vstring(),
)
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(options.NTHREADS),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary=cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm=cms.untracked.string('LZMA'),
    compressionLevel=cms.untracked.int32(9),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string('NANOAODSIM'),
        filterName=cms.untracked.string(''),
    ),
    fileName=cms.untracked.string(options.OUTFILE),
    outputCommands=process.NANOAODSIMEventContent.outputCommands,
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GLOBALTAG, '')

process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon
process = nanoAOD_customizeCommon(process)

# --- custom H -> Ds(*)W phi-pi candidate table ---
process.load('EDAnalyzers.Reconstruction.nano_phipimu_cff')
from EDAnalyzers.Reconstruction.nano_phipimu_cff import customise
if options.TUPLIZER == "gentag":
    process.phipimuCandGenTag.sameSignPi = cms.bool(bool(options.SS))
else:
    process.phipimuCand.sameSignPi = cms.bool(bool(options.SS))
    if options.SAMPLE == "DsstW":
        process.phipimuCand.channel = cms.string("DsstW")
    process.phipimuCand.doTruth = cms.bool(options.SAMPLE in ("DsW", "DsstW"))
process = customise(process, tuplizer=options.TUPLIZER)
from EDAnalyzers.Reconstruction.nano_phipimu_cff import set_btag_wp
process = set_btag_wp(process, options.YEAR)  # era ParticleNet WPs, required

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

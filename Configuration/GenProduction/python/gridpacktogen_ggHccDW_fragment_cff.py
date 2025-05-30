import FWCore.ParameterSet.Config as cms

# link to datacards: 
# https://github.com/Soumyatifr/genproductions/blob/e77039bd5ac909398a612c12ce10a6efbc2c37d0/bin/Powheg/production/2017/13TeV/Higgs/HJ_MiNLO_NNLOPS_NNPDF31_13TeV/HJ_MiNLO_NNLOPS_NNPDF31_13TeV.input

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/UL/13TeV/powheg/V2/gg_H_quark-mass-effects_slc7_amd64_gcc700_CMSSW_10_6_30_patch1_NNPDF31_13TeV_UL_M125/v2/gg_H_quark-mass-effects_slc7_amd64_gcc700_CMSSW_10_6_30_patch1_NNPDF31_13TeV_UL_M125.tgz'),
    nEvents = cms.untracked.uint32(5000),
    generateConcurrently = cms.untracked.bool(True),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

generator = cms.EDFilter("Pythia8ConcurrentHadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,

        processParameters = cms.vstring(
        'POWHEG:nFinal = 1',
        '25:onMode = off',
        '25:addChannel = 1  1.00   103   -24   431',
	    '25:onIfMatch = -24 431', # custom addition based on twiki page, see notes
        '25:m0 = 125.0',
        '24:onMode = off',
        '24:onIfMatch = -13 14',
	    # custom decay of Ds meson:
	    '431:onMode = off',
	    '431:onIfMatch = 333 211',
	    '333:onMode = off',
	    '333:onIfMatch = 321 -321',
        ),

        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8PSweightsSettings',
                                    'pythia8PowhegEmissionVetoSettings',
                                    'processParameters'
                                    )
        )
                         )

ProductionFilterSequence = cms.Sequence(generator)

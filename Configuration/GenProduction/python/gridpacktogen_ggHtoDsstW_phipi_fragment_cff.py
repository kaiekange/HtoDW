import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/RunIII/13p6TeV/slc7_amd64_gcc10/powheg/V2/gg_H_quark-mass-effects_mwindow1d0_slc7_amd64_gcc10_CMSSW_12_4_8.tgz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    generateConcurrently = cms.untracked.bool(False),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13600.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,
        processParameters = cms.vstring(
            'POWHEG:nFinal = 1',
            '25:onMode = off',
            '25:addChannel = 1  1.00   103   -24   433',
            '25:addChannel = 1  1.00   103   24   -433',
            '25:m0 = 125.0',
            '433:onMode = off',
            '24:onMode = off',
            '24:onIfMatch = -11 12',
            '24:onIfMatch = -13 14',
            '24:onIfMatch = -15 16'
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8PSweightsSettings',
                                    'pythia8PowhegEmissionVetoSettings',
                                    'processParameters',
        )
    ),

    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            verbose = cms.untracked.int32(1),
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            list_forced_decays = cms.vstring('myD_s*+','myD_s*-'),
            operates_on_particles = cms.vint32(433,-433),
            convertPythiaCodes = cms.untracked.bool(True),
            user_decay_embedded= cms.vstring(
"""
Alias myD_s*+ D_s*+
Alias myD_s*- D_s*-
Alias myD_s+ D_s+
Alias myD_s- D_s-
Alias myphi phi
Alias myK+ K+
Alias myK- K-
Alias mypi+ pi+
Alias mypi- pi-
Alias mypi0 pi0
Alias mygamma gamma
ChargeConj myD_s*+ myD_s*-
ChargeConj myD_s+ myD_s-
ChargeConj myphi myphi
ChargeConj myK+ myK-
ChargeConj mypi+ mypi-
ChargeConj mypi0 mypi0
ChargeConj mygamma mygamma

Decay myD_s*+
0.942000000 myD_s+ mygamma VSP_PWAVE;
0.058000000 myD_s+ mypi0 VSS;
Enddecay
CDecay myD_s*-

Decay myD_s+
1.0 myphi mypi+ SVS;
Enddecay
CDecay myD_s-

Decay myphi
1.0 myK+ myK- VSS;
Enddecay

Decay mypi0
1.0 mygamma mygamma PHSP;
Enddecay

End
"""
            ),
        ),
        parameterSets = cms.vstring('EvtGen130'),
    )
)

ProductionFilterSequence = cms.Sequence(generator)

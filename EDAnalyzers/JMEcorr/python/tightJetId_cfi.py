import FWCore.ParameterSet.Config as cms

tightJetId = cms.EDProducer(
    "PatJetIDValueMapProducer",
    filterParams=cms.PSet(
        version=cms.string("RUN3PUPPIruns2022FGruns2023CD"),
        quality=cms.string("TIGHTLEPVETO"),
    ),
    src=cms.InputTag("jecAppliedJetProducer", "CorrectedAK4"),
)

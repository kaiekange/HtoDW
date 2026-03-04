import FWCore.ParameterSet.Config as cms

myMiniAodAnalyzer = cms.EDAnalyzer("MyMiniAodAnalyzer",
    jetAK4Src = cms.InputTag("slimmedJetsPuppi"),
    jetAK8Src = cms.InputTag(""),   # set to "" to skip AK8
    metSrc    = cms.InputTag("slimmedMETsPuppi")       # or "slimmedMETsPuppi"
)


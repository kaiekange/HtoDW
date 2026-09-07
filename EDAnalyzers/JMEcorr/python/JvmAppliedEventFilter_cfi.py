
import FWCore.ParameterSet.Config as cms

jvmAppliedEventFilter = cms.EDFilter(
    "JvmAppliedEventFilter",
    Jets = cms.PSet(
        SourcesAK4 = cms.InputTag("slimmedJetsPuppi"),
        minPt   = cms.double(15.0),
        maxEta  = cms.double(4.0),
        # No usable default: Year is a knob on THREE modules (both JEC producers and
        # jvmAppliedEventFilter) and setting only some of them is silently wrong, not an
        # error. Each plugin throws on an empty Year; set it via
        #   from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra
        #   applyJmeEra(process, year="2022Pre", isData=False)
        Year    = cms.string(""),
        # False = keep every event and publish the decision as the "JvmPass" product, which
        # H_phipimu stores as the jet_vetomap_pass branch; the cut then happens offline.
        # True = the historical behaviour, drop the event here.
        #
        # This lives as a cfi default (unlike Year, which is deliberately forced per-pset):
        # Year's correct value differs per sample, so a default is silently wrong for some,
        # whereas Veto has one correct value for the whole analysis.
        #
        # Changed to False on 2026-08-17. Everything produced before that date was vetoed at
        # this step and has no jet_vetomap_pass branch, so old and new tuples must not be
        # mixed until every sample has been reproduced. It also makes nEventsGen_ a pre-veto
        # count, i.e. usable as a normalisation denominator, which it was not before.
        Veto    = cms.bool(False),
        JvmConfig = cms.FileInPath("EDAnalyzers/JMEcorr/data/JvmConfig.json"),
    )
    
)

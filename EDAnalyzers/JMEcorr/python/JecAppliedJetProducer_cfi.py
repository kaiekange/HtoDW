import FWCore.ParameterSet.Config as cms

jecAppliedJetProducer = cms.EDProducer(
    "JecAppliedJetProducer",
    isDebug = cms.bool(False),  # set True to dump [JERC DEBUG] logs

    Jets = cms.PSet(
        srcAK4 = cms.InputTag("slimmedJetsPuppi"),          # pat::JetCollection
        rho    = cms.InputTag("fixedGridRhoFastjetAll"),
        # No usable default: Year is a knob on THREE modules (both JEC producers and
        # jvmAppliedEventFilter) and setting only some of them is silently wrong, not an
        # error. Each plugin throws on an empty Year; set it via
        #   from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra
        #   applyJmeEra(process, year="2022Pre", isData=False)
        Year   = cms.string(""),
        IsData = cms.bool(False),

        # Optional era ("" => None)
        Era    = cms.string(""),

        # --- Choose one of: "Nominal", "JES", "JER"
        SystKind    = cms.string("Nominal"),

        # If SystKind == "JES"
        JesSystName = cms.string("AbsoluteStat"),  # correction set key
        JesSystVar  = cms.string("Up"),            # "Up" | "Down"

        # If SystKind == "JER"
        JerVar      = cms.string("nom"),           # "nom" | "up" | "down"
        JerRegion   = cms.PSet(                    # optional gate
            etaMin = cms.double(0.0),
            etaMax = cms.double(999.0),
            ptMin  = cms.double(0.0),
            ptMax  = cms.double(1.0e9),
        ),
        JecConfig = cms.FileInPath("EDAnalyzers/JMEcorr/data/JecConfigAK4.json"),
        JerToolConfig = cms.FileInPath("EDAnalyzers/JMEcorr/data/jer_smear.json.gz"),
    ),
)

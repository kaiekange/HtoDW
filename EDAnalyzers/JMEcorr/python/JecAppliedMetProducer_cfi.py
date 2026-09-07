import FWCore.ParameterSet.Config as cms

jecAppliedMetProducer = cms.EDProducer(
    "JecAppliedMetProducer",
    isDebug = cms.bool(False),

    Jets = cms.PSet(
        srcAK4 = cms.InputTag("slimmedJetsPuppi"),           # pat::JetCollection
        rho    = cms.InputTag("fixedGridRhoFastjetAll"),
        # No usable default: Year is a knob on THREE modules (both JEC producers and
        # jvmAppliedEventFilter) and setting only some of them is silently wrong, not an
        # error. Each plugin throws on an empty Year; set it via
        #   from EDAnalyzers.JMEcorr.jmeEra import applyJmeEra
        #   applyJmeEra(process, year="2022Pre", isData=False)
        Year   = cms.string(""),
        IsData = cms.bool(False),
        Era    = cms.string(""),

        # One of: "Nominal", "JES", "JER"
        SystKind    = cms.string("Nominal"),

        # If SystKind == "JES"
        JesSystName = cms.string("AbsoluteStat"),  # correction set key (empty => off)
        JesSystVar  = cms.string("Up"),            # "Up" | "Down"

        # If SystKind == "JER"
        JerVar      = cms.string("nom"),           # "nom" | "up" | "down"
        JerRegion   = cms.PSet(                    # optional region gate
            etaMin = cms.double(0.0),
            etaMax = cms.double(999.),
            ptMin  = cms.double(0.0),
            ptMax  = cms.double(1e9),
        ),

        # Optional: name of userFloat on jets storing muon subtraction fraction (0..1)
        MuonSubtrUserFloat = cms.string(""),       # e.g. "muonSubtrFactor"; empty => ignore
        JecConfig = cms.FileInPath("EDAnalyzers/JMEcorr/data/JecConfigAK4.json"),
        JerToolConfig = cms.FileInPath("EDAnalyzers/JMEcorr/data/jer_smear.json.gz"),
    ),

    # MET sources
    srcMET_PF    = cms.InputTag("slimmedMETs"),
    srcMET_PUPPI = cms.InputTag("slimmedMETsPuppi"),
)


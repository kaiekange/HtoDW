#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <memory>
#include <iostream>

#include "EDAnalyzers/JMEcorr/interface/JvmConfigReader.h"
#include "EDAnalyzers/JMEcorr/interface/JvmApplication.h"

class JvmAppliedEventFilter : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit JvmAppliedEventFilter(const edm::ParameterSet& iConfig);
  ~JvmAppliedEventFilter();

private:
  edm::ParameterSet                    configParamsJets_;
  edm::InputTag                        edmTagJetAK4_;
  std::string                          year_;
  bool                                 veto_;
  edm::EDGetTokenT<pat::JetCollection> edmTokenJetAK4_;
  edm::FileInPath                      jvmConfigPath_;
  JvmApplication::VetoChecker          jvmChecker_;

  bool filter(edm::Event&, const edm::EventSetup&) override;
};

JvmAppliedEventFilter::JvmAppliedEventFilter(const edm::ParameterSet& iConfig)
  : configParamsJets_(iConfig.getParameter<edm::ParameterSet>("Jets")),
    edmTagJetAK4_(configParamsJets_.getParameter<edm::InputTag>("SourcesAK4")),
    year_(configParamsJets_.getParameter<std::string>("Year")),
    veto_(configParamsJets_.getParameter<bool>("Veto")),
    jvmConfigPath_(configParamsJets_.getParameter<edm::FileInPath>("JvmConfig")),
    jvmChecker_() {

    // Year has no usable default in the cfi. This module is the one where a wrong Year is
    // fully silent -- the veto map is an (eta,phi) map with no run-number dimension, so a
    // preEE map applied to postEE data produces plausible numbers and no error anywhere.
    if (year_.empty()) {
        throw edm::Exception(edm::errors::Configuration)
            << "JvmAppliedEventFilter: 'Year' is empty. Call applyJmeEra(process, year, isData, era) "
               "from EDAnalyzers.JMEcorr.jmeEra in your pset -- it sets Year on all three JME modules.";
    }

    // Configure the JVM config singleton with the path from edm::FileInPath
    auto& jvmCfg = JvmConfigReader::JvmConfig::defaultInstance();
    jvmCfg.setConfigPath(jvmConfigPath_.fullPath(), /*reload=*/true);

    // Get the JVM for this year from the JSON config
    JvmConfigReader::Jvm jvm = jvmCfg.getJvmForYear(year_);
    if (jvm.use) {
        jvmChecker_.setJvm(jvm.ref, jvm.key);
    } else {
        // A year that is spelled correctly but absent from JvmConfig.json would otherwise
        // disable the veto map with only a line of stdout in a job log of many thousands.
        throw edm::Exception(edm::errors::Configuration)
            << "JvmAppliedEventFilter: no JVM entry for year \"" << year_
            << "\" in config \"" << jvmConfigPath_.fullPath()
            << "\". Refusing to run with the jet veto map silently disabled.";
    }

    if (!veto_) {
        std::cout << "[JvmAppliedEventFilter] Veto=False: no event is dropped. The decision "
                     "is published as the \"JvmPass\" product instead, for the analyzer to "
                     "store per event.\n";
    }

    // Declare consumes using the EDFilter base API
    edmTokenJetAK4_ = consumes<pat::JetCollection>(edmTagJetAK4_);

    // Published whatever Veto is set to, so the decision stays available to downstream
    // modules even when this filter is also dropping on it. Int rather than bool: a bare
    // bool event product needs an edm::Wrapper<bool> dictionary, which this package does
    // not declare.
    produces<int>("JvmPass");
}

JvmAppliedEventFilter::~JvmAppliedEventFilter() {}

bool JvmAppliedEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    bool passEvent = true;
    try {
        edm::Handle<pat::JetCollection> iJets;
        iEvent.getByToken(edmTokenJetAK4_, iJets);

        if (iJets.isValid()) {
            for (const auto& jIt : *iJets) {
                if (jvmChecker_.checkJetInVetoRegion(
                        jIt.eta(),
                        jIt.phi(),
                        jIt.pt(),
                        /* jetId: NanoAOD TightLepVeto not always available in MiniAOD */
                        6,
                        jIt.chargedEmEnergyFraction(),
                        jIt.neutralEmEnergyFraction())) {

                    passEvent = false;
                    break; // No need to loop over the remaining jets.
                }
            }
        }
    } catch (std::exception const& e) {
        std::cout << "[JvmAppliedEventFilter] exception in filter -> "
                  << e.what() << std::endl;
        // Be conservative on exception: do not veto the event. Falls through rather than
        // returning early, so the product below is always put -- an event missing it would
        // make the analyzer throw ProductNotFound.
        passEvent = true;
    }

    iEvent.put(std::make_unique<int>(passEvent ? 1 : 0), "JvmPass");

    // Veto=False keeps every event and leaves the decision to the analysis, which reads
    // the product above. Veto=True is the historical behaviour: the event is dropped here
    // and never reaches the analyzer.
    return veto_ ? passEvent : true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JvmAppliedEventFilter);


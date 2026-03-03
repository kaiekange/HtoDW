#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <memory>
#include <iostream>

#include "EDAnalyzers/Reconstruction/interface/JvmConfigReader.h"
#include "EDAnalyzers/Reconstruction/interface/JvmApplication.h"

class JvmAppliedEventFilter : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
    explicit JvmAppliedEventFilter(const edm::ParameterSet &iConfig);
    ~JvmAppliedEventFilter();

private:
    std::string year_ = "2022Pre";

    edm::EDGetTokenT<pat::JetCollection> srcJetsToken_;

    std::string jvmConfigPath_ = "EDAnalyzers/Reconstruction/json/JvmConfig.json";

    JvmApplication::VetoChecker jvmChecker_;

    Bool_t filter(edm::Event &, const edm::EventSetup &) override;
};

JvmAppliedEventFilter::JvmAppliedEventFilter(const edm::ParameterSet &iConfig) : jvmChecker_()
{
    // Configure the JVM config singleton with the path from edm::FileInPath
    auto &jvmCfg = JvmConfigReader::JvmConfig::defaultInstance();
    jvmCfg.setConfigPath(jvmConfigPath_, /*reload=*/true);

    // Get the JVM for this year from the JSON config
    JvmConfigReader::Jvm jvm = jvmCfg.getJvmForYear(year_);

    if (jvm.use)
        jvmChecker_.setJvm(jvm.ref, jvm.key);
    else
        std::cout << "[JvmAppliedEventFilter] No JVM entry for year \"" << year_ << "\" in config \"" << jvmConfigPath_ << "\". Filter will not veto any events.\n";

    // Declare consumes using the EDFilter base API
    edm::InputTag srcJetsTag_("slimmedJetsPuppi");
    srcJetsToken_ = consumes<pat::JetCollection>(srcJetsTag_);
}

JvmAppliedEventFilter::~JvmAppliedEventFilter() {}

Bool_t JvmAppliedEventFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{

    Bool_t passEvent = true;
    try
    {
        edm::Handle<pat::JetCollection> srcJetsHandle;
        iEvent.getByToken(srcJetsToken_, srcJetsHandle);

        if (srcJetsHandle.isValid())
        {
            for (const auto &jIt : *srcJetsHandle)
            {
                if (jvmChecker_.checkJetInVetoRegion(jIt.eta(), jIt.phi(), jIt.pt(), 6, jIt.chargedEmEnergyFraction(), jIt.neutralEmEnergyFraction()))
                {
                    passEvent = false;
                    break; // No need to loop over the remaining jets.
                }
            }
        }
    }
    catch (std::exception const &e)
    {
        std::cout << "[JvmAppliedEventFilter] exception in filter -> " << e.what() << std::endl;
        // Be conservative on exception: do not veto the event.
        return true;
    }
    return passEvent;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JvmAppliedEventFilter);

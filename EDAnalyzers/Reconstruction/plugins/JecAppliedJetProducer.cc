#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/Handle.h"

#include <memory>
#include <string>
#include <optional>
#include <vector>
#include <iostream>

#include "EDAnalyzers/Reconstruction/interface/JecConfigReader.h"
#include "EDAnalyzers/Reconstruction/interface/JecApplication.h"

class JecAppliedJetProducer : public edm::stream::EDProducer<>
{
public:
	explicit JecAppliedJetProducer(const edm::ParameterSet &iConfig);
	~JecAppliedJetProducer() override = default;

	void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

private:
	Bool_t isDebug_{false};
	Bool_t isData_{false};
	std::string year_ = "2022Pre";
	std::string era_ = "";

	// --- Tokens
	edm::EDGetTokenT<pat::JetCollection> srcJetsToken_;
	edm::EDGetTokenT<Double_t> rhoToken_;

	// --- Config / correction backends
	edm::FileInPath jecConfigPath_ = "EDAnalyzers/Reconstruction/json/JecConfigAK4.json";
	edm::FileInPath jerToolConfigPath_ = "EDAnalyzers/Reconstruction/json/jer_smear.json.gz";

	// Pointer to the global JEC config singleton
	JecConfigReader::JecConfig *jecCfg_{nullptr};

	// JES syst correction ref (MC only, cached once)
	std::optional<correction::Correction::Ref> jesSystRef_;

	// Helpers
	static Double_t computeRawFactorFromMiniAOD(const pat::Jet &j)
	{
		const Double_t ptCorr = j.pt();
		const Double_t ptRaw = j.correctedP4("Uncorrected").pt(); // MiniAOD stored
		if (ptCorr <= 0.0)
			return 0.0;
		const Double_t rf = 1.0 - (ptRaw / ptCorr);
		// clamp to sane range
		return std::min(std::max(rf, 0.0), 1.0);
	}

	static inline void setP4Scaled(pat::Jet &j, Double_t scale)
	{
		// Scale full p4 (JEC/JER are p4 scalings)
		const auto p4 = j.p4();
		const auto p4Scaled = reco::Particle::LorentzVector(
			p4.px() * scale, p4.py() * scale, p4.pz() * scale, p4.energy() * scale);
		j.setP4(p4Scaled);
	}

	static inline void addUserFloats(pat::Jet &j,
									 Double_t jes,
									 Double_t jesSyst,
									 Double_t jer,
									 Double_t total)
	{
		j.addUserFloat("JERC:jesFactor", jes);
		j.addUserFloat("JERC:jesSystFactor", jesSyst);
		j.addUserFloat("JERC:jerFactor", jer);
		j.addUserFloat("JERC:totalFactor", total);
	}
};

JecAppliedJetProducer::JecAppliedJetProducer(const edm::ParameterSet &iConfig)
{
	edm::InputTag srcJetsTag_("slimmedJetsPuppi");
	srcJetsToken_ = consumes<pat::JetCollection>(srcJetsTag_);

	edm::InputTag rhoTag_("fixedGridRhoFastjetAll");
	rhoToken_ = consumes<Double_t>(rhoTag_);

	// Era (needed for data path)
	if (jetsPset_.existsAs<std::string>("Era"))
		era_ = jetsPset_.getParameter<std::string>("Era");
	if (isData_ && era_.empty())
		throw edm::Exception(edm::errors::Configuration) << "JecAppliedJetProducer: For data (IsData=true) you must provide a non-empty 'Era' in the Jets PSet.";

	// Use the resolved filesystem path from edm::FileInPath and configure
	// the singleton defaults *before* first construction.
	JecConfigReader::ConfigPaths paths;
	paths.ak4 = jecConfigPath_.fullPath();

	JecConfigReader::JecConfig::setDefaultPaths(paths);
	JecConfigReader::JecConfig::setDefaultJerSmearPath(jerToolConfigPath_.fullPath());

	jecCfg_ = &JecConfigReader::JecConfig::defaultInstance();

	produces<pat::JetCollection>("CorrectedAK4");
}

void JecAppliedJetProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{

	// Inputs
	edm::Handle<pat::JetCollection> srcJetsHandle;
	iEvent.getByToken(srcJetsToken_, srcJetsHandle);

	edm::Handle<Double_t> rhoHandle;
	iEvent.getByToken(rhoToken_, rhoHandle);
	const Double_t rho = rhoHandle.isValid() ? *rhoHandle : 0.0;

	auto out = std::make_unique<pat::JetCollection>();
	out->reserve(srcJetsHandle->size());

	// Build systematics payload
	JecApplication::SystematicOptions systOpts{};

	// Per-event Applier (handles JEC/JER handles and run-based residual logic)
	JecApplication::Applier applier = [&]() -> JecApplication::Applier
	{
		if (isData_)
		{
			std::optional<Double_t> runNumber;
			if (JecApplication::requiresRunBasedResidual(year_))
				runNumber = static_cast<Double_t>(iEvent.id().run());
			return JecApplication::Applier::DataAK4(*jecCfg_, year_, era_, isDebug_, runNumber);
		}
		else
			return JecApplication::Applier::McAK4(*jecCfg_, year_, isDebug_);
	}();

	// Loop jets
	for (const auto &jet : *srcJetsHandle)
	{
		pat::Jet outJet = jet;

		const Double_t rawFactor = computeRawFactorFromMiniAOD(jet);

		JecApplication::JesInputs jesIn{
			jet.pt(),
			jet.eta(),
			jet.phi(),
			jet.jetArea(),
			rho,
			rawFactor};

		// --- JES nominal factor (via Applier) ---
		const Double_t jesFactor = applier.jesFactorNominal(jesIn);

		// Apply JES nominal
		setP4Scaled(outJet, jesFactor);

		// --- JER (MC only) ---
		double jerFactor = 1.0;
		if (!isData_)
		{
			// State after JES (and optional JES syst)
			JecApplication::JesInputs jAfter{
				outJet.pt(),
				outJet.eta(),
				outJet.phi(),
				jet.jetArea(),
				rho,
				0.0};

			// Gen-matching if available
			JecApplication::JerInputs jerIn{};
			jerIn.event = iEvent.id().event();
			jerIn.rho = rho;

			if (const auto *gj = jet.genJet())
			{
				jerIn.hasGen = true;
				jerIn.genPt = gj->pt();
				jerIn.genEta = gj->eta();
				jerIn.genPhi = gj->phi();
				jerIn.maxDr = 0.2; // AK4
			}

			jerFactor = applier.jerFactor(jAfter, jerIn, systOpts);
			setP4Scaled(outJet, jerFactor);
		}

		// Attach userFloats for downstream bookkeeping
		const Double_t total = jesFactor * jesSystFactor * jerFactor;
		addUserFloats(outJet, jesFactor, jesSystFactor, jerFactor, total);

		out->push_back(std::move(outJet));
	}

	iEvent.put(std::move(out), "CorrectedAK4");
}

// Register
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JecAppliedJetProducer);

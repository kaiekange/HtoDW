#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "EDAnalyzers/Reconstruction/interface/JecConfigReader.h"
#include "EDAnalyzers/Reconstruction/interface/JecApplication.h"

class JecAppliedMetProducer : public edm::stream::EDProducer<>
{
public:
	explicit JecAppliedMetProducer(const edm::ParameterSet &iConfig);
	~JecAppliedMetProducer() override = default;

	void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

private:
	Bool_t isDebug_{false};
	Bool_t isData_{false};
	std::string year_ = "2022Pre";
	std::string era_ = "";

	// Optional: read muon subtraction factor from jet userFloat (0..1). Empty => disabled.
	std::string muonSubtrUserFloat_ = "";

	// ---- Tokens ----
	edm::EDGetTokenT<pat::JetCollection> srcJetsToken_;
	edm::EDGetTokenT<Double_t> rhoToken_;
	edm::EDGetTokenT<pat::METCollection> metPFToken_;
	edm::EDGetTokenT<pat::METCollection> metPuppiToken_;

	// ---- Config / corrections (AK4) ----
	std::string jecConfigPath_ = "EDAnalyzers/Reconstruction/json/JecConfigAK4.json";
	std::string jerToolConfigPath_ = "EDAnalyzers/Reconstruction/json/jer_smear.json.gz";
	JecConfigReader::JecConfig *jecCfg_{nullptr};

	// JES syst correction ref (MC only, cached once)
	std::optional<correction::Correction::Ref> jesSystRef_;

	// Helpers
	static inline Double_t getRawPtFromMiniAOD(const pat::Jet &j)
	{
		return j.correctedP4("Uncorrected").pt(); // MiniAOD stores this
	}

	static inline Double_t getMuonSubtrFraction(const pat::Jet &j, const std::string &name)
	{
		if (name.empty())
			return 0.0;
		if (j.hasUserFloat(name))
		{
			const Double_t f = j.userFloat(name);
			if (std::isfinite(f))
				return std::max(0.0, std::min(1.0, f));
		}
		return 0.0;
	}

	static inline Double_t getMetRawPt(const pat::MET &met)
	{
		// Prefer uncorrected if available, else fall back with a warning.
		Double_t v = met.uncorPt();
		if (!(v >= 0.0))
		{
			std::cout << "[JecAppliedMetProducer] WARNING: met.uncorPt() invalid; using met.pt().\n";
			v = met.pt();
		}
		return v;
	}

	static inline Double_t getMetRawPhi(const pat::MET &met)
	{
		Double_t v = met.uncorPhi();
		if (!std::isfinite(v))
		{
			std::cout << "[JecAppliedMetProducer] WARNING: met.uncorPhi() invalid; using met.phi().\n";
			v = met.phi();
		}
		return v;
	}

	static inline pat::MET makeCorrectedMetFrom(const pat::MET &src, const TLorentzVector &p4corr)
	{
		pat::MET out = src; // copy PAT object (keeps extras)
		// Replace p4 with corrected pT/phi, zero eta/mass (MET mass=0, eta unused)
		const auto corr = reco::Particle::LorentzVector(p4corr.Px(), p4corr.Py(), 0.0, p4corr.Pt() /* E≈|pT| for massless */);
		out.setP4(corr);
		// Bookkeeping helpers (raw pt/phi will be set by caller)
		out.addUserFloat("JERC:MET_corr_pt", p4corr.Pt());
		out.addUserFloat("JERC:MET_corr_phi", p4corr.Phi());
		return out;
	}
};

JecAppliedMetProducer::JecAppliedMetProducer(const edm::ParameterSet &iConfig)
{
	edm::InputTag srcJetsTag_("slimmedJetsPuppi");
	srcJetsToken_ = consumes<pat::JetCollection>(srcJetsTag_);

	edm::InputTag rhoTag_("fixedGridRhoFastjetAll");
	rhoToken_ = consumes<Double_t>(rhoTag_);

	edm::InputTag srcMetPFTag_("slimmedMETs");
	metPFToken_ = consumes<pat::METCollection>(srcMetPFTag_);

	edm::InputTag srcMetPuppiTag_("slimmedMETsPuppi");
	metPuppiToken_ = consumes<pat::METCollection>(srcMetPuppiTag_);

	// Era (needed for data path)
	if (isData_ && era_.empty())
		throw edm::Exception(edm::errors::Configuration) << "JecAppliedMetProducer: For data (IsData=true) you must provide a non-empty 'Era' in the Jets PSet.";

	JecConfigReader::ConfigPaths paths;
	paths.ak4 = jecConfigPath_;

	JecConfigReader::JecConfig::setDefaultPaths(paths);
	JecConfigReader::JecConfig::setDefaultJerSmearPath(jerToolConfigPath_);

	// Get the singleton instance
	jecCfg_ = &JecConfigReader::JecConfig::defaultInstance();

	produces<pat::METCollection>("CorrectedMet");
}

void JecAppliedMetProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	// Inputs
	edm::Handle<pat::JetCollection> srcJetsHandle;
	iEvent.getByToken(srcJetsToken_, srcJetsHandle);

	edm::Handle<Double_t> rhoHandle;
	iEvent.getByToken(rhoToken_, rhoHandle);
	const Double_t rho = rhoHandle.isValid() ? *rhoHandle : 0.0;

	const Bool_t usePuppi = JecApplication::usesPuppiMet(year_);

	edm::Handle<pat::METCollection> metHandle;
	if (usePuppi)
		iEvent.getByToken(metPuppiToken_, metHandle);
	else
		iEvent.getByToken(metPFToken_, metHandle);

	if (!metHandle.isValid() || metHandle->empty())
		throw edm::Exception(edm::errors::ProductNotFound) << "JecAppliedMetProducer: MET collection is invalid or empty.";

	const pat::MET &metIn = metHandle->front();

	const Double_t metRawPt = getMetRawPt(metIn);
	const Double_t metRawPhi = getMetRawPhi(metIn);

	// Build jetsForMet (RAW AK4 BEFORE any modifications)
	std::vector<JecApplication::JetForMet> jetsForMet;
	std::vector<JecApplication::JerInputs> jersForMet;
	jetsForMet.reserve(srcJetsHandle->size());
	jersForMet.reserve(srcJetsHandle->size());

	for (const auto &j : *srcJetsHandle)
	{
		JecApplication::JetForMet v;
		v.phi = j.phi();
		v.eta = j.eta();
		v.area = j.jetArea();
		v.rawPt = getRawPtFromMiniAOD(j);
		v.muonSubtrFactor = getMuonSubtrFraction(j, muonSubtrUserFloat_);
		v.chEmEf = j.chargedEmEnergyFraction();
		v.neEmEf = j.neutralEmEnergyFraction();
		jetsForMet.push_back(v);

		JecApplication::JerInputs jerIn;
		jerIn.event = iEvent.id().event();
		jerIn.rho = rho;
		if (const auto *gj = j.genJet())
		{
			jerIn.hasGen = true;
			jerIn.genPt = gj->pt();
			jerIn.genEta = gj->eta();
			jerIn.genPhi = gj->phi();
			jerIn.maxDr = 0.2; // AK4
		}
		jersForMet.push_back(jerIn);
	}

	// Per-event Applier (same pattern as jet producer)
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

	// Systematics payload
	JecApplication::SystematicOptions systOpts{};

	// Call the shared Type-1 MET routine through Applier
	JecApplication::MetInputs metInputs{metRawPt, metRawPhi};
	const TLorentzVector p4CorrectedMET = applier.correctedMet(metInputs, jetsForMet, jersForMet, rho, systOpts);

	// Build and put output (copy original PAT MET; replace p4)
	auto out = std::make_unique<pat::METCollection>();
	pat::MET metOut = makeCorrectedMetFrom(metIn, p4CorrectedMET);
	// record which raw was used
	metOut.addUserFloat("JERC:MET_raw_pt_used", metRawPt);
	metOut.addUserFloat("JERC:MET_raw_phi_used", metRawPhi);
	out->push_back(std::move(metOut));

	iEvent.put(std::move(out), "CorrectedMet");
}

// Register module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JecAppliedMetProducer);

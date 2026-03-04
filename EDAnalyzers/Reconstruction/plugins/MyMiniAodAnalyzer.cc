// system
#include <memory>
#include <string>
#include <vector>
#include <cstdint>

// ROOT
#include "TTree.h"

// CMSSW framework (one::EDAnalyzer with shared TFileService)
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

// DataFormats (MiniAOD)
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class MyMiniAodAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
	explicit MyMiniAodAnalyzer(const edm::ParameterSet &);
	void analyze(const edm::Event &, const edm::EventSetup &) override;

private:
	void bookTree();

	// tokens
	edm::EDGetTokenT<std::vector<pat::Jet>> jetAK4Tok_;
	edm::EDGetTokenT<std::vector<pat::Jet>> jetAK8Tok_;
	edm::EDGetTokenT<std::vector<pat::MET>> metTok_;
	bool readAK4_{true}, readAK8_{true};

	// output (stable addresses for TTree)
	edm::Service<TFileService> tfs_;
	TTree *tree_{nullptr};

	// event-level
	UInt_t run_{0};
	ULong64_t event_{0};
	UInt_t luminosityBlock_{0};

	// jets/met
	std::vector<float> JetAK4_pt_;
	std::vector<float> JetAK8_pt_;
	float MET_pt_{0.f};
};

MyMiniAodAnalyzer::MyMiniAodAnalyzer(const edm::ParameterSet &ps)
{
	usesResource("TFileService");
	const auto jetAK4Src = ps.getParameter<edm::InputTag>("jetAK4Src");
	const auto jetAK8Src = ps.getParameter<edm::InputTag>("jetAK8Src");
	const auto metSrc = ps.getParameter<edm::InputTag>("metSrc");

	readAK4_ = !jetAK4Src.label().empty();
	readAK8_ = !jetAK8Src.label().empty();
	if (readAK4_)
		jetAK4Tok_ = consumes<std::vector<pat::Jet>>(jetAK4Src);
	if (readAK8_)
		jetAK8Tok_ = consumes<std::vector<pat::Jet>>(jetAK8Src);
	metTok_ = consumes<std::vector<pat::MET>>(metSrc);

	tree_ = tfs_->make<TTree>("Events", "Events");
	tree_->SetAutoSave(10000000000LL);
	tree_->SetAutoFlush(1000000);
	bookTree();
}

void MyMiniAodAnalyzer::bookTree()
{
	tree_->Branch("run", &run_, "run/i");
	tree_->Branch("luminosityBlock", &luminosityBlock_, "luminosityBlock/i");
	tree_->Branch("event", &event_, "event/l");

	tree_->Branch("JetAK4_pt", &JetAK4_pt_);
	tree_->Branch("JetAK8_pt", &JetAK8_pt_);
	tree_->Branch("MET_pt", &MET_pt_, "MET_pt/F");
}

void MyMiniAodAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &)
{
	// clear per event
	JetAK4_pt_.clear();
	JetAK8_pt_.clear();
	MET_pt_ = 0.f;

	// ids
	run_ = iEvent.id().run();
	luminosityBlock_ = iEvent.id().luminosityBlock();
	event_ = static_cast<ULong64_t>(iEvent.id().event());

	// AK4
	if (readAK4_)
	{
		edm::Handle<std::vector<pat::Jet>> h;
		iEvent.getByToken(jetAK4Tok_, h);
		if (h.isValid())
		{
			JetAK4_pt_.reserve(h->size());
			for (const auto &j : *h)
				JetAK4_pt_.push_back(j.pt());
		}
	}

	// AK8
	if (readAK8_)
	{
		edm::Handle<std::vector<pat::Jet>> h;
		iEvent.getByToken(jetAK8Tok_, h);
		if (h.isValid())
		{
			JetAK8_pt_.reserve(h->size());
			for (const auto &j : *h)
				JetAK8_pt_.push_back(j.pt());
		}
	}

	// MET (first entry)
	edm::Handle<std::vector<pat::MET>> hmet;
	iEvent.getByToken(metTok_, hmet);
	if (hmet.isValid() && !hmet->empty())
		MET_pt_ = hmet->front().pt();

	tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MyMiniAodAnalyzer);

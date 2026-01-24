#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TTree.h"
#include "TLorentzVector.h"

class MatchAnalyzer : public edm::one::EDAnalyzer<>
{
public:
	explicit MatchAnalyzer(const edm::ParameterSet &);
	~MatchAnalyzer() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	bool hasAncestor(const reco::GenParticle &gp, const int motherid, const int otherparticleid) const;
	void endJob() override;

	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;

	unsigned long nEventsProcessed_;
	unsigned long nEventsValidKp_;
	unsigned long nEventsValidKm_;
	unsigned long nEventsValidphi_;
	unsigned long nEventsValidpi_;
	unsigned long nEventsValidDs_;
	unsigned long nEventsValidmu_;

	const edm::Service<TFileService> fs;
	TTree * InfoTree;
};

MatchAnalyzer::MatchAnalyzer(const edm::ParameterSet &iConfig)
{
	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);
}

void MatchAnalyzer::beginJob()
{
	nEventsProcessed_ = 0;
	nEventsValidKp_ = 0;
	nEventsValidKm_ = 0;
	nEventsValidphi_ = 0;
	nEventsValidpi_ = 0;
	nEventsValidDs_ = 0;
	nEventsValidmu_ = 0;
	InfoTree = fs->make<TTree>("RunInfo", "RunInfo");
    InfoTree->Branch("nEventsProcessed_", &nEventsProcessed_);
    InfoTree->Branch("nEventsValidKp_", &nEventsValidKp_);
    InfoTree->Branch("nEventsValidKm_", &nEventsValidKm_);
	InfoTree->Branch("nEventsValidphi_", &nEventsValidphi_);
	InfoTree->Branch("nEventsValidpi_", &nEventsValidpi_);
	InfoTree->Branch("nEventsValidDs_", &nEventsValidDs_);
	InfoTree->Branch("nEventsValidmu_", &nEventsValidmu_);
}

void MatchAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
	iEvent.getByToken(associationToken_, associationHandle);

	nEventsProcessed_++;

	size_t idx_pi = -1;

	Int_t num_Kp = 0;
	Int_t num_Km = 0;
	Int_t num_pi = 0;
	Int_t num_mu = 0;

	for (size_t i = 0; i < PFCandHandle->size(); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (abs(pf.pdgId()) != 211)
			continue;
		if (!pf.hasTrackDetails())
			continue;

		pat::PackedCandidateRef packedPFCandRef(PFCandHandle, i);
		reco::GenParticleRef genRef = (*associationHandle)[packedPFCandRef];

		if (genRef.isNull())
			continue;

		const reco::GenParticle &gp = *genRef;

		if (abs(gp.pdgId()) == 211)
		{
			if (!hasAncestor(gp, 431, 333))
				continue;
			if (!hasAncestor(gp, 25, 24))
				continue;

			num_pi++;
			idx_pi = i;
		}
		else if (abs(gp.pdgId()) == 321)
		{
			if (!hasAncestor(gp, 333, 321))
				continue;
			if (!hasAncestor(gp, 431, 211))
				continue;
			if (!hasAncestor(gp, 25, 24))
				continue;

			if (pf.charge() == 1)
			{
				num_Kp++;
			}
			else if (pf.charge() == -1)
			{
				num_Km++;
			}
		}
	}

	if (num_Kp == 1)
		nEventsValidKp_++;
	if (num_Km == 1)
		nEventsValidKm_++;
	if (num_pi == 1)
		nEventsValidpi_++;


	if (num_Kp != 1)
		return;
	if (num_Km != 1)
		return;
	nEventsValidphi_++;
	if (num_pi != 1)
		return;
	nEventsValidDs_++;

	const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];

	for (size_t i = 0; i < MuonHandle->size(); i++)
	{
		const pat::Muon &muon = (*MuonHandle)[i];
		if (muon.charge() == pi_pf.charge())
			continue;

		if (!muon.genParticle())
			continue;

		const reco::GenParticle &gp = *(muon.genParticle());

		if (abs(gp.pdgId()) != 13)
			continue;
		if (!hasAncestor(gp, 24, 14))
			continue;
		if (!hasAncestor(gp, 25, 431))
			continue;

		num_mu++;
	}

	if (num_mu != 1)
		return;
	nEventsValidmu_++;
}

bool MatchAnalyzer::hasAncestor(const reco::GenParticle &gp, const int motherid, const int otherparticleid) const
{
	bool hasancestor = false;

	if (gp.numberOfMothers() == 0)
		return false;

	if (gp.numberOfMothers() == 1 && abs(gp.mother(0)->pdgId()) == motherid)
	{
		if (gp.mother(0)->numberOfDaughters() == 2)
		{
			if (abs(gp.mother(0)->daughter(0)->pdgId()) == otherparticleid || abs(gp.mother(0)->daughter(1)->pdgId()) == otherparticleid)
				hasancestor = true;
		}
	}

	if (hasancestor)
		return true;
	else
		return hasAncestor(*gp.motherRef(0), motherid, otherparticleid);
}

void MatchAnalyzer::endJob()
{
	std::cout << nEventsProcessed_ << ", Kp: " << nEventsValidKp_ << ", Km " << nEventsValidKm_ << ", pi" << nEventsValidpi_ << ", mu" << nEventsValidmu_ << std::endl;
	InfoTree->Fill();
}

void MatchAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MatchAnalyzer);

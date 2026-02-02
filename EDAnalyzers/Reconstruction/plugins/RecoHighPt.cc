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
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "TLorentzVector.h"
#include "TTree.h"

class RecoHighPt : public edm::one::EDAnalyzer<>
{
public:
	explicit RecoHighPt(const edm::ParameterSet &);
	~RecoHighPt() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;

	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

	const Double_t Mass_Ds = 1.96835;
	const Double_t Mass_phi = 1.019460;
	const Double_t Mass_K = 0.493677;
	const Double_t Mass_pi = 0.13957039;
	const Double_t Mass_mu = 0.1056583755;
	unsigned long nEventsProcessed_;
	unsigned long nEventsKpMatch_;
	unsigned long nEventsKpNonMatch_;
	unsigned long nEventsKmMatch_;
	unsigned long nEventsKmNonMatch_;
	unsigned long nEventspiMatch_;
	unsigned long nEventspiNonMatch_;
	unsigned long nEventsphiMatch_;
	unsigned long nEventsphiNonMatch_;
	unsigned long nEventsDsMatch_;
	unsigned long nEventsDsNonMatch_;
	unsigned long nEventsMuMatch_;
	unsigned long nEventsMuNonMatch_;
	unsigned long nEventsMatch_;
	unsigned long nEventsNonMatch_;

	KalmanVertexFitter KVFitter;
};

RecoHighPt::RecoHighPt(const edm::ParameterSet &iConfig)
	: KVFitter(true)
{
	edm::InputTag TriggerBitsTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);
}

void RecoHighPt::beginJob()
{
	nEventsProcessed_ = 0;
	nEventsKpMatch_ = 0;
	nEventsKpNonMatch_ = 0;
	nEventsKmMatch_ = 0;
	nEventsKmNonMatch_ = 0;
	nEventspiMatch_ = 0;
	nEventspiNonMatch_ = 0;
	nEventsphiMatch_ = 0;
	nEventsphiNonMatch_ = 0;
	nEventsDsMatch_ = 0;
	nEventsDsNonMatch_ = 0;
	nEventsMuMatch_ = 0;
	nEventsMuNonMatch_ = 0;
	nEventsMatch_ = 0;
	nEventsNonMatch_ = 0;
}

void RecoHighPt::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);

	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
	iEvent.getByToken(associationToken_, associationHandle);

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	// nEventsProcessed_++;

	Int_t Kp_match_idx = -777;
	Int_t Km_match_idx = -777;
	Int_t pi_match_idx = -777;
	Int_t mu_match_idx = -777;

	Float_t Kp_deltaR = 1.0;
	Float_t Km_deltaR = 1.0;
	Float_t pi_deltaR = 1.0;
	Float_t mu_deltaR = 1.0;

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		const pat::PackedCandidateRef pfRef(PFCandHandle, i);
		const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
		if (genRef.isNull())
			continue;
		const reco::GenParticle &gp = *genRef;
		Float_t deltaR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());
		if (deltaR > 0.05)
			continue;
		Float_t deltapt = abs((gp.pt() - pf.pt()) / gp.pt());
		if (deltapt > 0.3)
			continue;

		if ((gp.pdgId() == 321))
		{
			if (hasAncestor(gp, 333, 321, -321) && hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < Kp_deltaR)
				{
					Kp_match_idx = i;
					Kp_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == -321))
		{
			if (hasAncestor(gp, 333, 321, -321) && hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < Km_deltaR)
				{
					Km_match_idx = i;
					Km_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == 211))
		{
			if (hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < pi_deltaR)
				{
					pi_match_idx = i;
					pi_deltaR = deltaR;
				}
			}
		}
	}

	for (Int_t i = 0; i < static_cast<Int_t>(MuonHandle->size()); i++)
	{
		const pat::Muon &muon = (*MuonHandle)[i];
		if (!muon.genParticle())
			continue;
		const reco::GenParticle &gp = *(muon.genParticle());
		Float_t deltaR = reco::deltaR(gp.eta(), gp.phi(), muon.eta(), muon.phi());
		if (deltaR > 0.05)
			continue;
		Float_t deltapt = abs((gp.pt() - muon.pt()) / gp.pt());
		if (deltapt > 0.3)
			continue;

		if ((gp.pdgId() == 13) && hasAncestor(gp, -24, 13, -14) && hasAncestor(gp, 25, 431, -24))
		{
			if (deltaR < mu_deltaR)
			{
				mu_match_idx = i;
				mu_deltaR = deltaR;
			}
		}
	}

	std::vector<Int_t> idx_pf;

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (!pf.hasTrackDetails())
			continue;
		if (!pf.trackHighPurity())
			continue;
		if (pf.pt() < 0.5)
			continue;
		if (pf.p() < 1)
			continue;
		if (abs(pf.pdgId()) != 211)
			continue;

		idx_pf.push_back(i);
	}

	if (idx_pf.size() < 3)
		return;

	Int_t best_Kp_idx = -1;
	Int_t best_Km_idx = -1;
	Int_t best_pi_idx = -1;
	Float_t best_Ds_pt = -1.0;

	Int_t Ds_pass = 0;

	for (Int_t iKp = 0; iKp < static_cast<Int_t>(idx_pf.size()); iKp++)
	{
		const Int_t idx_Kp = idx_pf[iKp];

		const pat::PackedCandidate &cnd_Kp_pf = (*PFCandHandle)[idx_Kp];
		if (cnd_Kp_pf.charge() != 1)
			continue;

		for (Int_t iKm = 0; iKm < static_cast<Int_t>(idx_pf.size()); iKm++)
		{
			if (iKp == iKm)
				continue;
			const Int_t idx_Km = idx_pf[iKm];

			const pat::PackedCandidate &cnd_Km_pf = (*PFCandHandle)[idx_Km];
			if (cnd_Km_pf.charge() != -1)
				continue;

			if (reco::deltaR(cnd_Kp_pf.eta(), cnd_Kp_pf.phi(), cnd_Km_pf.eta(), cnd_Km_pf.phi()) > 0.1)
				continue;
			if (abs(cnd_Kp_pf.pt() - cnd_Km_pf.pt()) > 20)
				continue;

			std::vector<reco::TransientTrack> cnd_phi_ttrk = {
				(*ttBuilder).build(cnd_Kp_pf.pseudoTrack()),
				(*ttBuilder).build(cnd_Km_pf.pseudoTrack())};

			TransientVertex cnd_phi_tvtx = KVFitter.vertex(cnd_phi_ttrk);
			if (!cnd_phi_tvtx.isValid())
				continue;
			if (!cnd_phi_tvtx.hasRefittedTracks())
				continue;
			if (cnd_phi_tvtx.normalisedChiSquared() < 0)
				continue;
			if (cnd_phi_tvtx.normalisedChiSquared() > 10)
				continue;

			reco::Vertex cnd_phi_vtx = reco::Vertex(cnd_phi_tvtx);
			std::vector<reco::TransientTrack> cnd_phiFit_ttrk = cnd_phi_tvtx.refittedTracks();
			TLorentzVector cnd_phiFit_Kp_tlv, cnd_phiFit_Km_tlv, cnd_phiFit_phi_tlv;
			cnd_phiFit_Kp_tlv.SetXYZM(cnd_phiFit_ttrk[0].track().px(), cnd_phiFit_ttrk[0].track().py(), cnd_phiFit_ttrk[0].track().pz(), Mass_K);
			cnd_phiFit_Km_tlv.SetXYZM(cnd_phiFit_ttrk[1].track().px(), cnd_phiFit_ttrk[1].track().py(), cnd_phiFit_ttrk[1].track().pz(), Mass_K);
			cnd_phiFit_phi_tlv = cnd_phiFit_Kp_tlv + cnd_phiFit_Km_tlv;
			if (reco::deltaR(cnd_phiFit_Kp_tlv.Eta(), cnd_phiFit_Kp_tlv.Phi(), cnd_phiFit_phi_tlv.Eta(), cnd_phiFit_phi_tlv.Phi()) > 0.05)
				continue;
			if (reco::deltaR(cnd_phiFit_Km_tlv.Eta(), cnd_phiFit_Km_tlv.Phi(), cnd_phiFit_phi_tlv.Eta(), cnd_phiFit_phi_tlv.Phi()) > 0.05)
				continue;
			if ((cnd_phiFit_phi_tlv.M() - Mass_phi) > 0.1)
				continue;
			if (cnd_phiFit_phi_tlv.Pt() < 2)
				continue;

			for (Int_t ipi = 0; ipi < static_cast<Int_t>(idx_pf.size()); ipi++)
			{
				if ((iKp == ipi) || (iKm == ipi))
					continue;
				const Int_t idx_pi = idx_pf[ipi];
				const pat::PackedCandidate &cnd_pi_pf = (*PFCandHandle)[idx_pi];

				if (reco::deltaR(cnd_Kp_pf.eta(), cnd_Kp_pf.phi(), cnd_pi_pf.eta(), cnd_pi_pf.phi()) > 0.4)
					continue;
				if (reco::deltaR(cnd_Km_pf.eta(), cnd_Km_pf.phi(), cnd_pi_pf.eta(), cnd_pi_pf.phi()) > 0.4)
					continue;
				if (reco::deltaR(cnd_phiFit_phi_tlv.Eta(), cnd_phiFit_phi_tlv.Phi(), cnd_pi_pf.eta(), cnd_pi_pf.phi()) > 0.4)
					continue;

				std::vector<reco::TransientTrack> cnd_Ds_ttrk = {
					cnd_phiFit_ttrk[0],
					cnd_phiFit_ttrk[1],
					(*ttBuilder).build(cnd_pi_pf.pseudoTrack())};

				TransientVertex cnd_Ds_tvtx = KVFitter.vertex(cnd_Ds_ttrk);
				if (!cnd_Ds_tvtx.isValid())
					continue;
				if (!cnd_Ds_tvtx.hasRefittedTracks())
					continue;
				if (cnd_Ds_tvtx.normalisedChiSquared() < 0)
					continue;
				if (cnd_Ds_tvtx.normalisedChiSquared() > 10)
					continue;
				reco::Vertex cnd_Ds_vtx = reco::Vertex(cnd_Ds_tvtx);

				std::vector<reco::TransientTrack> cnd_DsFit_ttrk = cnd_Ds_tvtx.refittedTracks();
				TLorentzVector cnd_Kp_tlv, cnd_Km_tlv, cnd_pi_tlv, cnd_phi_tlv, cnd_Ds_tlv;
				cnd_Kp_tlv.SetXYZM(cnd_DsFit_ttrk[0].track().px(), cnd_DsFit_ttrk[0].track().py(), cnd_DsFit_ttrk[0].track().pz(), Mass_K);
				cnd_Km_tlv.SetXYZM(cnd_DsFit_ttrk[1].track().px(), cnd_DsFit_ttrk[1].track().py(), cnd_DsFit_ttrk[1].track().pz(), Mass_K);
				cnd_pi_tlv.SetXYZM(cnd_DsFit_ttrk[2].track().px(), cnd_DsFit_ttrk[2].track().py(), cnd_DsFit_ttrk[2].track().pz(), Mass_pi);
				cnd_phi_tlv = cnd_Kp_tlv + cnd_Km_tlv;
				cnd_Ds_tlv = cnd_Kp_tlv + cnd_Km_tlv + cnd_pi_tlv;

				if (reco::deltaR(cnd_Kp_tlv.Eta(), cnd_Kp_tlv.Phi(), cnd_Ds_tlv.Eta(), cnd_Ds_tlv.Phi()) > 0.15)
					continue;
				if (reco::deltaR(cnd_Km_tlv.Eta(), cnd_Km_tlv.Phi(), cnd_Ds_tlv.Eta(), cnd_Ds_tlv.Phi()) > 0.15)
					continue;
				if (reco::deltaR(cnd_pi_tlv.Eta(), cnd_pi_tlv.Phi(), cnd_Ds_tlv.Eta(), cnd_Ds_tlv.Phi()) > 0.4)
					continue;
				if (reco::deltaR(cnd_phi_tlv.Eta(), cnd_phi_tlv.Phi(), cnd_Ds_tlv.Eta(), cnd_Ds_tlv.Phi()) > 0.15)
					continue;

				if ((cnd_phi_tlv.M() - Mass_phi) > 0.03)
					continue;
				if (cnd_phi_tlv.Pt() < 3)
					continue;
				if ((cnd_Ds_tlv.M() - Mass_Ds) > 0.15)
					continue;
				if (cnd_Ds_tlv.Pt() < 7)
					continue;

				GlobalPoint Ds_pos(cnd_Ds_vtx.x(), cnd_Ds_vtx.y(), cnd_Ds_vtx.z());
				GlobalVector displacement = Ds_pos - primvtx_pos;
				GlobalVector PDirec(cnd_Ds_tlv.Px(), cnd_Ds_tlv.Py(), cnd_Ds_tlv.Pz());
				VertexDistance3D vtxdist3D;
				Measurement1D Ds_FD_3D = vtxdist3D.signedDistance(primvtx, cnd_Ds_vtx, PDirec);
				if (abs(Ds_FD_3D.value()) > 4)
					continue;
				Float_t direction_angle = displacement.dot(PDirec) / (displacement.mag() * PDirec.mag());
				if(direction_angle < 0)
					continue;

				if ((idx_Kp == Kp_match_idx) && (idx_Km == Km_match_idx) && (idx_pi == pi_match_idx))
					Ds_pass++;

				if (cnd_Ds_tlv.Pt() > best_Ds_pt)
				{
					best_Kp_idx = idx_Kp;
					best_Km_idx = idx_Km;
					best_pi_idx = idx_pi;
					best_Ds_pt = cnd_Ds_tlv.Pt();
				}
			}
		}
	}

	if (best_Kp_idx == -1)
		return;
	if (best_Km_idx == -1)
		return;
	if (best_pi_idx == -1)
		return;

	if (Ds_pass > 0 && !((best_Kp_idx == Kp_match_idx) && (best_Km_idx == Km_match_idx) && (best_pi_idx == pi_match_idx)))
		std::cout << "no pass" << std::endl;

	const pat::PackedCandidate &best_pi_pf = (*PFCandHandle)[best_pi_idx];

	Int_t best_mu_idx = -1;
	Float_t best_mu_pt = -1;
	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &cnd_muon = (*MuonHandle)[idx_mu];
		if (cnd_muon.charge() == best_pi_pf.charge())
			continue;
		if (cnd_muon.pt() < 7)
			continue;
		// if (cnd_muon.p() < 1)
		// 	continue;
		if (!cnd_muon.isTightMuon(primvtx))
			continue;

		if (cnd_muon.pt() > best_mu_pt)
		{
			best_mu_idx = idx_mu;
			best_mu_pt = cnd_muon.pt();
		}
	}

	if (best_mu_idx == -1)
		return;

	Bool_t Kp_match = (best_Kp_idx == Kp_match_idx);
	Bool_t Km_match = (best_Km_idx == Km_match_idx);
	Bool_t pi_match = (best_pi_idx == pi_match_idx);
	Bool_t mu_match = (best_mu_idx == mu_match_idx);

	if (Kp_match)
		nEventsKpMatch_++;
	else
		nEventsKpNonMatch_++;
	if (Km_match)
		nEventsKmMatch_++;
	else
		nEventsKmNonMatch_++;
	if (pi_match)
		nEventspiMatch_++;
	else
		nEventspiNonMatch_++;
	if (Kp_match && Km_match)
		nEventsphiMatch_++;
	else
		nEventsphiNonMatch_++;
	if (Kp_match && Km_match && pi_match)
		nEventsDsMatch_++;
	else
		nEventsDsNonMatch_++;
	if (mu_match)
		nEventsMuMatch_++;
	else
		nEventsMuNonMatch_++;
	if (Kp_match && Km_match && pi_match && mu_match)
		nEventsMatch_++;
	else
		nEventsNonMatch_++;
}

Bool_t RecoHighPt::hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
{
	if (gp.numberOfMothers() == 0)
		return false;

	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother)
			continue;

		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 2)
		{
			Bool_t daughter1_flag = false;
			Bool_t daughter2_flag = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id)
					daughter1_flag = true;
				if (mother->daughter(j)->pdgId() == daughter2id)
					daughter2_flag = true;
			}

			if (daughter1_flag && daughter2_flag)
				return true;
		}

		if (hasAncestor(*(gp.motherRef(i)), motherid, daughter1id, daughter2id))
			return true;
	}

	return false;
}

void RecoHighPt::endJob()
{
	std::cout << "Total: " << nEventsProcessed_ << std::endl;
	std::cout << "Kp match: " << nEventsKpMatch_ << ", nonmatch: " << nEventsKpNonMatch_ << std::endl;
	std::cout << "Km match: " << nEventsKmMatch_ << ", nonmatch: " << nEventsKmNonMatch_ << std::endl;
	std::cout << "pi match: " << nEventspiMatch_ << ", nonmatch: " << nEventspiNonMatch_ << std::endl;
	std::cout << "phi match: " << nEventsphiMatch_ << ", nonmatch: " << nEventsphiNonMatch_ << std::endl;
	std::cout << "Ds match: " << nEventsDsMatch_ << ", nonmatch: " << nEventsDsNonMatch_ << std::endl;
	std::cout << "Mu match: " << nEventsMuMatch_ << ", nonmatch: " << nEventsMuNonMatch_ << std::endl;
	std::cout << "H match: " << nEventsMatch_ << ", nonmatch: " << nEventsNonMatch_ << std::endl;
}

void RecoHighPt::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoHighPt);

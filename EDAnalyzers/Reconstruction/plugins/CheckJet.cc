#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <string>
#include <optional>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
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

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "EDAnalyzers/Reconstruction/interface/FullRecoTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class CheckJet : public edm::one::EDAnalyzer<>
{
public:
	explicit CheckJet(const edm::ParameterSet &);
	~CheckJet() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<pat::METCollection> MetToken_;
	edm::EDGetTokenT<pat::JetCollection> JetToken_;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

	KinematicParticleFactoryFromTransientTrack factory;
	KinematicParticleVertexFitter KPVFitter;
	KinematicConstrainedVertexFitter KCVFitter;

	const Float_t Mass_Ds = 1.96835;
	const Float_t Mass_phi = 1.019460;
	const Float_t Mass_K = 0.493677;
	const Float_t Mass_pi = 0.13957039;
	const Float_t Mass_mu = 0.1056583755;

	Float_t chi = 0;
	Float_t ndf = 0;
	Float_t Mass_K_sigma = 0.000013;
	Float_t Mass_pi_sigma = 0.00000017;

	TTree *atree;
	FullRecoTree *ftree;

	const edm::Service<TFileService> fs;
	TTree *outtree;

	std::vector<Float_t> jet_pt;
	Int_t num_jets_pt30;
};

CheckJet::CheckJet(const edm::ParameterSet &iConfig)
{
	edm::InputTag TriggerBitsTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	// edm::InputTag MetTag_("slimmedMETsPuppi");
	// MetToken_ = consumes<pat::METCollection>(MetTag_);

	// edm::InputTag JetTag_("slimmedJetsPuppi");
	// JetToken_ = consumes<pat::JetCollection>(JetTag_);

	edm::InputTag MetTag_("jecAppliedMetProducer", "CorrectedMet");
	MetToken_ = consumes<pat::METCollection>(MetTag_);

	edm::InputTag JetTag_("jecAppliedJetProducer", "CorrectedAK4");
	JetToken_ = consumes<pat::JetCollection>(JetTag_);

	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);
}

void CheckJet::beginJob()
{
	ftree = new FullRecoTree(atree);

	outtree = fs->make<TTree>("Jets", "Jets");
	outtree->Branch("jet_pt", &jet_pt);
	outtree->Branch("num_jets_pt30", &num_jets_pt30);
}

void CheckJet::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	edm::Handle<pat::METCollection> MetHandle;
	iEvent.getByToken(MetToken_, MetHandle);

	edm::Handle<pat::JetCollection> JetHandle;
	iEvent.getByToken(JetToken_, JetHandle);

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	ftree->Init();

	//==================================================//
	//--------------------------------------------------//
	// primary vertex
	//--------------------------------------------------//
	//==================================================//
	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_cov(primvtx.covariance());

	//==================================================//
	//--------------------------------------------------//
	// Triggers
	//--------------------------------------------------//
	//==================================================//
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);
	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		Bool_t pass = triggerHandle->accept(i);
		if (trigName.Contains("HLT_IsoMu24_v"))
			ftree->trig_HLT_IsoMu24 = pass;
		else if (trigName.Contains("HLT_Mu50_v"))
			ftree->trig_HLT_Mu50 = pass;
	}
	if (!(ftree->trig_HLT_IsoMu24) && !(ftree->trig_HLT_Mu50))
		return;

	//==================================================//
	//--------------------------------------------------//
	// Ds candidate
	//--------------------------------------------------//
	//==================================================//

	// prepare hadron candidates
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
	RefCountedKinematicTree best_DsTree;

	//--------------------------------------------------//
	// loop for K+
	//--------------------------------------------------//
	for (Int_t iKp = 0; iKp < static_cast<Int_t>(idx_pf.size()); iKp++)
	{
		const Int_t idx_Kp = idx_pf[iKp];

		const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[idx_Kp];
		// select on charge and pt
		if (Kp_pf.charge() != 1)
			continue;
		if (Kp_pf.pt() < 1)
			continue;

		TLorentzVector Kp_tlv;
		Kp_tlv.SetXYZM(Kp_pf.px(), Kp_pf.py(), Kp_pf.pz(), Mass_K);

		reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.bestTrack());

		//--------------------------------------------------//
		// loop for K-
		//--------------------------------------------------//
		for (Int_t iKm = 0; iKm < static_cast<Int_t>(idx_pf.size()); iKm++)
		{
			if (iKp == iKm)
				continue;
			const Int_t idx_Km = idx_pf[iKm];

			const pat::PackedCandidate &Km_pf = (*PFCandHandle)[idx_Km];
			// select on charge and pt
			if (Km_pf.charge() != -1)
				continue;
			if (Km_pf.pt() < 1)
				continue;

			// preliminary cut on phi
			TLorentzVector Km_tlv;
			Km_tlv.SetXYZM(Km_pf.px(), Km_pf.py(), Km_pf.pz(), Mass_K);
			if ((Kp_tlv + Km_tlv).Pt() < 2.5)
				continue;
			if (abs((Kp_tlv + Km_tlv).M() - Mass_phi) > 0.1)
				continue;

			reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.bestTrack());

			//--------------------------------------------------//
			// loop for pi
			//--------------------------------------------------//
			for (Int_t ipi = 0; ipi < static_cast<Int_t>(idx_pf.size()); ipi++)
			{
				if ((iKp == ipi) || (iKm == ipi))
					continue;
				const Int_t idx_pi = idx_pf[ipi];
				const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];

				// preliminary cut on Ds
				TLorentzVector pi_tlv;
				pi_tlv.SetXYZM(pi_pf.px(), pi_pf.py(), pi_pf.pz(), Mass_pi);
				if ((Kp_tlv + Km_tlv + pi_tlv).Pt() < 5)
					continue;
				if (abs((Kp_tlv + Km_tlv + pi_tlv).M() - Mass_Ds) > 0.3)
					continue;

				reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.bestTrack());

				//--------------------------------------------------//
				// Kinematic Fit on Ds
				//--------------------------------------------------//
				std::vector<RefCountedKinematicParticle> DsDaughters;
				DsDaughters.push_back(factory.particle(Kp_tt, Mass_K, chi, ndf, Mass_K_sigma));
				DsDaughters.push_back(factory.particle(Km_tt, Mass_K, chi, ndf, Mass_K_sigma));
				DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi, ndf, Mass_pi_sigma));
				RefCountedKinematicTree DsTree = KPVFitter.fit(DsDaughters);

				if (!DsTree->isValid())
					continue;
				DsTree->movePointerToTheTop();

				RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
				if (!Ds_vertex->vertexIsValid())
					continue;
				GlobalPoint Ds_vtx_pos = Ds_vertex->position();

				// select fit chi2 probability
				if (TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom()) < 0.01)
					continue;

				RefCountedKinematicParticle Ds = DsTree->currentParticle();
				GlobalVector Ds_mtm = Ds->currentState().globalMomentum();
				DsTree->movePointerToTheFirstChild();
				RefCountedKinematicParticle Kp = DsTree->currentParticle();
				GlobalVector Kp_mtm = Kp->currentState().globalMomentum();
				DsTree->movePointerToTheNextChild();
				RefCountedKinematicParticle Km = DsTree->currentParticle();
				GlobalVector Km_mtm = Km->currentState().globalMomentum();

				TLorentzVector final_Kp_tlv, final_Km_tlv;
				final_Kp_tlv.SetXYZM(Kp_mtm.x(), Kp_mtm.y(), Kp_mtm.z(), Mass_K);
				final_Km_tlv.SetXYZM(Km_mtm.x(), Km_mtm.y(), Km_mtm.z(), Mass_K);

				// select on phi and Ds candidates
				if (abs((final_Kp_tlv + final_Km_tlv).M() - Mass_phi) > 0.03)
					continue;
				if ((final_Kp_tlv + final_Km_tlv).Pt() < 3)
					continue;
				if (abs(Ds->currentState().mass() - Mass_Ds) > 0.15)
					continue;
				if (Ds_mtm.perp() < 7)
					continue;

				// select dira
				GlobalVector displacement = Ds_vtx_pos - primvtx_pos;
				if (displacement.dot(Ds_mtm) / (displacement.mag() * Ds_mtm.mag()) < 0.9)
					continue;

				//--------------------------------------------------//
				// best Ds candidate: highest pt
				//--------------------------------------------------//
				if (Ds_mtm.perp() > best_Ds_pt)
				{
					best_Kp_idx = idx_Kp;
					best_Km_idx = idx_Km;
					best_pi_idx = idx_pi;
					best_DsTree = DsTree;
					best_Ds_pt = Ds_mtm.perp();
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

	const pat::PackedCandidate &best_pi_pf = (*PFCandHandle)[best_pi_idx];

	//==================================================//
	//--------------------------------------------------//
	// Muon candidate
	//--------------------------------------------------//
	//==================================================//
	Int_t best_mu_idx = -1;
	Float_t best_mu_pt = -1;

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (muon.charge() == best_pi_pf.charge())
			continue;
		if (muon.pt() < 24)
			continue;
		if (!muon.isTightMuon(primvtx))
			continue;

		//--------------------------------------------------//
		// best mu candidate: highest pt
		//--------------------------------------------------//
		if (muon.pt() > best_mu_pt)
		{
			best_mu_idx = idx_mu;
			best_mu_pt = muon.pt();
		}
	}

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (!muon.isTightMuon(primvtx))
			continue;

		// if (muon.charge() == best_pi_pf.charge())
		// 	ftree->num_ss_mu_tight++;
		// else if (muon.charge() == -best_pi_pf.charge())
		// 	ftree->num_os_mu_tight++;
	}

	if (best_mu_idx == -1)
		return;

	// //==================================================//
	// //--------------------------------------------------//
	// // MET
	// //--------------------------------------------------//
	// //==================================================//
	// const pat::MET &Met = MetHandle->front();

	jet_pt.clear();
	// Float_t deltaR_Ds_Jet = 100.;
	num_jets_pt30 = 0;
	for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
	{
		const pat::Jet &Jet = (*JetHandle)[idx_jet];
		jet_pt.push_back(Jet.pt());
		
		if (Jet.pt() < 30)
			continue;
		num_jets_pt30++;
		// Float_t tmp_deltaR_Ds_Jet = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Jet.eta(), Jet.phi());

		// if (tmp_deltaR_Ds_Jet < deltaR_Ds_Jet)
		// 	deltaR_Ds_Jet = tmp_deltaR_Ds_Jet;
	}

	outtree->Fill();
}

Measurement1D CheckJet::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x();
	vDiff[1] = diff.y();
	vDiff[2] = 0.;

	Double_t dist = sqrt(pow(diff.x(), 2) + pow(diff.y(), 2));

	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = 0;
	if (dist != 0)
		err = sqrt(err2) / dist;

	return Measurement1D(dist, err);
}

Measurement1D CheckJet::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x();
	vDiff[1] = diff.y();
	vDiff[2] = diff.z();

	Double_t dist = diff.mag();

	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = 0.;
	if (dist != 0)
		err = sqrt(err2) / dist;

	return Measurement1D(dist, err);
}

void CheckJet::endJob()
{
}

void CheckJet::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CheckJet);

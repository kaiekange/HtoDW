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
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

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

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "EDAnalyzers/GenAnalyzer/interface/ThreeFitters_Tree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class ThreeFitters : public edm::one::EDAnalyzer<>
{
public:
	explicit ThreeFitters(const edm::ParameterSet &);
	~ThreeFitters() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const;
	Bool_t TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;
	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;

	KinematicParticleFactoryFromTransientTrack factory;
	KinematicParticleVertexFitter KPVFitter;

	KalmanVertexFitter KLMFitter{true};
	AdaptiveVertexFitter AVFFitter;

	const Float_t Mass_Ds = 1.96835;
	const Float_t Mass_phi = 1.019460;
	const Float_t Mass_K = 0.493677;
	const Float_t Mass_pi = 0.13957039;
	const Float_t Mass_mu = 0.1056583755;

	Float_t chi = 0;
	Float_t ndf = 0;
	Float_t Mass_K_sigma = 0.000013;
	Float_t Mass_pi_sigma = 0.00000017;

	const edm::Service<TFileService> fs;
	ThreeFitters_Tree *ftree;
};

ThreeFitters::ThreeFitters(const edm::ParameterSet &iConfig)
{
	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);
}

void ThreeFitters::beginJob()
{
	ftree = new ThreeFitters_Tree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
}

void ThreeFitters::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	ftree->Init();

	Int_t Gen_H_idx = 0;
	Int_t Gen_Ds_idx = 0;
	Int_t Gen_W_idx = 0;
	Int_t Gen_mu_idx = 0;
	Int_t Gen_Kp_idx = 0;
	Int_t Gen_Km_idx = 0;
	Int_t Gen_pi_idx = 0;

	Int_t num_Gen_H = 0;
	Int_t num_Gen_Ds = 0;
	Int_t num_Gen_W = 0;
	Int_t num_Gen_mu = 0;
	Int_t num_Gen_e = 0;
	Int_t num_Gen_Kp = 0;
	Int_t num_Gen_Km = 0;
	Int_t num_Gen_pi = 0;

	//==================================================//
	//--------------------------------------------------//
	// Gen
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];
		if (!gp.isHardProcess())
			continue;
		if (gp.pdgId() == 25)
		{
			Gen_H_idx = i;
			num_Gen_H++;
		}
		else if ((abs(gp.pdgId()) == 431))
		{
			Gen_Ds_idx = i;
			num_Gen_Ds++;
		}
		else if (abs(gp.pdgId()) == 24)
		{
			Gen_W_idx = i;
			num_Gen_W++;
		}
		else if (abs(gp.pdgId()) == 13)
		{
			Gen_mu_idx = i;
			num_Gen_mu++;
		}
		else if (abs(gp.pdgId()) == 11)
		{
			num_Gen_e++;
		}
	}
	if (!((num_Gen_H == 1) && (num_Gen_Ds == 1) && (num_Gen_W == 1) && (num_Gen_mu == 1) && (num_Gen_e == 0)))
		return;
	const reco::GenParticle &H_gp = (*GenPartHandle)[Gen_H_idx];
	const reco::GenParticle &Ds_gp = (*GenPartHandle)[Gen_Ds_idx];
	const reco::GenParticle &W_gp = (*GenPartHandle)[Gen_W_idx];
	const reco::GenParticle &mu_gp = (*GenPartHandle)[Gen_mu_idx];

	const Int_t Ds_charge = Ds_gp.charge();

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];
		if (gp.status() != 1)
			continue;
		if (gp.pdgId() == (Ds_charge * 321))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), 321, -321, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				Gen_Kp_idx = i;
				num_Gen_Kp++;
			}
		}
		else if (gp.pdgId() == (-Ds_charge * 321))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), 321, -321, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				Gen_Km_idx = i;
				num_Gen_Km++;
			}
		}
		else if (gp.pdgId() == (Ds_charge * 211))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), 321, -321, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				Gen_pi_idx = i;
				num_Gen_pi++;
			}
		}
	}
	if (!((num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1)))
		return;
	const reco::GenParticle &Kp_gp = (*GenPartHandle)[Gen_Kp_idx];
	const reco::GenParticle &Km_gp = (*GenPartHandle)[Gen_Km_idx];
	const reco::GenParticle &pi_gp = (*GenPartHandle)[Gen_pi_idx];

	ftree->Gen_PV_vx = H_gp.vx();
	ftree->Gen_PV_vy = H_gp.vy();
	ftree->Gen_PV_vz = H_gp.vz();

	ftree->Gen_SV_vx = Ds_gp.vx();
	ftree->Gen_SV_vy = Ds_gp.vy();
	ftree->Gen_SV_vz = Ds_gp.vz();

	ftree->Gen_Ds_eta = Ds_gp.eta();
	ftree->Gen_Ds_phi = Ds_gp.phi();
	ftree->Gen_Ds_px = Ds_gp.px();
	ftree->Gen_Ds_py = Ds_gp.py();
	ftree->Gen_Ds_pz = Ds_gp.pz();
	ftree->Gen_Ds_pt = Ds_gp.pt();
	ftree->Gen_Ds_p = Ds_gp.p();
	ftree->Gen_Ds_mass = Ds_gp.mass();

	//==================================================//
	//--------------------------------------------------//
	// Gen Match
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
	iEvent.getByToken(associationToken_, associationHandle);

	Int_t Kp_match_idx = -1;
	Int_t Km_match_idx = -1;
	Int_t pi_match_idx = -1;

	Float_t Kp_deltaR = 1.0;
	Float_t Km_deltaR = 1.0;
	Float_t pi_deltaR = 1.0;

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

		if ((gp.pdgId() == Kp_gp.pdgId()))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), Kp_gp.pdgId(), Km_gp.pdgId(), pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				if (deltaR < Kp_deltaR)
				{
					Kp_match_idx = i;
					Kp_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == Km_gp.pdgId()))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), Kp_gp.pdgId(), Km_gp.pdgId(), pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				if (deltaR < Km_deltaR)
				{
					Km_match_idx = i;
					Km_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == pi_gp.pdgId()))
		{
			if (ThreeBodyDecay(gp, Ds_gp.pdgId(), Kp_gp.pdgId(), Km_gp.pdgId(), pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{
				if (deltaR < pi_deltaR)
				{
					pi_match_idx = i;
					pi_deltaR = deltaR;
				}
			}
		}
	}
	if ((Kp_match_idx < 0) || (Km_match_idx < 0) || (pi_match_idx < 0))
		return;
	const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[Kp_match_idx];
	const pat::PackedCandidate &Km_pf = (*PFCandHandle)[Km_match_idx];
	const pat::PackedCandidate &pi_pf = (*PFCandHandle)[pi_match_idx];

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	Int_t mu_match_idx = -1;
	Float_t mu_deltaR = 1.0;

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
		if (gp.pdgId() != mu_gp.pdgId())
			continue;

		if (TwoBodyDecay(gp, W_gp.pdgId(), mu_gp.pdgId(), mu_gp.charge() * 14) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
		{
			if (deltaR < mu_deltaR)
			{
				mu_match_idx = i;
				mu_deltaR = deltaR;
			}
		}
	}
	if (mu_match_idx < 0)
		return;
	const pat::Muon &mu = (*MuonHandle)[mu_match_idx];

	//==================================================//
	//--------------------------------------------------//
	// primary vertex
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);

	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_cov(primvtx.covariance());
	ftree->primvtx_isValid = primvtx.isValid();
	ftree->primvtx_isFake = primvtx.isFake();
	ftree->primvtx_chi2 = primvtx.chi2();
	ftree->primvtx_ndof = primvtx.ndof();
	ftree->primvtx_chi2ndof = primvtx.chi2() / primvtx.ndof();
	ftree->primvtx_vx = primvtx.x();
	ftree->primvtx_vy = primvtx.y();
	ftree->primvtx_vz = primvtx.z();
	ftree->primvtx_vxerr = primvtx.xError();
	ftree->primvtx_vyerr = primvtx.yError();
	ftree->primvtx_vzerr = primvtx.zError();

	//==================================================//
	//--------------------------------------------------//
	// Reconstruction
	//--------------------------------------------------//
	//==================================================//
	reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.pseudoTrack());
	reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.pseudoTrack());
	reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.pseudoTrack());

	//--------------------------------------------------//
	// Kinematic Vertex Fit on Ds
	//--------------------------------------------------//
	std::vector<RefCountedKinematicParticle> KVF_DsDaughters;
	KVF_DsDaughters.push_back(factory.particle(Kp_tt, Mass_K, chi, ndf, Mass_K_sigma));
	KVF_DsDaughters.push_back(factory.particle(Km_tt, Mass_K, chi, ndf, Mass_K_sigma));
	KVF_DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi, ndf, Mass_pi_sigma));
	RefCountedKinematicTree KVF_DsTree = KPVFitter.fit(KVF_DsDaughters);
	if (!KVF_DsTree->isValid())
		return;
	KVF_DsTree->movePointerToTheTop();
	RefCountedKinematicVertex KVF_Ds_vertex = KVF_DsTree->currentDecayVertex();
	if (!KVF_Ds_vertex->vertexIsValid())
		return;
	GlobalPoint KVF_Ds_vtx_pos = KVF_Ds_vertex->position();
	GlobalError KVF_Ds_vtx_cov = KVF_Ds_vertex->error();

	// KVF_DsTree particles: Ds
	RefCountedKinematicParticle KVF_Ds = KVF_DsTree->currentParticle();
	GlobalVector KVF_Ds_mtm = KVF_Ds->currentState().globalMomentum();
	// KVF_DsTree particles: K+
	KVF_DsTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle KVF_Kp = KVF_DsTree->currentParticle();
	GlobalVector KVF_Kp_mtm = KVF_Kp->currentState().globalMomentum();
	// KVF_DsTree particles: K-
	KVF_DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle KVF_Km = KVF_DsTree->currentParticle();
	GlobalVector KVF_Km_mtm = KVF_Km->currentState().globalMomentum();
	// KVF_DsTree particles: pi
	KVF_DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle KVF_pi = KVF_DsTree->currentParticle();
	GlobalVector KVF_pi_mtm = KVF_pi->currentState().globalMomentum();

	TLorentzVector KVF_Kp_tlv, KVF_Km_tlv, KVF_pi_tlv, KVF_Ds_tlv, KVF_KK_tlv, KVF_Kppi_tlv, KVF_Kmpi_tlv;
	KVF_Kp_tlv.SetXYZM(KVF_Kp_mtm.x(), KVF_Kp_mtm.y(), KVF_Kp_mtm.z(), Mass_K);
	KVF_Km_tlv.SetXYZM(KVF_Km_mtm.x(), KVF_Km_mtm.y(), KVF_Km_mtm.z(), Mass_K);
	KVF_pi_tlv.SetXYZM(KVF_pi_mtm.x(), KVF_pi_mtm.y(), KVF_pi_mtm.z(), Mass_pi);
	KVF_KK_tlv = KVF_Kp_tlv + KVF_Km_tlv;
	KVF_Kppi_tlv = KVF_Kp_tlv + KVF_pi_tlv;
	KVF_Kmpi_tlv = KVF_Km_tlv + KVF_pi_tlv;

	GlobalVector KVF_displacement = KVF_Ds_vtx_pos - primvtx_pos;
	Measurement1D KVF_Ds_mtmDs_FD_XY = distanceXY(primvtx_pos, primvtx_cov, KVF_Ds_vtx_pos, KVF_Ds_vtx_cov);
	Measurement1D KVF_Ds_mtmDs_FD_3D = distance3D(primvtx_pos, primvtx_cov, KVF_Ds_vtx_pos, KVF_Ds_vtx_cov);

	ftree->KVF_KK_mass = KVF_KK_tlv.M();
	ftree->KVF_Kppi_mass = KVF_Kppi_tlv.M();
	ftree->KVF_Kmpi_mass = KVF_Kmpi_tlv.M();

	ftree->KVF_Ds_eta = KVF_Ds_mtm.eta();
	ftree->KVF_Ds_phi = KVF_Ds_mtm.phi();
	ftree->KVF_Ds_px = KVF_Ds_mtm.x();
	ftree->KVF_Ds_py = KVF_Ds_mtm.y();
	ftree->KVF_Ds_pz = KVF_Ds_mtm.z();
	ftree->KVF_Ds_pt = KVF_Ds_mtm.perp();
	ftree->KVF_Ds_p = KVF_Ds_mtm.mag();
	ftree->KVF_Ds_mass = KVF_Ds->currentState().mass();

	ftree->KVF_Ds_vx = KVF_Ds_vtx_pos.x();
	ftree->KVF_Ds_vy = KVF_Ds_vtx_pos.y();
	ftree->KVF_Ds_vz = KVF_Ds_vtx_pos.z();
	ftree->KVF_Ds_vxerr = KVF_Ds_vtx_cov.cxx();
	ftree->KVF_Ds_vyerr = KVF_Ds_vtx_cov.cyy();
	ftree->KVF_Ds_vzerr = KVF_Ds_vtx_cov.czz();
	ftree->KVF_Ds_chi2 = KVF_Ds_vertex->chiSquared();
	ftree->KVF_Ds_ndf = KVF_Ds_vertex->degreesOfFreedom();
	ftree->KVF_Ds_chi2ndf = KVF_Ds_vertex->chiSquared() / KVF_Ds_vertex->degreesOfFreedom();
	ftree->KVF_Ds_prob = TMath::Prob(KVF_Ds_vertex->chiSquared(), KVF_Ds_vertex->degreesOfFreedom());

	ftree->KVF_dira = KVF_displacement.dot(KVF_Ds_mtm) / (KVF_displacement.mag() * KVF_Ds_mtm.mag());
	ftree->KVF_Ds_FDxy = KVF_Ds_mtmDs_FD_XY.value();
	ftree->KVF_Ds_FDxyerr = KVF_Ds_mtmDs_FD_XY.error();
	ftree->KVF_Ds_FDxysig = KVF_Ds_mtmDs_FD_XY.significance();
	ftree->KVF_Ds_FD = KVF_Ds_mtmDs_FD_3D.value();
	ftree->KVF_Ds_FDerr = KVF_Ds_mtmDs_FD_3D.error();
	ftree->KVF_Ds_FDsig = KVF_Ds_mtmDs_FD_3D.significance();

	//--------------------------------------------------//
	// Kalman Vertex Fit on Ds
	//--------------------------------------------------//
	std::vector<reco::TransientTrack> KLM_DsDaughters = {
		(*ttBuilder).build(Kp_pf.pseudoTrack()),
		(*ttBuilder).build(Km_pf.pseudoTrack()),
		(*ttBuilder).build(pi_pf.pseudoTrack())
	};
	TransientVertex KLM_Ds_vertex = KLMFitter.vertex(KLM_DsDaughters);

	if (!KLM_Ds_vertex.isValid())
		return;
	if (!KLM_Ds_vertex.hasRefittedTracks())
		return;
	GlobalPoint KLM_Ds_vtx_pos = KLM_Ds_vertex.position();
	GlobalError KLM_Ds_vtx_cov = KLM_Ds_vertex.positionError();

	std::vector<reco::TransientTrack> KLM_RFTracks = KLM_Ds_vertex.refittedTracks();

	TLorentzVector KLM_Kp_tlv, KLM_Km_tlv, KLM_pi_tlv, KLM_Ds_tlv, KLM_KK_tlv, KLM_Kppi_tlv, KLM_Kmpi_tlv;
	KLM_Kp_tlv.SetXYZM(KLM_RFTracks[0].track().px(), KLM_RFTracks[0].track().py(), KLM_RFTracks[0].track().pz(), Mass_K);
	KLM_Km_tlv.SetXYZM(KLM_RFTracks[1].track().px(), KLM_RFTracks[1].track().py(), KLM_RFTracks[1].track().pz(), Mass_K);
	KLM_pi_tlv.SetXYZM(KLM_RFTracks[2].track().px(), KLM_RFTracks[2].track().py(), KLM_RFTracks[2].track().pz(), Mass_pi);
	KLM_KK_tlv = KLM_Kp_tlv + KLM_Km_tlv;
	KLM_Kppi_tlv = KLM_Kp_tlv + KLM_pi_tlv;
	KLM_Kmpi_tlv = KLM_Km_tlv + KLM_pi_tlv;
	KLM_Ds_tlv = KLM_Kp_tlv + KLM_Km_tlv + KLM_pi_tlv;

	GlobalVector KLM_displacement = KLM_Ds_vtx_pos - primvtx_pos;
	GlobalVector KLM_Ds_mtm(KLM_Ds_tlv.Px(), KLM_Ds_tlv.Py(), KLM_Ds_tlv.Pz());
	Measurement1D KLM_Ds_mtmDs_FD_XY = distanceXY(primvtx_pos, primvtx_cov, KLM_Ds_vtx_pos, KLM_Ds_vtx_cov);
	Measurement1D KLM_Ds_mtmDs_FD_3D = distance3D(primvtx_pos, primvtx_cov, KLM_Ds_vtx_pos, KLM_Ds_vtx_cov);

	ftree->KLM_KK_mass = KLM_KK_tlv.M();
	ftree->KLM_Kppi_mass = KLM_Kppi_tlv.M();
	ftree->KLM_Kmpi_mass = KLM_Kmpi_tlv.M();

	ftree->KLM_Ds_eta = KLM_Ds_tlv.Eta();
	ftree->KLM_Ds_phi = KLM_Ds_tlv.Phi();
	ftree->KLM_Ds_px = KLM_Ds_tlv.Px();
	ftree->KLM_Ds_py = KLM_Ds_tlv.Py();
	ftree->KLM_Ds_pz = KLM_Ds_tlv.Pz();
	ftree->KLM_Ds_pt = KLM_Ds_tlv.Pt();
	ftree->KLM_Ds_p = KLM_Ds_tlv.P();
	ftree->KLM_Ds_mass = KLM_Ds_tlv.M();;

	ftree->KLM_Ds_vx = KLM_Ds_vtx_pos.x();
	ftree->KLM_Ds_vy = KLM_Ds_vtx_pos.y();
	ftree->KLM_Ds_vz = KLM_Ds_vtx_pos.z();
	ftree->KLM_Ds_vxerr = KLM_Ds_vtx_cov.cxx();
	ftree->KLM_Ds_vyerr = KLM_Ds_vtx_cov.cyy();
	ftree->KLM_Ds_vzerr = KLM_Ds_vtx_cov.czz();
	ftree->KLM_Ds_chi2 = KLM_Ds_vertex.totalChiSquared();
	ftree->KLM_Ds_ndf = KLM_Ds_vertex.degreesOfFreedom();
	ftree->KLM_Ds_chi2ndf = KLM_Ds_vertex.totalChiSquared() / KLM_Ds_vertex.degreesOfFreedom();
	ftree->KLM_Ds_prob = TMath::Prob(KLM_Ds_vertex.totalChiSquared(), KLM_Ds_vertex.degreesOfFreedom());

	ftree->KLM_dira = KLM_displacement.dot(KLM_Ds_mtm) / (KLM_displacement.mag() * KLM_Ds_mtm.mag());
	ftree->KLM_Ds_FDxy = KLM_Ds_mtmDs_FD_XY.value();
	ftree->KLM_Ds_FDxyerr = KLM_Ds_mtmDs_FD_XY.error();
	ftree->KLM_Ds_FDxysig = KLM_Ds_mtmDs_FD_XY.significance();
	ftree->KLM_Ds_FD = KLM_Ds_mtmDs_FD_3D.value();
	ftree->KLM_Ds_FDerr = KLM_Ds_mtmDs_FD_3D.error();
	ftree->KLM_Ds_FDsig = KLM_Ds_mtmDs_FD_3D.significance();

	//--------------------------------------------------//
	// Adaptive Vertex Fitter on Ds
	//--------------------------------------------------//
	std::vector<reco::TransientTrack> AVF_DsDaughters = {
		(*ttBuilder).build(Kp_pf.pseudoTrack()),
		(*ttBuilder).build(Km_pf.pseudoTrack()),
		(*ttBuilder).build(pi_pf.pseudoTrack())
	};
	TransientVertex AVF_Ds_vertex = AVFFitter.vertex(AVF_DsDaughters);

	if (!AVF_Ds_vertex.isValid())
		return;   
	GlobalPoint AVF_Ds_vtx_pos = AVF_Ds_vertex.position();
	GlobalError AVF_Ds_vtx_cov = AVF_Ds_vertex.positionError();

	TLorentzVector AVF_Kp_tlv, AVF_Km_tlv, AVF_pi_tlv, AVF_Ds_tlv, AVF_KK_tlv, AVF_Kppi_tlv, AVF_Kmpi_tlv;
	AVF_Kp_tlv.SetXYZM(Kp_pf.px(), Kp_pf.py(), Kp_pf.pz(), Mass_K);
	AVF_Km_tlv.SetXYZM(Km_pf.px(), Km_pf.py(), Km_pf.pz(), Mass_K);
	AVF_pi_tlv.SetXYZM(pi_pf.px(), pi_pf.py(), pi_pf.pz(), Mass_pi);
	AVF_KK_tlv = AVF_Kp_tlv + AVF_Km_tlv;
	AVF_Kppi_tlv = AVF_Kp_tlv + AVF_pi_tlv;
	AVF_Kmpi_tlv = AVF_Km_tlv + AVF_pi_tlv;
	AVF_Ds_tlv = AVF_Kp_tlv + AVF_Km_tlv + AVF_pi_tlv;

	GlobalVector AVF_displacement = AVF_Ds_vtx_pos - primvtx_pos;
	GlobalVector AVF_Ds_mtm(AVF_Ds_tlv.Px(), AVF_Ds_tlv.Py(), AVF_Ds_tlv.Pz());
	Measurement1D AVF_Ds_mtmDs_FD_XY = distanceXY(primvtx_pos, primvtx_cov, AVF_Ds_vtx_pos, AVF_Ds_vtx_cov);
	Measurement1D AVF_Ds_mtmDs_FD_3D = distance3D(primvtx_pos, primvtx_cov, AVF_Ds_vtx_pos, AVF_Ds_vtx_cov);

	ftree->AVF_KK_mass = AVF_KK_tlv.M();
	ftree->AVF_Kppi_mass = AVF_Kppi_tlv.M();
	ftree->AVF_Kmpi_mass = AVF_Kmpi_tlv.M();

	ftree->AVF_Ds_eta = AVF_Ds_tlv.Eta();
	ftree->AVF_Ds_phi = AVF_Ds_tlv.Phi();
	ftree->AVF_Ds_px = AVF_Ds_tlv.Px();
	ftree->AVF_Ds_py = AVF_Ds_tlv.Py();
	ftree->AVF_Ds_pz = AVF_Ds_tlv.Pz();
	ftree->AVF_Ds_pt = AVF_Ds_tlv.Pt();
	ftree->AVF_Ds_p = AVF_Ds_tlv.P();
	ftree->AVF_Ds_mass = AVF_Ds_tlv.M();;

	ftree->AVF_Ds_vx = AVF_Ds_vtx_pos.x();
	ftree->AVF_Ds_vy = AVF_Ds_vtx_pos.y();
	ftree->AVF_Ds_vz = AVF_Ds_vtx_pos.z();
	ftree->AVF_Ds_vxerr = AVF_Ds_vtx_cov.cxx();
	ftree->AVF_Ds_vyerr = AVF_Ds_vtx_cov.cyy();
	ftree->AVF_Ds_vzerr = AVF_Ds_vtx_cov.czz();
	ftree->AVF_Ds_chi2 = AVF_Ds_vertex.totalChiSquared();
	ftree->AVF_Ds_ndf = AVF_Ds_vertex.degreesOfFreedom();
	ftree->AVF_Ds_chi2ndf = AVF_Ds_vertex.totalChiSquared() / AVF_Ds_vertex.degreesOfFreedom();
	ftree->AVF_Ds_prob = TMath::Prob(AVF_Ds_vertex.totalChiSquared(), AVF_Ds_vertex.degreesOfFreedom());

	ftree->AVF_dira = AVF_displacement.dot(AVF_Ds_mtm) / (AVF_displacement.mag() * AVF_Ds_mtm.mag());
	ftree->AVF_Ds_FDxy = AVF_Ds_mtmDs_FD_XY.value();
	ftree->AVF_Ds_FDxyerr = AVF_Ds_mtmDs_FD_XY.error();
	ftree->AVF_Ds_FDxysig = AVF_Ds_mtmDs_FD_XY.significance();
	ftree->AVF_Ds_FD = AVF_Ds_mtmDs_FD_3D.value();
	ftree->AVF_Ds_FDerr = AVF_Ds_mtmDs_FD_3D.error();
	ftree->AVF_Ds_FDsig = AVF_Ds_mtmDs_FD_3D.significance();

	ftree->tree->Fill();
}

void ThreeFitters::endJob()
{
}

Bool_t ThreeFitters::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
{
	if (gp.numberOfMothers() == 0)
		return false;

	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother)
			continue;

		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 3)
		{
			Bool_t daughter1_flag = false;
			Bool_t daughter2_flag = false;
			Bool_t daughter3_flag = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id)
					daughter1_flag = true;
				if (mother->daughter(j)->pdgId() == daughter2id)
					daughter2_flag = true;
				if (mother->daughter(j)->pdgId() == daughter3id)
					daughter3_flag = true;
			}

			if (daughter1_flag && daughter2_flag && daughter3_flag)
				return true;
		}

		if (ThreeBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id, daughter3id))
			return true;
	}

	return false;
}

Bool_t ThreeFitters::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

		if (TwoBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id))
			return true;
	}

	return false;
}

Measurement1D ThreeFitters::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D ThreeFitters::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void ThreeFitters::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ThreeFitters);

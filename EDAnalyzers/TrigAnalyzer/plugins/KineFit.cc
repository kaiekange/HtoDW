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

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "EDAnalyzers/TrigAnalyzer/interface/KineTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class KineFit : public edm::one::EDAnalyzer<>
{
public:
	explicit KineFit(const edm::ParameterSet &);
	~KineFit() override {};

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

	const Float_t Mass_Ds = 1.96835;
	const Float_t Mass_phi = 1.019460;
	const Float_t Mass_K = 0.493677;
	const Float_t Mass_pi = 0.13957039;
	const Float_t Mass_mu = 0.1056583755;

	// std::map<TString, Int_t> triggerCounts_;
	// std::map<TString, Int_t> triggerTotal_;

	const edm::Service<TFileService> fs;
	KineTree *ftree;

	KinematicParticleFactoryFromTransientTrack factory;
	KinematicParticleVertexFitter KPVFitter;
	KinematicParticleFitter CSFitter;
	KinematicConstrainedVertexFitter KCVFitter;
};

KineFit::KineFit(const edm::ParameterSet &iConfig)
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

void KineFit::beginJob()
{
	ftree = new KineTree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
}

void KineFit::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
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

	ftree->Init();

	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;

	// ftree->primvtx_isValid = primvtx.isValid();
	// ftree->primvtx_isFake = primvtx.isFake();
	// ftree->primvtx_chi2 = primvtx.chi2();
	// ftree->primvtx_ndof = primvtx.ndof();
	// ftree->primvtx_chi2ndof = primvtx.chi2() / primvtx.ndof();
	// ftree->primvtx_vx = primvtx.x();
	// ftree->primvtx_vy = primvtx.y();
	// ftree->primvtx_vz = primvtx.z();
	// ftree->primvtx_vxerr = primvtx.xError();
	// ftree->primvtx_vyerr = primvtx.yError();
	// ftree->primvtx_vzerr = primvtx.zError();
	// GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	// GlobalError primvtx_cov(primvtx.covariance());

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

	if (Kp_match_idx < 0)
		return;
	if (Km_match_idx < 0)
		return;
	if (pi_match_idx < 0)
		return;
	if (mu_match_idx < 0)
		return;

	const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[Kp_match_idx];
	const pat::PackedCandidate &Km_pf = (*PFCandHandle)[Km_match_idx];
	const pat::PackedCandidate &pi_pf = (*PFCandHandle)[pi_match_idx];
	const pat::Muon &muon = (*MuonHandle)[mu_match_idx];

	Float_t chi2 = 0.0;
	Float_t ndf = 0.0;
	Float_t Mass_Ds_sigma = 0.00007;
	Float_t Mass_phi_sigma = 0.000016;
	Float_t Mass_K_sigma = 0.000013;
	Float_t Mass_pi_sigma = 0.00000017;
	Float_t Mass_mu_sigma = 0.0000000023;

	/*--------------------------------------------------
		build transient tracks
	--------------------------------------------------*/
	reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.pseudoTrack());
	reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.pseudoTrack());
	reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.pseudoTrack());

	// /*--------------------------------------------------
	// 	Kinematic fit on phi
	// --------------------------------------------------*/
	// std::vector<RefCountedKinematicParticle> phiDaughters;
	// phiDaughters.push_back(factory.particle(Kp_tt, Mass_K, chi2, ndf, Mass_K_sigma));
	// phiDaughters.push_back(factory.particle(Km_tt, Mass_K, chi2, ndf, Mass_K_sigma));

	// RefCountedKinematicTree phiTree = KPVFitter.fit(phiDaughters);
	// if (!phiTree->isValid())
	// 	return;
	// phiTree->movePointerToTheTop();
	// RefCountedKinematicParticle phi = phiTree->currentParticle();
	// ftree->phi_mass = phi->currentState().mass();

	// /*--------------------------------------------------
	// 	Mass constraint on phi
	// --------------------------------------------------*/
	// KinematicConstraint *phi_mass_constraint = new MassKinematicConstraint(Mass_phi, Mass_phi_sigma);
	// phiTree = CSFitter.fit(phi_mass_constraint, phiTree);
	// if (!phiTree->isValid())
	// 	return;
	// phiTree->movePointerToTheTop();
	// RefCountedKinematicParticle phi_refitted = phiTree->currentParticle();
	// RefCountedKinematicVertex phi_vertex_refitted = phiTree->currentDecayVertex();
	// if (!phi_vertex_refitted->vertexIsValid())
	// 	return;

	// ftree->phi_mass_refitted = phi_refitted->currentState().mass();
	// TLorentzVector phi_tlv_refitted;
	// phi_tlv_refitted.SetXYZM(phi_refitted->currentState().globalMomentum().x(), phi_refitted->currentState().globalMomentum().y(), phi_refitted->currentState().globalMomentum().z(), phi_refitted->currentState().mass());
	// ftree->phi_px = phi_tlv_refitted.Px();
	// ftree->phi_py = phi_tlv_refitted.Py();
	// ftree->phi_pz = phi_tlv_refitted.Pz();
	// ftree->phi_pt = phi_tlv_refitted.Pt();
	// ftree->phi_p = phi_tlv_refitted.P();
	// ftree->phi_eta = phi_tlv_refitted.Eta();
	// ftree->phi_phi = phi_tlv_refitted.Phi();

	// ftree->phi_vx = phi_vertex_refitted->position().x();
	// ftree->phi_vy = phi_vertex_refitted->position().y();
	// ftree->phi_vz = phi_vertex_refitted->position().z();
	// ftree->phi_chi2 = phi_vertex_refitted->chiSquared();
	// ftree->phi_ndf = phi_vertex_refitted->degreesOfFreedom();
	// ftree->phi_chi2ndf = phi_vertex_refitted->chiSquared() / phi_vertex_refitted->degreesOfFreedom();
	// ftree->phi_prob = TMath::Prob(ftree->phi_chi2, ftree->phi_ndf);

	// /*--------------------------------------------------
	// 	Kinematic fit on Ds
	// --------------------------------------------------*/
	// std::vector<RefCountedKinematicParticle> DsDaughters;
	// DsDaughters.push_back(phi_refitted);
	// DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi2, ndf, Mass_pi_sigma));
	// RefCountedKinematicTree DsTree = KPVFitter.fit(DsDaughters);
	// if (!DsTree->isValid())
	// 	return;
	// DsTree->movePointerToTheTop();
	// RefCountedKinematicParticle Ds = DsTree->currentParticle();
	// RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
	// if (!Ds_vertex->vertexIsValid())
	// 	return;

	// ftree->Ds_mass = Ds->currentState().mass();
	// TLorentzVector Ds_tlv;
	// Ds_tlv.SetXYZM(Ds->currentState().globalMomentum().x(), Ds->currentState().globalMomentum().y(), Ds->currentState().globalMomentum().z(), Ds->currentState().mass());
	// ftree->Ds_px = Ds_tlv.Px();
	// ftree->Ds_py = Ds_tlv.Py();
	// ftree->Ds_pz = Ds_tlv.Pz();
	// ftree->Ds_pt = Ds_tlv.Pt();
	// ftree->Ds_p = Ds_tlv.P();
	// ftree->Ds_eta = Ds_tlv.Eta();
	// ftree->Ds_phi = Ds_tlv.Phi();

	// ftree->Ds_vx = Ds_vertex->position().x();
	// ftree->Ds_vy = Ds_vertex->position().y();
	// ftree->Ds_vz = Ds_vertex->position().z();
	// ftree->Ds_chi2 = Ds_vertex->chiSquared();
	// ftree->Ds_ndf = Ds_vertex->degreesOfFreedom();
	// ftree->Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
	// ftree->Ds_prob = TMath::Prob(ftree->Ds_chi2, ftree->Ds_ndf);

	ParticleMass phi_mass = 1.019460;
	std::vector<RefCountedKinematicParticle> DsDaughters;
	DsDaughters.push_back(factory.particle(Kp_tt, Mass_K, chi2, ndf, Mass_K_sigma));
	DsDaughters.push_back(factory.particle(Km_tt, Mass_K, chi2, ndf, Mass_K_sigma));
	DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi2, ndf, Mass_pi_sigma));

	MultiTrackKinematicConstraint *phi_mass_constraint = new TwoTrackMassKinematicConstraint(phi_mass);
	RefCountedKinematicTree DsTree = KCVFitter.fit(DsDaughters, phi_mass_constraint);
	if (!DsTree->isValid())
		return;
	DsTree->movePointerToTheTop();
	RefCountedKinematicParticle Ds = DsTree->currentParticle();
	RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
	if (!Ds_vertex->vertexIsValid())
		return;
	ftree->Ds_mass = Ds->currentState().mass();
	TLorentzVector Ds_tlv;
	Ds_tlv.SetXYZM(Ds->currentState().globalMomentum().x(), Ds->currentState().globalMomentum().y(), Ds->currentState().globalMomentum().z(), Ds->currentState().mass());
	ftree->Ds_px = Ds_tlv.Px();
	ftree->Ds_py = Ds_tlv.Py();
	ftree->Ds_pz = Ds_tlv.Pz();
	ftree->Ds_pt = Ds_tlv.Pt();
	ftree->Ds_p = Ds_tlv.P();
	ftree->Ds_eta = Ds_tlv.Eta();
	ftree->Ds_phi = Ds_tlv.Phi();

	ftree->Ds_vx = Ds_vertex->position().x();
	ftree->Ds_vy = Ds_vertex->position().y();
	ftree->Ds_vz = Ds_vertex->position().z();
	ftree->Ds_chi2 = Ds_vertex->chiSquared();
	ftree->Ds_ndf = Ds_vertex->degreesOfFreedom();
	ftree->Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
	ftree->Ds_prob = TMath::Prob(ftree->Ds_chi2, ftree->Ds_ndf);

	ftree->tree->Fill();
}

Bool_t KineFit::hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

void KineFit::endJob()
{
	// std::cout << nEventsProcessed_ << ", Kp: " << nEventsValidKp_ << ", Km " << nEventsValidKm_ << ", pi" << nEventsValidpi_ << ", mu" << nEventsValidmu_ << std::endl;

	// edm::LogPrint("KineFit") << "\n========== Trigger Statistics ==========";
	// edm::LogPrint("KineFit") << std::setw(60) << "Trigger Name"
	// 						   << std::setw(15) << "Fired"
	// 						   << std::setw(15) << "Total"
	// 						   << std::setw(15) << "Efficiency";
	// edm::LogPrint("KineFit") << std::string(105, '-');

	// std::vector<std::pair<TString, double>> triggerEfficiencies;

	// for (const auto &entry : triggerCounts_)
	// {
	// 	const TString trigName = entry.first;
	// 	Int_t fired = entry.second;
	// 	Int_t total = triggerTotal_[trigName];
	// 	double efficiency = (total > 0) ? (double)fired / total * 100.0 : 0.0;

	// 	triggerEfficiencies.push_back({trigName, efficiency});
	// }

	// std::sort(triggerEfficiencies.begin(), triggerEfficiencies.end(),
	// 		  [](const auto &a, const auto &b)
	// 		  { return a.second > b.second; });

	// for (const auto &item : triggerEfficiencies)
	// {
	// 	const TString &trigName = item.first;
	// 	Int_t fired = triggerCounts_[trigName];
	// 	Int_t total = triggerTotal_[trigName];
	// 	double efficiency = item.second;

	// 	edm::LogPrint("KineFit") << std::setw(80) << trigName
	// 							   << std::setw(15) << fired
	// 							   << std::setw(15) << total
	// 							   << std::setw(14) << std::fixed << std::setprecision(2)
	// 							   << efficiency << "%";
	// }

	// edm::LogPrint("KineFit") << "========================================\n";
}

void KineFit::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(KineFit);

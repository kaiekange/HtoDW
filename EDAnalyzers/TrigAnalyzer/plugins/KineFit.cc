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
	GlobalPoint primvtx_vtx(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_vtxerr(primvtx.covariance());

	Int_t Kp_match_idx = null;
	Int_t Km_match_idx = null;
	Int_t pi_match_idx = null;
	Int_t mu_match_idx = null;

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
	// const pat::Muon &muon = (*MuonHandle)[mu_match_idx];

	Float_t chi = 0.0;
	Float_t ndf = 0.0;
	Float_t Mass_K_sigma = 0.000013;
	Float_t Mass_pi_sigma = 0.00000017;
	Float_t Mass_mu_sigma = 0.0000000023;

	//--------------------------------------------------//
	// build transient tracks
	//--------------------------------------------------//
	reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.bestTrack());
	reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.bestTrack());
	reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.bestTrack());

	//--------------------------------------------------//
	// Kinematic fit on phi
	//--------------------------------------------------//
	std::vector<RefCountedKinematicParticle> phiDaughters;
	phiDaughters.push_back(factory.particle(Kp_tt, Mass_K, chi, ndf, Mass_K_sigma));
	phiDaughters.push_back(factory.particle(Km_tt, Mass_K, chi, ndf, Mass_K_sigma));
	RefCountedKinematicTree phiTree = KPVFitter.fit(phiDaughters);
	if (!phiTree->isValid())
		return;
	phiTree->movePointerToTheTop();

	// phiTree vertex
	RefCountedKinematicVertex phi_vertex = phiTree->currentDecayVertex();
	if (!phi_vertex->vertexIsValid())
		return;
	GlobalPoint phi_vtx = phi_vertex->position();
	GlobalError phi_vtxerr = phi_vertex->error();

	// phiTree particles: phi
	RefCountedKinematicParticle phiTree_phi = phiTree->currentParticle();
	GlobalVector phiTree_phi_mtm = phiTree_phi->currentState().globalMomentum();
	// phiTree particles: K+
	phiTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle phiTree_Kp = phiTree->currentParticle();
	GlobalVector phiTree_Kp_mtm = phiTree_Kp->currentState().globalMomentum();
	// phiTree particles: K-
	phiTree->movePointerToTheNextChild();
	RefCountedKinematicParticle phiTree_Km = phiTree->currentParticle();
	GlobalVector phiTree_Km_mtm = phiTree_Km->currentState().globalMomentum();

	//--------------------------------------------------//
	// Kinematic fit on Ds
	//--------------------------------------------------//
	std::vector<RefCountedKinematicParticle> DsDaughters;
	DsDaughters.push_back(phiTree_phi);
	DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi, ndf, Mass_pi_sigma));
	RefCountedKinematicTree DsTree = KPVFitter.fit(DsDaughters);
	if (!DsTree->isValid())
		return;
	DsTree->movePointerToTheTop();

	// DsTree vertex
	RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
	if (!Ds_vertex->vertexIsValid())
		return;
	GlobalPoint Ds_vtx = Ds_vertex->position();
	GlobalError Ds_vtxerr = Ds_vertex->error();

	// DsTree particles: Ds
	RefCountedKinematicParticle DsTree_Ds = DsTree->currentParticle();
	GlobalVector DsTree_Ds_mtm = DsTree_Ds->currentState().globalMomentum();
	const reco::TransientTrack DsTree_Ds_tt = DsTree_Ds->refittedTransientTrack();
	// DsTree particles: phi
	DsTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle DsTree_phi = DsTree->currentParticle();
	GlobalVector DsTree_phi_mtm = DsTree_phi->currentState().globalMomentum();
	const reco::TransientTrack DsTree_phi_tt = DsTree_phi->refittedTransientTrack();
	// DsTree particles: K+
	DsTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle DsTree_Kp = DsTree->currentParticle();
	GlobalVector DsTree_Kp_mtm = DsTree_Kp->currentState().globalMomentum();
	const reco::TransientTrack DsTree_Kp_tt = DsTree_Kp->refittedTransientTrack();
	// DsTree particles: K-
	DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle DsTree_Km = DsTree->currentParticle();
	GlobalVector DsTree_Km_mtm = DsTree_Km->currentState().globalMomentum();
	const reco::TransientTrack DsTree_Km_tt = DsTree_Km->refittedTransientTrack();
	// DsTree particles: pi
	DsTree->movePointerToTheTop();
	DsTree->movePointerToTheFirstChild();
	DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle DsTree_pi = DsTree->currentParticle();
	GlobalVector DsTree_pi_mtm = DsTree_pi->currentState().globalMomentum();
	const reco::TransientTrack DsTree_pi_tt = DsTree_pi->refittedTransientTrack();

	std::pair<Bool_t, Measurement1D> Kp_ipresult = IPTools::absoluteImpactParameter3D(DsTree_Kp_tt, primvtx);
	std::pair<Bool_t, Measurement1D> Km_ipresult = IPTools::absoluteImpactParameter3D(DsTree_Km_tt, primvtx);
	std::pair<Bool_t, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(DsTree_pi_tt, primvtx);
	std::pair<Bool_t, Measurement1D> phi_ipresult = IPTools::absoluteImpactParameter3D(DsTree_phi_tt, primvtx);
	std::pair<Bool_t, Measurement1D> Ds_ipresult = IPTools::absoluteImpactParameter3D(DsTree_Ds_tt, primvtx);

	GlobalVector displacement = Ds_vtx - primvtx_vtx;

	// VertexDistanceXY vtxdistXY;
	// Measurement1D Ds_FDXY = vtxdistXY.distance(primvtx_vtx, primvtx_vtxerr, Ds_vtx, Ds_vtxerr);
	// VertexDistance3D vtxdist3D;
	// Measurement1D Ds_FD_3D = vtxdist3D.distance(primvtx_vtx, primvtx_vtxerr, Ds_vtx, Ds_vtxerr);

	Float_t IsoR03_sumChargedHadronPt = 0;
	Float_t IsoR03_sumNeutralHadronEt = 0;
	Float_t IsoR03_sumPhotonEt = 0;
	Float_t IsoR03_sumPUPt = 0;
	Float_t IsoR04_sumChargedHadronPt = 0;
	Float_t IsoR04_sumNeutralHadronEt = 0;
	Float_t IsoR04_sumPhotonEt = 0;
	Float_t IsoR04_sumPUPt = 0;

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		if (i == Kp_match_idx)
			continue;
		if (i == Km_match_idx)
			continue;
		if (i == pi_match_idx)
			continue;
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];

		if (reco::deltaR(DsTree_Ds_mtm.eta(), DsTree_Ds_mtm.phi(), pf.eta(), pf.phi()) < 0.3)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				IsoR03_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				IsoR03_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				IsoR03_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				IsoR03_sumPUPt += pf.pt();
		}
		if (reco::deltaR(DsTree_Ds_mtm.eta(), DsTree_Ds_mtm.phi(), pf.eta(), pf.phi()) < 0.4)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				IsoR04_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				IsoR04_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				IsoR04_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				IsoR04_sumPUPt += pf.pt();
		}
	}
	Float_t IsoR03_PFIso = IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumNeutralHadronEt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt);
	Float_t IsoR03_RelIso = (IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumNeutralHadronEt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt)) / DsTree_Ds_mtm.perp();
	Float_t IsoR04_PFIso = IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumNeutralHadronEt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt);
	Float_t IsoR04_RelIso = (IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumNeutralHadronEt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt)) / DsTree_Ds_mtm.perp();

	//--------------------------------------------------//
	// Store in tree
	//--------------------------------------------------//
	ftree->Kp_charge = Kp_pf.charge();
	ftree->Kp_eta = Kp_pf.eta();
	ftree->Kp_phi = Kp_pf.phi();
	ftree->Kp_px = Kp_pf.px();
	ftree->Kp_py = Kp_pf.py();
	ftree->Kp_pz = Kp_pf.pz();
	ftree->Kp_pt = Kp_pf.pt();
	ftree->Kp_p = Kp_pf.p();

	ftree->Km_charge = Km_pf.charge();
	ftree->Km_eta = Km_pf.eta();
	ftree->Km_phi = Km_pf.phi();
	ftree->Km_px = Km_pf.px();
	ftree->Km_py = Km_pf.py();
	ftree->Km_pz = Km_pf.pz();
	ftree->Km_pt = Km_pf.pt();
	ftree->Km_p = Km_pf.p();

	ftree->pi_charge = pi_pf.charge();
	ftree->pi_eta = pi_pf.eta();
	ftree->pi_phi = pi_pf.phi();
	ftree->pi_px = pi_pf.px();
	ftree->pi_py = pi_pf.py();
	ftree->pi_pz = pi_pf.pz();
	ftree->pi_pt = pi_pf.pt();
	ftree->pi_p = pi_pf.p();

	ftree->phiFit_Kp_charge = Kp_pf.charge();
	ftree->phiFit_Kp_eta = phiTree_Kp_mtm.eta();
	ftree->phiFit_Kp_phi = phiTree_Kp_mtm.phi();
	ftree->phiFit_Kp_px = phiTree_Kp_mtm.x();
	ftree->phiFit_Kp_py = phiTree_Kp_mtm.y();
	ftree->phiFit_Kp_pz = phiTree_Kp_mtm.z();
	ftree->phiFit_Kp_pt = phiTree_Kp_mtm.perp();
	ftree->phiFit_Kp_p = phiTree_Kp_mtm.mag();
	ftree->phiFit_Kp_mass = phiTree_Kp->currentState().mass();

	ftree->phiFit_Km_charge = Km_pf.charge();
	ftree->phiFit_Km_eta = phiTree_Km_mtm.eta();
	ftree->phiFit_Km_phi = phiTree_Km_mtm.phi();
	ftree->phiFit_Km_px = phiTree_Km_mtm.x();
	ftree->phiFit_Km_py = phiTree_Km_mtm.y();
	ftree->phiFit_Km_pz = phiTree_Km_mtm.z();
	ftree->phiFit_Km_pt = phiTree_Km_mtm.perp();
	ftree->phiFit_Km_p = phiTree_Km_mtm.mag();
	ftree->phiFit_Km_mass = phiTree_Km->currentState().mass();

	ftree->phiFit_phi_charge = Kp_pf.charge() + Km_pf.charge();
	ftree->phiFit_phi_eta = phiTree_phi_mtm.eta();
	ftree->phiFit_phi_phi = phiTree_phi_mtm.phi();
	ftree->phiFit_phi_px = phiTree_phi_mtm.x();
	ftree->phiFit_phi_py = phiTree_phi_mtm.y();
	ftree->phiFit_phi_pz = phiTree_phi_mtm.z();
	ftree->phiFit_phi_pt = phiTree_phi_mtm.perp();
	ftree->phiFit_phi_p = phiTree_phi_mtm.mag();
	ftree->phiFit_phi_mass = phiTree_phi->currentState().mass();

	ftree->DsFit_Kp_charge = Kp_pf.charge();
	ftree->DsFit_Kp_eta = DsTree_Kp_mtm.eta();
	ftree->DsFit_Kp_phi = DsTree_Kp_mtm.phi();
	ftree->DsFit_Kp_px = DsTree_Kp_mtm.x();
	ftree->DsFit_Kp_py = DsTree_Kp_mtm.y();
	ftree->DsFit_Kp_pz = DsTree_Kp_mtm.z();
	ftree->DsFit_Kp_pt = DsTree_Kp_mtm.perp();
	ftree->DsFit_Kp_p = DsTree_Kp_mtm.mag();
	ftree->DsFit_Kp_mass = DsTree_Kp->currentState().mass();

	ftree->DsFit_Km_charge = Km_pf.charge();
	ftree->DsFit_Km_eta = DsTree_Km_mtm.eta();
	ftree->DsFit_Km_phi = DsTree_Km_mtm.phi();
	ftree->DsFit_Km_px = DsTree_Km_mtm.x();
	ftree->DsFit_Km_py = DsTree_Km_mtm.y();
	ftree->DsFit_Km_pz = DsTree_Km_mtm.z();
	ftree->DsFit_Km_pt = DsTree_Km_mtm.perp();
	ftree->DsFit_Km_p = DsTree_Km_mtm.mag();
	ftree->DsFit_Km_mass = DsTree_Km->currentState().mass();

	ftree->DsFit_pi_charge = pi_pf.charge();
	ftree->DsFit_pi_eta = DsTree_pi_mtm.eta();
	ftree->DsFit_pi_phi = DsTree_pi_mtm.phi();
	ftree->DsFit_pi_px = DsTree_pi_mtm.x();
	ftree->DsFit_pi_py = DsTree_pi_mtm.y();
	ftree->DsFit_pi_pz = DsTree_pi_mtm.z();
	ftree->DsFit_pi_pt = DsTree_pi_mtm.perp();
	ftree->DsFit_pi_p = DsTree_pi_mtm.mag();
	ftree->DsFit_pi_mass = DsTree_pi->currentState().mass();

	ftree->DsFit_phi_charge = Kp_pf.charge() + Km_pf.charge();
	ftree->DsFit_phi_eta = DsTree_phi_mtm.eta();
	ftree->DsFit_phi_phi = DsTree_phi_mtm.phi();
	ftree->DsFit_phi_px = DsTree_phi_mtm.x();
	ftree->DsFit_phi_py = DsTree_phi_mtm.y();
	ftree->DsFit_phi_pz = DsTree_phi_mtm.z();
	ftree->DsFit_phi_pt = DsTree_phi_mtm.perp();
	ftree->DsFit_phi_p = DsTree_phi_mtm.mag();
	ftree->DsFit_phi_mass = DsTree_phi->currentState().mass();

	ftree->DsFit_Ds_charge = Kp_pf.charge() + Km_pf.charge() + pi_pf.charge();
	ftree->DsFit_Ds_eta = DsTree_Ds_mtm.eta();
	ftree->DsFit_Ds_phi = DsTree_Ds_mtm.phi();
	ftree->DsFit_Ds_px = DsTree_Ds_mtm.x();
	ftree->DsFit_Ds_py = DsTree_Ds_mtm.y();
	ftree->DsFit_Ds_pz = DsTree_Ds_mtm.z();
	ftree->DsFit_Ds_pt = DsTree_Ds_mtm.perp();
	ftree->DsFit_Ds_p = DsTree_Ds_mtm.mag();
	ftree->DsFit_Ds_mass = DsTree_Ds->currentState().mass();

	ftree->phi_vx = phi_vtx.x();
	ftree->phi_vy = phi_vtx.y();
	ftree->phi_vz = phi_vtx.z();
	ftree->phi_vxerr = phi_vtxerr.cxx();
	ftree->phi_vyerr = phi_vtxerr.cyy();
	ftree->phi_vzerr = phi_vtxerr.czz();
	ftree->phi_chi2 = phi_vertex->chiSquared();
	ftree->phi_ndf = phi_vertex->degreesOfFreedom();
	ftree->phi_chi2ndf = phi_vertex->chiSquared() / phi_vertex->degreesOfFreedom();
	ftree->phi_prob = TMath::Prob(phi_vertex->chiSquared(), phi_vertex->degreesOfFreedom());

	ftree->Ds_vx = Ds_vtx.x();
	ftree->Ds_vy = Ds_vtx.y();
	ftree->Ds_vz = Ds_vtx.z();
	ftree->Ds_vxerr = Ds_vtxerr.cxx();
	ftree->Ds_vyerr = Ds_vtxerr.cyy();
	ftree->Ds_vzerr = Ds_vtxerr.czz();
	ftree->Ds_chi2 = Ds_vertex->chiSquared();
	ftree->Ds_ndf = Ds_vertex->degreesOfFreedom();
	ftree->Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
	ftree->Ds_prob = TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom());

	ftree->dxy_phi_Ds = std::hypot(phi_vtx.x() - Ds_vtx.x(), phi_vtx.y() - Ds_vtx.y());
	ftree->dz_phi_Ds = abs(phi_vtx.z() - Ds_vtx.z());

	ftree->dR_Kp_Km = reco::deltaR(Kp_pf, Km_pf);
	ftree->dR_Kp_pi = reco::deltaR(Kp_pf, pi_pf);
	ftree->dR_Km_pi = reco::deltaR(Km_pf, pi_pf);

	ftree->phiFit_dR_Kp_Km = reco::deltaR(phiTree_Kp_mtm, phiTree_Km_mtm);
	ftree->phiFit_dR_Kp_pi = reco::deltaR(phiTree_Kp_mtm.eta(), phiTree_Kp_mtm.phi(), pi_pf.eta(), pi_pf.phi());
	ftree->phiFit_dR_Kp_phi = reco::deltaR(phiTree_Kp_mtm, phiTree_phi_mtm);
	ftree->phiFit_dR_Km_pi = reco::deltaR(phiTree_Km_mtm.eta(), phiTree_Km_mtm.phi(), pi_pf.eta(), pi_pf.phi());
	ftree->phiFit_dR_Km_phi = reco::deltaR(phiTree_Km_mtm, phiTree_phi_mtm);
	ftree->phiFit_dR_pi_phi = reco::deltaR(phiTree_phi_mtm.eta(), phiTree_phi_mtm.phi(), pi_pf.eta(), pi_pf.phi());

	ftree->DsFit_dR_Kp_Km = reco::deltaR(DsTree_Kp_mtm, DsTree_Km_mtm);
	ftree->DsFit_dR_Kp_pi = reco::deltaR(DsTree_Kp_mtm, DsTree_pi_mtm);
	ftree->DsFit_dR_Kp_phi = reco::deltaR(DsTree_Kp_mtm, DsTree_phi_mtm);
	ftree->DsFit_dR_Kp_Ds = reco::deltaR(DsTree_Kp_mtm, DsTree_Ds_mtm);
	ftree->DsFit_dR_Km_pi = reco::deltaR(DsTree_Km_mtm, DsTree_pi_mtm);
	ftree->DsFit_dR_Km_phi = reco::deltaR(DsTree_Km_mtm, DsTree_phi_mtm);
	ftree->DsFit_dR_Km_Ds = reco::deltaR(DsTree_Km_mtm, DsTree_Ds_mtm);
	ftree->DsFit_dR_pi_phi = reco::deltaR(DsTree_pi_mtm, DsTree_phi_mtm);
	ftree->DsFit_dR_pi_Ds = reco::deltaR(DsTree_pi_mtm, DsTree_Ds_mtm);
	ftree->DsFit_dR_phi_Ds = reco::deltaR(DsTree_phi_mtm, DsTree_Ds_mtm);

	if (Kp_ipresult.first)
	{
		ftree->Kp_ip = Kp_ipresult.second.value();
		ftree->Kp_iperr = Kp_ipresult.second.error();
		ftree->Kp_ipsig = Kp_ipresult.second.significance();
	}
	if (Km_ipresult.first)
	{
		ftree->Km_ip = Km_ipresult.second.value();
		ftree->Km_iperr = Km_ipresult.second.error();
		ftree->Km_ipsig = Km_ipresult.second.significance();
	}
	if (pi_ipresult.first)
	{
		ftree->pi_ip = pi_ipresult.second.value();
		ftree->pi_iperr = pi_ipresult.second.error();
		ftree->pi_ipsig = pi_ipresult.second.significance();
	}
	if (phi_ipresult.first)
	{
		ftree->phi_ip = phi_ipresult.second.value();
		ftree->phi_iperr = phi_ipresult.second.error();
		ftree->phi_ipsig = phi_ipresult.second.significance();
	}
	if (Ds_ipresult.first)
	{
		ftree->Ds_ip = Ds_ipresult.second.value();
		ftree->Ds_iperr = Ds_ipresult.second.error();
		ftree->Ds_ipsig = Ds_ipresult.second.significance();
	}

	ftree->dira = displacement.dot(DsTree_Ds_mtm) / (displacement.mag() * DsTree_Ds_mtm.mag());

	// ftree->Ds_FDxy = Ds_FDXY.value();
	// ftree->Ds_FDxyerr = Ds_FDXY.error();
	// ftree->Ds_FDxysig = Ds_FDXY.significance();

	// ftree->Ds_FD = Ds_FD_3D.value();
	// ftree->Ds_FDerr = Ds_FD_3D.error();
	// ftree->Ds_FDsig = Ds_FD_3D.significance();

	ftree->Ds_IsoR03_sumChargedHadronPt = IsoR03_sumChargedHadronPt;
	ftree->Ds_IsoR03_sumNeutralHadronEt = IsoR03_sumNeutralHadronEt;
	ftree->Ds_IsoR03_sumPhotonEt = IsoR03_sumPhotonEt;
	ftree->Ds_IsoR03_sumPUPt = IsoR03_sumPUPt;
	ftree->Ds_IsoR03_PFIso = IsoR03_PFIso;
	ftree->Ds_IsoR03_RelIso = IsoR03_RelIso;

	ftree->Ds_IsoR04_sumChargedHadronPt = IsoR04_sumChargedHadronPt;
	ftree->Ds_IsoR04_sumNeutralHadronEt = IsoR04_sumNeutralHadronEt;
	ftree->Ds_IsoR04_sumPhotonEt = IsoR04_sumPhotonEt;
	ftree->Ds_IsoR04_sumPUPt = IsoR04_sumPUPt;
	ftree->Ds_IsoR04_PFIso = IsoR04_PFIso;
	ftree->Ds_IsoR04_RelIso = IsoR04_RelIso;

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

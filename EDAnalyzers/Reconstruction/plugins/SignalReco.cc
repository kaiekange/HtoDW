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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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

#include "EDAnalyzers/Reconstruction/interface/SignalRecoTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class SignalReco : public edm::one::EDAnalyzer<>
{
public:
	explicit SignalReco(const edm::ParameterSet &);
	~SignalReco() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;
	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<pat::ElectronCollection> ElectronToken_;
	edm::EDGetTokenT<pat::METCollection> MetToken_;
	edm::EDGetTokenT<pat::JetCollection> JetToken_;
	edm::EDGetTokenT<edm::ValueMap<Int_t>> tightJetIdToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
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

	const edm::Service<TFileService> fs;
	SignalTree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsTriggered_;
	unsigned long nEventsReconstructed_;
	unsigned long nEventsMatch_;
};

SignalReco::SignalReco(const edm::ParameterSet &iConfig)
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

	edm::InputTag ElectronTag_("slimmedElectrons");
	ElectronToken_ = consumes<pat::ElectronCollection>(ElectronTag_);

	// edm::InputTag MetTag_("slimmedMETsPuppi");
	// MetToken_ = consumes<pat::METCollection>(MetTag_);

	// edm::InputTag JetTag_("slimmedJetsPuppi");
	// JetToken_ = consumes<pat::JetCollection>(JetTag_);

	edm::InputTag MetTag_("jecAppliedMetProducer", "CorrectedMet");
	MetToken_ = consumes<pat::METCollection>(MetTag_);

	edm::InputTag JetTag_("jecAppliedJetProducer", "CorrectedAK4");
	JetToken_ = consumes<pat::JetCollection>(JetTag_);

	edm::InputTag tightJetIdTag_("tightJetId");
	tightJetIdToken_ = consumes<edm::ValueMap<Int_t>>(tightJetIdTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);
}

void SignalReco::beginJob()
{
	ftree = new SignalTree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
	InfoTree = fs->make<TTree>("RunInfo", "RunInfo");
	InfoTree->Branch("nEventsProcessed_", &nEventsProcessed_);
	InfoTree->Branch("nEventsTriggered_", &nEventsTriggered_);
	InfoTree->Branch("nEventsReconstructed_", &nEventsReconstructed_);
	InfoTree->Branch("nEventsMatch_", &nEventsMatch_);

	nEventsProcessed_ = 0;
	nEventsTriggered_ = 0;
	nEventsReconstructed_ = 0;
	nEventsMatch_ = 0;
}

void SignalReco::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

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

	edm::Handle<pat::ElectronCollection> ElectronHandle;
	iEvent.getByToken(ElectronToken_, ElectronHandle);

	edm::Handle<pat::METCollection> MetHandle;
	iEvent.getByToken(MetToken_, MetHandle);

	edm::Handle<pat::JetCollection> JetHandle;
	iEvent.getByToken(JetToken_, JetHandle);

	edm::Handle<edm::ValueMap<Int_t>> tightJetIdHandle;
	iEvent.getByToken(tightJetIdToken_, tightJetIdHandle);

	edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
	iEvent.getByToken(associationToken_, associationHandle);

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
	nEventsTriggered_++;
	//==================================================//
	//--------------------------------------------------//
	// Gen match
	//--------------------------------------------------//
	//==================================================//
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

	const pat::PackedCandidate &best_Kp_pf = (*PFCandHandle)[best_Kp_idx];
	const pat::PackedCandidate &best_Km_pf = (*PFCandHandle)[best_Km_idx];
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
	if (best_mu_idx == -1)
		return;

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];

		if (muon.charge() == best_pi_pf.charge())
		{
			if (muon.isTightMuon(primvtx))
				ftree->num_ss_mu_tight++;
			if (muon.isMediumMuon())
				ftree->num_ss_mu_medium++;
			if (muon.isLooseMuon())
				ftree->num_ss_mu_loose++;
			if (muon.isSoftMuon(primvtx))
				ftree->num_ss_mu_soft++;
			if (muon.isHighPtMuon(primvtx))
				ftree->num_ss_mu_highpt++;

			if (muon.isTightMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_ss_mu_tightglobal++;
			if (muon.isMediumMuon() && muon.isGlobalMuon())
				ftree->num_ss_mu_mediumglobal++;
			if (muon.isLooseMuon() && muon.isGlobalMuon())
				ftree->num_ss_mu_looseglobal++;
			if (muon.isSoftMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_ss_mu_softglobal++;
			if (muon.isHighPtMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_ss_mu_highptglobal++;
		}
		else if (muon.charge() == -best_pi_pf.charge())
		{
			if (muon.isTightMuon(primvtx))
				ftree->num_os_mu_tight++;
			if (muon.isMediumMuon())
				ftree->num_os_mu_medium++;
			if (muon.isLooseMuon())
				ftree->num_os_mu_loose++;
			if (muon.isSoftMuon(primvtx))
				ftree->num_os_mu_soft++;
			if (muon.isHighPtMuon(primvtx))
				ftree->num_os_mu_highpt++;

			if (muon.isTightMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_os_mu_tightglobal++;
			if (muon.isMediumMuon() && muon.isGlobalMuon())
				ftree->num_os_mu_mediumglobal++;
			if (muon.isLooseMuon() && muon.isGlobalMuon())
				ftree->num_os_mu_looseglobal++;
			if (muon.isSoftMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_os_mu_softglobal++;
			if (muon.isHighPtMuon(primvtx) && muon.isGlobalMuon())
				ftree->num_os_mu_highptglobal++;
		}
	}

	for (Int_t idx_el = 0; idx_el < static_cast<Int_t>(ElectronHandle->size()); idx_el++)
	{
		const pat::Electron &electron = (*ElectronHandle)[idx_el];
		if (electron.charge() == best_pi_pf.charge())
		{
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-loose"))
				ftree->num_ss_el_cutloose++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-medium"))
				ftree->num_ss_el_cutmedium++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
				ftree->num_ss_el_cuttight++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-veto"))
				ftree->num_ss_el_cutveto++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp80"))
				ftree->num_ss_el_mvaiso80++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp90"))
				ftree->num_ss_el_mvaiso90++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-noIso-V1-wp80"))
				ftree->num_ss_el_mvanoiso80++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-noIso-V1-wp90"))
				ftree->num_ss_el_mvanoiso90++;
		}
		else if (electron.charge() == -best_pi_pf.charge())
		{
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-loose"))
				ftree->num_os_el_cutloose++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-medium"))
				ftree->num_os_el_cutmedium++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
				ftree->num_os_el_cuttight++;
			if (electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-veto"))
				ftree->num_os_el_cutveto++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp80"))
				ftree->num_os_el_mvaiso80++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp90"))
				ftree->num_os_el_mvaiso90++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-noIso-V1-wp80"))
				ftree->num_os_el_mvanoiso80++;
			if (electron.electronID("mvaEleID-RunIIIWinter22-noIso-V1-wp90"))
				ftree->num_os_el_mvanoiso90++;
		}
	}

	const pat::Muon &muon = (*MuonHandle)[best_mu_idx];

	//==================================================//
	//--------------------------------------------------//
	// MET
	//--------------------------------------------------//
	//==================================================//
	const pat::MET &Met = MetHandle->front();

	nEventsReconstructed_++;

	if ((best_Kp_idx == Kp_match_idx) && (best_Km_idx == Km_match_idx) && (best_pi_idx == pi_match_idx) && (best_mu_idx == mu_match_idx))
		nEventsMatch_++;

	//==================================================//
	//--------------------------------------------------//
	// Reconstruction
	//--------------------------------------------------//
	//==================================================//
	best_DsTree->movePointerToTheTop();
	RefCountedKinematicVertex Ds_vertex = best_DsTree->currentDecayVertex();
	GlobalPoint Ds_vtx_pos = Ds_vertex->position();
	GlobalError Ds_vtx_cov = Ds_vertex->error();

	// DsTree particles: Ds
	RefCountedKinematicParticle Ds = best_DsTree->currentParticle();
	GlobalVector Ds_mtm = Ds->currentState().globalMomentum();
	const reco::TransientTrack Ds_rett = Ds->refittedTransientTrack();
	// DsTree particles: K+
	best_DsTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle Kp = best_DsTree->currentParticle();
	GlobalVector Kp_mtm = Kp->currentState().globalMomentum();
	const reco::TransientTrack Kp_rett = Kp->refittedTransientTrack();
	// DsTree particles: K-
	best_DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle Km = best_DsTree->currentParticle();
	GlobalVector Km_mtm = Km->currentState().globalMomentum();
	const reco::TransientTrack Km_rett = Km->refittedTransientTrack();
	// DsTree particles: pi
	best_DsTree->movePointerToTheNextChild();
	RefCountedKinematicParticle pi = best_DsTree->currentParticle();
	GlobalVector pi_mtm = pi->currentState().globalMomentum();
	const reco::TransientTrack pi_rett = pi->refittedTransientTrack();

	// Construct phi from
	GlobalVector phi_mtm = Kp_mtm + Km_mtm;
	TLorentzVector Kp_tlv, Km_tlv, phi_tlv;
	Kp_tlv.SetXYZM(Kp_mtm.x(), Kp_mtm.y(), Kp_mtm.z(), Mass_K);
	Km_tlv.SetXYZM(Km_mtm.x(), Km_mtm.y(), Km_mtm.z(), Mass_K);
	phi_tlv = Kp_tlv + Km_tlv;

	//--------------------------------------------------//
	// Kinematic variables
	//--------------------------------------------------//
	// IP
	std::pair<Bool_t, Measurement1D> Kp_ipresult = IPTools::absoluteImpactParameter3D(Kp_rett, primvtx);
	std::pair<Bool_t, Measurement1D> Km_ipresult = IPTools::absoluteImpactParameter3D(Km_rett, primvtx);
	std::pair<Bool_t, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(pi_rett, primvtx);
	std::pair<Bool_t, Measurement1D> Ds_ipresult = IPTools::absoluteImpactParameter3D(Ds_rett, primvtx);

	// PV-SV separation
	GlobalVector displacement = Ds_vtx_pos - primvtx_pos;

	Measurement1D Ds_FD_XY = distanceXY(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);
	Measurement1D Ds_FD_3D = distance3D(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);

	// Ds Isolation
	Float_t Ds_IsoR03_sumChargedHadronPt = 0;
	Float_t Ds_IsoR03_sumNeutralHadronEt = 0;
	Float_t Ds_IsoR03_sumPhotonEt = 0;
	Float_t Ds_IsoR03_sumPUPt = 0;
	Float_t Ds_IsoR04_sumChargedHadronPt = 0;
	Float_t Ds_IsoR04_sumNeutralHadronEt = 0;
	Float_t Ds_IsoR04_sumPhotonEt = 0;
	Float_t Ds_IsoR04_sumPUPt = 0;
	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		if (i == best_Kp_idx)
			continue;
		if (i == best_Km_idx)
			continue;
		if (i == best_pi_idx)
			continue;
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];

		if (reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), pf.eta(), pf.phi()) < 0.3)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				Ds_IsoR03_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				Ds_IsoR03_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				Ds_IsoR03_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				Ds_IsoR03_sumPUPt += pf.pt();
		}
		if (reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), pf.eta(), pf.phi()) < 0.4)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				Ds_IsoR04_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				Ds_IsoR04_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				Ds_IsoR04_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				Ds_IsoR04_sumPUPt += pf.pt();
		}
	}
	Float_t Ds_IsoR03_PFIso = Ds_IsoR03_sumChargedHadronPt + std::max(0.0, Ds_IsoR03_sumNeutralHadronEt + Ds_IsoR03_sumPhotonEt - 0.5 * Ds_IsoR03_sumPUPt);
	Float_t Ds_IsoR03_RelIso = (Ds_IsoR03_sumChargedHadronPt + std::max(0.0, Ds_IsoR03_sumNeutralHadronEt + Ds_IsoR03_sumPhotonEt - 0.5 * Ds_IsoR03_sumPUPt)) / Ds_mtm.perp();
	Float_t Ds_IsoR04_PFIso = Ds_IsoR04_sumChargedHadronPt + std::max(0.0, Ds_IsoR04_sumNeutralHadronEt + Ds_IsoR04_sumPhotonEt - 0.5 * Ds_IsoR04_sumPUPt);
	Float_t Ds_IsoR04_RelIso = (Ds_IsoR04_sumChargedHadronPt + std::max(0.0, Ds_IsoR04_sumNeutralHadronEt + Ds_IsoR04_sumPhotonEt - 0.5 * Ds_IsoR04_sumPUPt)) / Ds_mtm.perp();

	Float_t deltaR_Ds_Jet = 100.;
	Int_t num_jets_pt20 = 0;
	Int_t num_jets_pt30 = 0;
	Int_t num_jets_pt40 = 0;
	Int_t num_jets_pt50 = 0;
	Int_t num_jets_pt60 = 0;
	Int_t num_jets_pt70 = 0;
	Int_t num_jets_pt80 = 0;
	Int_t num_jets_pt90 = 0;
	Int_t num_jets_pt100 = 0;
	std::vector<Float_t> jet_pt;
	jet_pt.clear();
	for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
	{
		Int_t tightJetId = (*tightJetIdHandle)[edm::Ref<pat::JetCollection>(JetHandle, idx_jet)];
		if (tightJetId == 0)
			continue;

		const pat::Jet &Jet = (*JetHandle)[idx_jet];
		jet_pt.push_back(Jet.pt());

		if (Jet.pt() > 20)
			num_jets_pt20++;
		if (Jet.pt() > 30)
			num_jets_pt30++;
		if (Jet.pt() > 40)
			num_jets_pt40++;
		if (Jet.pt() > 50)
			num_jets_pt50++;
		if (Jet.pt() > 60)
			num_jets_pt60++;
		if (Jet.pt() > 70)
			num_jets_pt70++;
		if (Jet.pt() > 80)
			num_jets_pt80++;
		if (Jet.pt() > 90)
			num_jets_pt90++;
		if (Jet.pt() > 100)
			num_jets_pt100++;

		if(Jet.pt() > 30)
		{
			Float_t tmp_deltaR_Ds_Jet = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Jet.eta(), Jet.phi());
			if (tmp_deltaR_Ds_Jet < deltaR_Ds_Jet)
				deltaR_Ds_Jet = tmp_deltaR_Ds_Jet;
		}
	}


	//--------------------------------------------------//
	// muon candidate
	//--------------------------------------------------//
	std::pair<Bool_t, Measurement1D> mu_ipresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(muon.muonBestTrack()), primvtx);

	Float_t mu_IsoR03_sumChargedHadronPt = muon.pfIsolationR03().sumChargedHadronPt;
	Float_t mu_IsoR03_sumNeutralHadronEt = muon.pfIsolationR03().sumNeutralHadronEt;
	Float_t mu_IsoR03_sumPhotonEt = muon.pfIsolationR03().sumPhotonEt;
	Float_t mu_IsoR03_sumPUPt = muon.pfIsolationR03().sumPUPt;
	Float_t mu_IsoR03_PFIso = muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt);
	Float_t mu_IsoR03_RelIso = (muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt)) / muon.pt();

	Float_t mu_IsoR04_sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
	Float_t mu_IsoR04_sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
	Float_t mu_IsoR04_sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
	Float_t mu_IsoR04_sumPUPt = muon.pfIsolationR04().sumPUPt;
	Float_t mu_IsoR04_PFIso = muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt);
	Float_t mu_IsoR04_RelIso = (muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

	//==================================================//
	//--------------------------------------------------//
	// Storage
	//--------------------------------------------------//
	//==================================================//
	// refitted info from Kinematic Fit
	ftree->Kp_charge = best_Kp_pf.charge();
	ftree->Kp_eta = Kp_mtm.eta();
	ftree->Kp_phi = Kp_mtm.phi();
	ftree->Kp_px = Kp_mtm.x();
	ftree->Kp_py = Kp_mtm.y();
	ftree->Kp_pz = Kp_mtm.z();
	ftree->Kp_pt = Kp_mtm.perp();
	ftree->Kp_p = Kp_mtm.mag();
	ftree->Kp_mass = Kp->currentState().mass();

	ftree->Km_charge = best_Km_pf.charge();
	ftree->Km_eta = Km_mtm.eta();
	ftree->Km_phi = Km_mtm.phi();
	ftree->Km_px = Km_mtm.x();
	ftree->Km_py = Km_mtm.y();
	ftree->Km_pz = Km_mtm.z();
	ftree->Km_pt = Km_mtm.perp();
	ftree->Km_p = Km_mtm.mag();
	ftree->Km_mass = Km->currentState().mass();

	ftree->pi_charge = best_pi_pf.charge();
	ftree->pi_eta = pi_mtm.eta();
	ftree->pi_phi = pi_mtm.phi();
	ftree->pi_px = pi_mtm.x();
	ftree->pi_py = pi_mtm.y();
	ftree->pi_pz = pi_mtm.z();
	ftree->pi_pt = pi_mtm.perp();
	ftree->pi_p = pi_mtm.mag();
	ftree->pi_mass = pi->currentState().mass();

	ftree->phi_charge = best_Kp_pf.charge() + best_Km_pf.charge();
	ftree->phi_eta = phi_mtm.eta();
	ftree->phi_phi = phi_mtm.phi();
	ftree->phi_px = phi_mtm.x();
	ftree->phi_py = phi_mtm.y();
	ftree->phi_pz = phi_mtm.z();
	ftree->phi_pt = phi_mtm.perp();
	ftree->phi_p = phi_mtm.mag();
	ftree->phi_mass = phi_tlv.M();

	ftree->Ds_charge = best_Kp_pf.charge() + best_Km_pf.charge() + best_pi_pf.charge();
	ftree->Ds_eta = Ds_mtm.eta();
	ftree->Ds_phi = Ds_mtm.phi();
	ftree->Ds_px = Ds_mtm.x();
	ftree->Ds_py = Ds_mtm.y();
	ftree->Ds_pz = Ds_mtm.z();
	ftree->Ds_pt = Ds_mtm.perp();
	ftree->Ds_p = Ds_mtm.mag();
	ftree->Ds_mass = Ds->currentState().mass();

	// refitted vertex
	ftree->Ds_vx = Ds_vtx_pos.x();
	ftree->Ds_vy = Ds_vtx_pos.y();
	ftree->Ds_vz = Ds_vtx_pos.z();
	ftree->Ds_vxerr = Ds_vtx_cov.cxx();
	ftree->Ds_vyerr = Ds_vtx_cov.cyy();
	ftree->Ds_vzerr = Ds_vtx_cov.czz();
	ftree->Ds_chi2 = Ds_vertex->chiSquared();
	ftree->Ds_ndf = Ds_vertex->degreesOfFreedom();
	ftree->Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
	ftree->Ds_prob = TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom());

	// angular distance
	ftree->dR_Kp_Km = reco::deltaR(Kp_mtm, Km_mtm);
	ftree->dR_Kp_pi = reco::deltaR(Kp_mtm, pi_mtm);
	ftree->dR_Kp_phi = reco::deltaR(Kp_mtm, phi_mtm);
	ftree->dR_Kp_Ds = reco::deltaR(Kp_mtm, Ds_mtm);
	ftree->dR_Km_pi = reco::deltaR(Km_mtm, pi_mtm);
	ftree->dR_Km_phi = reco::deltaR(Km_mtm, phi_mtm);
	ftree->dR_Km_Ds = reco::deltaR(Km_mtm, Ds_mtm);
	ftree->dR_pi_phi = reco::deltaR(pi_mtm, phi_mtm);
	ftree->dR_pi_Ds = reco::deltaR(pi_mtm, Ds_mtm);
	ftree->dR_phi_Ds = reco::deltaR(phi_mtm, Ds_mtm);

	// IP
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
	if (Ds_ipresult.first)
	{
		ftree->Ds_ip = Ds_ipresult.second.value();
		ftree->Ds_iperr = Ds_ipresult.second.error();
		ftree->Ds_ipsig = Ds_ipresult.second.significance();
	}

	// PV-SV variable
	ftree->dira = displacement.dot(Ds_mtm) / (displacement.mag() * Ds_mtm.mag());

	ftree->Ds_FDxy = Ds_FD_XY.value();
	ftree->Ds_FDxyerr = Ds_FD_XY.error();
	ftree->Ds_FDxysig = Ds_FD_XY.significance();

	ftree->Ds_FD = Ds_FD_3D.value();
	ftree->Ds_FDerr = Ds_FD_3D.error();
	ftree->Ds_FDsig = Ds_FD_3D.significance();

	// Ds ioslation
	ftree->Ds_IsoR03_sumChargedHadronPt = Ds_IsoR03_sumChargedHadronPt;
	ftree->Ds_IsoR03_sumNeutralHadronEt = Ds_IsoR03_sumNeutralHadronEt;
	ftree->Ds_IsoR03_sumPhotonEt = Ds_IsoR03_sumPhotonEt;
	ftree->Ds_IsoR03_sumPUPt = Ds_IsoR03_sumPUPt;
	ftree->Ds_IsoR03_PFIso = Ds_IsoR03_PFIso;
	ftree->Ds_IsoR03_RelIso = Ds_IsoR03_RelIso;

	ftree->Ds_IsoR04_sumChargedHadronPt = Ds_IsoR04_sumChargedHadronPt;
	ftree->Ds_IsoR04_sumNeutralHadronEt = Ds_IsoR04_sumNeutralHadronEt;
	ftree->Ds_IsoR04_sumPhotonEt = Ds_IsoR04_sumPhotonEt;
	ftree->Ds_IsoR04_sumPUPt = Ds_IsoR04_sumPUPt;
	ftree->Ds_IsoR04_PFIso = Ds_IsoR04_PFIso;
	ftree->Ds_IsoR04_RelIso = Ds_IsoR04_RelIso;

	// muon kinematics
	ftree->mu_charge = muon.charge();
	ftree->mu_eta = muon.eta();
	ftree->mu_phi = muon.phi();
	ftree->mu_px = muon.px();
	ftree->mu_py = muon.py();
	ftree->mu_pz = muon.pz();
	ftree->mu_pt = muon.pt();
	ftree->mu_p = muon.p();
	ftree->mu_isHighPt = muon.isHighPtMuon(primvtx);
	ftree->mu_isLoose = muon.isLooseMuon();
	ftree->mu_isMedium = muon.isMediumMuon();
	ftree->mu_isSoft = muon.isSoftMuon(primvtx);
	ftree->mu_isTight = muon.isTightMuon(primvtx);
	ftree->mu_isPF = muon.isPFMuon();
	ftree->mu_isTracker = muon.isTrackerMuon();
	ftree->mu_isGlobal = muon.isGlobalMuon();

	if (mu_ipresult.first)
	{
		ftree->mu_ip = mu_ipresult.second.value();
		ftree->mu_iperr = mu_ipresult.second.error();
		ftree->mu_ipsig = mu_ipresult.second.significance();
	}

	ftree->dR_Kp_mu = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), muon.eta(), muon.phi());
	ftree->dR_Km_mu = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), muon.eta(), muon.phi());
	ftree->dR_pi_mu = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), muon.eta(), muon.phi());
	ftree->dR_phi_mu = reco::deltaR(phi_mtm.eta(), phi_mtm.phi(), muon.eta(), muon.phi());
	ftree->dR_Ds_mu = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), muon.eta(), muon.phi());

	// muon isolation
	ftree->mu_IsoR03_sumChargedHadronPt = mu_IsoR03_sumChargedHadronPt;
	ftree->mu_IsoR03_sumNeutralHadronEt = mu_IsoR03_sumNeutralHadronEt;
	ftree->mu_IsoR03_sumPhotonEt = mu_IsoR03_sumPhotonEt;
	ftree->mu_IsoR03_sumPUPt = mu_IsoR03_sumPUPt;
	ftree->mu_IsoR03_PFIso = mu_IsoR03_PFIso;
	ftree->mu_IsoR03_RelIso = mu_IsoR03_RelIso;

	ftree->mu_IsoR04_sumChargedHadronPt = mu_IsoR04_sumChargedHadronPt;
	ftree->mu_IsoR04_sumNeutralHadronEt = mu_IsoR04_sumNeutralHadronEt;
	ftree->mu_IsoR04_sumPhotonEt = mu_IsoR04_sumPhotonEt;
	ftree->mu_IsoR04_sumPUPt = mu_IsoR04_sumPUPt;
	ftree->mu_IsoR04_PFIso = mu_IsoR04_PFIso;
	ftree->mu_IsoR04_RelIso = mu_IsoR04_RelIso;

	// MET kinematics
	ftree->met_sig = Met.metSignificance();
	ftree->met_phi = Met.phi();
	ftree->met_px = Met.px();
	ftree->met_py = Met.py();
	ftree->met_pt = Met.pt();
	ftree->met_sumEt = Met.sumEt();

	// Jets
	ftree->num_jets = JetHandle->size();
	ftree->num_jets_pt20 = num_jets_pt20;
	ftree->num_jets_pt30 = num_jets_pt30;
	ftree->num_jets_pt40 = num_jets_pt40;
	ftree->num_jets_pt50 = num_jets_pt50;
	ftree->num_jets_pt60 = num_jets_pt60;
	ftree->num_jets_pt70 = num_jets_pt70;
	ftree->num_jets_pt80 = num_jets_pt80;
	ftree->num_jets_pt90 = num_jets_pt90;
	ftree->num_jets_pt100 = num_jets_pt100;
	ftree->jet_pt = jet_pt;
	ftree->deltaR_Ds_Jet = deltaR_Ds_Jet;
	
	ftree->Kp_match = (best_Kp_idx == Kp_match_idx);
	ftree->Km_match = (best_Km_idx == Km_match_idx);
	ftree->pi_match = (best_pi_idx == pi_match_idx);
	ftree->mu_match = (best_mu_idx == mu_match_idx);
	ftree->event_match = (best_Kp_idx == Kp_match_idx) && (best_Km_idx == Km_match_idx) && (best_pi_idx == pi_match_idx) && (best_mu_idx == mu_match_idx);

	ftree->tree->Fill();
}

Bool_t SignalReco::hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D SignalReco::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D SignalReco::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void SignalReco::endJob()
{
	InfoTree->Fill();
}

void SignalReco::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(SignalReco);

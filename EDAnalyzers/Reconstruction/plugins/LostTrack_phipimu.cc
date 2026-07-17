#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <string>

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "EDAnalyzers/Reconstruction/interface/LostTrack_phipimu_Tree.h"

class LostTrack_phipimu : public edm::one::EDAnalyzer<>
{
public:
	explicit LostTrack_phipimu(const edm::ParameterSet &);
	~LostTrack_phipimu() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const;
	Bool_t TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;
	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	// Charged-hadron pool cuts of H_phipimu: the idx_pf selection, plus the
	// per-leg pt cut applied inside the candidate loops (2 for K, 1 for pi).
	Bool_t PassPool(const pat::PackedCandidate &c, const Float_t ptmin) const
	{
		return c.hasTrackDetails() && c.trackHighPurity() && (abs(c.pdgId()) == 211) && (c.pt() > ptmin);
	}

	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> LostTrackToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<std::vector<reco::GenMET>> GenMETToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> ltAssociationToken_;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

	bool doTruth_;
	std::string channel_;

	const Float_t Mass_K    = 0.493677;
	const Float_t Mass_pi   = 0.13957039;
	const Float_t Mass_mu   = 0.1056583755;
	const Float_t Mass_phi   = 1.01946;
	const Float_t Mass_Ds   = 1.96835;
	const Float_t Mass_Dsst = 2.1066;

	Float_t Mass_K_sigma   = 0.000016;
	Float_t Mass_pi_sigma  = 0.000016;

	LostTrack_phipimu_Tree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsTriggered_;
	unsigned long nEventsLepton_;
	unsigned long nEventsGen_;
};

LostTrack_phipimu::LostTrack_phipimu(const edm::ParameterSet &iConfig)
{
	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag LostTrackTag_("lostTracks");
	LostTrackToken_ = consumes<pat::PackedCandidateCollection>(LostTrackTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

	edm::InputTag triggerTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(triggerTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag GenMETTag_("genMetTrue");
	GenMETToken_ = consumes<std::vector<reco::GenMET>>(GenMETTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::InputTag ltAssotiationTag_("lostTracksToGenAssociation");
	ltAssociationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(ltAssotiationTag_);

	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));

	doTruth_ = iConfig.getParameter<bool>("doTruth");
	channel_ = doTruth_ ? iConfig.getParameter<std::string>("channel") : "";
}

void LostTrack_phipimu::beginJob()
{
	edm::Service<TFileService> fs;
	ftree = new LostTrack_phipimu_Tree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();

	InfoTree = fs->make<TTree>("Info", "Info");
	InfoTree->Branch("nEventsProcessed_", &nEventsProcessed_);
	InfoTree->Branch("nEventsTriggered_", &nEventsTriggered_);
	InfoTree->Branch("nEventsLepton_", &nEventsLepton_);
	InfoTree->Branch("nEventsGen_", &nEventsGen_);

	nEventsProcessed_ = 0;
	nEventsTriggered_ = 0;
	nEventsLepton_ = 0;
	nEventsGen_ = 0;
}

void LostTrack_phipimu::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	ftree->Init();

	ftree->event_id = iEvent.id().event();
	ftree->run_id   = iEvent.id().run();
	ftree->lumi_id  = iEvent.luminosityBlock();

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::PackedCandidateCollection> LostTrackHandle;
	iEvent.getByToken(LostTrackToken_, LostTrackHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	Int_t Kp_match_idx    = -1;
	Int_t Km_match_idx    = -1;
	Int_t pi_match_idx    = -1;
	Int_t gamma_match_idx = -1;
	Int_t mu_match_idx    = -1;

	Int_t Kp_lt_match_idx = -1;
	Int_t Km_lt_match_idx = -1;
	Int_t pi_lt_match_idx = -1;

	//==================================================//
	// Gen (DsW channel)
	//==================================================//
	if (doTruth_ && channel_ == "DsW")
	{
		Int_t Gen_H_idx = 0, Gen_Ds_idx = 0, Gen_W_idx = 0, Gen_mu_idx = 0;
		Int_t Gen_Kp_idx = 0, Gen_Km_idx = 0, Gen_pi_idx = 0;
		Int_t num_Gen_H = 0, num_Gen_Ds = 0, num_Gen_W = 0, num_Gen_mu = 0, num_Gen_e = 0;
		Int_t num_Gen_Kp = 0, num_Gen_Km = 0, num_Gen_pi = 0;

		edm::Handle<reco::GenParticleCollection> GenPartHandle;
		iEvent.getByToken(GenPartToken_, GenPartHandle);

		for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
		{
			const reco::GenParticle &gp = (*GenPartHandle)[i];
			if (!gp.isHardProcess()) continue;
			if      (gp.pdgId() == 25)             { Gen_H_idx = i;  num_Gen_H++;   }
			else if (abs(gp.pdgId()) == 431)        { Gen_Ds_idx = i; num_Gen_Ds++;  }
			else if (abs(gp.pdgId()) == 24)         { Gen_W_idx = i;  num_Gen_W++;   }
			else if (abs(gp.pdgId()) == 13)         { Gen_mu_idx = i; num_Gen_mu++;  }
			else if (abs(gp.pdgId()) == 11)         { num_Gen_e++;                   }
		}
		if (!((num_Gen_H == 1) && (num_Gen_Ds == 1) && (num_Gen_W == 1) && (num_Gen_mu == 1) && (num_Gen_e == 0)))
			return;
		nEventsLepton_++;

		const reco::GenParticle &H_gp  = (*GenPartHandle)[Gen_H_idx];
		const reco::GenParticle &Ds_gp = (*GenPartHandle)[Gen_Ds_idx];
		const reco::GenParticle &W_gp  = (*GenPartHandle)[Gen_W_idx];
		const reco::GenParticle &mu_gp = (*GenPartHandle)[Gen_mu_idx];
		const Int_t Ds_charge = Ds_gp.charge();

		for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
		{
			const reco::GenParticle &gp = (*GenPartHandle)[i];
			if (gp.status() != 1) continue;
			if (gp.pdgId() == (Ds_charge * 321))
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ Gen_Kp_idx = i; num_Gen_Kp++; }
			}
			else if (gp.pdgId() == (-Ds_charge * 321))
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ Gen_Km_idx = i; num_Gen_Km++; }
			}
			else if (gp.pdgId() == (Ds_charge * 211))
			{
				if (TwoBodyDecay(gp, Ds_gp.pdgId(), 333, Ds_charge * 211) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ Gen_pi_idx = i; num_Gen_pi++; }
			}
		}
		if (!((num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1)))
			return;

		const reco::GenParticle &Kp_gp = (*GenPartHandle)[Gen_Kp_idx];
		const reco::GenParticle &Km_gp = (*GenPartHandle)[Gen_Km_idx];
		const reco::GenParticle &pi_gp = (*GenPartHandle)[Gen_pi_idx];
		nEventsGen_++;

		TLorentzVector Gen_Ds_tlv, Gen_mu_tlv, Gen_Dsmu_tlv;
		Gen_Ds_tlv.SetXYZM(Ds_gp.px(), Ds_gp.py(), Ds_gp.pz(), Mass_Ds);
		Gen_mu_tlv.SetXYZM(mu_gp.px(), mu_gp.py(), mu_gp.pz(), Mass_mu);
		Gen_Dsmu_tlv = Gen_Ds_tlv + Gen_mu_tlv;

		edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
		iEvent.getByToken(GenMETToken_, GenMETHandle);
		if (GenMETHandle->empty()) return;
		const reco::MET &GenMET = GenMETHandle->front();

		ftree->Gen_PV_vx = H_gp.vx();  ftree->Gen_PV_vy = H_gp.vy();  ftree->Gen_PV_vz = H_gp.vz();
		ftree->Gen_SV_vx = Kp_gp.vx(); ftree->Gen_SV_vy = Kp_gp.vy(); ftree->Gen_SV_vz = Kp_gp.vz();

		ftree->Gen_H_charge = H_gp.charge(); ftree->Gen_H_eta = H_gp.eta(); ftree->Gen_H_phi = H_gp.phi();
		ftree->Gen_H_px = H_gp.px(); ftree->Gen_H_py = H_gp.py(); ftree->Gen_H_pz = H_gp.pz();
		ftree->Gen_H_pt = H_gp.pt(); ftree->Gen_H_p = H_gp.p(); ftree->Gen_H_mass = H_gp.mass();

		ftree->Gen_Ds_charge = Ds_gp.charge(); ftree->Gen_Ds_eta = Ds_gp.eta(); ftree->Gen_Ds_phi = Ds_gp.phi();
		ftree->Gen_Ds_px = Ds_gp.px(); ftree->Gen_Ds_py = Ds_gp.py(); ftree->Gen_Ds_pz = Ds_gp.pz();
		ftree->Gen_Ds_pt = Ds_gp.pt(); ftree->Gen_Ds_p = Ds_gp.p(); ftree->Gen_Ds_mass = Ds_gp.mass();

		ftree->Gen_W_charge = W_gp.charge(); ftree->Gen_W_eta = W_gp.eta(); ftree->Gen_W_phi = W_gp.phi();
		ftree->Gen_W_px = W_gp.px(); ftree->Gen_W_py = W_gp.py(); ftree->Gen_W_pz = W_gp.pz();
		ftree->Gen_W_pt = W_gp.pt(); ftree->Gen_W_p = W_gp.p(); ftree->Gen_W_mass = W_gp.mass();

		ftree->Gen_Kp_charge = Kp_gp.charge(); ftree->Gen_Kp_eta = Kp_gp.eta(); ftree->Gen_Kp_phi = Kp_gp.phi();
		ftree->Gen_Kp_px = Kp_gp.px(); ftree->Gen_Kp_py = Kp_gp.py(); ftree->Gen_Kp_pz = Kp_gp.pz();
		ftree->Gen_Kp_pt = Kp_gp.pt(); ftree->Gen_Kp_p = Kp_gp.p(); ftree->Gen_Kp_mass = Kp_gp.mass();

		ftree->Gen_Km_charge = Km_gp.charge(); ftree->Gen_Km_eta = Km_gp.eta(); ftree->Gen_Km_phi = Km_gp.phi();
		ftree->Gen_Km_px = Km_gp.px(); ftree->Gen_Km_py = Km_gp.py(); ftree->Gen_Km_pz = Km_gp.pz();
		ftree->Gen_Km_pt = Km_gp.pt(); ftree->Gen_Km_p = Km_gp.p(); ftree->Gen_Km_mass = Km_gp.mass();

		ftree->Gen_pi_charge = pi_gp.charge(); ftree->Gen_pi_eta = pi_gp.eta(); ftree->Gen_pi_phi = pi_gp.phi();
		ftree->Gen_pi_px = pi_gp.px(); ftree->Gen_pi_py = pi_gp.py(); ftree->Gen_pi_pz = pi_gp.pz();
		ftree->Gen_pi_pt = pi_gp.pt(); ftree->Gen_pi_p = pi_gp.p(); ftree->Gen_pi_mass = pi_gp.mass();

		ftree->Gen_mu_charge = mu_gp.charge(); ftree->Gen_mu_eta = mu_gp.eta(); ftree->Gen_mu_phi = mu_gp.phi();
		ftree->Gen_mu_px = mu_gp.px(); ftree->Gen_mu_py = mu_gp.py(); ftree->Gen_mu_pz = mu_gp.pz();
		ftree->Gen_mu_pt = mu_gp.pt(); ftree->Gen_mu_p = mu_gp.p(); ftree->Gen_mu_mass = mu_gp.mass();

		ftree->Gen_Dsmu_charge = Ds_gp.charge() + mu_gp.charge();
		ftree->Gen_Dsmu_eta = Gen_Dsmu_tlv.Eta(); ftree->Gen_Dsmu_phi = Gen_Dsmu_tlv.Phi();
		ftree->Gen_Dsmu_px = Gen_Dsmu_tlv.Px(); ftree->Gen_Dsmu_py = Gen_Dsmu_tlv.Py(); ftree->Gen_Dsmu_pz = Gen_Dsmu_tlv.Pz();
		ftree->Gen_Dsmu_pt = Gen_Dsmu_tlv.Pt(); ftree->Gen_Dsmu_p = Gen_Dsmu_tlv.P(); ftree->Gen_Dsmu_mass = Gen_Dsmu_tlv.M();

		ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_gp, Km_gp); ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_gp, pi_gp);
		ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp); ftree->Gen_dR_Km_pi = reco::deltaR(Km_gp, pi_gp);
		ftree->Gen_dR_Km_Ds = reco::deltaR(Km_gp, Ds_gp); ftree->Gen_dR_pi_Ds = reco::deltaR(pi_gp, Ds_gp);
		ftree->Gen_dR_Kp_mu = reco::deltaR(Kp_gp, mu_gp); ftree->Gen_dR_Km_mu = reco::deltaR(Km_gp, mu_gp);
		ftree->Gen_dR_pi_mu = reco::deltaR(pi_gp, mu_gp); ftree->Gen_dR_Ds_mu = reco::deltaR(Ds_gp, mu_gp);

		ftree->Gen_met_phi = GenMET.phi(); ftree->Gen_met_px = GenMET.px(); ftree->Gen_met_py = GenMET.py();
		ftree->Gen_met_pt = GenMET.pt(); ftree->Gen_met_sumEt = GenMET.sumEt();

		// Gen match
		edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
		iEvent.getByToken(associationToken_, associationHandle);

		Float_t Kp_deltaR = 1.0, Km_deltaR = 1.0, pi_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
		{
			const pat::PackedCandidate &pf = (*PFCandHandle)[i];
			const pat::PackedCandidateRef pfRef(PFCandHandle, i);
			const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - pf.pt()) / gp.pt()) > 0.3) continue;

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_deltaR) { Kp_match_idx = i; Kp_deltaR = dR; } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_deltaR) { Km_match_idx = i; Km_deltaR = dR; } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_deltaR) { pi_match_idx = i; pi_deltaR = dR; } }
			}
		}

		// Gen match to lostTracks (same criteria, lostTracks collection)
		edm::Handle<edm::Association<reco::GenParticleCollection>> ltAssociationHandle;
		iEvent.getByToken(ltAssociationToken_, ltAssociationHandle);

		Float_t Kp_lt_deltaR = 1.0, Km_lt_deltaR = 1.0, pi_lt_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(LostTrackHandle->size()); i++)
		{
			const pat::PackedCandidate &lt = (*LostTrackHandle)[i];
			const pat::PackedCandidateRef ltRef(LostTrackHandle, i);
			const reco::GenParticleRef genRef = (*ltAssociationHandle)[ltRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), lt.eta(), lt.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - lt.pt()) / gp.pt()) > 0.3) continue;

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_lt_deltaR) { Kp_lt_match_idx = i; Kp_lt_deltaR = dR; } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_lt_deltaR) { Km_lt_match_idx = i; Km_lt_deltaR = dR; } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_lt_deltaR) { pi_lt_match_idx = i; pi_lt_deltaR = dR; } }
			}
		}

		Float_t mu_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(MuonHandle->size()); i++)
		{
			const pat::Muon &muon = (*MuonHandle)[i];
			if (!muon.genParticle()) continue;
			const reco::GenParticle &gp = *(muon.genParticle());
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), muon.eta(), muon.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - muon.pt()) / gp.pt()) > 0.3) continue;
			if (gp.pdgId() != mu_gp.pdgId()) continue;
			if (TwoBodyDecay(gp, W_gp.pdgId(), mu_gp.pdgId(), mu_gp.charge() * 14) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
			{ if (dR < mu_deltaR) { mu_match_idx = i; mu_deltaR = dR; } }
		}
	} // DsW gen

	//==================================================//
	// Gen (DsstW channel)
	//==================================================//
	else if (doTruth_ && channel_ == "DsstW")
	{
		Int_t Gen_H_idx = 0, Gen_Dsst_idx = 0, Gen_W_idx = 0, Gen_mu_idx = 0;
		Int_t Gen_Ds_idx = 0, Gen_gamma_idx = 0, Gen_Kp_idx = 0, Gen_Km_idx = 0, Gen_pi_idx = 0;
		Int_t num_Gen_H = 0, num_Gen_Dsst = 0, num_Gen_W = 0, num_Gen_mu = 0, num_Gen_e = 0;
		Int_t num_Gen_Ds = 0, num_Gen_gamma = 0, num_Gen_Kp = 0, num_Gen_Km = 0, num_Gen_pi = 0;

		edm::Handle<reco::GenParticleCollection> GenPartHandle;
		iEvent.getByToken(GenPartToken_, GenPartHandle);

		for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
		{
			const reco::GenParticle &gp = (*GenPartHandle)[i];
			if (!gp.isHardProcess()) continue;
			if      (gp.pdgId() == 25)             { Gen_H_idx = i;    num_Gen_H++;    }
			else if (abs(gp.pdgId()) == 433)        { Gen_Dsst_idx = i; num_Gen_Dsst++; }
			else if (abs(gp.pdgId()) == 24)         { Gen_W_idx = i;    num_Gen_W++;    }
			else if (abs(gp.pdgId()) == 13)         { Gen_mu_idx = i;   num_Gen_mu++;   }
			else if (abs(gp.pdgId()) == 11)         { num_Gen_e++;                      }
		}
		if (!((num_Gen_H == 1) && (num_Gen_Dsst == 1) && (num_Gen_W == 1) && (num_Gen_mu == 1) && (num_Gen_e == 0)))
			return;
		nEventsLepton_++;

		const reco::GenParticle &H_gp    = (*GenPartHandle)[Gen_H_idx];
		const reco::GenParticle &Dsst_gp = (*GenPartHandle)[Gen_Dsst_idx];
		const reco::GenParticle &W_gp    = (*GenPartHandle)[Gen_W_idx];
		const reco::GenParticle &mu_gp   = (*GenPartHandle)[Gen_mu_idx];
		const Int_t Dsst_charge = Dsst_gp.charge();

		for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
		{
			const reco::GenParticle &gp = (*GenPartHandle)[i];
			if ((gp.status() == 2) && gp.pdgId() == (Dsst_charge * 431))
			{
				if (TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()) && (gp.numberOfDaughters() > 1))
				{
					Bool_t d1 = false, d2 = false;
					for (size_t j = 0; j < gp.numberOfDaughters(); j++)
					{
						if (gp.daughter(j)->pdgId() == 333)               d1 = true;
						if (gp.daughter(j)->pdgId() == (Dsst_charge*211)) d2 = true;
					}
					if (d1 && d2) { Gen_Ds_idx = i; num_Gen_Ds++; }
				}
			}
			else if (gp.status() == 1)
			{
				if (gp.pdgId() == (Dsst_charge * 321))
				{
					if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					{ Gen_Kp_idx = i; num_Gen_Kp++; }
				}
				else if (gp.pdgId() == (-Dsst_charge * 321))
				{
					if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					{ Gen_Km_idx = i; num_Gen_Km++; }
				}
				else if (gp.pdgId() == (Dsst_charge * 211))
				{
					if (TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					{ Gen_pi_idx = i; num_Gen_pi++; }
				}
				else if (gp.pdgId() == 22)
				{
					if ((gp.mother(0)->pdgId() == Dsst_gp.pdgId()) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					{ Gen_gamma_idx = i; num_Gen_gamma++; }
				}
			}
		}
		if (!((num_Gen_Ds == 1) && (num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1) && (num_Gen_gamma == 1)))
			return;

		const reco::GenParticle &Ds_gp    = (*GenPartHandle)[Gen_Ds_idx];
		const reco::GenParticle &Kp_gp    = (*GenPartHandle)[Gen_Kp_idx];
		const reco::GenParticle &Km_gp    = (*GenPartHandle)[Gen_Km_idx];
		const reco::GenParticle &pi_gp    = (*GenPartHandle)[Gen_pi_idx];
		const reco::GenParticle &gamma_gp = (*GenPartHandle)[Gen_gamma_idx];
		nEventsGen_++;

		TLorentzVector Gen_Dsst_tlv, Gen_Ds_tlv, Gen_mu_tlv, Gen_Dsmu_tlv, Gen_Dsstmu_tlv;
		Gen_Dsst_tlv.SetXYZM(Dsst_gp.px(), Dsst_gp.py(), Dsst_gp.pz(), Mass_Dsst);
		Gen_Ds_tlv.SetXYZM(Ds_gp.px(), Ds_gp.py(), Ds_gp.pz(), Mass_Ds);
		Gen_mu_tlv.SetXYZM(mu_gp.px(), mu_gp.py(), mu_gp.pz(), Mass_mu);
		Gen_Dsmu_tlv  = Gen_Ds_tlv   + Gen_mu_tlv;
		Gen_Dsstmu_tlv = Gen_Dsst_tlv + Gen_mu_tlv;

		edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
		iEvent.getByToken(GenMETToken_, GenMETHandle);
		if (GenMETHandle->empty()) return;
		const reco::MET &GenMET = GenMETHandle->front();

		ftree->Gen_PV_vx = H_gp.vx();   ftree->Gen_PV_vy = H_gp.vy();   ftree->Gen_PV_vz = H_gp.vz();
		ftree->Gen_SV_vx = Kp_gp.vx();  ftree->Gen_SV_vy = Kp_gp.vy();  ftree->Gen_SV_vz = Kp_gp.vz();

		ftree->Gen_H_charge = H_gp.charge(); ftree->Gen_H_eta = H_gp.eta(); ftree->Gen_H_phi = H_gp.phi();
		ftree->Gen_H_px = H_gp.px(); ftree->Gen_H_py = H_gp.py(); ftree->Gen_H_pz = H_gp.pz();
		ftree->Gen_H_pt = H_gp.pt(); ftree->Gen_H_p = H_gp.p(); ftree->Gen_H_mass = H_gp.mass();

		ftree->Gen_Dsst_charge = Dsst_gp.charge(); ftree->Gen_Dsst_eta = Dsst_gp.eta(); ftree->Gen_Dsst_phi = Dsst_gp.phi();
		ftree->Gen_Dsst_px = Dsst_gp.px(); ftree->Gen_Dsst_py = Dsst_gp.py(); ftree->Gen_Dsst_pz = Dsst_gp.pz();
		ftree->Gen_Dsst_pt = Dsst_gp.pt(); ftree->Gen_Dsst_p = Dsst_gp.p(); ftree->Gen_Dsst_mass = Dsst_gp.mass();

		ftree->Gen_Ds_charge = Ds_gp.charge(); ftree->Gen_Ds_eta = Ds_gp.eta(); ftree->Gen_Ds_phi = Ds_gp.phi();
		ftree->Gen_Ds_px = Ds_gp.px(); ftree->Gen_Ds_py = Ds_gp.py(); ftree->Gen_Ds_pz = Ds_gp.pz();
		ftree->Gen_Ds_pt = Ds_gp.pt(); ftree->Gen_Ds_p = Ds_gp.p(); ftree->Gen_Ds_mass = Ds_gp.mass();

		ftree->Gen_gamma_charge = gamma_gp.charge(); ftree->Gen_gamma_eta = gamma_gp.eta(); ftree->Gen_gamma_phi = gamma_gp.phi();
		ftree->Gen_gamma_px = gamma_gp.px(); ftree->Gen_gamma_py = gamma_gp.py(); ftree->Gen_gamma_pz = gamma_gp.pz();
		ftree->Gen_gamma_pt = gamma_gp.pt(); ftree->Gen_gamma_p = gamma_gp.p(); ftree->Gen_gamma_mass = gamma_gp.mass();

		ftree->Gen_W_charge = W_gp.charge(); ftree->Gen_W_eta = W_gp.eta(); ftree->Gen_W_phi = W_gp.phi();
		ftree->Gen_W_px = W_gp.px(); ftree->Gen_W_py = W_gp.py(); ftree->Gen_W_pz = W_gp.pz();
		ftree->Gen_W_pt = W_gp.pt(); ftree->Gen_W_p = W_gp.p(); ftree->Gen_W_mass = W_gp.mass();

		ftree->Gen_Kp_charge = Kp_gp.charge(); ftree->Gen_Kp_eta = Kp_gp.eta(); ftree->Gen_Kp_phi = Kp_gp.phi();
		ftree->Gen_Kp_px = Kp_gp.px(); ftree->Gen_Kp_py = Kp_gp.py(); ftree->Gen_Kp_pz = Kp_gp.pz();
		ftree->Gen_Kp_pt = Kp_gp.pt(); ftree->Gen_Kp_p = Kp_gp.p(); ftree->Gen_Kp_mass = Kp_gp.mass();

		ftree->Gen_Km_charge = Km_gp.charge(); ftree->Gen_Km_eta = Km_gp.eta(); ftree->Gen_Km_phi = Km_gp.phi();
		ftree->Gen_Km_px = Km_gp.px(); ftree->Gen_Km_py = Km_gp.py(); ftree->Gen_Km_pz = Km_gp.pz();
		ftree->Gen_Km_pt = Km_gp.pt(); ftree->Gen_Km_p = Km_gp.p(); ftree->Gen_Km_mass = Km_gp.mass();

		ftree->Gen_pi_charge = pi_gp.charge(); ftree->Gen_pi_eta = pi_gp.eta(); ftree->Gen_pi_phi = pi_gp.phi();
		ftree->Gen_pi_px = pi_gp.px(); ftree->Gen_pi_py = pi_gp.py(); ftree->Gen_pi_pz = pi_gp.pz();
		ftree->Gen_pi_pt = pi_gp.pt(); ftree->Gen_pi_p = pi_gp.p(); ftree->Gen_pi_mass = pi_gp.mass();

		ftree->Gen_mu_charge = mu_gp.charge(); ftree->Gen_mu_eta = mu_gp.eta(); ftree->Gen_mu_phi = mu_gp.phi();
		ftree->Gen_mu_px = mu_gp.px(); ftree->Gen_mu_py = mu_gp.py(); ftree->Gen_mu_pz = mu_gp.pz();
		ftree->Gen_mu_pt = mu_gp.pt(); ftree->Gen_mu_p = mu_gp.p(); ftree->Gen_mu_mass = mu_gp.mass();

		ftree->Gen_Dsmu_charge = Ds_gp.charge() + mu_gp.charge();
		ftree->Gen_Dsmu_eta = Gen_Dsmu_tlv.Eta(); ftree->Gen_Dsmu_phi = Gen_Dsmu_tlv.Phi();
		ftree->Gen_Dsmu_px = Gen_Dsmu_tlv.Px(); ftree->Gen_Dsmu_py = Gen_Dsmu_tlv.Py(); ftree->Gen_Dsmu_pz = Gen_Dsmu_tlv.Pz();
		ftree->Gen_Dsmu_pt = Gen_Dsmu_tlv.Pt(); ftree->Gen_Dsmu_p = Gen_Dsmu_tlv.P(); ftree->Gen_Dsmu_mass = Gen_Dsmu_tlv.M();

		ftree->Gen_Dsstmu_charge = Dsst_gp.charge() + mu_gp.charge();
		ftree->Gen_Dsstmu_eta = Gen_Dsstmu_tlv.Eta(); ftree->Gen_Dsstmu_phi = Gen_Dsstmu_tlv.Phi();
		ftree->Gen_Dsstmu_px = Gen_Dsstmu_tlv.Px(); ftree->Gen_Dsstmu_py = Gen_Dsstmu_tlv.Py(); ftree->Gen_Dsstmu_pz = Gen_Dsstmu_tlv.Pz();
		ftree->Gen_Dsstmu_pt = Gen_Dsstmu_tlv.Pt(); ftree->Gen_Dsstmu_p = Gen_Dsstmu_tlv.P(); ftree->Gen_Dsstmu_mass = Gen_Dsstmu_tlv.M();

		ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_gp, Km_gp);     ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_gp, pi_gp);
		ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp);     ftree->Gen_dR_Km_pi = reco::deltaR(Km_gp, pi_gp);
		ftree->Gen_dR_Km_Ds = reco::deltaR(Km_gp, Ds_gp);     ftree->Gen_dR_pi_Ds = reco::deltaR(pi_gp, Ds_gp);
		ftree->Gen_dR_Kp_gamma = reco::deltaR(Kp_gp, gamma_gp);   ftree->Gen_dR_Kp_Dsst = reco::deltaR(Kp_gp, Dsst_gp);
		ftree->Gen_dR_Km_gamma = reco::deltaR(Km_gp, gamma_gp);   ftree->Gen_dR_Km_Dsst = reco::deltaR(Km_gp, Dsst_gp);
		ftree->Gen_dR_pi_gamma = reco::deltaR(pi_gp, gamma_gp);   ftree->Gen_dR_pi_Dsst = reco::deltaR(pi_gp, Dsst_gp);
		ftree->Gen_dR_Ds_gamma = reco::deltaR(Ds_gp, gamma_gp);   ftree->Gen_dR_Ds_Dsst = reco::deltaR(Ds_gp, Dsst_gp);
		ftree->Gen_dR_gamma_Dsst = reco::deltaR(gamma_gp, Dsst_gp);
		ftree->Gen_dR_Kp_mu = reco::deltaR(Kp_gp, mu_gp);     ftree->Gen_dR_Km_mu = reco::deltaR(Km_gp, mu_gp);
		ftree->Gen_dR_pi_mu = reco::deltaR(pi_gp, mu_gp);     ftree->Gen_dR_Ds_mu = reco::deltaR(Ds_gp, mu_gp);
		ftree->Gen_dR_gamma_mu = reco::deltaR(gamma_gp, mu_gp); ftree->Gen_dR_Dsst_mu = reco::deltaR(Dsst_gp, mu_gp);

		ftree->Gen_met_phi = GenMET.phi(); ftree->Gen_met_px = GenMET.px(); ftree->Gen_met_py = GenMET.py();
		ftree->Gen_met_pt = GenMET.pt(); ftree->Gen_met_sumEt = GenMET.sumEt();

		// Gen match
		edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
		iEvent.getByToken(associationToken_, associationHandle);

		Float_t Kp_deltaR = 1.0, Km_deltaR = 1.0, pi_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
		{
			const pat::PackedCandidate &pf = (*PFCandHandle)[i];
			const pat::PackedCandidateRef pfRef(PFCandHandle, i);
			const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - pf.pt()) / gp.pt()) > 0.3) continue;

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_deltaR) { Kp_match_idx = i; Kp_deltaR = dR; } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_deltaR) { Km_match_idx = i; Km_deltaR = dR; } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_deltaR) { pi_match_idx = i; pi_deltaR = dR; } }
			}
		}

		// Gen match to lostTracks (same criteria, lostTracks collection)
		edm::Handle<edm::Association<reco::GenParticleCollection>> ltAssociationHandle;
		iEvent.getByToken(ltAssociationToken_, ltAssociationHandle);

		Float_t Kp_lt_deltaR = 1.0, Km_lt_deltaR = 1.0, pi_lt_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(LostTrackHandle->size()); i++)
		{
			const pat::PackedCandidate &lt = (*LostTrackHandle)[i];
			const pat::PackedCandidateRef ltRef(LostTrackHandle, i);
			const reco::GenParticleRef genRef = (*ltAssociationHandle)[ltRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), lt.eta(), lt.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - lt.pt()) / gp.pt()) > 0.3) continue;

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_lt_deltaR) { Kp_lt_match_idx = i; Kp_lt_deltaR = dR; } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_lt_deltaR) { Km_lt_match_idx = i; Km_lt_deltaR = dR; } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_lt_deltaR) { pi_lt_match_idx = i; pi_lt_deltaR = dR; } }
			}
		}

		Float_t gamma_deltaR = 0.02;
		for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
		{
			const pat::PackedCandidate &pf = (*PFCandHandle)[i];
			if (pf.pdgId() != 22) continue;
			Float_t dR = reco::deltaR(pf.eta(), pf.phi(), gamma_gp.eta(), gamma_gp.phi());
			if (std::abs((gamma_gp.pt() - pf.pt()) / gamma_gp.pt()) > 0.5) continue;
			if (dR < gamma_deltaR) { gamma_match_idx = i; gamma_deltaR = dR; }
		}

		Float_t mu_deltaR = 1.0;
		for (Int_t i = 0; i < static_cast<Int_t>(MuonHandle->size()); i++)
		{
			const pat::Muon &muon = (*MuonHandle)[i];
			if (!muon.genParticle()) continue;
			const reco::GenParticle &gp = *(muon.genParticle());
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), muon.eta(), muon.phi());
			if (dR > 0.05) continue;
			if (std::abs((gp.pt() - muon.pt()) / gp.pt()) > 0.3) continue;
			if (gp.pdgId() != mu_gp.pdgId()) continue;
			if (TwoBodyDecay(gp, W_gp.pdgId(), mu_gp.pdgId(), mu_gp.charge() * 14) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
			{ if (dR < mu_deltaR) { mu_match_idx = i; mu_deltaR = dR; } }
		}
	} // DsstW gen

	//==================================================//
	// Triggers
	//==================================================//
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	Bool_t trig_HLT_IsoMu24 = false;
	Bool_t trig_HLT_Mu50    = false;
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);
	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		Bool_t pass = triggerHandle->accept(i);
		if (trigName.Contains("HLT_IsoMu24_v")) trig_HLT_IsoMu24 = pass;
		else if (trigName.Contains("HLT_Mu50_v")) trig_HLT_Mu50 = pass;
	}
	if (!trig_HLT_IsoMu24 && !trig_HLT_Mu50)
		return;
	nEventsTriggered_++;

	ftree->trig_HLT_IsoMu24 = trig_HLT_IsoMu24;
	ftree->trig_HLT_Mu50    = trig_HLT_Mu50;

	//==================================================//
	// Primary vertex
	//==================================================//
	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);
	if (primvtxHandle->empty()) return;
	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid()) return;
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_cov(primvtx.covariance());

	ftree->primvtx_isValid    = primvtx.isValid();
	ftree->primvtx_isFake     = primvtx.isFake();
	ftree->primvtx_chi2       = primvtx.chi2();
	ftree->primvtx_ndof       = primvtx.ndof();
	ftree->primvtx_chi2ndof   = primvtx.chi2() / primvtx.ndof();
	ftree->primvtx_vx         = primvtx.x();
	ftree->primvtx_vy         = primvtx.y();
	ftree->primvtx_vz         = primvtx.z();
	ftree->primvtx_vxerr      = primvtx.xError();
	ftree->primvtx_vyerr      = primvtx.yError();
	ftree->primvtx_vzerr      = primvtx.zError();

	//==================================================//
	// Muon
	//==================================================//
	Int_t best_mu_idx = -1;
	Float_t best_mu_pt = -1;

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (muon.pt() < 24) continue;
		if (!muon.isTightMuon(primvtx)) continue;
		if (muon.pt() > best_mu_pt) { best_mu_idx = idx_mu; best_mu_pt = muon.pt(); }
	}
	if (best_mu_idx == -1) return;
	const pat::Muon &mu = (*MuonHandle)[best_mu_idx];
	TLorentzVector mu_tlv;
	mu_tlv.SetXYZM(mu.px(), mu.py(), mu.pz(), Mass_mu);

	std::pair<Bool_t, Measurement1D> mu_ipresult   = IPTools::absoluteImpactParameter3D((*ttBuilder).build(mu.muonBestTrack()), primvtx);
	std::pair<Bool_t, Measurement1D> mu_ipxyresult = IPTools::absoluteTransverseImpactParameter((*ttBuilder).build(mu.muonBestTrack()), primvtx);

	Float_t mu_IsoR03_PFIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5f * mu.pfIsolationR03().sumPUPt);
	Float_t mu_IsoR04_PFIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5f * mu.pfIsolationR04().sumPUPt);

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (muon.charge() == mu.charge()  && muon.isTightMuon(primvtx)) ftree->num_rs_mu++;
		else if (muon.charge() == -mu.charge() && muon.isTightMuon(primvtx)) ftree->num_ws_mu++;
	}

	ftree->mu_charge   = mu.charge();
	ftree->mu_eta      = mu.eta();   ftree->mu_phi  = mu.phi();
	ftree->mu_px       = mu.px();    ftree->mu_py   = mu.py();   ftree->mu_pz = mu.pz();
	ftree->mu_pt       = mu.pt();    ftree->mu_p    = mu.p();
	ftree->mu_isHighPt = mu.isHighPtMuon(primvtx);
	ftree->mu_isLoose  = mu.isLooseMuon();
	ftree->mu_isMedium = mu.isMediumMuon();
	ftree->mu_isSoft   = mu.isSoftMuon(primvtx);
	ftree->mu_isTight  = mu.isTightMuon(primvtx);
	ftree->mu_isPF     = mu.isPFMuon();
	ftree->mu_isTracker = mu.isTrackerMuon();
	ftree->mu_isGlobal  = mu.isGlobalMuon();
	if (mu_ipresult.first)
	{
		ftree->mu_ip    = mu_ipresult.second.value();
		ftree->mu_iperr = mu_ipresult.second.error();
		ftree->mu_ipsig = mu_ipresult.second.significance();
	}
	if (mu_ipxyresult.first)
	{
		ftree->mu_ipxy    = mu_ipxyresult.second.value();
		ftree->mu_ipxyerr = mu_ipxyresult.second.error();
		ftree->mu_ipxysig = mu_ipxyresult.second.significance();
	}
	ftree->mu_IsoR03_sumChargedHadronPt = mu.pfIsolationR03().sumChargedHadronPt;
	ftree->mu_IsoR03_sumNeutralHadronEt = mu.pfIsolationR03().sumNeutralHadronEt;
	ftree->mu_IsoR03_sumPhotonEt        = mu.pfIsolationR03().sumPhotonEt;
	ftree->mu_IsoR03_sumPUPt            = mu.pfIsolationR03().sumPUPt;
	ftree->mu_IsoR03_PFIso              = mu_IsoR03_PFIso;
	ftree->mu_IsoR03_RelIso             = mu_IsoR03_PFIso / mu.pt();
	ftree->mu_IsoR04_sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
	ftree->mu_IsoR04_sumNeutralHadronEt = mu.pfIsolationR04().sumNeutralHadronEt;
	ftree->mu_IsoR04_sumPhotonEt        = mu.pfIsolationR04().sumPhotonEt;
	ftree->mu_IsoR04_sumPUPt            = mu.pfIsolationR04().sumPUPt;
	ftree->mu_IsoR04_PFIso              = mu_IsoR04_PFIso;
	ftree->mu_IsoR04_RelIso             = mu_IsoR04_PFIso / mu.pt();
	ftree->mu_match = (best_mu_idx == mu_match_idx);

	//==================================================//
	// lostTrack recovery
	//==================================================//
	ftree->Kp_pf_match    = (Kp_match_idx >= 0);
	ftree->Km_pf_match    = (Km_match_idx >= 0);
	ftree->pi_pf_match    = (pi_match_idx >= 0);
	ftree->gamma_pf_match = (gamma_match_idx >= 0);

	ftree->Kp_lt_match = (Kp_lt_match_idx >= 0);
	ftree->Km_lt_match = (Km_lt_match_idx >= 0);
	ftree->pi_lt_match = (pi_lt_match_idx >= 0);

	if (Kp_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Kp_match_idx];
		ftree->Kp_pf_pool = PassPool(pf, 2);
		ftree->Kp_pf_dR  = reco::deltaR(ftree->Gen_Kp_eta, ftree->Gen_Kp_phi, pf.eta(), pf.phi());
		ftree->Kp_pf_dpt = (pf.pt() - ftree->Gen_Kp_pt) / ftree->Gen_Kp_pt;
	}
	if (Km_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Km_match_idx];
		ftree->Km_pf_pool = PassPool(pf, 2);
		ftree->Km_pf_dR  = reco::deltaR(ftree->Gen_Km_eta, ftree->Gen_Km_phi, pf.eta(), pf.phi());
		ftree->Km_pf_dpt = (pf.pt() - ftree->Gen_Km_pt) / ftree->Gen_Km_pt;
	}
	if (pi_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[pi_match_idx];
		ftree->pi_pf_pool = PassPool(pf, 1);
		ftree->pi_pf_dR  = reco::deltaR(ftree->Gen_pi_eta, ftree->Gen_pi_phi, pf.eta(), pf.phi());
		ftree->pi_pf_dpt = (pf.pt() - ftree->Gen_pi_pt) / ftree->Gen_pi_pt;
	}

	if (Kp_lt_match_idx >= 0)
	{
		const pat::PackedCandidate &lt = (*LostTrackHandle)[Kp_lt_match_idx];
		ftree->Kp_lt_pool = PassPool(lt, 2);
		ftree->Kp_lt_pt = lt.pt(); ftree->Kp_lt_eta = lt.eta();
		ftree->Kp_lt_dR  = reco::deltaR(ftree->Gen_Kp_eta, ftree->Gen_Kp_phi, lt.eta(), lt.phi());
		ftree->Kp_lt_dpt = (lt.pt() - ftree->Gen_Kp_pt) / ftree->Gen_Kp_pt;
	}
	if (Km_lt_match_idx >= 0)
	{
		const pat::PackedCandidate &lt = (*LostTrackHandle)[Km_lt_match_idx];
		ftree->Km_lt_pool = PassPool(lt, 2);
		ftree->Km_lt_pt = lt.pt(); ftree->Km_lt_eta = lt.eta();
		ftree->Km_lt_dR  = reco::deltaR(ftree->Gen_Km_eta, ftree->Gen_Km_phi, lt.eta(), lt.phi());
		ftree->Km_lt_dpt = (lt.pt() - ftree->Gen_Km_pt) / ftree->Gen_Km_pt;
	}
	if (pi_lt_match_idx >= 0)
	{
		const pat::PackedCandidate &lt = (*LostTrackHandle)[pi_lt_match_idx];
		ftree->pi_lt_pool = PassPool(lt, 1);
		ftree->pi_lt_pt = lt.pt(); ftree->pi_lt_eta = lt.eta();
		ftree->pi_lt_dR  = reco::deltaR(ftree->Gen_pi_eta, ftree->Gen_pi_phi, lt.eta(), lt.phi());
		ftree->pi_lt_dpt = (lt.pt() - ftree->Gen_pi_pt) / ftree->Gen_pi_pt;
	}

	ftree->tree->Fill();
}
void LostTrack_phipimu::endJob()
{
	InfoTree->Fill();
}

Bool_t LostTrack_phipimu::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
{
	if (gp.numberOfMothers() == 0) return false;
	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother) continue;
		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 3)
		{
			Bool_t f1 = false, f2 = false, f3 = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id) f1 = true;
				if (mother->daughter(j)->pdgId() == daughter2id) f2 = true;
				if (mother->daughter(j)->pdgId() == daughter3id) f3 = true;
			}
			if (f1 && f2 && f3) return true;
		}
		if (ThreeBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id, daughter3id)) return true;
	}
	return false;
}

Bool_t LostTrack_phipimu::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
{
	if (gp.numberOfMothers() == 0) return false;
	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother) continue;
		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 2)
		{
			Bool_t f1 = false, f2 = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id) f1 = true;
				if (mother->daughter(j)->pdgId() == daughter2id) f2 = true;
			}
			if (f1 && f2) return true;
		}
		if (TwoBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id)) return true;
	}
	return false;
}

Measurement1D LostTrack_phipimu::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x(); vDiff[1] = diff.y(); vDiff[2] = 0.;
	Double_t dist = sqrt(pow(diff.x(), 2) + pow(diff.y(), 2));
	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = (dist != 0) ? sqrt(err2) / dist : 0;
	return Measurement1D(dist, err);
}

Measurement1D LostTrack_phipimu::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x(); vDiff[1] = diff.y(); vDiff[2] = diff.z();
	Double_t dist = diff.mag();
	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = (dist != 0) ? sqrt(err2) / dist : 0;
	return Measurement1D(dist, err);
}

void LostTrack_phipimu::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.add<bool>("doTruth", true);
	desc.add<std::string>("channel", "DsW");
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LostTrack_phipimu);

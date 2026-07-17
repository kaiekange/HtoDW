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

#include "EDAnalyzers/Reconstruction/interface/GeoMatch_phipimu_Tree.h"

class GeoMatch_phipimu : public edm::one::EDAnalyzer<>
{
public:
	explicit GeoMatch_phipimu(const edm::ParameterSet &);
	~GeoMatch_phipimu() override {};

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
	// per-daughter pt cut applied inside the candidate loops (2 for K, 1 for pi).
	Bool_t PassPool(const pat::PackedCandidate &c, const Float_t ptmin) const
	{
		return c.hasTrackDetails() && c.trackHighPurity() && (abs(c.pdgId()) == 211) && (c.pt() > ptmin);
	}

	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<std::vector<reco::GenMET>> GenMETToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
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

	// geometric matching criteria. Same values association is judged
	// at, so the two are compared at the same working point.
	const Float_t geo_dR_cut  = 0.05;
	const Float_t geo_dpt_cut = 0.3;

	GeoMatch_phipimu_Tree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsTriggered_;
	unsigned long nEventsLepton_;
	unsigned long nEventsGen_;
};

GeoMatch_phipimu::GeoMatch_phipimu(const edm::ParameterSet &iConfig)
{
	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

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

	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));

	doTruth_ = iConfig.getParameter<bool>("doTruth");
	channel_ = doTruth_ ? iConfig.getParameter<std::string>("channel") : "";
}

void GeoMatch_phipimu::beginJob()
{
	edm::Service<TFileService> fs;
	ftree = new GeoMatch_phipimu_Tree(fs->make<TTree>("Events", "Events"));
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

void GeoMatch_phipimu::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
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

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	// The cross-tag below fills outside the per-channel gen blocks, so it needs
	// its own handle on the association map.
	edm::Handle<edm::Association<reco::GenParticleCollection>> xassocHandle;
	if (doTruth_) iEvent.getByToken(associationToken_, xassocHandle);

	Int_t Kp_match_idx    = -1;
	Int_t Km_match_idx    = -1;
	Int_t pi_match_idx    = -1;
	Int_t gamma_match_idx = -1;
	Int_t mu_match_idx    = -1;

	// genRef.key() of association's pick, for the offline compare
	// against Gen_<dau>_idx that the cross-tag rests on.
	Int_t Kp_match_key = -1;
	Int_t Km_match_key = -1;
	Int_t pi_match_key = -1;

	Int_t Kp_geo_match_idx = -1;
	Int_t Km_geo_match_idx = -1;
	Int_t pi_geo_match_idx = -1;

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

		ftree->Gen_Kp_idx = Gen_Kp_idx; ftree->Gen_Km_idx = Gen_Km_idx; ftree->Gen_pi_idx = Gen_pi_idx;

		// Gen match: association. Which gen particle a candidate is stays
		// decided by the map plus the ancestry chain, exactly as in H_phipimu.
		// The dR < 0.05 / |dpt| < 0.3 gate is gone: the map, not geometry,
		// identifies the particle, so the gate is applied offline instead and the
		// working point stays scannable. The dR seeds are 999 rather than 1.0 so
		// nothing is truncated implicitly now that the gate is off.
		edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
		iEvent.getByToken(associationToken_, associationHandle);

		Float_t Kp_deltaR = 999.0, Km_deltaR = 999.0, pi_deltaR = 999.0;
		for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
		{
			const pat::PackedCandidate &pf = (*PFCandHandle)[i];
			const pat::PackedCandidateRef pfRef(PFCandHandle, i);
			const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_deltaR) { Kp_match_idx = i; Kp_deltaR = dR; Kp_match_key = static_cast<Int_t>(genRef.key()); } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, Kp_gp.pdgId(), Km_gp.pdgId()) && TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_deltaR) { Km_match_idx = i; Km_deltaR = dR; Km_match_key = static_cast<Int_t>(genRef.key()); } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Ds_gp.pdgId(), 333, pi_gp.pdgId()) && TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_deltaR) { pi_match_idx = i; pi_deltaR = dR; pi_match_key = static_cast<Int_t>(genRef.key()); } }
			}
		}

		// Gen match: geometric. Gen-driven, no association map: for each gen
		// daughter, keep every charged-hadron PF of that daughter's charge which
		// passes BOTH dR < 0.05 and |dpt| < 0.3, then take the smallest dR among
		// the survivors. Same shape as the gamma matching above: the dpt cut sits
		// INSIDE the loop, so a soft track that happens to be nearest cannot block
		// a good candidate slightly further out. The three tests are deliberately
		// not else-if: K+ and pi carry the same charge, so one PF must be offered
		// to both independently. A collision shows up as equal indices and is
		// resolved offline.
		{
			Float_t Kp_geo_deltaR = geo_dR_cut, Km_geo_deltaR = geo_dR_cut, pi_geo_deltaR = geo_dR_cut;
			for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
			{
				const pat::PackedCandidate &pf = (*PFCandHandle)[i];
				if (abs(pf.pdgId()) != 211) continue;   // PF calls every charged hadron +-211
				if (pf.charge() == Kp_gp.charge())
				{
					Float_t dR = reco::deltaR(Kp_gp.eta(), Kp_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((Kp_gp.pt() - pf.pt()) / Kp_gp.pt()) <= geo_dpt_cut)
						if (dR < Kp_geo_deltaR) { Kp_geo_match_idx = i; Kp_geo_deltaR = dR; }
				}
				if (pf.charge() == Km_gp.charge())
				{
					Float_t dR = reco::deltaR(Km_gp.eta(), Km_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((Km_gp.pt() - pf.pt()) / Km_gp.pt()) <= geo_dpt_cut)
						if (dR < Km_geo_deltaR) { Km_geo_match_idx = i; Km_geo_deltaR = dR; }
				}
				if (pf.charge() == pi_gp.charge())
				{
					Float_t dR = reco::deltaR(pi_gp.eta(), pi_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((pi_gp.pt() - pf.pt()) / pi_gp.pt()) <= geo_dpt_cut)
						if (dR < pi_geo_deltaR) { pi_geo_match_idx = i; pi_geo_deltaR = dR; }
				}
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

		ftree->Gen_Kp_idx = Gen_Kp_idx; ftree->Gen_Km_idx = Gen_Km_idx; ftree->Gen_pi_idx = Gen_pi_idx;

		// Gen match: association. See the DsW block for why the dR / dpt gate
		// is gone and the key is recorded; only the ancestry chain differs here.
		edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
		iEvent.getByToken(associationToken_, associationHandle);

		Float_t Kp_deltaR = 999.0, Km_deltaR = 999.0, pi_deltaR = 999.0;
		for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
		{
			const pat::PackedCandidate &pf = (*PFCandHandle)[i];
			const pat::PackedCandidateRef pfRef(PFCandHandle, i);
			const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
			if (genRef.isNull()) continue;
			const reco::GenParticle &gp = *genRef;
			Float_t dR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());

			if (gp.pdgId() == Kp_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Kp_deltaR) { Kp_match_idx = i; Kp_deltaR = dR; Kp_match_key = static_cast<Int_t>(genRef.key()); } }
			}
			else if (gp.pdgId() == Km_gp.pdgId())
			{
				if (TwoBodyDecay(gp, 333, 321, -321) && TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < Km_deltaR) { Km_match_idx = i; Km_deltaR = dR; Km_match_key = static_cast<Int_t>(genRef.key()); } }
			}
			else if (gp.pdgId() == pi_gp.pdgId())
			{
				if (TwoBodyDecay(gp, Dsst_charge * 431, 333, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{ if (dR < pi_deltaR) { pi_match_idx = i; pi_deltaR = dR; pi_match_key = static_cast<Int_t>(genRef.key()); } }
			}
		}

		// Gen match: geometric. Identical to the DsW block: both cuts inside
		// the loop, then the smallest dR among the survivors.
		{
			Float_t Kp_geo_deltaR = geo_dR_cut, Km_geo_deltaR = geo_dR_cut, pi_geo_deltaR = geo_dR_cut;
			for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
			{
				const pat::PackedCandidate &pf = (*PFCandHandle)[i];
				if (abs(pf.pdgId()) != 211) continue;   // PF calls every charged hadron +-211
				if (pf.charge() == Kp_gp.charge())
				{
					Float_t dR = reco::deltaR(Kp_gp.eta(), Kp_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((Kp_gp.pt() - pf.pt()) / Kp_gp.pt()) <= geo_dpt_cut)
						if (dR < Kp_geo_deltaR) { Kp_geo_match_idx = i; Kp_geo_deltaR = dR; }
				}
				if (pf.charge() == Km_gp.charge())
				{
					Float_t dR = reco::deltaR(Km_gp.eta(), Km_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((Km_gp.pt() - pf.pt()) / Km_gp.pt()) <= geo_dpt_cut)
						if (dR < Km_geo_deltaR) { Km_geo_match_idx = i; Km_geo_deltaR = dR; }
				}
				if (pf.charge() == pi_gp.charge())
				{
					Float_t dR = reco::deltaR(pi_gp.eta(), pi_gp.phi(), pf.eta(), pf.phi());
					if (std::abs((pi_gp.pt() - pf.pt()) / pi_gp.pt()) <= geo_dpt_cut)
						if (dR < pi_geo_deltaR) { pi_geo_match_idx = i; pi_geo_deltaR = dR; }
				}
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
	// association vs geometric truth-matching
	//==================================================//
	ftree->Kp_pf_match    = (Kp_match_idx >= 0);
	ftree->Km_pf_match    = (Km_match_idx >= 0);
	ftree->pi_pf_match    = (pi_match_idx >= 0);
	ftree->gamma_pf_match = (gamma_match_idx >= 0);

	ftree->Kp_geo_match = (Kp_geo_match_idx >= 0);
	ftree->Km_geo_match = (Km_geo_match_idx >= 0);
	ftree->pi_geo_match = (pi_geo_match_idx >= 0);

	if (Kp_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Kp_match_idx];
		ftree->Kp_pf_idx  = Kp_match_idx;
		ftree->Kp_pf_key  = Kp_match_key;
		ftree->Kp_pf_pool = PassPool(pf, 2);
		ftree->Kp_pf_dR  = reco::deltaR(ftree->Gen_Kp_eta, ftree->Gen_Kp_phi, pf.eta(), pf.phi());
		ftree->Kp_pf_dpt = (pf.pt() - ftree->Gen_Kp_pt) / ftree->Gen_Kp_pt;
	}
	if (Km_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Km_match_idx];
		ftree->Km_pf_idx  = Km_match_idx;
		ftree->Km_pf_key  = Km_match_key;
		ftree->Km_pf_pool = PassPool(pf, 2);
		ftree->Km_pf_dR  = reco::deltaR(ftree->Gen_Km_eta, ftree->Gen_Km_phi, pf.eta(), pf.phi());
		ftree->Km_pf_dpt = (pf.pt() - ftree->Gen_Km_pt) / ftree->Gen_Km_pt;
	}
	if (pi_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[pi_match_idx];
		ftree->pi_pf_idx  = pi_match_idx;
		ftree->pi_pf_key  = pi_match_key;
		ftree->pi_pf_pool = PassPool(pf, 1);
		ftree->pi_pf_dR  = reco::deltaR(ftree->Gen_pi_eta, ftree->Gen_pi_phi, pf.eta(), pf.phi());
		ftree->pi_pf_dpt = (pf.pt() - ftree->Gen_pi_pt) / ftree->Gen_pi_pt;
	}

	// The geometric, plus the cross-tag: for the PF geometry picked, look up
	// its own association. isNull leaves key/pdgId at their null sentinel, which
	// offline reads as the association-null recovery bucket.
	if (Kp_geo_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Kp_geo_match_idx];
		ftree->Kp_geo_idx  = Kp_geo_match_idx;
		ftree->Kp_geo_pool = PassPool(pf, 2);
		ftree->Kp_geo_dR  = reco::deltaR(ftree->Gen_Kp_eta, ftree->Gen_Kp_phi, pf.eta(), pf.phi());
		ftree->Kp_geo_dpt = (pf.pt() - ftree->Gen_Kp_pt) / ftree->Gen_Kp_pt;
		ftree->Kp_geo_charge = pf.charge();
		ftree->Kp_geo_px = pf.px(); ftree->Kp_geo_py = pf.py(); ftree->Kp_geo_pz = pf.pz();
		ftree->Kp_geo_pt = pf.pt(); ftree->Kp_geo_eta = pf.eta();

		const pat::PackedCandidateRef pfRef(PFCandHandle, Kp_geo_match_idx);
		const reco::GenParticleRef genRef = (*xassocHandle)[pfRef];
		if (!genRef.isNull())
		{
			ftree->Kp_geo_assoc_key   = static_cast<Int_t>(genRef.key());
			ftree->Kp_geo_assoc_pdgId = genRef->pdgId();
		}
	}
	if (Km_geo_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[Km_geo_match_idx];
		ftree->Km_geo_idx  = Km_geo_match_idx;
		ftree->Km_geo_pool = PassPool(pf, 2);
		ftree->Km_geo_dR  = reco::deltaR(ftree->Gen_Km_eta, ftree->Gen_Km_phi, pf.eta(), pf.phi());
		ftree->Km_geo_dpt = (pf.pt() - ftree->Gen_Km_pt) / ftree->Gen_Km_pt;
		ftree->Km_geo_charge = pf.charge();
		ftree->Km_geo_px = pf.px(); ftree->Km_geo_py = pf.py(); ftree->Km_geo_pz = pf.pz();
		ftree->Km_geo_pt = pf.pt(); ftree->Km_geo_eta = pf.eta();

		const pat::PackedCandidateRef pfRef(PFCandHandle, Km_geo_match_idx);
		const reco::GenParticleRef genRef = (*xassocHandle)[pfRef];
		if (!genRef.isNull())
		{
			ftree->Km_geo_assoc_key   = static_cast<Int_t>(genRef.key());
			ftree->Km_geo_assoc_pdgId = genRef->pdgId();
		}
	}
	if (pi_geo_match_idx >= 0)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[pi_geo_match_idx];
		ftree->pi_geo_idx  = pi_geo_match_idx;
		ftree->pi_geo_pool = PassPool(pf, 1);
		ftree->pi_geo_dR  = reco::deltaR(ftree->Gen_pi_eta, ftree->Gen_pi_phi, pf.eta(), pf.phi());
		ftree->pi_geo_dpt = (pf.pt() - ftree->Gen_pi_pt) / ftree->Gen_pi_pt;
		ftree->pi_geo_charge = pf.charge();
		ftree->pi_geo_px = pf.px(); ftree->pi_geo_py = pf.py(); ftree->pi_geo_pz = pf.pz();
		ftree->pi_geo_pt = pf.pt(); ftree->pi_geo_eta = pf.eta();

		const pat::PackedCandidateRef pfRef(PFCandHandle, pi_geo_match_idx);
		const reco::GenParticleRef genRef = (*xassocHandle)[pfRef];
		if (!genRef.isNull())
		{
			ftree->pi_geo_assoc_key   = static_cast<Int_t>(genRef.key());
			ftree->pi_geo_assoc_pdgId = genRef->pdgId();
		}
	}

	ftree->tree->Fill();
}
void GeoMatch_phipimu::endJob()
{
	InfoTree->Fill();
}

Bool_t GeoMatch_phipimu::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
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

Bool_t GeoMatch_phipimu::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D GeoMatch_phipimu::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D GeoMatch_phipimu::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void GeoMatch_phipimu::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.add<bool>("doTruth", true);
	desc.add<std::string>("channel", "DsW");
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GeoMatch_phipimu);

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

#include "EDAnalyzers/Reconstruction/interface/RecoAll_DsstW_KKpimu_Tree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class RecoAll_DsstW_KKpimu : public edm::one::EDAnalyzer<>
{
public:
	explicit RecoAll_DsstW_KKpimu(const edm::ParameterSet &);
	~RecoAll_DsstW_KKpimu() override {};

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
	edm::EDGetTokenT<std::vector<reco::GenMET>> GenMETToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<pat::METCollection> MetToken_;
	edm::EDGetTokenT<pat::JetCollection> JetToken_;
	edm::EDGetTokenT<edm::ValueMap<Int_t>> tightJetIdToken_;

	KinematicParticleFactoryFromTransientTrack factory;
	KinematicParticleVertexFitter KPVFitter;

	const Float_t Mass_Dsst = 2.1066;
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
	RecoAll_DsstW_KKpimu_Tree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsLepton_;
	unsigned long nEventsGen_;
	unsigned long nEventsTriggered_;
};

RecoAll_DsstW_KKpimu::RecoAll_DsstW_KKpimu(const edm::ParameterSet &iConfig)
{
	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag GenMETTag_("genMetTrue");
	GenMETToken_ = consumes<std::vector<reco::GenMET>>(GenMETTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag TriggerBitsTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

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
}

void RecoAll_DsstW_KKpimu::beginJob()
{
	ftree = new RecoAll_DsstW_KKpimu_Tree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
	InfoTree = fs->make<TTree>("RunInfo", "RunInfo");

	InfoTree->Branch("nEventsProcessed_", &nEventsProcessed_);
	InfoTree->Branch("nEventsLepton_", &nEventsLepton_);
	InfoTree->Branch("nEventsGen_", &nEventsGen_);
	InfoTree->Branch("nEventsTriggered_", &nEventsTriggered_);

	nEventsProcessed_ = 0;
	nEventsLepton_ = 0;
	nEventsGen_ = 0;
	nEventsTriggered_ = 0;
}

void RecoAll_DsstW_KKpimu::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	ftree->Init();

	ftree->event_id = iEvent.id().event();
	ftree->run_id = iEvent.id().run();
	ftree->lumi_id = iEvent.luminosityBlock();

	Int_t Gen_H_idx = 0;
	Int_t Gen_Dsst_idx = 0;
	Int_t Gen_W_idx = 0;
	Int_t Gen_mu_idx = 0;
	Int_t Gen_Ds_idx = 0;
	Int_t Gen_gamma_idx = 0;
	Int_t Gen_Kp_idx = 0;
	Int_t Gen_Km_idx = 0;
	Int_t Gen_pi_idx = 0;

	Int_t num_Gen_H = 0;
	Int_t num_Gen_Dsst = 0;
	Int_t num_Gen_W = 0;
	Int_t num_Gen_mu = 0;
	Int_t num_Gen_e = 0;
	Int_t num_Gen_Ds = 0;
	Int_t num_Gen_gamma = 0;
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
		else if ((abs(gp.pdgId()) == 433))
		{
			Gen_Dsst_idx = i;
			num_Gen_Dsst++;
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
	if (!((num_Gen_H == 1) && (num_Gen_Dsst == 1) && (num_Gen_W == 1) && (num_Gen_mu == 1) && (num_Gen_e == 0)))
		return;
	nEventsLepton_++;

	const reco::GenParticle &H_gp = (*GenPartHandle)[Gen_H_idx];
	const reco::GenParticle &Dsst_gp = (*GenPartHandle)[Gen_Dsst_idx];
	const reco::GenParticle &W_gp = (*GenPartHandle)[Gen_W_idx];
	const reco::GenParticle &mu_gp = (*GenPartHandle)[Gen_mu_idx];

	const Int_t Dsst_charge = Dsst_gp.charge();

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];
		if ((gp.status() == 2) && gp.pdgId() == (Dsst_charge * 431))
		{
			if (TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()) && (gp.numberOfDaughters() > 2))
			{
				Bool_t daughter1_flag = false;
				Bool_t daughter2_flag = false;
				Bool_t daughter3_flag = false;
				for (size_t j = 0; j < gp.numberOfDaughters(); j++)
				{
					if (gp.daughter(j)->pdgId() == 321)
						daughter1_flag = true;
					if (gp.daughter(j)->pdgId() == -321)
						daughter2_flag = true;
					if (gp.daughter(j)->pdgId() == (Dsst_charge * 211))
						daughter3_flag = true;
				}
				if (daughter1_flag && daughter2_flag && daughter3_flag)
				{
					Gen_Ds_idx = i;
					num_Gen_Ds++;
				}
			}
		}
		else if (gp.status() == 1)
		{
			if (gp.pdgId() == (Dsst_charge * 321))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{
					Gen_Kp_idx = i;
					num_Gen_Kp++;
				}
			}
			else if (gp.pdgId() == (-Dsst_charge * 321))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{
					Gen_Km_idx = i;
					num_Gen_Km++;
				}
			}
			else if (gp.pdgId() == (Dsst_charge * 211))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{
					Gen_pi_idx = i;
					num_Gen_pi++;
				}
			}
			else if (gp.pdgId() == 22)
			{
				if ((gp.mother(0)->pdgId() == Dsst_gp.pdgId()) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
				{
					Gen_gamma_idx = i;
					num_Gen_gamma++;
				}
			}
		}
	}

	if (!((num_Gen_Ds == 1) && (num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1) && (num_Gen_gamma == 1)))
		return;

	const reco::GenParticle &Ds_gp = (*GenPartHandle)[Gen_Ds_idx];
	const reco::GenParticle &Kp_gp = (*GenPartHandle)[Gen_Kp_idx];
	const reco::GenParticle &Km_gp = (*GenPartHandle)[Gen_Km_idx];
	const reco::GenParticle &pi_gp = (*GenPartHandle)[Gen_pi_idx];
	const reco::GenParticle &gamma_gp = (*GenPartHandle)[Gen_gamma_idx];
	nEventsGen_++;

	TLorentzVector Gen_vis_tlv, Gen_Dsmu_tlv, Gen_Dsst_tlv, Gen_Ds_tlv, Gen_mu_tlv;
	Gen_Dsst_tlv.SetXYZM(Dsst_gp.px(), Dsst_gp.py(), Dsst_gp.pz(), Mass_Dsst);
	Gen_Ds_tlv.SetXYZM(Ds_gp.px(), Ds_gp.py(), Ds_gp.pz(), Mass_Ds);
	Gen_mu_tlv.SetXYZM(mu_gp.px(), mu_gp.py(), mu_gp.pz(), Mass_mu);
	Gen_Dsmu_tlv = Gen_Ds_tlv + Gen_mu_tlv;
	Gen_vis_tlv = Gen_Dsst_tlv + Gen_mu_tlv;

	edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
	iEvent.getByToken(GenMETToken_, GenMETHandle);
	if (GenMETHandle->empty())
		return;
	const reco::MET &GenMET = GenMETHandle->front();

	ftree->Gen_PV_vx = H_gp.vx();
	ftree->Gen_PV_vy = H_gp.vy();
	ftree->Gen_PV_vz = H_gp.vz();

	ftree->Gen_SV_vx = Ds_gp.vx();
	ftree->Gen_SV_vy = Ds_gp.vy();
	ftree->Gen_SV_vz = Ds_gp.vz();

	ftree->Gen_H_charge = H_gp.charge();
	ftree->Gen_H_eta = H_gp.eta();
	ftree->Gen_H_phi = H_gp.phi();
	ftree->Gen_H_px = H_gp.px();
	ftree->Gen_H_py = H_gp.py();
	ftree->Gen_H_pz = H_gp.pz();
	ftree->Gen_H_pt = H_gp.pt();
	ftree->Gen_H_p = H_gp.p();
	ftree->Gen_H_mass = H_gp.mass();

	ftree->Gen_Dsst_charge = Dsst_gp.charge();
	ftree->Gen_Dsst_eta = Dsst_gp.eta();
	ftree->Gen_Dsst_phi = Dsst_gp.phi();
	ftree->Gen_Dsst_px = Dsst_gp.px();
	ftree->Gen_Dsst_py = Dsst_gp.py();
	ftree->Gen_Dsst_pz = Dsst_gp.pz();
	ftree->Gen_Dsst_pt = Dsst_gp.pt();
	ftree->Gen_Dsst_p = Dsst_gp.p();
	ftree->Gen_Dsst_mass = Dsst_gp.mass();

	ftree->Gen_W_charge = W_gp.charge();
	ftree->Gen_W_eta = W_gp.eta();
	ftree->Gen_W_phi = W_gp.phi();
	ftree->Gen_W_px = W_gp.px();
	ftree->Gen_W_py = W_gp.py();
	ftree->Gen_W_pz = W_gp.pz();
	ftree->Gen_W_pt = W_gp.pt();
	ftree->Gen_W_p = W_gp.p();
	ftree->Gen_W_mass = W_gp.mass();

	ftree->Gen_Ds_charge = Ds_gp.charge();
	ftree->Gen_Ds_eta = Ds_gp.eta();
	ftree->Gen_Ds_phi = Ds_gp.phi();
	ftree->Gen_Ds_px = Ds_gp.px();
	ftree->Gen_Ds_py = Ds_gp.py();
	ftree->Gen_Ds_pz = Ds_gp.pz();
	ftree->Gen_Ds_pt = Ds_gp.pt();
	ftree->Gen_Ds_p = Ds_gp.p();
	ftree->Gen_Ds_mass = Ds_gp.mass();

	ftree->Gen_gamma_charge = gamma_gp.charge();
	ftree->Gen_gamma_eta = gamma_gp.eta();
	ftree->Gen_gamma_phi = gamma_gp.phi();
	ftree->Gen_gamma_px = gamma_gp.px();
	ftree->Gen_gamma_py = gamma_gp.py();
	ftree->Gen_gamma_pz = gamma_gp.pz();
	ftree->Gen_gamma_pt = gamma_gp.pt();
	ftree->Gen_gamma_p = gamma_gp.p();
	ftree->Gen_gamma_mass = gamma_gp.mass();

	ftree->Gen_Kp_charge = Kp_gp.charge();
	ftree->Gen_Kp_eta = Kp_gp.eta();
	ftree->Gen_Kp_phi = Kp_gp.phi();
	ftree->Gen_Kp_px = Kp_gp.px();
	ftree->Gen_Kp_py = Kp_gp.py();
	ftree->Gen_Kp_pz = Kp_gp.pz();
	ftree->Gen_Kp_pt = Kp_gp.pt();
	ftree->Gen_Kp_p = Kp_gp.p();
	ftree->Gen_Kp_mass = Kp_gp.mass();

	ftree->Gen_Km_charge = Km_gp.charge();
	ftree->Gen_Km_eta = Km_gp.eta();
	ftree->Gen_Km_phi = Km_gp.phi();
	ftree->Gen_Km_px = Km_gp.px();
	ftree->Gen_Km_py = Km_gp.py();
	ftree->Gen_Km_pz = Km_gp.pz();
	ftree->Gen_Km_pt = Km_gp.pt();
	ftree->Gen_Km_p = Km_gp.p();
	ftree->Gen_Km_mass = Km_gp.mass();

	ftree->Gen_pi_charge = pi_gp.charge();
	ftree->Gen_pi_eta = pi_gp.eta();
	ftree->Gen_pi_phi = pi_gp.phi();
	ftree->Gen_pi_px = pi_gp.px();
	ftree->Gen_pi_py = pi_gp.py();
	ftree->Gen_pi_pz = pi_gp.pz();
	ftree->Gen_pi_pt = pi_gp.pt();
	ftree->Gen_pi_p = pi_gp.p();
	ftree->Gen_pi_mass = pi_gp.mass();

	ftree->Gen_mu_charge = mu_gp.charge();
	ftree->Gen_mu_eta = mu_gp.eta();
	ftree->Gen_mu_phi = mu_gp.phi();
	ftree->Gen_mu_px = mu_gp.px();
	ftree->Gen_mu_py = mu_gp.py();
	ftree->Gen_mu_pz = mu_gp.pz();
	ftree->Gen_mu_pt = mu_gp.pt();
	ftree->Gen_mu_p = mu_gp.p();
	ftree->Gen_mu_mass = mu_gp.mass();

	ftree->Gen_Dsmu_charge = Ds_gp.charge() + mu_gp.charge();
	ftree->Gen_Dsmu_eta = Gen_Dsmu_tlv.Eta();
	ftree->Gen_Dsmu_phi = Gen_Dsmu_tlv.Phi();
	ftree->Gen_Dsmu_px = Gen_Dsmu_tlv.Px();
	ftree->Gen_Dsmu_py = Gen_Dsmu_tlv.Py();
	ftree->Gen_Dsmu_pz = Gen_Dsmu_tlv.Pz();
	ftree->Gen_Dsmu_pt = Gen_Dsmu_tlv.Pt();
	ftree->Gen_Dsmu_p = Gen_Dsmu_tlv.P();
	ftree->Gen_Dsmu_mass = Gen_Dsmu_tlv.M();

	ftree->Gen_vis_charge = Dsst_gp.charge() + mu_gp.charge();
	ftree->Gen_vis_eta = Gen_vis_tlv.Eta();
	ftree->Gen_vis_phi = Gen_vis_tlv.Phi();
	ftree->Gen_vis_px = Gen_vis_tlv.Px();
	ftree->Gen_vis_py = Gen_vis_tlv.Py();
	ftree->Gen_vis_pz = Gen_vis_tlv.Pz();
	ftree->Gen_vis_pt = Gen_vis_tlv.Pt();
	ftree->Gen_vis_p = Gen_vis_tlv.P();
	ftree->Gen_vis_mass = Gen_vis_tlv.M();

	ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_gp, Km_gp);
	ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_gp, pi_gp);
	ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp);
	ftree->Gen_dR_Km_pi = reco::deltaR(Km_gp, pi_gp);
	ftree->Gen_dR_Km_Ds = reco::deltaR(Km_gp, Ds_gp);
	ftree->Gen_dR_pi_Ds = reco::deltaR(pi_gp, Ds_gp);

	ftree->Gen_dR_Kp_gamma = reco::deltaR(Kp_gp, gamma_gp);
	ftree->Gen_dR_Kp_Dsst = reco::deltaR(Kp_gp, Dsst_gp);
	ftree->Gen_dR_Km_gamma = reco::deltaR(Km_gp, gamma_gp);
	ftree->Gen_dR_Km_Dsst = reco::deltaR(Km_gp, Dsst_gp);
	ftree->Gen_dR_pi_gamma = reco::deltaR(pi_gp, gamma_gp);
	ftree->Gen_dR_pi_Dsst = reco::deltaR(pi_gp, Dsst_gp);
	ftree->Gen_dR_Ds_gamma = reco::deltaR(Ds_gp, gamma_gp);
	ftree->Gen_dR_Ds_Dsst = reco::deltaR(Ds_gp, Dsst_gp);
	ftree->Gen_dR_gamma_Dsst = reco::deltaR(gamma_gp, Dsst_gp);

	ftree->Gen_dR_Kp_mu = reco::deltaR(Kp_gp, mu_gp);
	ftree->Gen_dR_Km_mu = reco::deltaR(Km_gp, mu_gp);
	ftree->Gen_dR_pi_mu = reco::deltaR(pi_gp, mu_gp);
	ftree->Gen_dR_Ds_mu = reco::deltaR(Ds_gp, mu_gp);
	ftree->Gen_dR_gamma_mu = reco::deltaR(gamma_gp, mu_gp);
	ftree->Gen_dR_Dsst_mu = reco::deltaR(Dsst_gp, mu_gp);

	ftree->Gen_met_phi = GenMET.phi();
	ftree->Gen_met_px = GenMET.px();
	ftree->Gen_met_py = GenMET.py();
	ftree->Gen_met_pt = GenMET.pt();
	ftree->Gen_met_sumEt = GenMET.sumEt();

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
			if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
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
			if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
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
			if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
			{
				if (deltaR < pi_deltaR)
				{
					pi_match_idx = i;
					pi_deltaR = deltaR;
				}
			}
		}
	}

	Int_t gamma_match_idx = -1;
	Float_t gamma_deltaR = 0.02;

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (pf.pdgId() != 22)
			continue;
		Float_t deltaR = reco::deltaR(pf.eta(), pf.phi(), gamma_gp.eta(), gamma_gp.phi());
		Float_t deltapt = abs((gamma_gp.pt() - pf.pt()) / gamma_gp.pt());
		if(deltapt > 0.5)
			continue;
		if (deltaR < gamma_deltaR)
		{
			gamma_match_idx = i;
			gamma_deltaR = deltaR;
		}
	}

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

		if (TwoBodyDecay(gp, W_gp.pdgId(), mu_gp.pdgId(), mu_gp.charge() * 14) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
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
	// Triggers
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	Bool_t trig_HLT_IsoMu24 = false;
	Bool_t trig_HLT_Mu50 = false;
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);
	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		Bool_t pass = triggerHandle->accept(i);
		if (trigName.Contains("HLT_IsoMu24_v"))
			trig_HLT_IsoMu24 = pass;
		else if (trigName.Contains("HLT_Mu50_v"))
			trig_HLT_Mu50 = pass;
	}
	if (!(trig_HLT_IsoMu24) && !(trig_HLT_Mu50))
		return;
	nEventsTriggered_++;
	ftree->trig_HLT_IsoMu24 = trig_HLT_IsoMu24;
	ftree->trig_HLT_Mu50 = trig_HLT_Mu50;

	//==================================================//
	//--------------------------------------------------//
	// primary vertex
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);

	if (primvtxHandle->empty())
		return;
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
	// MET
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<pat::METCollection> MetHandle;
	iEvent.getByToken(MetToken_, MetHandle);
	if (MetHandle->empty())
		return;
	const pat::MET &Met = MetHandle->front();

	// MET kinematics
	ftree->met_sig = Met.metSignificance();
	ftree->met_phi = Met.phi();
	ftree->met_px = Met.px();
	ftree->met_py = Met.py();
	ftree->met_pt = Met.pt();
	ftree->met_sumEt = Met.sumEt();

	//==================================================//
	//--------------------------------------------------//
	// Jet
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<pat::JetCollection> JetHandle;
	iEvent.getByToken(JetToken_, JetHandle);

	edm::Handle<edm::ValueMap<Int_t>> tightJetIdHandle;
	iEvent.getByToken(tightJetIdToken_, tightJetIdHandle);

	for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
	{
		Int_t tightJetId = (*tightJetIdHandle)[edm::Ref<pat::JetCollection>(JetHandle, idx_jet)];
		if (tightJetId == 0)
			continue;
		ftree->num_jets++;
		const pat::Jet &Jet = (*JetHandle)[idx_jet];
		if (Jet.pt() > 20)
			ftree->num_jets_pt20++;
		if (Jet.pt() > 30)
			ftree->num_jets_pt30++;
		if (Jet.pt() > 40)
			ftree->num_jets_pt40++;
		if (Jet.pt() > 50)
			ftree->num_jets_pt50++;
		if (Jet.pt() > 60)
			ftree->num_jets_pt60++;
		if (Jet.pt() > 70)
			ftree->num_jets_pt70++;
		if (Jet.pt() > 80)
			ftree->num_jets_pt80++;
		if (Jet.pt() > 90)
			ftree->num_jets_pt90++;
		if (Jet.pt() > 100)
			ftree->num_jets_pt100++;
	}

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
	const pat::Muon &mu = (*MuonHandle)[best_mu_idx];
	TLorentzVector mu_tlv;
	mu_tlv.SetXYZM(mu.px(), mu.py(), mu.pz(), Mass_mu);

	std::pair<Bool_t, Measurement1D> mu_ipresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(mu.muonBestTrack()), primvtx);
	std::pair<Bool_t, Measurement1D> mu_ipxyresult = IPTools::absoluteTransverseImpactParameter((*ttBuilder).build(mu.muonBestTrack()), primvtx);

	Float_t mu_IsoR03_PFIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5f *mu.pfIsolationR03().sumPUPt);
	Float_t mu_IsoR04_PFIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5f *mu.pfIsolationR04().sumPUPt);

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (muon.charge() == mu.charge())
		{
			if (muon.isTightMuon(primvtx))
				ftree->num_rs_mu++;
		}
		else if (muon.charge() == -mu.charge())
		{
			if (muon.isTightMuon(primvtx))
				ftree->num_ws_mu++;
		}
	}
	
	Bool_t mu_match = (best_mu_idx == mu_match_idx);

	//--------------------------------------------------//
	// mu storage
	//--------------------------------------------------//
	ftree->mu_charge = mu.charge();
	ftree->mu_eta = mu.eta();
	ftree->mu_phi = mu.phi();
	ftree->mu_px = mu.px();
	ftree->mu_py = mu.py();
	ftree->mu_pz = mu.pz();
	ftree->mu_pt = mu.pt();
	ftree->mu_p = mu.p();
	ftree->mu_isHighPt = mu.isHighPtMuon(primvtx);
	ftree->mu_isLoose = mu.isLooseMuon();
	ftree->mu_isMedium = mu.isMediumMuon();
	ftree->mu_isSoft = mu.isSoftMuon(primvtx);
	ftree->mu_isTight = mu.isTightMuon(primvtx);
	ftree->mu_isPF = mu.isPFMuon();
	ftree->mu_isTracker = mu.isTrackerMuon();
	ftree->mu_isGlobal = mu.isGlobalMuon();

	if (mu_ipresult.first)
	{
		ftree->mu_ip = mu_ipresult.second.value();
		ftree->mu_iperr = mu_ipresult.second.error();
		ftree->mu_ipsig = mu_ipresult.second.significance();
	}
	if (mu_ipxyresult.first)
	{
		ftree->mu_ipxy = mu_ipxyresult.second.value();
		ftree->mu_ipxyerr = mu_ipxyresult.second.error();
		ftree->mu_ipxysig = mu_ipxyresult.second.significance();
	}

	ftree->mu_IsoR03_sumChargedHadronPt = mu.pfIsolationR03().sumChargedHadronPt;
	ftree->mu_IsoR03_sumNeutralHadronEt = mu.pfIsolationR03().sumNeutralHadronEt;
	ftree->mu_IsoR03_sumPhotonEt = mu.pfIsolationR03().sumPhotonEt;
	ftree->mu_IsoR03_sumPUPt = mu.pfIsolationR03().sumPUPt;
	ftree->mu_IsoR03_PFIso = mu_IsoR03_PFIso;
	ftree->mu_IsoR03_RelIso = mu_IsoR03_PFIso / mu.pt();

	ftree->mu_IsoR04_sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
	ftree->mu_IsoR04_sumNeutralHadronEt = mu.pfIsolationR04().sumNeutralHadronEt;
	ftree->mu_IsoR04_sumPhotonEt = mu.pfIsolationR04().sumPhotonEt;
	ftree->mu_IsoR04_sumPUPt = mu.pfIsolationR04().sumPUPt;
	ftree->mu_IsoR04_PFIso = mu_IsoR04_PFIso;
	ftree->mu_IsoR04_RelIso = mu_IsoR04_PFIso / mu.pt();

	ftree->mu_match = mu_match;

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
		if (pf.pt() < 1)
			continue;
		if (abs(pf.pdgId()) != 211)
			continue;

		idx_pf.push_back(i);
	}

	if (idx_pf.size() < 3)
		return;
	
	// prepare gamma candidates
	std::vector<Int_t> idx_photon;
	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (pf.pt() < 0.5)
			continue;
		if (abs(pf.pdgId()) != 22)
			continue;

		idx_photon.push_back(i);
	}

	if (idx_photon.size() < 1)
		return;

	//--------------------------------------------------//
	// loop for K+
	//--------------------------------------------------//
	for (Int_t iKp = 0; iKp < static_cast<Int_t>(idx_pf.size()); iKp++)
	{
		const Int_t idx_Kp = idx_pf[iKp];

		const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[idx_Kp];
		// select on charge and pt
		if (Kp_pf.charge() == mu.charge())
			continue;
		if (Kp_pf.pt() < 2)
			continue;

		TLorentzVector Kp_p4;
		Kp_p4.SetXYZM(Kp_pf.px(), Kp_pf.py(), Kp_pf.pz(), Mass_K);

		reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.bestTrack());

		Int_t Kp_charge = Kp_pf.charge();
		Int_t Kp_PVassoc = Kp_pf.pvAssociationQuality();

		Bool_t Kp_match = (idx_Kp == Kp_match_idx);

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
			if (Km_pf.charge() != mu.charge())
				continue;
			if (Km_pf.pt() < 2)
				continue;

			TLorentzVector Km_p4;
			Km_p4.SetXYZM(Km_pf.px(), Km_pf.py(), Km_pf.pz(), Mass_K);

			reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.bestTrack());

			Int_t Km_charge = Km_pf.charge();
			Int_t Km_PVassoc = Km_pf.pvAssociationQuality();

			Bool_t Km_match = (idx_Km == Km_match_idx);

			//--------------------------------------------------//
			// loop for pi
			//--------------------------------------------------//
			for (Int_t ipi = 0; ipi < static_cast<Int_t>(idx_pf.size()); ipi++)
			{
				if ((iKp == ipi) || (iKm == ipi))
					continue;
				const Int_t idx_pi = idx_pf[ipi];
				const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];

				if (pi_pf.charge() == mu.charge())
					continue;

				// preliminary cut on Ds
				TLorentzVector pi_p4;
				pi_p4.SetXYZM(pi_pf.px(), pi_pf.py(), pi_pf.pz(), Mass_pi);
				if ((Kp_p4 + Km_p4 + pi_p4).Pt() < 5)
					continue;
				if (abs((Kp_p4 + Km_p4 + pi_p4).M() - Mass_Ds) > 0.3)
					continue;

				reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.bestTrack());

				Int_t pi_charge = pi_pf.charge();
				Int_t pi_PVassoc = pi_pf.pvAssociationQuality();

				Bool_t pi_match = (idx_pi == pi_match_idx);

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
				GlobalError Ds_vtx_cov = Ds_vertex->error();

				// select fit chi2 probability
				if (TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom()) < 0.001)
					continue;

				// DsTree particles: Ds
				RefCountedKinematicParticle Ds = DsTree->currentParticle();
				GlobalVector Ds_mtm = Ds->currentState().globalMomentum();
				const reco::TransientTrack Ds_rett = Ds->refittedTransientTrack();

				// select on Ds candidates
				if (abs(Ds->currentState().mass() - Mass_Ds) > 0.15)
					continue;
				if (Ds_mtm.perp() < 7)
					continue;

				// select dira
				GlobalVector displacement = Ds_vtx_pos - primvtx_pos;
				Float_t dira = displacement.dot(Ds_mtm) / (displacement.mag() * Ds_mtm.mag());
				if (dira < 0.99)
					continue;

				TLorentzVector Ds_tlv;
				Ds_tlv.SetXYZM(Ds_mtm.x(), Ds_mtm.y(), Ds_mtm.z(), Ds->currentState().mass());

				// DsTree particles: K+
				DsTree->movePointerToTheFirstChild();
				RefCountedKinematicParticle Kp = DsTree->currentParticle();
				GlobalVector Kp_mtm = Kp->currentState().globalMomentum();
				const reco::TransientTrack Kp_rett = Kp->refittedTransientTrack();
				// DsTree particles: K-
				DsTree->movePointerToTheNextChild();
				RefCountedKinematicParticle Km = DsTree->currentParticle();
				GlobalVector Km_mtm = Km->currentState().globalMomentum();
				const reco::TransientTrack Km_rett = Km->refittedTransientTrack();
				// DsTree particles: pi
				DsTree->movePointerToTheNextChild();
				RefCountedKinematicParticle pi = DsTree->currentParticle();
				GlobalVector pi_mtm = pi->currentState().globalMomentum();
				const reco::TransientTrack pi_rett = pi->refittedTransientTrack();

				// Construct resonances
				TLorentzVector Kp_tlv, Km_tlv, pi_tlv, KK_tlv, Kppi_tlv, Kmpi_tlv, Dsmu_tlv;
				Kp_tlv.SetXYZM(Kp_mtm.x(), Kp_mtm.y(), Kp_mtm.z(), Mass_K);
				Km_tlv.SetXYZM(Km_mtm.x(), Km_mtm.y(), Km_mtm.z(), Mass_K);
				pi_tlv.SetXYZM(pi_mtm.x(), pi_mtm.y(), pi_mtm.z(), Mass_pi);
				KK_tlv = Kp_tlv + Km_tlv;
				Kppi_tlv = Kp_tlv + pi_tlv;
				Kmpi_tlv = Km_tlv + pi_tlv;
				Dsmu_tlv = Ds_tlv + mu_tlv;

				Bool_t Ds_match = (idx_Kp == Kp_match_idx) && (idx_Km == Km_match_idx) && (idx_pi == pi_match_idx);
				
				//--------------------------------------------------//
				// Kinematic variables
				//--------------------------------------------------//
				Float_t Kp_eta = Kp_mtm.eta();
				Float_t Kp_phi = Kp_mtm.phi();
				Float_t Kp_px = Kp_mtm.x();
				Float_t Kp_py = Kp_mtm.y();
				Float_t Kp_pz = Kp_mtm.z();
				Float_t Kp_pt = Kp_mtm.perp();
				Float_t Kp_p = Kp_mtm.mag();
				Float_t Kp_mass = Kp->currentState().mass();

				Float_t Km_eta = Km_mtm.eta();
				Float_t Km_phi = Km_mtm.phi();
				Float_t Km_px = Km_mtm.x();
				Float_t Km_py = Km_mtm.y();
				Float_t Km_pz = Km_mtm.z();
				Float_t Km_pt = Km_mtm.perp();
				Float_t Km_p = Km_mtm.mag();
				Float_t Km_mass = Km->currentState().mass();

				Float_t pi_eta = pi_mtm.eta();
				Float_t pi_phi = pi_mtm.phi();
				Float_t pi_px = pi_mtm.x();
				Float_t pi_py = pi_mtm.y();
				Float_t pi_pz = pi_mtm.z();
				Float_t pi_pt = pi_mtm.perp();
				Float_t pi_p = pi_mtm.mag();
				Float_t pi_mass = pi->currentState().mass();

				Int_t KK_charge = Kp_pf.charge() + Km_pf.charge();
				Float_t KK_eta = KK_tlv.Eta();
				Float_t KK_phi = KK_tlv.Phi();
				Float_t KK_px = KK_tlv.Px();
				Float_t KK_py = KK_tlv.Py();
				Float_t KK_pz = KK_tlv.Pz();
				Float_t KK_pt = KK_tlv.Pt();
				Float_t KK_p = KK_tlv.P();
				Float_t KK_mass = KK_tlv.M();

				Int_t Kppi_charge = Kp_pf.charge() + pi_pf.charge();
				Float_t Kppi_eta = Kppi_tlv.Eta();
				Float_t Kppi_phi = Kppi_tlv.Phi();
				Float_t Kppi_px = Kppi_tlv.Px();
				Float_t Kppi_py = Kppi_tlv.Py();
				Float_t Kppi_pz = Kppi_tlv.Pz();
				Float_t Kppi_pt = Kppi_tlv.Pt();
				Float_t Kppi_p = Kppi_tlv.P();
				Float_t Kppi_mass = Kppi_tlv.M();

				Int_t Kmpi_charge = Km_pf.charge() + pi_pf.charge();
				Float_t Kmpi_eta = Kmpi_tlv.Eta();
				Float_t Kmpi_phi = Kmpi_tlv.Phi();
				Float_t Kmpi_px = Kmpi_tlv.Px();
				Float_t Kmpi_py = Kmpi_tlv.Py();
				Float_t Kmpi_pz = Kmpi_tlv.Pz();
				Float_t Kmpi_pt = Kmpi_tlv.Pt();
				Float_t Kmpi_p = Kmpi_tlv.P();
				Float_t Kmpi_mass = Kmpi_tlv.M();

				Int_t Ds_charge = Kp_pf.charge() + Km_pf.charge() + pi_pf.charge();
				Float_t Ds_eta = Ds_mtm.eta();
				Float_t Ds_phi = Ds_mtm.phi();
				Float_t Ds_px = Ds_mtm.x();
				Float_t Ds_py = Ds_mtm.y();
				Float_t Ds_pz = Ds_mtm.z();
				Float_t Ds_pt = Ds_mtm.perp();
				Float_t Ds_p = Ds_mtm.mag();
				Float_t Ds_mass = Ds->currentState().mass();

				Int_t Dsmu_charge = Kp_pf.charge() + Km_pf.charge() + pi_pf.charge() + mu.charge();
				Float_t Dsmu_eta = Dsmu_tlv.Eta();
				Float_t Dsmu_phi = Dsmu_tlv.Phi();
				Float_t Dsmu_px = Dsmu_tlv.Px();
				Float_t Dsmu_py = Dsmu_tlv.Py();
				Float_t Dsmu_pz = Dsmu_tlv.Pz();
				Float_t Dsmu_pt = Dsmu_tlv.Pt();
				Float_t Dsmu_p = Dsmu_tlv.P();
				Float_t Dsmu_mass = Dsmu_tlv.M();

				Float_t Ds_vx = Ds_vtx_pos.x();
				Float_t Ds_vy = Ds_vtx_pos.y();
				Float_t Ds_vz = Ds_vtx_pos.z();
				Float_t Ds_vxerr = Ds_vtx_cov.cxx();
				Float_t Ds_vyerr = Ds_vtx_cov.cyy();
				Float_t Ds_vzerr = Ds_vtx_cov.czz();
				Float_t Ds_chi2 = Ds_vertex->chiSquared();
				Int_t Ds_ndf = Ds_vertex->degreesOfFreedom();
				Float_t Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
				Float_t Ds_prob = TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom());

				// IP
				std::pair<Bool_t, Measurement1D> Kp_ipresult = IPTools::absoluteImpactParameter3D(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipresult = IPTools::absoluteImpactParameter3D(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipresult = IPTools::absoluteImpactParameter3D(Ds_rett, primvtx);

				std::pair<Bool_t, Measurement1D> Kp_ipxyresult = IPTools::absoluteTransverseImpactParameter(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipxyresult = IPTools::absoluteTransverseImpactParameter(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipxyresult = IPTools::absoluteTransverseImpactParameter(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipxyresult = IPTools::absoluteTransverseImpactParameter(Ds_rett, primvtx);

				Float_t Kp_ip = null;
				Float_t Kp_iperr = null;
				Float_t Kp_ipsig = null;
				if (Kp_ipresult.first)
				{
					Kp_ip = Kp_ipresult.second.value();
					Kp_iperr = Kp_ipresult.second.error();
					Kp_ipsig = Kp_ipresult.second.significance();
				}
				Float_t Kp_ipxy = null;
				Float_t Kp_ipxyerr = null;
				Float_t Kp_ipxysig = null;
				if (Kp_ipxyresult.first)
				{
					Kp_ipxy = Kp_ipxyresult.second.value();
					Kp_ipxyerr = Kp_ipxyresult.second.error();
					Kp_ipxysig = Kp_ipxyresult.second.significance();
				}
				Float_t Km_ip = null;
				Float_t Km_iperr = null;
				Float_t Km_ipsig = null;
				if (Km_ipresult.first)
				{
					Km_ip = Km_ipresult.second.value();
					Km_iperr = Km_ipresult.second.error();
					Km_ipsig = Km_ipresult.second.significance();
				}
				Float_t Km_ipxy = null;
				Float_t Km_ipxyerr = null;
				Float_t Km_ipxysig = null;
				if (Km_ipxyresult.first)
				{
					Km_ipxy = Km_ipxyresult.second.value();
					Km_ipxyerr = Km_ipxyresult.second.error();
					Km_ipxysig = Km_ipxyresult.second.significance();
				}
				Float_t pi_ip = null;
				Float_t pi_iperr = null;
				Float_t pi_ipsig = null;
				if (pi_ipresult.first)
				{
					pi_ip = pi_ipresult.second.value();
					pi_iperr = pi_ipresult.second.error();
					pi_ipsig = pi_ipresult.second.significance();
				}
				Float_t pi_ipxy = null;
				Float_t pi_ipxyerr = null;
				Float_t pi_ipxysig = null;
				if (pi_ipxyresult.first)
				{
					pi_ipxy = pi_ipxyresult.second.value();
					pi_ipxyerr = pi_ipxyresult.second.error();
					pi_ipxysig = pi_ipxyresult.second.significance();
				}
				Float_t Ds_ip = null;
				Float_t Ds_iperr = null;
				Float_t Ds_ipsig = null;
				if (Ds_ipresult.first)
				{
					Ds_ip = Ds_ipresult.second.value();
					Ds_iperr = Ds_ipresult.second.error();
					Ds_ipsig = Ds_ipresult.second.significance();
				}
				Float_t Ds_ipxy = null;
				Float_t Ds_ipxyerr = null;
				Float_t Ds_ipxysig = null;
				if (Ds_ipxyresult.first)
				{
					Ds_ipxy = Ds_ipxyresult.second.value();
					Ds_ipxyerr = Ds_ipxyresult.second.error();
					Ds_ipxysig = Ds_ipxyresult.second.significance();
				}

				// PV-SV separation
				Measurement1D Ds_FD_XY = distanceXY(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);
				Measurement1D Ds_FD_3D = distance3D(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);

				Float_t Ds_FDxy = Ds_FD_XY.value();
				Float_t Ds_FDxyerr = Ds_FD_XY.error();
				Float_t Ds_FDxysig = Ds_FD_XY.significance();

				Float_t Ds_FD = Ds_FD_3D.value();
				Float_t Ds_FDerr = Ds_FD_3D.error();
				Float_t Ds_FDsig = Ds_FD_3D.significance();

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
					if (i == idx_Kp)
						continue;
					if (i == idx_Km)
						continue;
					if (i == idx_pi)
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
				Float_t Ds_IsoR03_PFIso = Ds_IsoR03_sumChargedHadronPt + std::max(0.0f, Ds_IsoR03_sumNeutralHadronEt + Ds_IsoR03_sumPhotonEt - 0.5f *Ds_IsoR03_sumPUPt);
				Float_t Ds_IsoR03_RelIso = Ds_IsoR03_PFIso / Ds_mtm.perp();
				Float_t Ds_IsoR04_PFIso = Ds_IsoR04_sumChargedHadronPt + std::max(0.0f, Ds_IsoR04_sumNeutralHadronEt + Ds_IsoR04_sumPhotonEt - 0.5f *Ds_IsoR04_sumPUPt);
				Float_t Ds_IsoR04_RelIso = Ds_IsoR04_PFIso / Ds_mtm.perp();

				Float_t dR_Kp_Km = reco::deltaR(Kp_mtm, Km_mtm);
				Float_t dR_Kp_pi = reco::deltaR(Kp_mtm, pi_mtm);
				Float_t dR_Kp_Ds = reco::deltaR(Kp_mtm, Ds_mtm);
				Float_t dR_Km_pi = reco::deltaR(Km_mtm, pi_mtm);
				Float_t dR_Km_Ds = reco::deltaR(Km_mtm, Ds_mtm);
				Float_t dR_pi_Ds = reco::deltaR(pi_mtm, Ds_mtm);

				Float_t dR_Kp_mu = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), mu.eta(), mu.phi());
				Float_t dR_Km_mu = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), mu.eta(), mu.phi());
				Float_t dR_pi_mu = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), mu.eta(), mu.phi());
				Float_t dR_Ds_mu = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), mu.eta(), mu.phi());

				//--------------------------------------------------//
				// loop for gamma
				//--------------------------------------------------//
				for (Int_t igamma = 0; igamma < static_cast<Int_t>(idx_photon.size()); igamma++)
				{
					const Int_t idx_gamma = idx_photon[igamma];
					const pat::PackedCandidate &gamma_pf = (*PFCandHandle)[idx_gamma];
					Int_t gamma_charge = gamma_pf.charge();
					Int_t gamma_PVassoc = gamma_pf.pvAssociationQuality();

					Bool_t gamma_match = (idx_gamma == gamma_match_idx);
					Bool_t Dsst_match = Ds_match && gamma_match;
					Bool_t event_match = Dsst_match && mu_match;

					// preliminary cut on Ds
					TLorentzVector gamma_tlv, Dsst_tlv, vis_tlv;
					gamma_tlv.SetXYZM(gamma_pf.px(), gamma_pf.py(), gamma_pf.pz(), 0.0);
					Dsst_tlv = Ds_tlv + gamma_tlv;
					vis_tlv = Dsst_tlv + mu_tlv;
					if (abs(Dsst_tlv.M() - Mass_Dsst) > 0.3)
						continue;
					if (vis_tlv.M() > 125)
						continue;

					
					Float_t gamma_eta = gamma_tlv.Eta();
					Float_t gamma_phi = gamma_tlv.Phi();
					Float_t gamma_px = gamma_tlv.Px();
					Float_t gamma_py = gamma_tlv.Py();
					Float_t gamma_pz = gamma_tlv.Pz();
					Float_t gamma_pt = gamma_tlv.Pt();
					Float_t gamma_p = gamma_tlv.P();
					Float_t gamma_mass = gamma_tlv.M();

					Int_t Dsst_charge = Kp_pf.charge() + Km_pf.charge() + pi_pf.charge();
					Float_t Dsst_eta = Dsst_tlv.Eta();
					Float_t Dsst_phi = Dsst_tlv.Phi();
					Float_t Dsst_px = Dsst_tlv.Px();
					Float_t Dsst_py = Dsst_tlv.Py();
					Float_t Dsst_pz = Dsst_tlv.Pz();
					Float_t Dsst_pt = Dsst_tlv.Pt();
					Float_t Dsst_p = Dsst_tlv.P();
					Float_t Dsst_mass = Dsst_tlv.M();

					Int_t vis_charge = Kp_pf.charge() + Km_pf.charge() + pi_pf.charge() + mu.charge();
					Float_t vis_eta = vis_tlv.Eta();
					Float_t vis_phi = vis_tlv.Phi();
					Float_t vis_px = vis_tlv.Px();
					Float_t vis_py = vis_tlv.Py();
					Float_t vis_pz = vis_tlv.Pz();
					Float_t vis_pt = vis_tlv.Pt();
					Float_t vis_p = vis_tlv.P();
					Float_t vis_mass = vis_tlv.M();

					// Dsst Isolation
					Float_t Dsst_IsoR03_sumChargedHadronPt = 0;
					Float_t Dsst_IsoR03_sumNeutralHadronEt = 0;
					Float_t Dsst_IsoR03_sumPhotonEt = 0;
					Float_t Dsst_IsoR03_sumPUPt = 0;
					Float_t Dsst_IsoR04_sumChargedHadronPt = 0;
					Float_t Dsst_IsoR04_sumNeutralHadronEt = 0;
					Float_t Dsst_IsoR04_sumPhotonEt = 0;
					Float_t Dsst_IsoR04_sumPUPt = 0;
					for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
					{
						if (i == idx_Kp)
							continue;
						if (i == idx_Km)
							continue;
						if (i == idx_pi)
							continue;
						if (i == idx_gamma)
							continue;
						const pat::PackedCandidate &pf = (*PFCandHandle)[i];

						if (reco::deltaR(Dsst_tlv.Eta(), Dsst_tlv.Phi(), pf.eta(), pf.phi()) < 0.3)
						{
							if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
								Dsst_IsoR03_sumChargedHadronPt += pf.pt();
							else if (pf.pdgId() == 130)
								Dsst_IsoR03_sumNeutralHadronEt += pf.et();
							else if (pf.pdgId() == 22)
								Dsst_IsoR03_sumPhotonEt += pf.et();
							if (pf.fromPV() <= 1)
								Dsst_IsoR03_sumPUPt += pf.pt();
						}
						if (reco::deltaR(Dsst_tlv.Eta(), Dsst_tlv.Phi(), pf.eta(), pf.phi()) < 0.4)
						{
							if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
								Dsst_IsoR04_sumChargedHadronPt += pf.pt();
							else if (pf.pdgId() == 130)
								Dsst_IsoR04_sumNeutralHadronEt += pf.et();
							else if (pf.pdgId() == 22)
								Dsst_IsoR04_sumPhotonEt += pf.et();
							if (pf.fromPV() <= 1)
								Dsst_IsoR04_sumPUPt += pf.pt();
						}
					}
					Float_t Dsst_IsoR03_PFIso = Dsst_IsoR03_sumChargedHadronPt + std::max(0.0f, Dsst_IsoR03_sumNeutralHadronEt + Dsst_IsoR03_sumPhotonEt - 0.5f *Dsst_IsoR03_sumPUPt);
					Float_t Dsst_IsoR03_RelIso = Dsst_IsoR03_PFIso / Dsst_tlv.Pt();
					Float_t Dsst_IsoR04_PFIso = Dsst_IsoR04_sumChargedHadronPt + std::max(0.0f, Dsst_IsoR04_sumNeutralHadronEt + Dsst_IsoR04_sumPhotonEt - 0.5f *Dsst_IsoR04_sumPUPt);
					Float_t Dsst_IsoR04_RelIso = Dsst_IsoR04_PFIso / Dsst_tlv.Pt();

					Float_t dR_Kp_gamma = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					Float_t dR_Kp_Dsst = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					Float_t dR_Km_gamma = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					Float_t dR_Km_Dsst = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					Float_t dR_pi_gamma = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					Float_t dR_pi_Dsst = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					Float_t dR_Ds_gamma = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					Float_t dR_Ds_Dsst = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					Float_t dR_gamma_Dsst = reco::deltaR(gamma_pf.eta(), gamma_pf.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());

					Float_t dR_gamma_mu = reco::deltaR(gamma_pf.eta(), gamma_pf.phi(), mu.eta(), mu.phi());
					Float_t dR_Dsst_mu = reco::deltaR(Dsst_tlv.Eta(), Dsst_tlv.Phi(), mu.eta(), mu.phi());

					// refitted info from Kinematic Fit
					ftree->Kp_charge = Kp_charge;
					ftree->Kp_eta = Kp_eta;
					ftree->Kp_phi = Kp_phi;
					ftree->Kp_px = Kp_px;
					ftree->Kp_py = Kp_py;
					ftree->Kp_pz = Kp_pz;
					ftree->Kp_pt = Kp_pt;
					ftree->Kp_p = Kp_p;
					ftree->Kp_mass = Kp_mass;
					ftree->Kp_PVassoc = Kp_PVassoc;

					ftree->Km_charge = Km_charge;
					ftree->Km_eta = Km_eta;
					ftree->Km_phi = Km_phi;
					ftree->Km_px = Km_px;
					ftree->Km_py = Km_py;
					ftree->Km_pz = Km_pz;
					ftree->Km_pt = Km_pt;
					ftree->Km_p = Km_p;
					ftree->Km_mass = Km_mass;
					ftree->Km_PVassoc = Km_PVassoc;

					ftree->pi_charge = pi_charge;
					ftree->pi_eta = pi_eta;
					ftree->pi_phi = pi_phi;
					ftree->pi_px = pi_px;
					ftree->pi_py = pi_py;
					ftree->pi_pz = pi_pz;
					ftree->pi_pt = pi_pt;
					ftree->pi_p = pi_p;
					ftree->pi_mass = pi_mass;
					ftree->pi_PVassoc = pi_PVassoc;

					ftree->KK_charge = KK_charge;
					ftree->KK_eta = KK_eta;
					ftree->KK_phi = KK_phi;
					ftree->KK_px = KK_px;
					ftree->KK_py = KK_py;
					ftree->KK_pz = KK_pz;
					ftree->KK_pt = KK_pt;
					ftree->KK_p = KK_p;
					ftree->KK_mass = KK_mass;

					ftree->Kppi_charge = Kppi_charge;
					ftree->Kppi_eta = Kppi_eta;
					ftree->Kppi_phi = Kppi_phi;
					ftree->Kppi_px = Kppi_px;
					ftree->Kppi_py = Kppi_py;
					ftree->Kppi_pz = Kppi_pz;
					ftree->Kppi_pt = Kppi_pt;
					ftree->Kppi_p = Kppi_p;
					ftree->Kppi_mass = Kppi_mass;

					ftree->Kmpi_charge = Kmpi_charge;
					ftree->Kmpi_eta = Kmpi_eta;
					ftree->Kmpi_phi = Kmpi_phi;
					ftree->Kmpi_px = Kmpi_px;
					ftree->Kmpi_py = Kmpi_py;
					ftree->Kmpi_pz = Kmpi_pz;
					ftree->Kmpi_pt = Kmpi_pt;
					ftree->Kmpi_p = Kmpi_p;
					ftree->Kmpi_mass = Kmpi_mass;

					ftree->Ds_charge = Ds_charge;
					ftree->Ds_eta = Ds_eta;
					ftree->Ds_phi = Ds_phi;
					ftree->Ds_px = Ds_px;
					ftree->Ds_py = Ds_py;
					ftree->Ds_pz = Ds_pz;
					ftree->Ds_pt = Ds_pt;
					ftree->Ds_p = Ds_p;
					ftree->Ds_mass = Ds_mass;

					ftree->Dsmu_charge = Dsmu_charge;
					ftree->Dsmu_eta = Dsmu_eta;
					ftree->Dsmu_phi = Dsmu_phi;
					ftree->Dsmu_px = Dsmu_px;
					ftree->Dsmu_py = Dsmu_py;
					ftree->Dsmu_pz = Dsmu_pz;
					ftree->Dsmu_pt = Dsmu_pt;
					ftree->Dsmu_p = Dsmu_p;
					ftree->Dsmu_mass = Dsmu_mass;

					ftree->gamma_charge = gamma_charge;
					ftree->gamma_eta = gamma_eta;
					ftree->gamma_phi = gamma_phi;
					ftree->gamma_px = gamma_px;
					ftree->gamma_py = gamma_py;
					ftree->gamma_pz = gamma_pz;
					ftree->gamma_pt = gamma_pt;
					ftree->gamma_p = gamma_p;
					ftree->gamma_mass = gamma_mass;
					ftree->gamma_PVassoc = gamma_PVassoc;

					ftree->Dsst_charge = Dsst_charge;
					ftree->Dsst_eta = Dsst_eta;
					ftree->Dsst_phi = Dsst_phi;
					ftree->Dsst_px = Dsst_px;
					ftree->Dsst_py = Dsst_py;
					ftree->Dsst_pz = Dsst_pz;
					ftree->Dsst_pt = Dsst_pt;
					ftree->Dsst_p = Dsst_p;
					ftree->Dsst_mass = Dsst_mass;

					ftree->vis_charge = vis_charge;
					ftree->vis_eta = vis_eta;
					ftree->vis_phi = vis_phi;
					ftree->vis_px = vis_px;
					ftree->vis_py = vis_py;
					ftree->vis_pz = vis_pz;
					ftree->vis_pt = vis_pt;
					ftree->vis_p = vis_p;
					ftree->vis_mass = vis_mass;

					// refitted vertex
					ftree->Ds_vx = Ds_vx;
					ftree->Ds_vy = Ds_vy;
					ftree->Ds_vz = Ds_vz;
					ftree->Ds_vxerr = Ds_vxerr;
					ftree->Ds_vyerr = Ds_vyerr;
					ftree->Ds_vzerr = Ds_vzerr;
					ftree->Ds_chi2 = Ds_chi2;
					ftree->Ds_ndf = Ds_ndf;
					ftree->Ds_chi2ndf = Ds_chi2ndf;
					ftree->Ds_prob = Ds_prob;

					// angular distance
					ftree->dR_Kp_Km = dR_Kp_Km;
					ftree->dR_Kp_pi = dR_Kp_pi;
					ftree->dR_Kp_Ds = dR_Kp_Ds;
					ftree->dR_Km_pi = dR_Km_pi;
					ftree->dR_Km_Ds = dR_Km_Ds;
					ftree->dR_pi_Ds = dR_pi_Ds;

					ftree->dR_Kp_gamma = dR_Kp_gamma;
					ftree->dR_Kp_Dsst = dR_Kp_Dsst;
					ftree->dR_Km_gamma = dR_Km_gamma;
					ftree->dR_Km_Dsst = dR_Km_Dsst;
					ftree->dR_pi_gamma = dR_pi_gamma;
					ftree->dR_pi_Dsst = dR_pi_Dsst;
					ftree->dR_Ds_gamma = dR_Ds_gamma;
					ftree->dR_Ds_Dsst = dR_Ds_Dsst;
					ftree->dR_gamma_Dsst = dR_gamma_Dsst;

					ftree->dR_Kp_mu = dR_Kp_mu;
					ftree->dR_Km_mu = dR_Km_mu;
					ftree->dR_pi_mu = dR_pi_mu;
					ftree->dR_Ds_mu = dR_Ds_mu;
					ftree->dR_gamma_mu = dR_gamma_mu;
					ftree->dR_Dsst_mu = dR_Dsst_mu;

					// IP
					ftree->Kp_ip = Kp_ip;
					ftree->Kp_iperr = Kp_iperr;
					ftree->Kp_ipsig = Kp_ipsig;
					ftree->Kp_ipxy = Kp_ipxy;
					ftree->Kp_ipxyerr = Kp_ipxyerr;
					ftree->Kp_ipxysig = Kp_ipxysig;

					ftree->Km_ip = Km_ip;
					ftree->Km_iperr = Km_iperr;
					ftree->Km_ipsig = Km_ipsig;
					ftree->Km_ipxy = Km_ipxy;
					ftree->Km_ipxyerr = Km_ipxyerr;
					ftree->Km_ipxysig = Km_ipxysig;

					ftree->pi_ip = pi_ip;
					ftree->pi_iperr = pi_iperr;
					ftree->pi_ipsig = pi_ipsig;
					ftree->pi_ipxy = pi_ipxy;
					ftree->pi_ipxyerr = pi_ipxyerr;
					ftree->pi_ipxysig = pi_ipxysig;

					ftree->Ds_ip = Ds_ip;
					ftree->Ds_iperr = Ds_iperr;
					ftree->Ds_ipsig = Ds_ipsig;
					ftree->Ds_ipxy = Ds_ipxy;
					ftree->Ds_ipxyerr = Ds_ipxyerr;
					ftree->Ds_ipxysig = Ds_ipxysig;

					// PV-SV variable
					ftree->dira = dira;

					ftree->Ds_FDxy = Ds_FDxy;
					ftree->Ds_FDxyerr = Ds_FDxyerr;
					ftree->Ds_FDxysig = Ds_FDxysig;

					ftree->Ds_FD = Ds_FD;
					ftree->Ds_FDerr = Ds_FDerr;
					ftree->Ds_FDsig = Ds_FDsig;

					// Ds isolation
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

					ftree->Dsst_IsoR03_sumChargedHadronPt = Dsst_IsoR03_sumChargedHadronPt;
					ftree->Dsst_IsoR03_sumNeutralHadronEt = Dsst_IsoR03_sumNeutralHadronEt;
					ftree->Dsst_IsoR03_sumPhotonEt = Dsst_IsoR03_sumPhotonEt;
					ftree->Dsst_IsoR03_sumPUPt = Dsst_IsoR03_sumPUPt;
					ftree->Dsst_IsoR03_PFIso = Dsst_IsoR03_PFIso;
					ftree->Dsst_IsoR03_RelIso = Dsst_IsoR03_RelIso;

					ftree->Dsst_IsoR04_sumChargedHadronPt = Dsst_IsoR04_sumChargedHadronPt;
					ftree->Dsst_IsoR04_sumNeutralHadronEt = Dsst_IsoR04_sumNeutralHadronEt;
					ftree->Dsst_IsoR04_sumPhotonEt = Dsst_IsoR04_sumPhotonEt;
					ftree->Dsst_IsoR04_sumPUPt = Dsst_IsoR04_sumPUPt;
					ftree->Dsst_IsoR04_PFIso = Dsst_IsoR04_PFIso;
					ftree->Dsst_IsoR04_RelIso = Dsst_IsoR04_RelIso;

					ftree->Kp_match = Kp_match;
					ftree->Km_match = Km_match;
					ftree->pi_match = pi_match;
					ftree->Ds_match = Ds_match;
					ftree->gamma_match = gamma_match;
					ftree->Dsst_match = Dsst_match;
					ftree->event_match = event_match;
					ftree->tree->Fill();
				}
			}
		}
	}
}

void RecoAll_DsstW_KKpimu::endJob()
{
	InfoTree->Fill();
}

Bool_t RecoAll_DsstW_KKpimu::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
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

Bool_t RecoAll_DsstW_KKpimu::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D RecoAll_DsstW_KKpimu::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D RecoAll_DsstW_KKpimu::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void RecoAll_DsstW_KKpimu::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoAll_DsstW_KKpimu);

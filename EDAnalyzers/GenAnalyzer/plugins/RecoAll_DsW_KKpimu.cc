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

#include "EDAnalyzers/GenAnalyzer/interface/RecoAll_DsW_KKpimu_Tree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class RecoAll_DsW_KKpimu : public edm::one::EDAnalyzer<>
{
public:
	explicit RecoAll_DsW_KKpimu(const edm::ParameterSet &);
	~RecoAll_DsW_KKpimu() override {};

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
	edm::EDGetTokenT<reco::GenJetCollection> GenJetToken_;
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
	RecoAll_DsW_KKpimu_Tree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsLepton_;
	unsigned long nEventsGen_;
	unsigned long nEventsTriggered_;
};

RecoAll_DsW_KKpimu::RecoAll_DsW_KKpimu(const edm::ParameterSet &iConfig)
{
	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag GenJetTag_("slimmedGenJets");
	GenJetToken_ = consumes<reco::GenJetCollection>(GenJetTag_);

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

void RecoAll_DsW_KKpimu::beginJob()
{
	ftree = new RecoAll_DsW_KKpimu_Tree(fs->make<TTree>("Events", "Events"));
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

void RecoAll_DsW_KKpimu::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

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
	nEventsLepton_++;

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
	nEventsGen_++;

	TLorentzVector Gen_vis_tlv, Gen_Ds_tlv, Gen_mu_tlv;
	Gen_Ds_tlv.SetXYZM(Ds_gp.px(), Ds_gp.py(), Ds_gp.pz(), Mass_Ds);
	Gen_mu_tlv.SetXYZM(mu_gp.px(), mu_gp.py(), mu_gp.pz(), Mass_mu);
	Gen_vis_tlv = Gen_Ds_tlv + Gen_mu_tlv;

	edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
	iEvent.getByToken(GenMETToken_, GenMETHandle);
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

	ftree->Gen_Ds_charge = Ds_gp.charge();
	ftree->Gen_Ds_eta = Ds_gp.eta();
	ftree->Gen_Ds_phi = Ds_gp.phi();
	ftree->Gen_Ds_px = Ds_gp.px();
	ftree->Gen_Ds_py = Ds_gp.py();
	ftree->Gen_Ds_pz = Ds_gp.pz();
	ftree->Gen_Ds_pt = Ds_gp.pt();
	ftree->Gen_Ds_p = Ds_gp.p();
	ftree->Gen_Ds_mass = Ds_gp.mass();

	ftree->Gen_W_charge = W_gp.charge();
	ftree->Gen_W_eta = W_gp.eta();
	ftree->Gen_W_phi = W_gp.phi();
	ftree->Gen_W_px = W_gp.px();
	ftree->Gen_W_py = W_gp.py();
	ftree->Gen_W_pz = W_gp.pz();
	ftree->Gen_W_pt = W_gp.pt();
	ftree->Gen_W_p = W_gp.p();
	ftree->Gen_W_mass = W_gp.mass();

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

	ftree->Gen_vis_charge = Ds_gp.charge() + mu_gp.charge();
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
	ftree->Gen_dR_Kp_mu = reco::deltaR(Kp_gp, mu_gp);
	ftree->Gen_dR_Km_mu = reco::deltaR(Km_gp, mu_gp);
	ftree->Gen_dR_pi_mu = reco::deltaR(pi_gp, mu_gp);
	ftree->Gen_dR_Ds_mu = reco::deltaR(Ds_gp, mu_gp);

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

	Float_t mu_IsoR03_PFIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0.0, mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5 * mu.pfIsolationR03().sumPUPt);
	Float_t mu_IsoR04_PFIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0.0, mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumPUPt);

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

	ftree->mu_match = (best_mu_idx == mu_match_idx);

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
				if ((Ds_tlv + mu_tlv).M() > 125)
					continue;

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
				TLorentzVector Kp_tlv, Km_tlv, pi_tlv, KK_tlv, Kppi_tlv, Kmpi_tlv, vis_tlv;
				Kp_tlv.SetXYZM(Kp_mtm.x(), Kp_mtm.y(), Kp_mtm.z(), Mass_K);
				Km_tlv.SetXYZM(Km_mtm.x(), Km_mtm.y(), Km_mtm.z(), Mass_K);
				pi_tlv.SetXYZM(pi_mtm.x(), pi_mtm.y(), pi_mtm.z(), Mass_pi);
				KK_tlv = Kp_tlv + Km_tlv;
				Kppi_tlv = Kp_tlv + pi_tlv;
				Kmpi_tlv = Km_tlv + pi_tlv;
				vis_tlv = Ds_tlv + mu_tlv;

				//--------------------------------------------------//
				// Kinematic variables
				//--------------------------------------------------//
				// IP
				std::pair<Bool_t, Measurement1D> Kp_ipresult = IPTools::absoluteImpactParameter3D(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipresult = IPTools::absoluteImpactParameter3D(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipresult = IPTools::absoluteImpactParameter3D(Ds_rett, primvtx);

				std::pair<Bool_t, Measurement1D> Kp_ipxyresult = IPTools::absoluteTransverseImpactParameter(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipxyresult = IPTools::absoluteTransverseImpactParameter(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipxyresult = IPTools::absoluteTransverseImpactParameter(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipxyresult = IPTools::absoluteTransverseImpactParameter(Ds_rett, primvtx);

				// PV-SV separation
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
				Float_t Ds_IsoR03_PFIso = Ds_IsoR03_sumChargedHadronPt + std::max(0.0, Ds_IsoR03_sumNeutralHadronEt + Ds_IsoR03_sumPhotonEt - 0.5 * Ds_IsoR03_sumPUPt);
				Float_t Ds_IsoR03_RelIso = (Ds_IsoR03_sumChargedHadronPt + std::max(0.0, Ds_IsoR03_sumNeutralHadronEt + Ds_IsoR03_sumPhotonEt - 0.5 * Ds_IsoR03_sumPUPt)) / Ds_mtm.perp();
				Float_t Ds_IsoR04_PFIso = Ds_IsoR04_sumChargedHadronPt + std::max(0.0, Ds_IsoR04_sumNeutralHadronEt + Ds_IsoR04_sumPhotonEt - 0.5 * Ds_IsoR04_sumPUPt);
				Float_t Ds_IsoR04_RelIso = (Ds_IsoR04_sumChargedHadronPt + std::max(0.0, Ds_IsoR04_sumNeutralHadronEt + Ds_IsoR04_sumPhotonEt - 0.5 * Ds_IsoR04_sumPUPt)) / Ds_mtm.perp();

				// refitted info from Kinematic Fit
				ftree->Kp_charge.push_back(Kp_pf.charge());
				ftree->Kp_eta.push_back(Kp_mtm.eta());
				ftree->Kp_phi.push_back(Kp_mtm.phi());
				ftree->Kp_px.push_back(Kp_mtm.x());
				ftree->Kp_py.push_back(Kp_mtm.y());
				ftree->Kp_pz.push_back(Kp_mtm.z());
				ftree->Kp_pt.push_back(Kp_mtm.perp());
				ftree->Kp_p.push_back(Kp_mtm.mag());
				ftree->Kp_mass.push_back(Kp->currentState().mass());
				ftree->Kp_PVassoc.push_back(Kp_pf.pvAssociationQuality());

				ftree->Km_charge.push_back(Km_pf.charge());
				ftree->Km_eta.push_back(Km_mtm.eta());
				ftree->Km_phi.push_back(Km_mtm.phi());
				ftree->Km_px.push_back(Km_mtm.x());
				ftree->Km_py.push_back(Km_mtm.y());
				ftree->Km_pz.push_back(Km_mtm.z());
				ftree->Km_pt.push_back(Km_mtm.perp());
				ftree->Km_p.push_back(Km_mtm.mag());
				ftree->Km_mass.push_back(Km->currentState().mass());
				ftree->Km_PVassoc.push_back(Km_pf.pvAssociationQuality());

				ftree->pi_charge.push_back(pi_pf.charge());
				ftree->pi_eta.push_back(pi_mtm.eta());
				ftree->pi_phi.push_back(pi_mtm.phi());
				ftree->pi_px.push_back(pi_mtm.x());
				ftree->pi_py.push_back(pi_mtm.y());
				ftree->pi_pz.push_back(pi_mtm.z());
				ftree->pi_pt.push_back(pi_mtm.perp());
				ftree->pi_p.push_back(pi_mtm.mag());
				ftree->pi_mass.push_back(pi->currentState().mass());
				ftree->pi_PVassoc.push_back(pi_pf.pvAssociationQuality());

				ftree->KK_charge.push_back(Kp_pf.charge() + Km_pf.charge());
				ftree->KK_eta.push_back(KK_tlv.Eta());
				ftree->KK_phi.push_back(KK_tlv.Phi());
				ftree->KK_px.push_back(KK_tlv.Px());
				ftree->KK_py.push_back(KK_tlv.Py());
				ftree->KK_pz.push_back(KK_tlv.Pz());
				ftree->KK_pt.push_back(KK_tlv.Pt());
				ftree->KK_p.push_back(KK_tlv.P());
				ftree->KK_mass.push_back(KK_tlv.M());

				ftree->Kppi_charge.push_back(Kp_pf.charge() + pi_pf.charge());
				ftree->Kppi_eta.push_back(Kppi_tlv.Eta());
				ftree->Kppi_phi.push_back(Kppi_tlv.Phi());
				ftree->Kppi_px.push_back(Kppi_tlv.Px());
				ftree->Kppi_py.push_back(Kppi_tlv.Py());
				ftree->Kppi_pz.push_back(Kppi_tlv.Pz());
				ftree->Kppi_pt.push_back(Kppi_tlv.Pt());
				ftree->Kppi_p.push_back(Kppi_tlv.P());
				ftree->Kppi_mass.push_back(Kppi_tlv.M());

				ftree->Kmpi_charge.push_back(Km_pf.charge() + pi_pf.charge());
				ftree->Kmpi_eta.push_back(Kmpi_tlv.Eta());
				ftree->Kmpi_phi.push_back(Kmpi_tlv.Phi());
				ftree->Kmpi_px.push_back(Kmpi_tlv.Px());
				ftree->Kmpi_py.push_back(Kmpi_tlv.Py());
				ftree->Kmpi_pz.push_back(Kmpi_tlv.Pz());
				ftree->Kmpi_pt.push_back(Kmpi_tlv.Pt());
				ftree->Kmpi_p.push_back(Kmpi_tlv.P());
				ftree->Kmpi_mass.push_back(Kmpi_tlv.M());

				ftree->Ds_charge.push_back(Kp_pf.charge() + Km_pf.charge() + pi_pf.charge());
				ftree->Ds_eta.push_back(Ds_mtm.eta());
				ftree->Ds_phi.push_back(Ds_mtm.phi());
				ftree->Ds_px.push_back(Ds_mtm.x());
				ftree->Ds_py.push_back(Ds_mtm.y());
				ftree->Ds_pz.push_back(Ds_mtm.z());
				ftree->Ds_pt.push_back(Ds_mtm.perp());
				ftree->Ds_p.push_back(Ds_mtm.mag());
				ftree->Ds_mass.push_back(Ds->currentState().mass());

				ftree->vis_charge.push_back(Kp_pf.charge() + Km_pf.charge() + pi_pf.charge() + mu.charge());
				ftree->vis_eta.push_back(vis_tlv.Eta());
				ftree->vis_phi.push_back(vis_tlv.Phi());
				ftree->vis_px.push_back(vis_tlv.Px());
				ftree->vis_py.push_back(vis_tlv.Py());
				ftree->vis_pz.push_back(vis_tlv.Pz());
				ftree->vis_pt.push_back(vis_tlv.Pt());
				ftree->vis_p.push_back(vis_tlv.P());
				ftree->vis_mass.push_back(vis_tlv.M());

				// refitted vertex
				ftree->Ds_vx.push_back(Ds_vtx_pos.x());
				ftree->Ds_vy.push_back(Ds_vtx_pos.y());
				ftree->Ds_vz.push_back(Ds_vtx_pos.z());
				ftree->Ds_vxerr.push_back(Ds_vtx_cov.cxx());
				ftree->Ds_vyerr.push_back(Ds_vtx_cov.cyy());
				ftree->Ds_vzerr.push_back(Ds_vtx_cov.czz());
				ftree->Ds_chi2.push_back(Ds_vertex->chiSquared());
				ftree->Ds_ndf.push_back(Ds_vertex->degreesOfFreedom());
				ftree->Ds_chi2ndf.push_back(Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom());
				ftree->Ds_prob.push_back(TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom()));

				// angular distance
				ftree->dR_Kp_Km.push_back(reco::deltaR(Kp_mtm, Km_mtm));
				ftree->dR_Kp_pi.push_back(reco::deltaR(Kp_mtm, pi_mtm));
				ftree->dR_Kp_Ds.push_back(reco::deltaR(Kp_mtm, Ds_mtm));
				ftree->dR_Km_pi.push_back(reco::deltaR(Km_mtm, pi_mtm));
				ftree->dR_Km_Ds.push_back(reco::deltaR(Km_mtm, Ds_mtm));
				ftree->dR_pi_Ds.push_back(reco::deltaR(pi_mtm, Ds_mtm));
				ftree->dR_Kp_mu.push_back(reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), mu.eta(), mu.phi()));
				ftree->dR_Km_mu.push_back(reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), mu.eta(), mu.phi()));
				ftree->dR_pi_mu.push_back(reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), mu.eta(), mu.phi()));
				ftree->dR_Ds_mu.push_back(reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), mu.eta(), mu.phi()));

				// IP
				if (Kp_ipresult.first)
				{
					ftree->Kp_ip.push_back(Kp_ipresult.second.value());
					ftree->Kp_iperr.push_back(Kp_ipresult.second.error());
					ftree->Kp_ipsig.push_back(Kp_ipresult.second.significance());
				}
				else
				{
					ftree->Kp_ip.push_back(null);
					ftree->Kp_iperr.push_back(null);
					ftree->Kp_ipsig.push_back(null);
				}
				if (Kp_ipxyresult.first)
				{
					ftree->Kp_ipxy.push_back(Kp_ipxyresult.second.value());
					ftree->Kp_ipxyerr.push_back(Kp_ipxyresult.second.error());
					ftree->Kp_ipxysig.push_back(Kp_ipxyresult.second.significance());
				}
				else
				{
					ftree->Kp_ipxy.push_back(null);
					ftree->Kp_ipxyerr.push_back(null);
					ftree->Kp_ipxysig.push_back(null);
				}
				if (Km_ipresult.first)
				{
					ftree->Km_ip.push_back(Km_ipresult.second.value());
					ftree->Km_iperr.push_back(Km_ipresult.second.error());
					ftree->Km_ipsig.push_back(Km_ipresult.second.significance());
				}
				else
				{
					ftree->Km_ip.push_back(null);
					ftree->Km_iperr.push_back(null);
					ftree->Km_ipsig.push_back(null);
				}
				if (Km_ipxyresult.first)
				{
					ftree->Km_ipxy.push_back(Km_ipxyresult.second.value());
					ftree->Km_ipxyerr.push_back(Km_ipxyresult.second.error());
					ftree->Km_ipxysig.push_back(Km_ipxyresult.second.significance());
				}
				else
				{
					ftree->Km_ipxy.push_back(null);
					ftree->Km_ipxyerr.push_back(null);
					ftree->Km_ipxysig.push_back(null);
				}
				if (pi_ipresult.first)
				{
					ftree->pi_ip.push_back(pi_ipresult.second.value());
					ftree->pi_iperr.push_back(pi_ipresult.second.error());
					ftree->pi_ipsig.push_back(pi_ipresult.second.significance());
				}
				else
				{
					ftree->pi_ip.push_back(null);
					ftree->pi_iperr.push_back(null);
					ftree->pi_ipsig.push_back(null);
				}
				if (pi_ipxyresult.first)
				{
					ftree->pi_ipxy.push_back(pi_ipxyresult.second.value());
					ftree->pi_ipxyerr.push_back(pi_ipxyresult.second.error());
					ftree->pi_ipxysig.push_back(pi_ipxyresult.second.significance());
				}
				else
				{
					ftree->pi_ipxy.push_back(null);
					ftree->pi_ipxyerr.push_back(null);
					ftree->pi_ipxysig.push_back(null);
				}
				if (Ds_ipresult.first)
				{
					ftree->Ds_ip.push_back(Ds_ipresult.second.value());
					ftree->Ds_iperr.push_back(Ds_ipresult.second.error());
					ftree->Ds_ipsig.push_back(Ds_ipresult.second.significance());
				}
				else
				{
					ftree->Ds_ip.push_back(null);
					ftree->Ds_iperr.push_back(null);
					ftree->Ds_ipsig.push_back(null);
				}
				if (Ds_ipxyresult.first)
				{
					ftree->Ds_ipxy.push_back(Ds_ipxyresult.second.value());
					ftree->Ds_ipxyerr.push_back(Ds_ipxyresult.second.error());
					ftree->Ds_ipxysig.push_back(Ds_ipxyresult.second.significance());
				}
				else
				{
					ftree->Ds_ipxy.push_back(null);
					ftree->Ds_ipxyerr.push_back(null);
					ftree->Ds_ipxysig.push_back(null);
				}

				// PV-SV variable
				ftree->dira.push_back(dira);

				ftree->Ds_FDxy.push_back(Ds_FD_XY.value());
				ftree->Ds_FDxyerr.push_back(Ds_FD_XY.error());
				ftree->Ds_FDxysig.push_back(Ds_FD_XY.significance());

				ftree->Ds_FD.push_back(Ds_FD_3D.value());
				ftree->Ds_FDerr.push_back(Ds_FD_3D.error());
				ftree->Ds_FDsig.push_back(Ds_FD_3D.significance());

				// Ds ioslation
				ftree->Ds_IsoR03_sumChargedHadronPt.push_back(Ds_IsoR03_sumChargedHadronPt);
				ftree->Ds_IsoR03_sumNeutralHadronEt.push_back(Ds_IsoR03_sumNeutralHadronEt);
				ftree->Ds_IsoR03_sumPhotonEt.push_back(Ds_IsoR03_sumPhotonEt);
				ftree->Ds_IsoR03_sumPUPt.push_back(Ds_IsoR03_sumPUPt);
				ftree->Ds_IsoR03_PFIso.push_back(Ds_IsoR03_PFIso);
				ftree->Ds_IsoR03_RelIso.push_back(Ds_IsoR03_RelIso);

				ftree->Ds_IsoR04_sumChargedHadronPt.push_back(Ds_IsoR04_sumChargedHadronPt);
				ftree->Ds_IsoR04_sumNeutralHadronEt.push_back(Ds_IsoR04_sumNeutralHadronEt);
				ftree->Ds_IsoR04_sumPhotonEt.push_back(Ds_IsoR04_sumPhotonEt);
				ftree->Ds_IsoR04_sumPUPt.push_back(Ds_IsoR04_sumPUPt);
				ftree->Ds_IsoR04_PFIso.push_back(Ds_IsoR04_PFIso);
				ftree->Ds_IsoR04_RelIso.push_back(Ds_IsoR04_RelIso);

				ftree->Kp_match.push_back((idx_Kp == Kp_match_idx));
				ftree->Km_match.push_back((idx_Km == Km_match_idx));
				ftree->pi_match.push_back((idx_pi == pi_match_idx));
				ftree->Ds_match.push_back((idx_Kp == Kp_match_idx) && (idx_Km == Km_match_idx) && (idx_pi == pi_match_idx));
				ftree->event_match.push_back((idx_Kp == Kp_match_idx) && (idx_Km == Km_match_idx) && (idx_pi == pi_match_idx) && ftree->mu_match);
			}
		}
	}

	ftree->nCandidates = ftree->Ds_match.size();
	ftree->tree->Fill();
}

void RecoAll_DsW_KKpimu::endJob()
{
	InfoTree->Fill();
}

Bool_t RecoAll_DsW_KKpimu::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
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

Bool_t RecoAll_DsW_KKpimu::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D RecoAll_DsW_KKpimu::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D RecoAll_DsW_KKpimu::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void RecoAll_DsW_KKpimu::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoAll_DsW_KKpimu);

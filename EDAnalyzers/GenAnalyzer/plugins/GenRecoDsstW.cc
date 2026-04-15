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

#include "EDAnalyzers/GenAnalyzer/interface/GenRecoDsstWTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class GenRecoDsstW : public edm::one::EDAnalyzer<>
{
public:
	explicit GenRecoDsstW(const edm::ParameterSet &);
	~GenRecoDsstW() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const;
	Bool_t TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;
	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<reco::GenJetCollection> GenJetToken_;
	edm::EDGetTokenT<std::vector<reco::GenMET>> GenMETToken_;

	const edm::Service<TFileService> fs;
	GenRecoDsstWTree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsTriggered_;
	unsigned long nEventsReconstructed_;
	unsigned long nEventsMatch_;
};

GenRecoDsstW::GenRecoDsstW(const edm::ParameterSet &iConfig)
{
	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag GenJetTag_("slimmedGenJets");
	GenJetToken_ = consumes<reco::GenJetCollection>(GenJetTag_);

	edm::InputTag GenMETTag_("genMetTrue");
	GenMETToken_ = consumes<std::vector<reco::GenMET>>(GenMETTag_);
}

void GenRecoDsstW::beginJob()
{
	ftree = new GenRecoDsstWTree(fs->make<TTree>("Events", "Events"));
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

void GenRecoDsstW::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	edm::Handle<reco::GenJetCollection> GenJetHandle;
	iEvent.getByToken(GenJetToken_, GenJetHandle);

	edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
	iEvent.getByToken(GenMETToken_, GenMETHandle);

	ftree->Init();

	Int_t Gen_H_idx = 0;
	Int_t Gen_Dsst_idx = 0;
	Int_t Gen_W_idx = 0;
	Int_t Gen_mu_idx = 0;
	Int_t Gen_e_idx = 0;
	Int_t Gen_l_idx = 0;
	Int_t Gen_Ds_idx = 0;
	Int_t Gen_Kp_idx = 0;
	Int_t Gen_Km_idx = 0;
	Int_t Gen_pi_idx = 0;
	Int_t Gen_gamma_idx = 0;

	Int_t num_Gen_H = 0;
	Int_t num_Gen_Dsst = 0;
	Int_t num_Gen_W = 0;
	Int_t num_Gen_mu = 0;
	Int_t num_Gen_e = 0;
	Int_t num_Gen_l = 0;
	Int_t num_Gen_Ds = 0;
	Int_t num_Gen_Kp = 0;
	Int_t num_Gen_Km = 0;
	Int_t num_Gen_pi = 0;
	Int_t num_Gen_gamma = 0;

	Bool_t MU_EVENT = false;

	//==================================================//
	//--------------------------------------------------//
	// Gen
	//--------------------------------------------------//
	//==================================================//
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
			Gen_e_idx = i;
			num_Gen_e++;
		}
	}

	num_Gen_l = num_Gen_mu + num_Gen_e;
	if (!((num_Gen_H == 1) && (num_Gen_Dsst == 1) && (num_Gen_W == 1) && (num_Gen_l == 1)))
		return;

	MU_EVENT = (num_Gen_mu == 1);
	ftree->MU_EVENT = MU_EVENT;

	if (MU_EVENT)
		Gen_l_idx = Gen_mu_idx;
	else
		Gen_l_idx = Gen_e_idx;

	const reco::GenParticle &H_gp = (*GenPartHandle)[Gen_H_idx];
	const reco::GenParticle &Dsst_gp = (*GenPartHandle)[Gen_Dsst_idx];
	const reco::GenParticle &W_gp = (*GenPartHandle)[Gen_W_idx];
	const reco::GenParticle &l_gp = (*GenPartHandle)[Gen_l_idx];

	const Int_t Dsst_charge = Dsst_gp.charge();

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];

		if (gp.status() == 2)
		{
			if (gp.pdgId() == (Dsst_charge * 431))
			{
				if (TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
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
				if (TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
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

	TLorentzVector Gen_Dsst_p4;
	Gen_Dsst_p4.SetXYZM(Dsst_gp.px(), Dsst_gp.py(), Dsst_gp.pz(), Dsst_gp.mass());
	TLorentzVector Gen_l_p4;
	Gen_l_p4.SetXYZM(l_gp.px(), l_gp.py(), l_gp.pz(), l_gp.mass());
	TLorentzVector Gen_vis_p4 = Gen_Dsst_p4 + Gen_l_p4;

	ftree->Gen_PV_vx = H_gp.vx();
	ftree->Gen_PV_vy = H_gp.vy();
	ftree->Gen_PV_vz = H_gp.vz();

	ftree->Gen_SV_vx = Dsst_gp.vx();
	ftree->Gen_SV_vy = Dsst_gp.vy();
	ftree->Gen_SV_vz = Dsst_gp.vz();

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

	ftree->Gen_l_charge = l_gp.charge();
	ftree->Gen_l_eta = l_gp.eta();
	ftree->Gen_l_phi = l_gp.phi();
	ftree->Gen_l_px = l_gp.px();
	ftree->Gen_l_py = l_gp.py();
	ftree->Gen_l_pz = l_gp.pz();
	ftree->Gen_l_pt = l_gp.pt();
	ftree->Gen_l_p = l_gp.p();
	ftree->Gen_l_mass = l_gp.mass();

	ftree->Gen_gamma_charge = gamma_gp.charge();
	ftree->Gen_gamma_eta = gamma_gp.eta();
	ftree->Gen_gamma_phi = gamma_gp.phi();
	ftree->Gen_gamma_px = gamma_gp.px();
	ftree->Gen_gamma_py = gamma_gp.py();
	ftree->Gen_gamma_pz = gamma_gp.pz();
	ftree->Gen_gamma_pt = gamma_gp.pt();
	ftree->Gen_gamma_p = gamma_gp.p();
	ftree->Gen_gamma_mass = gamma_gp.mass();

	ftree->Gen_vis_charge = Dsst_gp.charge() + l_gp.charge();
	ftree->Gen_vis_eta = Gen_vis_p4.Eta();
	ftree->Gen_vis_phi = Gen_vis_p4.Phi();
	ftree->Gen_vis_px = Gen_vis_p4.Px();
	ftree->Gen_vis_py = Gen_vis_p4.Py();
	ftree->Gen_vis_pz = Gen_vis_p4.Pz();
	ftree->Gen_vis_pt = Gen_vis_p4.Pt();
	ftree->Gen_vis_p = Gen_vis_p4.P();
	ftree->Gen_vis_mass = Gen_vis_p4.M();

	ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_gp, Km_gp);
	ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_gp, pi_gp);
	ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp);
	ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp);
	ftree->Gen_dR_Km_pi = reco::deltaR(Km_gp, pi_gp);
	ftree->Gen_dR_Km_Ds = reco::deltaR(Km_gp, Ds_gp);
	ftree->Gen_dR_pi_Ds = reco::deltaR(pi_gp, Ds_gp);
	ftree->Gen_dR_Kp_l = reco::deltaR(Kp_gp, l_gp);
	ftree->Gen_dR_Km_l = reco::deltaR(Km_gp, l_gp);
	ftree->Gen_dR_pi_l = reco::deltaR(pi_gp, l_gp);
	ftree->Gen_dR_Ds_l = reco::deltaR(Ds_gp, l_gp);

	Int_t Gen_num_jets_pt20 = 0;
	Int_t Gen_num_jets_pt30 = 0;
	Int_t Gen_num_jets_pt40 = 0;
	Int_t Gen_num_jets_pt50 = 0;
	Int_t Gen_num_jets_pt60 = 0;
	Int_t Gen_num_jets_pt70 = 0;
	Int_t Gen_num_jets_pt80 = 0;
	Int_t Gen_num_jets_pt90 = 0;
	Int_t Gen_num_jets_pt100 = 0;
	for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(GenJetHandle->size()); idx_jet++)
	{
		const reco::Jet &Jet = (*GenJetHandle)[idx_jet];
		if (Jet.pt() > 20)
			Gen_num_jets_pt20++;
		if (Jet.pt() > 30)
			Gen_num_jets_pt30++;
		if (Jet.pt() > 40)
			Gen_num_jets_pt40++;
		if (Jet.pt() > 50)
			Gen_num_jets_pt50++;
		if (Jet.pt() > 60)
			Gen_num_jets_pt60++;
		if (Jet.pt() > 70)
			Gen_num_jets_pt70++;
		if (Jet.pt() > 80)
			Gen_num_jets_pt80++;
		if (Jet.pt() > 90)
			Gen_num_jets_pt90++;
		if (Jet.pt() > 100)
			Gen_num_jets_pt100++;
	}
	ftree->Gen_num_jets = static_cast<Int_t>(GenJetHandle->size());
	ftree->Gen_num_jets_pt20 = Gen_num_jets_pt20;
	ftree->Gen_num_jets_pt30 = Gen_num_jets_pt30;
	ftree->Gen_num_jets_pt40 = Gen_num_jets_pt40;
	ftree->Gen_num_jets_pt50 = Gen_num_jets_pt50;
	ftree->Gen_num_jets_pt60 = Gen_num_jets_pt60;
	ftree->Gen_num_jets_pt70 = Gen_num_jets_pt70;
	ftree->Gen_num_jets_pt80 = Gen_num_jets_pt80;
	ftree->Gen_num_jets_pt90 = Gen_num_jets_pt90;
	ftree->Gen_num_jets_pt100 = Gen_num_jets_pt100;

	const reco::MET &GenMET = GenMETHandle->front();
	ftree->met_phi = GenMET.phi();
	ftree->met_px = GenMET.px();
	ftree->met_py = GenMET.py();
	ftree->met_pt = GenMET.pt();
	ftree->met_sumEt = GenMET.sumEt();

	ftree->tree->Fill();
}

Bool_t GenRecoDsstW::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
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

Bool_t GenRecoDsstW::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D GenRecoDsstW::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D GenRecoDsstW::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void GenRecoDsstW::endJob()
{
	InfoTree->Fill();
}

void GenRecoDsstW::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenRecoDsstW);

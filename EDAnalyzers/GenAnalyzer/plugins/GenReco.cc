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

#include "EDAnalyzers/GenAnalyzer/interface/GenRecoTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class GenReco : public edm::one::EDAnalyzer<>
{
public:
	explicit GenReco(const edm::ParameterSet &);
	~GenReco() override {};

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
	GenRecoTree *ftree;
	TTree *InfoTree;

	unsigned long nEventsProcessed_;
	unsigned long nEventsTriggered_;
	unsigned long nEventsReconstructed_;
	unsigned long nEventsMatch_;
};

GenReco::GenReco(const edm::ParameterSet &iConfig)
{
	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag GenJetTag_("slimmedGenJets");
	GenJetToken_ = consumes<reco::GenJetCollection>(GenJetTag_);

	edm::InputTag GenMETTag_("genMetTrue");
	GenMETToken_ = consumes<std::vector<reco::GenMET>>(GenMETTag_);
}

void GenReco::beginJob()
{
	ftree = new GenRecoTree(fs->make<TTree>("Events", "Events"));
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

void GenReco::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
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
	Int_t Gen_Ds_idx = 0;
	Int_t Gen_W_idx = 0;
	Int_t Gen_mu_idx = 0;
	Int_t Gen_e_idx = 0;
	Int_t Gen_l_idx = 0;
	Int_t Gen_Kp_idx = 0;
	Int_t Gen_Km_idx = 0;
	Int_t Gen_pi_idx = 0;

	Int_t num_Gen_H = 0;
	Int_t num_Gen_Ds = 0;
	Int_t num_Gen_W = 0;
	Int_t num_Gen_mu = 0;
	Int_t num_Gen_e = 0;
	Int_t num_Gen_l = 0;
	Int_t num_Gen_Kp = 0;
	Int_t num_Gen_Km = 0;
	Int_t num_Gen_pi = 0;

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
			Gen_e_idx = i;
			num_Gen_e++;
		}
	}

	num_Gen_l = num_Gen_mu + num_Gen_e;
	if (!((num_Gen_H == 1) && (num_Gen_Ds == 1) && (num_Gen_W == 1) && (num_Gen_l == 1)))
		return;

	MU_EVENT = (num_Gen_mu == 1);
	if (MU_EVENT)
		Gen_l_idx = Gen_mu_idx;
	else
		Gen_l_idx = Gen_e_idx;

	const reco::GenParticle &H_gp = (*GenPartHandle)[Gen_H_idx];
	const reco::GenParticle &Ds_gp = (*GenPartHandle)[Gen_Ds_idx];
	const reco::GenParticle &W_gp = (*GenPartHandle)[Gen_W_idx];
	const reco::GenParticle &l_gp = (*GenPartHandle)[Gen_l_idx];

	// const reco::Candidate *dtr = &Ds_gp;
	// while (dtr)
	// {
	// 	std::cout << "pdgId: " << dtr->pdgId() << " nDaughters: " << dtr->numberOfDaughters() << std::endl;
	// 	if (dtr->numberOfDaughters() == 0)
	// 	{
	// 		std::cout << "end" << std::endl;
	// 		break;
	// 	}
	// 	dtr = dtr->daughter(0);
	// }

	// const reco::Candidate* current = &Ds_gp;
    // while (true)
    // {
    //     bool foundSame = false;
	// 	const reco::Candidate* d = current->daughter(0);
	// 	if (d && d->pdgId() == current->pdgId())
	// 	{
	// 		current = d;
	// 		foundSame = true;
	// 	}
    //     if (!foundSame) break;
    // }

	// for (size_t i = 0; i < current->numberOfDaughters(); ++i)
	// {
	// 	std::cout << "daughters: " << current->daughter(i)->pdgId() << ", ";
	// }
	// std::cout << "end" << std::endl;

	const Int_t Ds_charge = Ds_gp.charge();

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];

		// if (!gp.isLastCopy())
		// 	continue;
		// if(gp.numberOfDaughters() != 0)
		// 	continue;
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

		// if (gp.pdgId() == (Ds_charge * 321))
		// {
		// 	Bool_t dec3 = ThreeBodyDecay(gp, Ds_gp.pdgId(), 321, -321, Ds_charge * 211);
		// 	Bool_t dec2 = TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId());
		// 	if(dec3 && !(dec2))
		// 		std::cout << "2 body fails" << std::endl;
		// 	else if(!(dec3) && dec2)
		// 		std::cout << "3 body fails" << std::endl;
		// }
	}

	// if (num_Gen_Kp == 0 && num_Gen_Km == 0 && num_Gen_pi == 0)
	// {
	// 	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	// 	{
	// 		const reco::GenParticle &gp = (*GenPartHandle)[i];

	// 		if (gp.pdgId() == (Ds_charge * 321))
	// 		{
	// 			Bool_t dec3 = ThreeBodyDecay(gp, Ds_gp.pdgId(), 321, -321, Ds_charge * 211);
	// 			Bool_t dec2 = TwoBodyDecay(gp, 25, Ds_gp.pdgId(), W_gp.pdgId());
	// 			if (dec3 || dec2)
	// 			{
	// 				const reco::Candidate *m = &gp;
	// 				while (m)
	// 				{
	// 					std::cout << "pdgId: " << m->pdgId()
	// 							  << "  status: " << m->status()
	// 							  << "  nMothers: " << m->numberOfMothers()
	// 							  << std::endl;
	// 					if (abs(m->pdgId()) == 431)
	// 					{
	// 						std::cout << "D_s daughters: ";
	// 						for (Int_t xxx = 0; xxx < static_cast<Int_t>(m->numberOfDaughters()); xxx++)
	// 						{
	// 							std::cout << m->daughter(xxx)->pdgId() << " ";
	// 						}
	// 						std::cout << "end of D_s daughters" << std::endl;
	// 					}
	// 					else if (abs(m->pdgId()) == 25)
	// 					{
	// 						std::cout << "H daughters: ";
	// 						for (Int_t xxx = 0; xxx < static_cast<Int_t>(m->numberOfDaughters()); xxx++)
	// 						{
	// 							std::cout << m->daughter(xxx)->pdgId() << " ";
	// 						}
	// 						std::cout << "end of H daughters" << std::endl;
	// 					}
	// 					if (m->numberOfMothers() == 0)
	// 					{
	// 						std::cout << "end" << std::endl;
	// 						break;
	// 					}
	// 					m = m->mother(0);
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	ftree->num_Gen_H = num_Gen_H;
	ftree->num_Gen_Ds = num_Gen_Ds;
	ftree->num_Gen_W = num_Gen_W;
	ftree->num_Gen_l = num_Gen_l;
	ftree->num_Gen_Kp = num_Gen_Kp;
	ftree->num_Gen_Km = num_Gen_Km;
	ftree->num_Gen_pi = num_Gen_pi;
	ftree->MU_EVENT = MU_EVENT;

	// // if (!((num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1)))
	// // 	return;

	const reco::GenParticle &Kp_gp = (*GenPartHandle)[Gen_Kp_idx];
	const reco::GenParticle &Km_gp = (*GenPartHandle)[Gen_Km_idx];
	const reco::GenParticle &pi_gp = (*GenPartHandle)[Gen_pi_idx];

	// TLorentzVector Gen_Ds_p4;
	// Gen_Ds_p4.SetXYZM(Ds_gp.px(), Ds_gp.py(), Ds_gp.pz(), Ds_gp.mass());
	// TLorentzVector Gen_mu_p4;
	// Gen_mu_p4.SetXYZM(mu_gp.px(), mu_gp.py(), mu_gp.pz(), mu_gp.mass());
	// TLorentzVector Gen_vis_p4 = Gen_Ds_p4 + Gen_mu_p4;

	// ftree->Gen_PV_vx = H_gp.vx();
	// ftree->Gen_PV_vy = H_gp.vy();
	// ftree->Gen_PV_vz = H_gp.vz();

	// ftree->Gen_SV_vx = Ds_gp.vx();
	// ftree->Gen_SV_vy = Ds_gp.vy();
	// ftree->Gen_SV_vz = Ds_gp.vz();

	// ftree->Gen_H_charge = H_gp.charge();
	// ftree->Gen_H_eta = H_gp.eta();
	// ftree->Gen_H_phi = H_gp.phi();
	// ftree->Gen_H_px = H_gp.px();
	// ftree->Gen_H_py = H_gp.py();
	// ftree->Gen_H_pz = H_gp.pz();
	// ftree->Gen_H_pt = H_gp.pt();
	// ftree->Gen_H_p = H_gp.p();
	// ftree->Gen_H_mass = H_gp.mass();

	// ftree->Gen_Ds_charge = Ds_gp.charge();
	// ftree->Gen_Ds_eta = Ds_gp.eta();
	// ftree->Gen_Ds_phi = Ds_gp.phi();
	// ftree->Gen_Ds_px = Ds_gp.px();
	// ftree->Gen_Ds_py = Ds_gp.py();
	// ftree->Gen_Ds_pz = Ds_gp.pz();
	// ftree->Gen_Ds_pt = Ds_gp.pt();
	// ftree->Gen_Ds_p = Ds_gp.p();
	// ftree->Gen_Ds_mass = Ds_gp.mass();

	// ftree->Gen_W_charge = W_gp.charge();
	// ftree->Gen_W_eta = W_gp.eta();
	// ftree->Gen_W_phi = W_gp.phi();
	// ftree->Gen_W_px = W_gp.px();
	// ftree->Gen_W_py = W_gp.py();
	// ftree->Gen_W_pz = W_gp.pz();
	// ftree->Gen_W_pt = W_gp.pt();
	// ftree->Gen_W_p = W_gp.p();
	// ftree->Gen_W_mass = W_gp.mass();

	// ftree->Gen_Kp_charge = Kp_gp.charge();
	// ftree->Gen_Kp_eta = Kp_gp.eta();
	// ftree->Gen_Kp_phi = Kp_gp.phi();
	// ftree->Gen_Kp_px = Kp_gp.px();
	// ftree->Gen_Kp_py = Kp_gp.py();
	// ftree->Gen_Kp_pz = Kp_gp.pz();
	// ftree->Gen_Kp_pt = Kp_gp.pt();
	// ftree->Gen_Kp_p = Kp_gp.p();
	// ftree->Gen_Kp_mass = Kp_gp.mass();

	// ftree->Gen_Km_charge = Km_gp.charge();
	// ftree->Gen_Km_eta = Km_gp.eta();
	// ftree->Gen_Km_phi = Km_gp.phi();
	// ftree->Gen_Km_px = Km_gp.px();
	// ftree->Gen_Km_py = Km_gp.py();
	// ftree->Gen_Km_pz = Km_gp.pz();
	// ftree->Gen_Km_pt = Km_gp.pt();
	// ftree->Gen_Km_p = Km_gp.p();
	// ftree->Gen_Km_mass = Km_gp.mass();

	// ftree->Gen_pi_charge = pi_gp.charge();
	// ftree->Gen_pi_eta = pi_gp.eta();
	// ftree->Gen_pi_phi = pi_gp.phi();
	// ftree->Gen_pi_px = pi_gp.px();
	// ftree->Gen_pi_py = pi_gp.py();
	// ftree->Gen_pi_pz = pi_gp.pz();
	// ftree->Gen_pi_pt = pi_gp.pt();
	// ftree->Gen_pi_p = pi_gp.p();
	// ftree->Gen_pi_mass = pi_gp.mass();

	// ftree->Gen_mu_charge = mu_gp.charge();
	// ftree->Gen_mu_eta = mu_gp.eta();
	// ftree->Gen_mu_phi = mu_gp.phi();
	// ftree->Gen_mu_px = mu_gp.px();
	// ftree->Gen_mu_py = mu_gp.py();
	// ftree->Gen_mu_pz = mu_gp.pz();
	// ftree->Gen_mu_pt = mu_gp.pt();
	// ftree->Gen_mu_p = mu_gp.p();
	// ftree->Gen_mu_mass = mu_gp.mass();

	// ftree->Gen_vis_charge = Ds_gp.charge() + mu_gp.charge();
	// ftree->Gen_vis_eta = Gen_vis_p4.Eta();
	// ftree->Gen_vis_phi = Gen_vis_p4.Phi();
	// ftree->Gen_vis_px = Gen_vis_p4.Px();
	// ftree->Gen_vis_py = Gen_vis_p4.Py();
	// ftree->Gen_vis_pz = Gen_vis_p4.Pz();
	// ftree->Gen_vis_pt = Gen_vis_p4.Pt();
	// ftree->Gen_vis_p = Gen_vis_p4.P();
	// ftree->Gen_vis_mass = Gen_vis_p4.M();

	// ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_gp, Km_gp);
	// ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_gp, pi_gp);
	// ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_gp, Ds_gp);
	// ftree->Gen_dR_Km_pi = reco::deltaR(Km_gp, pi_gp);
	// ftree->Gen_dR_Km_Ds = reco::deltaR(Km_gp, Ds_gp);
	// ftree->Gen_dR_pi_Ds = reco::deltaR(pi_gp, Ds_gp);
	// ftree->Gen_dR_Kp_mu = reco::deltaR(Kp_gp, mu_gp);
	// ftree->Gen_dR_Km_mu = reco::deltaR(Km_gp, mu_gp);
	// ftree->Gen_dR_pi_mu = reco::deltaR(pi_gp, mu_gp);
	// ftree->Gen_dR_Ds_mu = reco::deltaR(Ds_gp, mu_gp);

	// Int_t Gen_num_jets_pt20 = 0;
	// Int_t Gen_num_jets_pt30 = 0;
	// Int_t Gen_num_jets_pt40 = 0;
	// Int_t Gen_num_jets_pt50 = 0;
	// Int_t Gen_num_jets_pt60 = 0;
	// Int_t Gen_num_jets_pt70 = 0;
	// Int_t Gen_num_jets_pt80 = 0;
	// Int_t Gen_num_jets_pt90 = 0;
	// Int_t Gen_num_jets_pt100 = 0;
	// for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(GenJetHandle->size()); idx_jet++)
	// {
	// 	const reco::Jet &Jet = (*GenJetHandle)[idx_jet];
	// 	if (Jet.pt() > 20)
	// 		Gen_num_jets_pt20++;
	// 	if (Jet.pt() > 30)
	// 		Gen_num_jets_pt30++;
	// 	if (Jet.pt() > 40)
	// 		Gen_num_jets_pt40++;
	// 	if (Jet.pt() > 50)
	// 		Gen_num_jets_pt50++;
	// 	if (Jet.pt() > 60)
	// 		Gen_num_jets_pt60++;
	// 	if (Jet.pt() > 70)
	// 		Gen_num_jets_pt70++;
	// 	if (Jet.pt() > 80)
	// 		Gen_num_jets_pt80++;
	// 	if (Jet.pt() > 90)
	// 		Gen_num_jets_pt90++;
	// 	if (Jet.pt() > 100)
	// 		Gen_num_jets_pt100++;
	// }
	// ftree->Gen_num_jets = static_cast<Int_t>(GenJetHandle->size());
	// ftree->Gen_num_jets_pt20 = Gen_num_jets_pt20;
	// ftree->Gen_num_jets_pt30 = Gen_num_jets_pt30;
	// ftree->Gen_num_jets_pt40 = Gen_num_jets_pt40;
	// ftree->Gen_num_jets_pt50 = Gen_num_jets_pt50;
	// ftree->Gen_num_jets_pt60 = Gen_num_jets_pt60;
	// ftree->Gen_num_jets_pt70 = Gen_num_jets_pt70;
	// ftree->Gen_num_jets_pt80 = Gen_num_jets_pt80;
	// ftree->Gen_num_jets_pt90 = Gen_num_jets_pt90;
	// ftree->Gen_num_jets_pt100 = Gen_num_jets_pt100;

	// const reco::MET &GenMET = GenMETHandle->front();
	// ftree->met_phi = GenMET.phi();
	// ftree->met_px = GenMET.px();
	// ftree->met_py = GenMET.py();
	// ftree->met_pt = GenMET.pt();
	// ftree->met_sumEt = GenMET.sumEt();

	ftree->tree->Fill();
}

Bool_t GenReco::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
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

Bool_t GenReco::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

Measurement1D GenReco::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

Measurement1D GenReco::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
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

void GenReco::endJob()
{
	InfoTree->Fill();
}

void GenReco::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenReco);

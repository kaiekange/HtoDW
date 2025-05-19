// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      PreliminaryCut
//
/**\class PreliminaryCut PreliminaryCut.cc EDAnalyzer/GenParticleAnalyzer/plugins/PreliminaryCut.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  cmsusers user
//         Created:  Fri, 15 Nov 2024 15:09:59 GMT
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

#include "EDAnalyzers/GenParticleAnalyzer/interface/PFRecoTree.h"

class PreliminaryCut : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit PreliminaryCut(const edm::ParameterSet&);
        ~PreliminaryCut();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;

        edm::Service<TFileService> fs;
        PFRecoTree *ftree;

        std::vector<int> idx_Kp_vec;
        std::vector<int> idx_Km_vec;
        std::vector<int> idx_pi_vec;

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

PreliminaryCut::PreliminaryCut(const edm::ParameterSet& iConfig) :
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    ftree = new PFRecoTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

PreliminaryCut::~PreliminaryCut() {}

void PreliminaryCut::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    idx_Kp_vec.clear();
    idx_Km_vec.clear();
    idx_pi_vec.clear();

    ftree->Init();

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
    KalmanVertexFitter fitter(true);
    
    for(size_t i=0; i<packedPF->size(); i++){
        const auto& pf = (*packedPF)[i];

        if(!(pf.trackHighPurity())) continue;
        if(!(pf.hasTrackDetails())) continue;
        if(pf.pt() < 0.5) continue;
        if(pf.p() < 1) continue;

        if(pf.pdgId() == 211){
            idx_Kp_vec.push_back(i);
            idx_pi_vec.push_back(i);
        } else if(pf.pdgId() == -211){
            idx_Km_vec.push_back(i);
        }
    }

    for(size_t i=0; i<idx_Kp_vec.size(); i++){

        const auto& Kp_PF = (*packedPF)[idx_Kp_vec[i]];

        ftree->Kp_Reset();
        ftree->Kp_ETA = Kp_PF.eta();
        ftree->Kp_PHI = Kp_PF.phi();
        ftree->Kp_ORIVX_X = Kp_PF.vx();
        ftree->Kp_ORIVX_Y = Kp_PF.vy();
        ftree->Kp_ORIVX_Z = Kp_PF.vz();
        ftree->Kp_P = Kp_PF.p();
        ftree->Kp_PT = Kp_PF.pt();
        ftree->Kp_PX = Kp_PF.px();
        ftree->Kp_PY = Kp_PF.py();
        ftree->Kp_PZ = Kp_PF.pz();
        math::XYZPoint Kp_originvertex(Kp_PF.vx(), Kp_PF.vy(), Kp_PF.vz());
        TLorentzVector Kp_P4; Kp_P4.SetXYZM(Kp_PF.px(), Kp_PF.py(), Kp_PF.pz(), Mass_K);

        for(size_t j=0; j<idx_Km_vec.size(); j++){
            const auto& Km_PF = (*packedPF)[idx_Km_vec[j]];

            ftree->Km_Reset();
            ftree->Km_ETA = Km_PF.eta();
            ftree->Km_PHI = Km_PF.phi();
            ftree->Km_ORIVX_X = Km_PF.vx();
            ftree->Km_ORIVX_Y = Km_PF.vy();
            ftree->Km_ORIVX_Z = Km_PF.vz();
            ftree->Km_P = Km_PF.p();
            ftree->Km_PT = Km_PF.pt();
            ftree->Km_PX = Km_PF.px();
            ftree->Km_PY = Km_PF.py();
            ftree->Km_PZ = Km_PF.pz();
            math::XYZPoint Km_originvertex(Km_PF.vx(), Km_PF.vy(), Km_PF.vz());
            TLorentzVector Km_P4; Km_P4.SetXYZM(Km_PF.px(), Km_PF.py(), Km_PF.pz(), Mass_K);
            
            ftree->dR_Kp_Km = reco::deltaR(Kp_PF, Km_PF);
            ftree->d_Kp_Km = (Kp_originvertex - Km_originvertex).R();
            
            if( abs(ftree->Kp_PT - ftree->Km_PT) > 20) continue;
            if( ftree->dR_Kp_Km > 0.15 ) continue;
            if( ftree->d_Kp_Km > 0.2) continue;

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*theB).build(Kp_PF.pseudoTrack()),
                (*theB).build(Km_PF.pseudoTrack())
            };
            TransientVertex phi_Vertex = fitter.vertex(phi_Tracks);

            if( !(phi_Vertex.isValid()) ) continue;
            if( !(phi_Vertex.hasRefittedTracks()) ) continue;
           
            ftree->phi_Reset();

            std::vector<reco::TransientTrack> phiFit_Tracks = phi_Vertex.refittedTracks();

            TLorentzVector phiFit_Kp_P4;
            phiFit_Kp_P4.SetXYZM(phiFit_Tracks[0].track().px(), phiFit_Tracks[0].track().py(), phiFit_Tracks[0].track().pz(), Mass_K);
            ftree->phiFit_Kp_ETA = phiFit_Kp_P4.Eta();
            ftree->phiFit_Kp_PHI = phiFit_Kp_P4.Phi();
            ftree->phiFit_Kp_P = phiFit_Kp_P4.P();
            ftree->phiFit_Kp_PT = phiFit_Kp_P4.Pt();
            ftree->phiFit_Kp_PX = phiFit_Kp_P4.Px();
            ftree->phiFit_Kp_PY = phiFit_Kp_P4.Py();
            ftree->phiFit_Kp_PZ = phiFit_Kp_P4.Pz();

            TLorentzVector phiFit_Km_P4;
            phiFit_Km_P4.SetXYZM(phiFit_Tracks[1].track().px(), phiFit_Tracks[1].track().py(), phiFit_Tracks[1].track().pz(), Mass_K);
            ftree->phiFit_Km_ETA = phiFit_Km_P4.Eta();
            ftree->phiFit_Km_PHI = phiFit_Km_P4.Phi();
            ftree->phiFit_Km_P = phiFit_Km_P4.P();
            ftree->phiFit_Km_PT = phiFit_Km_P4.Pt();
            ftree->phiFit_Km_PX = phiFit_Km_P4.Px();
            ftree->phiFit_Km_PY = phiFit_Km_P4.Py();
            ftree->phiFit_Km_PZ = phiFit_Km_P4.Pz();

            TLorentzVector phi_P4 = phiFit_Kp_P4 + phiFit_Km_P4;
            ftree->phi_ETA = phi_P4.Eta();
            ftree->phi_PHI = phi_P4.Phi();
            ftree->phi_P = phi_P4.P();
            ftree->phi_PT = phi_P4.Pt();
            ftree->phi_PX = phi_P4.Px();
            ftree->phi_PY = phi_P4.Py();
            ftree->phi_PZ = phi_P4.Pz();
            ftree->phi_CHI2 = phi_Vertex.totalChiSquared();
            ftree->phi_M = phi_P4.M(); 
            ftree->phi_ENDVX_X = phi_Vertex.position().x();
            ftree->phi_ENDVX_Y = phi_Vertex.position().y();
            ftree->phi_ENDVX_Z = phi_Vertex.position().z();
            math::XYZPoint phi_endvertex(phi_Vertex.position().x(), phi_Vertex.position().y(), phi_Vertex.position().z());
            
            ftree->Kp_PP = Kp_P4.Vect().Pt(phi_P4.Vect());
            ftree->Kp_PL = Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            ftree->Km_PP = Km_P4.Vect().Pt(phi_P4.Vect());
            ftree->Km_PL = Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            
            ftree->phiFit_Kp_PP = phiFit_Kp_P4.Vect().Pt(phi_P4.Vect());
            ftree->phiFit_Kp_PL = phiFit_Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            ftree->phiFit_Km_PP = phiFit_Km_P4.Vect().Pt(phi_P4.Vect());
            ftree->phiFit_Km_PL = phiFit_Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            
            ftree->dR_Kp_phi = reco::deltaR(ftree->Kp_ETA, ftree->Kp_PHI, ftree->phi_ETA, ftree->phi_PHI);
            ftree->dR_Km_phi = reco::deltaR(ftree->Km_ETA, ftree->Km_PHI, ftree->phi_ETA, ftree->phi_PHI);
            ftree->d_Kp_phi = (Kp_originvertex - phi_endvertex).R();
            ftree->d_Km_phi = (Km_originvertex - phi_endvertex).R();

            if( ftree->phi_CHI2 < 0 ) continue;
            if( ftree->phi_CHI2 > 10 ) continue;

            ftree->num_reco_phi++;
            for(size_t k=0; k<idx_pi_vec.size(); k++){

                if( idx_pi_vec[k] == idx_Kp_vec[i] ) continue;
            
                const auto& pi_PF = (*packedPF)[idx_pi_vec[k]];
                ftree->pi_Reset();
                
                ftree->pi_ETA = pi_PF.eta();
                ftree->pi_PHI = pi_PF.phi();
                ftree->pi_ORIVX_X = pi_PF.vx();
                ftree->pi_ORIVX_Y = pi_PF.vy();
                ftree->pi_ORIVX_Z = pi_PF.vz();
                ftree->pi_P = pi_PF.p();
                ftree->pi_PT = pi_PF.pt();
                ftree->pi_PX = pi_PF.px();
                ftree->pi_PY = pi_PF.py();
                ftree->pi_PZ = pi_PF.pz();
                math::XYZPoint pi_originvertex(pi_PF.vx(), pi_PF.vy(), pi_PF.vz());
                TLorentzVector pi_P4; pi_P4.SetXYZM(pi_PF.px(), pi_PF.py(), pi_PF.pz(), Mass_pi);
                
                ftree->dR_Kp_pi = reco::deltaR(ftree->Kp_ETA, ftree->Kp_PHI, ftree->pi_ETA, ftree->pi_PHI);
                ftree->dR_Km_pi = reco::deltaR(ftree->Km_ETA, ftree->Km_PHI, ftree->pi_ETA, ftree->pi_PHI);
                ftree->dR_phi_pi = reco::deltaR(ftree->pi_ETA, ftree->pi_PHI, ftree->phi_ETA, ftree->phi_PHI);
                ftree->d_Kp_pi = (Kp_originvertex - pi_originvertex).R();
                ftree->d_Km_pi = (Km_originvertex - pi_originvertex).R();
                ftree->d_phi_pi = (phi_endvertex - pi_originvertex).R();
               
                std::cout << 0 << std::endl; 
                if( ftree->dR_Kp_pi > 0.6 ) continue;
                std::cout << 1 << std::endl; 
                if( ftree->dR_Km_pi > 0.6 ) continue;
                std::cout << 2 << std::endl; 
                if( ftree->dR_phi_pi > 1 ) continue;

                std::cout << 3 << std::endl; 
                if( ftree->d_Kp_pi > 0.4 ) continue;
                std::cout << 4 << std::endl; 
                if( ftree->d_Km_pi > 0.4 ) continue;
                std::cout << 5 << std::endl; 
                if( ftree->d_phi_pi > 6 ) continue;

                std::cout << 6 << std::endl; 
                std::vector<reco::TransientTrack> Ds_Tracks = {
                    phiFit_Tracks[0],
                    phiFit_Tracks[1],
                    (*theB).build(pi_PF.pseudoTrack())
                };

                TransientVertex Ds_Vertex = fitter.vertex(Ds_Tracks);

                std::cout << 7 << std::endl; 
                if( !(Ds_Vertex.isValid()) ) continue;
                std::cout << 8 << std::endl; 
                if( !(Ds_Vertex.hasRefittedTracks()) ) continue;

                ftree->Ds_Reset();

                std::vector<reco::TransientTrack> DsFit_Tracks = Ds_Vertex.refittedTracks();
                
                TLorentzVector DsFit_Kp_P4;
                DsFit_Kp_P4.SetXYZM(DsFit_Tracks[0].track().px(), DsFit_Tracks[0].track().py(), DsFit_Tracks[0].track().pz(), Mass_K);
                ftree->DsFit_Kp_ETA = DsFit_Kp_P4.Eta();
                ftree->DsFit_Kp_PHI = DsFit_Kp_P4.Phi();
                ftree->DsFit_Kp_P = DsFit_Kp_P4.P();
                ftree->DsFit_Kp_PT = DsFit_Kp_P4.Pt();
                ftree->DsFit_Kp_PX = DsFit_Kp_P4.Px();
                ftree->DsFit_Kp_PY = DsFit_Kp_P4.Py();
                ftree->DsFit_Kp_PZ = DsFit_Kp_P4.Pz();

                TLorentzVector DsFit_Km_P4;
                DsFit_Km_P4.SetXYZM(DsFit_Tracks[1].track().px(), DsFit_Tracks[1].track().py(), DsFit_Tracks[1].track().pz(), Mass_K);
                ftree->DsFit_Km_ETA = DsFit_Km_P4.Eta();
                ftree->DsFit_Km_PHI = DsFit_Km_P4.Phi();
                ftree->DsFit_Km_P = DsFit_Km_P4.P();
                ftree->DsFit_Km_PT = DsFit_Km_P4.Pt();
                ftree->DsFit_Km_PX = DsFit_Km_P4.Px();
                ftree->DsFit_Km_PY = DsFit_Km_P4.Py();
                ftree->DsFit_Km_PZ = DsFit_Km_P4.Pz();

                TLorentzVector DsFit_phi_P4 = DsFit_Kp_P4 + DsFit_Km_P4;
                ftree->DsFit_phi_ETA = DsFit_phi_P4.Eta();
                ftree->DsFit_phi_PHI = DsFit_phi_P4.Phi();
                ftree->DsFit_phi_P = DsFit_phi_P4.P();
                ftree->DsFit_phi_PT = DsFit_phi_P4.Pt();
                ftree->DsFit_phi_PX = DsFit_phi_P4.Px();
                ftree->DsFit_phi_PY = DsFit_phi_P4.Py();
                ftree->DsFit_phi_PZ = DsFit_phi_P4.Pz();
                ftree->DsFit_phi_M = DsFit_phi_P4.M();

                ftree->DsFit_Kp_PP = DsFit_Kp_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Kp_PL = DsFit_Kp_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();
                ftree->DsFit_Km_PP = DsFit_Km_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Km_PL = DsFit_Km_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();

                TLorentzVector DsFit_pi_P4;
                DsFit_pi_P4.SetXYZM(DsFit_Tracks[2].track().px(), DsFit_Tracks[2].track().py(), DsFit_Tracks[2].track().pz(), Mass_pi);
                ftree->DsFit_pi_ETA = DsFit_pi_P4.Eta();
                ftree->DsFit_pi_PHI = DsFit_pi_P4.Phi();
                ftree->DsFit_pi_P = DsFit_pi_P4.P();
                ftree->DsFit_pi_PT = DsFit_pi_P4.Pt();
                ftree->DsFit_pi_PX = DsFit_pi_P4.Px();
                ftree->DsFit_pi_PY = DsFit_pi_P4.Py();
                ftree->DsFit_pi_PZ = DsFit_pi_P4.Pz();

                TLorentzVector Ds_P4 = DsFit_phi_P4 + DsFit_pi_P4;
                ftree->Ds_ETA = Ds_P4.Eta();
                ftree->Ds_PHI = Ds_P4.Phi();
                ftree->Ds_P = Ds_P4.P();
                ftree->Ds_PT = Ds_P4.Pt();
                ftree->Ds_PX = Ds_P4.Px();
                ftree->Ds_PY = Ds_P4.Py();
                ftree->Ds_PZ = Ds_P4.Pz();
                ftree->Ds_CHI2 = Ds_Vertex.totalChiSquared();
                ftree->Ds_M = Ds_P4.M(); 
                ftree->Ds_ENDVX_X = Ds_Vertex.position().x();
                ftree->Ds_ENDVX_Y = Ds_Vertex.position().y();
                ftree->Ds_ENDVX_Z = Ds_Vertex.position().z();
                math::XYZPoint Ds_endvertex(Ds_Vertex.position().x(), Ds_Vertex.position().y(), Ds_Vertex.position().z());

                TLorentzVector DsFit_Mconstraint_phi_P4;
                DsFit_Mconstraint_phi_P4.SetXYZM(DsFit_phi_P4.Px(), DsFit_phi_P4.Py(), DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector DsFit_Mconstraint_Ds_P4 = DsFit_Mconstraint_phi_P4 + DsFit_pi_P4;
                ftree->DsFit_Mconstraint_Ds_M = DsFit_Mconstraint_Ds_P4.M();

                ftree->phi_PP = phi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->phi_PL = phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
                ftree->pi_PP = pi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->pi_PL = pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
            
                ftree->DsFit_phi_PP = DsFit_phi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->DsFit_phi_PL = DsFit_phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
                ftree->DsFit_pi_PP = DsFit_pi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->DsFit_pi_PL = DsFit_pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
            
                ftree->DsFit_Kp_PP = DsFit_Kp_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Kp_PL = DsFit_Kp_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();
                ftree->DsFit_Km_PP = DsFit_Km_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Km_PL = DsFit_Km_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();
        
                ftree->dR_Kp_Ds = reco::deltaR(ftree->Kp_ETA, ftree->Kp_PHI, ftree->Ds_ETA, ftree->Ds_PHI);
                ftree->dR_Km_Ds = reco::deltaR(ftree->Km_ETA, ftree->Km_PHI, ftree->Ds_ETA, ftree->Ds_PHI);
                ftree->dR_phi_Ds = reco::deltaR(ftree->phi_ETA, ftree->phi_PHI, ftree->Ds_ETA, ftree->Ds_PHI);
                ftree->dR_pi_Ds = reco::deltaR(ftree->pi_ETA, ftree->pi_PHI, ftree->Ds_ETA, ftree->Ds_PHI);
                ftree->d_Kp_Ds = (Kp_originvertex - Ds_endvertex).R();
                ftree->d_Km_Ds = (Km_originvertex - Ds_endvertex).R();
                ftree->d_phi_Ds = (phi_endvertex - Ds_endvertex).R();
                ftree->d_pi_Ds = (pi_originvertex - Ds_endvertex).R();
    
                std::cout << 599 << std::endl; 
                if( ftree->Ds_CHI2 > 25 ) continue;
                std::cout << 7232 << std::endl; 
                if( ftree->Ds_CHI2 < 0 ) continue;
                std::cout << 48263 << std::endl; 
                if( ftree->dR_phi_Ds > 0.4 ) continue;
                std::cout << 7723923 << std::endl; 
                if( ftree->d_phi_Ds > 4 ) continue;
                std::cout << 7823924 << std::endl; 

                ftree->num_reco_Ds++;

                ftree->Fill_Vector(); 
            }
        }
    }

    ftree->tree->Fill();

    return;
}

void PreliminaryCut::beginJob() {}

void PreliminaryCut::endJob() {}

void PreliminaryCut::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("PreliminaryCut", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PreliminaryCut);

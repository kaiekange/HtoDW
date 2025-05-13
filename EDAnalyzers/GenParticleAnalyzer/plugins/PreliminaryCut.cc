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
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


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

        void Kp_CLEAR();
        void Km_CLEAR();
        void pi_CLEAR();
        void phi_CLEAR();
        void Ds_CLEAR();

        edm::Service<TFileService> fs;
        TTree *tree_{nullptr};

        std::vector<int> idx_Kp_vec;
        std::vector<int> idx_Km_vec;
        std::vector<int> idx_pi_vec;

        double Kp_ETA;
        double Kp_PHI;
        double Kp_ORIVX_X;
        double Kp_ORIVX_Y;
        double Kp_ORIVX_Z;
        double Kp_P;
        double Kp_PT;
        double Kp_PX;
        double Kp_PY;
        double Kp_PZ;
        double Kp_PP;
        double Kp_PL;

        double Km_ETA;
        double Km_PHI;
        double Km_ORIVX_X;
        double Km_ORIVX_Y;
        double Km_ORIVX_Z;
        double Km_P;
        double Km_PT;
        double Km_PX;
        double Km_PY;
        double Km_PZ;
        double Km_PP;
        double Km_PL;

        double pi_ETA;
        double pi_PHI;
        double pi_ORIVX_X;
        double pi_ORIVX_Y;
        double pi_ORIVX_Z;
        double pi_P;
        double pi_PT;
        double pi_PX;
        double pi_PY;
        double pi_PZ;
        double pi_PP;
        double pi_PL;

        double phi_ETA;
        double phi_PHI;
        double phi_P;
        double phi_PT;
        double phi_PX;
        double phi_PY;
        double phi_PZ;
        double phi_PP;
        double phi_PL;
        double phi_CHI2;
        double phi_M;
        double phi_ENDVX_X;
        double phi_ENDVX_Y;
        double phi_ENDVX_Z;

        double Ds_ETA;
        double Ds_PHI;
        double Ds_P;
        double Ds_PT;
        double Ds_PX;
        double Ds_PY;
        double Ds_PZ;
        double Ds_CHI2;
        double Ds_M;
        double Ds_ENDVX_X;
        double Ds_ENDVX_Y;
        double Ds_ENDVX_Z;

        double dR_Kp_Km;
        double dR_Kp_pi;
        double dR_Kp_phi;
        double dR_Kp_Ds;
        double dR_Km_pi;
        double dR_Km_phi;
        double dR_Km_Ds;
        double dR_pi_phi;
        double dR_pi_Ds;
        double dR_phi_Ds;
        
        double d_Kp_Km;
        double d_Kp_pi;
        double d_Kp_phi;
        double d_Kp_Ds;
        double d_Km_pi;
        double d_Km_phi;
        double d_Km_Ds;
        double d_pi_phi;
        double d_pi_Ds;
        double d_phi_Ds;

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

PreliminaryCut::PreliminaryCut(const edm::ParameterSet& iConfig) :
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
}

PreliminaryCut::~PreliminaryCut() {}

void PreliminaryCut::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    idx_Kp_vec.clear();
    idx_Km_vec.clear();
    idx_pi_vec.clear();

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

        int idx_Kp = idx_Kp_vec[i];

        Kp_CLEAR();
        Kp_ETA = (*packedPF)[idx_Kp].eta();
        Kp_PHI = (*packedPF)[idx_Kp].phi();
        Kp_ORIVX_X = (*packedPF)[idx_Kp].vx();
        Kp_ORIVX_Y = (*packedPF)[idx_Kp].vy();
        Kp_ORIVX_Z = (*packedPF)[idx_Kp].vz();
        Kp_P = (*packedPF)[idx_Kp].p();
        Kp_PT = (*packedPF)[idx_Kp].pt();
        Kp_PX = (*packedPF)[idx_Kp].px();
        Kp_PY = (*packedPF)[idx_Kp].py();
        Kp_PZ = (*packedPF)[idx_Kp].pz();
        math::XYZPoint Kp_originvertex(Kp_ORIVX_X, Kp_ORIVX_Y, Kp_ORIVX_Z);

        for(size_t j=0; j<idx_Km_vec.size(); j++){
            int idx_Km = idx_Km_vec[j];

            Km_CLEAR();
            Km_ETA = (*packedPF)[idx_Km].eta();
            Km_PHI = (*packedPF)[idx_Km].phi();
            Km_ORIVX_X = (*packedPF)[idx_Km].vx();
            Km_ORIVX_Y = (*packedPF)[idx_Km].vy();
            Km_ORIVX_Z = (*packedPF)[idx_Km].vz();
            Km_P = (*packedPF)[idx_Km].p();
            Km_PT = (*packedPF)[idx_Km].pt();
            Km_PX = (*packedPF)[idx_Km].px();
            Km_PY = (*packedPF)[idx_Km].py();
            Km_PZ = (*packedPF)[idx_Km].pz();
            math::XYZPoint Km_originvertex(Km_ORIVX_X, Km_ORIVX_Y, Km_ORIVX_Z);
            
            dR_Kp_Km = reco::deltaR((*packedPF)[idx_Kp], (*packedPF)[idx_Km]);
            d_Kp_Km = (Kp_originvertex - Km_originvertex).R();
            
            if( abs(Kp_PT - Km_PT) > 20) continue;
            if(dR_Kp_Km > 0.15) continue;
            if(d_Kp_Km > 0.2) continue;

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*theB).build((*packedPF)[idx_Kp].pseudoTrack()),
                (*theB).build((*packedPF)[idx_Km].pseudoTrack())
            };
            TransientVertex phi_Vertex = fitter.vertex(phi_Tracks);

            if(!(phi_Vertex.isValid())) continue;
            if(!(phi_Vertex.hasRefittedTracks())) continue;
            
            std::vector<reco::TransientTrack> refitted_phi_Tracks = phi_Vertex.refittedTracks();

            TLorentzVector p4Kp; p4Kp.SetXYZM(refitted_phi_Tracks[0].track().px(), refitted_phi_Tracks[0].track().py(), refitted_phi_Tracks[0].track().pz(), Mass_K);
            TLorentzVector p4Km; p4Km.SetXYZM(refitted_phi_Tracks[1].track().px(), refitted_phi_Tracks[1].track().py(), refitted_phi_Tracks[1].track().pz(), Mass_K);
            TLorentzVector p4phi = p4Kp + p4Km;
            
            phi_CLEAR();
            phi_ETA = p4phi.Eta();
            phi_PHI = p4phi.Phi();
            phi_P = p4phi.P();
            phi_PT = p4phi.Pt();
            phi_PX = p4phi.Px();
            phi_PY = p4phi.Py();
            phi_PZ = p4phi.Pz();
            phi_CHI2 = phi_Vertex.totalChiSquared();
            phi_M = p4phi.M(); 
            phi_ENDVX_X = phi_Vertex.position().x();
            phi_ENDVX_Y = phi_Vertex.position().y();
            phi_ENDVX_Z = phi_Vertex.position().z();
            math::XYZPoint phi_endvertex(phi_ENDVX_X, phi_ENDVX_Y, phi_ENDVX_Z);
            
            Kp_PP = p4Kp.Vect().Pt(p4phi.Vect());
            Kp_PL = p4Kp.Vect().Dot(p4phi.Vect())/phi_P;
            Km_PP = p4Km.Vect().Pt(p4phi.Vect());
            Km_PL = p4Km.Vect().Dot(p4phi.Vect())/phi_P;
            
            dR_Kp_phi = reco::deltaR(Kp_ETA, Kp_PHI, phi_ETA, phi_PHI);
            dR_Km_phi = reco::deltaR(Km_ETA, Km_PHI, phi_ETA, phi_PHI);
            d_Kp_phi = (Kp_originvertex - phi_endvertex).R();
            d_Km_phi = (Km_originvertex - phi_endvertex).R();

            if(phi_CHI2 > 10) continue;
            if(phi_CHI2 < 0) continue;

            for(size_t k=0; k<idx_pi_vec.size(); k++){
                int idx_pi = idx_pi_vec[k];
                if(idx_pi == idx_Kp) continue;
               
                pi_CLEAR(); 
                pi_ETA = (*packedPF)[idx_pi].eta();
                pi_PHI = (*packedPF)[idx_pi].phi();
                pi_ORIVX_X = (*packedPF)[idx_pi].vx();
                pi_ORIVX_Y = (*packedPF)[idx_pi].vy();
                pi_ORIVX_Z = (*packedPF)[idx_pi].vz();
                pi_P = (*packedPF)[idx_pi].p();
                pi_PT = (*packedPF)[idx_pi].pt();
                pi_PX = (*packedPF)[idx_pi].px();
                pi_PY = (*packedPF)[idx_pi].py();
                pi_PZ = (*packedPF)[idx_pi].pz();
                math::XYZPoint pi_originvertex(pi_ORIVX_X, pi_ORIVX_Y, pi_ORIVX_Z);
                
                dR_Kp_pi = reco::deltaR(Kp_ETA, Kp_PHI, pi_ETA, pi_PHI);
                dR_Km_pi = reco::deltaR(Km_ETA, Km_PHI, pi_ETA, pi_PHI);
                dR_pi_phi = reco::deltaR(pi_ETA, pi_PHI, phi_ETA, phi_PHI);
                d_Kp_pi = (Kp_originvertex - pi_originvertex).R();
                d_Km_pi = (Km_originvertex - pi_originvertex).R();
                d_pi_phi = (pi_originvertex - phi_endvertex).R();
                
                if(dR_Kp_pi > 0.6) continue;
                if(dR_Km_pi > 0.6) continue;
                if(dR_pi_phi > 1) continue;

                if(d_Kp_pi > 0.4) continue;
                if(d_Km_pi > 0.4) continue;
                if(d_pi_phi > 6) continue;

                std::vector<reco::TransientTrack> Ds_Tracks = {
                    refitted_phi_Tracks[0],
                    refitted_phi_Tracks[1],
                    (*theB).build((*packedPF)[idx_pi].pseudoTrack())
                };

                TransientVertex Ds_Vertex = fitter.vertex(Ds_Tracks);

                if(!(Ds_Vertex.isValid())) continue;
                if(!(Ds_Vertex.hasRefittedTracks())) continue;

                std::vector<reco::TransientTrack> refitted_Ds_Tracks = Ds_Vertex.refittedTracks();
                
                TLorentzVector refitted_p4Kp; refitted_p4Kp.SetXYZM(refitted_Ds_Tracks[0].track().px(), refitted_Ds_Tracks[0].track().py(), refitted_Ds_Tracks[0].track().pz(), Mass_K);
                TLorentzVector refitted_p4Km; refitted_p4Km.SetXYZM(refitted_Ds_Tracks[1].track().px(), refitted_Ds_Tracks[1].track().py(), refitted_Ds_Tracks[1].track().pz(), Mass_K);
                TLorentzVector refitted_p4phi; refitted_p4phi = refitted_p4Kp + refitted_p4Km;

                TLorentzVector p4pi; p4pi.SetXYZM(refitted_Ds_Tracks[2].track().px(), refitted_Ds_Tracks[2].track().py(), refitted_Ds_Tracks[2].track().pz(), Mass_pi);
                TLorentzVector p4Ds; p4Ds = refitted_p4Kp + refitted_p4Km + p4pi;

                Ds_CLEAR();
                Ds_ETA = p4Ds.Eta();
                Ds_PHI = p4Ds.Phi();
                Ds_P = p4Ds.P();
                Ds_PT = p4Ds.Pt();
                Ds_PX = p4Ds.Px();
                Ds_PY = p4Ds.Py();
                Ds_PZ = p4Ds.Pz();
                Ds_CHI2 = Ds_Vertex.totalChiSquared();
                Ds_M = p4Ds.M(); 
                Ds_ENDVX_X = Ds_Vertex.position().x();
                Ds_ENDVX_Y = Ds_Vertex.position().y();
                Ds_ENDVX_Z = Ds_Vertex.position().z();
                math::XYZPoint Ds_endvertex(Ds_ENDVX_X, Ds_ENDVX_Y, Ds_ENDVX_Z);
                
                pi_PP = p4pi.Vect().Pt(p4Ds.Vect());
                pi_PL = p4pi.Vect().Dot(p4Ds.Vect())/Ds_P;
                phi_PP = refitted_p4phi.Vect().Pt(p4Ds.Vect());
                phi_PL = refitted_p4phi.Vect().Dot(p4Ds.Vect())/Ds_P;
                
                dR_Kp_Ds = reco::deltaR(Kp_ETA, Kp_PHI, Ds_ETA, Ds_PHI);
                dR_Km_Ds = reco::deltaR(Km_ETA, Km_PHI, Ds_ETA, Ds_PHI);
                dR_pi_Ds = reco::deltaR(pi_ETA, pi_PHI, Ds_ETA, Ds_PHI);
                dR_phi_Ds = reco::deltaR(phi_ETA, phi_PHI, Ds_ETA, Ds_PHI);
                d_Kp_Ds = (Kp_originvertex - Ds_endvertex).R();
                d_Km_Ds = (Km_originvertex - Ds_endvertex).R();
                d_pi_Ds = (pi_originvertex - Ds_endvertex).R();
                d_phi_Ds = (phi_endvertex - Ds_endvertex).R();

                if(Ds_CHI2 > 25) continue;
                if(Ds_CHI2 < 0) continue;
                if(dR_phi_Ds > 0.4) continue;
                if(d_phi_Ds > 4) continue;
                
                tree_->Fill();
            }
        }
    }

    return;
}

void PreliminaryCut::Kp_CLEAR()
{
    Kp_ETA = -777;
    Kp_PHI = -777;
    Kp_ORIVX_X = -777;
    Kp_ORIVX_Y = -777;
    Kp_ORIVX_Z = -777;
    Kp_P = -777;
    Kp_PT = -777;
    Kp_PX = -777;
    Kp_PY = -777;
    Kp_PZ = -777;
}
void PreliminaryCut::Km_CLEAR()
{
    Km_ETA = -777;
    Km_PHI = -777;
    Km_ORIVX_X = -777;
    Km_ORIVX_Y = -777;
    Km_ORIVX_Z = -777;
    Km_P = -777;
    Km_PT = -777;
    Km_PX = -777;
    Km_PY = -777;
    Km_PZ = -777;
    dR_Kp_Km = -777;
    d_Kp_Km = -777;
}
void PreliminaryCut::phi_CLEAR()
{
    phi_ETA = -777;
    phi_PHI = -777;
    phi_P = -777;
    phi_PT = -777;
    phi_PX = -777;
    phi_PY = -777;
    phi_PZ = -777;
    phi_CHI2 = -777;
    phi_M = -777;
    phi_ENDVX_X = -777;
    phi_ENDVX_Y = -777;
    phi_ENDVX_Z = -777;
    Kp_PP = -777;
    Kp_PL = -777;
    Km_PP = -777;
    Km_PL = -777;
    dR_Kp_phi = -777;
    dR_Km_phi = -777;
    d_Kp_phi = -777;
    d_Km_phi = -777;
}
void PreliminaryCut::pi_CLEAR()
{
    pi_ETA = -777;
    pi_PHI = -777;
    pi_ORIVX_X = -777;
    pi_ORIVX_Y = -777;
    pi_ORIVX_Z = -777;
    pi_P = -777;
    pi_PT = -777;
    pi_PX = -777;
    pi_PY = -777;
    pi_PZ = -777;
    dR_Kp_pi = -777;
    dR_Km_pi = -777;
    dR_pi_phi = -777;
    d_Kp_pi = -777;
    d_Km_pi = -777;
    d_pi_phi = -777;
}
void PreliminaryCut::Ds_CLEAR()
{
    Ds_ETA = -777;
    Ds_PHI = -777;
    Ds_P = -777;
    Ds_PT = -777;
    Ds_PX = -777;
    Ds_PY = -777;
    Ds_PZ = -777;
    Ds_CHI2 = -777;
    Ds_M = -777;
    Ds_ENDVX_X = -777;
    Ds_ENDVX_Y = -777;
    Ds_ENDVX_Z = -777;
    pi_PP = -777;
    pi_PL = -777;
    phi_PP = -777;
    phi_PL = -777;
    dR_Kp_Ds = -777;
    dR_Km_Ds = -777;
    dR_pi_Ds = -777;
    dR_phi_Ds = -777;
    d_Kp_Ds = -777;
    d_Km_Ds = -777;
    d_pi_Ds = -777;
    d_phi_Ds = -777;
}

void PreliminaryCut::beginJob()
{
    tree_ = fs->make<TTree>("Events", "Events");

    tree_->Branch("Kp_ETA", &Kp_ETA);
    tree_->Branch("Kp_PHI", &Kp_PHI);
    tree_->Branch("Kp_ORIVX_X", &Kp_ORIVX_X);
    tree_->Branch("Kp_ORIVX_Y", &Kp_ORIVX_Y);
    tree_->Branch("Kp_ORIVX_Z", &Kp_ORIVX_Z);
    tree_->Branch("Kp_P", &Kp_P);
    tree_->Branch("Kp_PT", &Kp_PT);
    tree_->Branch("Kp_PX", &Kp_PX);
    tree_->Branch("Kp_PY", &Kp_PY);
    tree_->Branch("Kp_PZ", &Kp_PZ);
    tree_->Branch("Kp_PP", &Kp_PP);
    tree_->Branch("Kp_PL", &Kp_PL);

    tree_->Branch("Km_ETA", &Km_ETA);
    tree_->Branch("Km_PHI", &Km_PHI);
    tree_->Branch("Km_ORIVX_X", &Km_ORIVX_X);
    tree_->Branch("Km_ORIVX_Y", &Km_ORIVX_Y);
    tree_->Branch("Km_ORIVX_Z", &Km_ORIVX_Z);
    tree_->Branch("Km_P", &Km_P);
    tree_->Branch("Km_PT", &Km_PT);
    tree_->Branch("Km_PX", &Km_PX);
    tree_->Branch("Km_PY", &Km_PY);
    tree_->Branch("Km_PZ", &Km_PZ);
    tree_->Branch("Km_PP", &Km_PP);
    tree_->Branch("Km_PL", &Km_PL);

    tree_->Branch("pi_ETA", &pi_ETA);
    tree_->Branch("pi_PHI", &pi_PHI);
    tree_->Branch("pi_ORIVX_X", &pi_ORIVX_X);
    tree_->Branch("pi_ORIVX_Y", &pi_ORIVX_Y);
    tree_->Branch("pi_ORIVX_Z", &pi_ORIVX_Z);
    tree_->Branch("pi_P", &pi_P);
    tree_->Branch("pi_PT", &pi_PT);
    tree_->Branch("pi_PX", &pi_PX);
    tree_->Branch("pi_PY", &pi_PY);
    tree_->Branch("pi_PZ", &pi_PZ);
    tree_->Branch("pi_PP", &pi_PP);
    tree_->Branch("pi_PL", &pi_PL);

    tree_->Branch("phi_ETA", &phi_ETA);
    tree_->Branch("phi_PHI", &phi_PHI);
    tree_->Branch("phi_P", &phi_P);
    tree_->Branch("phi_PT", &phi_PT);
    tree_->Branch("phi_PX", &phi_PX);
    tree_->Branch("phi_PY", &phi_PY);
    tree_->Branch("phi_PZ", &phi_PZ);
    tree_->Branch("phi_PP", &phi_PP);
    tree_->Branch("phi_PL", &phi_PL);
    tree_->Branch("phi_CHI2", &phi_CHI2);
    tree_->Branch("phi_M", &phi_M);
    tree_->Branch("phi_ENDVX_X", &phi_ENDVX_X);
    tree_->Branch("phi_ENDVX_Y", &phi_ENDVX_Y);
    tree_->Branch("phi_ENDVX_Z", &phi_ENDVX_Z);

    tree_->Branch("Ds_ETA", &Ds_ETA);
    tree_->Branch("Ds_PHI", &Ds_PHI);
    tree_->Branch("Ds_P", &Ds_P);
    tree_->Branch("Ds_PT", &Ds_PT);
    tree_->Branch("Ds_PX", &Ds_PX);
    tree_->Branch("Ds_PY", &Ds_PY);
    tree_->Branch("Ds_PZ", &Ds_PZ);
    tree_->Branch("Ds_CHI2", &Ds_CHI2);
    tree_->Branch("Ds_M", &Ds_M);
    tree_->Branch("Ds_ENDVX_X", &Ds_ENDVX_X);
    tree_->Branch("Ds_ENDVX_Y", &Ds_ENDVX_Y);
    tree_->Branch("Ds_ENDVX_Z", &Ds_ENDVX_Z);

    tree_->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree_->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree_->Branch("dR_Kp_phi", &dR_Kp_phi);
    tree_->Branch("dR_Kp_Ds", &dR_Kp_Ds);
    tree_->Branch("dR_Km_pi", &dR_Km_pi);
    tree_->Branch("dR_Km_phi", &dR_Km_phi);
    tree_->Branch("dR_Km_Ds", &dR_Km_Ds);
    tree_->Branch("dR_pi_phi", &dR_pi_phi);
    tree_->Branch("dR_pi_Ds", &dR_pi_Ds);
    tree_->Branch("dR_phi_Ds", &dR_phi_Ds);

    tree_->Branch("d_Kp_Km", &d_Kp_Km);
    tree_->Branch("d_Kp_pi", &d_Kp_pi);
    tree_->Branch("d_Kp_phi", &d_Kp_phi);
    tree_->Branch("d_Kp_Ds", &d_Kp_Ds);
    tree_->Branch("d_Km_pi", &d_Km_pi);
    tree_->Branch("d_Km_phi", &d_Km_phi);
    tree_->Branch("d_Km_Ds", &d_Km_Ds);
    tree_->Branch("d_pi_phi", &d_pi_phi);
    tree_->Branch("d_pi_Ds", &d_pi_Ds);
    tree_->Branch("d_phi_Ds", &d_phi_Ds);
}

void PreliminaryCut::endJob() {}

void PreliminaryCut::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("PreliminaryCut", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PreliminaryCut);

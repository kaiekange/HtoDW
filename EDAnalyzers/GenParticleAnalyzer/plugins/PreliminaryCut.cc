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

        void vector_clear();

        edm::Service<TFileService> fs;
        TTree *tree_{nullptr};

        std::vector<int> idx_Kp_vec;
        std::vector<int> idx_Km_vec;
        std::vector<int> idx_pi_vec;

        int n_phi;
        int n_Ds;

        double pt_Kp;
        double pt_Km;
        double pt_pi;
        double p_Kp;
        double p_Km;
        double p_pi;

        double dR_Kp_Km;
        double dR_Kp_pi;
        double dR_Km_pi;
        double d_Kp_Km;
        double d_Kp_pi;
        double d_Km_pi;
        
        double chi2_phi;
        double pt_phi;
        double p_phi;
        double m_phi;
        double dR_phi_pi;
        double d_phi_pi;
    
        double chi2_Ds;
        double pt_Ds;
        double p_Ds;
        double m_Ds;
        double dR_phi_Ds;
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

    n_phi = 0;
    n_Ds = 0;

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

        pt_Kp = -777;
        p_Kp = -777;

        pt_Kp = (*packedPF)[idx_Kp].pt();
        p_Kp = (*packedPF)[idx_Kp].p();

        for(size_t j=0; j<idx_Km_vec.size(); j++){
            int idx_Km = idx_Km_vec[j];

            pt_Km = -777;
            p_Km = -777;
            dR_Kp_Km = -777;
            d_Kp_Km = -777;
            chi2_phi = -777;
            pt_phi = -777;
            p_phi = -777;
            m_phi = -777;
            
            pt_Km = (*packedPF)[idx_Km].pt();
            p_Km = (*packedPF)[idx_Km].p();

            if( abs(pt_Kp - pt_Km) > 20) continue;
            
            dR_Kp_Km = reco::deltaR((*packedPF)[idx_Kp], (*packedPF)[idx_Km]);
            if(dR_Kp_Km > 0.15) continue;
            d_Kp_Km = ((*packedPF)[idx_Kp].vertex() - (*packedPF)[idx_Km].vertex()).R();
            if(d_Kp_Km > 0.2) continue;
            
            std::vector<reco::TransientTrack> phiTracks = {
                (*theB).build((*packedPF)[idx_Kp].pseudoTrack()),
                (*theB).build((*packedPF)[idx_Km].pseudoTrack())
            };
            TransientVertex phiVertex = fitter.vertex(phiTracks);

            if(!(phiVertex.isValid())) continue;

            chi2_phi = phiVertex.totalChiSquared();
            if(chi2_phi > 10) continue;

            if(!(phiVertex.hasRefittedTracks())) continue;
            
            std::vector<reco::TransientTrack> refitted_phiTracks = phiVertex.refittedTracks();

            TLorentzVector p4Kp; p4Kp.SetXYZM(refitted_phiTracks[0].track().px(), refitted_phiTracks[0].track().py(), refitted_phiTracks[0].track().pz(), Mass_K);
            TLorentzVector p4Km; p4Km.SetXYZM(refitted_phiTracks[1].track().px(), refitted_phiTracks[1].track().py(), refitted_phiTracks[1].track().pz(), Mass_K);
            TLorentzVector p4phi; p4phi = p4Kp + p4Km;
            
            pt_phi = p4phi.Pt();
            p_phi = p4phi.P();
            m_phi = p4phi.M();

            n_phi++;

            for(size_t k=0; k<idx_pi_vec.size(); k++){
                
                pt_pi = -777;
                p_pi = -777;
                dR_Kp_pi = -777;
                d_Kp_pi = -777;
                dR_Km_pi = -777;
                d_Km_pi = -777;
                dR_phi_pi = -777;
                d_phi_pi = -777;
                dR_phi_Ds = -777;
                d_phi_Ds = -777;
                chi2_Ds = -777;
                pt_Ds = -777;
                p_Ds = -777;
                m_Ds = -777;
                
                int idx_pi = idx_pi_vec[k];
                if(idx_pi == idx_Kp) continue;
                
                pt_pi = (*packedPF)[idx_pi].pt();
                p_pi = (*packedPF)[idx_pi].p();

                dR_Kp_pi = reco::deltaR((*packedPF)[idx_Kp], (*packedPF)[idx_pi]);
                if(dR_Kp_pi > 0.6) continue;
                d_Kp_pi = ((*packedPF)[idx_Kp].vertex() - (*packedPF)[idx_pi].vertex()).R();
                if(d_Kp_pi > 0.4) continue;
                dR_Km_pi = reco::deltaR((*packedPF)[idx_Km], (*packedPF)[idx_pi]);
                if(dR_Km_pi > 0.6) continue;
                d_Km_pi = ((*packedPF)[idx_Km].vertex() - (*packedPF)[idx_pi].vertex()).R();
                if(d_Km_pi > 0.4) continue;

                dR_phi_pi = reco::deltaR(p4phi.Eta(), p4phi.Phi(), (*packedPF)[idx_pi].eta(), (*packedPF)[idx_pi].phi());
                if(dR_phi_pi > 1) continue;
                math::XYZPoint vtx_phi(phiVertex.position().x(), phiVertex.position().y(), phiVertex.position().z());
                math::XYZPoint vtx_pi = (*packedPF)[idx_pi].vertex();
                d_phi_pi = (vtx_phi - vtx_pi).R();
                if(d_phi_pi > 6) continue;


                std::vector<reco::TransientTrack> DsTracks = {
                    refitted_phiTracks[0],
                    refitted_phiTracks[1],
                    (*theB).build((*packedPF)[idx_pi].pseudoTrack())
                };

                TransientVertex DsVertex = fitter.vertex(DsTracks);

                if(!(DsVertex.isValid())) continue;
                chi2_Ds = DsVertex.totalChiSquared();
                if(chi2_Ds > 25) continue;
                if(!(DsVertex.hasRefittedTracks())) continue;

                std::vector<reco::TransientTrack> refitted_DsTracks = DsVertex.refittedTracks();
                
                TLorentzVector refitted_p4Kp; refitted_p4Kp.SetXYZM(refitted_DsTracks[0].track().px(), refitted_DsTracks[0].track().py(), refitted_DsTracks[0].track().pz(), Mass_K);
                TLorentzVector refitted_p4Km; refitted_p4Km.SetXYZM(refitted_DsTracks[1].track().px(), refitted_DsTracks[1].track().py(), refitted_DsTracks[1].track().pz(), Mass_K);
                TLorentzVector refitted_p4phi; refitted_p4phi = refitted_p4Kp + refitted_p4Km;

                TLorentzVector p4pi; p4pi.SetXYZM(refitted_DsTracks[2].track().px(), refitted_DsTracks[2].track().py(), refitted_DsTracks[2].track().pz(), Mass_pi);
                TLorentzVector p4Ds; p4Ds = refitted_p4Kp + refitted_p4Km + p4pi;

                dR_phi_Ds = reco::deltaR(refitted_p4phi.Eta(), refitted_p4phi.Phi(), p4Ds.Eta(), p4Ds.Phi());
                if(dR_phi_Ds > 0.4) continue;
                math::XYZPoint vtx_Ds(DsVertex.position().x(), DsVertex.position().y(), DsVertex.position().z());
                d_phi_Ds = (vtx_phi - vtx_Ds).R();
                if(d_phi_Ds > 4) continue;
                
                pt_Ds = p4Ds.Pt();
                p_Ds = p4Ds.P();
                m_Ds = p4Ds.M();
                
                n_Ds++;
                tree_->Fill();
            }
        }
    }

    return;
}

void PreliminaryCut::beginJob()
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("Events", "Events");
    
    tree_->Branch("pt_Kp", &pt_Kp);
    tree_->Branch("pt_Km", &pt_Km);
    tree_->Branch("pt_pi", &pt_pi);
    tree_->Branch("p_Kp", &p_Kp);
    tree_->Branch("p_Km", &p_Km);
    tree_->Branch("p_pi", &p_pi);

    tree_->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree_->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree_->Branch("dR_Km_pi", &dR_Km_pi);
    tree_->Branch("d_Kp_Km", &d_Kp_Km);
    tree_->Branch("d_Kp_pi", &d_Kp_pi);
    tree_->Branch("d_Km_pi", &d_Km_pi);

    tree_->Branch("chi2_phi", &chi2_phi);
    tree_->Branch("pt_phi", &pt_phi);
    tree_->Branch("p_phi", &p_phi);
    tree_->Branch("m_phi", &m_phi);
    tree_->Branch("dR_phi_pi", &dR_phi_pi);
    tree_->Branch("d_phi_pi", &d_phi_pi);
    tree_->Branch("n_phi", &n_phi);

    tree_->Branch("chi2_Ds", &chi2_Ds);
    tree_->Branch("pt_Ds", &pt_Ds);
    tree_->Branch("p_Ds", &p_Ds);
    tree_->Branch("m_Ds", &m_Ds);
    tree_->Branch("dR_phi_Ds", &dR_phi_Ds);
    tree_->Branch("d_phi_Ds", &d_phi_Ds);
    tree_->Branch("n_Ds", &n_Ds);
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

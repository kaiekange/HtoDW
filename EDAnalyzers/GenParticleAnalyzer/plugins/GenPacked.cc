// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenPacked
//
/**\class GenPacked GenPacked.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenPacked.cc

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

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


class GenPacked : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenPacked(const edm::ParameterSet&);
        ~GenPacked();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;
        float getDeltaR(float eta1, float phi1, float eta2, float phi2);

        edm::Service<TFileService> myfile;
        TTree *mytree{};
        std::vector<double> all_pt, all_p;
        std::vector<double> Gen_p_Kp, Gen_p_Km, Gen_p_pi;
        std::vector<double> Gen_pt_Kp, Gen_pt_Km, Gen_pt_pi;
        std::vector<double> Gen_dR_Kp_Km, Gen_dR_Kp_pi, Gen_dR_Km_pi;
        std::vector<double> p_Kp, p_Km, p_pi;
        std::vector<double> pt_Kp, pt_Km, pt_pi;
        std::vector<double> dR_Kp_Km, dR_Kp_pi, dR_Km_pi;
        std::vector<double> chi2_phi, m_phi, chi2_Ds, m_Ds;
        std::vector<double> d_Kp_Km, d_Kp_pi, d_Km_pi, d_phi_pi, d_phi_Ds;
        std::vector<float> all_dR_Kp_vec, all_dR_Km_vec, all_dR_pi_vec;
        std::vector<float> dR_Kp_vec, dR_Km_vec, dR_pi_vec;

        const double mK = 0.493677;
        const double mpi = 0.13957039;
};

GenPacked::GenPacked(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    mytree = myfile->make<TTree>("Events", "Events");
    mytree->Branch("Gen_pt_Kp", &Gen_pt_Kp);
    mytree->Branch("Gen_pt_Km", &Gen_pt_Km);
    mytree->Branch("Gen_pt_pi", &Gen_pt_pi);
    mytree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
    mytree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
    mytree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);
    mytree->Branch("pt_Kp", &pt_Kp);
    mytree->Branch("pt_Km", &pt_Km);
    mytree->Branch("pt_pi", &pt_pi);
    mytree->Branch("dR_Kp_Km", &dR_Kp_Km);
    mytree->Branch("dR_Kp_pi", &dR_Kp_pi);
    mytree->Branch("dR_Km_pi", &dR_Km_pi);
    mytree->Branch("chi2_phi", &chi2_phi);
    mytree->Branch("chi2_Ds", &chi2_Ds);
    mytree->Branch("m_phi", &m_phi);
    mytree->Branch("m_Ds", &m_Ds);
    mytree->Branch("d_Kp_Km", &d_Kp_Km);
    mytree->Branch("d_Kp_pi", &d_Kp_pi);
    mytree->Branch("d_Km_pi", &d_Km_pi);
    mytree->Branch("d_phi_pi", &d_phi_pi);
    mytree->Branch("d_phi_Ds", &d_phi_Ds);
    /* mytree->Branch("all_dR_Kp_vec", &all_dR_Kp_vec); */
    /* mytree->Branch("all_dR_Km_vec", &all_dR_Km_vec); */
    /* mytree->Branch("all_dR_pi_vec", &all_dR_pi_vec); */
    /* mytree->Branch("dR_Kp_vec", &dR_Kp_vec); */
    /* mytree->Branch("dR_Km_vec", &dR_Km_vec); */
    /* mytree->Branch("dR_pi_vec", &dR_pi_vec); */
}

GenPacked::~GenPacked() {}

void GenPacked::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    all_pt.clear();
    all_p.clear();
    Gen_p_Kp.clear();
    Gen_p_Km.clear();
    Gen_p_pi.clear();
    Gen_pt_Kp.clear();
    Gen_pt_Km.clear();
    Gen_pt_pi.clear();
    Gen_dR_Kp_Km.clear();
    Gen_dR_Kp_pi.clear();
    Gen_dR_Km_pi.clear();
    p_Kp.clear();
    p_Km.clear();
    p_pi.clear();
    pt_Kp.clear();
    pt_Km.clear();
    pt_pi.clear();
    dR_Kp_Km.clear();
    dR_Kp_pi.clear();
    dR_Km_pi.clear();
    chi2_phi.clear();
    chi2_Ds.clear();
    m_phi.clear();
    m_Ds.clear();
    d_Kp_Km.clear();
    d_Kp_pi.clear();
    d_Km_pi.clear();
    d_phi_pi.clear();
    d_phi_Ds.clear();
    all_dR_Kp_vec.clear();
    all_dR_Km_vec.clear();
    all_dR_pi_vec.clear();
    dR_Kp_vec.clear();
    dR_Km_vec.clear();
    dR_pi_vec.clear();

    edm::Handle<reco::GenParticleCollection> prunedGen;
    iEvent.getByToken(prunedGenToken, prunedGen);

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    float phi_Kp=100, eta_Kp=100, phi_Km=100, eta_Km=100, phi_pi=100, eta_pi=100;
    int ctr_Kp=0, ctr_Km=0, ctr_pi=0; 

    for(size_t i=0; i<prunedGen->size(); ++i){
        if( (*prunedGen)[i].pdgId()==321 ){
            if(hasAncestor((*prunedGen)[i], 333, -321) && hasAncestor((*prunedGen)[i], 431, 211) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_Kp++;
                phi_Kp = (*prunedGen)[i].phi();
                eta_Kp = (*prunedGen)[i].eta();
                Gen_pt_Kp.push_back((*prunedGen)[i].pt());
                Gen_p_Kp.push_back((*prunedGen)[i].p());
            }
        }
        else if( (*prunedGen)[i].pdgId()==-321 ){
            if(hasAncestor((*prunedGen)[i], 333, 321) && hasAncestor((*prunedGen)[i], 431, 211) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_Km++;
                phi_Km = (*prunedGen)[i].phi();
                eta_Km = (*prunedGen)[i].eta();
                Gen_pt_Km.push_back((*prunedGen)[i].pt());
                Gen_p_Km.push_back((*prunedGen)[i].p());
            }
        }
        else if( (*prunedGen)[i].pdgId()==211 ){
            if(hasAncestor((*prunedGen)[i], 431, 333) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_pi++;
                phi_pi = (*prunedGen)[i].phi();
                eta_pi = (*prunedGen)[i].eta();
                Gen_pt_pi.push_back((*prunedGen)[i].pt());
                Gen_p_pi.push_back((*prunedGen)[i].p());
            }
        }
    }

    Gen_dR_Kp_Km.push_back(getDeltaR(eta_Kp, phi_Kp, eta_Km, phi_Km));
    Gen_dR_Kp_pi.push_back(getDeltaR(eta_Kp, phi_Kp, eta_pi, phi_pi));
    Gen_dR_Km_pi.push_back(getDeltaR(eta_Km, phi_Km, eta_pi, phi_pi));

    if(ctr_Kp==1 && ctr_Km==1 && ctr_pi==1){

        std::vector<int> idx_Kp_vec, idx_Km_vec, idx_pi_vec;

        for(size_t i=0; i<packedPF->size(); ++i){

            all_pt.push_back((*packedPF)[i].pt());
            all_p.push_back((*packedPF)[i].p());

            float dr_Kp = getDeltaR(eta_Kp, phi_Kp, (*packedPF)[i].eta(), (*packedPF)[i].phi());
            float dr_Km = getDeltaR(eta_Km, phi_Km, (*packedPF)[i].eta(), (*packedPF)[i].phi());
            float dr_pi = getDeltaR(eta_pi, phi_pi, (*packedPF)[i].eta(), (*packedPF)[i].phi());

            all_dR_Kp_vec.push_back(dr_Kp);
            all_dR_Km_vec.push_back(dr_Km);
            all_dR_pi_vec.push_back(dr_pi);

            if((*packedPF)[i].pdgId()==211){
                dR_Kp_vec.push_back(dr_Kp);
                idx_Kp_vec.push_back(i);
                dR_pi_vec.push_back(dr_pi);
                idx_pi_vec.push_back(i);
            }
            else if((*packedPF)[i].pdgId()==-211){
                dR_Km_vec.push_back(dr_Km);
                idx_Km_vec.push_back(i);
            }
        }

        if(dR_Kp_vec.size()>0 && dR_Km_vec.size()>0 && dR_pi_vec.size()>0){

            auto min_ele_Kp = std::min_element(dR_Kp_vec.begin(), dR_Kp_vec.end());
            auto min_ele_Km = std::min_element(dR_Km_vec.begin(), dR_Km_vec.end());
            auto min_ele_pi = std::min_element(dR_pi_vec.begin(), dR_pi_vec.end());

            int idx_min_ele_Kp = std::distance(dR_Kp_vec.begin(), min_ele_Kp);
            int idx_min_ele_Km = std::distance(dR_Km_vec.begin(), min_ele_Km);
            int idx_min_ele_pi = std::distance(dR_pi_vec.begin(), min_ele_pi);

            int idx_Kp = idx_Kp_vec[idx_min_ele_Kp]; 
            int idx_Km = idx_Km_vec[idx_min_ele_Km]; 
            int idx_pi = idx_pi_vec[idx_min_ele_pi]; 

            if(idx_Kp == idx_pi){
                float min_dR = 100;
                if(*min_ele_Kp > *min_ele_pi){
                    for(int i=0; i<int(dR_Kp_vec.size()) && i!=idx_Kp; i++){
                        if(dR_Kp_vec[i] < min_dR){
                            idx_Kp = idx_Kp_vec[i];
                            min_dR = dR_Kp_vec[i];
                        }
                    }
                }
                else{
                    for(int i=0; i<int(dR_pi_vec.size()) && i!=idx_pi; i++){
                        if(dR_pi_vec[i] < min_dR){
                            idx_pi = idx_pi_vec[i];
                            min_dR = dR_pi_vec[i];
                        }
                    }
                }
            }

            if(all_dR_Kp_vec[idx_Kp]<0.03 && all_dR_Km_vec[idx_Km]<0.03 && (*packedPF)[idx_Kp].hasTrackDetails() && (*packedPF)[idx_Km].hasTrackDetails()){

                pt_Kp.push_back((*packedPF)[idx_Kp].pt());
                p_Kp.push_back((*packedPF)[idx_Kp].p());
                pt_Km.push_back((*packedPF)[idx_Km].pt());
                p_Km.push_back((*packedPF)[idx_Km].p());
                dR_Kp_Km.push_back(getDeltaR((*packedPF)[idx_Kp].eta(), (*packedPF)[idx_Kp].phi(), (*packedPF)[idx_Km].eta(), (*packedPF)[idx_Km].phi()));

                edm::ESHandle<TransientTrackBuilder> theB;
                iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
                KalmanVertexFitter fitter(true);

                math::XYZPoint vtx_Kp = (*packedPF)[idx_Kp].vertex();
                math::XYZPoint vtx_Km = (*packedPF)[idx_Km].vertex();
                d_Kp_Km.push_back((vtx_Kp-vtx_Km).R());

                std::vector<reco::TransientTrack> ttrk_phi;
                reco::Track trk_Kp = (*packedPF)[idx_Kp].pseudoTrack();
                reco::Track trk_Km = (*packedPF)[idx_Km].pseudoTrack();
                reco::TransientTrack ttrk_Kp = (*theB).build(trk_Kp);
                reco::TransientTrack ttrk_Km = (*theB).build(trk_Km);
                ttrk_phi.push_back(ttrk_Kp);
                ttrk_phi.push_back(ttrk_Km);
                TransientVertex tvx_phi;
                tvx_phi = fitter.vertex(ttrk_phi);

                if(tvx_phi.isValid()){
                    if(tvx_phi.hasRefittedTracks()){
                        std::vector<reco::TransientTrack> rf_ttrk_phi = tvx_phi.refittedTracks();
                        reco::Track rf_phi_trk_Kp = rf_ttrk_phi[0].track();
                        reco::Track rf_phi_trk_Km = rf_ttrk_phi[1].track();

                        TLorentzVector rf_phi_p4_Kp, rf_phi_p4_Km;
                        rf_phi_p4_Kp.SetXYZM(rf_phi_trk_Kp.px(), rf_phi_trk_Kp.py(), rf_phi_trk_Kp.pz(), mK);
                        rf_phi_p4_Km.SetXYZM(rf_phi_trk_Km.px(), rf_phi_trk_Km.py(), rf_phi_trk_Km.pz(), mK);

                        chi2_phi.push_back(tvx_phi.totalChiSquared());
                        m_phi.push_back((rf_phi_p4_Kp + rf_phi_p4_Km).M());

                        if(all_dR_pi_vec[idx_pi]<0.03 && (*packedPF)[idx_pi].hasTrackDetails()){

                            pt_pi.push_back((*packedPF)[idx_pi].pt());
                            p_pi.push_back((*packedPF)[idx_pi].p());
                            dR_Kp_pi.push_back(getDeltaR((*packedPF)[idx_Kp].eta(), (*packedPF)[idx_Kp].phi(), (*packedPF)[idx_pi].eta(), (*packedPF)[idx_pi].phi()));
                            dR_Km_pi.push_back(getDeltaR((*packedPF)[idx_Km].eta(), (*packedPF)[idx_Km].phi(), (*packedPF)[idx_pi].eta(), (*packedPF)[idx_pi].phi()));
                            
                            math::XYZPoint vtx_phi(tvx_phi.position().x(), tvx_phi.position().y(), tvx_phi.position().z());
                            math::XYZPoint vtx_pi = (*packedPF)[idx_pi].vertex();
                            d_phi_pi.push_back((vtx_phi-vtx_pi).R());
                            d_Kp_pi.push_back((vtx_Kp-vtx_pi).R());
                            d_Km_pi.push_back((vtx_Km-vtx_pi).R());

                            std::vector<reco::TransientTrack> ttrk_Ds;
                            reco::Track trk_pi = (*packedPF)[idx_pi].pseudoTrack();
                            reco::TransientTrack ttrk_pi = (*theB).build(trk_pi);
                            ttrk_Ds.push_back(rf_ttrk_phi[0]);
                            ttrk_Ds.push_back(rf_ttrk_phi[1]);
                            ttrk_Ds.push_back(ttrk_pi);
                            TransientVertex tvx_Ds;
                            tvx_Ds = fitter.vertex(ttrk_Ds);

                            if(tvx_Ds.isValid()){
                                if(tvx_Ds.hasRefittedTracks()){
                                    math::XYZPoint vtx_Ds(tvx_Ds.position().x(), tvx_Ds.position().y(), tvx_Ds.position().z());
                                    d_phi_Ds.push_back((vtx_phi-vtx_Ds).R());

                                    std::vector<reco::TransientTrack> rf_ttrk_Ds = tvx_Ds.refittedTracks();
                                    reco::Track rf_Ds_trk_Kp = rf_ttrk_Ds[0].track();
                                    reco::Track rf_Ds_trk_Km = rf_ttrk_Ds[1].track();
                                    reco::Track rf_Ds_trk_pi = rf_ttrk_Ds[2].track();

                                    TLorentzVector rf_Ds_p4_Kp, rf_Ds_p4_Km, rf_Ds_p4_pi;
                                    rf_Ds_p4_Kp.SetXYZM(rf_Ds_trk_Kp.px(), rf_Ds_trk_Kp.py(), rf_Ds_trk_Kp.pz(), mK);
                                    rf_Ds_p4_Km.SetXYZM(rf_Ds_trk_Km.px(), rf_Ds_trk_Km.py(), rf_Ds_trk_Km.pz(), mK);
                                    rf_Ds_p4_pi.SetXYZM(rf_Ds_trk_pi.px(), rf_Ds_trk_pi.py(), rf_Ds_trk_pi.pz(), mpi);

                                    chi2_Ds.push_back(tvx_Ds.totalChiSquared());
                                    m_Ds.push_back((rf_Ds_p4_Kp + rf_Ds_p4_Km + rf_Ds_p4_pi).M());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    mytree->Fill();
    return;
}

float GenPacked::getDeltaR(float eta1, float phi1, float eta2, float phi2)
{      
    float DeltaPhi = fabs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

bool GenPacked::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const{
    bool hasancestor = false;

    if(gp.numberOfMothers()==0) return false;

    if(gp.numberOfMothers()==1 && gp.mother(0)->pdgId()==motherid){
        if(gp.mother(0)->numberOfDaughters()==2){
            if(gp.mother(0)->daughter(0)->pdgId()==otherparticleid || gp.mother(0)->daughter(1)->pdgId()==otherparticleid) hasancestor = true;
        }
    }

    if (hasancestor) return true;
    else return hasAncestor(*gp.motherRef(0), motherid, otherparticleid);
}

void GenPacked::beginJob() {}

void GenPacked::endJob() {}

void GenPacked::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPacked);

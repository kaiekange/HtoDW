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
        std::vector<float> all_dR_Kp_vec, all_dR_Km_vec, all_dR_pi_vec, all_dR_mu_vec;
        std::vector<float> dR_Kp_vec, dR_Km_vec, dR_pi_vec, dR_mu_vec;
};

GenPacked::GenPacked(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    mytree = myfile->make<TTree>("Events", "Events");
    mytree->Branch("all_dR_Kp_vec", &all_dR_Kp_vec);
    mytree->Branch("all_dR_Km_vec", &all_dR_Km_vec);
    mytree->Branch("all_dR_pi_vec", &all_dR_pi_vec);
    mytree->Branch("all_dR_mu_vec", &all_dR_mu_vec);
    mytree->Branch("dR_Kp_vec", &dR_Kp_vec);
    mytree->Branch("dR_Km_vec", &dR_Km_vec);
    mytree->Branch("dR_pi_vec", &dR_pi_vec);
    mytree->Branch("dR_mu_vec", &dR_mu_vec);
}

GenPacked::~GenPacked() {}

void GenPacked::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    all_dR_Kp_vec.clear();
    all_dR_Km_vec.clear();
    all_dR_pi_vec.clear();
    all_dR_mu_vec.clear();
    dR_Kp_vec.clear();
    dR_Km_vec.clear();
    dR_pi_vec.clear();
    dR_mu_vec.clear();

    edm::Handle<reco::GenParticleCollection> prunedGen;
    iEvent.getByToken(prunedGenToken, prunedGen);

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    float Kp_phi=100, Kp_eta=100, Km_phi=100, Km_eta=100, pi_phi=100, pi_eta=100, mu_phi=100, mu_eta=100;
    int ctr_Kp=0, ctr_Km=0, ctr_pi=0, ctr_mu=0; 

    for(size_t i=0; i<prunedGen->size(); ++i){
        if( (*prunedGen)[i].pdgId()==321 ){
            if(hasAncestor((*prunedGen)[i], 333, -321) && hasAncestor((*prunedGen)[i], 431, 211) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_Kp++;
                Kp_phi = (*prunedGen)[i].phi();
                Kp_eta = (*prunedGen)[i].eta();
            }
        }
        else if( (*prunedGen)[i].pdgId()==-321 ){
            if(hasAncestor((*prunedGen)[i], 333, 321) && hasAncestor((*prunedGen)[i], 431, 211) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_Km++;
                Km_phi = (*prunedGen)[i].phi();
                Km_eta = (*prunedGen)[i].eta();
            }
        }
        else if( (*prunedGen)[i].pdgId()==211 ){
            if(hasAncestor((*prunedGen)[i], 431, 333) && hasAncestor((*prunedGen)[i], 25, -24)){
                ctr_pi++;
                pi_phi = (*prunedGen)[i].phi();
                pi_eta = (*prunedGen)[i].eta();
            }
        }
        else if( (*prunedGen)[i].pdgId()==13 && (*prunedGen)[i].isHardProcess() ){
            if(hasAncestor((*prunedGen)[i], -24, -14) && hasAncestor((*prunedGen)[i], 25, 431)){
                ctr_mu++;
                mu_phi = (*prunedGen)[i].phi();
                mu_eta = (*prunedGen)[i].eta();
            }
        }
    }

    if((ctr_Kp==1) && (ctr_Km==1) && (ctr_pi==1) && (ctr_mu==1)){

        std::vector<int> idx_Kp_vec, idx_Km_vec, idx_pi_vec, idx_mu_vec;

        for(size_t i=0; i<packedPF->size(); ++i){

            float dr_Kp = getDeltaR(Kp_eta, Kp_phi, (*packedPF)[i].eta(), (*packedPF)[i].phi());
            float dr_Km = getDeltaR(Km_eta, Km_phi, (*packedPF)[i].eta(), (*packedPF)[i].phi());
            float dr_pi = getDeltaR(pi_eta, pi_phi, (*packedPF)[i].eta(), (*packedPF)[i].phi());
            float dr_mu = getDeltaR(mu_eta, mu_phi, (*packedPF)[i].eta(), (*packedPF)[i].phi());

            all_dR_Kp_vec.push_back(dr_Kp);
            all_dR_Km_vec.push_back(dr_Km);
            all_dR_pi_vec.push_back(dr_pi);
            all_dR_mu_vec.push_back(dr_mu);

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
            else if((*packedPF)[i].pdgId()==13){
                dR_mu_vec.push_back(dr_mu);
                idx_mu_vec.push_back(i);
            }

        }


        if(dR_Kp_vec.size()>0 && dR_Km_vec.size()>0 && dR_pi_vec.size()>0 && dR_mu_vec.size()>0){

            auto min_ele_Kp = std::min_element(dR_Kp_vec.begin(), dR_Kp_vec.end());
            auto min_ele_Km = std::min_element(dR_Km_vec.begin(), dR_Km_vec.end());
            auto min_ele_pi = std::min_element(dR_pi_vec.begin(), dR_pi_vec.end());
            auto min_ele_mu = std::min_element(dR_mu_vec.begin(), dR_mu_vec.end());

            int idx_min_ele_Kp = std::distance(dR_Kp_vec.begin(), min_ele_Kp);
            int idx_min_ele_Km = std::distance(dR_Km_vec.begin(), min_ele_Km);
            int idx_min_ele_pi = std::distance(dR_pi_vec.begin(), min_ele_pi);
            int idx_min_ele_mu = std::distance(dR_mu_vec.begin(), min_ele_mu);

            int idx_Kp = idx_Kp_vec[idx_min_ele_Kp]; 
            int idx_Km = idx_Km_vec[idx_min_ele_Km]; 
            int idx_pi = idx_pi_vec[idx_min_ele_pi]; 
            int idx_mu = idx_mu_vec[idx_min_ele_mu]; 
            
            /* std::cout << all_dR_Kp_vec[idx_Kp] << " " << all_dR_Km_vec[idx_Km] << " " << all_dR_pi_vec[idx_pi] << " " << all_dR_mu_vec[idx_mu] << std::endl; */ 
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
                /* std::cout << "check here" << std::endl; */
                /* std::cout << all_dR_Kp_vec[idx_Kp] << " " << all_dR_Km_vec[idx_Km] << " " << all_dR_pi_vec[idx_pi] << " " << all_dR_mu_vec[idx_mu] << std::endl; */ 
            }


        }
        mytree->Fill();
    }


    return;
}

float GenPacked::getDeltaR(float eta1, float phi1, float eta2, float phi2)
{      
    float DeltaPhi = fabs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

/* bool GenPacked::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const{ */
/*     bool hasancestor = false; */

/*     if(gp.numberOfMothers()!=1) return false; */
/*     if(gp.mother(0)->numberOfDaughters()!=2) return false; */
/*     if( gp.mother(0)->pdgId()==motherid && (gp.mother(0)->daughter(0)->pdgId()==otherparticleid || gp.mother(0)->daughter(1)->pdgId()==otherparticleid) ) hasancestor = true; */

/*     if (hasancestor) return true; */
/*     else return hasAncestor(*gp.motherRef(0), motherid, otherparticleid); */
/* } */

bool GenPacked::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const{
    bool hasancestor = false;

    if(gp.numberOfMothers()==0) return false;

    if(gp.numberOfMothers()==1 && gp.mother(0)->pdgId()==motherid){
        if(gp.mother(0)->numberOfDaughters()==2){
            if(gp.mother(0)->daughter(0)->pdgId()==otherparticleid || gp.mother(0)->daughter(1)->pdgId()==otherparticleid) hasancestor = true;
        }
    }

    /* if (hasancestor && motherid == 25 && otherparticleid == 431){ */
    /*     std::cout << gp.mother(0)->daughter(0)->numberOfDaughters() << std::endl; */
    /*     for(size_t i=0; i<gp.mother(0)->daughter(0)->numberOfDaughters(); ++i) std::cout << gp.mother(0)->daughter(0)->daughter(i)->pdgId() << std::endl; */
    /*     std::cout << gp.mother(0)->daughter(1)->numberOfDaughters() << std::endl; */
    /*     for(size_t i=0; i<gp.mother(0)->daughter(1)->numberOfDaughters(); ++i) std::cout << gp.mother(0)->daughter(1)->daughter(i)->pdgId() << std::endl; */
    /* } */
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

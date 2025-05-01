// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenInfo
//
/**\class GenInfo GenInfo.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenInfo.cc

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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TFile.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class GenInfo : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenInfo(const edm::ParameterSet&);
        ~GenInfo();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken;

        edm::Service<TFileService> myfile;
        TTree *mytree{};

        std::vector<double> H_x;
        std::vector<double> H_y;
        std::vector<double> H_z;
        std::vector<double> H_chi2;
        std::vector<double> Ds_x;
        std::vector<double> Ds_y;
        std::vector<double> Ds_z;
        std::vector<double> Ds_chi2;
        std::vector<double> W_x;
        std::vector<double> W_y;
        std::vector<double> W_z;
        std::vector<double> W_chi2;
        std::vector<double> phi_x;
        std::vector<double> phi_y;
        std::vector<double> phi_z;
        std::vector<double> phi_chi2;
        std::vector<double> Kp_x;
        std::vector<double> Kp_y;
        std::vector<double> Kp_z;
        std::vector<double> Kp_chi2;
        std::vector<double> Km_x;
        std::vector<double> Km_y;
        std::vector<double> Km_z;
        std::vector<double> Km_chi2;
        std::vector<double> pi_x;
        std::vector<double> pi_y;
        std::vector<double> pi_z;
        std::vector<double> pi_chi2;
        std::vector<double> mu_x;
        std::vector<double> mu_y;
        std::vector<double> mu_z;
        std::vector<double> mu_chi2;
        std::vector<double> nu_x;
        std::vector<double> nu_y;
        std::vector<double> nu_z;
        std::vector<double> nu_chi2;
};

GenInfo::GenInfo(const edm::ParameterSet& iConfig) :
    GenPartToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")))
{
    mytree = myfile->make<TTree>("Events", "Events");
    mytree->Branch("H_x", &H_x);
    mytree->Branch("H_y", &H_y);
    mytree->Branch("H_z", &H_z);
    mytree->Branch("H_chi2", &H_chi2);
    mytree->Branch("Ds_x", &Ds_x);
    mytree->Branch("Ds_y", &Ds_y);
    mytree->Branch("Ds_z", &Ds_z);
    mytree->Branch("Ds_chi2", &Ds_chi2);
    mytree->Branch("W_x", &W_x);
    mytree->Branch("W_y", &W_y);
    mytree->Branch("W_z", &W_z);
    mytree->Branch("W_chi2", &W_chi2);
    mytree->Branch("phi_x", &phi_x);
    mytree->Branch("phi_y", &phi_y);
    mytree->Branch("phi_z", &phi_z);
    mytree->Branch("phi_chi2", &phi_chi2);
    mytree->Branch("Kp_x", &Kp_x);
    mytree->Branch("Kp_y", &Kp_y);
    mytree->Branch("Kp_z", &Kp_z);
    mytree->Branch("Kp_chi2", &Kp_chi2);
    mytree->Branch("Km_x", &Km_x);
    mytree->Branch("Km_y", &Km_y);
    mytree->Branch("Km_z", &Km_z);
    mytree->Branch("Km_chi2", &Km_chi2);
    mytree->Branch("pi_x", &pi_x);
    mytree->Branch("pi_y", &pi_y);
    mytree->Branch("pi_z", &pi_z);
    mytree->Branch("pi_chi2", &pi_chi2);
    mytree->Branch("mu_x", &mu_x);
    mytree->Branch("mu_y", &mu_y);
    mytree->Branch("mu_z", &mu_z);
    mytree->Branch("mu_chi2", &mu_chi2);
    mytree->Branch("nu_x", &nu_x);
    mytree->Branch("nu_y", &nu_y);
    mytree->Branch("nu_z", &nu_z);
    mytree->Branch("nu_chi2", &nu_chi2);
}

GenInfo::~GenInfo() {}

void GenInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(GenPartToken, pruned);

    H_x.clear();    
    H_y.clear();    
    H_z.clear();    
    H_chi2.clear();    
    Ds_x.clear();
    Ds_y.clear();
    Ds_z.clear();
    Ds_chi2.clear();
    W_x.clear();
    W_y.clear();
    W_z.clear();
    W_chi2.clear();
    phi_x.clear();
    phi_y.clear();
    phi_z.clear();
    phi_chi2.clear();
    Kp_x.clear();
    Kp_y.clear();
    Kp_z.clear();
    Kp_chi2.clear();
    Km_x.clear();
    Km_y.clear();
    Km_z.clear();
    Km_chi2.clear();
    pi_x.clear();
    pi_y.clear();
    pi_z.clear();
    pi_chi2.clear();
    mu_x.clear();
    mu_y.clear();
    mu_z.clear();
    mu_chi2.clear();
    nu_x.clear();
    nu_y.clear();
    nu_z.clear();
    nu_chi2.clear();

    reco::GenParticleCollection::const_iterator it;

    for(it = pruned->begin(); it != pruned->end(); it++){

        if(it->pdgId()==25){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==431 && pdg1==-24) || (pdg0==-24 && pdg1==431) );
                if(correctDaughter) {
                    H_x.push_back(it->vx());
                    H_y.push_back(it->vy());
                    H_z.push_back(it->vz());
                    H_chi2.push_back(it->vertexChi2());
                }
            }
        }

        else if(it->pdgId()==431){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==333 && pdg1==211) || (pdg0==211 && pdg1==333) );
                bool correctAncestor = hasAncestor(*it, 25, -24);
                if(correctAncestor && correctDaughter) {
                    Ds_x.push_back(it->vx());
                    Ds_y.push_back(it->vy());
                    Ds_z.push_back(it->vz());
                    Ds_chi2.push_back(it->vertexChi2());
                }
            }
        }

        else if(it->pdgId()==-24){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==13 && pdg1==-14) || (pdg0==-14 && pdg1==13) );
                bool correctAncestor = hasAncestor(*it, 25, 431);
                if(correctAncestor && correctDaughter){
                    W_x.push_back(it->vx());
                    W_y.push_back(it->vy());
                    W_z.push_back(it->vz());
                    W_chi2.push_back(it->vertexChi2());
                }
            }
        }

        else if(it->pdgId()==333){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==321 && pdg1==-321) || (pdg0==-321 && pdg1==321) );
                bool HAncestor = hasAncestor(*it, 25, -24);
                bool DsAncestor = hasAncestor(*it, 431, 211);
                if(HAncestor && DsAncestor && correctDaughter){
                    phi_x.push_back(it->vx());
                    phi_y.push_back(it->vy());
                    phi_z.push_back(it->vz());
                    phi_chi2.push_back(it->vertexChi2());
                }
            }
        }

        else if(it->pdgId()==211){
            bool HAncestor = hasAncestor(*it, 25, -24);
            bool DsAncestor = false;
            if( it->mother(0)->pdgId()==431 ){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==333) DsAncestor = true;
                }
            }
            if(HAncestor && DsAncestor){
                pi_x.push_back(it->vx());
                pi_y.push_back(it->vy());
                pi_z.push_back(it->vz());
                pi_chi2.push_back(it->vertexChi2());
            } 
        }

        else if(it->pdgId()==321){
            bool HAncestor = hasAncestor(*it, 25, -24);
            bool DsAncestor = hasAncestor(*it, 431, 211);
            bool phiAncestor = false;
            if(it->mother(0)->pdgId()==333){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==-321) phiAncestor = true;
                }
            }
            if(HAncestor && DsAncestor && phiAncestor){
                Kp_x.push_back(it->vx());
                Kp_y.push_back(it->vy());
                Kp_z.push_back(it->vz());
                Kp_chi2.push_back(it->vertexChi2());
            }
        }

        else if(it->pdgId()==-321){
            bool HAncestor = hasAncestor(*it, 25, -24);
            bool DsAncestor = hasAncestor(*it, 431, 211);
            bool phiAncestor = false;
            if(it->mother(0)->pdgId()==333){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==321) phiAncestor = true;
                }
            }
            if(HAncestor && DsAncestor && phiAncestor){
                Km_x.push_back(it->vx());
                Km_y.push_back(it->vy());
                Km_z.push_back(it->vz());
                Km_chi2.push_back(it->vertexChi2());
            }
        }

        else if(it->pdgId()==13){
            bool HAncestor = hasAncestor(*it, 25, 431);
            bool WAncestor = false;
            if(it->mother(0)->pdgId()==-24){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==-14) WAncestor = true;
                }
            }
            if(HAncestor && WAncestor){
                mu_x.push_back(it->vx());
                mu_y.push_back(it->vy());
                mu_z.push_back(it->vz());
                mu_chi2.push_back(it->vertexChi2());
            } 
        }

        else if(it->pdgId()==-14){
            bool HAncestor = hasAncestor(*it, 25, 431);
            bool WAncestor = false;
            if(it->mother(0)->pdgId()==-24){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==13) WAncestor = true;
                }
            }
            if(HAncestor && WAncestor){
                nu_x.push_back(it->vx());
                nu_y.push_back(it->vy());
                nu_z.push_back(it->vz());
                nu_chi2.push_back(it->vertexChi2());
            }
        }
    }

    mytree->Fill();
    return;
}

bool GenInfo::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const{
    bool hasancestor = false;

    if(gp.numberOfMothers()==0) return false;

    if(gp.numberOfMothers()==1){
        if(gp.mother(0)->pdgId()==motherid){
            for(size_t i=0; i<gp.mother(0)->numberOfDaughters(); i++){
                if(gp.mother(0)->daughter(i)->pdgId()==otherparticleid) hasancestor = true;
            }
        }
    }

    if (hasancestor) return true;
    else return hasAncestor(*gp.motherRef(0), motherid, otherparticleid);
}

void GenInfo::beginJob() {}

void GenInfo::endJob() {}

void GenInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenInfo);

// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenParticleMatch
//
/**\class GenParticleMatch GenParticleMatch.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenParticleMatch.cc

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


class GenParticleMatch : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenParticleMatch(const edm::ParameterSet&);
        ~GenParticleMatch();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken;

        edm::Service<TFileService> myfile;
        TTree *mytree{};

        std::vector<reco::GenParticle> match_part;
};

GenParticleMatch::GenParticleMatch(const edm::ParameterSet& iConfig) :
    GenPartToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")))
{
    mytree = myfile->make<TTree>("Events", "Events");
    mytree->Branch("match_part", &match_part);
}

GenParticleMatch::~GenParticleMatch() {}

void GenParticleMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(GenPartToken, pruned);

    match_part.clear();    

    reco::GenParticleCollection::const_iterator it;

    for(it = pruned->begin(); it != pruned->end(); it++){

        if(it->pdgId()==25){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==431 && pdg1==-24) || (pdg0==-24 && pdg1==431) );
                if(correctDaughter) match_part.push_back(*it);
            }
        }

        else if(it->pdgId()==431){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==333 && pdg1==211) || (pdg0==211 && pdg1==333) );
                bool correctAncestor = hasAncestor(*it, 25, -24);
                if(correctAncestor && correctDaughter) match_part.push_back(*it);
            }
        }

        else if(it->pdgId()==-24){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==13 && pdg1==-14) || (pdg0==-14 && pdg1==13) );
                bool correctAncestor = hasAncestor(*it, 25, 431);
                if(correctAncestor && correctDaughter) match_part.push_back(*it);
            }
        }

        else if(it->pdgId()==333){
            if( it->numberOfDaughters()==2 ){
                int pdg0 = it->daughter(0)->pdgId();
                int pdg1 = it->daughter(1)->pdgId();
                bool correctDaughter = ( (pdg0==321 && pdg1==-321) || (pdg0==-321 && pdg1==321) );
                bool HAncestor = hasAncestor(*it, 25, -24);
                bool DsAncestor = hasAncestor(*it, 431, 211);
                if(HAncestor && DsAncestor && correctDaughter) match_part.push_back(*it);
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
            if(HAncestor && DsAncestor) match_part.push_back(*it);
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
            if(HAncestor && DsAncestor && phiAncestor) match_part.push_back(*it);
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
            if(HAncestor && DsAncestor && phiAncestor) match_part.push_back(*it);
        }

        else if(it->pdgId()==13){
            bool HAncestor = hasAncestor(*it, 25, 431);
            bool WAncestor = false;
            if(it->mother(0)->pdgId()==-24){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==-14) WAncestor = true;
                }
            }
            if(HAncestor && WAncestor) match_part.push_back(*it);
        }

        else if(it->pdgId()==-14){
            bool HAncestor = hasAncestor(*it, 25, 431);
            bool WAncestor = false;
            if(it->mother(0)->pdgId()==-24){
                for(size_t i=0; i<it->mother(0)->numberOfDaughters(); i++){
                    if(it->mother(0)->daughter(i)->pdgId()==13) WAncestor = true;
                }
            }
            if(HAncestor && WAncestor) match_part.push_back(*it);
        }
    }

    mytree->Fill();
    return;
}

bool GenParticleMatch::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const{
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

void GenParticleMatch::beginJob() {}

void GenParticleMatch::endJob() {}

void GenParticleMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleMatch);

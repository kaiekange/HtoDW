// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenParticleAnalyzer
//
/**\class GenParticleAnalyzer GenParticleAnalyzer.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenParticleAnalyzer.cc

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


class GenParticleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenParticleAnalyzer(const edm::ParameterSet&);
        ~GenParticleAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken;

        TTree *mytree;

        int mynumCandidates;
        std::vector<int> mystatus;
        std::vector<int> mypdgId;
        std::vector<int> mynummother;
        std::vector<int> mynumdaughter;
        std::vector<std::vector<int>> mymotherId;
        std::vector<std::vector<int>> mydaughterId;
};

GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig) :
    GenPartToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")))
{
    edm::Service<TFileService> myfile;
    mytree = myfile->make<TTree>("Events", "Events");

    mytree->Branch("numCandidates", &mynumCandidates);
    mytree->Branch("status", &mystatus);
    mytree->Branch("pdgId", &mypdgId);
    mytree->Branch("nummother", &mynummother);
    mytree->Branch("numdaughter", &mynumdaughter);
    mytree->Branch("motherId", &mymotherId);
    mytree->Branch("daughterId", &mydaughterId);

}

GenParticleAnalyzer::~GenParticleAnalyzer() {}


void GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    mynumCandidates = 0;
    mystatus.clear();
    mypdgId.clear();
    mynummother.clear();
    mynumdaughter.clear();
    mymotherId.clear();
    mydaughterId.clear();

    Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(GenPartToken, pruned);

    reco::GenParticleCollection::const_iterator it;
   
    for (it = pruned->begin(); it != pruned->end(); it++){
        mynumCandidates++;
        mystatus.push_back(it->status());
        mypdgId.push_back(it->pdgId());
        mynummother.push_back(it->numberOfMothers());
        mynumdaughter.push_back(it->numberOfDaughters());
        std::vector<int> motherId;
        for(size_t i=0; i<it->numberOfMothers(); i++){
            motherId.push_back(it->mother(i)->pdgId());
        }
        mymotherId.push_back(motherId);
        std::vector<int> daughterId;
        for(size_t i=0; i<it->numberOfDaughters(); i++){
            daughterId.push_back(it->daughter(i)->pdgId());
        }
        mydaughterId.push_back(daughterId);   
    }

    mytree->Fill();
    return;
}

void GenParticleAnalyzer::beginJob() {}

void GenParticleAnalyzer::endJob() {}

void GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);

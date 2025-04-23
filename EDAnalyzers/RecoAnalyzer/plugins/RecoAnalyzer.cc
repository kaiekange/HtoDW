// -*- C++ -*-
//
// Package:    EDAnalyzers/RecoAnalyzer
// Class:      RecoAnalyzer
//
/**\class RecoAnalyzer RecoAnalyzer.cc EDAnalyzers/RecoAnalyzer/plugins/RecoAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  cmsusers user
//         Created:  Tue, 22 Apr 2025 14:40:41 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/iterator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class RecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit RecoAnalyzer(const edm::ParameterSet&);
        ~RecoAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
        edm::Service<TFileService> outfile;
        TTree *mytree;

        std::vector<float> mass;
        /* std::vector<float> m_phi; */
        /* std::vector<float> dv; */
};

RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig)
    :
        packedToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packed")))

{
    mytree = outfile->make<TTree>("Events", "Events");
    mytree->Branch("mass", &mass);
    /* mytree->Branch("m_phi", &m_phi); */ 
    /* mytree->Branch("dv", &dv); */ 
}

RecoAnalyzer::~RecoAnalyzer(){}

void RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    mass.clear();
    int ctr = 0;
    /* m_phi.clear(); */
    /* dv.clear(); */

    edm::Handle<pat::PackedCandidateCollection> packed;
    iEvent.getByToken(packedToken_, packed);

    std::vector<int> mycharge;
    std::vector<TLorentzVector> myp4;
    std::vector<ROOT::Math::XYZVector> myvertex;

    for(const pat::PackedCandidate &pfcandidate : *packed){
        if(pfcandidate.trackHighPurity() && (abs(pfcandidate.charge()) == 1) && (pfcandidate.pt()>1)){
            mass.push_back(pfcandidate.mass());
            if (pfcandidate.mass()>0.49 && pfcandidate.mass()<0.5) {
                ctr++;
                std::cout << pfcandidate.mass() << std::endl;
            }
            /* mycharge.push_back(pfcandidate.charge()); */
            /* myp4.push_back(TLorentzVector(pfcandidate.px(), pfcandidate.py(), pfcandidate.pz(), pfcandidate.energy())); */
            /* ROOT::Math::XYZVector pfvertex(pfcandidate.vertex().x(), pfcandidate.vertex().y(), pfcandidate.vertex().z()); */ 
            /* myvertex.push_back(pfvertex); */
        }
    }

    std::cout << "number of kaon " << ctr << std::endl; 

/*     int partsize = mycharge.size(); */

/*     std::cout << partsize << std::endl; */

/*     if (partsize>3) { */
/*         for(int i=0; i<partsize-1; i++){ */
/*             for(int j=i+1; j<partsize; j++){ */
/*                 if(mycharge[i] == mycharge[j]) continue; */
/*                 if((myvertex[i] - myvertex[j]).R() > 0.1) continue; */
/*                 m_phi.push_back((myp4[i] + myp4[j]).M()); */
/*             } */
/*         } */
/*     } */

    mytree->Fill();
}


void RecoAnalyzer::beginJob(){}

void RecoAnalyzer::endJob(){}

void RecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoAnalyzer);

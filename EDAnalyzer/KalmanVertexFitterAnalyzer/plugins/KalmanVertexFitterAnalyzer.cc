// -*- C++ -*-
//
// Package:    EDAnalyzer/KalmanVertexFitterAnalyzer
// Class:      KalmanVertexFitterAnalyzer
//
/**\class KalmanVertexFitterAnalyzer KalmanVertexFitterAnalyzer.cc EDAnalyzer/KalmanVertexFitterAnalyzer/plugins/KalmanVertexFitterAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  cmsusers user
//         Created:  Fri, 20 Dec 2024 13:26:47 GMT
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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/iterator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class KalmanVertexFitterAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit KalmanVertexFitterAnalyzer(const edm::ParameterSet&);
        ~KalmanVertexFitterAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<pat::PackedCandidateCollection> token_packed;
        edm::ParameterSet kvfPSet;
};

KalmanVertexFitterAnalyzer::KalmanVertexFitterAnalyzer(const edm::ParameterSet& iConfig)
    : token_packed(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packed"))) 
{
    kvfPSet = iConfig.getParameter<edm::ParameterSet>("KVFParameters");
}

KalmanVertexFitterAnalyzer::~KalmanVertexFitterAnalyzer() {}

void KalmanVertexFitterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::PackedCandidateCollection> packed;
    iEvent.getByToken(token_packed, packed);

    if(packed.isValid()){
        std::vector<reco::TransientTrack> tks;

        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

        for (const pat::PackedCandidate &pfcandidate : *packed){
            if(pfcandidate.hasTrackDetails()){
                reco::Track track = pfcandidate.pseudoTrack();
                reco::TransientTrack trajectory = (*theB).build(track);
                tks.push_back(trajectory);
            }
        }

        KalmanVertexFitter fitter(true);
        TransientVertex myTransientVertex = fitter.vertex(tks);

        if(myTransientVertex.isValid()){
            std::cout << myTransientVertex.position().x() << std::endl; 
        }
    }

    else {
        std::cout << "not valid" << std::endl;
    }

  }


void KalmanVertexFitterAnalyzer::beginJob(){}

void KalmanVertexFitterAnalyzer::endJob(){}

void KalmanVertexFitterAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(KalmanVertexFitterAnalyzer);

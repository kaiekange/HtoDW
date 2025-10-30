// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      MuonGen
//
/**\class MuonGen MuonGen.cc EDAnalyzer/GenParticleAnalyzer/plugins/MuonGen.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  cmsusers user
//         Created:  Fri, 15 Nov 2024 15:09:59 GMT
//
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "Compression.h"

class MuonGen : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit MuonGen(const edm::ParameterSet& iConfig);
        ~MuonGen();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedToken;
        edm::EDGetTokenT<pat::MuonCollection> muonToken;

        /* const edm::Service<TFileService> fs; */
        /* MuonTree *ftree; */

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        int num_hasGen;

        std::vector<float> dR_mu;

        std::vector<int> Gen_mu_idx;
        std::vector<int> Gen_nu_idx;
        std::vector<int> Gen_W_idx;

        std::vector<int> idx_mu_vec;
};

MuonGen::MuonGen(const edm::ParameterSet& iConfig) :
    prunedToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPart"))),
    muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
{
}

MuonGen::~MuonGen() {}

void MuonGen::beginJob()
{
    num_hasGen = 0;
    /* TFile & f = fs->file(); */
    /* f.SetCompressionAlgorithm(ROOT::kZLIB); */
    /* ftree = new MuonTree(fs->make<TTree>("Events", "Events")); */
    /* ftree->CreateBranches(); */
}

void MuonGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    int num_Gen_mu = 0;
    int num_Gen_nu = 0;
    int num_Gen_W = 0;

    dR_mu.clear();

    Gen_mu_idx.clear();
    Gen_nu_idx.clear();
    Gen_W_idx.clear();

    edm::Handle<reco::GenParticleCollection> prunedHandle;
    iEvent.getByToken(prunedToken, prunedHandle);

    edm::Handle<pat::MuonCollection> muonHandle;
    iEvent.getByToken(muonToken, muonHandle);

    edm::ESHandle<MagneticField> magField;
    iSetup.get<IdealMagneticFieldRecord>().get(magField);

    edm::ESHandle<TransientTrackBuilder> ttBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttBuilder);

    for(size_t i=0; i<prunedHandle->size(); i++){
        const auto& gp = (*prunedHandle)[i];
        const int pdg = gp.pdgId();
        if( pdg==13 && hasAncestor(gp, -24, -14) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_mu_idx.push_back(i);
            num_Gen_mu++;
        } else if( pdg==-14 && hasAncestor(gp, -24, 13) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_nu_idx.push_back(i);
            num_Gen_nu++;
        } else if( pdg==-24 && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_W_idx.push_back(i);
            num_Gen_W++;
        }
    }

    if(num_Gen_mu!=1) return;
    if(num_Gen_nu!=1) return;
    if(num_Gen_W!=1) return;

    const auto& mu_GP = (*prunedHandle)[Gen_mu_idx[0]];
/*     const auto& nu_GP = (*prunedHandle)[Gen_nu_idx[0]]; */
/*     const auto& W_GP = (*prunedHandle)[Gen_W_idx[0]]; */

    for(size_t i=0; i<muonHandle->size(); i++){
        const auto& muon = (*muonHandle)[i];

        if (!muon.genParticle()) continue;
        const reco::GenParticle& genMuon = *(muon.genParticle());

        /* if(genMuon.isHardProcess()){ */
        /*     if(hasAncestor(genMuon, -24, -14)) std::cout << "mother is W" << std::endl; */
        /*     num_hasGen++; */
        /* } */
        if(hasAncestor(genMuon, -24, -14) && hasAncestor(genMuon, 25, 431)){
            std::cout << reco::deltaR(muon.eta(), muon.phi(), mu_GP.eta(), mu_GP.phi()) << std::endl; 
            num_hasGen++;
        }
    }

    return;
}

bool MuonGen::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
{
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

void MuonGen::endJob()
{
    std::cout << "\n\n\n num_has_Gen " << num_hasGen << "\n\n\n" << std::endl;
    edm::LogInfo("MuonGen") << "Processed all events.";
}

void MuonGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonGen);

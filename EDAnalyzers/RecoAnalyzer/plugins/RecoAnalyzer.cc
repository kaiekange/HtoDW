// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      RecoAnalyzer
//
/**\class RecoAnalyzer RecoAnalyzer.cc EDAnalyzer/GenParticleAnalyzer/plugins/RecoAnalyzer.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "Compression.h"

#include "EDAnalyzers/RecoAnalyzer/interface/RecoTree.h"
#include "EDAnalyzers/RecoAnalyzer/interface/VertexReProducer.h"

class RecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit RecoAnalyzer(const edm::ParameterSet& iConfig);
        ~RecoAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken;
        edm::EDGetTokenT<pat::MuonCollection> muonToken;
        edm::EDGetTokenT<reco::VertexCollection> pvToken;
        edm::EDGetTokenT<reco::BeamSpot> bsToken;

        KalmanVertexFitter KVFitter;
        VertexReProducer *revertex;

        const edm::Service<TFileService> fs;
        RecoTree *ftree;

        std::vector<int> idx_Kp_vec;
        std::vector<int> idx_Km_vec;
        std::vector<int> idx_pi_vec;
        std::vector<int> idx_mu_vec;

        const double Mass_Ds = 1.96835;
        const double Mass_phi = 1.019460;
        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_mu = 0.1056583755;

        unsigned long nEventsProcessed_;
};

RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig) :
    packedToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    pvToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primvtx"))),
    bsToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    KVFitter(true)
{
    revertex = new VertexReProducer(iConfig);
}

RecoAnalyzer::~RecoAnalyzer() {}

void RecoAnalyzer::beginJob()
{
    TFile & f = fs->file();
    f.SetCompressionAlgorithm(ROOT::kZLIB);
    ftree = new RecoTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
    nEventsProcessed_ = 0;
}

void RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    ++nEventsProcessed_;
    ftree->Init();

    edm::Handle<pat::PackedCandidateCollection> packedHandle;
    iEvent.getByToken(packedToken, packedHandle);
    edm::Handle<pat::MuonCollection> muonHandle;
    iEvent.getByToken(muonToken, muonHandle);
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(pvToken, pvHandle);
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken, bsHandle);

    edm::ESHandle<MagneticField> magField;
    iSetup.get<IdealMagneticFieldRecord>().get(magField);
    edm::ESHandle<TransientTrackBuilder> ttBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttBuilder);

    ftree->beamspotInfo.fillAll(*bsHandle );

    ftree->primvtxInfo.clearAll();
    const reco::Vertex& primvtx = pvHandle->front();
    if(primvtx.isValid()) ftree->primvtxInfo.fillAll( primvtx );
    GlobalPoint primvtx_pos( primvtx.x(), primvtx.y(), primvtx.z() );
    GlobalError primvtx_cov( primvtx.covariance(0,0), primvtx.covariance(0,1), primvtx.covariance(0,2), primvtx.covariance(1,1), primvtx.covariance(1,2), primvtx.covariance(2,2) );

    reco::Track::CovarianceMatrix dummy_track_cov;
    for (unsigned int ic = 0; ic < 5; ++ic) {
        for (unsigned int jc = 0; jc < 5; ++jc) {
            if(ic == jc) dummy_track_cov(ic,jc) = 1e-4;
            else dummy_track_cov(ic,jc) = 0.0;
        }
    }

    idx_Kp_vec.clear();
    idx_Km_vec.clear();
    idx_pi_vec.clear();

    for(size_t i=0; i<packedHandle->size(); i++){
        const auto& pf = (*packedHandle)[i];

        if( !(pf.trackHighPurity()) ) continue;
        if( !(pf.hasTrackDetails()) ) continue;
        if( pf.pt() < 0.5 ) continue;
        if( pf.p() < 1 ) continue;
        if( abs(pf.pdgId()) != 211 ) continue;

        if(pf.charge() == 1){
            idx_Kp_vec.push_back(i);
            idx_pi_vec.push_back(i);
        }
        else if(pf.charge() == -1){
            idx_Km_vec.push_back(i);
        }

        /* if(nEventsProcessed_ == 1){ */
        /*     std::cout << pf.pseudoTrack().chi2() << ", " << pf.pseudoTrack().ndof() << std::endl; */
        /*     const reco::Track::CovarianceMatrix& cov = pf.pseudoTrack().covariance(); */

        /*     std::cout << "---------------------------\n"; */
        /*     for (int i = 0; i < cov.kRows; ++i) { */
        /*         for (int j = 0; j < cov.kCols; ++j) { */
        /*             std::cout << cov(i,j) << " "; */
        /*         } */
        /*         std::cout << "\n"; */
        /*     } */
        /*     std::cout << "---------------------------\n"; */
        /* } */
    }

    for(size_t i=0; i<idx_Kp_vec.size(); i++){
        ftree->Kp.clearAll();

        const auto& Kp_PF = (*packedHandle)[idx_Kp_vec[i]];
        ftree->Kp.fillAll( Kp_PF );
        TLorentzVector Kp_P4;
        Kp_P4.SetXYZM(Kp_PF.px(), Kp_PF.py(), Kp_PF.pz(), Mass_K);

        for(size_t j=0; j<idx_Km_vec.size(); j++){
            ftree->Km.clearAll();
            ftree->phi.clearAll();
            ftree->dR_Kp_Km  = null;
            ftree->dR_Kp_phi = null;
            ftree->dR_Km_phi = null;
            ftree->phiFit.clearAll();
            ftree->phiFit_Kp.clearAll();
            ftree->phiFit_Km.clearAll();
            ftree->phiFit_phi.clearAll();
            ftree->phiFit_dR_Kp_Km  = null;
            ftree->phiFit_dR_Kp_phi = null;
            ftree->phiFit_dR_Km_phi = null;

            const auto& Km_PF = (*packedHandle)[idx_Km_vec[j]];
            ftree->Km.fillAll( Km_PF );
            TLorentzVector Km_P4;
            Km_P4.SetXYZM(Km_PF.px(), Km_PF.py(), Km_PF.pz(), Mass_K);

            TLorentzVector phi_P4 = Kp_P4 + Km_P4;
            ftree->phi.fillAll( phi_P4 ); 

            ftree->dR_Kp_Km  = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Km_P4.Eta(),  Km_P4.Phi());
            ftree->dR_Kp_phi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());
            ftree->dR_Km_phi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());

            if( abs( Kp_PF.pt() - Km_PF.pt() ) > 20 ) continue;
            if( ftree->dR_Kp_Km > 0.1 ) continue;

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*ttBuilder).build(Kp_PF.pseudoTrack()),
                (*ttBuilder).build(Km_PF.pseudoTrack())
            };
            TransientVertex phi_Vertex = KVFitter.vertex(phi_Tracks);

            if( !(phi_Vertex.isValid()) ) continue;
            if( !(phi_Vertex.hasRefittedTracks()) ) continue;
            if( phi_Vertex.normalisedChiSquared() < 0 ) continue;
            if( phi_Vertex.normalisedChiSquared() > 10 ) continue;

            std::vector<reco::TransientTrack> phiFit_Tracks = phi_Vertex.refittedTracks();

            ftree->phiFit.fillAll( reco::Vertex(phi_Vertex) );

            TLorentzVector phiFit_Kp_P4;
            TLorentzVector phiFit_Km_P4;

            phiFit_Kp_P4.SetXYZM(phiFit_Tracks[0].track().px(), phiFit_Tracks[0].track().py(), phiFit_Tracks[0].track().pz(), Mass_K);
            phiFit_Km_P4.SetXYZM(phiFit_Tracks[1].track().px(), phiFit_Tracks[1].track().py(), phiFit_Tracks[1].track().pz(), Mass_K);

            TLorentzVector phiFit_phi_P4 = phiFit_Kp_P4 + phiFit_Km_P4;

            ftree->phiFit_Kp.fillAll(  phiFit_Kp_P4 ); 
            ftree->phiFit_Km.fillAll(  phiFit_Km_P4 ); 
            ftree->phiFit_phi.fillAll( phiFit_phi_P4 ); 

            ftree->phiFit_dR_Kp_Km  = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_Km_P4.Eta(),  phiFit_Km_P4.Phi());
            ftree->phiFit_dR_Kp_phi = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());
            ftree->phiFit_dR_Km_phi = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());

            if( ftree->phiFit_dR_Kp_phi > 0.05 ) continue;
            if( ftree->phiFit_dR_Km_phi > 0.05 ) continue;
            if( phiFit_phi_P4.M() < 0.99 ) continue;
            if( phiFit_phi_P4.M() > 1.05 ) continue;

            reco::Track phi_dummy_track(0.0, 0.0, reco::Track::Point(phi_Vertex.position().x(), phi_Vertex.position().y(), phi_Vertex.position().z()), reco::Track::Vector(phiFit_phi_P4.Px(), phiFit_phi_P4.Py(), phiFit_phi_P4.Pz()), 0, dummy_track_cov);

            ftree->num_reco_phi++;

            for(size_t k=0; k<idx_pi_vec.size(); k++){
                if( idx_pi_vec[k] == idx_Kp_vec[i] ) continue;
                ftree->pi.clearAll();
                ftree->Ds.clearAll();
                ftree->phiFit_pi.clearAll();
                ftree->phiFit_Ds.clearAll();
                ftree->dR_Kp_pi  = null;
                ftree->dR_Km_pi  = null;
                ftree->dR_pi_phi = null;
                ftree->dR_Kp_Ds  = null;
                ftree->dR_Km_Ds  = null;
                ftree->dR_pi_Ds  = null;
                ftree->dR_phi_Ds = null;
                ftree->phiFit_dR_Kp_pi  = null;
                ftree->phiFit_dR_Km_pi  = null;
                ftree->phiFit_dR_pi_phi = null;
                ftree->phiFit_dR_Kp_Ds  = null;
                ftree->phiFit_dR_Km_Ds  = null;
                ftree->phiFit_dR_pi_Ds  = null;
                ftree->phiFit_dR_phi_Ds = null;
                ftree->DsFit.clearAll();
                ftree->DsFit_Kp.clearAll();
                ftree->DsFit_Km.clearAll();
                ftree->DsFit_pi.clearAll();
                ftree->DsFit_phi.clearAll();
                ftree->DsFit_Ds.clearAll();
                ftree->Ds_primvtx_FD.clearAll();
                ftree->Kp_primvtx_ip.clearAll();
                ftree->Km_primvtx_ip.clearAll();
                ftree->pi_primvtx_ip.clearAll();
                ftree->phi_primvtx_ip.clearAll();
                ftree->Ds_primvtx_ip.clearAll();
                ftree->PVnoDs.clearAll();
                ftree->Ds_PVnoDs_FD.clearAll();
                ftree->Kp_PVnoDs_ip.clearAll();
                ftree->Km_PVnoDs_ip.clearAll();
                ftree->pi_PVnoDs_ip.clearAll();
                ftree->phi_PVnoDs_ip.clearAll();
                ftree->Ds_PVnoDs_ip.clearAll();
                ftree->PVwithDs.clearAll();
                ftree->Ds_PVwithDs_FD.clearAll();
                ftree->Kp_PVwithDs_ip.clearAll();
                ftree->Km_PVwithDs_ip.clearAll();
                ftree->pi_PVwithDs_ip.clearAll();
                ftree->phi_PVwithDs_ip.clearAll();
                ftree->Ds_PVwithDs_ip.clearAll();
                ftree->Ds_IsoR03.clearAll();
                ftree->Ds_IsoR04.clearAll();
                ftree->DsFit_dR_Kp_Km  = null;
                ftree->DsFit_dR_Kp_pi  = null;
                ftree->DsFit_dR_Km_pi  = null;
                ftree->DsFit_dR_Kp_phi = null;
                ftree->DsFit_dR_Km_phi = null;
                ftree->DsFit_dR_pi_phi = null;
                ftree->DsFit_dR_Kp_Ds  = null;
                ftree->DsFit_dR_Km_Ds  = null;
                ftree->DsFit_dR_pi_Ds  = null;
                ftree->DsFit_dR_phi_Ds = null;
                ftree->dxy_phi_Ds = null;
                ftree->dz_phi_Ds  = null;

                const auto& pi_PF = (*packedHandle)[idx_pi_vec[k]];
                ftree->pi.fillAll( pi_PF );

                TLorentzVector pi_P4;
                pi_P4.SetXYZM(pi_PF.px(), pi_PF.py(), pi_PF.pz(), Mass_pi);

                TLorentzVector Ds_P4 = Kp_P4 + Km_P4 + pi_P4;
                ftree->Ds.fillAll( Ds_P4 );

                ftree->dR_Kp_pi  = reco::deltaR(Kp_P4.Eta(),  Kp_P4.Phi(),  pi_P4.Eta(),  pi_P4.Phi());
                ftree->dR_Km_pi  = reco::deltaR(Km_P4.Eta(),  Km_P4.Phi(),  pi_P4.Eta(),  pi_P4.Phi());
                ftree->dR_pi_phi = reco::deltaR(pi_P4.Eta(),  pi_P4.Phi(),  phi_P4.Eta(), phi_P4.Phi());
                ftree->dR_Kp_Ds  = reco::deltaR(Kp_P4.Eta(),  Kp_P4.Phi(),  Ds_P4.Eta(),  Ds_P4.Phi());
                ftree->dR_Km_Ds  = reco::deltaR(Km_P4.Eta(),  Km_P4.Phi(),  Ds_P4.Eta(),  Ds_P4.Phi());
                ftree->dR_pi_Ds  = reco::deltaR(pi_P4.Eta(),  pi_P4.Phi(),  Ds_P4.Eta(),  Ds_P4.Phi());
                ftree->dR_phi_Ds = reco::deltaR(phi_P4.Eta(), phi_P4.Phi(), Ds_P4.Eta(),  Ds_P4.Phi());

                TLorentzVector phiFit_pi_P4 = pi_P4;
                ftree->phiFit_pi.fillAll( phiFit_pi_P4 );

                TLorentzVector phiFit_Ds_P4 = phiFit_Kp_P4 + phiFit_Km_P4 + phiFit_pi_P4;
                ftree->phiFit_Ds.fillAll( phiFit_Ds_P4 );

                ftree->phiFit_dR_Kp_pi  = reco::deltaR(phiFit_Kp_P4.Eta(),  phiFit_Kp_P4.Phi(),  phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi());
                ftree->phiFit_dR_Km_pi  = reco::deltaR(phiFit_Km_P4.Eta(),  phiFit_Km_P4.Phi(),  phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi());
                ftree->phiFit_dR_pi_phi = reco::deltaR(phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi(),  phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());
                ftree->phiFit_dR_Kp_Ds  = reco::deltaR(phiFit_Kp_P4.Eta(),  phiFit_Kp_P4.Phi(),  phiFit_Ds_P4.Eta(),  phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_Km_Ds  = reco::deltaR(phiFit_Km_P4.Eta(),  phiFit_Km_P4.Phi(),  phiFit_Ds_P4.Eta(),  phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_pi_Ds  = reco::deltaR(phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi(),  phiFit_Ds_P4.Eta(),  phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_phi_Ds = reco::deltaR(phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi(), phiFit_Ds_P4.Eta(),  phiFit_Ds_P4.Phi());

                if( ftree->dR_Kp_pi > 0.4 ) continue;
                if( ftree->dR_Km_pi > 0.4 ) continue;
                if( ftree->phiFit_dR_pi_phi > 0.4 ) continue;

                std::vector<reco::TransientTrack> Ds_Tracks = {
                    phiFit_Tracks[0],
                    phiFit_Tracks[1],
                    (*ttBuilder).build(pi_PF.pseudoTrack())
                };

                TransientVertex Ds_Vertex = KVFitter.vertex(Ds_Tracks);

                if( !(Ds_Vertex.isValid()) ) continue;
                if( !(Ds_Vertex.hasRefittedTracks()) ) continue;
                if( Ds_Vertex.normalisedChiSquared() < 0 ) continue;
                if( Ds_Vertex.normalisedChiSquared() > 10 ) continue;

                ftree->DsFit.fillAll( reco::Vertex(Ds_Vertex) );

                std::vector<reco::TransientTrack> DsFit_Tracks = Ds_Vertex.refittedTracks();

                TLorentzVector DsFit_Kp_P4;
                TLorentzVector DsFit_Km_P4;
                TLorentzVector DsFit_pi_P4;

                DsFit_Kp_P4.SetXYZM(DsFit_Tracks[0].track().px(), DsFit_Tracks[0].track().py(), DsFit_Tracks[0].track().pz(), Mass_K);
                DsFit_Km_P4.SetXYZM(DsFit_Tracks[1].track().px(), DsFit_Tracks[1].track().py(), DsFit_Tracks[1].track().pz(), Mass_K);
                DsFit_pi_P4.SetXYZM(DsFit_Tracks[2].track().px(), DsFit_Tracks[2].track().py(), DsFit_Tracks[2].track().pz(), Mass_pi);

                TLorentzVector DsFit_phi_P4 = DsFit_Kp_P4 + DsFit_Km_P4;
                TLorentzVector DsFit_Ds_P4  = DsFit_Kp_P4 + DsFit_Km_P4 + DsFit_pi_P4;

                ftree->DsFit_Kp.fillAll(  DsFit_Kp_P4  );
                ftree->DsFit_Km.fillAll(  DsFit_Km_P4  );
                ftree->DsFit_pi.fillAll(  DsFit_pi_P4  );
                ftree->DsFit_phi.fillAll( DsFit_phi_P4 );
                ftree->DsFit_Ds.fillAll(  DsFit_Ds_P4  );

                TLorentzVector DsFit_Mconstraint_phi_P4;
                DsFit_Mconstraint_phi_P4.SetXYZM(DsFit_phi_P4.Px(), DsFit_phi_P4.Py(), DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector DsFit_Mconstraint_Ds_P4 = DsFit_Mconstraint_phi_P4 + DsFit_pi_P4;
                ftree->DsFit_Mconstraint_Ds_invm = DsFit_Mconstraint_Ds_P4.M();

                ftree->DsFit_dR_Kp_Km  = reco::deltaR(DsFit_Kp_P4.Eta(),  DsFit_Kp_P4.Phi(),  DsFit_Km_P4.Eta(),  DsFit_Km_P4.Phi());
                ftree->DsFit_dR_Kp_phi = reco::deltaR(DsFit_Kp_P4.Eta(),  DsFit_Kp_P4.Phi(),  DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi());
                ftree->DsFit_dR_Km_phi = reco::deltaR(DsFit_Km_P4.Eta(),  DsFit_Km_P4.Phi(),  DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi());
                ftree->DsFit_dR_Kp_pi  = reco::deltaR(DsFit_Kp_P4.Eta(),  DsFit_Kp_P4.Phi(),  DsFit_pi_P4.Eta(),  DsFit_pi_P4.Phi());
                ftree->DsFit_dR_Km_pi  = reco::deltaR(DsFit_Km_P4.Eta(),  DsFit_Km_P4.Phi(),  DsFit_pi_P4.Eta(),  DsFit_pi_P4.Phi());
                ftree->DsFit_dR_pi_phi = reco::deltaR(DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi(), DsFit_pi_P4.Eta(),  DsFit_pi_P4.Phi());
                ftree->DsFit_dR_Kp_Ds  = reco::deltaR(DsFit_Kp_P4.Eta(),  DsFit_Kp_P4.Phi(),  DsFit_Ds_P4.Eta(),  DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_Km_Ds  = reco::deltaR(DsFit_Km_P4.Eta(),  DsFit_Km_P4.Phi(),  DsFit_Ds_P4.Eta(),  DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_phi_Ds = reco::deltaR(DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi(), DsFit_Ds_P4.Eta(),  DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_pi_Ds  = reco::deltaR(DsFit_pi_P4.Eta(),  DsFit_pi_P4.Phi(),  DsFit_Ds_P4.Eta(),  DsFit_Ds_P4.Phi());

                ftree->dxy_phi_Ds = sqrt( pow( phi_Vertex.position().x() - Ds_Vertex.position().x(), 2 ) + pow( phi_Vertex.position().y() - Ds_Vertex.position().y(), 2 ) );
                ftree->dz_phi_Ds  = abs( phi_Vertex.position().z() - Ds_Vertex.position().z() );

                if( ftree->DsFit_dR_Kp_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_Km_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_pi_Ds > 0.4 ) continue;
                if( ftree->DsFit_dR_phi_Ds > 0.15 ) continue;
                if( DsFit_Ds_P4.M() < 1.85) continue;
                if( DsFit_Ds_P4.M() > 2.1) continue;

                reco::Track Ds_dummy_track(0.0, 0.0, reco::Track::Point(Ds_Vertex.position().x(), Ds_Vertex.position().y(), Ds_Vertex.position().z()), reco::Track::Vector(DsFit_Ds_P4.Px(), DsFit_Ds_P4.Py(), DsFit_Ds_P4.Pz()), 0, dummy_track_cov);

                ftree->num_reco_Ds++;

                GlobalPoint Ds_pos = Ds_Vertex.position();
                GlobalError Ds_cov = Ds_Vertex.positionError();
                GlobalVector Ds_PDirec(DsFit_Ds_P4.Px(), DsFit_Ds_P4.Py(), DsFit_Ds_P4.Pz());

                ftree->Ds_primvtx_FD.fillAll(Ds_pos-primvtx_pos, Ds_cov+primvtx_cov, Ds_PDirec);
                ftree->Kp_primvtx_ip.fillAll(DsFit_Tracks[0], primvtx);
                ftree->Km_primvtx_ip.fillAll(DsFit_Tracks[1], primvtx);
                ftree->pi_primvtx_ip.fillAll(DsFit_Tracks[2], primvtx);
                ftree->phi_primvtx_ip.fillAll((*ttBuilder).build(phi_dummy_track), primvtx);
                ftree->Ds_primvtx_ip.fillAll( (*ttBuilder).build(Ds_dummy_track),  primvtx);

                reco::TrackCollection PVtracks_noDs;
                reco::TrackCollection PVtracks_withDs;

                for(size_t it=0; it<packedHandle->size(); it++){
                    if(it==static_cast<size_t>(idx_Kp_vec[i])) continue;
                    if(it==static_cast<size_t>(idx_Km_vec[j])) continue;
                    if(it==static_cast<size_t>(idx_pi_vec[k])) continue;
                    const auto& pf = (*packedHandle)[it];

                    if(reco::deltaR(DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.3){
                        if( abs(pf.pdgId()) == 211 && (pf.fromPV()>1) ) ftree->Ds_IsoR03.sumChargedHadronPt += pf.pt();
                        else if( pf.pdgId() == 130 ) ftree->Ds_IsoR03.sumNeutralHadronEt += pf.et();
                        else if( pf.pdgId() == 22 ) ftree->Ds_IsoR03.sumPhotonEt += pf.et();

                        if( pf.fromPV() <= 1 ) ftree->Ds_IsoR03.sumPUPt += pf.pt();
                    }

                    if(reco::deltaR(DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.4){
                        if( abs(pf.pdgId()) == 211 && (pf.fromPV()>1) ) ftree->Ds_IsoR04.sumChargedHadronPt += pf.pt();
                        else if( pf.pdgId() == 130 ) ftree->Ds_IsoR04.sumNeutralHadronEt += pf.et();
                        else if( pf.pdgId() == 22 ) ftree->Ds_IsoR04.sumPhotonEt += pf.et();

                        if( pf.fromPV() <= 1 ) ftree->Ds_IsoR04.sumPUPt += pf.pt();
                    }

                    if( ! pf.hasTrackDetails() ) continue;
                    PVtracks_noDs.push_back(pf.pseudoTrack());
                    PVtracks_withDs.push_back(pf.pseudoTrack());
                }
                PVtracks_withDs.push_back(Ds_dummy_track);

                ftree->Ds_IsoR03.PFIso = ( ftree->Ds_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->Ds_IsoR03.sumNeutralHadronEt + ftree->Ds_IsoR03.sumPhotonEt - 0.5*ftree->Ds_IsoR03.sumPUPt) ) / DsFit_Ds_P4.Pt();
                ftree->Ds_IsoR04.PFIso = ( ftree->Ds_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->Ds_IsoR04.sumNeutralHadronEt + ftree->Ds_IsoR04.sumPhotonEt - 0.5*ftree->Ds_IsoR04.sumPUPt) ) / DsFit_Ds_P4.Pt();

                std::vector<TransientVertex> pvs_noDs = revertex->makeVertices(PVtracks_noDs, *bsHandle, iSetup, "noBS");
                if( !(pvs_noDs.empty()) ){
                    reco::Vertex pv_noDs = reco::Vertex(pvs_noDs.front());
                    GlobalPoint PVnoDs_pos = pvs_noDs.front().position();
                    GlobalError PVnoDs_cov = pvs_noDs.front().positionError();

                    ftree->PVnoDs.fillAll( pv_noDs);
                    ftree->Ds_PVnoDs_FD.fillAll(Ds_pos-PVnoDs_pos, Ds_cov+PVnoDs_cov, Ds_PDirec);
                    ftree->Kp_PVnoDs_ip.fillAll(DsFit_Tracks[0], pv_noDs);
                    ftree->Km_PVnoDs_ip.fillAll(DsFit_Tracks[1], pv_noDs);
                    ftree->pi_PVnoDs_ip.fillAll(DsFit_Tracks[2], pv_noDs);
                    ftree->phi_PVnoDs_ip.fillAll((*ttBuilder).build(phi_dummy_track), pv_noDs);
                    ftree->Ds_PVnoDs_ip.fillAll( (*ttBuilder).build(Ds_dummy_track),  pv_noDs);
                }

                std::vector<TransientVertex> pvs_withDs = revertex->makeVertices(PVtracks_withDs, *bsHandle, iSetup, "noBS");
                if( !(pvs_withDs.empty()) ){
                    reco::Vertex pv_withDs = reco::Vertex(pvs_withDs.front());
                    GlobalPoint PVwithDs_pos = pvs_withDs.front().position();
                    GlobalError PVwithDs_cov = pvs_withDs.front().positionError();

                    ftree->PVwithDs.fillAll( pv_withDs);
                    ftree->Ds_PVwithDs_FD.fillAll(Ds_pos-PVwithDs_pos, Ds_cov+PVwithDs_cov, Ds_PDirec);
                    ftree->Kp_PVwithDs_ip.fillAll(DsFit_Tracks[0], pv_withDs);
                    ftree->Km_PVwithDs_ip.fillAll(DsFit_Tracks[1], pv_withDs);
                    ftree->pi_PVwithDs_ip.fillAll(DsFit_Tracks[2], pv_withDs);
                    ftree->phi_PVwithDs_ip.fillAll((*ttBuilder).build(phi_dummy_track), pv_withDs);
                    ftree->Ds_PVwithDs_ip.fillAll( (*ttBuilder).build(Ds_dummy_track),  pv_withDs);

                }

                ftree->FillDs();
            }
        }
    }

    double max_Ds_pt = null;
    int maxidx_Ds_pt = -1;

    for(int i=0; i<ftree->num_reco_Ds; i++){
        if(ftree->DsFit_Ds_vec.pt[i]>max_Ds_pt){
            max_Ds_pt = ftree->DsFit_Ds_vec.pt[i];
            maxidx_Ds_pt = i;
        }
    }
    if( maxidx_Ds_pt < 0 ) return;

    ftree->FillBestDs(maxidx_Ds_pt);

    double max_mu_pt = null;
    int maxidx_mu_pt = -1;

    for(size_t i=0; i<muonHandle->size(); i++){
        const auto& muon = (*muonHandle)[i];
        if(muon.charge() != -1) continue;
        if(muon.pt() < 0.5) continue;
        if(muon.p() < 1) continue;
        if(!muon.isTightMuon(primvtx)) continue;

        if(muon.pt()>max_mu_pt){
            max_mu_pt = muon.pt();
            maxidx_mu_pt = i;
        }
    }

    if( maxidx_mu_pt < 0 ) return;
    ftree->mu.clearAll();
    ftree->mu_primvtx_ip.clearAll();

    const auto& muon = (*muonHandle)[maxidx_mu_pt];

    ftree->mu.fillAll( muon, primvtx );
    ftree->mu_IsoR03.sumChargedHadronPt = muon.pfIsolationR03().sumChargedHadronPt;
    ftree->mu_IsoR03.sumNeutralHadronEt = muon.pfIsolationR03().sumNeutralHadronEt;
    ftree->mu_IsoR03.sumPhotonEt = muon.pfIsolationR03().sumPhotonEt;
    ftree->mu_IsoR03.sumPUPt = muon.pfIsolationR03().sumPUPt;
    ftree->mu_IsoR03.PFIso = (ftree->mu_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR03.sumNeutralHadronEt + ftree->mu_IsoR03.sumPhotonEt - 0.5*ftree->mu_IsoR03.sumPUPt)) / muon.pt();
    ftree->mu_IsoR04.sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
    ftree->mu_IsoR04.sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
    ftree->mu_IsoR04.sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
    ftree->mu_IsoR04.sumPUPt = muon.pfIsolationR04().sumPUPt;
    ftree->mu_IsoR04.PFIso = (ftree->mu_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR04.sumNeutralHadronEt + ftree->mu_IsoR04.sumPhotonEt - 0.5*ftree->mu_IsoR04.sumPUPt)) / muon.pt();
    ftree->mu_primvtx_ip.fillAll( (*ttBuilder).build(muon.muonBestTrack()), primvtx);

    ftree->FillBestmu();

    ftree->tree->Fill();

    return;
}

void RecoAnalyzer::endJob()
{
    edm::LogInfo("RecoAnalyzer") << "Processed all events.";
    edm::LogPrint("RecoTuplizer") << "Total number of events processed: " << nEventsProcessed_;
}

void RecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoAnalyzer);

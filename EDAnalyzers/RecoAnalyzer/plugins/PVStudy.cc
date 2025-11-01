// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      PVStudy
//
/**\class PVStudy PVStudy.cc EDAnalyzer/GenParticleAnalyzer/plugins/PVStudy.cc

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

#include "EDAnalyzers/RecoAnalyzer/interface/PVTree.h"
#include "EDAnalyzers/RecoAnalyzer/interface/VertexReProducer.h"

class PVStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit PVStudy(const edm::ParameterSet& iConfig);
        ~PVStudy();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken;
        edm::EDGetTokenT<pat::MuonCollection> muonToken;
        edm::EDGetTokenT<reco::VertexCollection> pvToken;
        edm::EDGetTokenT<reco::BeamSpot> bsToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        KalmanVertexFitter KVFitter;
        VertexReProducer *revertex;

        const edm::Service<TFileService> fs;
        PVTree *ftree;

        std::vector<int> Gen_Kp_idx;
        std::vector<int> Gen_Km_idx;
        std::vector<int> Gen_pi_idx;
        std::vector<int> Gen_phi_idx;
        std::vector<int> Gen_Ds_idx;
        std::vector<int> Gen_mu_idx;
        std::vector<int> Gen_nu_idx;
        std::vector<int> Gen_W_idx;
        std::vector<int> Gen_H_idx;

        std::vector<int> idx_Kp_vec;
        std::vector<int> idx_Km_vec;
        std::vector<int> idx_pi_vec;
        std::vector<int> idx_mu_vec;

        const double Mass_Ds = 1.96835;
        const double Mass_phi = 1.019460;
        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_mu = 0.1056583755;
};

PVStudy::PVStudy(const edm::ParameterSet& iConfig) :
    prunedToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPart"))),
    packedToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    pvToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primvtx"))),
    bsToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    KVFitter(true)
{
    revertex = new VertexReProducer(iConfig);
}

PVStudy::~PVStudy() {}

void PVStudy::beginJob()
{
    TFile & f = fs->file();
    f.SetCompressionAlgorithm(ROOT::kZLIB);
    ftree = new PVTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

void PVStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Gen_Kp_idx.clear();
    Gen_Km_idx.clear();
    Gen_pi_idx.clear();
    Gen_phi_idx.clear();
    Gen_Ds_idx.clear();
    Gen_mu_idx.clear();
    Gen_nu_idx.clear();
    Gen_W_idx.clear();
    Gen_H_idx.clear();

    ftree->Init();

    edm::Handle<reco::GenParticleCollection> prunedHandle;
    iEvent.getByToken(prunedToken, prunedHandle);
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

    for(size_t i=0; i<prunedHandle->size(); i++){
        const auto& gp = (*prunedHandle)[i];
        const int pdg = gp.pdgId();
        if( pdg==321 && hasAncestor(gp, 333, -321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_Kp_idx.push_back(i);
            ftree->num_Gen_Kp++;
        } else if( pdg==-321 && hasAncestor(gp, 333, 321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_Km_idx.push_back(i);
            ftree->num_Gen_Km++;
        } else if( pdg==211 && hasAncestor(gp, 431, 333) && hasAncestor(gp, 25, -24) ){
            Gen_pi_idx.push_back(i);
            ftree->num_Gen_pi++;
        } else if( pdg==333 && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) && gp.numberOfDaughters()==2 ){
            int pdg0 = gp.daughter(0)->pdgId();
            int pdg1 = gp.daughter(1)->pdgId();
            if( (pdg0 == 321 && pdg1 == -321) || (pdg0 == -321 && pdg1 == 321)){
                Gen_phi_idx.push_back(i);
                ftree->num_Gen_phi++;
            }
        } else if( pdg==431 && hasAncestor(gp, 25, -24) && gp.isHardProcess() ){
            Gen_Ds_idx.push_back(i);
            ftree->num_Gen_Ds++;
        } else if( pdg==13 && hasAncestor(gp, -24, -14) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_mu_idx.push_back(i);
            ftree->num_Gen_mu++;
        } else if( pdg==-14 && hasAncestor(gp, -24, 13) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_nu_idx.push_back(i);
            ftree->num_Gen_nu++;
        } else if( pdg==-24 && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_W_idx.push_back(i);
            ftree->num_Gen_W++;
        } else if( pdg==25 && gp.isHardProcess() ){
            Gen_H_idx.push_back(i);
            ftree->num_Gen_H++;
        }
    }

    ftree->Gen_Reset();

    if(ftree->num_Gen_Kp  != 1) return;
    if(ftree->num_Gen_Km  != 1) return;
    if(ftree->num_Gen_pi  != 1) return;
    if(ftree->num_Gen_phi != 1) return;
    if(ftree->num_Gen_Ds  != 1) return;
    if(ftree->num_Gen_mu  != 1) return;
    if(ftree->num_Gen_nu  != 1) return;
    if(ftree->num_Gen_W   != 1) return;
    if(ftree->num_Gen_H   != 1) return;

    const auto& Kp_GP  = (*prunedHandle)[Gen_Kp_idx[0]];
    const auto& Km_GP  = (*prunedHandle)[Gen_Km_idx[0]];
    const auto& pi_GP  = (*prunedHandle)[Gen_pi_idx[0]];
    const auto& phi_GP = (*prunedHandle)[Gen_phi_idx[0]];
    const auto& Ds_GP  = (*prunedHandle)[Gen_Ds_idx[0]];
    const auto& mu_GP  = (*prunedHandle)[Gen_mu_idx[0]];
    const auto& nu_GP  = (*prunedHandle)[Gen_nu_idx[0]];
    const auto& W_GP   = (*prunedHandle)[Gen_W_idx[0]];
    const auto& H_GP   = (*prunedHandle)[Gen_H_idx[0]];

    ftree->Gen_Kp.fillAll(Kp_GP);
    ftree->Gen_Km.fillAll(Km_GP);
    ftree->Gen_pi.fillAll(pi_GP);
    ftree->Gen_phi.fillAll(phi_GP);
    ftree->Gen_Ds.fillAll(Ds_GP);
    ftree->Gen_mu.fillAll(mu_GP);
    ftree->Gen_nu.fillAll(nu_GP);
    ftree->Gen_W.fillAll(W_GP);
    ftree->Gen_H.fillAll(H_GP);

    TVector3 Gen_Kp_P3 (Kp_GP.px(),  Kp_GP.py(),  Kp_GP.pz());
    TVector3 Gen_Km_P3 (Km_GP.px(),  Km_GP.py(),  Km_GP.pz());
    TVector3 Gen_pi_P3 (pi_GP.px(),  pi_GP.py(),  pi_GP.pz());
    TVector3 Gen_phi_P3(phi_GP.px(), phi_GP.py(), phi_GP.pz());
    TVector3 Gen_Ds_P3 (Ds_GP.px(),  Ds_GP.py(),  Ds_GP.pz());
    TVector3 Gen_mu_P3 (mu_GP.px(),  mu_GP.py(),  mu_GP.pz());
    TVector3 Gen_nu_P3 (nu_GP.px(),  nu_GP.py(),  nu_GP.pz());
    TVector3 Gen_W_P3  (W_GP.px(),   W_GP.py(),   W_GP.pz());
    TVector3 Gen_H_P3  (H_GP.px(),   H_GP.py(),   H_GP.pz());

    ftree->Gen_dR_Kp_Km  = reco::deltaR(Kp_GP,  Km_GP);
    ftree->Gen_dR_Kp_phi = reco::deltaR(Kp_GP,  phi_GP);
    ftree->Gen_dR_Km_phi = reco::deltaR(Km_GP,  phi_GP);
    ftree->Gen_dR_Kp_pi  = reco::deltaR(Kp_GP,  pi_GP);
    ftree->Gen_dR_Km_pi  = reco::deltaR(Km_GP,  pi_GP);
    ftree->Gen_dR_pi_phi = reco::deltaR(pi_GP,  phi_GP);
    ftree->Gen_dR_Kp_Ds  = reco::deltaR(Kp_GP,  Ds_GP);
    ftree->Gen_dR_Km_Ds  = reco::deltaR(Km_GP,  Ds_GP);
    ftree->Gen_dR_phi_Ds = reco::deltaR(phi_GP, Ds_GP);
    ftree->Gen_dR_pi_Ds  = reco::deltaR(pi_GP,  Ds_GP);
    ftree->Gen_dR_Kp_mu  = reco::deltaR(Kp_GP,  mu_GP);
    ftree->Gen_dR_Km_mu  = reco::deltaR(Km_GP,  mu_GP);
    ftree->Gen_dR_phi_mu = reco::deltaR(phi_GP, mu_GP);
    ftree->Gen_dR_pi_mu  = reco::deltaR(pi_GP,  mu_GP);
    ftree->Gen_dR_Ds_mu  = reco::deltaR(Ds_GP,  mu_GP);

    ftree->Gen_Ds_dx   = Kp_GP.vx() - H_GP.vx();
    ftree->Gen_Ds_dy   = Kp_GP.vy() - H_GP.vy();
    ftree->Gen_Ds_dz   = Kp_GP.vz() - H_GP.vz();
    ftree->Gen_Ds_FDxy = std::hypot( ftree->Gen_Ds_dx,   ftree->Gen_Ds_dy );
    ftree->Gen_Ds_FD   = std::hypot( ftree->Gen_Ds_FDxy, ftree->Gen_Ds_dz );

    ftree->beamspotInfo.fillAll(*bsHandle );

    ftree->primvtxInfo.clearAll();
    const reco::Vertex& primvtx = pvHandle->front();
    if(primvtx.isValid()) ftree->primvtxInfo.fillAll( primvtx );
    GlobalPoint primvtx_pos( primvtx.x(), primvtx.y(), primvtx.z() );
    GlobalError primvtx_cov( primvtx.covariance(0,0), primvtx.covariance(0,1), primvtx.covariance(0,2), primvtx.covariance(1,1), primvtx.covariance(1,2), primvtx.covariance(2,2) );

    ftree->Match_Reset();
    reco::Track::CovarianceMatrix dummy_track_cov;
    for (unsigned int ic = 0; ic < 5; ++ic) {
        dummy_track_cov(ic,ic) = 1.0;
    }

    struct MatchInfo {int index; float dr;};
    std::vector<MatchInfo> MatchKp, MatchKm, Matchpi, Matchmu;

    for(size_t i=0; i<packedHandle->size(); i++){

        const auto& pf = (*packedHandle)[i];
        if( !pf.trackHighPurity() ) continue;
        if( !pf.hasTrackDetails() ) continue;
        if( abs(pf.pdgId()) != 211 ) continue;

        float dpt_Kp = abs( pf.pt() - Kp_GP.pt() ) / Kp_GP.pt();
        float dpt_Km = abs( pf.pt() - Km_GP.pt() ) / Km_GP.pt();
        float dpt_pi = abs( pf.pt() - pi_GP.pt() ) / pi_GP.pt();

        float dr_Kp = reco::deltaR( pf.eta(), pf.phi(), Kp_GP.eta(), Kp_GP.phi() );
        float dr_Km = reco::deltaR( pf.eta(), pf.phi(), Km_GP.eta(), Km_GP.phi() );
        float dr_pi = reco::deltaR( pf.eta(), pf.phi(), pi_GP.eta(), pi_GP.phi() );

        if( pf.charge() == 1 ){
            if( (dpt_Kp < 0.3) && (dr_Kp < 0.1) ) MatchKp.push_back( {static_cast<int>(i), dr_Kp} );
            if( (dpt_pi < 0.3) && (dr_pi < 0.1) ) Matchpi.push_back( {static_cast<int>(i), dr_pi} );
        }
        else if( pf.charge() == -1 ){
            if( (dpt_Km < 0.3) && (dr_Km < 0.1) ) MatchKm.push_back( {static_cast<int>(i), dr_Km} );
        }
    }

    for(size_t i=0; i<muonHandle->size(); i++){

        const auto& muon = (*muonHandle)[i];
        if (!muon.genParticle()) continue;

        const reco::GenParticle& genMuon = *(muon.genParticle());
        if( !hasAncestor(genMuon, -24, -14) ) continue;
        if( !hasAncestor(genMuon, 25, 431) ) continue;

        float dr_mu = reco::deltaR(muon.eta(), muon.phi(), mu_GP.eta(), mu_GP.phi());
        Matchmu.push_back({static_cast<int>(i), dr_mu});
    }

    auto findBestMatch = [](const std::vector<MatchInfo>& matches) -> MatchInfo {
        auto it = std::min_element(matches.begin(), matches.end(),
                [](const auto& a, const auto& b) { return a.dr < b.dr; });
        return (it != matches.end()) ? *it : MatchInfo{-1, 100};
    };

    auto findBestMatchExcluding = [](const std::vector<MatchInfo>& matches, int excludeIndex) -> MatchInfo {
        MatchInfo best{-1, std::numeric_limits<float>::max()};
        for (const auto& m : matches) {
            if (m.index != excludeIndex && m.dr < best.dr) {
                best = m;
            }
        }
        return (best.dr < std::numeric_limits<float>::max()) ? best : MatchInfo{-1, 100};
    };

    auto bestKp = findBestMatch(MatchKp);
    auto bestKm = findBestMatch(MatchKm);
    auto bestpi = findBestMatch(Matchpi);
    auto bestmu = findBestMatch(Matchmu);

    if (bestKp.index == bestpi.index && bestKp.index != -1){
        if (bestKp.dr <= bestpi.dr) {
            bestpi = findBestMatchExcluding(Matchpi, bestKp.index);
        }
        else{
            const auto altKp = findBestMatchExcluding(MatchKp, bestpi.index);
            if (altKp.index != -1) {
                bestKp = altKp;
            }
            else {
                bestpi = findBestMatchExcluding(Matchpi, bestKp.index);
            }
        }
    }

    if(bestKp.dr < 0.1){
        ftree->num_match_Kp++;
        ftree->match_dR_Kp = bestKp.dr; 
        if(bestKp.dr < 0.03) ftree->num_tight_match_Kp++;
    }
    if(bestKm.dr < 0.1){
        ftree->num_match_Km++;
        ftree->match_dR_Km = bestKm.dr; 
        if(bestKm.dr < 0.03) ftree->num_tight_match_Km++;
    }
    if(bestpi.dr < 0.1){
        ftree->num_match_pi++;
        ftree->match_dR_pi = bestpi.dr; 
        if(bestpi.dr < 0.03) ftree->num_tight_match_pi++;
    }
    if(bestmu.dr < 100){
        ftree->num_match_mu++;
        ftree->match_dR_mu = bestmu.dr; 
    }

    if( (ftree->num_match_Kp > 0) && (ftree->num_match_Km > 0) && (ftree->num_match_pi > 0) ){

        const auto& match_Kp_PF = (*packedHandle)[bestKp.index];
        const auto& match_Km_PF = (*packedHandle)[bestKm.index];
        const auto& match_pi_PF = (*packedHandle)[bestpi.index];

        ftree->match_Kp.fillAll( match_Kp_PF );
        ftree->match_Km.fillAll( match_Km_PF );
        ftree->match_pi.fillAll( match_pi_PF );

        TLorentzVector match_Kp_P4;
        TLorentzVector match_Km_P4;
        TLorentzVector match_pi_P4;
        match_Kp_P4.SetXYZM(match_Kp_PF.px(), match_Kp_PF.py(), match_Kp_PF.pz(), Mass_K);
        match_Km_P4.SetXYZM(match_Km_PF.px(), match_Km_PF.py(), match_Km_PF.pz(), Mass_K);
        match_pi_P4.SetXYZM(match_pi_PF.px(), match_pi_PF.py(), match_pi_PF.pz(), Mass_pi);

        TLorentzVector match_phi_P4 = match_Kp_P4 + match_Km_P4;
        TLorentzVector match_Ds_P4 = match_Kp_P4 + match_Km_P4 + match_pi_P4;

        ftree->match_phi.fillAll( match_phi_P4 );
        ftree->match_Ds.fillAll(  match_Ds_P4  );

        ftree->match_dR_Kp_Km  = reco::deltaR(match_Kp_P4.Eta(),  match_Kp_P4.Phi(),  match_Km_P4.Eta(),  match_Km_P4.Phi());
        ftree->match_dR_Kp_phi = reco::deltaR(match_Kp_P4.Eta(),  match_Kp_P4.Phi(),  match_phi_P4.Eta(), match_phi_P4.Phi());
        ftree->match_dR_Km_phi = reco::deltaR(match_Km_P4.Eta(),  match_Km_P4.Phi(),  match_phi_P4.Eta(), match_phi_P4.Phi());
        ftree->match_dR_Kp_pi  = reco::deltaR(match_Kp_P4.Eta(),  match_Kp_P4.Phi(),  match_pi_P4.Eta(),  match_pi_P4.Phi());
        ftree->match_dR_Km_pi  = reco::deltaR(match_Km_P4.Eta(),  match_Km_P4.Phi(),  match_pi_P4.Eta(),  match_pi_P4.Phi());
        ftree->match_dR_pi_phi = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_pi_P4.Eta(),  match_pi_P4.Phi());
        ftree->match_dR_Kp_Ds  = reco::deltaR(match_Kp_P4.Eta(),  match_Kp_P4.Phi(),  match_Ds_P4.Eta(),  match_Ds_P4.Phi());
        ftree->match_dR_Km_Ds  = reco::deltaR(match_Km_P4.Eta(),  match_Km_P4.Phi(),  match_Ds_P4.Eta(),  match_Ds_P4.Phi());
        ftree->match_dR_phi_Ds = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_Ds_P4.Eta(),  match_Ds_P4.Phi());
        ftree->match_dR_pi_Ds  = reco::deltaR(match_pi_P4.Eta(),  match_pi_P4.Phi(),  match_Ds_P4.Eta(),  match_Ds_P4.Phi());

        std::vector<reco::TransientTrack> match_phi_Tracks = {
            (*ttBuilder).build(match_Kp_PF.pseudoTrack()),
            (*ttBuilder).build(match_Km_PF.pseudoTrack())
        };

        TransientVertex match_phi_Vertex = KVFitter.vertex(match_phi_Tracks);

        if( match_phi_Vertex.isValid() && match_phi_Vertex.hasRefittedTracks() ){

            ftree->match_phiFit.fillAll( reco::Vertex(match_phi_Vertex) );

            std::vector<reco::TransientTrack> match_phiFit_Tracks = match_phi_Vertex.refittedTracks();

            TLorentzVector match_phiFit_Kp_P4;
            TLorentzVector match_phiFit_Km_P4;
            match_phiFit_Kp_P4.SetXYZM(match_phiFit_Tracks[0].track().px(), match_phiFit_Tracks[0].track().py(), match_phiFit_Tracks[0].track().pz(), Mass_K);
            match_phiFit_Km_P4.SetXYZM(match_phiFit_Tracks[1].track().px(), match_phiFit_Tracks[1].track().py(), match_phiFit_Tracks[1].track().pz(), Mass_K);

            TLorentzVector match_phiFit_pi_P4  = match_pi_P4;
            TLorentzVector match_phiFit_phi_P4 = match_phiFit_Kp_P4 + match_phiFit_Km_P4;
            TLorentzVector match_phiFit_Ds_P4  = match_phiFit_Kp_P4 + match_phiFit_Km_P4 + match_phiFit_pi_P4;

            ftree->match_phiFit_Kp.fillAll(  match_phiFit_Kp_P4  );
            ftree->match_phiFit_Km.fillAll(  match_phiFit_Km_P4  );
            ftree->match_phiFit_pi.fillAll(  match_phiFit_pi_P4  );
            ftree->match_phiFit_phi.fillAll( match_phiFit_phi_P4 );
            ftree->match_phiFit_Ds.fillAll(  match_phiFit_Ds_P4  );

            ftree->match_phiFit_dR_Kp_Km  = reco::deltaR(match_phiFit_Kp_P4.Eta(),  match_phiFit_Kp_P4.Phi(),  match_phiFit_Km_P4.Eta(),  match_phiFit_Km_P4.Phi());
            ftree->match_phiFit_dR_Kp_phi = reco::deltaR(match_phiFit_Kp_P4.Eta(),  match_phiFit_Kp_P4.Phi(),  match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi());
            ftree->match_phiFit_dR_Km_phi = reco::deltaR(match_phiFit_Km_P4.Eta(),  match_phiFit_Km_P4.Phi(),  match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi());
            ftree->match_phiFit_dR_Kp_pi  = reco::deltaR(match_phiFit_Kp_P4.Eta(),  match_phiFit_Kp_P4.Phi(),  match_phiFit_pi_P4.Eta(),  match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_Km_pi  = reco::deltaR(match_phiFit_Km_P4.Eta(),  match_phiFit_Km_P4.Phi(),  match_phiFit_pi_P4.Eta(),  match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_pi_phi = reco::deltaR(match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi(), match_phiFit_pi_P4.Eta(),  match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_Kp_Ds  = reco::deltaR(match_phiFit_Kp_P4.Eta(),  match_phiFit_Kp_P4.Phi(),  match_phiFit_Ds_P4.Eta(),  match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_Km_Ds  = reco::deltaR(match_phiFit_Km_P4.Eta(),  match_phiFit_Km_P4.Phi(),  match_phiFit_Ds_P4.Eta(),  match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_phi_Ds = reco::deltaR(match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi(), match_phiFit_Ds_P4.Eta(),  match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_pi_Ds  = reco::deltaR(match_phiFit_pi_P4.Eta(),  match_phiFit_pi_P4.Phi(),  match_phiFit_Ds_P4.Eta(),  match_phiFit_Ds_P4.Phi());

            reco::Track match_phi_dummy_track(0.0, 0.0, reco::Track::Point(match_phi_Vertex.position().x(), match_phi_Vertex.position().y(), match_phi_Vertex.position().z()), reco::Track::Vector(match_phiFit_phi_P4.Px(), match_phiFit_phi_P4.Py(), match_phiFit_phi_P4.Pz()), 0, dummy_track_cov);

            std::vector<reco::TransientTrack> match_Ds_Tracks = {
                match_phiFit_Tracks[0],
                match_phiFit_Tracks[1],
                (*ttBuilder).build(match_pi_PF.pseudoTrack())
            };

            TransientVertex match_Ds_Vertex = KVFitter.vertex(match_Ds_Tracks);

            if( match_Ds_Vertex.isValid() && match_Ds_Vertex.hasRefittedTracks() ){

                ftree->match_DsFit.fillAll( reco::Vertex(match_Ds_Vertex) );

                std::vector<reco::TransientTrack> match_DsFit_Tracks = match_Ds_Vertex.refittedTracks();

                TLorentzVector match_DsFit_Kp_P4;
                TLorentzVector match_DsFit_Km_P4;
                TLorentzVector match_DsFit_pi_P4;
                match_DsFit_Kp_P4.SetXYZM(match_DsFit_Tracks[0].track().px(), match_DsFit_Tracks[0].track().py(), match_DsFit_Tracks[0].track().pz(), Mass_K);
                match_DsFit_Km_P4.SetXYZM(match_DsFit_Tracks[1].track().px(), match_DsFit_Tracks[1].track().py(), match_DsFit_Tracks[1].track().pz(), Mass_K);
                match_DsFit_pi_P4.SetXYZM(match_DsFit_Tracks[2].track().px(), match_DsFit_Tracks[2].track().py(), match_DsFit_Tracks[2].track().pz(), Mass_pi);

                TLorentzVector match_DsFit_phi_P4 = match_DsFit_Kp_P4 + match_DsFit_Km_P4;
                TLorentzVector match_DsFit_Ds_P4  = match_DsFit_Kp_P4 + match_DsFit_Km_P4 + match_DsFit_pi_P4;

                ftree->match_DsFit_Kp.fillAll(  match_DsFit_Kp_P4  );
                ftree->match_DsFit_Km.fillAll(  match_DsFit_Km_P4  );
                ftree->match_DsFit_pi.fillAll(  match_DsFit_pi_P4  );
                ftree->match_DsFit_phi.fillAll( match_DsFit_phi_P4 );
                ftree->match_DsFit_Ds.fillAll(  match_DsFit_Ds_P4  );

                TLorentzVector match_DsFit_Mconstraint_phi_P4;
                match_DsFit_Mconstraint_phi_P4.SetXYZM(match_DsFit_phi_P4.Px(), match_DsFit_phi_P4.Py(), match_DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector match_DsFit_Mconstraint_Ds_P4 = match_DsFit_Mconstraint_phi_P4 + match_DsFit_pi_P4;
                ftree->match_DsFit_Mconstraint_Ds_invm = match_DsFit_Mconstraint_Ds_P4.M();

                ftree->match_DsFit_dR_Kp_Km  = reco::deltaR(match_DsFit_Kp_P4.Eta(),  match_DsFit_Kp_P4.Phi(),  match_DsFit_Km_P4.Eta(),  match_DsFit_Km_P4.Phi());
                ftree->match_DsFit_dR_Kp_phi = reco::deltaR(match_DsFit_Kp_P4.Eta(),  match_DsFit_Kp_P4.Phi(),  match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi());
                ftree->match_DsFit_dR_Km_phi = reco::deltaR(match_DsFit_Km_P4.Eta(),  match_DsFit_Km_P4.Phi(),  match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi());
                ftree->match_DsFit_dR_Kp_pi  = reco::deltaR(match_DsFit_Kp_P4.Eta(),  match_DsFit_Kp_P4.Phi(),  match_DsFit_pi_P4.Eta(),  match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_Km_pi  = reco::deltaR(match_DsFit_Km_P4.Eta(),  match_DsFit_Km_P4.Phi(),  match_DsFit_pi_P4.Eta(),  match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_pi_phi = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_pi_P4.Eta(),  match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_Kp_Ds  = reco::deltaR(match_DsFit_Kp_P4.Eta(),  match_DsFit_Kp_P4.Phi(),  match_DsFit_Ds_P4.Eta(),  match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_Km_Ds  = reco::deltaR(match_DsFit_Km_P4.Eta(),  match_DsFit_Km_P4.Phi(),  match_DsFit_Ds_P4.Eta(),  match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_phi_Ds = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_Ds_P4.Eta(),  match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_pi_Ds  = reco::deltaR(match_DsFit_pi_P4.Eta(),  match_DsFit_pi_P4.Phi(),  match_DsFit_Ds_P4.Eta(),  match_DsFit_Ds_P4.Phi());

                ftree->match_dxy_phi_Ds = std::hypot( match_phi_Vertex.position().x() - match_Ds_Vertex.position().x(), match_phi_Vertex.position().y() - match_Ds_Vertex.position().y() );
                ftree->match_dz_phi_Ds  = abs( match_phi_Vertex.position().z() - match_Ds_Vertex.position().z() );

                reco::Track match_Ds_dummy_track(0.0, 0.0, reco::Track::Point(match_Ds_Vertex.position().x(), match_Ds_Vertex.position().y(), match_Ds_Vertex.position().z()), reco::Track::Vector(match_DsFit_Ds_P4.Px(), match_DsFit_Ds_P4.Py(), match_DsFit_Ds_P4.Pz()), 1, dummy_track_cov);

                GlobalPoint match_Ds_pos = match_Ds_Vertex.position();
                GlobalError match_Ds_cov = match_Ds_Vertex.positionError();
                GlobalVector match_Ds_PDirec(match_DsFit_Ds_P4.Px(), match_DsFit_Ds_P4.Py(), match_DsFit_Ds_P4.Pz());

                ftree->match_Ds_primvtx_FD.fillAll(match_Ds_pos-primvtx_pos, match_Ds_cov+primvtx_cov, match_Ds_PDirec);

                ftree->match_Kp_primvtx_ip.fillAll(match_DsFit_Tracks[0], primvtx);
                ftree->match_Km_primvtx_ip.fillAll(match_DsFit_Tracks[1], primvtx);
                ftree->match_pi_primvtx_ip.fillAll(match_DsFit_Tracks[2], primvtx);
                ftree->match_phi_primvtx_ip.fillAll((*ttBuilder).build(match_phi_dummy_track), primvtx);
                ftree->match_Ds_primvtx_ip.fillAll( (*ttBuilder).build(match_Ds_dummy_track),  primvtx);

                reco::TrackCollection match_PVtracks_noDs;
                reco::TrackCollection match_PVtracks_withDs;

                for(size_t it=0; it<packedHandle->size(); it++){
                    if(it==static_cast<size_t>(bestKp.index)) continue;
                    if(it==static_cast<size_t>(bestKm.index)) continue;
                    if(it==static_cast<size_t>(bestpi.index)) continue;
                    const auto& pf = (*packedHandle)[it];

                    if(reco::deltaR(match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.3){
                        if( abs(pf.pdgId()) == 211 && (pf.fromPV() > 1) ) ftree->match_Ds_IsoR03.sumChargedHadronPt += pf.pt();
                        else if( pf.pdgId() == 130 ) ftree->match_Ds_IsoR03.sumNeutralHadronPt += pf.pt();
                        else if( pf.pdgId() == 22 ) ftree->match_Ds_IsoR03.sumPhotonPt += pf.pt();

                        if( pf.fromPV() <= 1 ) ftree->match_Ds_IsoR03.sumPUPt += pf.pt();
                    }

                    if(reco::deltaR(match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.4){
                        if( abs(pf.pdgId()) == 211 && (pf.fromPV() > 1) ) ftree->match_Ds_IsoR04.sumChargedHadronPt += pf.pt();
                        else if( pf.pdgId() == 130 ) ftree->match_Ds_IsoR04.sumNeutralHadronPt += pf.pt();
                        else if( pf.pdgId() == 22 ) ftree->match_Ds_IsoR04.sumPhotonPt += pf.pt();

                        if( pf.fromPV() <= 1 ) ftree->match_Ds_IsoR04.sumPUPt += pf.pt();
                    }

                    if( ! pf.hasTrackDetails() ) continue;
                    match_PVtracks_noDs.push_back(pf.pseudoTrack());
                    match_PVtracks_withDs.push_back(pf.pseudoTrack());
                }
                match_PVtracks_withDs.push_back(match_Ds_dummy_track); 

                ftree->match_Ds_IsoR03.PFIso = ( ftree->match_Ds_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->match_Ds_IsoR03.sumNeutralHadronPt + ftree->match_Ds_IsoR03.sumPhotonPt - 0.5*ftree->match_Ds_IsoR03.sumPUPt) ) / match_DsFit_Ds_P4.Pt();
                ftree->match_Ds_IsoR04.PFIso = ( ftree->match_Ds_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->match_Ds_IsoR04.sumNeutralHadronPt + ftree->match_Ds_IsoR04.sumPhotonPt - 0.5*ftree->match_Ds_IsoR04.sumPUPt) ) / match_DsFit_Ds_P4.Pt();

                std::vector<TransientVertex> match_pvs_noDs = revertex->makeVertices(match_PVtracks_noDs, *bsHandle, iSetup, "noBS");
                if( !(match_pvs_noDs.empty()) ){
                    reco::Vertex match_pv_noDs = reco::Vertex(match_pvs_noDs.front());
                    ftree->match_PVnoDs.fillAll( match_pv_noDs );

                    GlobalPoint match_PVnoDs_pos = match_pvs_noDs.front().position();
                    GlobalError match_PVnoDs_cov = match_pvs_noDs.front().positionError();

                    ftree->match_Ds_PVnoDs_FD.fillAll(match_Ds_pos-match_PVnoDs_pos, match_Ds_cov+match_PVnoDs_cov, match_Ds_PDirec);

                    ftree->match_Kp_PVnoDs_ip.fillAll(match_DsFit_Tracks[0], match_pv_noDs);
                    ftree->match_Km_PVnoDs_ip.fillAll(match_DsFit_Tracks[1], match_pv_noDs);
                    ftree->match_pi_PVnoDs_ip.fillAll(match_DsFit_Tracks[2], match_pv_noDs);
                    ftree->match_phi_PVnoDs_ip.fillAll((*ttBuilder).build(match_phi_dummy_track), match_pv_noDs);
                    ftree->match_Ds_PVnoDs_ip.fillAll( (*ttBuilder).build(match_Ds_dummy_track),  match_pv_noDs);
                }

                std::vector<TransientVertex> match_pvs_withDs = revertex->makeVertices(match_PVtracks_withDs, *bsHandle, iSetup, "noBS");
                if( !(match_pvs_withDs.empty()) ){
                    reco::Vertex match_pv_withDs = reco::Vertex(match_pvs_withDs.front());
                    ftree->match_PVwithDs.fillAll( match_pv_withDs );

                    GlobalPoint match_PVwithDs_pos = match_pvs_withDs.front().position();
                    GlobalError match_PVwithDs_cov = match_pvs_withDs.front().positionError();

                    ftree->match_Ds_PVwithDs_FD.fillAll(match_Ds_pos-match_PVwithDs_pos, match_Ds_cov+match_PVwithDs_cov, match_Ds_PDirec);

                    ftree->match_Kp_PVwithDs_ip.fillAll(match_DsFit_Tracks[0], match_pv_withDs);
                    ftree->match_Km_PVwithDs_ip.fillAll(match_DsFit_Tracks[1], match_pv_withDs);
                    ftree->match_pi_PVwithDs_ip.fillAll(match_DsFit_Tracks[2], match_pv_withDs);
                    ftree->match_phi_PVwithDs_ip.fillAll((*ttBuilder).build(match_phi_dummy_track), match_pv_withDs);
                    ftree->match_Ds_PVwithDs_ip.fillAll( (*ttBuilder).build(match_Ds_dummy_track),  match_pv_withDs);
                }

                ftree->FillMatchDs();
            }
        }
    }

    ftree->match_mu.clearAll();
    ftree->match_mu_primvtx_ip.clearAll();
    if( ftree->num_match_mu > 0 ){

        const auto& match_muon = (*muonHandle)[bestmu.index];

        ftree->match_mu.fillAll( match_muon, primvtx );
        ftree->match_mu_IsoR03.sumChargedHadronPt = match_muon.pfIsolationR03().sumChargedHadronPt;
        ftree->match_mu_IsoR03.sumNeutralHadronPt = match_muon.pfIsolationR03().sumNeutralHadronEt;
        ftree->match_mu_IsoR03.sumPhotonPt = match_muon.pfIsolationR03().sumPhotonEt;
        ftree->match_mu_IsoR03.sumPUPt = match_muon.pfIsolationR03().sumPUPt;
        ftree->match_mu_IsoR03.PFIso = (ftree->match_mu_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->match_mu_IsoR03.sumNeutralHadronPt + ftree->match_mu_IsoR03.sumPhotonPt - 0.5*ftree->match_mu_IsoR03.sumPUPt)) / match_muon.pt();
        ftree->match_mu_IsoR04.sumChargedHadronPt = match_muon.pfIsolationR04().sumChargedHadronPt;
        ftree->match_mu_IsoR04.sumNeutralHadronPt = match_muon.pfIsolationR04().sumNeutralHadronEt;
        ftree->match_mu_IsoR04.sumPhotonPt = match_muon.pfIsolationR04().sumPhotonEt;
        ftree->match_mu_IsoR04.sumPUPt = match_muon.pfIsolationR04().sumPUPt;
        ftree->match_mu_IsoR04.PFIso = (ftree->match_mu_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->match_mu_IsoR04.sumNeutralHadronPt + ftree->match_mu_IsoR04.sumPhotonPt - 0.5*ftree->match_mu_IsoR04.sumPUPt)) / match_muon.pt();
        ftree->match_mu_primvtx_ip.fillAll( (*ttBuilder).build(match_muon.muonBestTrack()), primvtx);

        ftree->FillMatchmu();
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
                ftree->Kp_match        = false;
                ftree->Km_match        = false;
                ftree->pi_match        = false;
                ftree->match_entry     = false;
                ftree->non_match_entry = false;

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
                        else if( pf.pdgId() == 130 ) ftree->Ds_IsoR03.sumNeutralHadronPt += pf.pt();
                        else if( pf.pdgId() == 22 ) ftree->Ds_IsoR03.sumPhotonPt += pf.pt();

                        if( pf.fromPV() <= 1 ) ftree->Ds_IsoR03.sumPUPt += pf.pt();
                    }

                    if(reco::deltaR(DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.4){
                        if( abs(pf.pdgId()) == 211 && (pf.fromPV()>1) ) ftree->Ds_IsoR04.sumChargedHadronPt += pf.pt();
                        else if( pf.pdgId() == 130 ) ftree->Ds_IsoR04.sumNeutralHadronPt += pf.pt();
                        else if( pf.pdgId() == 22 ) ftree->Ds_IsoR04.sumPhotonPt += pf.pt();

                        if( pf.fromPV() <= 1 ) ftree->Ds_IsoR04.sumPUPt += pf.pt();
                    }

                    if( ! pf.hasTrackDetails() ) continue;
                    PVtracks_noDs.push_back(pf.pseudoTrack());
                    PVtracks_withDs.push_back(pf.pseudoTrack());
                }
                PVtracks_withDs.push_back(Ds_dummy_track);

                ftree->Ds_IsoR03.PFIso = ( ftree->Ds_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->Ds_IsoR03.sumNeutralHadronPt + ftree->Ds_IsoR03.sumPhotonPt - 0.5*ftree->Ds_IsoR03.sumPUPt) ) / DsFit_Ds_P4.Pt();
                ftree->Ds_IsoR04.PFIso = ( ftree->Ds_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->Ds_IsoR04.sumNeutralHadronPt + ftree->Ds_IsoR04.sumPhotonPt - 0.5*ftree->Ds_IsoR04.sumPUPt) ) / DsFit_Ds_P4.Pt();

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

                if( idx_Kp_vec[i] == bestKp.index ) ftree->Kp_match = true;
                else ftree->Kp_match = false;
                if( idx_Km_vec[j] == bestKm.index ) ftree->Km_match = true;
                else ftree->Km_match = false;
                if( idx_pi_vec[k] == bestpi.index ) ftree->pi_match = true;
                else ftree->pi_match = false;

                if(ftree->Kp_match && ftree->Km_match && ftree->pi_match) ftree->match_entry = true;
                else ftree->match_entry = false;
                if(!ftree->Kp_match && !ftree->Km_match && !ftree->pi_match) ftree->non_match_entry = true;
                else ftree->non_match_entry = false;

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
    if(maxidx_Ds_pt > -1) {
        ftree->FillBestDs(maxidx_Ds_pt);
    }


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

    if( maxidx_mu_pt > -1 ){
        ftree->mu.clearAll();
        ftree->mu_primvtx_ip.clearAll();
        ftree->best_mu_match = false;

        const auto& muon = (*muonHandle)[maxidx_mu_pt];

        ftree->mu.fillAll( muon, primvtx );
        ftree->mu_IsoR03.sumChargedHadronPt = muon.pfIsolationR03().sumChargedHadronPt;
        ftree->mu_IsoR03.sumNeutralHadronPt = muon.pfIsolationR03().sumNeutralHadronEt;
        ftree->mu_IsoR03.sumPhotonPt = muon.pfIsolationR03().sumPhotonEt;
        ftree->mu_IsoR03.sumPUPt = muon.pfIsolationR03().sumPUPt;
        ftree->mu_IsoR03.PFIso = (ftree->mu_IsoR03.sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR03.sumNeutralHadronPt + ftree->mu_IsoR03.sumPhotonPt - 0.5*ftree->mu_IsoR03.sumPUPt)) / muon.pt();
        ftree->mu_IsoR04.sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
        ftree->mu_IsoR04.sumNeutralHadronPt = muon.pfIsolationR04().sumNeutralHadronEt;
        ftree->mu_IsoR04.sumPhotonPt = muon.pfIsolationR04().sumPhotonEt;
        ftree->mu_IsoR04.sumPUPt = muon.pfIsolationR04().sumPUPt;
        ftree->mu_IsoR04.PFIso = (ftree->mu_IsoR04.sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR04.sumNeutralHadronPt + ftree->mu_IsoR04.sumPhotonPt - 0.5*ftree->mu_IsoR04.sumPUPt)) / muon.pt();
        ftree->mu_primvtx_ip.fillAll( (*ttBuilder).build(muon.muonBestTrack()), primvtx);
        
        if( maxidx_mu_pt == bestmu.index ) ftree->best_mu_match = true;
        else ftree->best_mu_match = false;

        ftree->FillBestmu();
    }

    ftree->tree->Fill();

    return;
}

bool PVStudy::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
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

void PVStudy::endJob()
{
    edm::LogInfo("PVStudy") << "Processed all events.";
}

void PVStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PVStudy);

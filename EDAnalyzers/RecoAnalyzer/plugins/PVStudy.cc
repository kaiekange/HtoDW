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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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
        edm::EDGetTokenT<reco::VertexCollection> pvToken;
        edm::EDGetTokenT<reco::BeamSpot> bsToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        KalmanVertexFitter KVFitter;
        VertexReProducer *revertex;

        const edm::Service<TFileService> fs;
        PVTree *ftree;

        std::vector<float> all_dR_Kp;
        std::vector<float> all_dR_Km;
        std::vector<float> all_dR_pi;
        std::vector<float> dR_Kp;
        std::vector<float> dR_Km;
        std::vector<float> dR_pi;

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

        const double Mass_Ds = 1.96835;
        const double Mass_phi = 1.019460;
        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
};

PVStudy::PVStudy(const edm::ParameterSet& iConfig) :
    prunedToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPart"))),
    packedToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
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
    all_dR_Kp.clear();
    all_dR_Km.clear();
    all_dR_pi.clear();
    dR_Kp.clear();
    dR_Km.clear();
    dR_pi.clear();

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

    if(ftree->num_Gen_Kp!=1) return;
    if(ftree->num_Gen_Km!=1) return;
    if(ftree->num_Gen_pi!=1) return;
    if(ftree->num_Gen_phi!=1) return;
    if(ftree->num_Gen_Ds!=1) return;
    if(ftree->num_Gen_mu!=1) return;
    if(ftree->num_Gen_nu!=1) return;
    if(ftree->num_Gen_W!=1) return;
    if(ftree->num_Gen_H!=1) return;

    const auto& Kp_GP = (*prunedHandle)[Gen_Kp_idx[0]];
    const auto& Km_GP = (*prunedHandle)[Gen_Km_idx[0]];
    const auto& pi_GP = (*prunedHandle)[Gen_pi_idx[0]];
    const auto& phi_GP = (*prunedHandle)[Gen_phi_idx[0]];
    const auto& Ds_GP = (*prunedHandle)[Gen_Ds_idx[0]];
    const auto& mu_GP = (*prunedHandle)[Gen_mu_idx[0]];
    const auto& nu_GP = (*prunedHandle)[Gen_nu_idx[0]];
    const auto& W_GP = (*prunedHandle)[Gen_W_idx[0]];
    const auto& H_GP = (*prunedHandle)[Gen_H_idx[0]];

    ftree->Gen_Kp_eta = Kp_GP.eta();
    ftree->Gen_Kp_phi = Kp_GP.phi();
    ftree->Gen_Kp_vx = Kp_GP.vx();
    ftree->Gen_Kp_vy = Kp_GP.vy();
    ftree->Gen_Kp_vz = Kp_GP.vz();
    ftree->Gen_Kp_p = Kp_GP.p();
    ftree->Gen_Kp_pt = Kp_GP.pt();
    ftree->Gen_Kp_px = Kp_GP.px();
    ftree->Gen_Kp_py = Kp_GP.py();
    ftree->Gen_Kp_pz = Kp_GP.pz();
    TVector3 Gen_Kp_P3(Kp_GP.px(), Kp_GP.py(), Kp_GP.pz());

    ftree->Gen_Km_eta = Km_GP.eta();
    ftree->Gen_Km_phi = Km_GP.phi();
    ftree->Gen_Km_vx = Km_GP.vx();
    ftree->Gen_Km_vy = Km_GP.vy();
    ftree->Gen_Km_vz = Km_GP.vz();
    ftree->Gen_Km_p = Km_GP.p();
    ftree->Gen_Km_pt = Km_GP.pt();
    ftree->Gen_Km_px = Km_GP.px();
    ftree->Gen_Km_py = Km_GP.py();
    ftree->Gen_Km_pz = Km_GP.pz();
    TVector3 Gen_Km_P3(Km_GP.px(), Km_GP.py(), Km_GP.pz());

    ftree->Gen_pi_eta = pi_GP.eta();
    ftree->Gen_pi_phi = pi_GP.phi();
    ftree->Gen_pi_vx = pi_GP.vx();
    ftree->Gen_pi_vy = pi_GP.vy();
    ftree->Gen_pi_vz = pi_GP.vz();
    ftree->Gen_pi_p = pi_GP.p();
    ftree->Gen_pi_pt = pi_GP.pt();
    ftree->Gen_pi_px = pi_GP.px();
    ftree->Gen_pi_py = pi_GP.py();
    ftree->Gen_pi_pz = pi_GP.pz();
    TVector3 Gen_pi_P3(pi_GP.px(), pi_GP.py(), pi_GP.pz());

    ftree->Gen_phi_eta = phi_GP.eta();
    ftree->Gen_phi_phi = phi_GP.phi();
    ftree->Gen_phi_vx = phi_GP.vx();
    ftree->Gen_phi_vy = phi_GP.vy();
    ftree->Gen_phi_vz = phi_GP.vz();
    ftree->Gen_phi_p = phi_GP.p();
    ftree->Gen_phi_pt = phi_GP.pt();
    ftree->Gen_phi_px = phi_GP.px();
    ftree->Gen_phi_py = phi_GP.py();
    ftree->Gen_phi_pz = phi_GP.pz();
    TVector3 Gen_phi_P3(phi_GP.px(), phi_GP.py(), phi_GP.pz());

    ftree->Gen_Ds_eta = Ds_GP.eta();
    ftree->Gen_Ds_phi = Ds_GP.phi();
    ftree->Gen_Ds_vx = Ds_GP.vx();
    ftree->Gen_Ds_vy = Ds_GP.vy();
    ftree->Gen_Ds_vz = Ds_GP.vz();
    ftree->Gen_Ds_p = Ds_GP.p();
    ftree->Gen_Ds_pt = Ds_GP.pt();
    ftree->Gen_Ds_px = Ds_GP.px();
    ftree->Gen_Ds_py = Ds_GP.py();
    ftree->Gen_Ds_pz = Ds_GP.pz();
    TVector3 Gen_Ds_P3(Ds_GP.px(), Ds_GP.py(), Ds_GP.pz());

    ftree->Gen_mu_eta = mu_GP.eta();
    ftree->Gen_mu_phi = mu_GP.phi();
    ftree->Gen_mu_vx = mu_GP.vx();
    ftree->Gen_mu_vy = mu_GP.vy();
    ftree->Gen_mu_vz = mu_GP.vz();
    ftree->Gen_mu_p = mu_GP.p();
    ftree->Gen_mu_pt = mu_GP.pt();
    ftree->Gen_mu_px = mu_GP.px();
    ftree->Gen_mu_py = mu_GP.py();
    ftree->Gen_mu_pz = mu_GP.pz();
    TVector3 Gen_mu_P3(mu_GP.px(), mu_GP.py(), mu_GP.pz());

    ftree->Gen_nu_eta = nu_GP.eta();
    ftree->Gen_nu_phi = nu_GP.phi();
    ftree->Gen_nu_vx = nu_GP.vx();
    ftree->Gen_nu_vy = nu_GP.vy();
    ftree->Gen_nu_vz = nu_GP.vz();
    ftree->Gen_nu_p = nu_GP.p();
    ftree->Gen_nu_pt = nu_GP.pt();
    ftree->Gen_nu_px = nu_GP.px();
    ftree->Gen_nu_py = nu_GP.py();
    ftree->Gen_nu_pz = nu_GP.pz();
    TVector3 Gen_nu_P3(nu_GP.px(), nu_GP.py(), nu_GP.pz());

    ftree->Gen_W_eta = W_GP.eta();
    ftree->Gen_W_phi = W_GP.phi();
    ftree->Gen_W_vx = W_GP.vx();
    ftree->Gen_W_vy = W_GP.vy();
    ftree->Gen_W_vz = W_GP.vz();
    ftree->Gen_W_p = W_GP.p();
    ftree->Gen_W_pt = W_GP.pt();
    ftree->Gen_W_px = W_GP.px();
    ftree->Gen_W_py = W_GP.py();
    ftree->Gen_W_pz = W_GP.pz();
    TVector3 Gen_W_P3(W_GP.px(), W_GP.py(), W_GP.pz());

    ftree->Gen_H_eta = H_GP.eta();
    ftree->Gen_H_phi = H_GP.phi();
    ftree->Gen_H_vx = H_GP.vx();
    ftree->Gen_H_vy = H_GP.vy();
    ftree->Gen_H_vz = H_GP.vz();
    ftree->Gen_H_p = H_GP.p();
    ftree->Gen_H_pt = H_GP.pt();
    ftree->Gen_H_px = H_GP.px();
    ftree->Gen_H_py = H_GP.py();
    ftree->Gen_H_pz = H_GP.pz();
    TVector3 Gen_H_P3(H_GP.px(), H_GP.py(), H_GP.pz());

    ftree->Gen_Kp_pp = Gen_Kp_P3.Pt(Gen_phi_P3);
    ftree->Gen_Kp_pl = Gen_Kp_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag();
    ftree->Gen_Km_pp = Gen_Km_P3.Pt(Gen_phi_P3);
    ftree->Gen_Km_pl = Gen_Km_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag();
    ftree->Gen_pi_pp = Gen_pi_P3.Pt(Gen_Ds_P3);
    ftree->Gen_pi_pl = Gen_pi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag();
    ftree->Gen_phi_pp = Gen_phi_P3.Pt(Gen_Ds_P3);
    ftree->Gen_phi_pl = Gen_phi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag();

    ftree->Gen_dR_Kp_Km = reco::deltaR(Kp_GP, Km_GP);
    ftree->Gen_dR_Kp_phi = reco::deltaR(Kp_GP, phi_GP);
    ftree->Gen_dR_Km_phi = reco::deltaR(Km_GP, phi_GP);
    ftree->Gen_dR_Kp_pi = reco::deltaR(Kp_GP, pi_GP);
    ftree->Gen_dR_Km_pi = reco::deltaR(Km_GP, pi_GP);
    ftree->Gen_dR_pi_phi = reco::deltaR(phi_GP, pi_GP);
    ftree->Gen_dR_Kp_Ds = reco::deltaR(Kp_GP, Ds_GP);
    ftree->Gen_dR_Km_Ds = reco::deltaR(Km_GP, Ds_GP);
    ftree->Gen_dR_phi_Ds = reco::deltaR(phi_GP, Ds_GP);
    ftree->Gen_dR_pi_Ds = reco::deltaR(pi_GP, Ds_GP);

    ftree->Gen_dxy_Kp_Km = sqrt(pow(ftree->Gen_Kp_vx-ftree->Gen_Km_vx,2) + pow(ftree->Gen_Kp_vy-ftree->Gen_Km_vy,2));
    ftree->Gen_dxy_Kp_pi = sqrt(pow(ftree->Gen_Kp_vx-ftree->Gen_pi_vx,2) + pow(ftree->Gen_Kp_vy-ftree->Gen_pi_vy,2));
    ftree->Gen_dxy_Km_pi = sqrt(pow(ftree->Gen_Km_vx-ftree->Gen_pi_vx,2) + pow(ftree->Gen_Km_vy-ftree->Gen_pi_vy,2));
    ftree->Gen_dz_Kp_Km = abs(ftree->Gen_Kp_vz-ftree->Gen_Km_vz);
    ftree->Gen_dz_Kp_pi = abs(ftree->Gen_Kp_vz-ftree->Gen_pi_vz);
    ftree->Gen_dz_Km_pi = abs(ftree->Gen_Km_vz-ftree->Gen_pi_vz);

    ftree->BS_Reset();
    ftree->BS_type = bsHandle->type();
    ftree->BS_x0 = bsHandle->x0();
    ftree->BS_y0 = bsHandle->y0();
    ftree->BS_z0 = bsHandle->z0();
    ftree->BS_sigmaZ = bsHandle->sigmaZ();
    ftree->BS_dxdz = bsHandle->dxdz();
    ftree->BS_dydz = bsHandle->dydz();
    ftree->BS_BWX = bsHandle->BeamWidthX();
    ftree->BS_BWY = bsHandle->BeamWidthY();
    ftree->BS_x0err = bsHandle->x0Error();
    ftree->BS_y0err = bsHandle->y0Error();
    ftree->BS_z0err = bsHandle->z0Error();
    ftree->BS_sigmaZ0err = bsHandle->sigmaZ0Error();
    ftree->BS_dxdzerr = bsHandle->dxdzError();
    ftree->BS_dydzerr = bsHandle->dydzError();
    ftree->BS_BWXerr = bsHandle->BeamWidthXError();
    ftree->BS_BWYerr = bsHandle->BeamWidthYError();
    ftree->BS_emitX = bsHandle->emittanceX();
    ftree->BS_emitY = bsHandle->emittanceY();
    ftree->BS_betaStar = bsHandle->betaStar();

    ftree->PV_Reset();
    const reco::Vertex& primaryvertex = pvHandle->front();
    ftree->PV_vx = primaryvertex.x();
    ftree->PV_vy = primaryvertex.y();
    ftree->PV_vz = primaryvertex.z();
    ftree->PV_vxerr = primaryvertex.xError();
    ftree->PV_vyerr = primaryvertex.yError();
    ftree->PV_vzerr = primaryvertex.zError();
    GlobalPoint primaryvertex_pos(primaryvertex.x(), primaryvertex.y(), primaryvertex.z());
    GlobalError primaryvertex_cov(
        primaryvertex.covariance(0,0), primaryvertex.covariance(0,1), primaryvertex.covariance(0,2),
        primaryvertex.covariance(1,1), primaryvertex.covariance(1,2),
        primaryvertex.covariance(2,2)
    );

    ftree->Match_Reset();

    struct MatchInfo {int index; float dr;};
    std::vector<MatchInfo> MatchKp, MatchKm, Matchpi;

    for(size_t i=0; i<packedHandle->size(); i++){
        const auto& pf = (*packedHandle)[i];

        const float dr_Kp = reco::deltaR(pf.eta(), pf.phi(), Kp_GP.eta(), Kp_GP.phi());
        const float dr_Km = reco::deltaR(pf.eta(), pf.phi(), Km_GP.eta(), Km_GP.phi());
        const float dr_pi = reco::deltaR(pf.eta(), pf.phi(), pi_GP.eta(), pi_GP.phi());

        all_dR_Kp.push_back(dr_Kp);
        all_dR_Km.push_back(dr_Km);
        all_dR_pi.push_back(dr_pi);

        if(pf.pdgId() == 211 && pf.trackHighPurity()){
            dR_Kp.push_back(dr_Kp);
            dR_pi.push_back(dr_pi);
            MatchKp.push_back({static_cast<int>(i), dr_Kp});
            Matchpi.push_back({static_cast<int>(i), dr_pi});
        } else if(pf.pdgId() == -211 && pf.trackHighPurity()){
            dR_Km.push_back(dr_Km);
            MatchKm.push_back({static_cast<int>(i), dr_Km});
        }
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

    if (bestKp.index == bestpi.index && bestKp.index != -1){
        if (bestKp.dr <= bestpi.dr) {
            bestpi = findBestMatchExcluding(Matchpi, bestKp.index);
        } else{
            const auto altKp = findBestMatchExcluding(MatchKp, bestpi.index);
            if (altKp.index != -1) {
                bestKp = altKp;
            } else {
                bestpi = findBestMatchExcluding(Matchpi, bestKp.index);
            }
        }
    }

    if(bestKp.dr < 0.03){
        ftree->num_match_Kp++;
        ftree->match_Kp_idx = bestKp.index;
    }
    if(bestKm.dr < 0.03){
        ftree->num_match_Km++;
        ftree->match_Km_idx = bestKm.index;
    }
    if(bestpi.dr < 0.03){
        ftree->num_match_pi++;
        ftree->match_pi_idx = bestpi.index;
    }

    if( bestKp.dr < 0.03 && bestKm.dr < 0.03 && bestpi.dr < 0.03 && (*packedHandle)[bestKp.index].hasTrackDetails() && (*packedHandle)[bestKm.index].hasTrackDetails() && (*packedHandle)[bestpi.index].hasTrackDetails() && (*packedHandle)[bestKp.index].trackHighPurity() && (*packedHandle)[bestKm.index].trackHighPurity() && (*packedHandle)[bestpi.index].trackHighPurity() ){

        const auto& match_Kp_PF = (*packedHandle)[bestKp.index];
        const auto& match_Km_PF = (*packedHandle)[bestKm.index];
        const auto& match_pi_PF = (*packedHandle)[bestpi.index];

        // Original Kp
        ftree->match_Kp_isIsolatedChargedHadron = match_Kp_PF.isIsolatedChargedHadron();
        ftree->match_Kp_charge = match_Kp_PF.charge();
        ftree->match_Kp_eta = match_Kp_PF.eta();
        ftree->match_Kp_phi = match_Kp_PF.phi();
        ftree->match_Kp_vx = match_Kp_PF.vx();
        ftree->match_Kp_vy = match_Kp_PF.vy();
        ftree->match_Kp_vz = match_Kp_PF.vz();
        ftree->match_Kp_p = match_Kp_PF.p();
        ftree->match_Kp_pt = match_Kp_PF.pt();
        ftree->match_Kp_px = match_Kp_PF.px();
        ftree->match_Kp_py = match_Kp_PF.py();
        ftree->match_Kp_pz = match_Kp_PF.pz();
        TLorentzVector match_Kp_P4;
        match_Kp_P4.SetXYZM(match_Kp_PF.px(), match_Kp_PF.py(), match_Kp_PF.pz(), Mass_K);

        // Original Km
        ftree->match_Km_isIsolatedChargedHadron = match_Km_PF.isIsolatedChargedHadron();
        ftree->match_Km_charge = match_Km_PF.charge();
        ftree->match_Km_eta = match_Km_PF.eta();
        ftree->match_Km_phi = match_Km_PF.phi();
        ftree->match_Km_vx = match_Km_PF.vx();
        ftree->match_Km_vy = match_Km_PF.vy();
        ftree->match_Km_vz = match_Km_PF.vz();
        ftree->match_Km_p = match_Km_PF.p();
        ftree->match_Km_pt = match_Km_PF.pt();
        ftree->match_Km_px = match_Km_PF.px();
        ftree->match_Km_py = match_Km_PF.py();
        ftree->match_Km_pz = match_Km_PF.pz();
        TLorentzVector match_Km_P4;
        match_Km_P4.SetXYZM(match_Km_PF.px(), match_Km_PF.py(), match_Km_PF.pz(), Mass_K);

        // Original pi
        ftree->match_pi_isIsolatedChargedHadron = match_pi_PF.isIsolatedChargedHadron();
        ftree->match_pi_charge = match_pi_PF.charge();
        ftree->match_pi_eta = match_pi_PF.eta();
        ftree->match_pi_phi = match_pi_PF.phi();
        ftree->match_pi_vx = match_pi_PF.vx();
        ftree->match_pi_vy = match_pi_PF.vy();
        ftree->match_pi_vz = match_pi_PF.vz();
        ftree->match_pi_p = match_pi_PF.p();
        ftree->match_pi_pt = match_pi_PF.pt();
        ftree->match_pi_px = match_pi_PF.px();
        ftree->match_pi_py = match_pi_PF.py();
        ftree->match_pi_pz = match_pi_PF.pz();
        TLorentzVector match_pi_P4;
        match_pi_P4.SetXYZM(match_pi_PF.px(), match_pi_PF.py(), match_pi_PF.pz(), Mass_pi);

        // Original phi
        TLorentzVector match_phi_P4 = match_Kp_P4 + match_Km_P4;
        ftree->match_phi_eta = match_phi_P4.Eta();
        ftree->match_phi_phi = match_phi_P4.Phi();
        ftree->match_phi_p = match_phi_P4.P();
        ftree->match_phi_pt = match_phi_P4.Pt();
        ftree->match_phi_px = match_phi_P4.Px();
        ftree->match_phi_py = match_phi_P4.Py();
        ftree->match_phi_pz = match_phi_P4.Pz();
        ftree->match_phi_invm = match_phi_P4.M();

        // Original Ds
        TLorentzVector match_Ds_P4 = match_Kp_P4 + match_Km_P4 + match_pi_P4;
        ftree->match_Ds_eta = match_Ds_P4.Eta();
        ftree->match_Ds_phi = match_Ds_P4.Phi();
        ftree->match_Ds_p = match_Ds_P4.P();
        ftree->match_Ds_pt = match_Ds_P4.Pt();
        ftree->match_Ds_px = match_Ds_P4.Px();
        ftree->match_Ds_py = match_Ds_P4.Py();
        ftree->match_Ds_pz = match_Ds_P4.Pz();
        ftree->match_Ds_invm = match_Ds_P4.M();

        // some dxy dz
        ftree->match_dxy_Kp_Km = sqrt(pow(ftree->match_Kp_vx-ftree->match_Km_vx,2) + pow(ftree->match_Kp_vy-ftree->match_Km_vy,2));
        ftree->match_dxy_Kp_pi = sqrt(pow(ftree->match_Kp_vx-ftree->match_pi_vx,2) + pow(ftree->match_Kp_vy-ftree->match_pi_vy,2));
        ftree->match_dxy_Km_pi = sqrt(pow(ftree->match_Km_vx-ftree->match_pi_vx,2) + pow(ftree->match_Km_vy-ftree->match_pi_vy,2));

        ftree->match_dz_Kp_Km = abs(ftree->match_Kp_vz-ftree->match_Km_vz);
        ftree->match_dz_Kp_pi = abs(ftree->match_Kp_vz-ftree->match_pi_vz);
        ftree->match_dz_Km_pi = abs(ftree->match_Km_vz-ftree->match_pi_vz);

        // Original PP PL
        ftree->match_Kp_pp = match_Kp_P4.Vect().Pt(match_phi_P4.Vect());
        ftree->match_Kp_pl = match_Kp_P4.Vect().Dot(match_phi_P4.Vect())/match_phi_P4.P();
        ftree->match_Km_pp = match_Km_P4.Vect().Pt(match_phi_P4.Vect());
        ftree->match_Km_pl = match_Km_P4.Vect().Dot(match_phi_P4.Vect())/match_phi_P4.P();

        ftree->match_phi_pp = match_phi_P4.Vect().Pt(match_Ds_P4.Vect());
        ftree->match_phi_pl = match_phi_P4.Vect().Dot(match_Ds_P4.Vect())/match_Ds_P4.P();
        ftree->match_pi_pp = match_pi_P4.Vect().Pt(match_Ds_P4.Vect());
        ftree->match_pi_pl = match_pi_P4.Vect().Dot(match_Ds_P4.Vect())/match_Ds_P4.P();

        // Original dR
        ftree->match_dR_Kp_Km = reco::deltaR(match_Kp_P4.Eta(), match_Kp_P4.Phi(), match_Km_P4.Eta(), match_Km_P4.Phi());
        ftree->match_dR_Kp_phi = reco::deltaR(match_Kp_P4.Eta(), match_Kp_P4.Phi(), match_phi_P4.Eta(), match_phi_P4.Phi());
        ftree->match_dR_Km_phi = reco::deltaR(match_Km_P4.Eta(), match_Km_P4.Phi(), match_phi_P4.Eta(), match_phi_P4.Phi());
        ftree->match_dR_Kp_pi = reco::deltaR(match_Kp_P4.Eta(), match_Kp_P4.Phi(), match_pi_P4.Eta(), match_pi_P4.Phi());
        ftree->match_dR_Km_pi = reco::deltaR(match_Km_P4.Eta(), match_Km_P4.Phi(), match_pi_P4.Eta(), match_pi_P4.Phi());
        ftree->match_dR_pi_phi = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_pi_P4.Eta(), match_pi_P4.Phi());
        ftree->match_dR_Kp_Ds = reco::deltaR(match_Kp_P4.Eta(), match_Kp_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
        ftree->match_dR_Km_Ds = reco::deltaR(match_Km_P4.Eta(), match_Km_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
        ftree->match_dR_phi_Ds = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
        ftree->match_dR_pi_Ds = reco::deltaR(match_pi_P4.Eta(), match_pi_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());

        // phi fit
        std::vector<reco::TransientTrack> match_phi_Tracks = {
            (*ttBuilder).build(match_Kp_PF.pseudoTrack()),
            (*ttBuilder).build(match_Km_PF.pseudoTrack())
        };

        TransientVertex match_phi_Vertex = KVFitter.vertex(match_phi_Tracks);

        if( match_phi_Vertex.isValid() && match_phi_Vertex.hasRefittedTracks() ){

            // phi fit variables 
            ftree->match_phiFit_chi2 = match_phi_Vertex.totalChiSquared();
            ftree->match_phiFit_ndof = match_phi_Vertex.degreesOfFreedom();
            ftree->match_phiFit_chi2ndof = match_phi_Vertex.normalisedChiSquared();
            ftree->match_phiFit_vx = match_phi_Vertex.position().x();
            ftree->match_phiFit_vy = match_phi_Vertex.position().y();
            ftree->match_phiFit_vz = match_phi_Vertex.position().z();
            ftree->match_phiFit_vxerr = std::sqrt(match_phi_Vertex.positionError().cxx());
            ftree->match_phiFit_vyerr = std::sqrt(match_phi_Vertex.positionError().cyy());
            ftree->match_phiFit_vzerr = std::sqrt(match_phi_Vertex.positionError().czz());

            // some dxy dz
            ftree->match_dxy_Kp_phi = sqrt(pow(ftree->match_Kp_vx-ftree->match_phiFit_vx,2) + pow(ftree->match_Kp_vy-ftree->match_phiFit_vy,2));
            ftree->match_dxy_Km_phi = sqrt(pow(ftree->match_Km_vx-ftree->match_phiFit_vx,2) + pow(ftree->match_Km_vy-ftree->match_phiFit_vy,2));
            ftree->match_dxy_pi_phi = sqrt(pow(ftree->match_pi_vx-ftree->match_phiFit_vx,2) + pow(ftree->match_pi_vy-ftree->match_phiFit_vy,2));

            ftree->match_dz_Kp_phi = abs(ftree->match_Kp_vz-ftree->match_phiFit_vz);
            ftree->match_dz_Km_phi = abs(ftree->match_Km_vz-ftree->match_phiFit_vz);
            ftree->match_dz_pi_phi = abs(ftree->match_pi_vz-ftree->match_phiFit_vz);

            std::vector<reco::TransientTrack> match_phiFit_Tracks = match_phi_Vertex.refittedTracks();

            // phi fit Kp
            TLorentzVector match_phiFit_Kp_P4;
            match_phiFit_Kp_P4.SetXYZM(match_phiFit_Tracks[0].track().px(), match_phiFit_Tracks[0].track().py(), match_phiFit_Tracks[0].track().pz(), Mass_K);
            ftree->match_phiFit_Kp_eta = match_phiFit_Kp_P4.Eta();
            ftree->match_phiFit_Kp_phi = match_phiFit_Kp_P4.Phi();
            ftree->match_phiFit_Kp_p = match_phiFit_Kp_P4.P();
            ftree->match_phiFit_Kp_pt = match_phiFit_Kp_P4.Pt();
            ftree->match_phiFit_Kp_px = match_phiFit_Kp_P4.Px();
            ftree->match_phiFit_Kp_py = match_phiFit_Kp_P4.Py();
            ftree->match_phiFit_Kp_pz = match_phiFit_Kp_P4.Pz();

            // phi fit Km
            TLorentzVector match_phiFit_Km_P4;
            match_phiFit_Km_P4.SetXYZM(match_phiFit_Tracks[1].track().px(), match_phiFit_Tracks[1].track().py(), match_phiFit_Tracks[1].track().pz(), Mass_K);
            ftree->match_phiFit_Km_eta = match_phiFit_Km_P4.Eta();
            ftree->match_phiFit_Km_phi = match_phiFit_Km_P4.Phi();
            ftree->match_phiFit_Km_p = match_phiFit_Km_P4.P();
            ftree->match_phiFit_Km_pt = match_phiFit_Km_P4.Pt();
            ftree->match_phiFit_Km_px = match_phiFit_Km_P4.Px();
            ftree->match_phiFit_Km_py = match_phiFit_Km_P4.Py();
            ftree->match_phiFit_Km_pz = match_phiFit_Km_P4.Pz();

            // phi fit pi = orignal pi
            TLorentzVector match_phiFit_pi_P4 = match_pi_P4;
            ftree->match_phiFit_pi_eta = match_phiFit_pi_P4.Eta();
            ftree->match_phiFit_pi_phi = match_phiFit_pi_P4.Phi();
            ftree->match_phiFit_pi_p = match_phiFit_pi_P4.P();
            ftree->match_phiFit_pi_pt = match_phiFit_pi_P4.Pt();
            ftree->match_phiFit_pi_px = match_phiFit_pi_P4.Px();
            ftree->match_phiFit_pi_py = match_phiFit_pi_P4.Py();
            ftree->match_phiFit_pi_pz = match_phiFit_pi_P4.Pz();

            // phi fit phi
            TLorentzVector match_phiFit_phi_P4 = match_phiFit_Kp_P4 + match_phiFit_Km_P4;
            ftree->match_phiFit_phi_eta = match_phiFit_phi_P4.Eta();
            ftree->match_phiFit_phi_phi = match_phiFit_phi_P4.Phi();
            ftree->match_phiFit_phi_p = match_phiFit_phi_P4.P();
            ftree->match_phiFit_phi_pt = match_phiFit_phi_P4.Pt();
            ftree->match_phiFit_phi_px = match_phiFit_phi_P4.Px();
            ftree->match_phiFit_phi_py = match_phiFit_phi_P4.Py();
            ftree->match_phiFit_phi_pz = match_phiFit_phi_P4.Pz();
            ftree->match_phiFit_phi_invm = match_phiFit_phi_P4.M();

            // phi fit Ds
            TLorentzVector match_phiFit_Ds_P4 = match_phiFit_Kp_P4 + match_phiFit_Km_P4 + match_phiFit_pi_P4;
            ftree->match_phiFit_Ds_eta = match_phiFit_Ds_P4.Eta();
            ftree->match_phiFit_Ds_phi = match_phiFit_Ds_P4.Phi();
            ftree->match_phiFit_Ds_p = match_phiFit_Ds_P4.P();
            ftree->match_phiFit_Ds_pt = match_phiFit_Ds_P4.Pt();
            ftree->match_phiFit_Ds_px = match_phiFit_Ds_P4.Px();
            ftree->match_phiFit_Ds_py = match_phiFit_Ds_P4.Py();
            ftree->match_phiFit_Ds_pz = match_phiFit_Ds_P4.Pz();
            ftree->match_phiFit_Ds_invm = match_phiFit_Ds_P4.M();

            // phi fit PP PL 
            ftree->match_phiFit_Kp_pp = match_phiFit_Kp_P4.Vect().Pt(match_phiFit_phi_P4.Vect());
            ftree->match_phiFit_Kp_pl = match_phiFit_Kp_P4.Vect().Dot(match_phiFit_phi_P4.Vect())/match_phiFit_phi_P4.P();
            ftree->match_phiFit_Km_pp = match_phiFit_Km_P4.Vect().Pt(match_phiFit_phi_P4.Vect());
            ftree->match_phiFit_Km_pl = match_phiFit_Km_P4.Vect().Dot(match_phiFit_phi_P4.Vect())/match_phiFit_phi_P4.P();

            ftree->match_phiFit_phi_pp = match_phiFit_phi_P4.Vect().Pt(match_phiFit_Ds_P4.Vect());
            ftree->match_phiFit_phi_pl = match_phiFit_phi_P4.Vect().Dot(match_phiFit_Ds_P4.Vect())/match_phiFit_Ds_P4.P();
            ftree->match_phiFit_pi_pp = match_phiFit_pi_P4.Vect().Pt(match_phiFit_Ds_P4.Vect());
            ftree->match_phiFit_pi_pl = match_phiFit_pi_P4.Vect().Dot(match_phiFit_Ds_P4.Vect())/match_phiFit_Ds_P4.P();

            // phi fit dR
            ftree->match_phiFit_dR_Kp_Km = reco::deltaR(match_phiFit_Kp_P4.Eta(), match_phiFit_Kp_P4.Phi(), match_phiFit_Km_P4.Eta(), match_phiFit_Km_P4.Phi());
            ftree->match_phiFit_dR_Kp_phi = reco::deltaR(match_phiFit_Kp_P4.Eta(), match_phiFit_Kp_P4.Phi(), match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi());
            ftree->match_phiFit_dR_Km_phi = reco::deltaR(match_phiFit_Km_P4.Eta(), match_phiFit_Km_P4.Phi(), match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi());
            ftree->match_phiFit_dR_Kp_pi = reco::deltaR(match_phiFit_Kp_P4.Eta(), match_phiFit_Kp_P4.Phi(), match_phiFit_pi_P4.Eta(), match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_Km_pi = reco::deltaR(match_phiFit_Km_P4.Eta(), match_phiFit_Km_P4.Phi(), match_phiFit_pi_P4.Eta(), match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_pi_phi = reco::deltaR(match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi(), match_phiFit_pi_P4.Eta(), match_phiFit_pi_P4.Phi());
            ftree->match_phiFit_dR_Kp_Ds = reco::deltaR(match_phiFit_Kp_P4.Eta(), match_phiFit_Kp_P4.Phi(), match_phiFit_Ds_P4.Eta(), match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_Km_Ds = reco::deltaR(match_phiFit_Km_P4.Eta(), match_phiFit_Km_P4.Phi(), match_phiFit_Ds_P4.Eta(), match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_phi_Ds = reco::deltaR(match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi(), match_phiFit_Ds_P4.Eta(), match_phiFit_Ds_P4.Phi());
            ftree->match_phiFit_dR_pi_Ds = reco::deltaR(match_phiFit_pi_P4.Eta(), match_phiFit_pi_P4.Phi(), match_phiFit_Ds_P4.Eta(), match_phiFit_Ds_P4.Phi());

            ftree->match_alpha_phi = (ftree->match_phiFit_Kp_pl - ftree->match_phiFit_Km_pl) / (ftree->match_phiFit_Kp_pl + ftree->match_phiFit_Km_pl);
            ftree->match_beta_phi = std::sqrt(pow(ftree->match_phiFit_phi_p,2) / (pow(ftree->match_phiFit_phi_invm,2) + pow(ftree->match_phiFit_phi_p,2)));
            ftree->match_APvar_phi = pow(ftree->match_phiFit_Kp_pp,2) + pow(ftree->match_alpha_phi*ftree->match_beta_phi*ftree->match_phiFit_phi_invm,2)/4;

            // Ds fit
            std::vector<reco::TransientTrack> match_Ds_Tracks = {
                match_phiFit_Tracks[0],
                match_phiFit_Tracks[1],
                (*ttBuilder).build(match_pi_PF.pseudoTrack())
            };

            TransientVertex match_Ds_Vertex = KVFitter.vertex(match_Ds_Tracks);

            if( match_Ds_Vertex.isValid() && match_Ds_Vertex.hasRefittedTracks() ){

                // Ds fit variables 
                ftree->match_DsFit_chi2 = match_Ds_Vertex.totalChiSquared();
                ftree->match_DsFit_ndof = match_Ds_Vertex.degreesOfFreedom();
                ftree->match_DsFit_chi2ndof = match_Ds_Vertex.normalisedChiSquared();
                ftree->match_DsFit_vx = match_Ds_Vertex.position().x();
                ftree->match_DsFit_vy = match_Ds_Vertex.position().y();
                ftree->match_DsFit_vz = match_Ds_Vertex.position().z();
                ftree->match_DsFit_vxerr = std::sqrt(match_Ds_Vertex.positionError().cxx());
                ftree->match_DsFit_vyerr = std::sqrt(match_Ds_Vertex.positionError().cyy());
                ftree->match_DsFit_vzerr = std::sqrt(match_Ds_Vertex.positionError().czz());

                // some dxy dz
                ftree->match_dxy_Kp_Ds = sqrt(pow(ftree->match_Kp_vx-ftree->match_DsFit_vx,2) + pow(ftree->match_Kp_vy-ftree->match_DsFit_vy,2));
                ftree->match_dxy_Km_Ds = sqrt(pow(ftree->match_Km_vx-ftree->match_DsFit_vx,2) + pow(ftree->match_Km_vy-ftree->match_DsFit_vy,2));
                ftree->match_dxy_pi_Ds = sqrt(pow(ftree->match_pi_vx-ftree->match_DsFit_vx,2) + pow(ftree->match_pi_vy-ftree->match_DsFit_vy,2));
                ftree->match_dxy_phi_Ds = sqrt(pow(ftree->match_phiFit_vx-ftree->match_DsFit_vx,2) + pow(ftree->match_phiFit_vy-ftree->match_DsFit_vy,2));

                ftree->match_dz_Kp_Ds = abs(ftree->match_Kp_vz-ftree->match_DsFit_vz);
                ftree->match_dz_Km_Ds = abs(ftree->match_Km_vz-ftree->match_DsFit_vz);
                ftree->match_dz_pi_Ds = abs(ftree->match_pi_vz-ftree->match_DsFit_vz);
                ftree->match_dz_phi_Ds = abs(ftree->match_phiFit_vz-ftree->match_DsFit_vz);

                std::vector<reco::TransientTrack> match_DsFit_Tracks = match_Ds_Vertex.refittedTracks();

                // Ds fit Kp
                TLorentzVector match_DsFit_Kp_P4;
                match_DsFit_Kp_P4.SetXYZM(match_DsFit_Tracks[0].track().px(), match_DsFit_Tracks[0].track().py(), match_DsFit_Tracks[0].track().pz(), Mass_K);
                ftree->match_DsFit_Kp_eta = match_DsFit_Kp_P4.Eta();
                ftree->match_DsFit_Kp_phi = match_DsFit_Kp_P4.Phi();
                ftree->match_DsFit_Kp_p = match_DsFit_Kp_P4.P();
                ftree->match_DsFit_Kp_pt = match_DsFit_Kp_P4.Pt();
                ftree->match_DsFit_Kp_px = match_DsFit_Kp_P4.Px();
                ftree->match_DsFit_Kp_py = match_DsFit_Kp_P4.Py();
                ftree->match_DsFit_Kp_pz = match_DsFit_Kp_P4.Pz();

                // Ds fit Km
                TLorentzVector match_DsFit_Km_P4;
                match_DsFit_Km_P4.SetXYZM(match_DsFit_Tracks[1].track().px(), match_DsFit_Tracks[1].track().py(), match_DsFit_Tracks[1].track().pz(), Mass_K);
                ftree->match_DsFit_Km_eta = match_DsFit_Km_P4.Eta();
                ftree->match_DsFit_Km_phi = match_DsFit_Km_P4.Phi();
                ftree->match_DsFit_Km_p = match_DsFit_Km_P4.P();
                ftree->match_DsFit_Km_pt = match_DsFit_Km_P4.Pt();
                ftree->match_DsFit_Km_px = match_DsFit_Km_P4.Px();
                ftree->match_DsFit_Km_py = match_DsFit_Km_P4.Py();
                ftree->match_DsFit_Km_pz = match_DsFit_Km_P4.Pz();

                // Ds fit pi
                TLorentzVector match_DsFit_pi_P4;
                match_DsFit_pi_P4.SetXYZM(match_DsFit_Tracks[2].track().px(), match_DsFit_Tracks[2].track().py(), match_DsFit_Tracks[2].track().pz(), Mass_pi);
                ftree->match_DsFit_pi_eta = match_DsFit_pi_P4.Eta();
                ftree->match_DsFit_pi_phi = match_DsFit_pi_P4.Phi();
                ftree->match_DsFit_pi_p = match_DsFit_pi_P4.P();
                ftree->match_DsFit_pi_pt = match_DsFit_pi_P4.Pt();
                ftree->match_DsFit_pi_px = match_DsFit_pi_P4.Px();
                ftree->match_DsFit_pi_py = match_DsFit_pi_P4.Py();
                ftree->match_DsFit_pi_pz = match_DsFit_pi_P4.Pz();

                // Ds fit phi
                TLorentzVector match_DsFit_phi_P4 = match_DsFit_Kp_P4 + match_DsFit_Km_P4;
                ftree->match_DsFit_phi_eta = match_DsFit_phi_P4.Eta();
                ftree->match_DsFit_phi_phi = match_DsFit_phi_P4.Phi();
                ftree->match_DsFit_phi_p = match_DsFit_phi_P4.P();
                ftree->match_DsFit_phi_pt = match_DsFit_phi_P4.Pt();
                ftree->match_DsFit_phi_px = match_DsFit_phi_P4.Px();
                ftree->match_DsFit_phi_py = match_DsFit_phi_P4.Py();
                ftree->match_DsFit_phi_pz = match_DsFit_phi_P4.Pz();
                ftree->match_DsFit_phi_invm = match_DsFit_phi_P4.M();

                // Ds fit Ds
                TLorentzVector match_DsFit_Ds_P4 = match_DsFit_Kp_P4 + match_DsFit_Km_P4 + match_DsFit_pi_P4;
                ftree->match_DsFit_Ds_eta = match_DsFit_Ds_P4.Eta();
                ftree->match_DsFit_Ds_phi = match_DsFit_Ds_P4.Phi();
                ftree->match_DsFit_Ds_p = match_DsFit_Ds_P4.P();
                ftree->match_DsFit_Ds_pt = match_DsFit_Ds_P4.Pt();
                ftree->match_DsFit_Ds_px = match_DsFit_Ds_P4.Px();
                ftree->match_DsFit_Ds_py = match_DsFit_Ds_P4.Py();
                ftree->match_DsFit_Ds_pz = match_DsFit_Ds_P4.Pz();
                ftree->match_DsFit_Ds_invm = match_DsFit_Ds_P4.M();

                TLorentzVector match_DsFit_Mconstraint_phi_P4;
                match_DsFit_Mconstraint_phi_P4.SetXYZM(match_DsFit_phi_P4.Px(), match_DsFit_phi_P4.Py(), match_DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector match_DsFit_Mconstraint_Ds_P4 = match_DsFit_Mconstraint_phi_P4 + match_DsFit_pi_P4;
                ftree->match_DsFit_Mconstraint_Ds_invm = match_DsFit_Mconstraint_Ds_P4.M();

                // Ds fit PP PL 
                ftree->match_DsFit_Kp_pp = match_DsFit_Kp_P4.Vect().Pt(match_DsFit_phi_P4.Vect());
                ftree->match_DsFit_Kp_pl = match_DsFit_Kp_P4.Vect().Dot(match_DsFit_phi_P4.Vect())/match_DsFit_phi_P4.P();
                ftree->match_DsFit_Km_pp = match_DsFit_Km_P4.Vect().Pt(match_DsFit_phi_P4.Vect());
                ftree->match_DsFit_Km_pl = match_DsFit_Km_P4.Vect().Dot(match_DsFit_phi_P4.Vect())/match_DsFit_phi_P4.P();

                ftree->match_DsFit_phi_pp = match_DsFit_phi_P4.Vect().Pt(match_DsFit_Ds_P4.Vect());
                ftree->match_DsFit_phi_pl = match_DsFit_phi_P4.Vect().Dot(match_DsFit_Ds_P4.Vect())/match_DsFit_Ds_P4.P();
                ftree->match_DsFit_pi_pp = match_DsFit_pi_P4.Vect().Pt(match_DsFit_Ds_P4.Vect());
                ftree->match_DsFit_pi_pl = match_DsFit_pi_P4.Vect().Dot(match_DsFit_Ds_P4.Vect())/match_DsFit_Ds_P4.P();

                // Ds fit dR
                ftree->match_DsFit_dR_Kp_Km = reco::deltaR(match_DsFit_Kp_P4.Eta(), match_DsFit_Kp_P4.Phi(), match_DsFit_Km_P4.Eta(), match_DsFit_Km_P4.Phi());
                ftree->match_DsFit_dR_Kp_phi = reco::deltaR(match_DsFit_Kp_P4.Eta(), match_DsFit_Kp_P4.Phi(), match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi());
                ftree->match_DsFit_dR_Km_phi = reco::deltaR(match_DsFit_Km_P4.Eta(), match_DsFit_Km_P4.Phi(), match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi());
                ftree->match_DsFit_dR_Kp_pi = reco::deltaR(match_DsFit_Kp_P4.Eta(), match_DsFit_Kp_P4.Phi(), match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_Km_pi = reco::deltaR(match_DsFit_Km_P4.Eta(), match_DsFit_Km_P4.Phi(), match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_pi_phi = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi());
                ftree->match_DsFit_dR_Kp_Ds = reco::deltaR(match_DsFit_Kp_P4.Eta(), match_DsFit_Kp_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_Km_Ds = reco::deltaR(match_DsFit_Km_P4.Eta(), match_DsFit_Km_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_phi_Ds = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                ftree->match_DsFit_dR_pi_Ds = reco::deltaR(match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());

                double match_pst_Ds = sqrt((pow(ftree->match_DsFit_Ds_invm,4)+pow(ftree->match_phiFit_phi_invm,4)+pow(Mass_pi,4)-2*pow(ftree->match_DsFit_Ds_invm*ftree->match_phiFit_phi_invm,2)-2*pow(ftree->match_DsFit_Ds_invm*Mass_pi,2)-2*pow(ftree->match_phiFit_phi_invm*Mass_pi,2))/(4*pow(ftree->match_DsFit_Ds_invm,2)));
                double match_Estphi = std::sqrt(pow(ftree->match_phiFit_phi_invm,2) + pow(match_pst_Ds,2));
                double match_Estpi = std::sqrt(pow(Mass_pi,2) + pow(match_pst_Ds,2));
                ftree->match_alpha_Ds = (ftree->match_DsFit_phi_pl - ftree->match_DsFit_pi_pl) / (ftree->match_DsFit_phi_pl + ftree->match_DsFit_pi_pl);
                ftree->match_beta_Ds = std::sqrt(pow(ftree->match_DsFit_Ds_p,2) / (pow(ftree->match_DsFit_Ds_invm,2) + pow(ftree->match_DsFit_Ds_p,2)));
                ftree->match_APvar_Ds = pow(ftree->match_DsFit_phi_pp,2) + pow(ftree->match_beta_Ds*(match_Estphi-match_Estpi)-ftree->match_alpha_Ds*ftree->match_beta_Ds*(match_Estphi+match_Estpi),2)/4;

                GlobalPoint match_Ds_pos = match_Ds_Vertex.position();
                GlobalError match_Ds_cov = match_Ds_Vertex.positionError();
                GlobalVector match_Ds_PDirec(ftree->match_DsFit_Ds_px, ftree->match_DsFit_Ds_py, ftree->match_DsFit_Ds_pz);

                GlobalVector match_Ds_primvtx_pos = match_Ds_pos - primaryvertex_pos;
                ftree->match_Ds_primvtx_FDxy = std::hypot(match_Ds_primvtx_pos.x(), match_Ds_primvtx_pos.y());
                ftree->match_Ds_primvtx_FDz = std::abs(match_Ds_primvtx_pos.z());
                ftree->match_Ds_primvtx_FD = std::hypot(ftree->match_Ds_primvtx_FDxy, ftree->match_Ds_primvtx_FDz);
                GlobalError match_Ds_primvtx_coverr = match_Ds_cov + primaryvertex_cov;
                AlgebraicSymMatrix33 match_Ds_primvtx_cov = match_Ds_primvtx_coverr.matrix();

                ftree->match_Ds_primvtx_FDxyerr = std::sqrt(
                        match_Ds_primvtx_pos.x()*match_Ds_primvtx_pos.x()*match_Ds_primvtx_cov(0,0) +
                        match_Ds_primvtx_pos.y()*match_Ds_primvtx_pos.y()*match_Ds_primvtx_cov(1,1) +
                        2*match_Ds_primvtx_pos.x()*match_Ds_primvtx_pos.y()*match_Ds_primvtx_cov(0,1) ) / ftree->match_Ds_primvtx_FDxy;
                Measurement1D match_Ds_primvtx_FDxy_withErr(ftree->match_Ds_primvtx_FDxy, ftree->match_Ds_primvtx_FDxyerr);
                ftree->match_Ds_primvtx_FDxychi2 = match_Ds_primvtx_FDxy_withErr.significance();

                ftree->match_Ds_primvtx_FDzerr = std::sqrt( match_Ds_primvtx_cov(2,2) );
                Measurement1D match_Ds_primvtx_FDz_withErr(ftree->match_Ds_primvtx_FDz, ftree->match_Ds_primvtx_FDzerr);
                ftree->match_Ds_primvtx_FDzchi2 = match_Ds_primvtx_FDz_withErr.significance();

                ftree->match_Ds_primvtx_FDerr = std::sqrt(
                        match_Ds_primvtx_pos.x()*match_Ds_primvtx_pos.x()*match_Ds_primvtx_cov(0,0) +
                        match_Ds_primvtx_pos.y()*match_Ds_primvtx_pos.y()*match_Ds_primvtx_cov(1,1) +
                        match_Ds_primvtx_pos.z()*match_Ds_primvtx_pos.z()*match_Ds_primvtx_cov(2,2) +
                        2*match_Ds_primvtx_pos.x()*match_Ds_primvtx_pos.y()*match_Ds_primvtx_cov(0,1) +
                        2*match_Ds_primvtx_pos.x()*match_Ds_primvtx_pos.z()*match_Ds_primvtx_cov(0,2) +
                        2*match_Ds_primvtx_pos.y()*match_Ds_primvtx_pos.z()*match_Ds_primvtx_cov(1,2)) / ftree->match_Ds_primvtx_FD;
                Measurement1D match_Ds_primvtx_FD_withErr(ftree->match_Ds_primvtx_FD, ftree->match_Ds_primvtx_FDerr);
                ftree->match_Ds_primvtx_FDchi2 = match_Ds_primvtx_FD_withErr.significance();


                ftree->match_Ds_primvtx_dira = match_Ds_primvtx_pos.dot(match_Ds_PDirec) / ( match_Ds_primvtx_pos.mag() * match_Ds_PDirec.mag() );
                ftree->match_Ds_primvtx_dira_angle = std::acos(ftree->match_Ds_primvtx_dira);

                std::pair<bool, Measurement1D> match_Kp_primvtx_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[0], primaryvertex);
                if( match_Kp_primvtx_IPresult.first ){
                    ftree->match_Kp_primvtx_ip = match_Kp_primvtx_IPresult.second.value();
                    ftree->match_Kp_primvtx_iperr = match_Kp_primvtx_IPresult.second.error();
                    ftree->match_Kp_primvtx_ipchi2 = match_Kp_primvtx_IPresult.second.significance();
                }
                std::pair<bool, Measurement1D> match_Km_primvtx_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[1], primaryvertex);
                if( match_Km_primvtx_IPresult.first ){
                    ftree->match_Km_primvtx_ip = match_Km_primvtx_IPresult.second.value();
                    ftree->match_Km_primvtx_iperr = match_Km_primvtx_IPresult.second.error();
                    ftree->match_Km_primvtx_ipchi2 = match_Km_primvtx_IPresult.second.significance();
                }
                std::pair<bool, Measurement1D> match_pi_primvtx_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[2], primaryvertex);
                if( match_pi_primvtx_IPresult.first ){
                    ftree->match_pi_primvtx_ip = match_pi_primvtx_IPresult.second.value();
                    ftree->match_pi_primvtx_iperr = match_pi_primvtx_IPresult.second.error();
                    ftree->match_pi_primvtx_ipchi2 = match_pi_primvtx_IPresult.second.significance();
                }

                reco::TrackCollection match_PVtracks_noDs;

                for(size_t it=0; it<packedHandle->size(); it++){
                    if(it==static_cast<size_t>(ftree->match_Kp_idx)) continue;
                    if(it==static_cast<size_t>(ftree->match_Km_idx)) continue;
                    if(it==static_cast<size_t>(ftree->match_pi_idx)) continue;
                    const auto& pf = (*packedHandle)[it];
                    if(reco::deltaR(match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.3){
                        ftree->match_Ds_Iso_R0p3 += pf.pt(); 
                    }
                    if(reco::deltaR(match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.4){
                        ftree->match_Ds_Iso_R0p4 += pf.pt(); 
                    }
                    if( ! pf.hasTrackDetails() ) continue;
                    match_PVtracks_noDs.push_back(pf.pseudoTrack());
                }

                ftree->match_Ds_IsoRel_R0p3 = ftree->match_Ds_Iso_R0p3 / ftree->match_DsFit_Ds_pt;
                ftree->match_Ds_IsoRel_R0p4 = ftree->match_Ds_Iso_R0p4 / ftree->match_DsFit_Ds_pt;

                std::vector<TransientVertex> match_pvs_noDs_withBS = revertex->makeVertices(match_PVtracks_noDs, *bsHandle, iSetup, "WithBS");
                if( !(match_pvs_noDs_withBS.empty()) ){
                    reco::Vertex match_pv_noDs_withBS = reco::Vertex(match_pvs_noDs_withBS.front());

                    ftree->match_PV_noDs_withBS_IsValid = match_pv_noDs_withBS.isValid();
                    ftree->match_PV_noDs_withBS_IsFake = match_pv_noDs_withBS.isFake();
                    ftree->match_PV_noDs_withBS_chi2 = match_pv_noDs_withBS.chi2();
                    ftree->match_PV_noDs_withBS_ndof = match_pv_noDs_withBS.ndof();
                    ftree->match_PV_noDs_withBS_chi2ndof = match_pv_noDs_withBS.chi2()/match_pv_noDs_withBS.ndof();
                    ftree->match_PV_noDs_withBS_x = match_pv_noDs_withBS.x();
                    ftree->match_PV_noDs_withBS_y = match_pv_noDs_withBS.y();
                    ftree->match_PV_noDs_withBS_z = match_pv_noDs_withBS.z();
                    ftree->match_PV_noDs_withBS_xerr = match_pv_noDs_withBS.xError();
                    ftree->match_PV_noDs_withBS_yerr = match_pv_noDs_withBS.yError();
                    ftree->match_PV_noDs_withBS_zerr = match_pv_noDs_withBS.zError();
                }

                std::vector<TransientVertex> match_pvs_noDs_noBS = revertex->makeVertices(match_PVtracks_noDs, *bsHandle, iSetup, "noBS");
                if( !(match_pvs_noDs_noBS.empty()) ){
                    reco::Vertex match_pv_noDs_noBS = reco::Vertex(match_pvs_noDs_noBS.front());

                    ftree->match_PV_noDs_noBS_IsValid = match_pv_noDs_noBS.isValid();
                    ftree->match_PV_noDs_noBS_IsFake = match_pv_noDs_noBS.isFake();
                    ftree->match_PV_noDs_noBS_chi2 = match_pv_noDs_noBS.chi2();
                    ftree->match_PV_noDs_noBS_ndof = match_pv_noDs_noBS.ndof();
                    ftree->match_PV_noDs_noBS_chi2ndof = match_pv_noDs_noBS.chi2()/match_pv_noDs_noBS.ndof();
                    ftree->match_PV_noDs_noBS_x = match_pv_noDs_noBS.x();
                    ftree->match_PV_noDs_noBS_y = match_pv_noDs_noBS.y();
                    ftree->match_PV_noDs_noBS_z = match_pv_noDs_noBS.z();
                    ftree->match_PV_noDs_noBS_xerr = match_pv_noDs_noBS.xError();
                    ftree->match_PV_noDs_noBS_yerr = match_pv_noDs_noBS.yError();
                    ftree->match_PV_noDs_noBS_zerr = match_pv_noDs_noBS.zError();

                    GlobalPoint match_PV_noDs_noBS_pos = match_pvs_noDs_noBS.front().position();
                    GlobalError match_PV_noDs_noBS_cov = match_pvs_noDs_noBS.front().positionError();

                    GlobalVector match_Ds_PVnoDs_pos = match_Ds_pos - match_PV_noDs_noBS_pos;
                    ftree->match_Ds_PVnoDs_FDxy = std::hypot(match_Ds_PVnoDs_pos.x(), match_Ds_PVnoDs_pos.y());
                    ftree->match_Ds_PVnoDs_FDz = std::abs(match_Ds_PVnoDs_pos.z());
                    ftree->match_Ds_PVnoDs_FD = std::hypot(ftree->match_Ds_PVnoDs_FDxy, ftree->match_Ds_PVnoDs_FDz);
                    GlobalError match_Ds_PVnoDs_coverr = match_Ds_cov + match_PV_noDs_noBS_cov;
                    AlgebraicSymMatrix33 match_Ds_PVnoDs_cov = match_Ds_PVnoDs_coverr.matrix();

                    ftree->match_Ds_PVnoDs_FDxyerr = std::sqrt(
                            match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_cov(0,0) +
                            match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_cov(1,1) +
                            2*match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_cov(0,1) ) / ftree->match_Ds_PVnoDs_FDxy;
                    Measurement1D match_Ds_PVnoDs_FDxy_withErr(ftree->match_Ds_PVnoDs_FDxy, ftree->match_Ds_PVnoDs_FDxyerr);
                    ftree->match_Ds_PVnoDs_FDxychi2 = match_Ds_PVnoDs_FDxy_withErr.significance();

                    ftree->match_Ds_PVnoDs_FDzerr = std::sqrt( match_Ds_PVnoDs_cov(2,2) );
                    Measurement1D match_Ds_PVnoDs_FDz_withErr(ftree->match_Ds_PVnoDs_FDz, ftree->match_Ds_PVnoDs_FDzerr);
                    ftree->match_Ds_PVnoDs_FDzchi2 = match_Ds_PVnoDs_FDz_withErr.significance();

                    ftree->match_Ds_PVnoDs_FDerr = std::sqrt(
                            match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_cov(0,0) +
                            match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_cov(1,1) +
                            match_Ds_PVnoDs_pos.z()*match_Ds_PVnoDs_pos.z()*match_Ds_PVnoDs_cov(2,2) +
                            2*match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_cov(0,1) +
                            2*match_Ds_PVnoDs_pos.x()*match_Ds_PVnoDs_pos.z()*match_Ds_PVnoDs_cov(0,2) +
                            2*match_Ds_PVnoDs_pos.y()*match_Ds_PVnoDs_pos.z()*match_Ds_PVnoDs_cov(1,2)) / ftree->match_Ds_PVnoDs_FD;
                    Measurement1D match_Ds_PVnoDs_FD_withErr(ftree->match_Ds_PVnoDs_FD, ftree->match_Ds_PVnoDs_FDerr);
                    ftree->match_Ds_PVnoDs_FDchi2 = match_Ds_PVnoDs_FD_withErr.significance();

                    ftree->match_Ds_PVnoDs_dira = match_Ds_PVnoDs_pos.dot(match_Ds_PDirec) / ( match_Ds_PVnoDs_pos.mag() * match_Ds_PDirec.mag() );
                    ftree->match_Ds_PVnoDs_dira_angle = std::acos(ftree->match_Ds_PVnoDs_dira);

                    std::pair<bool, Measurement1D> match_Kp_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[0], match_pv_noDs_noBS);
                    if( match_Kp_PVnoDs_IPresult.first ){
                        ftree->match_Kp_PVnoDs_ip = match_Kp_PVnoDs_IPresult.second.value();
                        ftree->match_Kp_PVnoDs_iperr = match_Kp_PVnoDs_IPresult.second.error();
                        ftree->match_Kp_PVnoDs_ipchi2 = match_Kp_PVnoDs_IPresult.second.significance();
                    }
                    std::pair<bool, Measurement1D> match_Km_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[1], match_pv_noDs_noBS);
                    if( match_Km_PVnoDs_IPresult.first ){
                        ftree->match_Km_PVnoDs_ip = match_Km_PVnoDs_IPresult.second.value();
                        ftree->match_Km_PVnoDs_iperr = match_Km_PVnoDs_IPresult.second.error();
                        ftree->match_Km_PVnoDs_ipchi2 = match_Km_PVnoDs_IPresult.second.significance();
                    }
                    std::pair<bool, Measurement1D> match_pi_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(match_DsFit_Tracks[2], match_pv_noDs_noBS);
                    if( match_pi_PVnoDs_IPresult.first ){
                        ftree->match_pi_PVnoDs_ip = match_pi_PVnoDs_IPresult.second.value();
                        ftree->match_pi_PVnoDs_iperr = match_pi_PVnoDs_IPresult.second.error();
                        ftree->match_pi_PVnoDs_ipchi2 = match_pi_PVnoDs_IPresult.second.significance();
                    }
                }
                ftree->Match_Fill_Vector();
            }
        }
    }

    idx_Kp_vec.clear();
    idx_Km_vec.clear();
    idx_pi_vec.clear();

    for(size_t i=0; i<packedHandle->size(); i++){
        const auto& pf = (*packedHandle)[i];

        if(!(pf.trackHighPurity())) continue;
        if(!(pf.hasTrackDetails())) continue;
        if(pf.pt() < 0.5) continue;
        if(pf.p() < 1) continue;

        if(pf.pdgId() == 211){
            idx_Kp_vec.push_back(i);
            idx_pi_vec.push_back(i);
        } else if(pf.pdgId() == -211){
            idx_Km_vec.push_back(i);
        }
    }

    for(size_t i=0; i<idx_Kp_vec.size(); i++){

        const auto& Kp_PF = (*packedHandle)[idx_Kp_vec[i]];

        ftree->Kp_Reset();

        ftree->Kp_isIsolatedChargedHadron = Kp_PF.isIsolatedChargedHadron();
        ftree->Kp_charge                  = Kp_PF.charge();
        ftree->Kp_eta                     = Kp_PF.eta();
        ftree->Kp_phi                     = Kp_PF.phi();
        ftree->Kp_vx                      = Kp_PF.vx();
        ftree->Kp_vy                      = Kp_PF.vy();
        ftree->Kp_vz                      = Kp_PF.vz();
        ftree->Kp_p                       = Kp_PF.p();
        ftree->Kp_pt                      = Kp_PF.pt();
        ftree->Kp_px                      = Kp_PF.px();
        ftree->Kp_py                      = Kp_PF.py();
        ftree->Kp_pz                      = Kp_PF.pz();

        TLorentzVector Kp_P4;
        Kp_P4.SetXYZM(Kp_PF.px(), Kp_PF.py(), Kp_PF.pz(), Mass_K);

        for(size_t j=0; j<idx_Km_vec.size(); j++){

            const auto& Km_PF = (*packedHandle)[idx_Km_vec[j]];

            ftree->Km_Reset();

            ftree->Km_isIsolatedChargedHadron = Km_PF.isIsolatedChargedHadron();
            ftree->Km_charge = Km_PF.charge();
            ftree->Km_eta    = Km_PF.eta();
            ftree->Km_phi    = Km_PF.phi();
            ftree->Km_vx     = Km_PF.vx();
            ftree->Km_vy     = Km_PF.vy();
            ftree->Km_vz     = Km_PF.vz();
            ftree->Km_p      = Km_PF.p();
            ftree->Km_pt     = Km_PF.pt();
            ftree->Km_px     = Km_PF.px();
            ftree->Km_py     = Km_PF.py();
            ftree->Km_pz     = Km_PF.pz();

            TLorentzVector Km_P4;
            Km_P4.SetXYZM(Km_PF.px(), Km_PF.py(), Km_PF.pz(), Mass_K);


            TLorentzVector phi_P4 = Kp_P4 + Km_P4;
            ftree->phi_eta = phi_P4.Eta();
            ftree->phi_phi = phi_P4.Phi();
            ftree->phi_p = phi_P4.P();
            ftree->phi_pt = phi_P4.Pt();
            ftree->phi_px = phi_P4.Px();
            ftree->phi_py = phi_P4.Py();
            ftree->phi_pz = phi_P4.Pz();
            ftree->phi_invm = phi_P4.M();

            ftree->Kp_pp = Kp_P4.Vect().Pt(phi_P4.Vect());
            ftree->Kp_pl = Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            ftree->Km_pp = Km_P4.Vect().Pt(phi_P4.Vect());
            ftree->Km_pl = Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();

            ftree->dR_Kp_Km  = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Km_P4.Eta(),  Km_P4.Phi());
            ftree->dR_Kp_phi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());
            ftree->dR_Km_phi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());

            ftree->dxy_Kp_Km = sqrt( pow( Kp_PF.vx()-Km_PF.vx(), 2 ) + pow( Kp_PF.vy()-Km_PF.vy(), 2 ) );
            ftree->dz_Kp_Km = abs( Kp_PF.vz()-Km_PF.vz() );

            if( abs( Kp_PF.pt() - Km_PF.pt() ) > 20 ) continue;
            if( ftree->dR_Kp_Km > 0.1 ) continue;

            ftree->phi_Reset();
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

            ftree->phiFit_chi2     = phi_Vertex.totalChiSquared();
            ftree->phiFit_ndof     = phi_Vertex.degreesOfFreedom();
            ftree->phiFit_chi2ndof = phi_Vertex.normalisedChiSquared();
            ftree->phiFit_vx       = phi_Vertex.position().x();
            ftree->phiFit_vy       = phi_Vertex.position().y();
            ftree->phiFit_vz       = phi_Vertex.position().z();
            ftree->phiFit_vxerr    = std::sqrt(phi_Vertex.positionError().cxx());
            ftree->phiFit_vyerr    = std::sqrt(phi_Vertex.positionError().cyy());
            ftree->phiFit_vzerr    = std::sqrt(phi_Vertex.positionError().czz());

            TLorentzVector phiFit_Kp_P4;
            phiFit_Kp_P4.SetXYZM(phiFit_Tracks[0].track().px(), phiFit_Tracks[0].track().py(), phiFit_Tracks[0].track().pz(), Mass_K);
            ftree->phiFit_Kp_eta = phiFit_Kp_P4.Eta();
            ftree->phiFit_Kp_phi = phiFit_Kp_P4.Phi();
            ftree->phiFit_Kp_p   = phiFit_Kp_P4.P();
            ftree->phiFit_Kp_pt  = phiFit_Kp_P4.Pt();
            ftree->phiFit_Kp_px  = phiFit_Kp_P4.Px();
            ftree->phiFit_Kp_py  = phiFit_Kp_P4.Py();
            ftree->phiFit_Kp_pz  = phiFit_Kp_P4.Pz();

            TLorentzVector phiFit_Km_P4;
            phiFit_Km_P4.SetXYZM(phiFit_Tracks[1].track().px(), phiFit_Tracks[1].track().py(), phiFit_Tracks[1].track().pz(), Mass_K);
            ftree->phiFit_Km_eta = phiFit_Km_P4.Eta();
            ftree->phiFit_Km_phi = phiFit_Km_P4.Phi();
            ftree->phiFit_Km_p   = phiFit_Km_P4.P();
            ftree->phiFit_Km_pt  = phiFit_Km_P4.Pt();
            ftree->phiFit_Km_px  = phiFit_Km_P4.Px();
            ftree->phiFit_Km_py  = phiFit_Km_P4.Py();
            ftree->phiFit_Km_pz  = phiFit_Km_P4.Pz();

            TLorentzVector phiFit_phi_P4 = phiFit_Kp_P4 + phiFit_Km_P4;
            ftree->phiFit_phi_eta  = phiFit_phi_P4.Eta();
            ftree->phiFit_phi_phi  = phiFit_phi_P4.Phi();
            ftree->phiFit_phi_p    = phiFit_phi_P4.P();
            ftree->phiFit_phi_pt   = phiFit_phi_P4.Pt();
            ftree->phiFit_phi_px   = phiFit_phi_P4.Px();
            ftree->phiFit_phi_py   = phiFit_phi_P4.Py();
            ftree->phiFit_phi_pz   = phiFit_phi_P4.Pz();
            ftree->phiFit_phi_invm = phiFit_phi_P4.M();

            ftree->phiFit_Kp_pp = phiFit_Kp_P4.Vect().Pt(phiFit_phi_P4.Vect());
            ftree->phiFit_Kp_pl = phiFit_Kp_P4.Vect().Dot(phiFit_phi_P4.Vect())/phiFit_phi_P4.P();
            ftree->phiFit_Km_pp = phiFit_Km_P4.Vect().Pt(phiFit_phi_P4.Vect());
            ftree->phiFit_Km_pl = phiFit_Km_P4.Vect().Dot(phiFit_phi_P4.Vect())/phiFit_phi_P4.P();

            ftree->phiFit_dR_Kp_Km  = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_Km_P4.Eta(),  phiFit_Km_P4.Phi());
            ftree->phiFit_dR_Kp_phi = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());
            ftree->phiFit_dR_Km_phi = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());

            ftree->dxy_Kp_phi = sqrt( pow( Kp_PF.vx()-phi_Vertex.position().x(), 2 ) + pow( Kp_PF.vy()-phi_Vertex.position().y(), 2 ) );
            ftree->dxy_Km_phi = sqrt( pow( Km_PF.vx()-phi_Vertex.position().x(), 2 ) + pow( Km_PF.vy()-phi_Vertex.position().y(), 2 ) );

            ftree->dz_Kp_phi = abs( Kp_PF.vz()-phi_Vertex.position().z() );
            ftree->dz_Km_phi = abs( Km_PF.vz()-phi_Vertex.position().z() );

            if( ftree->phiFit_dR_Kp_phi > 0.05 ) continue;
            if( ftree->phiFit_dR_Km_phi > 0.05 ) continue;
            if( ftree->phiFit_phi_invm < 0.99 ) continue;
            if( ftree->phiFit_phi_invm > 1.05 ) continue;

            ftree->alpha_phi = (ftree->phiFit_Kp_pl - ftree->phiFit_Km_pl) / (ftree->phiFit_Kp_pl + ftree->phiFit_Km_pl);
            ftree->beta_phi = std::sqrt(pow(ftree->phiFit_phi_p,2) / (pow(ftree->phiFit_phi_invm,2) + pow(ftree->phiFit_phi_p,2)));
            ftree->APvar_phi = pow(ftree->phiFit_Kp_pp,2) + pow(ftree->alpha_phi*ftree->beta_phi*ftree->phiFit_phi_invm,2)/4;
            /* if(APvar_phi < 0.006) continue; */
            /* if(APvar_phi > 0.028) continue; */

            ftree->num_reco_phi++;

            for(size_t k=0; k<idx_pi_vec.size(); k++){

                if( idx_pi_vec[k] == idx_Kp_vec[i] ) continue;

                const auto& pi_PF = (*packedHandle)[idx_pi_vec[k]];
                ftree->pi_Reset();

                ftree->pi_isIsolatedChargedHadron = pi_PF.isIsolatedChargedHadron();
                ftree->pi_charge = pi_PF.charge();
                ftree->pi_eta    = pi_PF.eta();
                ftree->pi_phi    = pi_PF.phi();
                ftree->pi_vx     = pi_PF.vx();
                ftree->pi_vy     = pi_PF.vy();
                ftree->pi_vz     = pi_PF.vz();
                ftree->pi_p      = pi_PF.p();
                ftree->pi_pt     = pi_PF.pt();
                ftree->pi_px     = pi_PF.px();
                ftree->pi_py     = pi_PF.py();
                ftree->pi_pz     = pi_PF.pz();
                TLorentzVector pi_P4;
                pi_P4.SetXYZM(pi_PF.px(), pi_PF.py(), pi_PF.pz(), Mass_pi);

                ftree->dR_Kp_pi  = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), pi_P4.Eta(),  pi_P4.Phi());
                ftree->dR_Km_pi  = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), pi_P4.Eta(),  pi_P4.Phi());
                ftree->dR_pi_phi = reco::deltaR(pi_P4.Eta(), pi_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());

                TLorentzVector Ds_P4 = Kp_P4 + Km_P4 + pi_P4;
                ftree->Ds_eta  = Ds_P4.Eta();
                ftree->Ds_phi  = Ds_P4.Phi();
                ftree->Ds_p    = Ds_P4.P();
                ftree->Ds_pt   = Ds_P4.Pt();
                ftree->Ds_px   = Ds_P4.Px();
                ftree->Ds_py   = Ds_P4.Py();
                ftree->Ds_pz   = Ds_P4.Pz();
                ftree->Ds_invm = Ds_P4.M();

                ftree->dR_Kp_Ds  = reco::deltaR(Kp_P4.Eta(),  Kp_P4.Phi(),  Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_Km_Ds  = reco::deltaR(Km_P4.Eta(),  Km_P4.Phi(),  Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_pi_Ds  = reco::deltaR(pi_P4.Eta(),  pi_P4.Phi(),  Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_phi_Ds = reco::deltaR(phi_P4.Eta(), phi_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());

                ftree->pi_pp  = pi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->pi_pl  = pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
                ftree->phi_pp = phi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->phi_pl = phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();

                TLorentzVector phiFit_pi_P4 = pi_P4;
                ftree->phiFit_pi_eta = phiFit_pi_P4.Eta();
                ftree->phiFit_pi_phi = phiFit_pi_P4.Phi();
                ftree->phiFit_pi_p   = phiFit_pi_P4.P();
                ftree->phiFit_pi_pt  = phiFit_pi_P4.Pt();
                ftree->phiFit_pi_px  = phiFit_pi_P4.Px();
                ftree->phiFit_pi_py  = phiFit_pi_P4.Py();
                ftree->phiFit_pi_pz  = phiFit_pi_P4.Pz();

                ftree->phiFit_dR_Kp_pi  = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi());
                ftree->phiFit_dR_Km_pi  = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi());
                ftree->phiFit_dR_pi_phi = reco::deltaR(phiFit_pi_P4.Eta(), phiFit_pi_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());

                TLorentzVector phiFit_Ds_P4 = phiFit_Kp_P4 + phiFit_Km_P4 + phiFit_pi_P4;
                ftree->phiFit_Ds_eta  = phiFit_Ds_P4.Eta();
                ftree->phiFit_Ds_phi  = phiFit_Ds_P4.Phi();
                ftree->phiFit_Ds_p    = phiFit_Ds_P4.P();
                ftree->phiFit_Ds_pt   = phiFit_Ds_P4.Pt();
                ftree->phiFit_Ds_px   = phiFit_Ds_P4.Px();
                ftree->phiFit_Ds_py   = phiFit_Ds_P4.Py();
                ftree->phiFit_Ds_pz   = phiFit_Ds_P4.Pz();
                ftree->phiFit_Ds_invm = phiFit_Ds_P4.M();

                ftree->phiFit_pi_pp  = phiFit_pi_P4.Vect().Pt(phiFit_Ds_P4.Vect());
                ftree->phiFit_pi_pl  = phiFit_pi_P4.Vect().Dot(phiFit_Ds_P4.Vect())/phiFit_Ds_P4.P();
                ftree->phiFit_phi_pp = phiFit_phi_P4.Vect().Pt(phiFit_Ds_P4.Vect());
                ftree->phiFit_phi_pl = phiFit_phi_P4.Vect().Dot(phiFit_Ds_P4.Vect())/phiFit_Ds_P4.P();

                ftree->phiFit_dR_Kp_Ds  = reco::deltaR(phiFit_Kp_P4.Eta(),  phiFit_Kp_P4.Phi(),  phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_Km_Ds  = reco::deltaR(phiFit_Km_P4.Eta(),  phiFit_Km_P4.Phi(),  phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_pi_Ds  = reco::deltaR(phiFit_pi_P4.Eta(),  phiFit_pi_P4.Phi(),  phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_phi_Ds = reco::deltaR(phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi(), phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());

                if( ftree->dR_Kp_pi > 0.4 ) continue;
                if( ftree->dR_Km_pi > 0.4 ) continue;
                if( ftree->phiFit_dR_pi_phi > 0.4 ) continue;

                ftree->dxy_Kp_pi  = sqrt( pow( Kp_PF.vx()-pi_PF.vx(), 2 ) + pow( Kp_PF.vy()-pi_PF.vy(), 2 ) );
                ftree->dxy_Km_pi  = sqrt( pow( Km_PF.vx()-pi_PF.vx(), 2 ) + pow( Km_PF.vy()-pi_PF.vy(), 2 ) );
                ftree->dxy_pi_phi = sqrt( pow( pi_PF.vx()-phi_Vertex.position().x(), 2 ) + pow( pi_PF.vy()-phi_Vertex.position().y(), 2 ) );

                ftree->dz_Kp_pi  = abs( Kp_PF.vz()-pi_PF.vz() );
                ftree->dz_Km_pi  = abs( Km_PF.vz()-pi_PF.vz() );
                ftree->dz_pi_phi = abs( pi_PF.vz()-phi_Vertex.position().z() );

                ftree->Ds_Reset();
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

                ftree->DsFit_chi2     = Ds_Vertex.totalChiSquared();
                ftree->DsFit_ndof     = Ds_Vertex.degreesOfFreedom();
                ftree->DsFit_chi2ndof = Ds_Vertex.normalisedChiSquared();
                ftree->DsFit_vx       = Ds_Vertex.position().x();
                ftree->DsFit_vy       = Ds_Vertex.position().y();
                ftree->DsFit_vz       = Ds_Vertex.position().z();
                ftree->DsFit_vxerr    = std::sqrt(Ds_Vertex.positionError().cxx());
                ftree->DsFit_vyerr    = std::sqrt(Ds_Vertex.positionError().cyy());
                ftree->DsFit_vzerr    = std::sqrt(Ds_Vertex.positionError().czz());

                std::vector<reco::TransientTrack> DsFit_Tracks = Ds_Vertex.refittedTracks();

                TLorentzVector DsFit_Kp_P4;
                DsFit_Kp_P4.SetXYZM(DsFit_Tracks[0].track().px(), DsFit_Tracks[0].track().py(), DsFit_Tracks[0].track().pz(), Mass_K);
                ftree->DsFit_Kp_eta = DsFit_Kp_P4.Eta();
                ftree->DsFit_Kp_phi = DsFit_Kp_P4.Phi();
                ftree->DsFit_Kp_p   = DsFit_Kp_P4.P();
                ftree->DsFit_Kp_pt  = DsFit_Kp_P4.Pt();
                ftree->DsFit_Kp_px  = DsFit_Kp_P4.Px();
                ftree->DsFit_Kp_py  = DsFit_Kp_P4.Py();
                ftree->DsFit_Kp_pz  = DsFit_Kp_P4.Pz();

                // Ds fit Km
                TLorentzVector DsFit_Km_P4;
                DsFit_Km_P4.SetXYZM(DsFit_Tracks[1].track().px(), DsFit_Tracks[1].track().py(), DsFit_Tracks[1].track().pz(), Mass_K);
                ftree->DsFit_Km_eta = DsFit_Km_P4.Eta();
                ftree->DsFit_Km_phi = DsFit_Km_P4.Phi();
                ftree->DsFit_Km_p   = DsFit_Km_P4.P();
                ftree->DsFit_Km_pt  = DsFit_Km_P4.Pt();
                ftree->DsFit_Km_px  = DsFit_Km_P4.Px();
                ftree->DsFit_Km_py  = DsFit_Km_P4.Py();
                ftree->DsFit_Km_pz  = DsFit_Km_P4.Pz();

                // Ds fit pi
                TLorentzVector DsFit_pi_P4;
                DsFit_pi_P4.SetXYZM(DsFit_Tracks[2].track().px(), DsFit_Tracks[2].track().py(), DsFit_Tracks[2].track().pz(), Mass_pi);
                ftree->DsFit_pi_eta = DsFit_pi_P4.Eta();
                ftree->DsFit_pi_phi = DsFit_pi_P4.Phi();
                ftree->DsFit_pi_p   = DsFit_pi_P4.P();
                ftree->DsFit_pi_pt  = DsFit_pi_P4.Pt();
                ftree->DsFit_pi_px  = DsFit_pi_P4.Px();
                ftree->DsFit_pi_py  = DsFit_pi_P4.Py();
                ftree->DsFit_pi_pz  = DsFit_pi_P4.Pz();

                // Ds fit phi
                TLorentzVector DsFit_phi_P4 = DsFit_Kp_P4 + DsFit_Km_P4;
                ftree->DsFit_phi_eta  = DsFit_phi_P4.Eta();
                ftree->DsFit_phi_phi  = DsFit_phi_P4.Phi();
                ftree->DsFit_phi_p    = DsFit_phi_P4.P();
                ftree->DsFit_phi_pt   = DsFit_phi_P4.Pt();
                ftree->DsFit_phi_px   = DsFit_phi_P4.Px();
                ftree->DsFit_phi_py   = DsFit_phi_P4.Py();
                ftree->DsFit_phi_pz   = DsFit_phi_P4.Pz();
                ftree->DsFit_phi_invm = DsFit_phi_P4.M();

                // Ds fit Ds
                TLorentzVector DsFit_Ds_P4 = DsFit_Kp_P4 + DsFit_Km_P4 + DsFit_pi_P4;
                ftree->DsFit_Ds_eta  = DsFit_Ds_P4.Eta();
                ftree->DsFit_Ds_phi  = DsFit_Ds_P4.Phi();
                ftree->DsFit_Ds_p    = DsFit_Ds_P4.P();
                ftree->DsFit_Ds_pt   = DsFit_Ds_P4.Pt();
                ftree->DsFit_Ds_px   = DsFit_Ds_P4.Px();
                ftree->DsFit_Ds_py   = DsFit_Ds_P4.Py();
                ftree->DsFit_Ds_pz   = DsFit_Ds_P4.Pz();
                ftree->DsFit_Ds_invm = DsFit_Ds_P4.M();

                TLorentzVector DsFit_Mconstraint_phi_P4;
                DsFit_Mconstraint_phi_P4.SetXYZM(DsFit_phi_P4.Px(), DsFit_phi_P4.Py(), DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector DsFit_Mconstraint_Ds_P4 = DsFit_Mconstraint_phi_P4 + DsFit_pi_P4;
                ftree->DsFit_Mconstraint_Ds_invm = DsFit_Mconstraint_Ds_P4.M();

                ftree->DsFit_Kp_pp = DsFit_Kp_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Kp_pl = DsFit_Kp_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();
                ftree->DsFit_Km_pp = DsFit_Km_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Km_pl = DsFit_Km_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();

                ftree->DsFit_phi_pp = DsFit_phi_P4.Vect().Pt(DsFit_Ds_P4.Vect());
                ftree->DsFit_phi_pl = DsFit_phi_P4.Vect().Dot(DsFit_Ds_P4.Vect())/DsFit_Ds_P4.P();
                ftree->DsFit_pi_pp  = DsFit_pi_P4.Vect().Pt(DsFit_Ds_P4.Vect());
                ftree->DsFit_pi_pl  = DsFit_pi_P4.Vect().Dot(DsFit_Ds_P4.Vect())/DsFit_Ds_P4.P();

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

                ftree->dxy_Kp_Ds  = sqrt( pow( Kp_PF.vx()-Ds_Vertex.position().x(), 2 ) + pow( Kp_PF.vy()-Ds_Vertex.position().y(), 2 ) );
                ftree->dxy_Km_Ds  = sqrt( pow( Km_PF.vx()-Ds_Vertex.position().x(), 2 ) + pow( Km_PF.vy()-Ds_Vertex.position().y(), 2 ) );
                ftree->dxy_pi_Ds  = sqrt( pow( pi_PF.vx()-Ds_Vertex.position().x(), 2 ) + pow( pi_PF.vy()-Ds_Vertex.position().y(), 2 ) );
                ftree->dxy_phi_Ds = sqrt( pow( phi_Vertex.position().x()-Ds_Vertex.position().x(), 2 ) + pow( phi_Vertex.position().y()-Ds_Vertex.position().y(), 2 ) );

                ftree->dz_Kp_Ds  = abs( Kp_PF.vz()-Ds_Vertex.position().z() );
                ftree->dz_Km_Ds  = abs( Km_PF.vz()-Ds_Vertex.position().z() );
                ftree->dz_pi_Ds  = abs( pi_PF.vz()-Ds_Vertex.position().z() );
                ftree->dz_phi_Ds = abs( phi_Vertex.position().z()-Ds_Vertex.position().z() );

                if( ftree->DsFit_dR_Kp_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_Km_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_pi_Ds > 0.4 ) continue;
                if( ftree->DsFit_dR_phi_Ds > 0.15 ) continue;
                if( ftree->DsFit_Ds_invm < 1.85) continue;
                if( ftree->DsFit_Ds_invm > 2.1) continue;

                double pst_Ds = sqrt((pow(ftree->DsFit_Ds_invm,4)+pow(ftree->phiFit_phi_invm,4)+pow(Mass_pi,4)-2*pow(ftree->DsFit_Ds_invm*ftree->phiFit_phi_invm,2)-2*pow(ftree->DsFit_Ds_invm*Mass_pi,2)-2*pow(ftree->phiFit_phi_invm*Mass_pi,2))/(4*pow(ftree->DsFit_Ds_invm,2)));
                double Est_phi = std::sqrt(pow(ftree->phiFit_phi_invm,2) + pow(pst_Ds,2));
                double Est_pi = std::sqrt(pow(Mass_pi,2) + pow(pst_Ds,2));
                ftree->alpha_Ds = (ftree->DsFit_phi_pl - ftree->DsFit_pi_pl) / (ftree->DsFit_phi_pl + ftree->DsFit_pi_pl);
                ftree->beta_Ds = std::sqrt(pow(ftree->DsFit_Ds_p,2) / (pow(ftree->DsFit_Ds_invm,2) + pow(ftree->DsFit_Ds_p,2)));
                ftree->APvar_Ds = pow(ftree->DsFit_phi_pp,2) + pow(ftree->beta_Ds*(Est_phi-Est_pi)-ftree->alpha_Ds*ftree->beta_Ds*(Est_phi+Est_pi),2)/4;
                /* if(ftree->APvar_Ds < 0.4) continue; */
                /* if(ftree->APvar_Ds > 0.6) continue; */

                GlobalPoint Ds_pos = Ds_Vertex.position();
                GlobalError Ds_cov = Ds_Vertex.positionError();
                GlobalVector Ds_PDirec(ftree->DsFit_Ds_px, ftree->DsFit_Ds_py, ftree->DsFit_Ds_pz);

                GlobalVector Ds_primvtx_pos = Ds_pos - primaryvertex_pos;
                ftree->Ds_primvtx_FDxy = std::hypot(Ds_primvtx_pos.x(), Ds_primvtx_pos.y());
                ftree->Ds_primvtx_FDz = std::abs(Ds_primvtx_pos.z());
                ftree->Ds_primvtx_FD = std::hypot(ftree->Ds_primvtx_FDxy, ftree->Ds_primvtx_FDz);
                GlobalError Ds_primvtx_coverr = Ds_cov + primaryvertex_cov;
                AlgebraicSymMatrix33 Ds_primvtx_cov = Ds_primvtx_coverr.matrix();

                ftree->Ds_primvtx_FDxyerr = std::sqrt(
                        Ds_primvtx_pos.x()*Ds_primvtx_pos.x()*Ds_primvtx_cov(0,0) +
                        Ds_primvtx_pos.y()*Ds_primvtx_pos.y()*Ds_primvtx_cov(1,1) +
                        2*Ds_primvtx_pos.x()*Ds_primvtx_pos.y()*Ds_primvtx_cov(0,1) ) / ftree->Ds_primvtx_FDxy;
                Measurement1D Ds_primvtx_FDxy_withErr(ftree->Ds_primvtx_FDxy, ftree->Ds_primvtx_FDxyerr);
                ftree->Ds_primvtx_FDxychi2 = Ds_primvtx_FDxy_withErr.significance();

                ftree->Ds_primvtx_FDzerr = std::sqrt( Ds_primvtx_cov(2,2) );
                Measurement1D Ds_primvtx_FDz_withErr(ftree->Ds_primvtx_FDz, ftree->Ds_primvtx_FDzerr);
                ftree->Ds_primvtx_FDzchi2 = Ds_primvtx_FDz_withErr.significance();

                ftree->Ds_primvtx_FDerr = std::sqrt(
                        Ds_primvtx_pos.x()*Ds_primvtx_pos.x()*Ds_primvtx_cov(0,0) +
                        Ds_primvtx_pos.y()*Ds_primvtx_pos.y()*Ds_primvtx_cov(1,1) +
                        Ds_primvtx_pos.z()*Ds_primvtx_pos.z()*Ds_primvtx_cov(2,2) +
                        2*Ds_primvtx_pos.x()*Ds_primvtx_pos.y()*Ds_primvtx_cov(0,1) +
                        2*Ds_primvtx_pos.x()*Ds_primvtx_pos.z()*Ds_primvtx_cov(0,2) +
                        2*Ds_primvtx_pos.y()*Ds_primvtx_pos.z()*Ds_primvtx_cov(1,2)) / ftree->Ds_primvtx_FD;
                Measurement1D Ds_primvtx_FD_withErr(ftree->Ds_primvtx_FD, ftree->Ds_primvtx_FDerr);
                ftree->Ds_primvtx_FDchi2 = Ds_primvtx_FD_withErr.significance();


                ftree->Ds_primvtx_dira = Ds_primvtx_pos.dot(Ds_PDirec) / ( Ds_primvtx_pos.mag() * Ds_PDirec.mag() );
                ftree->Ds_primvtx_dira_angle = std::acos(ftree->Ds_primvtx_dira);

                std::pair<bool, Measurement1D> Kp_primvtx_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[0], primaryvertex);
                if( Kp_primvtx_IPresult.first ){
                    ftree->Kp_primvtx_ip = Kp_primvtx_IPresult.second.value();
                    ftree->Kp_primvtx_iperr = Kp_primvtx_IPresult.second.error();
                    ftree->Kp_primvtx_ipchi2 = Kp_primvtx_IPresult.second.significance();
                }
                std::pair<bool, Measurement1D> Km_primvtx_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[1], primaryvertex);
                if( Km_primvtx_IPresult.first ){
                    ftree->Km_primvtx_ip = Km_primvtx_IPresult.second.value();
                    ftree->Km_primvtx_iperr = Km_primvtx_IPresult.second.error();
                    ftree->Km_primvtx_ipchi2 = Km_primvtx_IPresult.second.significance();
                }
                std::pair<bool, Measurement1D> pi_primvtx_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[2], primaryvertex);
                if( pi_primvtx_IPresult.first ){
                    ftree->pi_primvtx_ip = pi_primvtx_IPresult.second.value();
                    ftree->pi_primvtx_iperr = pi_primvtx_IPresult.second.error();
                    ftree->pi_primvtx_ipchi2 = pi_primvtx_IPresult.second.significance();
                }

                ftree->num_reco_Ds++;



                ftree->PV_noDs_Reset();

                reco::TrackCollection PVtracks_noDs;
                for(size_t it=0; it<packedHandle->size(); it++){
                    if(it==static_cast<size_t>(idx_Kp_vec[i])) continue;
                    if(it==static_cast<size_t>(idx_Km_vec[j])) continue;
                    if(it==static_cast<size_t>(idx_pi_vec[k])) continue;
                    const auto& pf = (*packedHandle)[it];
                    if(reco::deltaR(DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.3){
                        ftree->Ds_Iso_R0p3 += pf.pt(); 
                    }
                    if(reco::deltaR(DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi(), pf.eta(), pf.phi()) < 0.4){
                        ftree->Ds_Iso_R0p4 += pf.pt(); 
                    }
                    if( ! pf.hasTrackDetails() ) continue;
                    PVtracks_noDs.push_back(pf.pseudoTrack());
                }

                ftree->Ds_IsoRel_R0p3 = ftree->Ds_Iso_R0p3 / ftree->DsFit_Ds_pt;
                ftree->Ds_IsoRel_R0p4 = ftree->Ds_Iso_R0p4 / ftree->DsFit_Ds_pt;

                std::vector<TransientVertex> pvs_noDs_withBS = revertex->makeVertices(PVtracks_noDs, *bsHandle, iSetup, "WithBS");
                if( !(pvs_noDs_withBS.empty()) ){
                    reco::Vertex pv_noDs_withBS = reco::Vertex(pvs_noDs_withBS.front());

                    ftree->PV_noDs_withBS_IsValid = pv_noDs_withBS.isValid();
                    ftree->PV_noDs_withBS_IsFake = pv_noDs_withBS.isFake();
                    ftree->PV_noDs_withBS_chi2 = pv_noDs_withBS.chi2();
                    ftree->PV_noDs_withBS_ndof = pv_noDs_withBS.ndof();
                    ftree->PV_noDs_withBS_chi2ndof = pv_noDs_withBS.chi2()/pv_noDs_withBS.ndof();
                    ftree->PV_noDs_withBS_x = pv_noDs_withBS.x();
                    ftree->PV_noDs_withBS_y = pv_noDs_withBS.y();
                    ftree->PV_noDs_withBS_z = pv_noDs_withBS.z();
                    ftree->PV_noDs_withBS_xerr = pv_noDs_withBS.xError();
                    ftree->PV_noDs_withBS_yerr = pv_noDs_withBS.yError();
                    ftree->PV_noDs_withBS_zerr = pv_noDs_withBS.zError();
                }

                std::vector<TransientVertex> pvs_noDs_noBS = revertex->makeVertices(PVtracks_noDs, *bsHandle, iSetup, "noBS");
                if( !(pvs_noDs_noBS.empty()) ){
                    reco::Vertex pv_noDs_noBS = reco::Vertex(pvs_noDs_noBS.front());

                    ftree->PV_noDs_noBS_IsValid = pv_noDs_noBS.isValid();
                    ftree->PV_noDs_noBS_IsFake = pv_noDs_noBS.isFake();
                    ftree->PV_noDs_noBS_chi2 = pv_noDs_noBS.chi2();
                    ftree->PV_noDs_noBS_ndof = pv_noDs_noBS.ndof();
                    ftree->PV_noDs_noBS_chi2ndof = pv_noDs_noBS.chi2()/pv_noDs_noBS.ndof();
                    ftree->PV_noDs_noBS_x = pv_noDs_noBS.x();
                    ftree->PV_noDs_noBS_y = pv_noDs_noBS.y();
                    ftree->PV_noDs_noBS_z = pv_noDs_noBS.z();
                    ftree->PV_noDs_noBS_xerr = pv_noDs_noBS.xError();
                    ftree->PV_noDs_noBS_yerr = pv_noDs_noBS.yError();
                    ftree->PV_noDs_noBS_zerr = pv_noDs_noBS.zError();

                    GlobalPoint PV_noDs_noBS_pos = pvs_noDs_noBS.front().position();
                    GlobalError PV_noDs_noBS_cov = pvs_noDs_noBS.front().positionError();

                    GlobalVector Ds_PVnoDs_pos = Ds_pos - PV_noDs_noBS_pos;
                    ftree->Ds_PVnoDs_FDxy = std::hypot(Ds_PVnoDs_pos.x(), Ds_PVnoDs_pos.y());
                    ftree->Ds_PVnoDs_FDz = std::abs(Ds_PVnoDs_pos.z());
                    ftree->Ds_PVnoDs_FD = std::hypot(ftree->Ds_PVnoDs_FDxy, ftree->Ds_PVnoDs_FDz);
                    GlobalError Ds_PVnoDs_coverr = Ds_cov + PV_noDs_noBS_cov;
                    AlgebraicSymMatrix33 Ds_PVnoDs_cov = Ds_PVnoDs_coverr.matrix();

                    ftree->Ds_PVnoDs_FDxyerr = std::sqrt(
                            Ds_PVnoDs_pos.x()*Ds_PVnoDs_pos.x()*Ds_PVnoDs_cov(0,0) +
                            Ds_PVnoDs_pos.y()*Ds_PVnoDs_pos.y()*Ds_PVnoDs_cov(1,1) +
                            2*Ds_PVnoDs_pos.x()*Ds_PVnoDs_pos.y()*Ds_PVnoDs_cov(0,1) ) / ftree->Ds_PVnoDs_FDxy;
                    Measurement1D Ds_PVnoDs_FDxy_withErr(ftree->Ds_PVnoDs_FDxy, ftree->Ds_PVnoDs_FDxyerr);
                    ftree->Ds_PVnoDs_FDxychi2 = Ds_PVnoDs_FDxy_withErr.significance();

                    ftree->Ds_PVnoDs_FDzerr = std::sqrt( Ds_PVnoDs_cov(2,2) );
                    Measurement1D Ds_PVnoDs_FDz_withErr(ftree->Ds_PVnoDs_FDz, ftree->Ds_PVnoDs_FDzerr);
                    ftree->Ds_PVnoDs_FDzchi2 = Ds_PVnoDs_FDz_withErr.significance();

                    ftree->Ds_PVnoDs_FDerr = std::sqrt(
                            Ds_PVnoDs_pos.x()*Ds_PVnoDs_pos.x()*Ds_PVnoDs_cov(0,0) +
                            Ds_PVnoDs_pos.y()*Ds_PVnoDs_pos.y()*Ds_PVnoDs_cov(1,1) +
                            Ds_PVnoDs_pos.z()*Ds_PVnoDs_pos.z()*Ds_PVnoDs_cov(2,2) +
                            2*Ds_PVnoDs_pos.x()*Ds_PVnoDs_pos.y()*Ds_PVnoDs_cov(0,1) +
                            2*Ds_PVnoDs_pos.x()*Ds_PVnoDs_pos.z()*Ds_PVnoDs_cov(0,2) +
                            2*Ds_PVnoDs_pos.y()*Ds_PVnoDs_pos.z()*Ds_PVnoDs_cov(1,2)) / ftree->Ds_PVnoDs_FD;
                    Measurement1D Ds_PVnoDs_FD_withErr(ftree->Ds_PVnoDs_FD, ftree->Ds_PVnoDs_FDerr);
                    ftree->Ds_PVnoDs_FDchi2 = Ds_PVnoDs_FD_withErr.significance();

                    ftree->Ds_PVnoDs_dira = Ds_PVnoDs_pos.dot(Ds_PDirec) / ( Ds_PVnoDs_pos.mag() * Ds_PDirec.mag() );
                    ftree->Ds_PVnoDs_dira_angle = std::acos(ftree->Ds_PVnoDs_dira);

                    std::pair<bool, Measurement1D> Kp_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[0], pv_noDs_noBS);
                    if( Kp_PVnoDs_IPresult.first ){
                        ftree->Kp_PVnoDs_ip = Kp_PVnoDs_IPresult.second.value();
                        ftree->Kp_PVnoDs_iperr = Kp_PVnoDs_IPresult.second.error();
                        ftree->Kp_PVnoDs_ipchi2 = Kp_PVnoDs_IPresult.second.significance();
                    }
                    std::pair<bool, Measurement1D> Km_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[1], pv_noDs_noBS);
                    if( Km_PVnoDs_IPresult.first ){
                        ftree->Km_PVnoDs_ip = Km_PVnoDs_IPresult.second.value();
                        ftree->Km_PVnoDs_iperr = Km_PVnoDs_IPresult.second.error();
                        ftree->Km_PVnoDs_ipchi2 = Km_PVnoDs_IPresult.second.significance();
                    }
                    std::pair<bool, Measurement1D> pi_PVnoDs_IPresult = IPTools::absoluteImpactParameter3D(DsFit_Tracks[2], pv_noDs_noBS);
                    if( pi_PVnoDs_IPresult.first ){
                        ftree->pi_PVnoDs_ip = pi_PVnoDs_IPresult.second.value();
                        ftree->pi_PVnoDs_iperr = pi_PVnoDs_IPresult.second.error();
                        ftree->pi_PVnoDs_ipchi2 = pi_PVnoDs_IPresult.second.significance();
                    }
                }

                if( idx_Kp_vec[i] == ftree->match_Kp_idx ) ftree->Kp_match = true;
                else ftree->Kp_match = false;
                if( idx_Km_vec[j] == ftree->match_Km_idx ) ftree->Km_match = true;
                else ftree->Km_match = false;
                if( idx_pi_vec[k] == ftree->match_pi_idx ) ftree->pi_match = true;
                else ftree->pi_match = false;

                if(ftree->Kp_match && ftree->Km_match && ftree->pi_match) ftree->match_entry = true;
                else ftree->match_entry = false;
                if(!ftree->Kp_match && !ftree->Km_match && !ftree->pi_match) ftree->non_match_entry = true;
                else ftree->non_match_entry = false;

                ftree->Fill_Vector();
            }
        }
    }

    int num_candidates = ftree->num_reco_Ds;

    double maxPT = 0;
    int maxidx = -1;

    for(int i=0; i<num_candidates; i++){
        if(ftree->DsFit_Ds_pt_vec[i]>maxPT){
            maxPT = ftree->DsFit_Ds_pt_vec[i];
            maxidx = i;
        }
    }
    if(maxidx > -1) {
        ftree->Best_Fill_Vector(maxidx);
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
    /* edm::ParameterSetDescription desc; */
    /* desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles")); */
    /* desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));; */
    /* descriptions.add("PVStudy", desc); */
}

//define this as a plug-in
DEFINE_FWK_MODULE(PVStudy);

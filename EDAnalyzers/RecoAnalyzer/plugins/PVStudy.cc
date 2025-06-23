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
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

#include "EDAnalyzers/GenParticleAnalyzer/interface/PVTree.h"


class PVStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit PVStudy(const edm::ParameterSet&);
        ~PVStudy();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        const edm::Service<TFileService> fs;
        PVTree *ftree;

        // dR of hadrons to Gen
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
        const double pst_phi = Mass_phi/2;
        const double pst_Ds = sqrt((pow(Mass_Ds,4)+pow(Mass_phi,4)+pow(Mass_pi,4)-2*pow(Mass_Ds*Mass_phi,2)-2*pow(Mass_Ds*Mass_pi,2)-2*pow(Mass_phi*Mass_pi,2))/(4*pow(Mass_Ds,2)));
        const double Est_phi = std::sqrt(pow(Mass_phi,2) + pow(pst_Ds,2));
        const double Est_pi = std::sqrt(pow(Mass_pi,2) + pow(pst_Ds,2));
};

PVStudy::PVStudy(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    ftree = new PVTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

PVStudy::~PVStudy() {}

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

    edm::Handle<reco::GenParticleCollection> prunedGen;
    iEvent.getByToken(prunedGenToken, prunedGen);

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

    KalmanVertexFitter KVFitter(true);
    AdaptiveVertexFitter AVFitter;

    for(size_t i=0; i<prunedGen->size(); i++){
        const auto& gp = (*prunedGen)[i];
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
    ftree->Match_Reset();

    if(ftree->num_Gen_Kp==1 && ftree->num_Gen_Km==1 && ftree->num_Gen_pi==1 && ftree->num_Gen_phi==1 && ftree->num_Gen_Ds==1 && ftree->num_Gen_mu==1 && ftree->num_Gen_nu==1 && ftree->num_Gen_W==1 && ftree->num_Gen_H==1){

        const auto& Kp_GP = (*prunedGen)[Gen_Kp_idx[0]];
        const auto& Km_GP = (*prunedGen)[Gen_Km_idx[0]];
        const auto& pi_GP = (*prunedGen)[Gen_pi_idx[0]];
        const auto& phi_GP = (*prunedGen)[Gen_phi_idx[0]];
        const auto& Ds_GP = (*prunedGen)[Gen_Ds_idx[0]];
        const auto& mu_GP = (*prunedGen)[Gen_mu_idx[0]];
        const auto& nu_GP = (*prunedGen)[Gen_nu_idx[0]];
        const auto& W_GP = (*prunedGen)[Gen_W_idx[0]];
        const auto& H_GP = (*prunedGen)[Gen_H_idx[0]];

        ftree->Gen_Kp_ETA = Kp_GP.eta();
        ftree->Gen_Kp_PHI = Kp_GP.phi();
        ftree->Gen_Kp_ORIVX_X = Kp_GP.vx();
        ftree->Gen_Kp_ORIVX_Y = Kp_GP.vy();
        ftree->Gen_Kp_ORIVX_Z = Kp_GP.vz();
        ftree->Gen_Kp_P = Kp_GP.p();
        ftree->Gen_Kp_PT = Kp_GP.pt();
        ftree->Gen_Kp_PX = Kp_GP.px();
        ftree->Gen_Kp_PY = Kp_GP.py();
        ftree->Gen_Kp_PZ = Kp_GP.pz();
        TVector3 Gen_Kp_P3(Kp_GP.px(), Kp_GP.py(), Kp_GP.pz());

        ftree->Gen_Km_ETA = Km_GP.eta();
        ftree->Gen_Km_PHI = Km_GP.phi();
        ftree->Gen_Km_ORIVX_X = Km_GP.vx();
        ftree->Gen_Km_ORIVX_Y = Km_GP.vy();
        ftree->Gen_Km_ORIVX_Z = Km_GP.vz();
        ftree->Gen_Km_P = Km_GP.p();
        ftree->Gen_Km_PT = Km_GP.pt();
        ftree->Gen_Km_PX = Km_GP.px();
        ftree->Gen_Km_PY = Km_GP.py();
        ftree->Gen_Km_PZ = Km_GP.pz();
        TVector3 Gen_Km_P3(Km_GP.px(), Km_GP.py(), Km_GP.pz());

        ftree->Gen_pi_ETA = pi_GP.eta();
        ftree->Gen_pi_PHI = pi_GP.phi();
        ftree->Gen_pi_ORIVX_X = pi_GP.vx();
        ftree->Gen_pi_ORIVX_Y = pi_GP.vy();
        ftree->Gen_pi_ORIVX_Z = pi_GP.vz();
        ftree->Gen_pi_P = pi_GP.p();
        ftree->Gen_pi_PT = pi_GP.pt();
        ftree->Gen_pi_PX = pi_GP.px();
        ftree->Gen_pi_PY = pi_GP.py();
        ftree->Gen_pi_PZ = pi_GP.pz();
        TVector3 Gen_pi_P3(pi_GP.px(), pi_GP.py(), pi_GP.pz());

        ftree->Gen_phi_ETA = phi_GP.eta();
        ftree->Gen_phi_PHI = phi_GP.phi();
        ftree->Gen_phi_ORIVX_X = phi_GP.vx();
        ftree->Gen_phi_ORIVX_Y = phi_GP.vy();
        ftree->Gen_phi_ORIVX_Z = phi_GP.vz();
        ftree->Gen_phi_P = phi_GP.p();
        ftree->Gen_phi_PT = phi_GP.pt();
        ftree->Gen_phi_PX = phi_GP.px();
        ftree->Gen_phi_PY = phi_GP.py();
        ftree->Gen_phi_PZ = phi_GP.pz();
        TVector3 Gen_phi_P3(phi_GP.px(), phi_GP.py(), phi_GP.pz());

        ftree->Gen_Ds_ETA = Ds_GP.eta();
        ftree->Gen_Ds_PHI = Ds_GP.phi();
        ftree->Gen_Ds_ORIVX_X = Ds_GP.vx();
        ftree->Gen_Ds_ORIVX_Y = Ds_GP.vy();
        ftree->Gen_Ds_ORIVX_Z = Ds_GP.vz();
        ftree->Gen_Ds_P = Ds_GP.p();
        ftree->Gen_Ds_PT = Ds_GP.pt();
        ftree->Gen_Ds_PX = Ds_GP.px();
        ftree->Gen_Ds_PY = Ds_GP.py();
        ftree->Gen_Ds_PZ = Ds_GP.pz();
        TVector3 Gen_Ds_P3(Ds_GP.px(), Ds_GP.py(), Ds_GP.pz());

        ftree->Gen_mu_ETA = mu_GP.eta();
        ftree->Gen_mu_PHI = mu_GP.phi();
        ftree->Gen_mu_ORIVX_X = mu_GP.vx();
        ftree->Gen_mu_ORIVX_Y = mu_GP.vy();
        ftree->Gen_mu_ORIVX_Z = mu_GP.vz();
        ftree->Gen_mu_P = mu_GP.p();
        ftree->Gen_mu_PT = mu_GP.pt();
        ftree->Gen_mu_PX = mu_GP.px();
        ftree->Gen_mu_PY = mu_GP.py();
        ftree->Gen_mu_PZ = mu_GP.pz();
        TVector3 Gen_mu_P3(mu_GP.px(), mu_GP.py(), mu_GP.pz());

        ftree->Gen_nu_ETA = nu_GP.eta();
        ftree->Gen_nu_PHI = nu_GP.phi();
        ftree->Gen_nu_ORIVX_X = nu_GP.vx();
        ftree->Gen_nu_ORIVX_Y = nu_GP.vy();
        ftree->Gen_nu_ORIVX_Z = nu_GP.vz();
        ftree->Gen_nu_P = nu_GP.p();
        ftree->Gen_nu_PT = nu_GP.pt();
        ftree->Gen_nu_PX = nu_GP.px();
        ftree->Gen_nu_PY = nu_GP.py();
        ftree->Gen_nu_PZ = nu_GP.pz();
        TVector3 Gen_nu_P3(nu_GP.px(), nu_GP.py(), nu_GP.pz());

        ftree->Gen_W_ETA = W_GP.eta();
        ftree->Gen_W_PHI = W_GP.phi();
        ftree->Gen_W_ORIVX_X = W_GP.vx();
        ftree->Gen_W_ORIVX_Y = W_GP.vy();
        ftree->Gen_W_ORIVX_Z = W_GP.vz();
        ftree->Gen_W_P = W_GP.p();
        ftree->Gen_W_PT = W_GP.pt();
        ftree->Gen_W_PX = W_GP.px();
        ftree->Gen_W_PY = W_GP.py();
        ftree->Gen_W_PZ = W_GP.pz();
        TVector3 Gen_W_P3(W_GP.px(), W_GP.py(), W_GP.pz());

        ftree->Gen_H_ETA = H_GP.eta();
        ftree->Gen_H_PHI = H_GP.phi();
        ftree->Gen_H_ORIVX_X = H_GP.vx();
        ftree->Gen_H_ORIVX_Y = H_GP.vy();
        ftree->Gen_H_ORIVX_Z = H_GP.vz();
        ftree->Gen_H_P = H_GP.p();
        ftree->Gen_H_PT = H_GP.pt();
        ftree->Gen_H_PX = H_GP.px();
        ftree->Gen_H_PY = H_GP.py();
        ftree->Gen_H_PZ = H_GP.pz();
        TVector3 Gen_H_P3(H_GP.px(), H_GP.py(), H_GP.pz());

        ftree->Gen_Kp_PP = Gen_Kp_P3.Pt(Gen_phi_P3);
        ftree->Gen_Kp_PL = Gen_Kp_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag();
        ftree->Gen_Km_PP = Gen_Km_P3.Pt(Gen_phi_P3);
        ftree->Gen_Km_PL = Gen_Km_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag();
        ftree->Gen_pi_PP = Gen_pi_P3.Pt(Gen_Ds_P3);
        ftree->Gen_pi_PL = Gen_pi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag();
        ftree->Gen_phi_PP = Gen_phi_P3.Pt(Gen_Ds_P3);
        ftree->Gen_phi_PL = Gen_phi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag();

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

        ftree->Gen_dxy_Kp_Km = sqrt(pow(ftree->Gen_Kp_ORIVX_X-ftree->Gen_Km_ORIVX_X,2) + pow(ftree->Gen_Kp_ORIVX_Y-ftree->Gen_Km_ORIVX_Y,2));
        ftree->Gen_dxy_Kp_pi = sqrt(pow(ftree->Gen_Kp_ORIVX_X-ftree->Gen_pi_ORIVX_X,2) + pow(ftree->Gen_Kp_ORIVX_Y-ftree->Gen_pi_ORIVX_Y,2));
        ftree->Gen_dxy_Km_pi = sqrt(pow(ftree->Gen_Km_ORIVX_X-ftree->Gen_pi_ORIVX_X,2) + pow(ftree->Gen_Km_ORIVX_Y-ftree->Gen_pi_ORIVX_Y,2));
        ftree->Gen_dz_Kp_Km = abs(ftree->Gen_Kp_ORIVX_Z-ftree->Gen_Km_ORIVX_Z);
        ftree->Gen_dz_Kp_pi = abs(ftree->Gen_Kp_ORIVX_Z-ftree->Gen_pi_ORIVX_Z);
        ftree->Gen_dz_Km_pi = abs(ftree->Gen_Km_ORIVX_Z-ftree->Gen_pi_ORIVX_Z);

        ftree->Gen_Fill_Vector();

        struct MatchInfo {int index; float dr;};
        std::vector<MatchInfo> MatchKp, MatchKm, Matchpi;

        for(size_t i=0; i<packedPF->size(); i++){
            const auto& pf = (*packedPF)[i];

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

        if( bestKp.dr < 0.03 && bestKm.dr < 0.03 && bestpi.dr < 0.03 && (*packedPF)[bestKp.index].hasTrackDetails() && (*packedPF)[bestKm.index].hasTrackDetails() && (*packedPF)[bestpi.index].hasTrackDetails() && (*packedPF)[bestKp.index].trackHighPurity() && (*packedPF)[bestKm.index].trackHighPurity() && (*packedPF)[bestpi.index].trackHighPurity() ){

            const auto& match_Kp_PF = (*packedPF)[bestKp.index];
            const auto& match_Km_PF = (*packedPF)[bestKm.index];
            const auto& match_pi_PF = (*packedPF)[bestpi.index];

            // Original Kp
            ftree->match_Kp_ETA = match_Kp_PF.eta();
            ftree->match_Kp_PHI = match_Kp_PF.phi();
            ftree->match_Kp_ORIVX_X = match_Kp_PF.vx();
            ftree->match_Kp_ORIVX_Y = match_Kp_PF.vy();
            ftree->match_Kp_ORIVX_Z = match_Kp_PF.vz();
            ftree->match_Kp_P = match_Kp_PF.p();
            ftree->match_Kp_PT = match_Kp_PF.pt();
            ftree->match_Kp_PX = match_Kp_PF.px();
            ftree->match_Kp_PY = match_Kp_PF.py();
            ftree->match_Kp_PZ = match_Kp_PF.pz();
            TLorentzVector match_Kp_P4;
            match_Kp_P4.SetXYZM(match_Kp_PF.px(), match_Kp_PF.py(), match_Kp_PF.pz(), Mass_K);

            // Original Km
            ftree->match_Km_ETA = match_Km_PF.eta();
            ftree->match_Km_PHI = match_Km_PF.phi();
            ftree->match_Km_ORIVX_X = match_Km_PF.vx();
            ftree->match_Km_ORIVX_Y = match_Km_PF.vy();
            ftree->match_Km_ORIVX_Z = match_Km_PF.vz();
            ftree->match_Km_P = match_Km_PF.p();
            ftree->match_Km_PT = match_Km_PF.pt();
            ftree->match_Km_PX = match_Km_PF.px();
            ftree->match_Km_PY = match_Km_PF.py();
            ftree->match_Km_PZ = match_Km_PF.pz();
            TLorentzVector match_Km_P4;
            match_Km_P4.SetXYZM(match_Km_PF.px(), match_Km_PF.py(), match_Km_PF.pz(), Mass_K);

            // Original pi
            ftree->match_pi_ETA = match_pi_PF.eta();
            ftree->match_pi_PHI = match_pi_PF.phi();
            ftree->match_pi_ORIVX_X = match_pi_PF.vx();
            ftree->match_pi_ORIVX_Y = match_pi_PF.vy();
            ftree->match_pi_ORIVX_Z = match_pi_PF.vz();
            ftree->match_pi_P = match_pi_PF.p();
            ftree->match_pi_PT = match_pi_PF.pt();
            ftree->match_pi_PX = match_pi_PF.px();
            ftree->match_pi_PY = match_pi_PF.py();
            ftree->match_pi_PZ = match_pi_PF.pz();
            TLorentzVector match_pi_P4;
            match_pi_P4.SetXYZM(match_pi_PF.px(), match_pi_PF.py(), match_pi_PF.pz(), Mass_pi);

            // Original phi
            TLorentzVector match_phi_P4 = match_Kp_P4 + match_Km_P4;
            ftree->match_phi_ETA = match_phi_P4.Eta();
            ftree->match_phi_PHI = match_phi_P4.Phi();
            ftree->match_phi_P = match_phi_P4.P();
            ftree->match_phi_PT = match_phi_P4.Pt();
            ftree->match_phi_PX = match_phi_P4.Px();
            ftree->match_phi_PY = match_phi_P4.Py();
            ftree->match_phi_PZ = match_phi_P4.Pz();
            ftree->match_phi_M = match_phi_P4.M();

            // Original Ds
            TLorentzVector match_Ds_P4 = match_Kp_P4 + match_Km_P4 + match_pi_P4;
            ftree->match_Ds_ETA = match_Ds_P4.Eta();
            ftree->match_Ds_PHI = match_Ds_P4.Phi();
            ftree->match_Ds_P = match_Ds_P4.P();
            ftree->match_Ds_PT = match_Ds_P4.Pt();
            ftree->match_Ds_PX = match_Ds_P4.Px();
            ftree->match_Ds_PY = match_Ds_P4.Py();
            ftree->match_Ds_PZ = match_Ds_P4.Pz();
            ftree->match_Ds_M = match_Ds_P4.M();

            // some dxy dz
            ftree->match_dxy_Kp_Km = sqrt(pow(ftree->match_Kp_ORIVX_X-ftree->match_Km_ORIVX_X,2) + pow(ftree->match_Kp_ORIVX_Y-ftree->match_Km_ORIVX_Y,2));
            ftree->match_dxy_Kp_pi = sqrt(pow(ftree->match_Kp_ORIVX_X-ftree->match_pi_ORIVX_X,2) + pow(ftree->match_Kp_ORIVX_Y-ftree->match_pi_ORIVX_Y,2));
            ftree->match_dxy_Km_pi = sqrt(pow(ftree->match_Km_ORIVX_X-ftree->match_pi_ORIVX_X,2) + pow(ftree->match_Km_ORIVX_Y-ftree->match_pi_ORIVX_Y,2));

            ftree->match_dz_Kp_Km = abs(ftree->match_Kp_ORIVX_Z-ftree->match_Km_ORIVX_Z);
            ftree->match_dz_Kp_pi = abs(ftree->match_Kp_ORIVX_Z-ftree->match_pi_ORIVX_Z);
            ftree->match_dz_Km_pi = abs(ftree->match_Km_ORIVX_Z-ftree->match_pi_ORIVX_Z);

            // Original PP PL
            ftree->match_Kp_PP = match_Kp_P4.Vect().Pt(match_phi_P4.Vect());
            ftree->match_Kp_PL = match_Kp_P4.Vect().Dot(match_phi_P4.Vect())/match_phi_P4.P();
            ftree->match_Km_PP = match_Km_P4.Vect().Pt(match_phi_P4.Vect());
            ftree->match_Km_PL = match_Km_P4.Vect().Dot(match_phi_P4.Vect())/match_phi_P4.P();

            ftree->match_phi_PP = match_phi_P4.Vect().Pt(match_Ds_P4.Vect());
            ftree->match_phi_PL = match_phi_P4.Vect().Dot(match_Ds_P4.Vect())/match_Ds_P4.P();
            ftree->match_pi_PP = match_pi_P4.Vect().Pt(match_Ds_P4.Vect());
            ftree->match_pi_PL = match_pi_P4.Vect().Dot(match_Ds_P4.Vect())/match_Ds_P4.P();

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
                (*theB).build(match_Kp_PF.pseudoTrack()),
                (*theB).build(match_Km_PF.pseudoTrack())
            };

            TransientVertex match_phi_Vertex = KVFitter.vertex(match_phi_Tracks);

            if( match_phi_Vertex.isValid() && match_phi_Vertex.hasRefittedTracks() ){

                // phi fit variables 
                ftree->match_phiFit_CHI2 = match_phi_Vertex.totalChiSquared();
                ftree->match_phiFit_NDOF = match_phi_Vertex.degreesOfFreedom();
                ftree->match_phiFit_CHI2NDOF = match_phi_Vertex.normalisedChiSquared();
                ftree->match_phiFit_ENDVX_X = match_phi_Vertex.position().x();
                ftree->match_phiFit_ENDVX_Y = match_phi_Vertex.position().y();
                ftree->match_phiFit_ENDVX_Z = match_phi_Vertex.position().z();
                ftree->match_phiFit_ENDVX_XERR = std::sqrt(match_phi_Vertex.positionError().cxx());
                ftree->match_phiFit_ENDVX_YERR = std::sqrt(match_phi_Vertex.positionError().cyy());
                ftree->match_phiFit_ENDVX_ZERR = std::sqrt(match_phi_Vertex.positionError().czz());

                // some dxy dz
                ftree->match_dxy_Kp_phi = sqrt(pow(ftree->match_Kp_ORIVX_X-ftree->match_phiFit_ENDVX_X,2) + pow(ftree->match_Kp_ORIVX_Y-ftree->match_phiFit_ENDVX_Y,2));
                ftree->match_dxy_Km_phi = sqrt(pow(ftree->match_Km_ORIVX_X-ftree->match_phiFit_ENDVX_X,2) + pow(ftree->match_Km_ORIVX_Y-ftree->match_phiFit_ENDVX_Y,2));
                ftree->match_dxy_pi_phi = sqrt(pow(ftree->match_pi_ORIVX_X-ftree->match_phiFit_ENDVX_X,2) + pow(ftree->match_pi_ORIVX_Y-ftree->match_phiFit_ENDVX_Y,2));

                ftree->match_dz_Kp_phi = abs(ftree->match_Kp_ORIVX_Z-ftree->match_phiFit_ENDVX_Z);
                ftree->match_dz_Km_phi = abs(ftree->match_Km_ORIVX_Z-ftree->match_phiFit_ENDVX_Z);
                ftree->match_dz_pi_phi = abs(ftree->match_pi_ORIVX_Z-ftree->match_phiFit_ENDVX_Z);

                std::vector<reco::TransientTrack> match_phiFit_Tracks = match_phi_Vertex.refittedTracks();

                // phi fit Kp
                TLorentzVector match_phiFit_Kp_P4;
                match_phiFit_Kp_P4.SetXYZM(match_phiFit_Tracks[0].track().px(), match_phiFit_Tracks[0].track().py(), match_phiFit_Tracks[0].track().pz(), Mass_K);
                ftree->match_phiFit_Kp_ETA = match_phiFit_Kp_P4.Eta();
                ftree->match_phiFit_Kp_PHI = match_phiFit_Kp_P4.Phi();
                ftree->match_phiFit_Kp_P = match_phiFit_Kp_P4.P();
                ftree->match_phiFit_Kp_PT = match_phiFit_Kp_P4.Pt();
                ftree->match_phiFit_Kp_PX = match_phiFit_Kp_P4.Px();
                ftree->match_phiFit_Kp_PY = match_phiFit_Kp_P4.Py();
                ftree->match_phiFit_Kp_PZ = match_phiFit_Kp_P4.Pz();

                // phi fit Km
                TLorentzVector match_phiFit_Km_P4;
                match_phiFit_Km_P4.SetXYZM(match_phiFit_Tracks[1].track().px(), match_phiFit_Tracks[1].track().py(), match_phiFit_Tracks[1].track().pz(), Mass_K);
                ftree->match_phiFit_Km_ETA = match_phiFit_Km_P4.Eta();
                ftree->match_phiFit_Km_PHI = match_phiFit_Km_P4.Phi();
                ftree->match_phiFit_Km_P = match_phiFit_Km_P4.P();
                ftree->match_phiFit_Km_PT = match_phiFit_Km_P4.Pt();
                ftree->match_phiFit_Km_PX = match_phiFit_Km_P4.Px();
                ftree->match_phiFit_Km_PY = match_phiFit_Km_P4.Py();
                ftree->match_phiFit_Km_PZ = match_phiFit_Km_P4.Pz();

                // phi fit pi = orignal pi
                TLorentzVector match_phiFit_pi_P4 = match_pi_P4;
                ftree->match_phiFit_pi_ETA = match_phiFit_pi_P4.Eta();
                ftree->match_phiFit_pi_PHI = match_phiFit_pi_P4.Phi();
                ftree->match_phiFit_pi_P = match_phiFit_pi_P4.P();
                ftree->match_phiFit_pi_PT = match_phiFit_pi_P4.Pt();
                ftree->match_phiFit_pi_PX = match_phiFit_pi_P4.Px();
                ftree->match_phiFit_pi_PY = match_phiFit_pi_P4.Py();
                ftree->match_phiFit_pi_PZ = match_phiFit_pi_P4.Pz();

                // phi fit phi
                TLorentzVector match_phiFit_phi_P4 = match_phiFit_Kp_P4 + match_phiFit_Km_P4;
                ftree->match_phiFit_phi_ETA = match_phiFit_phi_P4.Eta();
                ftree->match_phiFit_phi_PHI = match_phiFit_phi_P4.Phi();
                ftree->match_phiFit_phi_P = match_phiFit_phi_P4.P();
                ftree->match_phiFit_phi_PT = match_phiFit_phi_P4.Pt();
                ftree->match_phiFit_phi_PX = match_phiFit_phi_P4.Px();
                ftree->match_phiFit_phi_PY = match_phiFit_phi_P4.Py();
                ftree->match_phiFit_phi_PZ = match_phiFit_phi_P4.Pz();
                ftree->match_phiFit_phi_M = match_phiFit_phi_P4.M();

                // phi fit Ds
                TLorentzVector match_phiFit_Ds_P4 = match_phiFit_Kp_P4 + match_phiFit_Km_P4 + match_phiFit_pi_P4;
                ftree->match_phiFit_Ds_ETA = match_phiFit_Ds_P4.Eta();
                ftree->match_phiFit_Ds_PHI = match_phiFit_Ds_P4.Phi();
                ftree->match_phiFit_Ds_P = match_phiFit_Ds_P4.P();
                ftree->match_phiFit_Ds_PT = match_phiFit_Ds_P4.Pt();
                ftree->match_phiFit_Ds_PX = match_phiFit_Ds_P4.Px();
                ftree->match_phiFit_Ds_PY = match_phiFit_Ds_P4.Py();
                ftree->match_phiFit_Ds_PZ = match_phiFit_Ds_P4.Pz();
                ftree->match_phiFit_Ds_M = match_phiFit_Ds_P4.M();

                // phi fit PP PL 
                ftree->match_phiFit_Kp_PP = match_phiFit_Kp_P4.Vect().Pt(match_phiFit_phi_P4.Vect());
                ftree->match_phiFit_Kp_PL = match_phiFit_Kp_P4.Vect().Dot(match_phiFit_phi_P4.Vect())/match_phiFit_phi_P4.P();
                ftree->match_phiFit_Km_PP = match_phiFit_Km_P4.Vect().Pt(match_phiFit_phi_P4.Vect());
                ftree->match_phiFit_Km_PL = match_phiFit_Km_P4.Vect().Dot(match_phiFit_phi_P4.Vect())/match_phiFit_phi_P4.P();

                ftree->match_phiFit_phi_PP = match_phiFit_phi_P4.Vect().Pt(match_phiFit_Ds_P4.Vect());
                ftree->match_phiFit_phi_PL = match_phiFit_phi_P4.Vect().Dot(match_phiFit_Ds_P4.Vect())/match_phiFit_Ds_P4.P();
                ftree->match_phiFit_pi_PP = match_phiFit_pi_P4.Vect().Pt(match_phiFit_Ds_P4.Vect());
                ftree->match_phiFit_pi_PL = match_phiFit_pi_P4.Vect().Dot(match_phiFit_Ds_P4.Vect())/match_phiFit_Ds_P4.P();

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

                // Ds fit
                std::vector<reco::TransientTrack> match_Ds_Tracks = {
                    match_phiFit_Tracks[0],
                    match_phiFit_Tracks[1],
                    (*theB).build(match_pi_PF.pseudoTrack())
                };

                TransientVertex match_Ds_Vertex = KVFitter.vertex(match_Ds_Tracks);

                if( match_Ds_Vertex.isValid() && match_Ds_Vertex.hasRefittedTracks() ){

                    // Ds fit variables 
                    ftree->match_DsFit_CHI2 = match_Ds_Vertex.totalChiSquared();
                    ftree->match_DsFit_NDOF = match_Ds_Vertex.degreesOfFreedom();
                    ftree->match_DsFit_CHI2NDOF = match_Ds_Vertex.normalisedChiSquared();
                    ftree->match_DsFit_ENDVX_X = match_Ds_Vertex.position().x();
                    ftree->match_DsFit_ENDVX_Y = match_Ds_Vertex.position().y();
                    ftree->match_DsFit_ENDVX_Z = match_Ds_Vertex.position().z();
                    ftree->match_DsFit_ENDVX_XERR = std::sqrt(match_Ds_Vertex.positionError().cxx());
                    ftree->match_DsFit_ENDVX_YERR = std::sqrt(match_Ds_Vertex.positionError().cyy());
                    ftree->match_DsFit_ENDVX_ZERR = std::sqrt(match_Ds_Vertex.positionError().czz());

                    // some dxy dz
                    ftree->match_dxy_Kp_Ds = sqrt(pow(ftree->match_Kp_ORIVX_X-ftree->match_DsFit_ENDVX_X,2) + pow(ftree->match_Kp_ORIVX_Y-ftree->match_DsFit_ENDVX_Y,2));
                    ftree->match_dxy_Km_Ds = sqrt(pow(ftree->match_Km_ORIVX_X-ftree->match_DsFit_ENDVX_X,2) + pow(ftree->match_Km_ORIVX_Y-ftree->match_DsFit_ENDVX_Y,2));
                    ftree->match_dxy_pi_Ds = sqrt(pow(ftree->match_pi_ORIVX_X-ftree->match_DsFit_ENDVX_X,2) + pow(ftree->match_pi_ORIVX_Y-ftree->match_DsFit_ENDVX_Y,2));
                    ftree->match_dxy_phi_Ds = sqrt(pow(ftree->match_phiFit_ENDVX_X-ftree->match_DsFit_ENDVX_X,2) + pow(ftree->match_phiFit_ENDVX_Y-ftree->match_DsFit_ENDVX_Y,2));

                    ftree->match_dz_Kp_Ds = abs(ftree->match_Kp_ORIVX_Z-ftree->match_DsFit_ENDVX_Z);
                    ftree->match_dz_Km_Ds = abs(ftree->match_Km_ORIVX_Z-ftree->match_DsFit_ENDVX_Z);
                    ftree->match_dz_pi_Ds = abs(ftree->match_pi_ORIVX_Z-ftree->match_DsFit_ENDVX_Z);
                    ftree->match_dz_phi_Ds = abs(ftree->match_phiFit_ENDVX_Z-ftree->match_DsFit_ENDVX_Z);

                    std::vector<reco::TransientTrack> match_DsFit_Tracks = match_Ds_Vertex.refittedTracks();

                    // Ds fit Kp
                    TLorentzVector match_DsFit_Kp_P4;
                    match_DsFit_Kp_P4.SetXYZM(match_DsFit_Tracks[0].track().px(), match_DsFit_Tracks[0].track().py(), match_DsFit_Tracks[0].track().pz(), Mass_K);
                    ftree->match_DsFit_Kp_ETA = match_DsFit_Kp_P4.Eta();
                    ftree->match_DsFit_Kp_PHI = match_DsFit_Kp_P4.Phi();
                    ftree->match_DsFit_Kp_P = match_DsFit_Kp_P4.P();
                    ftree->match_DsFit_Kp_PT = match_DsFit_Kp_P4.Pt();
                    ftree->match_DsFit_Kp_PX = match_DsFit_Kp_P4.Px();
                    ftree->match_DsFit_Kp_PY = match_DsFit_Kp_P4.Py();
                    ftree->match_DsFit_Kp_PZ = match_DsFit_Kp_P4.Pz();

                    // Ds fit Km
                    TLorentzVector match_DsFit_Km_P4;
                    match_DsFit_Km_P4.SetXYZM(match_DsFit_Tracks[1].track().px(), match_DsFit_Tracks[1].track().py(), match_DsFit_Tracks[1].track().pz(), Mass_K);
                    ftree->match_DsFit_Km_ETA = match_DsFit_Km_P4.Eta();
                    ftree->match_DsFit_Km_PHI = match_DsFit_Km_P4.Phi();
                    ftree->match_DsFit_Km_P = match_DsFit_Km_P4.P();
                    ftree->match_DsFit_Km_PT = match_DsFit_Km_P4.Pt();
                    ftree->match_DsFit_Km_PX = match_DsFit_Km_P4.Px();
                    ftree->match_DsFit_Km_PY = match_DsFit_Km_P4.Py();
                    ftree->match_DsFit_Km_PZ = match_DsFit_Km_P4.Pz();

                    // Ds fit pi
                    TLorentzVector match_DsFit_pi_P4;
                    match_DsFit_pi_P4.SetXYZM(match_DsFit_Tracks[2].track().px(), match_DsFit_Tracks[2].track().py(), match_DsFit_Tracks[2].track().pz(), Mass_pi);
                    ftree->match_DsFit_pi_ETA = match_DsFit_pi_P4.Eta();
                    ftree->match_DsFit_pi_PHI = match_DsFit_pi_P4.Phi();
                    ftree->match_DsFit_pi_P = match_DsFit_pi_P4.P();
                    ftree->match_DsFit_pi_PT = match_DsFit_pi_P4.Pt();
                    ftree->match_DsFit_pi_PX = match_DsFit_pi_P4.Px();
                    ftree->match_DsFit_pi_PY = match_DsFit_pi_P4.Py();
                    ftree->match_DsFit_pi_PZ = match_DsFit_pi_P4.Pz();

                    // Ds fit phi
                    TLorentzVector match_DsFit_phi_P4 = match_DsFit_Kp_P4 + match_DsFit_Km_P4;
                    ftree->match_DsFit_phi_ETA = match_DsFit_phi_P4.Eta();
                    ftree->match_DsFit_phi_PHI = match_DsFit_phi_P4.Phi();
                    ftree->match_DsFit_phi_P = match_DsFit_phi_P4.P();
                    ftree->match_DsFit_phi_PT = match_DsFit_phi_P4.Pt();
                    ftree->match_DsFit_phi_PX = match_DsFit_phi_P4.Px();
                    ftree->match_DsFit_phi_PY = match_DsFit_phi_P4.Py();
                    ftree->match_DsFit_phi_PZ = match_DsFit_phi_P4.Pz();
                    ftree->match_DsFit_phi_M = match_DsFit_phi_P4.M();

                    // Ds fit Ds
                    TLorentzVector match_DsFit_Ds_P4 = match_DsFit_Kp_P4 + match_DsFit_Km_P4 + match_DsFit_pi_P4;
                    ftree->match_DsFit_Ds_ETA = match_DsFit_Ds_P4.Eta();
                    ftree->match_DsFit_Ds_PHI = match_DsFit_Ds_P4.Phi();
                    ftree->match_DsFit_Ds_P = match_DsFit_Ds_P4.P();
                    ftree->match_DsFit_Ds_PT = match_DsFit_Ds_P4.Pt();
                    ftree->match_DsFit_Ds_PX = match_DsFit_Ds_P4.Px();
                    ftree->match_DsFit_Ds_PY = match_DsFit_Ds_P4.Py();
                    ftree->match_DsFit_Ds_PZ = match_DsFit_Ds_P4.Pz();
                    ftree->match_DsFit_Ds_M = match_DsFit_Ds_P4.M();

                    TLorentzVector match_DsFit_Mconstraint_phi_P4;
                    match_DsFit_Mconstraint_phi_P4.SetXYZM(match_DsFit_phi_P4.Px(), match_DsFit_phi_P4.Py(), match_DsFit_phi_P4.Pz(), Mass_phi);
                    TLorentzVector match_DsFit_Mconstraint_Ds_P4 = match_DsFit_Mconstraint_phi_P4 + match_DsFit_pi_P4;
                    ftree->match_DsFit_Mconstraint_Ds_M = match_DsFit_Mconstraint_Ds_P4.M();

                    // Ds fit PP PL 
                    ftree->match_DsFit_Kp_PP = match_DsFit_Kp_P4.Vect().Pt(match_DsFit_phi_P4.Vect());
                    ftree->match_DsFit_Kp_PL = match_DsFit_Kp_P4.Vect().Dot(match_DsFit_phi_P4.Vect())/match_DsFit_phi_P4.P();
                    ftree->match_DsFit_Km_PP = match_DsFit_Km_P4.Vect().Pt(match_DsFit_phi_P4.Vect());
                    ftree->match_DsFit_Km_PL = match_DsFit_Km_P4.Vect().Dot(match_DsFit_phi_P4.Vect())/match_DsFit_phi_P4.P();

                    ftree->match_DsFit_phi_PP = match_DsFit_phi_P4.Vect().Pt(match_DsFit_Ds_P4.Vect());
                    ftree->match_DsFit_phi_PL = match_DsFit_phi_P4.Vect().Dot(match_DsFit_Ds_P4.Vect())/match_DsFit_Ds_P4.P();
                    ftree->match_DsFit_pi_PP = match_DsFit_pi_P4.Vect().Pt(match_DsFit_Ds_P4.Vect());
                    ftree->match_DsFit_pi_PL = match_DsFit_pi_P4.Vect().Dot(match_DsFit_Ds_P4.Vect())/match_DsFit_Ds_P4.P();

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

                    ftree->Match_Fill_Vector();


                    std::vector<reco::TransientTrack> match_PV_ttracks;

                    for(size_t i=0; i<packedPF->size(); i++){

                        if(int(i) ==  ftree->match_Kp_idx) continue;
                        if(int(i) ==  ftree->match_Km_idx) continue;
                        if(int(i) ==  ftree->match_pi_idx) continue;
                        
                        const auto& pf = (*packedPF)[i];

                        if( !(pf.trackHighPurity()) ) continue;
                        if( !(pf.hasTrackDetails()) ) continue;

                        if(pf.pt() > 1){
                            match_PV_ttracks.push_back( (*theB).build(pf.pseudoTrack()) );
                        }
                    }

                    if(match_PV_ttracks.size() > 2){
                        TransientVertex match_PV_tvertex = AVFitter.vertex(match_PV_ttracks);
                        if( match_PV_tvertex.isValid() ){ 
                            ftree->match_PV_CHI2 = match_PV_tvertex.totalChiSquared();
                            ftree->match_PV_NDOF = match_PV_tvertex.degreesOfFreedom();
                            ftree->match_PV_CHI2NDOF = match_PV_tvertex.normalisedChiSquared();
                            ftree->match_PV_X = match_PV_tvertex.position().x();
                            ftree->match_PV_Y = match_PV_tvertex.position().y();
                            ftree->match_PV_Z = match_PV_tvertex.position().z();
                            ftree->match_PV_XERR = std::sqrt(match_PV_tvertex.positionError().cxx());
                            ftree->match_PV_YERR = std::sqrt(match_PV_tvertex.positionError().cyy());
                            ftree->match_PV_ZERR = std::sqrt(match_PV_tvertex.positionError().czz());
                        }
                    }
                    
                    ftree->Match_Fill_PV();
                }
            }
        }
    } 

    idx_Kp_vec.clear();
    idx_Km_vec.clear();
    idx_pi_vec.clear();

    for(size_t i=0; i<packedPF->size(); i++){
        const auto& pf = (*packedPF)[i];

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

        const auto& Kp_PF = (*packedPF)[idx_Kp_vec[i]];

        ftree->Kp_Reset();

        ftree->Kp_ETA = Kp_PF.eta();
        ftree->Kp_PHI = Kp_PF.phi();
        ftree->Kp_ORIVX_X = Kp_PF.vx();
        ftree->Kp_ORIVX_Y = Kp_PF.vy();
        ftree->Kp_ORIVX_Z = Kp_PF.vz();
        ftree->Kp_P = Kp_PF.p();
        ftree->Kp_PT = Kp_PF.pt();
        ftree->Kp_PX = Kp_PF.px();
        ftree->Kp_PY = Kp_PF.py();
        ftree->Kp_PZ = Kp_PF.pz();
        TLorentzVector Kp_P4;
        Kp_P4.SetXYZM(Kp_PF.px(), Kp_PF.py(), Kp_PF.pz(), Mass_K);

        for(size_t j=0; j<idx_Km_vec.size(); j++){
            const auto& Km_PF = (*packedPF)[idx_Km_vec[j]];

            ftree->Km_Reset();
            ftree->Km_ETA = Km_PF.eta();
            ftree->Km_PHI = Km_PF.phi();
            ftree->Km_ORIVX_X = Km_PF.vx();
            ftree->Km_ORIVX_Y = Km_PF.vy();
            ftree->Km_ORIVX_Z = Km_PF.vz();
            ftree->Km_P = Km_PF.p();
            ftree->Km_PT = Km_PF.pt();
            ftree->Km_PX = Km_PF.px();
            ftree->Km_PY = Km_PF.py();
            ftree->Km_PZ = Km_PF.pz();
            TLorentzVector Km_P4;
            Km_P4.SetXYZM(Km_PF.px(), Km_PF.py(), Km_PF.pz(), Mass_K);

            TLorentzVector phi_P4 = Kp_P4 + Km_P4;
            ftree->phi_ETA = phi_P4.Eta();
            ftree->phi_PHI = phi_P4.Phi();
            ftree->phi_P = phi_P4.P();
            ftree->phi_PT = phi_P4.Pt();
            ftree->phi_PX = phi_P4.Px();
            ftree->phi_PY = phi_P4.Py();
            ftree->phi_PZ = phi_P4.Pz();
            ftree->phi_M = phi_P4.M();

            ftree->Kp_PP = Kp_P4.Vect().Pt(phi_P4.Vect());
            ftree->Kp_PL = Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();
            ftree->Km_PP = Km_P4.Vect().Pt(phi_P4.Vect());
            ftree->Km_PL = Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P();

            ftree->dR_Kp_Km = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Km_P4.Eta(), Km_P4.Phi());
            ftree->dR_Kp_phi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());
            ftree->dR_Km_phi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());

            ftree->dxy_Kp_Km = sqrt(pow(ftree->Kp_ORIVX_X-ftree->Km_ORIVX_X,2) + pow(ftree->Kp_ORIVX_Y-ftree->Km_ORIVX_Y,2));
            ftree->dz_Kp_Km = abs(ftree->Kp_ORIVX_Z-ftree->Km_ORIVX_Z);

            if( abs(ftree->Kp_PT - ftree->Km_PT) > 20) continue;
            /* if( ftree->dR_Kp_Km > 0.15 ) continue; */
            if( ftree->dR_Kp_Km > 0.1 ) continue;
            /* if( ftree->dxy_Kp_Km > 0.2) continue; */
            /* if( ftree->dz_Kp_Km > 0.2) continue; */

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*theB).build(Kp_PF.pseudoTrack()),
                (*theB).build(Km_PF.pseudoTrack())
            };
            TransientVertex phi_Vertex = KVFitter.vertex(phi_Tracks);

            if( !(phi_Vertex.isValid()) ) continue;
            if( !(phi_Vertex.hasRefittedTracks()) ) continue;

            std::vector<reco::TransientTrack> phiFit_Tracks = phi_Vertex.refittedTracks();

            ftree->phiFit_CHI2 = phi_Vertex.totalChiSquared();
            ftree->phiFit_NDOF = phi_Vertex.degreesOfFreedom();
            ftree->phiFit_CHI2NDOF = phi_Vertex.normalisedChiSquared();
            ftree->phiFit_ENDVX_X = phi_Vertex.position().x();
            ftree->phiFit_ENDVX_Y = phi_Vertex.position().y();
            ftree->phiFit_ENDVX_Z = phi_Vertex.position().z();
            ftree->phiFit_ENDVX_XERR = std::sqrt(phi_Vertex.positionError().cxx());
            ftree->phiFit_ENDVX_YERR = std::sqrt(phi_Vertex.positionError().cyy());
            ftree->phiFit_ENDVX_ZERR = std::sqrt(phi_Vertex.positionError().czz());

            TLorentzVector phiFit_Kp_P4;
            phiFit_Kp_P4.SetXYZM(phiFit_Tracks[0].track().px(), phiFit_Tracks[0].track().py(), phiFit_Tracks[0].track().pz(), Mass_K);
            ftree->phiFit_Kp_ETA = phiFit_Kp_P4.Eta();
            ftree->phiFit_Kp_PHI = phiFit_Kp_P4.Phi();
            ftree->phiFit_Kp_P = phiFit_Kp_P4.P();
            ftree->phiFit_Kp_PT = phiFit_Kp_P4.Pt();
            ftree->phiFit_Kp_PX = phiFit_Kp_P4.Px();
            ftree->phiFit_Kp_PY = phiFit_Kp_P4.Py();
            ftree->phiFit_Kp_PZ = phiFit_Kp_P4.Pz();

            TLorentzVector phiFit_Km_P4;
            phiFit_Km_P4.SetXYZM(phiFit_Tracks[1].track().px(), phiFit_Tracks[1].track().py(), phiFit_Tracks[1].track().pz(), Mass_K);
            ftree->phiFit_Km_ETA = phiFit_Km_P4.Eta();
            ftree->phiFit_Km_PHI = phiFit_Km_P4.Phi();
            ftree->phiFit_Km_P = phiFit_Km_P4.P();
            ftree->phiFit_Km_PT = phiFit_Km_P4.Pt();
            ftree->phiFit_Km_PX = phiFit_Km_P4.Px();
            ftree->phiFit_Km_PY = phiFit_Km_P4.Py();
            ftree->phiFit_Km_PZ = phiFit_Km_P4.Pz();

            TLorentzVector phiFit_phi_P4 = phiFit_Kp_P4 + phiFit_Km_P4;
            ftree->phiFit_phi_ETA = phiFit_phi_P4.Eta();
            ftree->phiFit_phi_PHI = phiFit_phi_P4.Phi();
            ftree->phiFit_phi_P = phiFit_phi_P4.P();
            ftree->phiFit_phi_PT = phiFit_phi_P4.Pt();
            ftree->phiFit_phi_PX = phiFit_phi_P4.Px();
            ftree->phiFit_phi_PY = phiFit_phi_P4.Py();
            ftree->phiFit_phi_PZ = phiFit_phi_P4.Pz();
            ftree->phiFit_phi_M = phiFit_phi_P4.M();

            ftree->phiFit_Kp_PP = phiFit_Kp_P4.Vect().Pt(phiFit_phi_P4.Vect());
            ftree->phiFit_Kp_PL = phiFit_Kp_P4.Vect().Dot(phiFit_phi_P4.Vect())/phiFit_phi_P4.P();
            ftree->phiFit_Km_PP = phiFit_Km_P4.Vect().Pt(phiFit_phi_P4.Vect());
            ftree->phiFit_Km_PL = phiFit_Km_P4.Vect().Dot(phiFit_phi_P4.Vect())/phiFit_phi_P4.P();

            ftree->phiFit_dR_Kp_Km = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi());
            ftree->phiFit_dR_Kp_phi = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());
            ftree->phiFit_dR_Km_phi = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());

            ftree->dxy_Kp_phi = sqrt(pow(ftree->Kp_ORIVX_X-ftree->phiFit_ENDVX_X,2) + pow(ftree->Kp_ORIVX_Y-ftree->phiFit_ENDVX_Y,2));
            ftree->dxy_Km_phi = sqrt(pow(ftree->Km_ORIVX_X-ftree->phiFit_ENDVX_X,2) + pow(ftree->Km_ORIVX_Y-ftree->phiFit_ENDVX_Y,2));

            ftree->dz_Kp_phi = abs(ftree->Kp_ORIVX_Z-ftree->phiFit_ENDVX_Z);
            ftree->dz_Km_phi = abs(ftree->Km_ORIVX_Z-ftree->phiFit_ENDVX_Z);

            if( ftree->phiFit_CHI2NDOF < 0 ) continue;
            if( ftree->phiFit_CHI2NDOF > 10 ) continue;
            ftree->num_reco_phi++;

            for(size_t k=0; k<idx_pi_vec.size(); k++){

                if( idx_pi_vec[k] == idx_Kp_vec[i] ) continue;

                const auto& pi_PF = (*packedPF)[idx_pi_vec[k]];
                ftree->pi_Reset();

                ftree->pi_ETA = pi_PF.eta();
                ftree->pi_PHI = pi_PF.phi();
                ftree->pi_ORIVX_X = pi_PF.vx();
                ftree->pi_ORIVX_Y = pi_PF.vy();
                ftree->pi_ORIVX_Z = pi_PF.vz();
                ftree->pi_P = pi_PF.p();
                ftree->pi_PT = pi_PF.pt();
                ftree->pi_PX = pi_PF.px();
                ftree->pi_PY = pi_PF.py();
                ftree->pi_PZ = pi_PF.pz();
                TLorentzVector pi_P4;
                pi_P4.SetXYZM(pi_PF.px(), pi_PF.py(), pi_PF.pz(), Mass_pi);

                ftree->dR_Kp_pi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), pi_P4.Eta(), pi_P4.Phi());
                ftree->dR_Km_pi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), pi_P4.Eta(), pi_P4.Phi());
                ftree->dR_pi_phi = reco::deltaR(pi_P4.Eta(), pi_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());

                TLorentzVector Ds_P4 = Kp_P4 + Km_P4 + pi_P4;
                ftree->Ds_ETA = Ds_P4.Eta();
                ftree->Ds_PHI = Ds_P4.Phi();
                ftree->Ds_P = Ds_P4.P();
                ftree->Ds_PT = Ds_P4.Pt();
                ftree->Ds_PX = Ds_P4.Px();
                ftree->Ds_PY = Ds_P4.Py();
                ftree->Ds_PZ = Ds_P4.Pz();
                ftree->Ds_M = Ds_P4.M();

                ftree->dR_Kp_Ds = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_Km_Ds = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_pi_Ds = reco::deltaR(pi_P4.Eta(), pi_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
                ftree->dR_phi_Ds = reco::deltaR(phi_P4.Eta(), phi_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());

                ftree->pi_PP = pi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->pi_PL = pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();
                ftree->phi_PP = phi_P4.Vect().Pt(Ds_P4.Vect());
                ftree->phi_PL = phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P();

                TLorentzVector phiFit_pi_P4 = pi_P4;
                ftree->phiFit_pi_ETA = phiFit_pi_P4.Eta();
                ftree->phiFit_pi_PHI = phiFit_pi_P4.Phi();
                ftree->phiFit_pi_P = phiFit_pi_P4.P();
                ftree->phiFit_pi_PT = phiFit_pi_P4.Pt();
                ftree->phiFit_pi_PX = phiFit_pi_P4.Px();
                ftree->phiFit_pi_PY = phiFit_pi_P4.Py();
                ftree->phiFit_pi_PZ = phiFit_pi_P4.Pz();

                ftree->phiFit_dR_Kp_pi = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_pi_P4.Eta(), phiFit_pi_P4.Phi());
                ftree->phiFit_dR_Km_pi = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_pi_P4.Eta(), phiFit_pi_P4.Phi());
                ftree->phiFit_dR_pi_phi = reco::deltaR(phiFit_pi_P4.Eta(), phiFit_pi_P4.Phi(), phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi());

                TLorentzVector phiFit_Ds_P4 = phiFit_Kp_P4 + phiFit_Km_P4 + phiFit_pi_P4;
                ftree->phiFit_Ds_ETA = phiFit_Ds_P4.Eta();
                ftree->phiFit_Ds_PHI = phiFit_Ds_P4.Phi();
                ftree->phiFit_Ds_P = phiFit_Ds_P4.P();
                ftree->phiFit_Ds_PT = phiFit_Ds_P4.Pt();
                ftree->phiFit_Ds_PX = phiFit_Ds_P4.Px();
                ftree->phiFit_Ds_PY = phiFit_Ds_P4.Py();
                ftree->phiFit_Ds_PZ = phiFit_Ds_P4.Pz();
                ftree->phiFit_Ds_M = phiFit_Ds_P4.M();

                ftree->phiFit_pi_PP = phiFit_pi_P4.Vect().Pt(phiFit_Ds_P4.Vect());
                ftree->phiFit_pi_PL = phiFit_pi_P4.Vect().Dot(phiFit_Ds_P4.Vect())/phiFit_Ds_P4.P();
                ftree->phiFit_phi_PP = phiFit_phi_P4.Vect().Pt(phiFit_Ds_P4.Vect());
                ftree->phiFit_phi_PL = phiFit_phi_P4.Vect().Dot(phiFit_Ds_P4.Vect())/phiFit_Ds_P4.P();

                ftree->phiFit_dR_Kp_Ds = reco::deltaR(phiFit_Kp_P4.Eta(), phiFit_Kp_P4.Phi(), phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_Km_Ds = reco::deltaR(phiFit_Km_P4.Eta(), phiFit_Km_P4.Phi(), phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_pi_Ds = reco::deltaR(phiFit_pi_P4.Eta(), phiFit_pi_P4.Phi(), phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());
                ftree->phiFit_dR_phi_Ds = reco::deltaR(phiFit_phi_P4.Eta(), phiFit_phi_P4.Phi(), phiFit_Ds_P4.Eta(), phiFit_Ds_P4.Phi());

                /* if( ftree->dR_Kp_pi > 0.6 ) continue; */
                /* if( ftree->dR_Km_pi > 0.6 ) continue; */
                if( ftree->dR_Kp_pi > 0.4 ) continue;
                if( ftree->dR_Km_pi > 0.4 ) continue;
                if( ftree->phiFit_dR_Kp_phi > 0.05 ) continue;
                if( ftree->phiFit_dR_Km_phi > 0.05 ) continue;
                if( ftree->phiFit_dR_pi_phi > 0.4 ) continue;
                if( ftree->phiFit_phi_M < 0.99 ) continue;
                if( ftree->phiFit_phi_M > 1.05 ) continue;

                ftree->dxy_Kp_pi = sqrt(pow(ftree->Kp_ORIVX_X-ftree->pi_ORIVX_X,2) + pow(ftree->Kp_ORIVX_Y-ftree->pi_ORIVX_Y,2));
                ftree->dxy_Km_pi = sqrt(pow(ftree->Km_ORIVX_X-ftree->pi_ORIVX_X,2) + pow(ftree->Km_ORIVX_Y-ftree->pi_ORIVX_Y,2));
                ftree->dxy_pi_phi = sqrt(pow(ftree->pi_ORIVX_X-ftree->phiFit_ENDVX_X,2) + pow(ftree->pi_ORIVX_Y-ftree->phiFit_ENDVX_Y,2));

                ftree->dz_Kp_pi = abs(ftree->Kp_ORIVX_Z-ftree->pi_ORIVX_Z);
                ftree->dz_Km_pi = abs(ftree->Km_ORIVX_Z-ftree->pi_ORIVX_Z);
                ftree->dz_pi_phi = abs(ftree->pi_ORIVX_Z-ftree->phiFit_ENDVX_Z);

                /* if( ftree->dxy_Kp_pi > 0.4 ) continue; */
                /* if( ftree->dz_Kp_pi > 0.4 ) continue; */

                /* if( ftree->dxy_Km_pi > 0.4 ) continue; */
                /* if( ftree->dz_Km_pi > 0.4 ) continue; */

                /* if( ftree->dxy_pi_phi > 6 ) continue; */
                /* if( ftree->dz_pi_phi > 6 ) continue; */

                std::vector<reco::TransientTrack> Ds_Tracks = {
                    phiFit_Tracks[0],
                    phiFit_Tracks[1],
                    (*theB).build(pi_PF.pseudoTrack())
                };

                TransientVertex Ds_Vertex = KVFitter.vertex(Ds_Tracks);

                if( !(Ds_Vertex.isValid()) ) continue;
                if( !(Ds_Vertex.hasRefittedTracks()) ) continue;

                ftree->DsFit_CHI2 = Ds_Vertex.totalChiSquared();
                ftree->DsFit_NDOF = Ds_Vertex.degreesOfFreedom();
                ftree->DsFit_CHI2NDOF = Ds_Vertex.normalisedChiSquared();
                ftree->DsFit_ENDVX_X = Ds_Vertex.position().x();
                ftree->DsFit_ENDVX_Y = Ds_Vertex.position().y();
                ftree->DsFit_ENDVX_Z = Ds_Vertex.position().z();
                ftree->DsFit_ENDVX_XERR = std::sqrt(Ds_Vertex.positionError().cxx());
                ftree->DsFit_ENDVX_YERR = std::sqrt(Ds_Vertex.positionError().cyy());
                ftree->DsFit_ENDVX_ZERR = std::sqrt(Ds_Vertex.positionError().czz());

                std::vector<reco::TransientTrack> DsFit_Tracks = Ds_Vertex.refittedTracks();

                // Ds fit Kp
                TLorentzVector DsFit_Kp_P4;
                DsFit_Kp_P4.SetXYZM(DsFit_Tracks[0].track().px(), DsFit_Tracks[0].track().py(), DsFit_Tracks[0].track().pz(), Mass_K);
                ftree->DsFit_Kp_ETA = DsFit_Kp_P4.Eta();
                ftree->DsFit_Kp_PHI = DsFit_Kp_P4.Phi();
                ftree->DsFit_Kp_P = DsFit_Kp_P4.P();
                ftree->DsFit_Kp_PT = DsFit_Kp_P4.Pt();
                ftree->DsFit_Kp_PX = DsFit_Kp_P4.Px();
                ftree->DsFit_Kp_PY = DsFit_Kp_P4.Py();
                ftree->DsFit_Kp_PZ = DsFit_Kp_P4.Pz();

                // Ds fit Km
                TLorentzVector DsFit_Km_P4;
                DsFit_Km_P4.SetXYZM(DsFit_Tracks[1].track().px(), DsFit_Tracks[1].track().py(), DsFit_Tracks[1].track().pz(), Mass_K);
                ftree->DsFit_Km_ETA = DsFit_Km_P4.Eta();
                ftree->DsFit_Km_PHI = DsFit_Km_P4.Phi();
                ftree->DsFit_Km_P = DsFit_Km_P4.P();
                ftree->DsFit_Km_PT = DsFit_Km_P4.Pt();
                ftree->DsFit_Km_PX = DsFit_Km_P4.Px();
                ftree->DsFit_Km_PY = DsFit_Km_P4.Py();
                ftree->DsFit_Km_PZ = DsFit_Km_P4.Pz();

                // Ds fit pi
                TLorentzVector DsFit_pi_P4;
                DsFit_pi_P4.SetXYZM(DsFit_Tracks[2].track().px(), DsFit_Tracks[2].track().py(), DsFit_Tracks[2].track().pz(), Mass_pi);
                ftree->DsFit_pi_ETA = DsFit_pi_P4.Eta();
                ftree->DsFit_pi_PHI = DsFit_pi_P4.Phi();
                ftree->DsFit_pi_P = DsFit_pi_P4.P();
                ftree->DsFit_pi_PT = DsFit_pi_P4.Pt();
                ftree->DsFit_pi_PX = DsFit_pi_P4.Px();
                ftree->DsFit_pi_PY = DsFit_pi_P4.Py();
                ftree->DsFit_pi_PZ = DsFit_pi_P4.Pz();

                // Ds fit phi
                TLorentzVector DsFit_phi_P4 = DsFit_Kp_P4 + DsFit_Km_P4;
                ftree->DsFit_phi_ETA = DsFit_phi_P4.Eta();
                ftree->DsFit_phi_PHI = DsFit_phi_P4.Phi();
                ftree->DsFit_phi_P = DsFit_phi_P4.P();
                ftree->DsFit_phi_PT = DsFit_phi_P4.Pt();
                ftree->DsFit_phi_PX = DsFit_phi_P4.Px();
                ftree->DsFit_phi_PY = DsFit_phi_P4.Py();
                ftree->DsFit_phi_PZ = DsFit_phi_P4.Pz();
                ftree->DsFit_phi_M = DsFit_phi_P4.M();

                // Ds fit Ds
                TLorentzVector DsFit_Ds_P4 = DsFit_Kp_P4 + DsFit_Km_P4 + DsFit_pi_P4;
                ftree->DsFit_Ds_ETA = DsFit_Ds_P4.Eta();
                ftree->DsFit_Ds_PHI = DsFit_Ds_P4.Phi();
                ftree->DsFit_Ds_P = DsFit_Ds_P4.P();
                ftree->DsFit_Ds_PT = DsFit_Ds_P4.Pt();
                ftree->DsFit_Ds_PX = DsFit_Ds_P4.Px();
                ftree->DsFit_Ds_PY = DsFit_Ds_P4.Py();
                ftree->DsFit_Ds_PZ = DsFit_Ds_P4.Pz();
                ftree->DsFit_Ds_M = DsFit_Ds_P4.M();

                TLorentzVector DsFit_Mconstraint_phi_P4;
                DsFit_Mconstraint_phi_P4.SetXYZM(DsFit_phi_P4.Px(), DsFit_phi_P4.Py(), DsFit_phi_P4.Pz(), Mass_phi);
                TLorentzVector DsFit_Mconstraint_Ds_P4 = DsFit_Mconstraint_phi_P4 + DsFit_pi_P4;
                ftree->DsFit_Mconstraint_Ds_M = DsFit_Mconstraint_Ds_P4.M();

                ftree->DsFit_Kp_PP = DsFit_Kp_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Kp_PL = DsFit_Kp_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();
                ftree->DsFit_Km_PP = DsFit_Km_P4.Vect().Pt(DsFit_phi_P4.Vect());
                ftree->DsFit_Km_PL = DsFit_Km_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P();

                ftree->DsFit_phi_PP = DsFit_phi_P4.Vect().Pt(DsFit_Ds_P4.Vect());
                ftree->DsFit_phi_PL = DsFit_phi_P4.Vect().Dot(DsFit_Ds_P4.Vect())/DsFit_Ds_P4.P();
                ftree->DsFit_pi_PP = DsFit_pi_P4.Vect().Pt(DsFit_Ds_P4.Vect());
                ftree->DsFit_pi_PL = DsFit_pi_P4.Vect().Dot(DsFit_Ds_P4.Vect())/DsFit_Ds_P4.P();

                ftree->DsFit_dR_Kp_Km = reco::deltaR(DsFit_Kp_P4.Eta(), DsFit_Kp_P4.Phi(), DsFit_Km_P4.Eta(), DsFit_Km_P4.Phi());
                ftree->DsFit_dR_Kp_phi = reco::deltaR(DsFit_Kp_P4.Eta(), DsFit_Kp_P4.Phi(), DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi());
                ftree->DsFit_dR_Km_phi = reco::deltaR(DsFit_Km_P4.Eta(), DsFit_Km_P4.Phi(), DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi());
                ftree->DsFit_dR_Kp_pi = reco::deltaR(DsFit_Kp_P4.Eta(), DsFit_Kp_P4.Phi(), DsFit_pi_P4.Eta(), DsFit_pi_P4.Phi());
                ftree->DsFit_dR_Km_pi = reco::deltaR(DsFit_Km_P4.Eta(), DsFit_Km_P4.Phi(), DsFit_pi_P4.Eta(), DsFit_pi_P4.Phi());
                ftree->DsFit_dR_pi_phi = reco::deltaR(DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi(), DsFit_pi_P4.Eta(), DsFit_pi_P4.Phi());
                ftree->DsFit_dR_Kp_Ds = reco::deltaR(DsFit_Kp_P4.Eta(), DsFit_Kp_P4.Phi(), DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_Km_Ds = reco::deltaR(DsFit_Km_P4.Eta(), DsFit_Km_P4.Phi(), DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_phi_Ds = reco::deltaR(DsFit_phi_P4.Eta(), DsFit_phi_P4.Phi(), DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi());
                ftree->DsFit_dR_pi_Ds = reco::deltaR(DsFit_pi_P4.Eta(), DsFit_pi_P4.Phi(), DsFit_Ds_P4.Eta(), DsFit_Ds_P4.Phi());

                ftree->dxy_Kp_Ds = sqrt(pow(ftree->Kp_ORIVX_X-ftree->DsFit_ENDVX_X,2) + pow(ftree->Kp_ORIVX_Y-ftree->DsFit_ENDVX_Y,2));
                ftree->dxy_Km_Ds = sqrt(pow(ftree->Km_ORIVX_X-ftree->DsFit_ENDVX_X,2) + pow(ftree->Km_ORIVX_Y-ftree->DsFit_ENDVX_Y,2));
                ftree->dxy_pi_Ds = sqrt(pow(ftree->pi_ORIVX_X-ftree->DsFit_ENDVX_X,2) + pow(ftree->pi_ORIVX_Y-ftree->DsFit_ENDVX_Y,2));
                ftree->dxy_phi_Ds = sqrt(pow(ftree->phiFit_ENDVX_X-ftree->DsFit_ENDVX_X,2) + pow(ftree->phiFit_ENDVX_Y-ftree->DsFit_ENDVX_Y,2));

                ftree->dz_Kp_Ds = abs(ftree->Kp_ORIVX_Z-ftree->DsFit_ENDVX_Z);
                ftree->dz_Km_Ds = abs(ftree->Km_ORIVX_Z-ftree->DsFit_ENDVX_Z);
                ftree->dz_pi_Ds = abs(ftree->pi_ORIVX_Z-ftree->DsFit_ENDVX_Z);
                ftree->dz_phi_Ds = abs(ftree->phiFit_ENDVX_Z-ftree->DsFit_ENDVX_Z);

                if( ftree->DsFit_dR_Kp_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_Km_Ds > 0.15 ) continue;
                if( ftree->DsFit_dR_pi_Ds > 0.4 ) continue;
                if( ftree->DsFit_dR_phi_Ds > 0.15 ) continue;
                /* if( ftree->dxy_phi_Ds > 4 ) continue; */
                /* if( ftree->dz_phi_Ds > 4 ) continue; */
                if( ftree->DsFit_CHI2NDOF < 0 ) continue;
                if( ftree->DsFit_CHI2NDOF > 10 ) continue;
                if( ftree->DsFit_Ds_M < 1.85) continue;
                if( ftree->DsFit_Ds_M > 2.1) continue;

                double alpha_phi = (ftree->phiFit_Kp_PL - ftree->phiFit_Km_PL) / (ftree->phiFit_Kp_PL + ftree->phiFit_Km_PL);
                double alpha_Ds = (ftree->DsFit_phi_PL - ftree->DsFit_pi_PL) / (ftree->DsFit_phi_PL + ftree->DsFit_pi_PL);

                double beta_phi = std::sqrt(pow(ftree->phiFit_phi_P,2) / (pow(Mass_phi,2) + pow(ftree->phiFit_phi_P,2)));
                double beta_Ds = std::sqrt(pow(ftree->DsFit_Ds_P,2) / (pow(Mass_Ds,2) + pow(ftree->DsFit_Ds_P,2)));

                double var_phi = pow(ftree->phiFit_Kp_PP,2) + pow(alpha_phi*beta_phi*Mass_phi,2)/4;
                double var_Ds = pow(ftree->DsFit_phi_PP,2) + pow(beta_Ds*(Est_phi-Est_pi)-alpha_Ds*beta_Ds*(Est_phi+Est_pi),2)/4;

                if(var_phi < 0.006) continue;
                if(var_phi > 0.028) continue;
                if(var_Ds < 0.4) continue;
                if(var_Ds > 0.6) continue;

                ftree->num_reco_Ds++;

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

                std::vector<reco::TransientTrack> PV_ttracks;
                
                for(size_t l=0; l<packedPF->size(); l++){

                    if(int(l) ==  idx_Kp_vec[i]) continue;
                    if(int(l) ==  idx_Km_vec[j]) continue;
                    if(int(l) ==  idx_pi_vec[k]) continue;

                    const auto& otherpf = (*packedPF)[l];

                    if( !(otherpf.trackHighPurity()) ) continue;
                    if( !(otherpf.hasTrackDetails()) ) continue;

                    if(otherpf.pt() > 1){
                        PV_ttracks.push_back( (*theB).build(otherpf.pseudoTrack()) );
                    }
                }

                if(PV_ttracks.size() > 1){
                    TransientVertex PV_tvertex = AVFitter.vertex(PV_ttracks);
                    if( PV_tvertex.isValid() ){ 
                        ftree->PV_CHI2 = PV_tvertex.totalChiSquared();
                        ftree->PV_NDOF = PV_tvertex.degreesOfFreedom();
                        ftree->PV_CHI2NDOF = PV_tvertex.normalisedChiSquared();
                        ftree->PV_X = PV_tvertex.position().x();
                        ftree->PV_Y = PV_tvertex.position().y();
                        ftree->PV_Z = PV_tvertex.position().z();
                        ftree->PV_XERR = std::sqrt(PV_tvertex.positionError().cxx());
                        ftree->PV_YERR = std::sqrt(PV_tvertex.positionError().cyy());
                        ftree->PV_ZERR = std::sqrt(PV_tvertex.positionError().czz());
                    }
                }

                ftree->Fill_PV();
            }
        }
    }

    int num_candidates = ftree->match_entry_vec.size();

    double maxPT = 0;
    int maxidx = -1;

    for(int i=0; i<num_candidates; i++){
        if(ftree->DsFit_Ds_PT_vec[i]>maxPT){
            maxPT = ftree->DsFit_Ds_PT_vec[i];
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

void PVStudy::beginJob() {}

void PVStudy::endJob() {}

void PVStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("PVStudy", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PVStudy);

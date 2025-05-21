// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      SelectionStudy
//
/**\class SelectionStudy SelectionStudy.cc EDAnalyzer/GenParticleAnalyzer/plugins/SelectionStudy.cc

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
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

#include "EDAnalyzers/GenParticleAnalyzer/interface/SSTree.h"


class SelectionStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit SelectionStudy(const edm::ParameterSet&);
        ~SelectionStudy();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        const edm::Service<TFileService> fs;
        SSTree *ftree;

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

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

SelectionStudy::SelectionStudy(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    ftree = new SSTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

SelectionStudy::~SelectionStudy() {}

void SelectionStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    ftree->Init();

    edm::Handle<reco::GenParticleCollection> prunedGen;
    iEvent.getByToken(prunedGenToken, prunedGen);

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

    KalmanVertexFitter fitter(true);

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
        }
    }

    ftree->Gen_Reset();
    ftree->Match_Reset();
    
    if(ftree->num_Gen_Kp==1 && ftree->num_Gen_Km==1 && ftree->num_Gen_pi==1 && ftree->num_Gen_phi==1 && ftree->num_Gen_Ds==1){

        const auto& Kp_GP = (*prunedGen)[Gen_Kp_idx[0]];
        const auto& Km_GP = (*prunedGen)[Gen_Km_idx[0]];
        const auto& pi_GP = (*prunedGen)[Gen_pi_idx[0]];
        const auto& phi_GP = (*prunedGen)[Gen_phi_idx[0]];
        const auto& Ds_GP = (*prunedGen)[Gen_Ds_idx[0]];

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
        ftree->Gen_dR_phi_pi = reco::deltaR(phi_GP, pi_GP);
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
        if(bestKp.dr < 0.03) ftree->num_match_Kp++;
        if(bestKm.dr < 0.03) ftree->num_match_Km++;
        if(bestpi.dr < 0.03) ftree->num_match_pi++;

        ftree->match_Kp_idx = bestKp.index;
        ftree->match_Km_idx = bestKm.index;
        ftree->match_pi_idx = bestpi.index;

        const auto& match_Kp_PF = (*packedPF)[bestKp.index];
        const auto& match_Km_PF = (*packedPF)[bestKm.index];
        const auto& match_pi_PF = (*packedPF)[bestpi.index];

        if( bestKp.dr < 0.03 && bestKm.dr < 0.03 && bestpi.dr < 0.03 && match_Kp_PF.hasTrackDetails() && match_Km_PF.hasTrackDetails() && match_pi_PF.hasTrackDetails() ){

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
            ftree->match_dR_phi_pi = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_pi_P4.Eta(), match_pi_P4.Phi());
            ftree->match_dR_Kp_Ds = reco::deltaR(match_Kp_P4.Eta(), match_Kp_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
            ftree->match_dR_Km_Ds = reco::deltaR(match_Km_P4.Eta(), match_Km_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
            ftree->match_dR_phi_Ds = reco::deltaR(match_phi_P4.Eta(), match_phi_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());
            ftree->match_dR_pi_Ds = reco::deltaR(match_pi_P4.Eta(), match_pi_P4.Phi(), match_Ds_P4.Eta(), match_Ds_P4.Phi());

            // phi fit
            std::vector<reco::TransientTrack> match_phi_Tracks = {
                (*theB).build(match_Kp_PF.pseudoTrack()),
                (*theB).build(match_Km_PF.pseudoTrack())
            };

            TransientVertex match_phi_Vertex = fitter.vertex(match_phi_Tracks);

            if( match_phi_Vertex.isValid() && match_phi_Vertex.hasRefittedTracks() ){
               
                // phi fit variables 
                ftree->match_phiFit_CHI2 = match_phi_Vertex.totalChiSquared();
                ftree->match_phiFit_NDOF = match_phi_Vertex.degreesOfFreedom();
                ftree->match_phiFit_CHI2NDOF = match_phi_Vertex.normalisedChiSquared();
                ftree->match_phiFit_ENDVX_X = match_phi_Vertex.position().x();
                ftree->match_phiFit_ENDVX_Y = match_phi_Vertex.position().y();
                ftree->match_phiFit_ENDVX_Z = match_phi_Vertex.position().z();

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
                ftree->match_phiFit_dR_phi_pi = reco::deltaR(match_phiFit_phi_P4.Eta(), match_phiFit_phi_P4.Phi(), match_phiFit_pi_P4.Eta(), match_phiFit_pi_P4.Phi());
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

                TransientVertex match_Ds_Vertex = fitter.vertex(match_Ds_Tracks);

                if( match_Ds_Vertex.isValid() && match_Ds_Vertex.hasRefittedTracks() ){
               
                    // Ds fit variables 
                    ftree->match_DsFit_CHI2 = match_Ds_Vertex.totalChiSquared();
                    ftree->match_DsFit_NDOF = match_Ds_Vertex.degreesOfFreedom();
                    ftree->match_DsFit_CHI2NDOF = match_Ds_Vertex.normalisedChiSquared();
                    ftree->match_DsFit_ENDVX_X = match_Ds_Vertex.position().x();
                    ftree->match_DsFit_ENDVX_Y = match_Ds_Vertex.position().y();
                    ftree->match_DsFit_ENDVX_Z = match_Ds_Vertex.position().z();

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
                    ftree->match_DsFit_dR_phi_pi = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi());
                    ftree->match_DsFit_dR_Kp_Ds = reco::deltaR(match_DsFit_Kp_P4.Eta(), match_DsFit_Kp_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                    ftree->match_DsFit_dR_Km_Ds = reco::deltaR(match_DsFit_Km_P4.Eta(), match_DsFit_Km_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                    ftree->match_DsFit_dR_phi_Ds = reco::deltaR(match_DsFit_phi_P4.Eta(), match_DsFit_phi_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                    ftree->match_DsFit_dR_pi_Ds = reco::deltaR(match_DsFit_pi_P4.Eta(), match_DsFit_pi_P4.Phi(), match_DsFit_Ds_P4.Eta(), match_DsFit_Ds_P4.Phi());
                
                    ftree->Match_Fill_Vector();
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

            TLorentzVector

            ftree->dR_Kp_Km = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Km_P4.Eta(), Km_P4.Phi());
            ftree->dR_Kp_phi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());
            ftree->dR_Km_phi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), phi_P4.Eta(), phi_P4.Phi());
            ftree->dR_Kp_pi = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), pi_P4.Eta(), pi_P4.Phi());
            ftree->dR_Km_pi = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), pi_P4.Eta(), pi_P4.Phi());
            ftree->dR_phi_pi = reco::deltaR(phi_P4.Eta(), phi_P4.Phi(), pi_P4.Eta(), pi_P4.Phi());
            ftree->dR_Kp_Ds = reco::deltaR(Kp_P4.Eta(), Kp_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
            ftree->dR_Km_Ds = reco::deltaR(Km_P4.Eta(), Km_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
            ftree->dR_phi_Ds = reco::deltaR(phi_P4.Eta(), phi_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());
            ftree->dR_pi_Ds = reco::deltaR(pi_P4.Eta(), pi_P4.Phi(), Ds_P4.Eta(), Ds_P4.Phi());


    ftree->tree->Fill();
    return;
}

bool SelectionStudy::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
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

void SelectionStudy::beginJob() {}

void SelectionStudy::endJob() {}

void SelectionStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("SelectionStudy", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectionStudy);

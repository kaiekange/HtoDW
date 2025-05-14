// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenTrackMatch
//
/**\class GenTrackMatch GenTrackMatch.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenTrackMatch.cc

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
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GenTrackMatch : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenTrackMatch(const edm::ParameterSet&);
        ~GenTrackMatch();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;
        void vector_clear();

        edm::Service<TFileService> fs;
        TTree *tree_{nullptr};

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

        // Gen particles info
        int num_Gen_Kp;
        int num_Gen_Km;
        int num_Gen_pi;
        int num_Gen_phi;
        int num_Gen_Ds;

        std::vector<double> Gen_Kp_ETA;
        std::vector<double> Gen_Kp_PHI;
        std::vector<double> Gen_Kp_ORIVX_X;
        std::vector<double> Gen_Kp_ORIVX_Y;
        std::vector<double> Gen_Kp_ORIVX_Z;
        std::vector<double> Gen_Kp_P;
        std::vector<double> Gen_Kp_PT;
        std::vector<double> Gen_Kp_PX;
        std::vector<double> Gen_Kp_PY;
        std::vector<double> Gen_Kp_PZ;
        std::vector<double> Gen_Kp_PP;
        std::vector<double> Gen_Kp_PL;

        std::vector<double> Gen_Km_ETA;
        std::vector<double> Gen_Km_PHI;
        std::vector<double> Gen_Km_ORIVX_X;
        std::vector<double> Gen_Km_ORIVX_Y;
        std::vector<double> Gen_Km_ORIVX_Z;
        std::vector<double> Gen_Km_P;
        std::vector<double> Gen_Km_PT;
        std::vector<double> Gen_Km_PX;
        std::vector<double> Gen_Km_PY;
        std::vector<double> Gen_Km_PZ;
        std::vector<double> Gen_Km_PP;
        std::vector<double> Gen_Km_PL;

        std::vector<double> Gen_pi_ETA;
        std::vector<double> Gen_pi_PHI;
        std::vector<double> Gen_pi_ORIVX_X;
        std::vector<double> Gen_pi_ORIVX_Y;
        std::vector<double> Gen_pi_ORIVX_Z;
        std::vector<double> Gen_pi_P;
        std::vector<double> Gen_pi_PT;
        std::vector<double> Gen_pi_PX;
        std::vector<double> Gen_pi_PY;
        std::vector<double> Gen_pi_PZ;
        std::vector<double> Gen_pi_PP;
        std::vector<double> Gen_pi_PL;

        std::vector<double> Gen_phi_ETA;
        std::vector<double> Gen_phi_PHI;
        std::vector<double> Gen_phi_ORIVX_X;
        std::vector<double> Gen_phi_ORIVX_Y;
        std::vector<double> Gen_phi_ORIVX_Z;
        std::vector<double> Gen_phi_P;
        std::vector<double> Gen_phi_PT;
        std::vector<double> Gen_phi_PX;
        std::vector<double> Gen_phi_PY;
        std::vector<double> Gen_phi_PZ;
        std::vector<double> Gen_phi_PP;
        std::vector<double> Gen_phi_PL;

        std::vector<double> Gen_Ds_ETA;
        std::vector<double> Gen_Ds_PHI;
        std::vector<double> Gen_Ds_ORIVX_X;
        std::vector<double> Gen_Ds_ORIVX_Y;
        std::vector<double> Gen_Ds_ORIVX_Z;
        std::vector<double> Gen_Ds_P;
        std::vector<double> Gen_Ds_PT;
        std::vector<double> Gen_Ds_PX;
        std::vector<double> Gen_Ds_PY;
        std::vector<double> Gen_Ds_PZ;

        std::vector<double> Gen_dR_Kp_Km;
        std::vector<double> Gen_dR_Kp_pi;
        std::vector<double> Gen_dR_Kp_phi;
        std::vector<double> Gen_dR_Kp_Ds;
        std::vector<double> Gen_dR_Km_pi;
        std::vector<double> Gen_dR_Km_phi;
        std::vector<double> Gen_dR_Km_Ds;
        std::vector<double> Gen_dR_pi_phi;
        std::vector<double> Gen_dR_pi_Ds;
        std::vector<double> Gen_dR_phi_Ds;

        std::vector<double> Gen_d_Kp_Km;
        std::vector<double> Gen_d_Kp_pi;
        std::vector<double> Gen_d_Km_pi;

        // matched particles info
        int num_match_Kp;
        int num_match_Km;
        int num_match_pi;
        std::vector<double> Kp_ETA;
        std::vector<double> Kp_PHI;
        std::vector<double> Kp_ORIVX_X;
        std::vector<double> Kp_ORIVX_Y;
        std::vector<double> Kp_ORIVX_Z;
        std::vector<double> Kp_P;
        std::vector<double> Kp_PT;
        std::vector<double> Kp_PX;
        std::vector<double> Kp_PY;
        std::vector<double> Kp_PZ;
        std::vector<double> Kp_PP;
        std::vector<double> Kp_PL;

        std::vector<double> Km_ETA;
        std::vector<double> Km_PHI;
        std::vector<double> Km_ORIVX_X;
        std::vector<double> Km_ORIVX_Y;
        std::vector<double> Km_ORIVX_Z;
        std::vector<double> Km_P;
        std::vector<double> Km_PT;
        std::vector<double> Km_PX;
        std::vector<double> Km_PY;
        std::vector<double> Km_PZ;
        std::vector<double> Km_PP;
        std::vector<double> Km_PL;

        std::vector<double> pi_ETA;
        std::vector<double> pi_PHI;
        std::vector<double> pi_ORIVX_X;
        std::vector<double> pi_ORIVX_Y;
        std::vector<double> pi_ORIVX_Z;
        std::vector<double> pi_P;
        std::vector<double> pi_PT;
        std::vector<double> pi_PX;
        std::vector<double> pi_PY;
        std::vector<double> pi_PZ;
        std::vector<double> pi_PP;
        std::vector<double> pi_PL;

        std::vector<double> phi_ETA;
        std::vector<double> phi_PHI;
        std::vector<double> phi_P;
        std::vector<double> phi_PT;
        std::vector<double> phi_PX;
        std::vector<double> phi_PY;
        std::vector<double> phi_PZ;
        std::vector<double> phi_PP;
        std::vector<double> phi_PL;
        std::vector<double> phi_CHI2;
        std::vector<double> phi_M;
        std::vector<double> phi_ENDVX_X;
        std::vector<double> phi_ENDVX_Y;
        std::vector<double> phi_ENDVX_Z;

        std::vector<double> Ds_ETA;
        std::vector<double> Ds_PHI;
        std::vector<double> Ds_P;
        std::vector<double> Ds_PT;
        std::vector<double> Ds_PX;
        std::vector<double> Ds_PY;
        std::vector<double> Ds_PZ;
        std::vector<double> Ds_CHI2;
        std::vector<double> Ds_M;
        std::vector<double> Ds_ENDVX_X;
        std::vector<double> Ds_ENDVX_Y;
        std::vector<double> Ds_ENDVX_Z;

        std::vector<double> dR_Kp_Km;
        std::vector<double> dR_Kp_pi;
        std::vector<double> dR_Kp_phi;
        std::vector<double> dR_Kp_Ds;
        std::vector<double> dR_Km_pi;
        std::vector<double> dR_Km_phi;
        std::vector<double> dR_Km_Ds;
        std::vector<double> dR_pi_phi;
        std::vector<double> dR_pi_Ds;
        std::vector<double> dR_phi_Ds;

        std::vector<double> d_Kp_Km;
        std::vector<double> d_Kp_pi;
        std::vector<double> d_Kp_phi;
        std::vector<double> d_Kp_Ds;
        std::vector<double> d_Km_pi;
        std::vector<double> d_Km_phi;
        std::vector<double> d_Km_Ds;
        std::vector<double> d_pi_phi;
        std::vector<double> d_pi_Ds;
        std::vector<double> d_phi_Ds;

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

GenTrackMatch::GenTrackMatch(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
}

GenTrackMatch::~GenTrackMatch() {}

void GenTrackMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    vector_clear();

    edm::Handle<reco::GenParticleCollection> prunedGen;
    iEvent.getByToken(prunedGenToken, prunedGen);

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

    TVector3 Gen_p3_Kp, Gen_p3_Km, Gen_p3_pi;

    for(size_t i=0; i<prunedGen->size(); i++){
        const auto& gp = (*prunedGen)[i];
        const int pdg = gp.pdgId();
        if( pdg==321 && hasAncestor(gp, 333, -321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_Kp_idx.push_back(i);
            num_Gen_Kp++;
        } else if( pdg==-321 && hasAncestor(gp, 333, 321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_Km_idx.push_back(i);
            num_Gen_Km++;
        } else if( pdg==211 && hasAncestor(gp, 431, 333) && hasAncestor(gp, 25, -24) ){
            Gen_pi_idx.push_back(i);
            num_Gen_pi++;
        } else if( pdg==333 && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) && gp.numberOfDaughters()==2 ){
            int pdg0 = gp.daughter(0)->pdgId();
            int pdg1 = gp.daughter(1)->pdgId();
            if( (pdg0 == 321 && pdg1 == -321) || (pdg0 == -321 && pdg1 == 321)){
                Gen_phi_idx.push_back(i);
                num_Gen_phi++;
            }
        } else if( pdg==431 && hasAncestor(gp, 25, -24) && gp.isHardProcess() ){
            Gen_Ds_idx.push_back(i);
            num_Gen_Ds++;
        }
    }

    if(num_Gen_Kp==1 && num_Gen_Km==1 && num_Gen_pi==1 && num_Gen_phi==1 && num_Gen_Ds==1){

        const auto& Kp_GP = (*prunedGen)[Gen_Kp_idx[0]];
        const auto& Km_GP = (*prunedGen)[Gen_Km_idx[0]];
        const auto& pi_GP = (*prunedGen)[Gen_pi_idx[0]];
        const auto& phi_GP = (*prunedGen)[Gen_phi_idx[0]];
        const auto& Ds_GP = (*prunedGen)[Gen_Ds_idx[0]];

        Gen_Kp_ETA.push_back(Kp_GP.eta());
        Gen_Kp_PHI.push_back(Kp_GP.phi());
        Gen_Kp_ORIVX_X.push_back(Kp_GP.vx());
        Gen_Kp_ORIVX_Y.push_back(Kp_GP.vy());
        Gen_Kp_ORIVX_Z.push_back(Kp_GP.vz());
        Gen_Kp_P.push_back(Kp_GP.p());
        Gen_Kp_PT.push_back(Kp_GP.pt());
        Gen_Kp_PX.push_back(Kp_GP.px());
        Gen_Kp_PY.push_back(Kp_GP.py());
        Gen_Kp_PZ.push_back(Kp_GP.pz());
        TVector3 Gen_Kp_P3(Gen_Kp_PX[0], Gen_Kp_PY[0], Gen_Kp_PZ[0]);

        Gen_Km_ETA.push_back(Km_GP.eta());
        Gen_Km_PHI.push_back(Km_GP.phi());
        Gen_Km_ORIVX_X.push_back(Km_GP.vx());
        Gen_Km_ORIVX_Y.push_back(Km_GP.vy());
        Gen_Km_ORIVX_Z.push_back(Km_GP.vz());
        Gen_Km_P.push_back(Km_GP.p());
        Gen_Km_PT.push_back(Km_GP.pt());
        Gen_Km_PX.push_back(Km_GP.px());
        Gen_Km_PY.push_back(Km_GP.py());
        Gen_Km_PZ.push_back(Km_GP.pz());
        TVector3 Gen_Km_P3(Gen_Km_PX[0], Gen_Km_PY[0], Gen_Km_PZ[0]);

        Gen_pi_ETA.push_back(pi_GP.eta());
        Gen_pi_PHI.push_back(pi_GP.phi());
        Gen_pi_ORIVX_X.push_back(pi_GP.vx());
        Gen_pi_ORIVX_Y.push_back(pi_GP.vy());
        Gen_pi_ORIVX_Z.push_back(pi_GP.vz());
        Gen_pi_P.push_back(pi_GP.p());
        Gen_pi_PT.push_back(pi_GP.pt());
        Gen_pi_PX.push_back(pi_GP.px());
        Gen_pi_PY.push_back(pi_GP.py());
        Gen_pi_PZ.push_back(pi_GP.pz());
        TVector3 Gen_pi_P3(Gen_pi_PX[0], Gen_pi_PY[0], Gen_pi_PZ[0]);

        Gen_phi_ETA.push_back(phi_GP.eta());
        Gen_phi_PHI.push_back(phi_GP.phi());
        Gen_phi_ORIVX_X.push_back(phi_GP.vx());
        Gen_phi_ORIVX_Y.push_back(phi_GP.vy());
        Gen_phi_ORIVX_Z.push_back(phi_GP.vz());
        Gen_phi_P.push_back(phi_GP.p());
        Gen_phi_PT.push_back(phi_GP.pt());
        Gen_phi_PX.push_back(phi_GP.px());
        Gen_phi_PY.push_back(phi_GP.py());
        Gen_phi_PZ.push_back(phi_GP.pz());
        TVector3 Gen_phi_P3(Gen_phi_PX[0], Gen_phi_PY[0], Gen_phi_PZ[0]);

        Gen_Ds_ETA.push_back(Ds_GP.eta());
        Gen_Ds_PHI.push_back(Ds_GP.phi());
        Gen_Ds_ORIVX_X.push_back(Ds_GP.vx());
        Gen_Ds_ORIVX_Y.push_back(Ds_GP.vy());
        Gen_Ds_ORIVX_Z.push_back(Ds_GP.vz());
        Gen_Ds_P.push_back(Ds_GP.p());
        Gen_Ds_PT.push_back(Ds_GP.pt());
        Gen_Ds_PX.push_back(Ds_GP.px());
        Gen_Ds_PY.push_back(Ds_GP.py());
        Gen_Ds_PZ.push_back(Ds_GP.pz());
        TVector3 Gen_Ds_P3(Gen_Ds_PX[0], Gen_Ds_PY[0], Gen_Ds_PZ[0]);

        Gen_Kp_PP.push_back(Gen_Kp_P3.Pt(Gen_phi_P3));
        Gen_Kp_PL.push_back(Gen_Kp_P3.Dot(Gen_phi_P3)/Gen_phi_P[0]);
        Gen_Km_PP.push_back(Gen_Km_P3.Pt(Gen_phi_P3));
        Gen_Km_PL.push_back(Gen_Km_P3.Dot(Gen_phi_P3)/Gen_phi_P[0]);
        Gen_pi_PP.push_back(Gen_pi_P3.Pt(Gen_Ds_P3));
        Gen_pi_PL.push_back(Gen_pi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P[0]);
        Gen_phi_PP.push_back(Gen_phi_P3.Pt(Gen_Ds_P3));
        Gen_phi_PL.push_back(Gen_phi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P[0]);
        
        Gen_dR_Kp_Km.push_back(reco::deltaR(Kp_GP, Km_GP));
        Gen_dR_Kp_pi.push_back(reco::deltaR(Kp_GP, pi_GP));
        Gen_dR_Kp_phi.push_back(reco::deltaR(Kp_GP, phi_GP));
        Gen_dR_Kp_Ds.push_back(reco::deltaR(Kp_GP, Ds_GP));
        Gen_dR_Km_pi.push_back(reco::deltaR(Km_GP, pi_GP));
        Gen_dR_Km_phi.push_back(reco::deltaR(Km_GP, phi_GP));
        Gen_dR_Km_Ds.push_back(reco::deltaR(Km_GP, Ds_GP));
        Gen_dR_pi_phi.push_back(reco::deltaR(pi_GP, phi_GP));
        Gen_dR_pi_Ds.push_back(reco::deltaR(pi_GP, Ds_GP));
        Gen_dR_phi_Ds.push_back(reco::deltaR(phi_GP, Ds_GP));

        Gen_d_Kp_Km.push_back((Kp_GP.vertex() - Km_GP.vertex()).R());
        Gen_d_Kp_pi.push_back((Kp_GP.vertex() - pi_GP.vertex()).R());
        Gen_d_Km_pi.push_back((Km_GP.vertex() - pi_GP.vertex()).R());

        struct MatchInfo {int index; float dr;};
        std::vector<MatchInfo> MatchKp, MatchKm, Matchpi;

        for(size_t i=0; i<packedPF->size(); i++){
            const auto& pf = (*packedPF)[i];

            const float dr_Kp = reco::deltaR(pf.eta(), pf.phi(), Gen_Kp_ETA[0], Gen_Kp_PHI[0]);
            const float dr_Km = reco::deltaR(pf.eta(), pf.phi(), Gen_Km_ETA[0], Gen_Km_PHI[0]);
            const float dr_pi = reco::deltaR(pf.eta(), pf.phi(), Gen_pi_ETA[0], Gen_pi_PHI[0]);

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
        if(bestKp.dr < 0.03) num_match_Kp++;
        if(bestKm.dr < 0.03) num_match_Km++;
        if(bestpi.dr < 0.03) num_match_pi++;

        if( bestKp.dr < 0.03 && bestKm.dr < 0.03 && bestpi.dr < 0.03 && (*packedPF)[bestKp.index].hasTrackDetails() && (*packedPF)[bestKm.index].hasTrackDetails() && (*packedPF)[bestpi.index].hasTrackDetails() ){

            const auto& Kp_PF = (*packedPF)[bestKp.index];
            const auto& Km_PF = (*packedPF)[bestKm.index];
            const auto& pi_PF = (*packedPF)[bestpi.index];

            Kp_ETA.push_back(Kp_PF.eta());
            Kp_PHI.push_back(Kp_PF.phi());
            Kp_ORIVX_X.push_back(Kp_PF.vx());
            Kp_ORIVX_Y.push_back(Kp_PF.vy());
            Kp_ORIVX_Z.push_back(Kp_PF.vz());
            Kp_P.push_back(Kp_PF.p());
            Kp_PT.push_back(Kp_PF.pt());
            Kp_PX.push_back(Kp_PF.px());
            Kp_PY.push_back(Kp_PF.py());
            Kp_PZ.push_back(Kp_PF.pz());
            math::XYZPoint Kp_originvertex(Kp_ORIVX_X[0], Kp_ORIVX_Y[0], Kp_ORIVX_Z[0]);

            Km_ETA.push_back(Km_PF.eta());
            Km_PHI.push_back(Km_PF.phi());
            Km_ORIVX_X.push_back(Km_PF.vx());
            Km_ORIVX_Y.push_back(Km_PF.vy());
            Km_ORIVX_Z.push_back(Km_PF.vz());
            Km_P.push_back(Km_PF.p());
            Km_PT.push_back(Km_PF.pt());
            Km_PX.push_back(Km_PF.px());
            Km_PY.push_back(Km_PF.py());
            Km_PZ.push_back(Km_PF.pz());
            math::XYZPoint Km_originvertex(Km_ORIVX_X[0], Km_ORIVX_Y[0], Km_ORIVX_Z[0]);

            pi_ETA.push_back(pi_PF.eta());
            pi_PHI.push_back(pi_PF.phi());
            pi_ORIVX_X.push_back(pi_PF.vx());
            pi_ORIVX_Y.push_back(pi_PF.vy());
            pi_ORIVX_Z.push_back(pi_PF.vz());
            pi_P.push_back(pi_PF.p());
            pi_PT.push_back(pi_PF.pt());
            pi_PX.push_back(pi_PF.px());
            pi_PY.push_back(pi_PF.py());
            pi_PZ.push_back(pi_PF.pz());
            math::XYZPoint pi_originvertex(pi_ORIVX_X[0], pi_ORIVX_Y[0], pi_ORIVX_Z[0]);

            dR_Kp_Km.push_back(reco::deltaR(Kp_ETA[0], Kp_PHI[0], Km_ETA[0], Km_PHI[0]));
            dR_Kp_pi.push_back(reco::deltaR(Kp_ETA[0], Kp_PHI[0], pi_ETA[0], pi_PHI[0]));
            dR_Km_pi.push_back(reco::deltaR(Km_ETA[0], Km_PHI[0], pi_ETA[0], pi_PHI[0]));
            d_Kp_Km.push_back((Kp_originvertex - Km_originvertex).R());
            d_Kp_pi.push_back((Kp_originvertex - pi_originvertex).R());
            d_Km_pi.push_back((Km_originvertex - pi_originvertex).R());

            KalmanVertexFitter fitter(true);

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*theB).build(Kp_PF.pseudoTrack()),
                (*theB).build(Km_PF.pseudoTrack())
            };

            TransientVertex phi_Vertex = fitter.vertex(phi_Tracks);
            
            if( phi_Vertex.isValid() && phi_Vertex.hasRefittedTracks() ){

                std::vector<reco::TransientTrack> refitted_phi_Tracks = phi_Vertex.refittedTracks();
                TLorentzVector p4Kp; p4Kp.SetXYZM(refitted_phi_Tracks[0].track().px(), refitted_phi_Tracks[0].track().py(), refitted_phi_Tracks[0].track().pz(), Mass_K);
                TLorentzVector p4Km; p4Km.SetXYZM(refitted_phi_Tracks[1].track().px(), refitted_phi_Tracks[1].track().py(), refitted_phi_Tracks[1].track().pz(), Mass_K);
                TLorentzVector p4phi = p4Kp + p4Km;

                phi_ETA.push_back(p4phi.Eta());
                phi_PHI.push_back(p4phi.Phi());
                phi_P.push_back(p4phi.P());
                phi_PT.push_back(p4phi.Pt());
                phi_PX.push_back(p4phi.Px());
                phi_PY.push_back(p4phi.Py());
                phi_PZ.push_back(p4phi.Pz());
                phi_CHI2.push_back(phi_Vertex.totalChiSquared());
                phi_M.push_back(p4phi.M());
                phi_ENDVX_X.push_back(phi_Vertex.position().x());
                phi_ENDVX_Y.push_back(phi_Vertex.position().y());
                phi_ENDVX_Z.push_back(phi_Vertex.position().z());
                math::XYZPoint phi_endvertex(phi_ENDVX_X[0], phi_ENDVX_Y[0], phi_ENDVX_Z[0]);

                Kp_PP.push_back(p4Kp.Vect().Pt(p4phi.Vect()));
                Kp_PL.push_back(p4Kp.Vect().Dot(p4phi.Vect())/phi_P[0]);
                Km_PP.push_back(p4Km.Vect().Pt(p4phi.Vect()));
                Km_PL.push_back(p4Km.Vect().Dot(p4phi.Vect())/phi_P[0]);

                dR_Kp_phi.push_back(reco::deltaR(Kp_ETA[0], Kp_PHI[0], phi_ETA[0], phi_PHI[0]));
                dR_Km_phi.push_back(reco::deltaR(Km_ETA[0], Km_PHI[0], phi_ETA[0], phi_PHI[0]));
                dR_pi_phi.push_back(reco::deltaR(pi_ETA[0], pi_PHI[0], phi_ETA[0], phi_PHI[0]));
                d_Kp_phi.push_back((Kp_originvertex - phi_endvertex).R());
                d_Km_phi.push_back((Km_originvertex - phi_endvertex).R());
                d_pi_phi.push_back((pi_originvertex - phi_endvertex).R());

                std::vector<reco::TransientTrack> Ds_Tracks = {
                    refitted_phi_Tracks[0],
                    refitted_phi_Tracks[1],
                    (*theB).build(pi_PF.pseudoTrack())
                };

                TransientVertex Ds_Vertex = fitter.vertex(Ds_Tracks);

                if( Ds_Vertex.isValid() && Ds_Vertex.hasRefittedTracks() ){
                    std::vector<reco::TransientTrack> refitted_Ds_Tracks = Ds_Vertex.refittedTracks();
                    TLorentzVector refitted_p4Kp; refitted_p4Kp.SetXYZM(refitted_Ds_Tracks[0].track().px(), refitted_Ds_Tracks[0].track().py(), refitted_Ds_Tracks[0].track().pz(), Mass_K);
                    TLorentzVector refitted_p4Km; refitted_p4Km.SetXYZM(refitted_Ds_Tracks[1].track().px(), refitted_Ds_Tracks[1].track().py(), refitted_Ds_Tracks[1].track().pz(), Mass_K);
                    TLorentzVector refitted_p4phi = refitted_p4Kp + refitted_p4Km;
                    TLorentzVector p4pi; p4pi.SetXYZM(refitted_Ds_Tracks[2].track().px(), refitted_Ds_Tracks[2].track().py(), refitted_Ds_Tracks[2].track().pz(), Mass_pi);
                    TLorentzVector p4Ds = refitted_p4phi + p4pi;

                    Ds_ETA.push_back(p4Ds.Eta());
                    Ds_PHI.push_back(p4Ds.Phi());
                    Ds_P.push_back(p4Ds.P());
                    Ds_PT.push_back(p4Ds.Pt());
                    Ds_PX.push_back(p4Ds.Px());
                    Ds_PY.push_back(p4Ds.Py());
                    Ds_PZ.push_back(p4Ds.Pz());
                    Ds_CHI2.push_back(Ds_Vertex.totalChiSquared());
                    Ds_M.push_back(p4Ds.M());
                    Ds_ENDVX_X.push_back(Ds_Vertex.position().x());
                    Ds_ENDVX_Y.push_back(Ds_Vertex.position().y());
                    Ds_ENDVX_Z.push_back(Ds_Vertex.position().z());
                    math::XYZPoint Ds_endvertex(Ds_ENDVX_X[0], Ds_ENDVX_Y[0], Ds_ENDVX_Z[0]);

                    pi_PP.push_back(p4pi.Vect().Pt(p4Ds.Vect()));
                    pi_PL.push_back(p4pi.Vect().Dot(p4Ds.Vect())/Ds_P[0]);
                    phi_PP.push_back(refitted_p4phi.Vect().Pt(p4Ds.Vect()));
                    phi_PL.push_back(refitted_p4phi.Vect().Dot(p4Ds.Vect())/Ds_P[0]);

                    dR_Kp_Ds.push_back(reco::deltaR(Kp_ETA[0], Kp_PHI[0], Ds_ETA[0], Ds_PHI[0]));
                    dR_Km_Ds.push_back(reco::deltaR(Km_ETA[0], Km_PHI[0], Ds_ETA[0], Ds_PHI[0]));
                    dR_pi_Ds.push_back(reco::deltaR(pi_ETA[0], pi_PHI[0], Ds_ETA[0], Ds_PHI[0]));
                    dR_phi_Ds.push_back(reco::deltaR(phi_ETA[0], phi_PHI[0], Ds_ETA[0], Ds_PHI[0]));
                    d_Kp_Ds.push_back((Kp_originvertex - Ds_endvertex).R());
                    d_Km_Ds.push_back((Km_originvertex - Ds_endvertex).R());
                    d_pi_Ds.push_back((pi_originvertex - Ds_endvertex).R());
                    d_phi_Ds.push_back((phi_endvertex - Ds_endvertex).R());
                }
            }
        }
    } 

    tree_->Fill();
    return;
}

bool GenTrackMatch::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
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

void GenTrackMatch::vector_clear()
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
    // Gen particles info
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    num_Gen_phi = 0;
    num_Gen_Ds = 0;

    Gen_Kp_ETA.clear();
    Gen_Kp_PHI.clear();
    Gen_Kp_ORIVX_X.clear();
    Gen_Kp_ORIVX_Y.clear();
    Gen_Kp_ORIVX_Z.clear();
    Gen_Kp_P.clear();
    Gen_Kp_PT.clear();
    Gen_Kp_PX.clear();
    Gen_Kp_PY.clear();
    Gen_Kp_PZ.clear();
    Gen_Kp_PP.clear();
    Gen_Kp_PL.clear();

    Gen_Km_ETA.clear();
    Gen_Km_PHI.clear();
    Gen_Km_ORIVX_X.clear();
    Gen_Km_ORIVX_Y.clear();
    Gen_Km_ORIVX_Z.clear();
    Gen_Km_P.clear();
    Gen_Km_PT.clear();
    Gen_Km_PX.clear();
    Gen_Km_PY.clear();
    Gen_Km_PZ.clear();
    Gen_Km_PP.clear();
    Gen_Km_PL.clear();

    Gen_pi_ETA.clear();
    Gen_pi_PHI.clear();
    Gen_pi_ORIVX_X.clear();
    Gen_pi_ORIVX_Y.clear();
    Gen_pi_ORIVX_Z.clear();
    Gen_pi_P.clear();
    Gen_pi_PT.clear();
    Gen_pi_PX.clear();
    Gen_pi_PY.clear();
    Gen_pi_PZ.clear();
    Gen_pi_PP.clear();
    Gen_pi_PL.clear();

    Gen_phi_ETA.clear();
    Gen_phi_PHI.clear();
    Gen_phi_ORIVX_X.clear();
    Gen_phi_ORIVX_Y.clear();
    Gen_phi_ORIVX_Z.clear();
    Gen_phi_P.clear();
    Gen_phi_PT.clear();
    Gen_phi_PX.clear();
    Gen_phi_PY.clear();
    Gen_phi_PZ.clear();
    Gen_phi_PP.clear();
    Gen_phi_PL.clear();

    Gen_Ds_ETA.clear();
    Gen_Ds_PHI.clear();
    Gen_Ds_ORIVX_X.clear();
    Gen_Ds_ORIVX_Y.clear();
    Gen_Ds_ORIVX_Z.clear();
    Gen_Ds_P.clear();
    Gen_Ds_PT.clear();
    Gen_Ds_PX.clear();
    Gen_Ds_PY.clear();
    Gen_Ds_PZ.clear();

    Gen_dR_Kp_Km.clear();
    Gen_dR_Kp_pi.clear();
    Gen_dR_Kp_phi.clear();
    Gen_dR_Kp_Ds.clear();
    Gen_dR_Km_pi.clear();
    Gen_dR_Km_phi.clear();
    Gen_dR_Km_Ds.clear();
    Gen_dR_pi_phi.clear();
    Gen_dR_pi_Ds.clear();
    Gen_dR_phi_Ds.clear();

    Gen_d_Kp_Km.clear();
    Gen_d_Kp_pi.clear();
    Gen_d_Km_pi.clear();

    // matched particles info
    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;
    Kp_ETA.clear();
    Kp_PHI.clear();
    Kp_ORIVX_X.clear();
    Kp_ORIVX_Y.clear();
    Kp_ORIVX_Z.clear();
    Kp_P.clear();
    Kp_PT.clear();
    Kp_PX.clear();
    Kp_PY.clear();
    Kp_PZ.clear();
    Kp_PP.clear();
    Kp_PL.clear();

    Km_ETA.clear();
    Km_PHI.clear();
    Km_ORIVX_X.clear();
    Km_ORIVX_Y.clear();
    Km_ORIVX_Z.clear();
    Km_P.clear();
    Km_PT.clear();
    Km_PX.clear();
    Km_PY.clear();
    Km_PZ.clear();
    Km_PP.clear();
    Km_PL.clear();

    pi_ETA.clear();
    pi_PHI.clear();
    pi_ORIVX_X.clear();
    pi_ORIVX_Y.clear();
    pi_ORIVX_Z.clear();
    pi_P.clear();
    pi_PT.clear();
    pi_PX.clear();
    pi_PY.clear();
    pi_PZ.clear();
    pi_PP.clear();
    pi_PL.clear();

    phi_ETA.clear();
    phi_PHI.clear();
    phi_P.clear();
    phi_PT.clear();
    phi_PX.clear();
    phi_PY.clear();
    phi_PZ.clear();
    phi_PP.clear();
    phi_PL.clear();
    phi_CHI2.clear();
    phi_M.clear();
    phi_ENDVX_X.clear();
    phi_ENDVX_Y.clear();
    phi_ENDVX_Z.clear();

    Ds_ETA.clear();
    Ds_PHI.clear();
    Ds_P.clear();
    Ds_PT.clear();
    Ds_PX.clear();
    Ds_PY.clear();
    Ds_PZ.clear();
    Ds_CHI2.clear();
    Ds_M.clear();
    Ds_ENDVX_X.clear();
    Ds_ENDVX_Y.clear();
    Ds_ENDVX_Z.clear();

    dR_Kp_Km.clear();
    dR_Kp_pi.clear();
    dR_Kp_phi.clear();
    dR_Kp_Ds.clear();
    dR_Km_pi.clear();
    dR_Km_phi.clear();
    dR_Km_Ds.clear();
    dR_pi_phi.clear();
    dR_pi_Ds.clear();
    dR_phi_Ds.clear();

    d_Kp_Km.clear();
    d_Kp_pi.clear();
    d_Kp_phi.clear();
    d_Kp_Ds.clear();
    d_Km_pi.clear();
    d_Km_phi.clear();
    d_Km_Ds.clear();
    d_pi_phi.clear();
    d_pi_Ds.clear();
    d_phi_Ds.clear();
}

void GenTrackMatch::beginJob()
{
    tree_ = fs->make<TTree>("Events", "Events");

    tree_->Branch("num_Gen_Kp", &num_Gen_Kp);
    tree_->Branch("num_Gen_Km", &num_Gen_Km);
    tree_->Branch("num_Gen_pi", &num_Gen_pi);
    tree_->Branch("num_Gen_phi", &num_Gen_phi);
    tree_->Branch("num_Gen_Ds", &num_Gen_Ds);

    tree_->Branch("Gen_Kp_ETA", &Gen_Kp_ETA);
    tree_->Branch("Gen_Kp_PHI", &Gen_Kp_PHI);
    tree_->Branch("Gen_Kp_ORIVX_X", &Gen_Kp_ORIVX_X);
    tree_->Branch("Gen_Kp_ORIVX_Y", &Gen_Kp_ORIVX_Y);
    tree_->Branch("Gen_Kp_ORIVX_Z", &Gen_Kp_ORIVX_Z);
    tree_->Branch("Gen_Kp_P", &Gen_Kp_P);
    tree_->Branch("Gen_Kp_PT", &Gen_Kp_PT);
    tree_->Branch("Gen_Kp_PX", &Gen_Kp_PX);
    tree_->Branch("Gen_Kp_PY", &Gen_Kp_PY);
    tree_->Branch("Gen_Kp_PZ", &Gen_Kp_PZ);
    tree_->Branch("Gen_Kp_PP", &Gen_Kp_PP);
    tree_->Branch("Gen_Kp_PL", &Gen_Kp_PL);

    tree_->Branch("Gen_Km_ETA", &Gen_Km_ETA);
    tree_->Branch("Gen_Km_PHI", &Gen_Km_PHI);
    tree_->Branch("Gen_Km_ORIVX_X", &Gen_Km_ORIVX_X);
    tree_->Branch("Gen_Km_ORIVX_Y", &Gen_Km_ORIVX_Y);
    tree_->Branch("Gen_Km_ORIVX_Z", &Gen_Km_ORIVX_Z);
    tree_->Branch("Gen_Km_P", &Gen_Km_P);
    tree_->Branch("Gen_Km_PT", &Gen_Km_PT);
    tree_->Branch("Gen_Km_PX", &Gen_Km_PX);
    tree_->Branch("Gen_Km_PY", &Gen_Km_PY);
    tree_->Branch("Gen_Km_PZ", &Gen_Km_PZ);
    tree_->Branch("Gen_Km_PP", &Gen_Km_PP);
    tree_->Branch("Gen_Km_PL", &Gen_Km_PL);

    tree_->Branch("Gen_pi_ETA", &Gen_pi_ETA);
    tree_->Branch("Gen_pi_PHI", &Gen_pi_PHI);
    tree_->Branch("Gen_pi_ORIVX_X", &Gen_pi_ORIVX_X);
    tree_->Branch("Gen_pi_ORIVX_Y", &Gen_pi_ORIVX_Y);
    tree_->Branch("Gen_pi_ORIVX_Z", &Gen_pi_ORIVX_Z);
    tree_->Branch("Gen_pi_P", &Gen_pi_P);
    tree_->Branch("Gen_pi_PT", &Gen_pi_PT);
    tree_->Branch("Gen_pi_PX", &Gen_pi_PX);
    tree_->Branch("Gen_pi_PY", &Gen_pi_PY);
    tree_->Branch("Gen_pi_PZ", &Gen_pi_PZ);
    tree_->Branch("Gen_pi_PP", &Gen_pi_PP);
    tree_->Branch("Gen_pi_PL", &Gen_pi_PL);

    tree_->Branch("Gen_phi_ETA", &Gen_phi_ETA);
    tree_->Branch("Gen_phi_PHI", &Gen_phi_PHI);
    tree_->Branch("Gen_phi_ORIVX_X", &Gen_phi_ORIVX_X);
    tree_->Branch("Gen_phi_ORIVX_Y", &Gen_phi_ORIVX_Y);
    tree_->Branch("Gen_phi_ORIVX_Z", &Gen_phi_ORIVX_Z);
    tree_->Branch("Gen_phi_P", &Gen_phi_P);
    tree_->Branch("Gen_phi_PT", &Gen_phi_PT);
    tree_->Branch("Gen_phi_PX", &Gen_phi_PX);
    tree_->Branch("Gen_phi_PY", &Gen_phi_PY);
    tree_->Branch("Gen_phi_PZ", &Gen_phi_PZ);
    tree_->Branch("Gen_phi_PP", &Gen_phi_PP);
    tree_->Branch("Gen_phi_PL", &Gen_phi_PL);

    tree_->Branch("Gen_Ds_ETA", &Gen_Ds_ETA);
    tree_->Branch("Gen_Ds_PHI", &Gen_Ds_PHI);
    tree_->Branch("Gen_Ds_ORIVX_X", &Gen_Ds_ORIVX_X);
    tree_->Branch("Gen_Ds_ORIVX_Y", &Gen_Ds_ORIVX_Y);
    tree_->Branch("Gen_Ds_ORIVX_Z", &Gen_Ds_ORIVX_Z);
    tree_->Branch("Gen_Ds_P", &Gen_Ds_P);
    tree_->Branch("Gen_Ds_PT", &Gen_Ds_PT);
    tree_->Branch("Gen_Ds_PX", &Gen_Ds_PX);
    tree_->Branch("Gen_Ds_PY", &Gen_Ds_PY);
    tree_->Branch("Gen_Ds_PZ", &Gen_Ds_PZ);

    tree_->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
    tree_->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
    tree_->Branch("Gen_dR_Kp_phi", &Gen_dR_Kp_phi);
    tree_->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds);
    tree_->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);
    tree_->Branch("Gen_dR_Km_phi", &Gen_dR_Km_phi);
    tree_->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds);
    tree_->Branch("Gen_dR_pi_phi", &Gen_dR_pi_phi);
    tree_->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds);
    tree_->Branch("Gen_dR_phi_Ds", &Gen_dR_phi_Ds);

    tree_->Branch("Gen_d_Kp_Km", &Gen_d_Kp_Km);
    tree_->Branch("Gen_d_Kp_pi", &Gen_d_Kp_pi);
    tree_->Branch("Gen_d_Km_pi", &Gen_d_Km_pi);

    tree_->Branch("num_match_Kp", &num_match_Kp);
    tree_->Branch("num_match_Km", &num_match_Km);
    tree_->Branch("num_match_pi", &num_match_pi);
    tree_->Branch("Kp_ETA", &Kp_ETA);
    tree_->Branch("Kp_PHI", &Kp_PHI);
    tree_->Branch("Kp_ORIVX_X", &Kp_ORIVX_X);
    tree_->Branch("Kp_ORIVX_Y", &Kp_ORIVX_Y);
    tree_->Branch("Kp_ORIVX_Z", &Kp_ORIVX_Z);
    tree_->Branch("Kp_P", &Kp_P);
    tree_->Branch("Kp_PT", &Kp_PT);
    tree_->Branch("Kp_PX", &Kp_PX);
    tree_->Branch("Kp_PY", &Kp_PY);
    tree_->Branch("Kp_PZ", &Kp_PZ);
    tree_->Branch("Kp_PP", &Kp_PP);
    tree_->Branch("Kp_PL", &Kp_PL);

    tree_->Branch("Km_ETA", &Km_ETA);
    tree_->Branch("Km_PHI", &Km_PHI);
    tree_->Branch("Km_ORIVX_X", &Km_ORIVX_X);
    tree_->Branch("Km_ORIVX_Y", &Km_ORIVX_Y);
    tree_->Branch("Km_ORIVX_Z", &Km_ORIVX_Z);
    tree_->Branch("Km_P", &Km_P);
    tree_->Branch("Km_PT", &Km_PT);
    tree_->Branch("Km_PX", &Km_PX);
    tree_->Branch("Km_PY", &Km_PY);
    tree_->Branch("Km_PZ", &Km_PZ);
    tree_->Branch("Km_PP", &Km_PP);
    tree_->Branch("Km_PL", &Km_PL);

    tree_->Branch("pi_ETA", &pi_ETA);
    tree_->Branch("pi_PHI", &pi_PHI);
    tree_->Branch("pi_ORIVX_X", &pi_ORIVX_X);
    tree_->Branch("pi_ORIVX_Y", &pi_ORIVX_Y);
    tree_->Branch("pi_ORIVX_Z", &pi_ORIVX_Z);
    tree_->Branch("pi_P", &pi_P);
    tree_->Branch("pi_PT", &pi_PT);
    tree_->Branch("pi_PX", &pi_PX);
    tree_->Branch("pi_PY", &pi_PY);
    tree_->Branch("pi_PZ", &pi_PZ);
    tree_->Branch("pi_PP", &pi_PP);
    tree_->Branch("pi_PL", &pi_PL);

    tree_->Branch("phi_ETA", &phi_ETA);
    tree_->Branch("phi_PHI", &phi_PHI);
    tree_->Branch("phi_P", &phi_P);
    tree_->Branch("phi_PT", &phi_PT);
    tree_->Branch("phi_PX", &phi_PX);
    tree_->Branch("phi_PY", &phi_PY);
    tree_->Branch("phi_PZ", &phi_PZ);
    tree_->Branch("phi_PP", &phi_PP);
    tree_->Branch("phi_PL", &phi_PL);
    tree_->Branch("phi_CHI2", &phi_CHI2);
    tree_->Branch("phi_M", &phi_M);
    tree_->Branch("phi_ENDVX_X", &phi_ENDVX_X);
    tree_->Branch("phi_ENDVX_Y", &phi_ENDVX_Y);
    tree_->Branch("phi_ENDVX_Z", &phi_ENDVX_Z);

    tree_->Branch("Ds_ETA", &Ds_ETA);
    tree_->Branch("Ds_PHI", &Ds_PHI);
    tree_->Branch("Ds_P", &Ds_P);
    tree_->Branch("Ds_PT", &Ds_PT);
    tree_->Branch("Ds_PX", &Ds_PX);
    tree_->Branch("Ds_PY", &Ds_PY);
    tree_->Branch("Ds_PZ", &Ds_PZ);
    tree_->Branch("Ds_CHI2", &Ds_CHI2);
    tree_->Branch("Ds_M", &Ds_M);
    tree_->Branch("Ds_ENDVX_X", &Ds_ENDVX_X);
    tree_->Branch("Ds_ENDVX_Y", &Ds_ENDVX_Y);
    tree_->Branch("Ds_ENDVX_Z", &Ds_ENDVX_Z);

    tree_->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree_->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree_->Branch("dR_Kp_phi", &dR_Kp_phi);
    tree_->Branch("dR_Kp_Ds", &dR_Kp_Ds);
    tree_->Branch("dR_Km_pi", &dR_Km_pi);
    tree_->Branch("dR_Km_phi", &dR_Km_phi);
    tree_->Branch("dR_Km_Ds", &dR_Km_Ds);
    tree_->Branch("dR_pi_phi", &dR_pi_phi);
    tree_->Branch("dR_pi_Ds", &dR_pi_Ds);
    tree_->Branch("dR_phi_Ds", &dR_phi_Ds);

    tree_->Branch("d_Kp_Km", &d_Kp_Km);
    tree_->Branch("d_Kp_pi", &d_Kp_pi);
    tree_->Branch("d_Kp_phi", &d_Kp_phi);
    tree_->Branch("d_Kp_Ds", &d_Kp_Ds);
    tree_->Branch("d_Km_pi", &d_Km_pi);
    tree_->Branch("d_Km_phi", &d_Km_phi);
    tree_->Branch("d_Km_Ds", &d_Km_Ds);
    tree_->Branch("d_pi_phi", &d_pi_phi);
    tree_->Branch("d_pi_Ds", &d_pi_Ds);
    tree_->Branch("d_phi_Ds", &d_phi_Ds);
}

void GenTrackMatch::endJob() {}

void GenTrackMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("GenTrackMatch", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenTrackMatch);

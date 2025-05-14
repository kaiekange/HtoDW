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

#include "EDAnalyzers/GenParticleAnalyzer/interface/GenTree.h"


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

        const edm::Service<TFileService> fs;
        GenTree *ftree;

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

GenTrackMatch::GenTrackMatch(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
    ftree = new GenTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

GenTrackMatch::~GenTrackMatch() {}

void GenTrackMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    if(ftree->num_Gen_Kp==1 && ftree->num_Gen_Km==1 && ftree->num_Gen_pi==1 && ftree->num_Gen_phi==1 && ftree->num_Gen_Ds==1){

        const auto& Kp_GP = (*prunedGen)[Gen_Kp_idx[0]];
        const auto& Km_GP = (*prunedGen)[Gen_Km_idx[0]];
        const auto& pi_GP = (*prunedGen)[Gen_pi_idx[0]];
        const auto& phi_GP = (*prunedGen)[Gen_phi_idx[0]];
        const auto& Ds_GP = (*prunedGen)[Gen_Ds_idx[0]];

        ftree->Gen_Kp_ETA.push_back(Kp_GP.eta());
        ftree->Gen_Kp_PHI.push_back(Kp_GP.phi());
        ftree->Gen_Kp_ORIVX_X.push_back(Kp_GP.vx());
        ftree->Gen_Kp_ORIVX_Y.push_back(Kp_GP.vy());
        ftree->Gen_Kp_ORIVX_Z.push_back(Kp_GP.vz());
        ftree->Gen_Kp_P.push_back(Kp_GP.p());
        ftree->Gen_Kp_PT.push_back(Kp_GP.pt());
        ftree->Gen_Kp_PX.push_back(Kp_GP.px());
        ftree->Gen_Kp_PY.push_back(Kp_GP.py());
        ftree->Gen_Kp_PZ.push_back(Kp_GP.pz());
        TVector3 Gen_Kp_P3(Kp_GP.px(), Kp_GP.py(), Kp_GP.pz());

        ftree->Gen_Km_ETA.push_back(Km_GP.eta());
        ftree->Gen_Km_PHI.push_back(Km_GP.phi());
        ftree->Gen_Km_ORIVX_X.push_back(Km_GP.vx());
        ftree->Gen_Km_ORIVX_Y.push_back(Km_GP.vy());
        ftree->Gen_Km_ORIVX_Z.push_back(Km_GP.vz());
        ftree->Gen_Km_P.push_back(Km_GP.p());
        ftree->Gen_Km_PT.push_back(Km_GP.pt());
        ftree->Gen_Km_PX.push_back(Km_GP.px());
        ftree->Gen_Km_PY.push_back(Km_GP.py());
        ftree->Gen_Km_PZ.push_back(Km_GP.pz());
        TVector3 Gen_Km_P3(Km_GP.px(), Km_GP.py(), Km_GP.pz());

        ftree->Gen_pi_ETA.push_back(pi_GP.eta());
        ftree->Gen_pi_PHI.push_back(pi_GP.phi());
        ftree->Gen_pi_ORIVX_X.push_back(pi_GP.vx());
        ftree->Gen_pi_ORIVX_Y.push_back(pi_GP.vy());
        ftree->Gen_pi_ORIVX_Z.push_back(pi_GP.vz());
        ftree->Gen_pi_P.push_back(pi_GP.p());
        ftree->Gen_pi_PT.push_back(pi_GP.pt());
        ftree->Gen_pi_PX.push_back(pi_GP.px());
        ftree->Gen_pi_PY.push_back(pi_GP.py());
        ftree->Gen_pi_PZ.push_back(pi_GP.pz());
        TVector3 Gen_pi_P3(pi_GP.px(), pi_GP.py(), pi_GP.pz());

        ftree->Gen_phi_ETA.push_back(phi_GP.eta());
        ftree->Gen_phi_PHI.push_back(phi_GP.phi());
        ftree->Gen_phi_ORIVX_X.push_back(phi_GP.vx());
        ftree->Gen_phi_ORIVX_Y.push_back(phi_GP.vy());
        ftree->Gen_phi_ORIVX_Z.push_back(phi_GP.vz());
        ftree->Gen_phi_P.push_back(phi_GP.p());
        ftree->Gen_phi_PT.push_back(phi_GP.pt());
        ftree->Gen_phi_PX.push_back(phi_GP.px());
        ftree->Gen_phi_PY.push_back(phi_GP.py());
        ftree->Gen_phi_PZ.push_back(phi_GP.pz());
        TVector3 Gen_phi_P3(phi_GP.px(), phi_GP.py(), phi_GP.pz());

        ftree->Gen_Ds_ETA.push_back(Ds_GP.eta());
        ftree->Gen_Ds_PHI.push_back(Ds_GP.phi());
        ftree->Gen_Ds_ORIVX_X.push_back(Ds_GP.vx());
        ftree->Gen_Ds_ORIVX_Y.push_back(Ds_GP.vy());
        ftree->Gen_Ds_ORIVX_Z.push_back(Ds_GP.vz());
        ftree->Gen_Ds_P.push_back(Ds_GP.p());
        ftree->Gen_Ds_PT.push_back(Ds_GP.pt());
        ftree->Gen_Ds_PX.push_back(Ds_GP.px());
        ftree->Gen_Ds_PY.push_back(Ds_GP.py());
        ftree->Gen_Ds_PZ.push_back(Ds_GP.pz());
        TVector3 Gen_Ds_P3(Ds_GP.px(), Ds_GP.py(), Ds_GP.pz());

        ftree->Gen_Kp_PP.push_back(Gen_Kp_P3.Pt(Gen_phi_P3));
        ftree->Gen_Kp_PL.push_back(Gen_Kp_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag());
        ftree->Gen_Km_PP.push_back(Gen_Km_P3.Pt(Gen_phi_P3));
        ftree->Gen_Km_PL.push_back(Gen_Km_P3.Dot(Gen_phi_P3)/Gen_phi_P3.Mag());
        ftree->Gen_pi_PP.push_back(Gen_pi_P3.Pt(Gen_Ds_P3));
        ftree->Gen_pi_PL.push_back(Gen_pi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag());
        ftree->Gen_phi_PP.push_back(Gen_phi_P3.Pt(Gen_Ds_P3));
        ftree->Gen_phi_PL.push_back(Gen_phi_P3.Dot(Gen_Ds_P3)/Gen_Ds_P3.Mag());
        
        ftree->Gen_dR_Kp_Km.push_back(reco::deltaR(Kp_GP, Km_GP));
        ftree->Gen_dR_Kp_phi.push_back(reco::deltaR(Kp_GP, phi_GP));
        ftree->Gen_dR_Km_phi.push_back(reco::deltaR(Km_GP, phi_GP));
        ftree->Gen_dR_Kp_pi.push_back(reco::deltaR(Kp_GP, pi_GP));
        ftree->Gen_dR_Km_pi.push_back(reco::deltaR(Km_GP, pi_GP));
        ftree->Gen_dR_phi_pi.push_back(reco::deltaR(phi_GP, pi_GP));
        ftree->Gen_dR_Kp_Ds.push_back(reco::deltaR(Kp_GP, Ds_GP));
        ftree->Gen_dR_Km_Ds.push_back(reco::deltaR(Km_GP, Ds_GP));
        ftree->Gen_dR_phi_Ds.push_back(reco::deltaR(phi_GP, Ds_GP));
        ftree->Gen_dR_pi_Ds.push_back(reco::deltaR(pi_GP, Ds_GP));

        ftree->Gen_d_Kp_Km.push_back((Kp_GP.vertex() - Km_GP.vertex()).R());
        ftree->Gen_d_Kp_pi.push_back((Kp_GP.vertex() - pi_GP.vertex()).R());
        ftree->Gen_d_Km_pi.push_back((Km_GP.vertex() - pi_GP.vertex()).R());

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

        if( bestKp.dr < 0.03 && bestKm.dr < 0.03 && bestpi.dr < 0.03 && (*packedPF)[bestKp.index].hasTrackDetails() && (*packedPF)[bestKm.index].hasTrackDetails() && (*packedPF)[bestpi.index].hasTrackDetails() ){

            const auto& Kp_PF = (*packedPF)[bestKp.index];
            const auto& Km_PF = (*packedPF)[bestKm.index];
            const auto& pi_PF = (*packedPF)[bestpi.index];

            ftree->Kp_ETA.push_back(Kp_PF.eta());
            ftree->Kp_PHI.push_back(Kp_PF.phi());
            ftree->Kp_ORIVX_X.push_back(Kp_PF.vx());
            ftree->Kp_ORIVX_Y.push_back(Kp_PF.vy());
            ftree->Kp_ORIVX_Z.push_back(Kp_PF.vz());
            ftree->Kp_P.push_back(Kp_PF.p());
            ftree->Kp_PT.push_back(Kp_PF.pt());
            ftree->Kp_PX.push_back(Kp_PF.px());
            ftree->Kp_PY.push_back(Kp_PF.py());
            ftree->Kp_PZ.push_back(Kp_PF.pz());
            math::XYZPoint Kp_originvertex(Kp_PF.vx(), Kp_PF.vy(), Kp_PF.vz());
            TLorentzVector Kp_P4; Kp_P4.SetXYZM(Kp_PF.px(), Kp_PF.py(), Kp_PF.pz(), Mass_K);

            ftree->Km_ETA.push_back(Km_PF.eta());
            ftree->Km_PHI.push_back(Km_PF.phi());
            ftree->Km_ORIVX_X.push_back(Km_PF.vx());
            ftree->Km_ORIVX_Y.push_back(Km_PF.vy());
            ftree->Km_ORIVX_Z.push_back(Km_PF.vz());
            ftree->Km_P.push_back(Km_PF.p());
            ftree->Km_PT.push_back(Km_PF.pt());
            ftree->Km_PX.push_back(Km_PF.px());
            ftree->Km_PY.push_back(Km_PF.py());
            ftree->Km_PZ.push_back(Km_PF.pz());
            math::XYZPoint Km_originvertex(Km_PF.vx(), Km_PF.vy(), Km_PF.vz());
            TLorentzVector Km_P4; Km_P4.SetXYZM(Km_PF.px(), Km_PF.py(), Km_PF.pz(), Mass_K);

            ftree->pi_ETA.push_back(pi_PF.eta());
            ftree->pi_PHI.push_back(pi_PF.phi());
            ftree->pi_ORIVX_X.push_back(pi_PF.vx());
            ftree->pi_ORIVX_Y.push_back(pi_PF.vy());
            ftree->pi_ORIVX_Z.push_back(pi_PF.vz());
            ftree->pi_P.push_back(pi_PF.p());
            ftree->pi_PT.push_back(pi_PF.pt());
            ftree->pi_PX.push_back(pi_PF.px());
            ftree->pi_PY.push_back(pi_PF.py());
            ftree->pi_PZ.push_back(pi_PF.pz());
            math::XYZPoint pi_originvertex(pi_PF.vx(), pi_PF.vy(), pi_PF.vz());
            TLorentzVector pi_P4; pi_P4.SetXYZM(pi_PF.px(), pi_PF.py(), pi_PF.pz(), Mass_pi);

            ftree->dR_Kp_Km.push_back(reco::deltaR(Kp_PF, Km_PF));
            ftree->dR_Kp_pi.push_back(reco::deltaR(Kp_PF, pi_PF));
            ftree->dR_Km_pi.push_back(reco::deltaR(Km_PF, pi_PF));
            ftree->d_Kp_Km.push_back((Kp_originvertex - Km_originvertex).R());
            ftree->d_Kp_pi.push_back((Kp_originvertex - pi_originvertex).R());
            ftree->d_Km_pi.push_back((Km_originvertex - pi_originvertex).R());

            KalmanVertexFitter fitter(true);

            std::vector<reco::TransientTrack> phi_Tracks = {
                (*theB).build(Kp_PF.pseudoTrack()),
                (*theB).build(Km_PF.pseudoTrack())
            };

            TransientVertex phi_Vertex = fitter.vertex(phi_Tracks);
            
            if( phi_Vertex.isValid() && phi_Vertex.hasRefittedTracks() ){

                std::vector<reco::TransientTrack> phiFit_Tracks = phi_Vertex.refittedTracks();
                TLorentzVector phiFit_Kp_P4;
                phiFit_Kp_P4.SetXYZM(phiFit_Tracks[0].track().px(), phiFit_Tracks[0].track().py(), phiFit_Tracks[0].track().pz(), Mass_K);
                ftree->phiFit_Kp_ETA.push_back(phiFit_Kp_P4.Eta());
                ftree->phiFit_Kp_PHI.push_back(phiFit_Kp_P4.Phi());
                ftree->phiFit_Kp_P.push_back(phiFit_Kp_P4.P());
                ftree->phiFit_Kp_PT.push_back(phiFit_Kp_P4.Pt());
                ftree->phiFit_Kp_PX.push_back(phiFit_Kp_P4.Px());
                ftree->phiFit_Kp_PY.push_back(phiFit_Kp_P4.Py());
                ftree->phiFit_Kp_PZ.push_back(phiFit_Kp_P4.Pz());

                TLorentzVector phiFit_Km_P4;
                phiFit_Km_P4.SetXYZM(phiFit_Tracks[1].track().px(), phiFit_Tracks[1].track().py(), phiFit_Tracks[1].track().pz(), Mass_K);
                ftree->phiFit_Km_ETA.push_back(phiFit_Km_P4.Eta());
                ftree->phiFit_Km_PHI.push_back(phiFit_Km_P4.Phi());
                ftree->phiFit_Km_P.push_back(phiFit_Km_P4.P());
                ftree->phiFit_Km_PT.push_back(phiFit_Km_P4.Pt());
                ftree->phiFit_Km_PX.push_back(phiFit_Km_P4.Px());
                ftree->phiFit_Km_PY.push_back(phiFit_Km_P4.Py());
                ftree->phiFit_Km_PZ.push_back(phiFit_Km_P4.Pz());

                TLorentzVector phi_P4 = phiFit_Kp_P4 + phiFit_Km_P4;
                ftree->phi_ETA.push_back(phi_P4.Eta());
                ftree->phi_PHI.push_back(phi_P4.Phi());
                ftree->phi_P.push_back(phi_P4.P());
                ftree->phi_PT.push_back(phi_P4.Pt());
                ftree->phi_PX.push_back(phi_P4.Px());
                ftree->phi_PY.push_back(phi_P4.Py());
                ftree->phi_PZ.push_back(phi_P4.Pz());
                ftree->phi_CHI2.push_back(phi_Vertex.totalChiSquared());
                ftree->phi_M.push_back(phi_P4.M());
                ftree->phi_ENDVX_X.push_back(phi_Vertex.position().x());
                ftree->phi_ENDVX_Y.push_back(phi_Vertex.position().y());
                ftree->phi_ENDVX_Z.push_back(phi_Vertex.position().z());
                math::XYZPoint phi_endvertex(phi_Vertex.position().x(), phi_Vertex.position().y(), phi_Vertex.position().z());

                ftree->Kp_PP.push_back(Kp_P4.Vect().Pt(phi_P4.Vect()));
                ftree->Kp_PL.push_back(Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P());
                ftree->Km_PP.push_back(Km_P4.Vect().Pt(phi_P4.Vect()));
                ftree->Km_PL.push_back(Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P());

                ftree->phiFit_Kp_PP.push_back(phiFit_Kp_P4.Vect().Pt(phi_P4.Vect()));
                ftree->phiFit_Kp_PL.push_back(phiFit_Kp_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P());
                ftree->phiFit_Km_PP.push_back(phiFit_Km_P4.Vect().Pt(phi_P4.Vect()));
                ftree->phiFit_Km_PL.push_back(phiFit_Km_P4.Vect().Dot(phi_P4.Vect())/phi_P4.P());

                ftree->dR_Kp_phi.push_back(reco::deltaR(ftree->Kp_ETA[0], ftree->Kp_PHI[0], ftree->phi_ETA[0], ftree->phi_PHI[0]));
                ftree->dR_Km_phi.push_back(reco::deltaR(ftree->Km_ETA[0], ftree->Km_PHI[0], ftree->phi_ETA[0], ftree->phi_PHI[0]));
                ftree->dR_phi_pi.push_back(reco::deltaR(ftree->phi_ETA[0], ftree->phi_PHI[0], ftree->pi_ETA[0], ftree->pi_PHI[0]));
                ftree->d_Kp_phi.push_back((Kp_originvertex - phi_endvertex).R());
                ftree->d_Km_phi.push_back((Km_originvertex - phi_endvertex).R());
                ftree->d_phi_pi.push_back((phi_endvertex - pi_originvertex).R());

                std::vector<reco::TransientTrack> Ds_Tracks = {
                    phiFit_Tracks[0],
                    phiFit_Tracks[1],
                    (*theB).build(pi_PF.pseudoTrack())
                };

                TransientVertex Ds_Vertex = fitter.vertex(Ds_Tracks);

                if( Ds_Vertex.isValid() && Ds_Vertex.hasRefittedTracks() ){
                    std::vector<reco::TransientTrack> DsFit_Tracks = Ds_Vertex.refittedTracks();

                    TLorentzVector DsFit_Kp_P4;
                    DsFit_Kp_P4.SetXYZM(DsFit_Tracks[0].track().px(), DsFit_Tracks[0].track().py(), DsFit_Tracks[0].track().pz(), Mass_K);
                    ftree->DsFit_Kp_ETA.push_back(DsFit_Kp_P4.Eta());
                    ftree->DsFit_Kp_PHI.push_back(DsFit_Kp_P4.Phi());
                    ftree->DsFit_Kp_P.push_back(DsFit_Kp_P4.P());
                    ftree->DsFit_Kp_PT.push_back(DsFit_Kp_P4.Pt());
                    ftree->DsFit_Kp_PX.push_back(DsFit_Kp_P4.Px());
                    ftree->DsFit_Kp_PY.push_back(DsFit_Kp_P4.Py());
                    ftree->DsFit_Kp_PZ.push_back(DsFit_Kp_P4.Pz());


                    TLorentzVector DsFit_Km_P4;
                    DsFit_Km_P4.SetXYZM(DsFit_Tracks[1].track().px(), DsFit_Tracks[1].track().py(), DsFit_Tracks[1].track().pz(), Mass_K);
                    ftree->DsFit_Km_ETA.push_back(DsFit_Km_P4.Eta());
                    ftree->DsFit_Km_PHI.push_back(DsFit_Km_P4.Phi());
                    ftree->DsFit_Km_P.push_back(DsFit_Km_P4.P());
                    ftree->DsFit_Km_PT.push_back(DsFit_Km_P4.Pt());
                    ftree->DsFit_Km_PX.push_back(DsFit_Km_P4.Px());
                    ftree->DsFit_Km_PY.push_back(DsFit_Km_P4.Py());
                    ftree->DsFit_Km_PZ.push_back(DsFit_Km_P4.Pz());

                    TLorentzVector DsFit_pi_P4;
                    DsFit_pi_P4.SetXYZM(DsFit_Tracks[2].track().px(), DsFit_Tracks[2].track().py(), DsFit_Tracks[2].track().pz(), Mass_pi);
                    ftree->DsFit_pi_ETA.push_back(DsFit_pi_P4.Eta());
                    ftree->DsFit_pi_PHI.push_back(DsFit_pi_P4.Phi());
                    ftree->DsFit_pi_P.push_back(DsFit_pi_P4.P());
                    ftree->DsFit_pi_PT.push_back(DsFit_pi_P4.Pt());
                    ftree->DsFit_pi_PX.push_back(DsFit_pi_P4.Px());
                    ftree->DsFit_pi_PY.push_back(DsFit_pi_P4.Py());
                    ftree->DsFit_pi_PZ.push_back(DsFit_pi_P4.Pz());

                    TLorentzVector DsFit_phi_P4 = DsFit_Kp_P4 + DsFit_Km_P4;
                    ftree->DsFit_phi_ETA.push_back(DsFit_phi_P4.Eta());
                    ftree->DsFit_phi_PHI.push_back(DsFit_phi_P4.Phi());
                    ftree->DsFit_phi_P.push_back(DsFit_phi_P4.P());
                    ftree->DsFit_phi_PT.push_back(DsFit_phi_P4.Pt());
                    ftree->DsFit_phi_PX.push_back(DsFit_phi_P4.Px());
                    ftree->DsFit_phi_PY.push_back(DsFit_phi_P4.Py());
                    ftree->DsFit_phi_PZ.push_back(DsFit_phi_P4.Pz());
                    ftree->DsFit_phi_M.push_back(DsFit_phi_P4.M());

                    TLorentzVector Ds_P4 = DsFit_phi_P4 + DsFit_pi_P4;
                    ftree->Ds_ETA.push_back(Ds_P4.Eta());
                    ftree->Ds_PHI.push_back(Ds_P4.Phi());
                    ftree->Ds_P.push_back(Ds_P4.P());
                    ftree->Ds_PT.push_back(Ds_P4.Pt());
                    ftree->Ds_PX.push_back(Ds_P4.Px());
                    ftree->Ds_PY.push_back(Ds_P4.Py());
                    ftree->Ds_PZ.push_back(Ds_P4.Pz());
                    ftree->Ds_CHI2.push_back(Ds_Vertex.totalChiSquared());
                    ftree->Ds_M.push_back(Ds_P4.M());
                    ftree->Ds_ENDVX_X.push_back(Ds_Vertex.position().x());
                    ftree->Ds_ENDVX_Y.push_back(Ds_Vertex.position().y());
                    ftree->Ds_ENDVX_Z.push_back(Ds_Vertex.position().z());
                    math::XYZPoint Ds_endvertex(Ds_Vertex.position().x(), Ds_Vertex.position().y(), Ds_Vertex.position().z());

                    TLorentzVector DsFit_Mconstraint_phi_P4;
                    DsFit_Mconstraint_phi_P4.SetXYZM(DsFit_phi_P4.Px(), DsFit_phi_P4.Py(), DsFit_phi_P4.Pz(), Mass_phi);
                    TLorentzVector DsFit_Mconstraint_Ds_P4 = DsFit_Mconstraint_phi_P4 + DsFit_pi_P4;
                    ftree->DsFit_Mconstraint_Ds_M.push_back(DsFit_Mconstraint_Ds_P4.M());
                    
                    ftree->pi_PP.push_back(pi_P4.Vect().Pt(Ds_P4.Vect()));
                    ftree->pi_PL.push_back(pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P());
                    ftree->phi_PP.push_back(phi_P4.Vect().Pt(Ds_P4.Vect()));
                    ftree->phi_PL.push_back(phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P());
                
                    ftree->DsFit_pi_PP.push_back(DsFit_pi_P4.Vect().Pt(Ds_P4.Vect()));
                    ftree->DsFit_pi_PL.push_back(DsFit_pi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P());
                    ftree->DsFit_phi_PP.push_back(DsFit_phi_P4.Vect().Pt(Ds_P4.Vect()));
                    ftree->DsFit_phi_PL.push_back(DsFit_phi_P4.Vect().Dot(Ds_P4.Vect())/Ds_P4.P());
                
                    ftree->DsFit_Kp_PP.push_back(DsFit_Kp_P4.Vect().Pt(DsFit_phi_P4.Vect()));
                    ftree->DsFit_Kp_PL.push_back(DsFit_Kp_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P());
                    ftree->DsFit_Km_PP.push_back(DsFit_Km_P4.Vect().Pt(DsFit_phi_P4.Vect()));
                    ftree->DsFit_Km_PL.push_back(DsFit_Km_P4.Vect().Dot(DsFit_phi_P4.Vect())/DsFit_phi_P4.P());

                    ftree->dR_Kp_Ds.push_back(reco::deltaR(ftree->Kp_ETA[0], ftree->Kp_PHI[0], ftree->Ds_ETA[0], ftree->Ds_PHI[0]));
                    ftree->dR_Km_Ds.push_back(reco::deltaR(ftree->Km_ETA[0], ftree->Km_PHI[0], ftree->Ds_ETA[0], ftree->Ds_PHI[0]));
                    ftree->dR_phi_Ds.push_back(reco::deltaR(ftree->phi_ETA[0], ftree->phi_PHI[0], ftree->Ds_ETA[0], ftree->Ds_PHI[0]));
                    ftree->dR_pi_Ds.push_back(reco::deltaR(ftree->pi_ETA[0], ftree->pi_PHI[0], ftree->Ds_ETA[0], ftree->Ds_PHI[0]));
                    ftree->d_Kp_Ds.push_back((Kp_originvertex - Ds_endvertex).R());
                    ftree->d_Km_Ds.push_back((Km_originvertex - Ds_endvertex).R());
                    ftree->d_phi_Ds.push_back((phi_endvertex - Ds_endvertex).R());
                    ftree->d_pi_Ds.push_back((pi_originvertex - Ds_endvertex).R());
                }
            }
        }
    } 

    ftree->tree->Fill();
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

void GenTrackMatch::beginJob() {}

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

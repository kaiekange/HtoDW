// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      GenPacked
//
/**\class GenPacked GenPacked.cc EDAnalyzer/GenParticleAnalyzer/plugins/GenPacked.cc

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

class GenPacked : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit GenPacked(const edm::ParameterSet&);
        ~GenPacked();

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

        std::vector<double> all_pt;
        std::vector<double> all_p;

        // dR of hadrons to Gen
        std::vector<float> all_dR_Kp;
        std::vector<float> all_dR_Km;
        std::vector<float> all_dR_pi;
        std::vector<float> dR_Kp;
        std::vector<float> dR_Km;
        std::vector<float> dR_pi;

        // Gen particles info
        int num_Gen_Kp;
        int num_Gen_Km;
        int num_Gen_pi;
        std::vector<double> Gen_p_Kp;
        std::vector<double> Gen_p_Km;
        std::vector<double> Gen_p_pi;
        std::vector<double> Gen_pt_Kp;
        std::vector<double> Gen_pt_Km;
        std::vector<double> Gen_pt_pi;
        std::vector<double> Gen_eta_Kp;
        std::vector<double> Gen_eta_Km;
        std::vector<double> Gen_eta_pi;
        std::vector<double> Gen_phi_Kp;
        std::vector<double> Gen_phi_Km;
        std::vector<double> Gen_phi_pi;

        // dR & d betweem Gen tracks
        std::vector<double> Gen_dR_Kp_Km;
        std::vector<double> Gen_dR_Kp_pi;
        std::vector<double> Gen_dR_Km_pi;

        // matched particles info
        int num_match_Kp;
        int num_match_Km;
        int num_match_pi;
        std::vector<double> p_Kp;
        std::vector<double> p_Km;
        std::vector<double> p_pi;
        std::vector<double> pt_Kp;
        std::vector<double> pt_Km;
        std::vector<double> pt_pi;
        std::vector<double> eta_Kp;
        std::vector<double> eta_Km;
        std::vector<double> eta_pi;
        std::vector<double> phi_Kp;
        std::vector<double> phi_Km;
        std::vector<double> phi_pi;

        // dR & d betweem matched tracks
        std::vector<double> dR_Kp_Km;
        std::vector<double> dR_Kp_pi;
        std::vector<double> dR_Km_pi;
        std::vector<double> d_Kp_Km;
        std::vector<double> d_Kp_pi;
        std::vector<double> d_Km_pi;

        // phi variables
        std::vector<double> chi2_phi;
        std::vector<double> pt_phi;
        std::vector<double> p_phi;
        std::vector<double> m_phi;
        std::vector<double> refitted_pt_phi;
        std::vector<double> refitted_p_phi;
        std::vector<double> refitted_m_phi;
        std::vector<double> dR_phi_pi;
        std::vector<double> d_phi_pi;

        // Ds variable
        std::vector<double> chi2_Ds;
        std::vector<double> pt_Ds;
        std::vector<double> p_Ds;
        std::vector<double> m_Ds;
        std::vector<double> refitted_pt_Ds;
        std::vector<double> refitted_p_Ds;
        std::vector<double> refitted_m_Ds;
        std::vector<double> refitted2_pt_Ds;
        std::vector<double> refitted2_p_Ds;
        std::vector<double> refitted2_m_Ds;
        std::vector<double> pt_Ds_phimass_constraint;
        std::vector<double> p_Ds_phimass_constraint;
        std::vector<double> m_Ds_phimass_constraint;
        std::vector<double> dR_phi_Ds;
        std::vector<double> d_phi_Ds;

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

GenPacked::GenPacked(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
}

GenPacked::~GenPacked() {}

void GenPacked::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    for(const auto& gp : *prunedGen){
        const int pdg = gp.pdgId();
        if( pdg==321 && hasAncestor(gp, 333, -321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_eta_Kp.push_back(gp.eta());
            Gen_phi_Kp.push_back(gp.phi());
            Gen_pt_Kp.push_back(gp.pt());
            Gen_p_Kp.push_back(gp.p());
            num_Gen_Kp++;
        } else if( pdg==-321 && hasAncestor(gp, 333, 321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24) ){
            Gen_eta_Km.push_back(gp.eta());
            Gen_phi_Km.push_back(gp.phi());
            Gen_pt_Km.push_back(gp.pt());
            Gen_p_Km.push_back(gp.p());
            num_Gen_Km++;
        } else if( pdg==211 && hasAncestor(gp, 431, 333) && hasAncestor(gp, 25, -24) ){
            Gen_eta_pi.push_back(gp.eta());
            Gen_phi_pi.push_back(gp.phi());
            Gen_pt_pi.push_back(gp.pt());
            Gen_p_pi.push_back(gp.p());
            num_Gen_pi++;
        }
    }

    if(num_Gen_Kp==1 && num_Gen_Km==1 && num_Gen_pi==1){

        Gen_dR_Kp_Km.push_back(reco::deltaR(Gen_eta_Kp[0], Gen_phi_Kp[0], Gen_eta_Km[0], Gen_phi_Km[0]));
        Gen_dR_Kp_pi.push_back(reco::deltaR(Gen_eta_Kp[0], Gen_phi_Kp[0], Gen_eta_pi[0], Gen_phi_pi[0]));
        Gen_dR_Km_pi.push_back(reco::deltaR(Gen_eta_Km[0], Gen_phi_Km[0], Gen_eta_pi[0], Gen_phi_pi[0]));

        struct MatchInfo {int index; float dr;};
        std::vector<MatchInfo> MatchKp, MatchKm, Matchpi;

        for(size_t i=0; i<packedPF->size(); i++){
            const auto& pf = (*packedPF)[i];

            const float dr_Kp = reco::deltaR(pf.eta(), pf.phi(), Gen_eta_Kp[0], Gen_phi_Kp[0]);
            const float dr_Km = reco::deltaR(pf.eta(), pf.phi(), Gen_eta_Km[0], Gen_phi_Km[0]);
            const float dr_pi = reco::deltaR(pf.eta(), pf.phi(), Gen_eta_pi[0], Gen_phi_pi[0]);

            all_dR_Kp.push_back(dr_Kp);
            all_dR_Km.push_back(dr_Km);
            all_dR_pi.push_back(dr_pi);
            all_pt.push_back(pf.pt());
            all_p.push_back(pf.p());

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

            eta_Kp.push_back((*packedPF)[bestKp.index].eta());
            phi_Kp.push_back((*packedPF)[bestKp.index].phi());
            pt_Kp.push_back((*packedPF)[bestKp.index].pt());
            p_Kp.push_back((*packedPF)[bestKp.index].p());

            eta_Km.push_back((*packedPF)[bestKm.index].eta());
            phi_Km.push_back((*packedPF)[bestKm.index].phi());
            pt_Km.push_back((*packedPF)[bestKm.index].pt());
            p_Km.push_back((*packedPF)[bestKm.index].p());

            eta_pi.push_back((*packedPF)[bestpi.index].eta());
            phi_pi.push_back((*packedPF)[bestpi.index].phi());
            pt_pi.push_back((*packedPF)[bestpi.index].pt());
            p_pi.push_back((*packedPF)[bestpi.index].p());

            dR_Kp_Km.push_back(reco::deltaR((*packedPF)[bestKp.index], (*packedPF)[bestKm.index]));
            dR_Kp_pi.push_back(reco::deltaR((*packedPF)[bestKp.index], (*packedPF)[bestpi.index]));
            dR_Km_pi.push_back(reco::deltaR((*packedPF)[bestKm.index], (*packedPF)[bestpi.index]));

            d_Kp_Km.push_back(((*packedPF)[bestKp.index].vertex() - (*packedPF)[bestKm.index].vertex()).R());
            d_Kp_pi.push_back(((*packedPF)[bestKp.index].vertex() - (*packedPF)[bestpi.index].vertex()).R());
            d_Km_pi.push_back(((*packedPF)[bestKm.index].vertex() - (*packedPF)[bestpi.index].vertex()).R());

            KalmanVertexFitter fitter(true);

            std::vector<reco::TransientTrack> phiTracks = {
                (*theB).build((*packedPF)[bestKp.index].pseudoTrack()),
                (*theB).build((*packedPF)[bestKm.index].pseudoTrack())
            };

            TransientVertex phiVertex = fitter.vertex(phiTracks);

            if( phiVertex.isValid() && phiVertex.hasRefittedTracks() ){
                chi2_phi.push_back(phiVertex.totalChiSquared());

                std::vector<reco::TransientTrack> refitted_phiTracks = phiVertex.refittedTracks();

                TLorentzVector p4Kp; p4Kp.SetXYZM(phiTracks[0].track().px(), phiTracks[0].track().py(), phiTracks[0].track().pz(), Mass_K);
                TLorentzVector p4Km; p4Km.SetXYZM(phiTracks[1].track().px(), phiTracks[1].track().py(), phiTracks[1].track().pz(), Mass_K);
                TLorentzVector p4phi; p4phi = p4Kp + p4Km;

                pt_phi.push_back(p4phi.Pt());
                p_phi.push_back(p4phi.P());
                m_phi.push_back(p4phi.M());

                TLorentzVector refitted_p4Kp; refitted_p4Kp.SetXYZM(refitted_phiTracks[0].track().px(), refitted_phiTracks[0].track().py(), refitted_phiTracks[0].track().pz(), Mass_K);
                TLorentzVector refitted_p4Km; refitted_p4Km.SetXYZM(refitted_phiTracks[1].track().px(), refitted_phiTracks[1].track().py(), refitted_phiTracks[1].track().pz(), Mass_K);
                TLorentzVector refitted_p4phi; refitted_p4phi = refitted_p4Kp + refitted_p4Km;

                refitted_pt_phi.push_back(refitted_p4phi.Pt());
                refitted_p_phi.push_back(refitted_p4phi.P());
                refitted_m_phi.push_back(refitted_p4phi.M());

                dR_phi_pi.push_back(reco::deltaR(refitted_p4phi.Eta(), refitted_p4phi.Phi(), (*packedPF)[bestpi.index].eta(), (*packedPF)[bestpi.index].phi()));
                math::XYZPoint vtx_phi(phiVertex.position().x(), phiVertex.position().y(), phiVertex.position().z());
                math::XYZPoint vtx_pi = (*packedPF)[bestpi.index].vertex();
                d_phi_pi.push_back((vtx_phi - vtx_pi).R());

                std::vector<reco::TransientTrack> DsTracks = {
                    refitted_phiTracks[0],
                    refitted_phiTracks[1],
                    (*theB).build((*packedPF)[bestpi.index].pseudoTrack())
                };

                TransientVertex DsVertex = fitter.vertex(DsTracks);

                if( DsVertex.isValid() && DsVertex.hasRefittedTracks() ){

                    chi2_Ds.push_back(DsVertex.totalChiSquared());

                    std::vector<reco::TransientTrack> refitted_DsTracks = DsVertex.refittedTracks();

                    TLorentzVector p4pi; p4pi.SetXYZM(DsTracks[2].track().px(), DsTracks[2].track().py(), DsTracks[2].track().pz(), Mass_pi);
                    TLorentzVector p4Ds; p4Ds = p4pi + p4phi;
                    pt_Ds.push_back(p4Ds.Pt());
                    p_Ds.push_back(p4Ds.P());
                    m_Ds.push_back(p4Ds.M());

                    TLorentzVector refitted_p4Ds; refitted_p4Ds = p4pi + refitted_p4phi;
                    refitted_pt_Ds.push_back(refitted_p4Ds.Pt());
                    refitted_p_Ds.push_back(refitted_p4Ds.P());
                    refitted_m_Ds.push_back(refitted_p4Ds.M());

                    TLorentzVector refitted_p4pi; refitted_p4pi.SetXYZM(refitted_DsTracks[2].track().px(), refitted_DsTracks[2].track().py(), refitted_DsTracks[2].track().pz(), Mass_pi);
                    TLorentzVector refitted2_p4Kp; refitted2_p4Kp.SetXYZM(refitted_DsTracks[0].track().px(), refitted_DsTracks[0].track().py(), refitted_DsTracks[0].track().pz(), Mass_K);
                    TLorentzVector refitted2_p4Km; refitted2_p4Km.SetXYZM(refitted_DsTracks[1].track().px(), refitted_DsTracks[1].track().py(), refitted_DsTracks[1].track().pz(), Mass_K);
                    TLorentzVector refitted2_p4Ds; refitted2_p4Ds = refitted_p4pi + refitted2_p4Kp + refitted2_p4Km;
                    refitted2_pt_Ds.push_back(refitted2_p4Ds.Pt());
                    refitted2_p_Ds.push_back(refitted2_p4Ds.P());
                    refitted2_m_Ds.push_back(refitted2_p4Ds.M());

                    dR_phi_Ds.push_back(reco::deltaR(refitted_p4phi.Eta(), refitted_p4phi.Phi(), refitted2_p4Ds.Eta(), refitted2_p4Ds.Phi()));
                    math::XYZPoint vtx_Ds(DsVertex.position().x(), DsVertex.position().y(), DsVertex.position().z());
                    d_phi_Ds.push_back((vtx_phi - vtx_Ds).R());

                    TLorentzVector p4phi_mass_constraint; p4phi_mass_constraint.SetXYZM(refitted_p4phi.Px(), refitted_p4phi.Py(), refitted_p4phi.Pz(), Mass_phi);
                    TLorentzVector p4Ds_phimass_constraint; p4Ds_phimass_constraint = refitted_p4pi + p4phi_mass_constraint;

                    pt_Ds_phimass_constraint.push_back(p4Ds_phimass_constraint.Pt());
                    p_Ds_phimass_constraint.push_back(p4Ds_phimass_constraint.P());
                    m_Ds_phimass_constraint.push_back(p4Ds_phimass_constraint.M());

                }
            }
        }
    } 

    tree_->Fill();
    return;
}

bool GenPacked::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
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

void GenPacked::vector_clear()
{
    all_pt.clear();
    all_p.clear();

    // dR of hadrons to Gen
    all_dR_Kp.clear();
    all_dR_Km.clear();
    all_dR_pi.clear();
    dR_Kp.clear();
    dR_Km.clear();
    dR_pi.clear();

    // Gen particles info
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    Gen_p_Kp.clear();
    Gen_p_Km.clear();
    Gen_p_pi.clear();
    Gen_pt_Kp.clear();
    Gen_pt_Km.clear();
    Gen_pt_pi.clear();
    Gen_eta_Kp.clear();
    Gen_eta_Km.clear();
    Gen_eta_pi.clear();
    Gen_phi_Kp.clear();
    Gen_phi_Km.clear();
    Gen_phi_pi.clear();

    // dR & d betweem Gen tracks
    Gen_dR_Kp_Km.clear();
    Gen_dR_Kp_pi.clear();
    Gen_dR_Km_pi.clear();

    // matched particles info
    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;
    p_Kp.clear();
    p_Km.clear();
    p_pi.clear();
    pt_Kp.clear();
    pt_Km.clear();
    pt_pi.clear();
    eta_Kp.clear();
    eta_Km.clear();
    eta_pi.clear();
    phi_Kp.clear();
    phi_Km.clear();
    phi_pi.clear();

    // dR & d betweem matched tracks
    dR_Kp_Km.clear();
    dR_Kp_pi.clear();
    dR_Km_pi.clear();
    d_Kp_Km.clear();
    d_Kp_pi.clear();
    d_Km_pi.clear();

    // phi variables
    chi2_phi.clear();
    pt_phi.clear();
    p_phi.clear();
    m_phi.clear();
    refitted_pt_phi.clear();
    refitted_p_phi.clear();
    refitted_m_phi.clear();
    dR_phi_pi.clear();
    d_phi_pi.clear();

    // Ds variable
    chi2_Ds.clear();
    pt_Ds.clear();
    p_Ds.clear();
    m_Ds.clear();
    refitted_pt_Ds.clear();
    refitted_p_Ds.clear();
    refitted_m_Ds.clear();
    refitted2_pt_Ds.clear();
    refitted2_p_Ds.clear();
    refitted2_m_Ds.clear();
    pt_Ds_phimass_constraint.clear();
    p_Ds_phimass_constraint.clear();
    m_Ds_phimass_constraint.clear();
    dR_phi_Ds.clear();
    d_phi_Ds.clear();
}

void GenPacked::beginJob()
{
    tree_ = fs->make<TTree>("Events", "Events");

    tree_->Branch("all_pt", &all_pt);
    tree_->Branch("all_p", &all_p);

    tree_->Branch("all_dR_Kp", &all_dR_Kp);
    tree_->Branch("all_dR_Km", &all_dR_Km);
    tree_->Branch("all_dR_pi", &all_dR_pi);
    tree_->Branch("dR_Kp", &dR_Kp);
    tree_->Branch("dR_Km", &dR_Km);
    tree_->Branch("dR_pi", &dR_pi);

    tree_->Branch("num_Gen_Kp", &num_Gen_Kp);
    tree_->Branch("num_Gen_Km", &num_Gen_Km);
    tree_->Branch("num_Gen_pi", &num_Gen_pi);
    tree_->Branch("Gen_p_Kp", &Gen_p_Kp);
    tree_->Branch("Gen_p_Km", &Gen_p_Km);
    tree_->Branch("Gen_p_pi", &Gen_p_pi);
    tree_->Branch("Gen_pt_Kp", &Gen_pt_Kp);
    tree_->Branch("Gen_pt_Km", &Gen_pt_Km);
    tree_->Branch("Gen_pt_pi", &Gen_pt_pi);
    tree_->Branch("Gen_eta_Kp", &Gen_eta_Kp);
    tree_->Branch("Gen_eta_Km", &Gen_eta_Km);
    tree_->Branch("Gen_eta_pi", &Gen_eta_pi);
    tree_->Branch("Gen_phi_Kp", &Gen_phi_Kp);
    tree_->Branch("Gen_phi_Km", &Gen_phi_Km);
    tree_->Branch("Gen_phi_pi", &Gen_phi_pi);

    tree_->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
    tree_->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
    tree_->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);

    tree_->Branch("num_match_Kp", &num_match_Kp);
    tree_->Branch("num_match_Km", &num_match_Km);
    tree_->Branch("num_match_pi", &num_match_pi);
    tree_->Branch("p_Kp", &p_Kp);
    tree_->Branch("p_Km", &p_Km);
    tree_->Branch("p_pi", &p_pi);
    tree_->Branch("pt_Kp", &pt_Kp);
    tree_->Branch("pt_Km", &pt_Km);
    tree_->Branch("pt_pi", &pt_pi);
    tree_->Branch("eta_Kp", &eta_Kp);
    tree_->Branch("eta_Km", &eta_Km);
    tree_->Branch("eta_pi", &eta_pi);
    tree_->Branch("phi_Kp", &phi_Kp);
    tree_->Branch("phi_Km", &phi_Km);
    tree_->Branch("phi_pi", &phi_pi);

    tree_->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree_->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree_->Branch("dR_Km_pi", &dR_Km_pi);
    tree_->Branch("d_Kp_Km", &d_Kp_Km);
    tree_->Branch("d_Kp_pi", &d_Kp_pi);
    tree_->Branch("d_Km_pi", &d_Km_pi);

    tree_->Branch("chi2_phi", &chi2_phi);
    tree_->Branch("pt_phi", &pt_phi);
    tree_->Branch("p_phi", &p_phi);
    tree_->Branch("m_phi", &m_phi);
    tree_->Branch("refitted_pt_phi", &refitted_pt_phi);
    tree_->Branch("refitted_p_phi", &refitted_p_phi);
    tree_->Branch("refitted_m_phi", &refitted_m_phi);
    tree_->Branch("dR_phi_pi", &dR_phi_pi);
    tree_->Branch("d_phi_pi", &d_phi_pi);

    tree_->Branch("chi2_Ds", &chi2_Ds);
    tree_->Branch("pt_Ds", &pt_Ds);
    tree_->Branch("p_Ds", &p_Ds);
    tree_->Branch("m_Ds", &m_Ds);
    tree_->Branch("refitted_pt_Ds", &refitted_pt_Ds);
    tree_->Branch("refitted_p_Ds", &refitted_p_Ds);
    tree_->Branch("refitted_m_Ds", &refitted_m_Ds);
    tree_->Branch("refitted2_pt_Ds", &refitted2_pt_Ds);
    tree_->Branch("refitted2_p_Ds", &refitted2_p_Ds);
    tree_->Branch("refitted2_m_Ds", &refitted2_m_Ds);
    tree_->Branch("pt_Ds_phimass_constraint", &pt_Ds_phimass_constraint);
    tree_->Branch("p_Ds_phimass_constraint", &p_Ds_phimass_constraint);
    tree_->Branch("m_Ds_phimass_constraint", &m_Ds_phimass_constraint);
    tree_->Branch("dR_phi_Ds", &dR_phi_Ds);
    tree_->Branch("d_phi_Ds", &d_phi_Ds);
}

void GenPacked::endJob() {}

void GenPacked::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
    desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));;
    descriptions.add("GenPacked", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPacked);

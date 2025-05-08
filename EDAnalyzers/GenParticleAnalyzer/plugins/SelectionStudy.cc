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


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
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
        float getDeltaR(float eta1, float phi1, float eta2, float phi2);
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

        // matched particles momenta
        std::vector<double> p_Kp;
        std::vector<double> p_Km;
        std::vector<double> p_pi;
        std::vector<double> pt_Kp;
        std::vector<double> pt_Km;
        std::vector<double> pt_pi;

        // dR betweem tracks 
        std::vector<double> dR_Kp_Km;
        std::vector<double> dR_Kp_pi;
        std::vector<double> dR_Km_pi;
        std::vector<double> d_Kp_Km;
        std::vector<double> d_Kp_pi;
        std::vector<double> d_Km_pi;

        // fit variables
        std::vector<double> chi2_phi;
        std::vector<double> m_phi;
        std::vector<double> refitted_m_phi;
        std::vector<double> chi2_Ds;
        std::vector<double> m_Ds;
        std::vector<double> refitted_m_Ds;
        std::vector<double> refitted2_m_Ds;

        // dR of resonances
        std::vector<double> d_phi_pi;
        std::vector<double> dR_phi_pi;
        std::vector<double> d_phi_Ds;
        std::vector<double> dR_phi_Ds;

        const double Mass_K = 0.493677;
        const double Mass_pi = 0.13957039;
        const double Mass_phi = 1.019461;
};

SelectionStudy::SelectionStudy(const edm::ParameterSet& iConfig) :
    prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates")))
{
}

SelectionStudy::~SelectionStudy() {}

void SelectionStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    float phi_Kp=100;
    float eta_Kp=100;
    float phi_Km=100;
    float eta_Km=100;
    float phi_pi=100;
    float eta_pi=100;

    int nKp=0;
    int nKm=0;
    int npi=0; 
    
    for (const auto& gp : *prunedGen){
        const int pdg = gp.pdgId();

        if(pdg==321 && hasAncestor(gp, 333, -321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24)){
            eta_Kp = gp.eta();
            phi_Kp = gp.phi();
            nKp++; 
        }

        else if(pdg==-321 && hasAncestor(gp, 333, 321) && hasAncestor(gp, 431, 211) && hasAncestor(gp, 25, -24)){
            eta_Km = gp.eta();
            phi_Km = gp.phi();
            nKm++;
        }

        else if(pdg==211 && hasAncestor(gp, 431, 333) && hasAncestor(gp, 25, -24)){
            eta_pi = gp.eta();
            phi_pi = gp.phi();
            npi++;
        }
    }

    if(nKp!=1 || nKm!=1 || npi!=1) return;

    struct MatchInfo {int index; float dr;};
    std::vector<MatchInfo> MatchKp, MatchKm, Matchpi;

    for(size_t i=0; i<packedPF->size(); i++){
        const auto& pf = (*packedPF)[i];

        const float dr_Kp = reco::deltaR(pf.eta(), pf.phi(), eta_Kp, phi_Kp);  
        const float dr_Km = reco::deltaR(pf.eta(), pf.phi(), eta_Km, phi_Km);  
        const float dr_pi = reco::deltaR(pf.eta(), pf.phi(), eta_pi, phi_pi);

        all_dR_Kp.push_back(dr_Kp);
        all_dR_Km.push_back(dr_Km);
        all_dR_pi.push_back(dr_pi);

        if(pf.pdgId() == 211 && pf.trackHighPurity()){
            MatchKp.push_back({static_cast<int>(i), dr_Kp});
            Matchpi.push_back({static_cast<int>(i), dr_pi});
        } else if(pf.pdgId() == -211 && pf.trackHighPurity()){
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

    if(bestKp.dr > 0.03) return;
    if(bestKm.dr > 0.03) return;
    if(bestpi.dr > 0.03) return;

    if(!(*packedPF)[bestKp.index].hasTrackDetails()) return;
    if(!(*packedPF)[bestKm.index].hasTrackDetails()) return;
    if(!(*packedPF)[bestpi.index].hasTrackDetails()) return;

    pt_Kp.push_back((*packedPF)[bestKp.index].pt());
    pt_Km.push_back((*packedPF)[bestKm.index].pt());
    pt_pi.push_back((*packedPF)[bestpi.index].pt());

    p_Kp.push_back((*packedPF)[bestKp.index].p());
    p_Km.push_back((*packedPF)[bestKm.index].p());
    p_pi.push_back((*packedPF)[bestpi.index].p());

    d_Kp_Km.push_back(((*packedPF)[bestKp.index].vertex() - (*packedPF)[bestKm.index].vertex()).R());
    d_Kp_pi.push_back(((*packedPF)[bestKp.index].vertex() - (*packedPF)[bestpi.index].vertex()).R());
    d_Km_pi.push_back(((*packedPF)[bestKm.index].vertex() - (*packedPF)[bestpi.index].vertex()).R());

    dR_Kp_Km.push_back(reco::deltaR((*packedPF)[bestKp.index], (*packedPF)[bestKm.index]));
    dR_Kp_pi.push_back(reco::deltaR((*packedPF)[bestKp.index], (*packedPF)[bestpi.index]));
    dR_Km_pi.push_back(reco::deltaR((*packedPF)[bestKm.index], (*packedPF)[bestpi.index]));

    KalmanVertexFitter fitter(true);

    std::vector<reco::TransientTrack> phiTracks = {
        (*theB).build((*packedPF)[bestKp.index].pseudoTrack()),
        (*theB).build((*packedPF)[bestKm.index].pseudoTrack())
    };

    TransientVertex phiVertex = fitter.vertex(phiTracks);

    if(!(phiVertex.isValid())) return;
    if(!(phiVertex.hasRefittedTracks())) return;

    chi2_phi.push_back(phiVertex.totalChiSquared());

    std::vector<reco::TransientTrack> refitted_phiTracks = phiVertex.refittedTracks();

    TLorentzVector p4Kp, p4Km, refitted_p4Kp, refitted_p4Km, p4phi, refitted_p4phi;
    p4Kp.SetXYZM(phiTracks[0].track().px(), phiTracks[0].track().py(), phiTracks[0].track().pz(), Mass_K);
    p4Km.SetXYZM(phiTracks[1].track().px(), phiTracks[1].track().py(), phiTracks[1].track().pz(), Mass_K);
    p4phi = p4Kp + p4Km;
    refitted_p4Kp.SetXYZM(refitted_phiTracks[0].track().px(), refitted_phiTracks[0].track().py(), refitted_phiTracks[0].track().pz(), Mass_K);
    refitted_p4Km.SetXYZM(refitted_phiTracks[1].track().px(), refitted_phiTracks[1].track().py(), refitted_phiTracks[1].track().pz(), Mass_K);
    refitted_p4phi = refitted_p4Kp + refitted_p4Km;
    m_phi.push_back(p4phi.M());
    refitted_m_phi.push_back(refitted_p4phi.M());

    math::XYZPoint vtx_phi(phiVertex.position().x(), phiVertex.position().y(), phiVertex.position().z());
    math::XYZPoint vtx_pi = (*packedPF)[bestpi.index].vertex();
    d_phi_pi.push_back((vtx_phi - vtx_pi).R());
    dR_phi_pi.push_back(reco::deltaR(refitted_p4phi.Eta(), refitted_p4phi.Phi(), (*packedPF)[bestpi.index].eta(), (*packedPF)[bestpi.index].phi()));

    std::vector<reco::TransientTrack> DsTracks = {
        refitted_phiTracks[0],
        refitted_phiTracks[1],
        (*theB).build((*packedPF)[bestpi.index].pseudoTrack())
    };

    TransientVertex DsVertex = fitter.vertex(DsTracks);

    if(!(DsVertex.isValid())) return;
    if(!(DsVertex.hasRefittedTracks())) return;

    chi2_Ds.push_back(DsVertex.totalChiSquared());

    std::vector<reco::TransientTrack> refitted_DsTracks = DsVertex.refittedTracks();

    TLorentzVector p4pi, refitted_p4pi, refitted2_p4Kp, refitted2_p4Km, p4Ds, refitted_p4Ds, refitted2_p4Ds;
    p4pi.SetXYZM(DsTracks[2].track().px(), DsTracks[2].track().py(), DsTracks[2].track().pz(), Mass_pi);
    p4Ds = p4pi + p4phi;
    m_Ds.push_back(p4Ds.M());

    refitted_p4Ds = p4pi + refitted_p4phi;
    refitted_m_Ds.push_back(refitted_p4Ds.M());

    refitted_p4pi.SetXYZM(refitted_DsTracks[2].track().px(), refitted_DsTracks[2].track().py(), refitted_DsTracks[2].track().pz(), Mass_pi); 
    refitted2_p4Kp.SetXYZM(refitted_DsTracks[0].track().px(), refitted_DsTracks[0].track().py(), refitted_DsTracks[0].track().pz(), Mass_K);
    refitted2_p4Km.SetXYZM(refitted_DsTracks[1].track().px(), refitted_DsTracks[1].track().py(), refitted_DsTracks[1].track().pz(), Mass_K);
    refitted2_p4Ds = refitted_p4pi + refitted2_p4Kp + refitted2_p4Km; 
    refitted2_m_Ds.push_back(refitted2_p4Ds.M());

    math::XYZPoint vtx_Ds(DsVertex.position().x(), DsVertex.position().y(), DsVertex.position().z());
    d_phi_Ds.push_back((vtx_phi - vtx_Ds).R());
    dR_phi_Ds.push_back(reco::deltaR(refitted_p4phi.Eta(), refitted_p4phi.Phi(), refitted2_p4Ds.Eta(), refitted2_p4Ds.Phi()));

    tree_->Fill();
    return;
}

void SelectionStudy::vector_clear()
{
    all_dR_Kp.clear();
    all_dR_Km.clear();
    all_dR_pi.clear();
    dR_Kp.clear();
    dR_Km.clear();
    dR_pi.clear();

    p_Kp.clear();
    p_Km.clear();
    p_pi.clear();
    pt_Kp.clear();
    pt_Km.clear();
    pt_pi.clear();

    dR_Kp_Km.clear();
    dR_Kp_pi.clear();
    dR_Km_pi.clear();
    d_Kp_Km.clear();
    d_Kp_pi.clear();
    d_Km_pi.clear();

    chi2_phi.clear();
    m_phi.clear();
    refitted_m_phi.clear();
    chi2_Ds.clear();
    m_Ds.clear();
    refitted_m_Ds.clear();
    refitted2_m_Ds.clear();

    d_phi_pi.clear();
    dR_phi_pi.clear();
    d_phi_Ds.clear();
    dR_phi_Ds.clear();
}

float SelectionStudy::getDeltaR(float eta1, float phi1, float eta2, float phi2)
{      
    float DeltaPhi = fabs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

bool SelectionStudy::hasAncestor(const reco::GenParticle & particle, const int motherid, const int otherparticleid) const
{
    bool hasancestor = false;

    if(particle.numberOfMothers()==0) return false;

    if(particle.numberOfMothers()==1 && particle.mother(0)->pdgId()==motherid){
        if(particle.mother(0)->numberOfDaughters()==2){
            if(particle.mother(0)->daughter(0)->pdgId()==otherparticleid || particle.mother(0)->daughter(1)->pdgId()==otherparticleid) hasancestor = true;
        }
    }

    if (hasancestor) return true;
    else return hasAncestor(*particle.motherRef(0), motherid, otherparticleid);
}

void SelectionStudy::beginJob()
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("Events", "Events");

    tree_->Branch("p_Kp", &p_Kp);
    tree_->Branch("p_Km", &p_Km);
    tree_->Branch("p_pi", &p_pi);
    tree_->Branch("pt_Kp", &pt_Kp);
    tree_->Branch("pt_Km", &pt_Km);
    tree_->Branch("pt_pi", &pt_pi);

    tree_->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree_->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree_->Branch("dR_Km_pi", &dR_Km_pi);
    tree_->Branch("d_Kp_Km", &d_Kp_Km);
    tree_->Branch("d_Kp_pi", &d_Kp_pi);
    tree_->Branch("d_Km_pi", &d_Km_pi);

    tree_->Branch("chi2_phi", &chi2_phi);
    tree_->Branch("m_phi", &m_phi);
    tree_->Branch("refitted_m_phi", &refitted_m_phi);
    tree_->Branch("chi2_Ds", &chi2_Ds);
    tree_->Branch("m_Ds", &m_Ds);
    tree_->Branch("refitted_m_Ds", &refitted_m_Ds);
    tree_->Branch("refitted2_m_Ds", &refitted2_m_Ds);

    tree_->Branch("d_phi_pi", &d_phi_pi);
    tree_->Branch("dR_phi_pi", &dR_phi_pi);
    tree_->Branch("d_phi_Ds", &d_phi_Ds);
    tree_->Branch("dR_phi_Ds", &dR_phi_Ds);
}

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

// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      MuonStudy
//
/**\class MuonStudy MuonStudy.cc EDAnalyzer/GenParticleAnalyzer/plugins/MuonStudy.cc

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

#include "EDAnalyzers/RecoAnalyzer/interface/MuonTree.h"
#include "EDAnalyzers/RecoAnalyzer/interface/VertexReProducer.h"

class MuonStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit MuonStudy(const edm::ParameterSet& iConfig);
        ~MuonStudy();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> prunedToken;
        edm::EDGetTokenT<pat::MuonCollection> muonToken;
        edm::EDGetTokenT<reco::VertexCollection> pvToken;
        edm::EDGetTokenT<reco::BeamSpot> bsToken;

        const edm::Service<TFileService> fs;
        MuonTree *ftree;

        bool hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid ) const;

        std::vector<float> dR_mu;

        std::vector<int> Gen_mu_idx;
        std::vector<int> Gen_nu_idx;
        std::vector<int> Gen_W_idx;

        std::vector<int> idx_mu_vec;
};

MuonStudy::MuonStudy(const edm::ParameterSet& iConfig) :
    prunedToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPart"))),
    muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    pvToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primvtx"))),
    bsToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")))
{
}

MuonStudy::~MuonStudy() {}

void MuonStudy::beginJob()
{
    TFile & f = fs->file();
    f.SetCompressionAlgorithm(ROOT::kZLIB);
    ftree = new MuonTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();
}

void MuonStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    dR_mu.clear();

    Gen_mu_idx.clear();
    Gen_nu_idx.clear();
    Gen_W_idx.clear();

    ftree->Init();

    edm::Handle<reco::GenParticleCollection> prunedHandle;
    iEvent.getByToken(prunedToken, prunedHandle);

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
        if( pdg==13 && hasAncestor(gp, -24, -14) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_mu_idx.push_back(i);
            ftree->num_Gen_mu++;
        } else if( pdg==-14 && hasAncestor(gp, -24, 13) && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_nu_idx.push_back(i);
            ftree->num_Gen_nu++;
        } else if( pdg==-24 && hasAncestor(gp, 25, 431) && gp.isHardProcess() ){
            Gen_W_idx.push_back(i);
            ftree->num_Gen_W++;
        }
    }

    ftree->Gen_Reset();

    if(ftree->num_Gen_mu!=1) return;
    if(ftree->num_Gen_nu!=1) return;
    if(ftree->num_Gen_W!=1) return;

    const auto& mu_GP = (*prunedHandle)[Gen_mu_idx[0]];
    const auto& nu_GP = (*prunedHandle)[Gen_nu_idx[0]];
    const auto& W_GP = (*prunedHandle)[Gen_W_idx[0]];

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
    std::vector<MatchInfo> Matchmu;

    for(size_t i=0; i<muonHandle->size(); i++){
        const auto& muon = (*muonHandle)[i];
        if (!muon.genParticle()) continue;
        const reco::GenParticle& genMuon = *(muon.genParticle());
        if(hasAncestor(genMuon, -24, -14) && hasAncestor(genMuon, 25, 431)){
            float dr_mu = reco::deltaR(muon.eta(), muon.phi(), mu_GP.eta(), mu_GP.phi());
            dR_mu.push_back(dr_mu);
            Matchmu.push_back({static_cast<int>(i), dr_mu});
            ftree->num_match_mu++;
        }
    }

    auto findBestMatch = [](const std::vector<MatchInfo>& matches) -> MatchInfo {
        auto it = std::min_element(matches.begin(), matches.end(),
                [](const auto& a, const auto& b) { return a.dr < b.dr; });
        return (it != matches.end()) ? *it : MatchInfo{-1, 100};
    };

    auto bestmu = findBestMatch(Matchmu);

    if( ftree->num_match_mu > 0 ){

        ftree->match_mu_idx = bestmu.index;

        const auto& match_muon = (*muonHandle)[bestmu.index];

        ftree->match_mu_dR = bestmu.dr;
        ftree->match_mu_charge = match_muon.charge();
        ftree->match_mu_eta = match_muon.eta();
        ftree->match_mu_phi = match_muon.phi();
        ftree->match_mu_vx = match_muon.vx();
        ftree->match_mu_vy = match_muon.vy();
        ftree->match_mu_vz = match_muon.vz();
        ftree->match_mu_p = match_muon.p();
        ftree->match_mu_pt = match_muon.pt();
        ftree->match_mu_px = match_muon.px();
        ftree->match_mu_py = match_muon.py();
        ftree->match_mu_pz = match_muon.pz();
        ftree->match_mu_isHighPt = match_muon.isHighPtMuon(primaryvertex);
        ftree->match_mu_isLoose = match_muon.isLooseMuon();
        ftree->match_mu_isMedium = match_muon.isMediumMuon();
        ftree->match_mu_isSoft = match_muon.isSoftMuon(primaryvertex);
        ftree->match_mu_isTight = match_muon.isTightMuon(primaryvertex);
        ftree->match_mu_isPF = match_muon.isPFMuon();
        ftree->match_mu_isTracker = match_muon.isTrackerMuon();
        ftree->match_mu_isGlobal = match_muon.isGlobalMuon();
        ftree->match_mu_IsoR03_sumChargedHadronPt = match_muon.pfIsolationR03().sumChargedHadronPt;
        ftree->match_mu_IsoR03_sumChargedParticlePt = match_muon.pfIsolationR03().sumChargedParticlePt;
        ftree->match_mu_IsoR03_sumNeutralHadronEt = match_muon.pfIsolationR03().sumNeutralHadronEt;
        ftree->match_mu_IsoR03_sumPhotonEt = match_muon.pfIsolationR03().sumPhotonEt;
        ftree->match_mu_IsoR03_sumPUPt = match_muon.pfIsolationR03().sumPUPt;
        ftree->match_mu_PFIsoR03 = (ftree->match_mu_IsoR03_sumChargedHadronPt + std::max(0.0, ftree->match_mu_IsoR03_sumNeutralHadronEt + ftree->match_mu_IsoR03_sumPhotonEt - 0.5*ftree->match_mu_IsoR03_sumPUPt))/ftree->match_mu_pt;
        ftree->match_mu_IsoR04_sumChargedHadronPt = match_muon.pfIsolationR04().sumChargedHadronPt;
        ftree->match_mu_IsoR04_sumChargedParticlePt = match_muon.pfIsolationR04().sumChargedParticlePt;
        ftree->match_mu_IsoR04_sumNeutralHadronEt = match_muon.pfIsolationR04().sumNeutralHadronEt;
        ftree->match_mu_IsoR04_sumPhotonEt = match_muon.pfIsolationR04().sumPhotonEt;
        ftree->match_mu_IsoR04_sumPUPt = match_muon.pfIsolationR04().sumPUPt;
        ftree->match_mu_PFIsoR04 = (ftree->match_mu_IsoR04_sumChargedHadronPt + std::max(0.0, ftree->match_mu_IsoR04_sumNeutralHadronEt + ftree->match_mu_IsoR04_sumPhotonEt - 0.5*ftree->match_mu_IsoR04_sumPUPt))/ftree->match_mu_pt;
        ftree->match_mu_primvtx_dxy = match_muon.muonBestTrack()->dxy(primaryvertex.position());
        ftree->match_mu_primvtx_dxyerr = match_muon.muonBestTrack()->dxyError();
        ftree->match_mu_primvtx_dz = match_muon.muonBestTrack()->dz(primaryvertex.position());
        ftree->match_mu_primvtx_dzerr = match_muon.muonBestTrack()->dzError();

        std::pair<bool, Measurement1D> match_mu_primvtx_IPresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(match_muon.muonBestTrack()), primaryvertex);
        if( match_mu_primvtx_IPresult.first ){
            ftree->match_mu_primvtx_ip = match_mu_primvtx_IPresult.second.value();
            ftree->match_mu_primvtx_iperr = match_mu_primvtx_IPresult.second.error();
            ftree->match_mu_primvtx_ipchi2 = match_mu_primvtx_IPresult.second.significance();
        }

        ftree->Match_Fill_Vector();
    }

    idx_mu_vec.clear();

    ftree->mu_Reset();

    for(size_t i=0; i<muonHandle->size(); i++){
        const auto& muon = (*muonHandle)[i];
        if(muon.charge() != -1) continue; 
        if(muon.pt() < 0.5) continue;
        if(muon.p() < 1) continue;

        ftree->mu_charge = muon.charge();
        ftree->mu_eta = muon.eta();
        ftree->mu_phi = muon.phi();
        ftree->mu_vx = muon.vx();
        ftree->mu_vy = muon.vy();
        ftree->mu_vz = muon.vz();
        ftree->mu_p = muon.p();
        ftree->mu_pt = muon.pt();
        ftree->mu_px = muon.px();
        ftree->mu_py = muon.py();
        ftree->mu_pz = muon.pz();
        ftree->mu_isHighPt = muon.isHighPtMuon(primaryvertex);
        ftree->mu_isLoose = muon.isLooseMuon();
        ftree->mu_isMedium = muon.isMediumMuon();
        ftree->mu_isSoft = muon.isSoftMuon(primaryvertex);
        ftree->mu_isTight = muon.isTightMuon(primaryvertex);
        ftree->mu_isPF = muon.isPFMuon();
        ftree->mu_isTracker = muon.isTrackerMuon();
        ftree->mu_isGlobal = muon.isGlobalMuon();
        ftree->mu_IsoR03_sumChargedHadronPt = muon.pfIsolationR03().sumChargedHadronPt;
        ftree->mu_IsoR03_sumChargedParticlePt = muon.pfIsolationR03().sumChargedParticlePt;
        ftree->mu_IsoR03_sumNeutralHadronEt = muon.pfIsolationR03().sumNeutralHadronEt;
        ftree->mu_IsoR03_sumPhotonEt = muon.pfIsolationR03().sumPhotonEt;
        ftree->mu_IsoR03_sumPUPt = muon.pfIsolationR03().sumPUPt;
        ftree->mu_PFIsoR03 = (ftree->mu_IsoR03_sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR03_sumNeutralHadronEt + ftree->mu_IsoR03_sumPhotonEt - 0.5*ftree->mu_IsoR03_sumPUPt))/ftree->mu_pt;
        ftree->mu_IsoR04_sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
        ftree->mu_IsoR04_sumChargedParticlePt = muon.pfIsolationR04().sumChargedParticlePt;
        ftree->mu_IsoR04_sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
        ftree->mu_IsoR04_sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
        ftree->mu_IsoR04_sumPUPt = muon.pfIsolationR04().sumPUPt;
        ftree->mu_PFIsoR04 = (ftree->mu_IsoR04_sumChargedHadronPt + std::max(0.0, ftree->mu_IsoR04_sumNeutralHadronEt + ftree->mu_IsoR04_sumPhotonEt - 0.5*ftree->mu_IsoR04_sumPUPt))/ftree->mu_pt;
        ftree->mu_primvtx_dxy = muon.muonBestTrack()->dxy(primaryvertex.position());
        ftree->mu_primvtx_dxyerr = muon.muonBestTrack()->dxyError();
        ftree->mu_primvtx_dz = muon.muonBestTrack()->dz(primaryvertex.position());
        ftree->mu_primvtx_dzerr = muon.muonBestTrack()->dzError();

        std::pair<bool, Measurement1D> mu_primvtx_IPresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(muon.muonBestTrack()), primaryvertex);
        if( mu_primvtx_IPresult.first ){
            ftree->mu_primvtx_ip = mu_primvtx_IPresult.second.value();
            ftree->mu_primvtx_iperr = mu_primvtx_IPresult.second.error();
            ftree->mu_primvtx_ipchi2 = mu_primvtx_IPresult.second.significance();
        }

        if(i == static_cast<size_t>(bestmu.index)) ftree->mu_match = true;
        else ftree->mu_match = false;

        ftree->num_mu++;

        ftree->Fill_Vector();
    }

    ftree->tree->Fill();

    return;
}

bool MuonStudy::hasAncestor(const reco::GenParticle & gp, const int motherid, const int otherparticleid) const
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

void MuonStudy::endJob()
{
    edm::LogInfo("MuonStudy") << "Processed all events.";
}

void MuonStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonStudy);

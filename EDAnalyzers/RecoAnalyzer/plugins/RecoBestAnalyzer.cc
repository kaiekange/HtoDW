// -*- C++ -*-
//
// Package:    EDAnalyzer/GenParticleAnalyzer
// Class:      RecoBestAnalyzer
//
/**\class RecoBestAnalyzer RecoBestAnalyzer.cc EDAnalyzer/GenParticleAnalyzer/plugins/RecoBestAnalyzer.cc

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

#include "EDAnalyzers/RecoAnalyzer/interface/RecoBestTree.h"
#include "EDAnalyzers/RecoAnalyzer/interface/VertexReProducer.h"


class RecoBestAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
        explicit RecoBestAnalyzer(const edm::ParameterSet& iConfig);
        ~RecoBestAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFToken;
        edm::EDGetTokenT<reco::BeamSpot> beamspotToken;

        VertexReProducer *revertex;

        const edm::Service<TFileService> fs;
        RecoBestTree *ftree;

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

RecoBestAnalyzer::RecoBestAnalyzer(const edm::ParameterSet& iConfig) :
    packedPFToken(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates"))),
    beamspotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")))
{
    ftree = new RecoBestTree(fs->make<TTree>("Events", "Events"));
    ftree->CreateBranches();

    revertex = new VertexReProducer(iConfig);
}

RecoBestAnalyzer::~RecoBestAnalyzer() {}

void RecoBestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    ftree->Init();

    edm::Handle<pat::PackedCandidateCollection> packedPF;
    iEvent.getByToken(packedPFToken, packedPF);

    edm::Handle<reco::BeamSpot> pvbeamspot;
    iEvent.getByToken(beamspotToken, pvbeamspot);

    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

    KalmanVertexFitter KVFitter(true);
    AdaptiveVertexFitter AVFitter;

    ftree->BS_Reset();
    ftree->BS_type = pvbeamspot->type();
    ftree->BS_X0 = pvbeamspot->x0();
    ftree->BS_Y0 = pvbeamspot->y0();
    ftree->BS_Z0 = pvbeamspot->z0();
    ftree->BS_SigmaZ = pvbeamspot->sigmaZ();
    ftree->BS_dXdZ = pvbeamspot->dxdz();
    ftree->BS_dYdZ = pvbeamspot->dydz();
    ftree->BS_BWX = pvbeamspot->BeamWidthX();
    ftree->BS_BWY = pvbeamspot->BeamWidthY();
    ftree->BS_X0ERR = pvbeamspot->x0Error();
    ftree->BS_Y0ERR = pvbeamspot->y0Error();
    ftree->BS_Z0ERR = pvbeamspot->z0Error();
    ftree->BS_SigmaZ0ERR = pvbeamspot->sigmaZ0Error();
    ftree->BS_dXdZERR = pvbeamspot->dxdzError();
    ftree->BS_dYdZERR = pvbeamspot->dydzError();
    ftree->BS_BWXERR = pvbeamspot->BeamWidthXError();
    ftree->BS_BWYERR = pvbeamspot->BeamWidthYError();
    ftree->BS_EmitX = pvbeamspot->emittanceX();
    ftree->BS_EmitY = pvbeamspot->emittanceY();
    ftree->BS_BetaStar = pvbeamspot->betaStar();

    ftree->PV_Reset();
    reco::TrackCollection tracks;
    for(size_t i=0; i<packedPF->size(); i++){
        const auto& pf = (*packedPF)[i];
        if( ! pf.hasTrackDetails() ) continue;
        tracks.push_back(pf.pseudoTrack());
    }
    std::vector<TransientVertex> pvs_withBS = revertex->makeVertices(tracks, *pvbeamspot, iSetup, "WithBS");
    if( !(pvs_withBS.empty()) ){
        reco::Vertex pv_withBS = reco::Vertex(pvs_withBS.front());

        ftree->PV_withBS_IsValid = pv_withBS.isValid();
        ftree->PV_withBS_IsFake = pv_withBS.isFake();
        ftree->PV_withBS_CHI2 = pv_withBS.chi2();
        ftree->PV_withBS_NDOF = pv_withBS.ndof();
        ftree->PV_withBS_CHI2NDOF = pv_withBS.chi2()/pv_withBS.ndof();
        ftree->PV_withBS_X = pv_withBS.x();
        ftree->PV_withBS_Y = pv_withBS.y();
        ftree->PV_withBS_Z = pv_withBS.z();
        ftree->PV_withBS_XERR = pv_withBS.xError();
        ftree->PV_withBS_YERR = pv_withBS.yError();
        ftree->PV_withBS_ZERR = pv_withBS.zError();

        ftree->PV_withBS_Fill_Vector();
    }

    std::vector<TransientVertex> pvs_noBS = revertex->makeVertices(tracks, *pvbeamspot, iSetup, "noBS");
    if( !(pvs_noBS.empty()) ){
        reco::Vertex pv_noBS = reco::Vertex(pvs_noBS.front());

        ftree->PV_noBS_IsValid = pv_noBS.isValid();
        ftree->PV_noBS_IsFake = pv_noBS.isFake();
        ftree->PV_noBS_CHI2 = pv_noBS.chi2();
        ftree->PV_noBS_NDOF = pv_noBS.ndof();
        ftree->PV_noBS_CHI2NDOF = pv_noBS.chi2()/pv_noBS.ndof();
        ftree->PV_noBS_X = pv_noBS.x();
        ftree->PV_noBS_Y = pv_noBS.y();
        ftree->PV_noBS_Z = pv_noBS.z();
        ftree->PV_noBS_XERR = pv_noBS.xError();
        ftree->PV_noBS_YERR = pv_noBS.yError();
        ftree->PV_noBS_ZERR = pv_noBS.zError();

        ftree->PV_noBS_Fill_Vector();
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
            if( ftree->dR_Kp_Km > 0.1 ) continue;

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

                ftree->Fill_Vector();
            }
        }
    }

    int num_candidates = ftree->num_reco_Ds;

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

void RecoBestAnalyzer::beginJob() {}

void RecoBestAnalyzer::endJob() {}

void RecoBestAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
    /* edm::ParameterSetDescription desc; */
    /* desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"));; */
    /* descriptions.add("RecoBestAnalyzer", desc); */
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoBestAnalyzer);

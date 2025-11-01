#ifndef PVTREE_H 
#define PVTREE_H

#include <TTree.h>
#include <vector>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TLorentzVector.h"
#include "Rtypes.h"

#define null -777


class PVTree 
{
    public:

        PVTree(TTree *tree_);
        TTree *tree;

        void CreateBranches();
        void Init();

        void Gen_Reset();

        struct GenInfo {
            double eta, phi, vx, vy, vz, p, pt, px, py, pz;
            void clearAll() {
                eta = null;
                phi = null;
                vx  = null;
                vy  = null;
                vz  = null;
                p   = null;
                pt  = null;
                px  = null;
                py  = null;
                pz  = null;
            }
            void fillAll( const reco::GenParticle& );
        };

        struct BSInfo {
            int type;
            double x0, y0, z0, sigmaZ, dxdz, dydz;
            double BWX, BWY, x0err, y0err, z0err, sigmaZ0err;
            double dxdzerr, dydzerr, BWXerr, BWYerr;
            double emitX, emitY, betaStar;
            void clearAll() {
                type       = null;
                x0         = null;
                y0         = null;
                z0         = null;
                sigmaZ     = null;
                dxdz       = null;
                dydz       = null;
                BWX        = null;
                BWY        = null;
                x0err      = null;
                y0err      = null;
                z0err      = null;
                sigmaZ0err = null;
                dxdzerr    = null;
                dydzerr    = null;
                BWXerr     = null;
                BWYerr     = null;
                emitX      = null;
                emitY      = null;
                betaStar   = null;
            }
            void fillAll( const reco::BeamSpot& );
        };

        struct VtxInfo {
            bool isValid, isFake;
            double chi2, ndof, chi2ndof;
            double vx, vy, vz, vxerr, vyerr, vzerr;
            void clearAll() {
                isValid  = false;
                isFake   = true;
                chi2     = null;
                ndof     = null;
                chi2ndof = null;
                vx       = null;
                vy       = null;
                vz       = null;
                vxerr    = null;
                vyerr    = null;
                vzerr    = null;
            }
            void fillAll( const reco::Vertex& );
        };
        struct VtxInfoVec {
            std::vector<bool> isValid, isFake;
            std::vector<double> chi2, ndof, chi2ndof;
            std::vector<double> vx, vy, vz, vxerr, vyerr, vzerr;
            void clearAll() {
                isValid.clear();
                isFake.clear();
                chi2.clear();
                ndof.clear();
                chi2ndof.clear();
                vx.clear();
                vy.clear();
                vz.clear();
                vxerr.clear();
                vyerr.clear();
                vzerr.clear();
            }
            void fillAll( const VtxInfo& info ) {
                isValid.push_back(info.isValid);
                isFake.push_back(info.isFake);
                chi2.push_back(info.chi2);
                ndof.push_back(info.ndof);
                chi2ndof.push_back(info.chi2ndof);
                vx.push_back(info.vx);
                vy.push_back(info.vy);
                vz.push_back(info.vz);
                vxerr.push_back(info.vxerr);
                vyerr.push_back(info.vyerr);
                vzerr.push_back(info.vzerr);
            }
            VtxInfo findIdx( int idx ) const {
                VtxInfo info;
                info.isValid  = isValid[idx];
                info.isFake   = isFake[idx];
                info.chi2     = chi2[idx];
                info.ndof     = ndof[idx];
                info.chi2ndof = chi2ndof[idx];
                info.vx       = vx[idx];
                info.vy       = vy[idx];
                info.vz       = vz[idx];
                info.vxerr    = vxerr[idx];
                info.vyerr    = vyerr[idx];
                info.vzerr    = vzerr[idx];
                return info;
            }
        };

        struct PFInfo {
            bool isIsolatedChargedHadron;
            int charge;
            double eta, phi, vx, vy, vz;
            double p, pt, px, py, pz;
            void clearAll() {
                isIsolatedChargedHadron = false;
                charge                  = null;
                eta                     = null;
                phi                     = null;
                vx                      = null;
                vy                      = null;
                vz                      = null;
                p                       = null;
                pt                      = null;
                px                      = null;
                py                      = null;
                pz                      = null;
            }
            void fillAll( const pat::PackedCandidate& );
        };
        struct PFInfoVec {
            std::vector<bool> isIsolatedChargedHadron;
            std::vector<int> charge;
            std::vector<double> eta, phi, vx, vy, vz;
            std::vector<double> p, pt, px, py, pz;
            void clearAll() {
                isIsolatedChargedHadron.clear();
                charge.clear();
                eta.clear();
                phi.clear();
                vx.clear();
                vy.clear();
                vz.clear();
                p.clear();
                pt.clear();
                px.clear();
                py.clear();
                pz.clear();
            }
            void fillAll( const PFInfo& info ) {
                isIsolatedChargedHadron.push_back(info.isIsolatedChargedHadron);
                charge.push_back(info.charge);
                eta.push_back(info.eta);
                phi.push_back(info.phi);
                vx.push_back(info.vx);
                vy.push_back(info.vy);
                vz.push_back(info.vz);
                p.push_back(info.p);
                pt.push_back(info.pt);
                px.push_back(info.px);
                py.push_back(info.py);
                pz.push_back(info.pz);
            }
            PFInfo findIdx(int idx) const {
                PFInfo info;
                info.isIsolatedChargedHadron = isIsolatedChargedHadron[idx];
                info.charge                  = charge[idx];
                info.eta                     = eta[idx];
                info.phi                     = phi[idx];
                info.vx                      = vx[idx];
                info.vy                      = vy[idx];
                info.vz                      = vz[idx];
                info.p                       = p[idx];
                info.pt                      = pt[idx];
                info.px                      = px[idx];
                info.py                      = py[idx];
                info.pz                      = pz[idx];
                return info;
            }
        };

        struct TLVInfo {
            double eta, phi, p, pt, px, py, pz, invm;
            void clearAll() {
                eta  = null;
                phi  = null;
                p    = null;
                pt   = null;
                px   = null;
                py   = null;
                pz   = null;
                invm = null;
            }
            void fillAll( const TLorentzVector& );
        };
        struct TLVInfoVec {
            std::vector<double> eta, phi, p, pt, px, py, pz, invm;
            void clearAll() {
                eta.clear();
                phi.clear();
                p.clear();
                pt.clear();
                px.clear();
                py.clear();
                pz.clear();
                invm.clear();
            }
            void fillAll( const TLVInfo& info ) {
                eta.push_back(info.eta);
                phi.push_back(info.phi);
                p.push_back(info.p);
                pt.push_back(info.pt);
                px.push_back(info.px);
                py.push_back(info.py);
                pz.push_back(info.pz);
                invm.push_back(info.invm);
            }
            TLVInfo findIdx(int idx) const {
                TLVInfo info;
                info.eta  = eta[idx];
                info.phi  = phi[idx];
                info.p    = p[idx];
                info.pt   = pt[idx];
                info.px   = px[idx];
                info.py   = py[idx];
                info.pz   = pz[idx];
                info.invm = invm[idx];
                return info;
            }
        };

        struct FDInfo {
            double FDxy, FDz, FD;
            double FDxyerr, FDzerr, FDerr;
            double FDxychi2, FDzchi2, FDchi2;
            double dira, dira_angle;
            void clearAll() {
                FDxy       = null;
                FDz        = null;
                FD         = null;
                FDxyerr    = null;
                FDzerr     = null;
                FDerr      = null;
                FDxychi2   = null;
                FDzchi2    = null;
                FDchi2     = null;
                dira       = null;
                dira_angle = null;
            }
            void fillAll( const GlobalVector&, const GlobalError&, const GlobalVector& );
        };
        struct FDInfoVec {
            std::vector<double> FDxy, FDz, FD;
            std::vector<double> FDxyerr, FDzerr, FDerr;
            std::vector<double> FDxychi2, FDzchi2, FDchi2;
            std::vector<double> dira, dira_angle;
            void clearAll() {
                FDxy.clear();
                FDz.clear();
                FD.clear();
                FDxyerr.clear();
                FDzerr.clear();
                FDerr.clear();
                FDxychi2.clear();
                FDzchi2.clear();
                FDchi2.clear();
                dira.clear();
                dira_angle.clear();
            }
            void fillAll( const FDInfo& info ) {
                FDxy.push_back(info.FDxy);
                FDz.push_back(info.FDz);
                FD.push_back(info.FD);
                FDxyerr.push_back(info.FDxyerr);
                FDzerr.push_back(info.FDzerr);
                FDerr.push_back(info.FDerr);
                FDxychi2.push_back(info.FDxychi2);
                FDzchi2.push_back(info.FDzchi2);
                FDchi2.push_back(info.FDchi2);
                dira.push_back(info.dira);
                dira_angle.push_back(info.dira_angle);
            }
            FDInfo findIdx(int idx) const {
                FDInfo info;
                info.FDxy       = FDxy[idx];
                info.FDz        = FDz[idx];
                info.FD         = FD[idx];
                info.FDxyerr    = FDxyerr[idx];
                info.FDzerr     = FDzerr[idx];
                info.FDerr      = FDerr[idx];
                info.FDxychi2   = FDxychi2[idx];
                info.FDzchi2    = FDzchi2[idx];
                info.FDchi2     = FDchi2[idx];
                info.dira       = dira[idx];
                info.dira_angle = dira_angle[idx];
                return info;
            }
        };

        struct IPInfo {
            double ip, iperr, ipchi2;
            void clearAll(){
                ip = null;
                iperr = null;
                ipchi2 = null;
            }
            void fillAll( const reco::TransientTrack&, const reco::Vertex& );
        };
        struct IPInfoVec {
            std::vector<double> ip, iperr, ipchi2;
            void clearAll() {
                ip.clear();
                iperr.clear();
                ipchi2.clear();
            }
            void fillAll( const IPInfo& info ){
                ip.push_back(info.ip);
                iperr.push_back(info.iperr);
                ipchi2.push_back(info.ipchi2);
            }
            IPInfo findIdx(int idx) const {
                IPInfo info;
                info.ip     = ip[idx];
                info.iperr  = iperr[idx];
                info.ipchi2 = ipchi2[idx];
                return info;
            }
        };

        struct IsoInfo {
            double sumChargedHadronPt, sumNeutralHadronPt, sumPhotonPt, sumPUPt, PFIso;
            void clearAll() {
                sumChargedHadronPt = 0.0;
                sumNeutralHadronPt = 0.0;
                sumPhotonPt        = 0.0;
                sumPUPt            = 0.0;
                PFIso              = 0.0;
            }
        };
        struct IsoInfoVec {
            std::vector<double> sumChargedHadronPt, sumNeutralHadronPt, sumPhotonPt, sumPUPt, PFIso;
            void clearAll() {
                sumChargedHadronPt.clear();
                sumNeutralHadronPt.clear();
                sumPhotonPt.clear();
                sumPUPt.clear();
                PFIso.clear();
            }
            void fillAll( const IsoInfo& info ){
                sumChargedHadronPt.push_back(info.sumChargedHadronPt);
                sumNeutralHadronPt.push_back(info.sumNeutralHadronPt);
                sumPhotonPt.push_back(info.sumPhotonPt);
                sumPUPt.push_back(info.sumPUPt);
                PFIso.push_back(info.PFIso);
            }
            IsoInfo findIdx(int idx) const {
                IsoInfo info;
                info.sumChargedHadronPt = sumChargedHadronPt[idx];
                info.sumNeutralHadronPt = sumNeutralHadronPt[idx];
                info.sumPhotonPt        = sumPhotonPt[idx];
                info.sumPUPt            = sumPUPt[idx];
                info.PFIso              = PFIso[idx];
                return info;
            }
        };

        struct MuonInfo {
            int charge;
            double eta, phi, vx, vy, vz;
            double p, pt, px, py, pz;
            double dxy, dxyerr, dz, dzerr;
            bool isHighPt, isLoose, isMedium;
            bool isSoft, isTight, isPF;
            bool isTracker, isGlobal;
            void clearAll(){
                charge    = null;
                eta       = null;
                phi       = null;
                vx        = null;
                vy        = null;
                vz        = null;
                p         = null;
                pt        = null;
                px        = null;
                py        = null;
                pz        = null;
                dxy       = null;
                dxyerr    = null;
                dz        = null;
                dzerr     = null;
                isHighPt  = false;
                isLoose   = false;
                isMedium  = false;
                isSoft    = false;
                isTight   = false;
                isPF      = false;
                isTracker = false;
                isGlobal  = false;
            }
            void fillAll( const pat::Muon&, const reco::Vertex& );
        };
        struct MuonInfoVec {
            std::vector<int> charge;
            std::vector<double> eta, phi, vx, vy, vz;
            std::vector<double> p, pt, px, py, pz;
            std::vector<double> dxy, dxyerr, dz, dzerr;
            std::vector<bool> isHighPt, isLoose, isMedium;
            std::vector<bool> isSoft, isTight, isPF;
            std::vector<bool> isTracker, isGlobal;
            void clearAll() {
                charge.clear();
                eta.clear();
                phi.clear();
                vx.clear();
                vy.clear();
                vz.clear();
                p.clear();
                pt.clear();
                px.clear();
                py.clear();
                pz.clear();
                dxy.clear();
                dxyerr.clear();
                dz.clear();
                dzerr.clear();
                isHighPt.clear();
                isLoose.clear();
                isMedium.clear();
                isSoft.clear();
                isTight.clear();
                isPF.clear();
                isTracker.clear();
                isGlobal.clear();
            }
            void fillAll( const MuonInfo& info ){
                charge.push_back(info.charge);
                eta.push_back(info.eta);
                phi.push_back(info.phi);
                vx.push_back(info.vx);
                vy.push_back(info.vy);
                vz.push_back(info.vz);
                p.push_back(info.p);
                pt.push_back(info.pt);
                px.push_back(info.px);
                py.push_back(info.py);
                pz.push_back(info.pz);
                dxy.push_back(info.dxy);
                dxyerr.push_back(info.dxyerr);
                dz.push_back(info.dz);
                dzerr.push_back(info.dzerr);
                isHighPt.push_back(info.isHighPt);
                isLoose.push_back(info.isLoose);
                isMedium.push_back(info.isMedium);
                isSoft.push_back(info.isSoft);
                isTight.push_back(info.isTight);
                isPF.push_back(info.isPF);
                isTracker.push_back(info.isTracker);
                isGlobal.push_back(info.isGlobal);
            }
            MuonInfo findIdx(int idx) const {
                MuonInfo info;
                info.charge    = charge[idx];
                info.eta       = eta[idx];
                info.phi       = phi[idx];
                info.vx        = vx[idx];
                info.vy        = vy[idx];
                info.vz        = vz[idx];
                info.p         = p[idx];
                info.pt        = pt[idx];
                info.px        = px[idx];
                info.py        = py[idx];
                info.pz        = pz[idx];
                info.dxy       = dxy[idx];
                info.dxyerr    = dxyerr[idx];
                info.dz        = dz[idx];
                info.dzerr     = dzerr[idx];
                info.isHighPt  = isHighPt[idx];
                info.isLoose   = isLoose[idx];
                info.isMedium  = isMedium[idx];
                info.isSoft    = isSoft[idx];
                info.isTight   = isTight[idx];
                info.isPF      = isPF[idx];
                info.isTracker = isTracker[idx];
                info.isGlobal  = isGlobal[idx];
                return info;
            }
        };


        void Match_Reset();
        void FillMatchDs();
        void FillMatchmu();

        void FillDs();
        void FillBestDs(int idxmax);
        void FillBestmu();

        int num_Gen_Kp, num_Gen_Km, num_Gen_pi, num_Gen_phi, num_Gen_Ds, num_Gen_mu, num_Gen_nu, num_Gen_W, num_Gen_H;
        GenInfo Gen_H, Gen_Ds, Gen_W, Gen_phi, Gen_Kp, Gen_Km, Gen_pi, Gen_mu, Gen_nu;

        double Gen_dR_Kp_Km;
        double Gen_dR_Kp_phi;
        double Gen_dR_Km_phi;
        double Gen_dR_Kp_pi; 
        double Gen_dR_Km_pi; 
        double Gen_dR_pi_phi; 
        double Gen_dR_Kp_Ds; 
        double Gen_dR_Km_Ds; 
        double Gen_dR_phi_Ds; 
        double Gen_dR_pi_Ds; 
        double Gen_dR_Kp_mu; 
        double Gen_dR_Km_mu; 
        double Gen_dR_phi_mu; 
        double Gen_dR_pi_mu; 
        double Gen_dR_Ds_mu;

        double Gen_Ds_dx;
        double Gen_Ds_dy;
        double Gen_Ds_dz;
        double Gen_Ds_FDxy;
        double Gen_Ds_FD;

        BSInfo beamspotInfo;
        VtxInfo primvtxInfo;

        // Matched particles
        int num_match_Kp;
        int num_match_Km;
        int num_match_pi;
        int num_match_mu;
        int num_tight_match_Kp;
        int num_tight_match_Km;
        int num_tight_match_pi;
        int match_dR_Kp;
        int match_dR_Km;
        int match_dR_pi;
        int match_dR_mu;

        PFInfo match_Kp;          PFInfoVec match_Kp_vec;
        PFInfo match_Km;          PFInfoVec match_Km_vec;
        PFInfo match_pi;          PFInfoVec match_pi_vec;
        TLVInfo match_phi;        TLVInfoVec match_phi_vec;
        TLVInfo match_Ds;         TLVInfoVec match_Ds_vec;
        TLVInfo match_phiFit_Kp;  TLVInfoVec match_phiFit_Kp_vec;
        TLVInfo match_phiFit_Km;  TLVInfoVec match_phiFit_Km_vec;
        TLVInfo match_phiFit_pi;  TLVInfoVec match_phiFit_pi_vec;
        TLVInfo match_phiFit_phi; TLVInfoVec match_phiFit_phi_vec;
        TLVInfo match_phiFit_Ds;  TLVInfoVec match_phiFit_Ds_vec;
        TLVInfo match_DsFit_Kp;   TLVInfoVec match_DsFit_Kp_vec;
        TLVInfo match_DsFit_Km;   TLVInfoVec match_DsFit_Km_vec;
        TLVInfo match_DsFit_pi;   TLVInfoVec match_DsFit_pi_vec;
        TLVInfo match_DsFit_phi;  TLVInfoVec match_DsFit_phi_vec;
        TLVInfo match_DsFit_Ds;   TLVInfoVec match_DsFit_Ds_vec;

        double match_dR_Kp_Km;  std::vector<double> match_dR_Kp_Km_vec;
        double match_dR_Kp_phi; std::vector<double> match_dR_Kp_phi_vec;
        double match_dR_Km_phi; std::vector<double> match_dR_Km_phi_vec;
        double match_dR_Kp_pi;  std::vector<double> match_dR_Kp_pi_vec;
        double match_dR_Km_pi;  std::vector<double> match_dR_Km_pi_vec;
        double match_dR_pi_phi; std::vector<double> match_dR_pi_phi_vec;
        double match_dR_Kp_Ds;  std::vector<double> match_dR_Kp_Ds_vec;
        double match_dR_Km_Ds;  std::vector<double> match_dR_Km_Ds_vec;
        double match_dR_phi_Ds; std::vector<double> match_dR_phi_Ds_vec;
        double match_dR_pi_Ds;  std::vector<double> match_dR_pi_Ds_vec;

        double match_phiFit_dR_Kp_Km;  std::vector<double> match_phiFit_dR_Kp_Km_vec;
        double match_phiFit_dR_Kp_phi; std::vector<double> match_phiFit_dR_Kp_phi_vec;
        double match_phiFit_dR_Km_phi; std::vector<double> match_phiFit_dR_Km_phi_vec;
        double match_phiFit_dR_Kp_pi;  std::vector<double> match_phiFit_dR_Kp_pi_vec;
        double match_phiFit_dR_Km_pi;  std::vector<double> match_phiFit_dR_Km_pi_vec;
        double match_phiFit_dR_pi_phi; std::vector<double> match_phiFit_dR_pi_phi_vec;
        double match_phiFit_dR_Kp_Ds;  std::vector<double> match_phiFit_dR_Kp_Ds_vec;
        double match_phiFit_dR_Km_Ds;  std::vector<double> match_phiFit_dR_Km_Ds_vec;
        double match_phiFit_dR_phi_Ds; std::vector<double> match_phiFit_dR_phi_Ds_vec;
        double match_phiFit_dR_pi_Ds;  std::vector<double> match_phiFit_dR_pi_Ds_vec;

        double match_DsFit_dR_Kp_Km;  std::vector<double> match_DsFit_dR_Kp_Km_vec;
        double match_DsFit_dR_Kp_phi; std::vector<double> match_DsFit_dR_Kp_phi_vec;
        double match_DsFit_dR_Km_phi; std::vector<double> match_DsFit_dR_Km_phi_vec;
        double match_DsFit_dR_Kp_pi;  std::vector<double> match_DsFit_dR_Kp_pi_vec;
        double match_DsFit_dR_Km_pi;  std::vector<double> match_DsFit_dR_Km_pi_vec;
        double match_DsFit_dR_pi_phi; std::vector<double> match_DsFit_dR_pi_phi_vec;
        double match_DsFit_dR_Kp_Ds;  std::vector<double> match_DsFit_dR_Kp_Ds_vec;
        double match_DsFit_dR_Km_Ds;  std::vector<double> match_DsFit_dR_Km_Ds_vec;
        double match_DsFit_dR_phi_Ds; std::vector<double> match_DsFit_dR_phi_Ds_vec;
        double match_DsFit_dR_pi_Ds;  std::vector<double> match_DsFit_dR_pi_Ds_vec;

        double match_dxy_phi_Ds; std::vector<double> match_dxy_phi_Ds_vec;
        double match_dz_phi_Ds;  std::vector<double> match_dz_phi_Ds_vec;

        double match_DsFit_Mconstraint_Ds_invm; std::vector<double> match_DsFit_Mconstraint_Ds_invm_vec;

        VtxInfo match_phiFit;   VtxInfoVec match_phiFit_vec;
        VtxInfo match_DsFit;    VtxInfoVec match_DsFit_vec;
        VtxInfo match_PVnoDs;   VtxInfoVec match_PVnoDs_vec;
        VtxInfo match_PVwithDs; VtxInfoVec match_PVwithDs_vec;

        FDInfo match_Ds_primvtx_FD;   FDInfoVec match_Ds_primvtx_FD_vec;
        IPInfo match_Kp_primvtx_ip;   IPInfoVec match_Kp_primvtx_ip_vec;
        IPInfo match_Km_primvtx_ip;   IPInfoVec match_Km_primvtx_ip_vec;
        IPInfo match_pi_primvtx_ip;   IPInfoVec match_pi_primvtx_ip_vec;
        IPInfo match_phi_primvtx_ip;  IPInfoVec match_phi_primvtx_ip_vec;
        IPInfo match_Ds_primvtx_ip;   IPInfoVec match_Ds_primvtx_ip_vec;
        FDInfo match_Ds_PVnoDs_FD;    FDInfoVec match_Ds_PVnoDs_FD_vec;
        IPInfo match_Kp_PVnoDs_ip;    IPInfoVec match_Kp_PVnoDs_ip_vec;
        IPInfo match_Km_PVnoDs_ip;    IPInfoVec match_Km_PVnoDs_ip_vec;
        IPInfo match_pi_PVnoDs_ip;    IPInfoVec match_pi_PVnoDs_ip_vec;
        IPInfo match_phi_PVnoDs_ip;   IPInfoVec match_phi_PVnoDs_ip_vec;
        IPInfo match_Ds_PVnoDs_ip;    IPInfoVec match_Ds_PVnoDs_ip_vec;
        FDInfo match_Ds_PVwithDs_FD;  FDInfoVec match_Ds_PVwithDs_FD_vec;
        IPInfo match_Kp_PVwithDs_ip;  IPInfoVec match_Kp_PVwithDs_ip_vec;
        IPInfo match_Km_PVwithDs_ip;  IPInfoVec match_Km_PVwithDs_ip_vec;
        IPInfo match_pi_PVwithDs_ip;  IPInfoVec match_pi_PVwithDs_ip_vec;
        IPInfo match_phi_PVwithDs_ip; IPInfoVec match_phi_PVwithDs_ip_vec;
        IPInfo match_Ds_PVwithDs_ip;  IPInfoVec match_Ds_PVwithDs_ip_vec;

        IsoInfo match_Ds_IsoR03; IsoInfoVec match_Ds_IsoR03_vec;
        IsoInfo match_Ds_IsoR04; IsoInfoVec match_Ds_IsoR04_vec;

        MuonInfo match_mu;          MuonInfoVec match_mu_vec;
        IsoInfo match_mu_IsoR03;    IsoInfoVec match_mu_IsoR03_vec;
        IsoInfo match_mu_IsoR04;    IsoInfoVec match_mu_IsoR04_vec;
        IPInfo match_mu_primvtx_ip; IPInfoVec match_mu_primvtx_ip_vec;

        int num_reco_phi;
        int num_reco_Ds;

        PFInfo Kp;          PFInfoVec Kp_vec;
        PFInfo Km;          PFInfoVec Km_vec;
        PFInfo pi;          PFInfoVec pi_vec;
        TLVInfo phi;        TLVInfoVec phi_vec;
        TLVInfo Ds;         TLVInfoVec Ds_vec;
        TLVInfo phiFit_Kp;  TLVInfoVec phiFit_Kp_vec;
        TLVInfo phiFit_Km;  TLVInfoVec phiFit_Km_vec;
        TLVInfo phiFit_pi;  TLVInfoVec phiFit_pi_vec;
        TLVInfo phiFit_phi; TLVInfoVec phiFit_phi_vec;
        TLVInfo phiFit_Ds;  TLVInfoVec phiFit_Ds_vec;
        TLVInfo DsFit_Kp;   TLVInfoVec DsFit_Kp_vec;
        TLVInfo DsFit_Km;   TLVInfoVec DsFit_Km_vec;
        TLVInfo DsFit_pi;   TLVInfoVec DsFit_pi_vec;
        TLVInfo DsFit_phi;  TLVInfoVec DsFit_phi_vec;
        TLVInfo DsFit_Ds;   TLVInfoVec DsFit_Ds_vec;

        double dR_Kp_Km;  std::vector<double> dR_Kp_Km_vec;
        double dR_Kp_phi; std::vector<double> dR_Kp_phi_vec;
        double dR_Km_phi; std::vector<double> dR_Km_phi_vec;
        double dR_Kp_pi;  std::vector<double> dR_Kp_pi_vec;
        double dR_Km_pi;  std::vector<double> dR_Km_pi_vec;
        double dR_pi_phi; std::vector<double> dR_pi_phi_vec;
        double dR_Kp_Ds;  std::vector<double> dR_Kp_Ds_vec;
        double dR_Km_Ds;  std::vector<double> dR_Km_Ds_vec;
        double dR_phi_Ds; std::vector<double> dR_phi_Ds_vec;
        double dR_pi_Ds;  std::vector<double> dR_pi_Ds_vec;

        double phiFit_dR_Kp_Km;  std::vector<double> phiFit_dR_Kp_Km_vec;
        double phiFit_dR_Kp_phi; std::vector<double> phiFit_dR_Kp_phi_vec;
        double phiFit_dR_Km_phi; std::vector<double> phiFit_dR_Km_phi_vec;
        double phiFit_dR_Kp_pi;  std::vector<double> phiFit_dR_Kp_pi_vec;
        double phiFit_dR_Km_pi;  std::vector<double> phiFit_dR_Km_pi_vec;
        double phiFit_dR_pi_phi; std::vector<double> phiFit_dR_pi_phi_vec;
        double phiFit_dR_Kp_Ds;  std::vector<double> phiFit_dR_Kp_Ds_vec;
        double phiFit_dR_Km_Ds;  std::vector<double> phiFit_dR_Km_Ds_vec;
        double phiFit_dR_phi_Ds; std::vector<double> phiFit_dR_phi_Ds_vec;
        double phiFit_dR_pi_Ds;  std::vector<double> phiFit_dR_pi_Ds_vec;

        double DsFit_dR_Kp_Km;  std::vector<double> DsFit_dR_Kp_Km_vec;
        double DsFit_dR_Kp_phi; std::vector<double> DsFit_dR_Kp_phi_vec;
        double DsFit_dR_Km_phi; std::vector<double> DsFit_dR_Km_phi_vec;
        double DsFit_dR_Kp_pi;  std::vector<double> DsFit_dR_Kp_pi_vec;
        double DsFit_dR_Km_pi;  std::vector<double> DsFit_dR_Km_pi_vec;
        double DsFit_dR_pi_phi; std::vector<double> DsFit_dR_pi_phi_vec;
        double DsFit_dR_Kp_Ds;  std::vector<double> DsFit_dR_Kp_Ds_vec;
        double DsFit_dR_Km_Ds;  std::vector<double> DsFit_dR_Km_Ds_vec;
        double DsFit_dR_phi_Ds; std::vector<double> DsFit_dR_phi_Ds_vec;
        double DsFit_dR_pi_Ds;  std::vector<double> DsFit_dR_pi_Ds_vec;

        double dxy_phi_Ds; std::vector<double> dxy_phi_Ds_vec;
        double dz_phi_Ds;  std::vector<double> dz_phi_Ds_vec;

        double DsFit_Mconstraint_Ds_invm; std::vector<double> DsFit_Mconstraint_Ds_invm_vec;

        VtxInfo phiFit;   VtxInfoVec phiFit_vec;
        VtxInfo DsFit;    VtxInfoVec DsFit_vec;
        VtxInfo PVnoDs;   VtxInfoVec PVnoDs_vec;
        VtxInfo PVwithDs; VtxInfoVec PVwithDs_vec;

        FDInfo Ds_primvtx_FD;   FDInfoVec Ds_primvtx_FD_vec;
        IPInfo Kp_primvtx_ip;   IPInfoVec Kp_primvtx_ip_vec;
        IPInfo Km_primvtx_ip;   IPInfoVec Km_primvtx_ip_vec;
        IPInfo pi_primvtx_ip;   IPInfoVec pi_primvtx_ip_vec;
        IPInfo phi_primvtx_ip;  IPInfoVec phi_primvtx_ip_vec;
        IPInfo Ds_primvtx_ip;   IPInfoVec Ds_primvtx_ip_vec;
        FDInfo Ds_PVnoDs_FD;    FDInfoVec Ds_PVnoDs_FD_vec;
        IPInfo Kp_PVnoDs_ip;    IPInfoVec Kp_PVnoDs_ip_vec;
        IPInfo Km_PVnoDs_ip;    IPInfoVec Km_PVnoDs_ip_vec;
        IPInfo pi_PVnoDs_ip;    IPInfoVec pi_PVnoDs_ip_vec;
        IPInfo phi_PVnoDs_ip;   IPInfoVec phi_PVnoDs_ip_vec;
        IPInfo Ds_PVnoDs_ip;    IPInfoVec Ds_PVnoDs_ip_vec;
        FDInfo Ds_PVwithDs_FD;  FDInfoVec Ds_PVwithDs_FD_vec;
        IPInfo Kp_PVwithDs_ip;  IPInfoVec Kp_PVwithDs_ip_vec;
        IPInfo Km_PVwithDs_ip;  IPInfoVec Km_PVwithDs_ip_vec;
        IPInfo pi_PVwithDs_ip;  IPInfoVec pi_PVwithDs_ip_vec;
        IPInfo phi_PVwithDs_ip; IPInfoVec phi_PVwithDs_ip_vec;
        IPInfo Ds_PVwithDs_ip;  IPInfoVec Ds_PVwithDs_ip_vec;

        IsoInfo Ds_IsoR03; IsoInfoVec Ds_IsoR03_vec;
        IsoInfo Ds_IsoR04; IsoInfoVec Ds_IsoR04_vec;

        bool Kp_match;        std::vector<bool> Kp_match_vec;
        bool Km_match;        std::vector<bool> Km_match_vec;
        bool pi_match;        std::vector<bool> pi_match_vec;
        bool match_entry;     std::vector<bool> match_entry_vec;
        bool non_match_entry; std::vector<bool> non_match_entry_vec;

        PFInfoVec best_Kp_vec;
        PFInfoVec best_Km_vec;
        PFInfoVec best_pi_vec;
        TLVInfoVec best_phi_vec;
        TLVInfoVec best_Ds_vec;
        TLVInfoVec best_phiFit_Kp_vec;
        TLVInfoVec best_phiFit_Km_vec;
        TLVInfoVec best_phiFit_pi_vec;
        TLVInfoVec best_phiFit_phi_vec;
        TLVInfoVec best_phiFit_Ds_vec;
        TLVInfoVec best_DsFit_Kp_vec;
        TLVInfoVec best_DsFit_Km_vec;
        TLVInfoVec best_DsFit_pi_vec;
        TLVInfoVec best_DsFit_phi_vec;
        TLVInfoVec best_DsFit_Ds_vec;

        std::vector<double> best_dR_Kp_Km_vec;
        std::vector<double> best_dR_Kp_phi_vec;
        std::vector<double> best_dR_Km_phi_vec;
        std::vector<double> best_dR_Kp_pi_vec;
        std::vector<double> best_dR_Km_pi_vec;
        std::vector<double> best_dR_pi_phi_vec;
        std::vector<double> best_dR_Kp_Ds_vec;
        std::vector<double> best_dR_Km_Ds_vec;
        std::vector<double> best_dR_phi_Ds_vec;
        std::vector<double> best_dR_pi_Ds_vec;

        std::vector<double> best_phiFit_dR_Kp_Km_vec;
        std::vector<double> best_phiFit_dR_Kp_phi_vec;
        std::vector<double> best_phiFit_dR_Km_phi_vec;
        std::vector<double> best_phiFit_dR_Kp_pi_vec;
        std::vector<double> best_phiFit_dR_Km_pi_vec;
        std::vector<double> best_phiFit_dR_pi_phi_vec;
        std::vector<double> best_phiFit_dR_Kp_Ds_vec;
        std::vector<double> best_phiFit_dR_Km_Ds_vec;
        std::vector<double> best_phiFit_dR_phi_Ds_vec;
        std::vector<double> best_phiFit_dR_pi_Ds_vec;

        std::vector<double> best_DsFit_dR_Kp_Km_vec;
        std::vector<double> best_DsFit_dR_Kp_phi_vec;
        std::vector<double> best_DsFit_dR_Km_phi_vec;
        std::vector<double> best_DsFit_dR_Kp_pi_vec;
        std::vector<double> best_DsFit_dR_Km_pi_vec;
        std::vector<double> best_DsFit_dR_pi_phi_vec;
        std::vector<double> best_DsFit_dR_Kp_Ds_vec;
        std::vector<double> best_DsFit_dR_Km_Ds_vec;
        std::vector<double> best_DsFit_dR_phi_Ds_vec;
        std::vector<double> best_DsFit_dR_pi_Ds_vec;

        std::vector<double> best_dxy_phi_Ds_vec;
        std::vector<double> best_dz_phi_Ds_vec;

        std::vector<double> best_DsFit_Mconstraint_Ds_invm_vec;

        VtxInfoVec best_phiFit_vec;
        VtxInfoVec best_DsFit_vec;
        VtxInfoVec best_PVnoDs_vec;
        VtxInfoVec best_PVwithDs_vec;

        FDInfoVec best_Ds_primvtx_FD_vec;
        IPInfoVec best_Kp_primvtx_ip_vec;
        IPInfoVec best_Km_primvtx_ip_vec;
        IPInfoVec best_pi_primvtx_ip_vec;
        IPInfoVec best_phi_primvtx_ip_vec;
        IPInfoVec best_Ds_primvtx_ip_vec;
        FDInfoVec best_Ds_PVnoDs_FD_vec;
        IPInfoVec best_Kp_PVnoDs_ip_vec;
        IPInfoVec best_Km_PVnoDs_ip_vec;
        IPInfoVec best_pi_PVnoDs_ip_vec;
        IPInfoVec best_phi_PVnoDs_ip_vec;
        IPInfoVec best_Ds_PVnoDs_ip_vec;
        FDInfoVec best_Ds_PVwithDs_FD_vec;
        IPInfoVec best_Kp_PVwithDs_ip_vec;
        IPInfoVec best_Km_PVwithDs_ip_vec;
        IPInfoVec best_pi_PVwithDs_ip_vec;
        IPInfoVec best_phi_PVwithDs_ip_vec;
        IPInfoVec best_Ds_PVwithDs_ip_vec;

        IsoInfoVec best_Ds_IsoR03_vec;
        IsoInfoVec best_Ds_IsoR04_vec;
        std::vector<bool> best_match_entry_vec;

        MuonInfo mu;          MuonInfoVec mu_vec;
        IsoInfo mu_IsoR03;    IsoInfoVec mu_IsoR03_vec;
        IsoInfo mu_IsoR04;    IsoInfoVec mu_IsoR04_vec;
        IPInfo mu_primvtx_ip; IPInfoVec mu_primvtx_ip_vec;

        bool best_mu_match;   std::vector<bool> best_mu_match_vec; 

};

#endif

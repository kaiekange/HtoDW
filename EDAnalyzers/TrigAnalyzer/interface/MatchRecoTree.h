#ifndef MATCHRECOTREE_H
#define MATCHRECOTREE_H

#include <TTree.h>

#define null -777

class MatchRecoTree
{
public:
    MatchRecoTree(TTree *tree_);
    TTree *tree;

    void CreateBranches();
    void Init();

    Bool_t primvtx_isValid;
    Bool_t primvtx_isFake;
    Float_t primvtx_chi2;
    Float_t primvtx_ndof;
    Float_t primvtx_chi2ndof;
    Float_t primvtx_vx;
    Float_t primvtx_vy;
    Float_t primvtx_vz;
    Float_t primvtx_vxerr;
    Float_t primvtx_vyerr;
    Float_t primvtx_vzerr;

    Int_t Kp_charge;
    Float_t Kp_eta;
    Float_t Kp_phi;
    Float_t Kp_px;
    Float_t Kp_py;
    Float_t Kp_pz;
    Float_t Kp_pt;
    Float_t Kp_p;

    Int_t Km_charge;
    Float_t Km_eta;
    Float_t Km_phi;
    Float_t Km_px;
    Float_t Km_py;
    Float_t Km_pz;
    Float_t Km_pt;
    Float_t Km_p;

    Int_t pi_charge;
    Float_t pi_eta;
    Float_t pi_phi;
    Float_t pi_px;
    Float_t pi_py;
    Float_t pi_pz;
    Float_t pi_pt;
    Float_t pi_p;

    Int_t phi_charge;
    Float_t phi_eta;
    Float_t phi_phi;
    Float_t phi_px;
    Float_t phi_py;
    Float_t phi_pz;
    Float_t phi_pt;
    Float_t phi_p;
    Float_t phi_m;

    Int_t Ds_charge;
    Float_t Ds_eta;
    Float_t Ds_phi;
    Float_t Ds_px;
    Float_t Ds_py;
    Float_t Ds_pz;
    Float_t Ds_pt;
    Float_t Ds_p;
    Float_t Ds_m;

    Bool_t phiFit_isValid;
    Bool_t phiFit_isFake;
    Float_t phiFit_chi2;
    Float_t phiFit_ndof;
    Float_t phiFit_chi2ndof;
    Float_t phiFit_vx;
    Float_t phiFit_vy;
    Float_t phiFit_vz;
    Float_t phiFit_vxerr;
    Float_t phiFit_vyerr;
    Float_t phiFit_vzerr;

    Bool_t DsFit_isValid;
    Bool_t DsFit_isFake;
    Float_t DsFit_chi2;
    Float_t DsFit_ndof;
    Float_t DsFit_chi2ndof;
    Float_t DsFit_vx;
    Float_t DsFit_vy;
    Float_t DsFit_vz;
    Float_t DsFit_vxerr;
    Float_t DsFit_vyerr;
    Float_t DsFit_vzerr;

    Float_t orig_dR_Kp_Km;
    Float_t orig_dR_Kp_pi;
    Float_t orig_dR_Km_pi;

    Float_t phiFit_dR_Kp_Km;
    Float_t phiFit_dR_Kp_pi;
    Float_t phiFit_dR_Kp_phi;
    Float_t phiFit_dR_Km_pi;
    Float_t phiFit_dR_Km_phi;
    Float_t phiFit_dR_pi_phi;

    Float_t dR_Kp_Km;
    Float_t dR_Kp_pi;
    Float_t dR_Kp_phi;
    Float_t dR_Kp_Ds;
    Float_t dR_Km_pi;
    Float_t dR_Km_phi;
    Float_t dR_Km_Ds;
    Float_t dR_pi_phi;
    Float_t dR_pi_Ds;
    Float_t dR_phi_Ds;
    Float_t dxy_phi_Ds;
    Float_t dz_phi_Ds;

    Float_t Kp_ip;
    Float_t Kp_iperr;
    Float_t Kp_ipchi2;
    Float_t Km_ip;
    Float_t Km_iperr;
    Float_t Km_ipchi2;
    Float_t pi_ip;
    Float_t pi_iperr;
    Float_t pi_ipchi2;
    Float_t Ds_FDxy;
    Float_t Ds_FDxyerr;
    Float_t Ds_FDxysig;
    Float_t Ds_FD;
    Float_t Ds_FDerr;
    Float_t Ds_FDsig;
    Float_t dira;
    Float_t dira_angle;

    Float_t Ds_IsoR03_sumChargedHadronPt;
    Float_t Ds_IsoR03_sumNeutralHadronEt;
    Float_t Ds_IsoR03_sumPhotonEt;
    Float_t Ds_IsoR03_sumPUPt;
    Float_t Ds_IsoR03_PFIso;
    Float_t Ds_IsoR03_RelIso;
    Float_t Ds_IsoR04_sumChargedHadronPt;
    Float_t Ds_IsoR04_sumNeutralHadronEt;
    Float_t Ds_IsoR04_sumPhotonEt;
    Float_t Ds_IsoR04_sumPUPt;
    Float_t Ds_IsoR04_PFIso;
    Float_t Ds_IsoR04_RelIso;

    Int_t mu_charge;
    Float_t mu_eta;
    Float_t mu_phi;
    Float_t mu_px;
    Float_t mu_py;
    Float_t mu_pz;
    Float_t mu_pt;
    Float_t mu_p;
    Bool_t mu_isHighPt;
    Bool_t mu_isLoose;
    Bool_t mu_isMedium;
    Bool_t mu_isSoft;
    Bool_t mu_isTight;
    Bool_t mu_isPF;
    Bool_t mu_isTracker;
    Bool_t mu_isGlobal;

    Float_t dR_Kp_mu;
    Float_t dR_Km_mu;
    Float_t dR_pi_mu;
    Float_t dR_phi_mu;
    Float_t dR_Ds_mu;

    Float_t mu_ip;
    Float_t mu_iperr;
    Float_t mu_ipchi2;

    Float_t mu_IsoR03_sumChargedHadronPt;
    Float_t mu_IsoR03_sumNeutralHadronEt;
    Float_t mu_IsoR03_sumPhotonEt;
    Float_t mu_IsoR03_sumPUPt;
    Float_t mu_IsoR03_PFIso;
    Float_t mu_IsoR03_RelIso;
    Float_t mu_IsoR04_sumChargedHadronPt;
    Float_t mu_IsoR04_sumNeutralHadronEt;
    Float_t mu_IsoR04_sumPhotonEt;
    Float_t mu_IsoR04_sumPUPt;
    Float_t mu_IsoR04_PFIso;
    Float_t mu_IsoR04_RelIso;
};
#endif
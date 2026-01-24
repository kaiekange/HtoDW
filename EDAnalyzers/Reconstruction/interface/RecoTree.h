#ifndef RECOTREE_H
#define RECOTREE_H

#include <TTree.h>

#define null -777

class RecoTree
{
public:
    RecoTree(TTree *tree_);
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

    std::vector<Int_t> K1_charge;
    std::vector<Float_t> K1_eta;
    std::vector<Float_t> K1_phi;
    std::vector<Float_t> K1_px;
    std::vector<Float_t> K1_py;
    std::vector<Float_t> K1_pz;
    std::vector<Float_t> K1_pt;
    std::vector<Float_t> K1_p;

    std::vector<Int_t> K2_charge;
    std::vector<Float_t> K2_eta;
    std::vector<Float_t> K2_phi;
    std::vector<Float_t> K2_px;
    std::vector<Float_t> K2_py;
    std::vector<Float_t> K2_pz;
    std::vector<Float_t> K2_pt;
    std::vector<Float_t> K2_p;

    std::vector<Int_t> pi_charge;
    std::vector<Float_t> pi_eta;
    std::vector<Float_t> pi_phi;
    std::vector<Float_t> pi_px;
    std::vector<Float_t> pi_py;
    std::vector<Float_t> pi_pz;
    std::vector<Float_t> pi_pt;
    std::vector<Float_t> pi_p;

    std::vector<Int_t> phi_charge;
    std::vector<Float_t> phi_eta;
    std::vector<Float_t> phi_phi;
    std::vector<Float_t> phi_px;
    std::vector<Float_t> phi_py;
    std::vector<Float_t> phi_pz;
    std::vector<Float_t> phi_pt;
    std::vector<Float_t> phi_p;
    std::vector<Float_t> phi_m;

    std::vector<Int_t> Ds_charge;
    std::vector<Float_t> Ds_eta;
    std::vector<Float_t> Ds_phi;
    std::vector<Float_t> Ds_px;
    std::vector<Float_t> Ds_py;
    std::vector<Float_t> Ds_pz;
    std::vector<Float_t> Ds_pt;
    std::vector<Float_t> Ds_p;
    std::vector<Float_t> Ds_m;
    std::vector<Float_t> Ds_m_phiM;

    std::vector<Int_t> mu_charge;
    std::vector<Float_t> mu_eta;
    std::vector<Float_t> mu_phi;
    std::vector<Float_t> mu_px;
    std::vector<Float_t> mu_py;
    std::vector<Float_t> mu_pz;
    std::vector<Float_t> mu_pt;
    std::vector<Float_t> mu_p;
    std::vector<Bool_t> mu_isHighPt;
    std::vector<Bool_t> mu_isLoose;
    std::vector<Bool_t> mu_isMedium;
    std::vector<Bool_t> mu_isSoft;
    std::vector<Bool_t> mu_isTight;
    std::vector<Bool_t> mu_isPF;
    std::vector<Bool_t> mu_isTracker;
    std::vector<Bool_t> mu_isGlobal;

    std::vector<Float_t> dR_K1_K2;
    std::vector<Float_t> dR_K1_pi;
    std::vector<Float_t> dR_K1_phi;
    std::vector<Float_t> dR_K1_Ds;
    std::vector<Float_t> dR_K1_mu;
    std::vector<Float_t> dR_K2_pi;
    std::vector<Float_t> dR_K2_phi;
    std::vector<Float_t> dR_K2_Ds;
    std::vector<Float_t> dR_K2_mu;
    std::vector<Float_t> dR_pi_phi;
    std::vector<Float_t> dR_pi_Ds;
    std::vector<Float_t> dR_pi_mu;
    std::vector<Float_t> dR_phi_Ds;
    std::vector<Float_t> dR_phi_mu;
    std::vector<Float_t> dR_Ds_mu;

    std::vector<Bool_t> phiFit_isValid;
    std::vector<Bool_t> phiFit_isFake;
    std::vector<Float_t> phiFit_chi2;
    std::vector<Float_t> phiFit_ndof;
    std::vector<Float_t> phiFit_chi2ndof;
    std::vector<Float_t> phiFit_vx;
    std::vector<Float_t> phiFit_vy;
    std::vector<Float_t> phiFit_vz;
    std::vector<Float_t> phiFit_vxerr;
    std::vector<Float_t> phiFit_vyerr;
    std::vector<Float_t> phiFit_vzerr;

    std::vector<Bool_t> DsFit_isValid;
    std::vector<Bool_t> DsFit_isFake;
    std::vector<Float_t> DsFit_chi2;
    std::vector<Float_t> DsFit_ndof;
    std::vector<Float_t> DsFit_chi2ndof;
    std::vector<Float_t> DsFit_vx;
    std::vector<Float_t> DsFit_vy;
    std::vector<Float_t> DsFit_vz;
    std::vector<Float_t> DsFit_vxerr;
    std::vector<Float_t> DsFit_vyerr;
    std::vector<Float_t> DsFit_vzerr;

    std::vector<Float_t> dxy_phi_Ds;
    std::vector<Float_t> dz_phi_Ds;

    std::vector<Float_t> K1_ip;
    std::vector<Float_t> K1_iperr;
    std::vector<Float_t> K1_ipchi2;
    std::vector<Float_t> K2_ip;
    std::vector<Float_t> K2_iperr;
    std::vector<Float_t> K2_ipchi2;
    std::vector<Float_t> pi_ip;
    std::vector<Float_t> pi_iperr;
    std::vector<Float_t> pi_ipchi2;
    std::vector<Float_t> mu_ip;
    std::vector<Float_t> mu_iperr;
    std::vector<Float_t> mu_ipchi2;

    std::vector<Float_t> dira;
    std::vector<Float_t> dira_angle;

    std::vector<Float_t> Ds_FDxy;
    std::vector<Float_t> Ds_FDxyerr;
    std::vector<Float_t> Ds_FDxysig;

    std::vector<Float_t> Ds_FD;
    std::vector<Float_t> Ds_FDerr;
    std::vector<Float_t> Ds_FDsig;

    std::vector<Float_t> Ds_IsoR03_sumChargedHadronPt;
    std::vector<Float_t> Ds_IsoR03_sumNeutralHadronEt;
    std::vector<Float_t> Ds_IsoR03_sumPhotonEt;
    std::vector<Float_t> Ds_IsoR03_sumPUPt;
    std::vector<Float_t> Ds_IsoR03_PFIso;
    std::vector<Float_t> Ds_IsoR03_RelIso;

    std::vector<Float_t> Ds_IsoR04_sumChargedHadronPt;
    std::vector<Float_t> Ds_IsoR04_sumNeutralHadronEt;
    std::vector<Float_t> Ds_IsoR04_sumPhotonEt;
    std::vector<Float_t> Ds_IsoR04_sumPUPt;
    std::vector<Float_t> Ds_IsoR04_PFIso;
    std::vector<Float_t> Ds_IsoR04_RelIso;

    std::vector<Float_t> mu_IsoR03_sumChargedHadronPt;
    std::vector<Float_t> mu_IsoR03_sumNeutralHadronEt;
    std::vector<Float_t> mu_IsoR03_sumPhotonEt;
    std::vector<Float_t> mu_IsoR03_sumPUPt;
    std::vector<Float_t> mu_IsoR03_PFIso;
    std::vector<Float_t> mu_IsoR03_RelIso;

    std::vector<Float_t> mu_IsoR04_sumChargedHadronPt;
    std::vector<Float_t> mu_IsoR04_sumNeutralHadronEt;
    std::vector<Float_t> mu_IsoR04_sumPhotonEt;
    std::vector<Float_t> mu_IsoR04_sumPUPt;
    std::vector<Float_t> mu_IsoR04_PFIso;
    std::vector<Float_t> mu_IsoR04_RelIso;

    std::vector<Bool_t> Kp_match;
    std::vector<Bool_t> Km_match;
    std::vector<Bool_t> pi_match;
    std::vector<Bool_t> mu_match;
};
#endif

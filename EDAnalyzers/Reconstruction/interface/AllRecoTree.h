#ifndef ALLAllRecoTree_H
#define ALLAllRecoTree_H

#include <TTree.h>

#define null -777

class AllRecoTree
{
public:
    AllRecoTree(TTree *tree_);
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

    std::vector<Int_t> Kp_charge;
    std::vector<Float_t> Kp_eta;
    std::vector<Float_t> Kp_phi;
    std::vector<Float_t> Kp_px;
    std::vector<Float_t> Kp_py;
    std::vector<Float_t> Kp_pz;
    std::vector<Float_t> Kp_pt;
    std::vector<Float_t> Kp_p;

    std::vector<Int_t> Km_charge;
    std::vector<Float_t> Km_eta;
    std::vector<Float_t> Km_phi;
    std::vector<Float_t> Km_px;
    std::vector<Float_t> Km_py;
    std::vector<Float_t> Km_pz;
    std::vector<Float_t> Km_pt;
    std::vector<Float_t> Km_p;

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

    std::vector<Float_t> orig_dR_Kp_Km;
    std::vector<Float_t> orig_dR_Kp_pi;
    std::vector<Float_t> orig_dR_Km_pi;

    std::vector<Float_t> phiFit_dR_Kp_Km;
    std::vector<Float_t> phiFit_dR_Kp_pi;
    std::vector<Float_t> phiFit_dR_Kp_phi;
    std::vector<Float_t> phiFit_dR_Km_pi;
    std::vector<Float_t> phiFit_dR_Km_phi;
    std::vector<Float_t> phiFit_dR_pi_phi;
    
    std::vector<Float_t> dR_Kp_Km;
    std::vector<Float_t> dR_Kp_pi;
    std::vector<Float_t> dR_Kp_phi;
    std::vector<Float_t> dR_Kp_Ds;
    std::vector<Float_t> dR_Km_pi;
    std::vector<Float_t> dR_Km_phi;
    std::vector<Float_t> dR_Km_Ds;
    std::vector<Float_t> dR_pi_phi;
    std::vector<Float_t> dR_pi_Ds;
    std::vector<Float_t> dR_phi_Ds;

    std::vector<Float_t> dxy_phi_Ds;
    std::vector<Float_t> dz_phi_Ds;

    std::vector<Float_t> Kp_ip;
    std::vector<Float_t> Kp_iperr;
    std::vector<Float_t> Kp_ipchi2;
    std::vector<Float_t> Km_ip;
    std::vector<Float_t> Km_iperr;
    std::vector<Float_t> Km_ipchi2;
    std::vector<Float_t> pi_ip;
    std::vector<Float_t> pi_iperr;
    std::vector<Float_t> pi_ipchi2;

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

    std::vector<Bool_t> Kp_match;
    std::vector<Bool_t> Km_match;
    std::vector<Bool_t> pi_match;

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

    std::vector<Float_t> mu_ip;
    std::vector<Float_t> mu_iperr;
    std::vector<Float_t> mu_ipchi2;

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

    std::vector<Bool_t> mu_match;
};
#endif

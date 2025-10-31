#include "EDAnalyzers/RecoAnalyzer/interface/PVTree.h"
#include <iostream>

PVTree::PVTree(TTree *tree_)
{
    tree = tree_;
}
void PVTree::Init()
{
    // Gen particles
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    num_Gen_phi = 0;
    num_Gen_Ds = 0;
    num_Gen_mu = 0;
    num_Gen_nu = 0;
    num_Gen_W = 0;
    num_Gen_H = 0;

    // Matched particles
    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;
    num_match_mu = 0;
    num_tight_match_Kp = 0;
    num_tight_match_Km = 0;
    num_tight_match_pi = 0;

    match_dR_Kp = null;
    match_dR_Km = null;
    match_dR_pi = null;
    match_dR_mu = null;

    // Original info
    match_Kp_isIsolatedChargedHadron_vec.clear();
    match_Kp_charge_vec.clear();
    match_Kp_eta_vec.clear();
    match_Kp_phi_vec.clear();
    match_Kp_vx_vec.clear();
    match_Kp_vy_vec.clear();
    match_Kp_vz_vec.clear();
    match_Kp_p_vec.clear();
    match_Kp_pt_vec.clear();
    match_Kp_px_vec.clear();
    match_Kp_py_vec.clear();
    match_Kp_pz_vec.clear();

    match_Km_isIsolatedChargedHadron_vec.clear();
    match_Km_charge_vec.clear();
    match_Km_eta_vec.clear();
    match_Km_phi_vec.clear();
    match_Km_vx_vec.clear();
    match_Km_vy_vec.clear();
    match_Km_vz_vec.clear();
    match_Km_p_vec.clear();
    match_Km_pt_vec.clear();
    match_Km_px_vec.clear();
    match_Km_py_vec.clear();
    match_Km_pz_vec.clear();

    match_pi_isIsolatedChargedHadron_vec.clear();
    match_pi_charge_vec.clear();
    match_pi_eta_vec.clear();
    match_pi_phi_vec.clear();
    match_pi_vx_vec.clear();
    match_pi_vy_vec.clear();
    match_pi_vz_vec.clear();
    match_pi_p_vec.clear();
    match_pi_pt_vec.clear();
    match_pi_px_vec.clear();
    match_pi_py_vec.clear();
    match_pi_pz_vec.clear();

    match_phi_eta_vec.clear();
    match_phi_phi_vec.clear();
    match_phi_p_vec.clear();
    match_phi_pt_vec.clear();
    match_phi_px_vec.clear();
    match_phi_py_vec.clear();
    match_phi_pz_vec.clear();
    match_phi_invm_vec.clear();

    match_Ds_eta_vec.clear();
    match_Ds_phi_vec.clear();
    match_Ds_p_vec.clear();
    match_Ds_pt_vec.clear();
    match_Ds_px_vec.clear();
    match_Ds_py_vec.clear();
    match_Ds_pz_vec.clear();
    match_Ds_invm_vec.clear();

    match_dR_Kp_Km_vec.clear();
    match_dR_Kp_phi_vec.clear();
    match_dR_Km_phi_vec.clear();
    match_dR_Kp_pi_vec.clear();
    match_dR_Km_pi_vec.clear();
    match_dR_pi_phi_vec.clear();
    match_dR_Kp_Ds_vec.clear();
    match_dR_Km_Ds_vec.clear();
    match_dR_phi_Ds_vec.clear();
    match_dR_pi_Ds_vec.clear();

    match_dxy_phi_Ds_vec.clear();
    match_dz_phi_Ds_vec.clear();

    // Fit on phi
    match_phiFit_chi2_vec.clear();
    match_phiFit_ndof_vec.clear();
    match_phiFit_chi2ndof_vec.clear();
    match_phiFit_vx_vec.clear();
    match_phiFit_vy_vec.clear();
    match_phiFit_vz_vec.clear();
    match_phiFit_vxerr_vec.clear();
    match_phiFit_vyerr_vec.clear();
    match_phiFit_vzerr_vec.clear();

    match_phiFit_Kp_eta_vec.clear();
    match_phiFit_Kp_phi_vec.clear();
    match_phiFit_Kp_p_vec.clear();
    match_phiFit_Kp_pt_vec.clear();
    match_phiFit_Kp_px_vec.clear();
    match_phiFit_Kp_py_vec.clear();
    match_phiFit_Kp_pz_vec.clear();

    match_phiFit_Km_eta_vec.clear();
    match_phiFit_Km_phi_vec.clear();
    match_phiFit_Km_p_vec.clear();
    match_phiFit_Km_pt_vec.clear();
    match_phiFit_Km_px_vec.clear();
    match_phiFit_Km_py_vec.clear();
    match_phiFit_Km_pz_vec.clear();

    match_phiFit_pi_eta_vec.clear();
    match_phiFit_pi_phi_vec.clear();
    match_phiFit_pi_p_vec.clear();
    match_phiFit_pi_pt_vec.clear();
    match_phiFit_pi_px_vec.clear();
    match_phiFit_pi_py_vec.clear();
    match_phiFit_pi_pz_vec.clear();

    match_phiFit_phi_eta_vec.clear();
    match_phiFit_phi_phi_vec.clear();
    match_phiFit_phi_p_vec.clear();
    match_phiFit_phi_pt_vec.clear();
    match_phiFit_phi_px_vec.clear();
    match_phiFit_phi_py_vec.clear();
    match_phiFit_phi_pz_vec.clear();
    match_phiFit_phi_invm_vec.clear();

    match_phiFit_Ds_eta_vec.clear();
    match_phiFit_Ds_phi_vec.clear();
    match_phiFit_Ds_p_vec.clear();
    match_phiFit_Ds_pt_vec.clear();
    match_phiFit_Ds_px_vec.clear();
    match_phiFit_Ds_py_vec.clear();
    match_phiFit_Ds_pz_vec.clear();
    match_phiFit_Ds_invm_vec.clear();

    match_phiFit_dR_Kp_Km_vec.clear();
    match_phiFit_dR_Kp_phi_vec.clear();
    match_phiFit_dR_Km_phi_vec.clear();
    match_phiFit_dR_Kp_pi_vec.clear();
    match_phiFit_dR_Km_pi_vec.clear();
    match_phiFit_dR_pi_phi_vec.clear();
    match_phiFit_dR_Kp_Ds_vec.clear();
    match_phiFit_dR_Km_Ds_vec.clear();
    match_phiFit_dR_phi_Ds_vec.clear();
    match_phiFit_dR_pi_Ds_vec.clear();

    // Fit on Ds 
    match_DsFit_chi2_vec.clear();
    match_DsFit_ndof_vec.clear();
    match_DsFit_chi2ndof_vec.clear();
    match_DsFit_vx_vec.clear();
    match_DsFit_vy_vec.clear();
    match_DsFit_vz_vec.clear();
    match_DsFit_vxerr_vec.clear();
    match_DsFit_vyerr_vec.clear();
    match_DsFit_vzerr_vec.clear();

    match_DsFit_Kp_eta_vec.clear();
    match_DsFit_Kp_phi_vec.clear();
    match_DsFit_Kp_p_vec.clear();
    match_DsFit_Kp_pt_vec.clear();
    match_DsFit_Kp_px_vec.clear();
    match_DsFit_Kp_py_vec.clear();
    match_DsFit_Kp_pz_vec.clear();

    match_DsFit_Km_eta_vec.clear();
    match_DsFit_Km_phi_vec.clear();
    match_DsFit_Km_p_vec.clear();
    match_DsFit_Km_pt_vec.clear();
    match_DsFit_Km_px_vec.clear();
    match_DsFit_Km_py_vec.clear();
    match_DsFit_Km_pz_vec.clear();

    match_DsFit_pi_eta_vec.clear();
    match_DsFit_pi_phi_vec.clear();
    match_DsFit_pi_p_vec.clear();
    match_DsFit_pi_pt_vec.clear();
    match_DsFit_pi_px_vec.clear();
    match_DsFit_pi_py_vec.clear();
    match_DsFit_pi_pz_vec.clear();

    match_DsFit_phi_eta_vec.clear();
    match_DsFit_phi_phi_vec.clear();
    match_DsFit_phi_p_vec.clear();
    match_DsFit_phi_pt_vec.clear();
    match_DsFit_phi_px_vec.clear();
    match_DsFit_phi_py_vec.clear();
    match_DsFit_phi_pz_vec.clear();
    match_DsFit_phi_invm_vec.clear();

    match_DsFit_Ds_eta_vec.clear();
    match_DsFit_Ds_phi_vec.clear();
    match_DsFit_Ds_p_vec.clear();
    match_DsFit_Ds_pt_vec.clear();
    match_DsFit_Ds_px_vec.clear();
    match_DsFit_Ds_py_vec.clear();
    match_DsFit_Ds_pz_vec.clear();
    match_DsFit_Ds_invm_vec.clear();

    match_DsFit_dR_Kp_Km_vec.clear();
    match_DsFit_dR_Kp_phi_vec.clear();
    match_DsFit_dR_Km_phi_vec.clear();
    match_DsFit_dR_Kp_pi_vec.clear();
    match_DsFit_dR_Km_pi_vec.clear();
    match_DsFit_dR_pi_phi_vec.clear();
    match_DsFit_dR_Kp_Ds_vec.clear();
    match_DsFit_dR_Km_Ds_vec.clear();
    match_DsFit_dR_phi_Ds_vec.clear();
    match_DsFit_dR_pi_Ds_vec.clear();

    match_DsFit_Mconstraint_Ds_invm_vec.clear();

    match_Ds_primvtx_FDxy_vec.clear();
    match_Ds_primvtx_FDz_vec.clear();
    match_Ds_primvtx_FD_vec.clear();
    match_Ds_primvtx_FDxyerr_vec.clear();
    match_Ds_primvtx_FDxychi2_vec.clear();
    match_Ds_primvtx_FDzerr_vec.clear();
    match_Ds_primvtx_FDzchi2_vec.clear();
    match_Ds_primvtx_FDerr_vec.clear();
    match_Ds_primvtx_FDchi2_vec.clear();
    match_Ds_primvtx_dira_vec.clear();
    match_Ds_primvtx_dira_angle_vec.clear();
    match_Kp_primvtx_ip_vec.clear();
    match_Kp_primvtx_iperr_vec.clear();
    match_Kp_primvtx_ipchi2_vec.clear();
    match_Km_primvtx_ip_vec.clear();
    match_Km_primvtx_iperr_vec.clear();
    match_Km_primvtx_ipchi2_vec.clear();
    match_pi_primvtx_ip_vec.clear();
    match_pi_primvtx_iperr_vec.clear();
    match_pi_primvtx_ipchi2_vec.clear();
    match_phi_primvtx_ip_vec.clear();
    match_phi_primvtx_iperr_vec.clear();
    match_phi_primvtx_ipchi2_vec.clear();
    match_Ds_primvtx_ip_vec.clear();
    match_Ds_primvtx_iperr_vec.clear();
    match_Ds_primvtx_ipchi2_vec.clear();

    match_Ds_IsoR03_sumChargedHadronPt_vec.clear();
    match_Ds_IsoR03_sumNeutralHadronPt_vec.clear();
    match_Ds_IsoR03_sumPhotonPt_vec.clear();
    match_Ds_IsoR03_sumPUPt_vec.clear();
    match_Ds_PFIsoR03_vec.clear();
    match_Ds_IsoR04_sumChargedHadronPt_vec.clear();
    match_Ds_IsoR04_sumNeutralHadronPt_vec.clear();
    match_Ds_IsoR04_sumPhotonPt_vec.clear();
    match_Ds_IsoR04_sumPUPt_vec.clear();
    match_Ds_PFIsoR04_vec.clear();

    match_PV_noDs_IsValid_vec.clear();
    match_PV_noDs_IsFake_vec.clear();
    match_PV_noDs_chi2_vec.clear();
    match_PV_noDs_ndof_vec.clear();
    match_PV_noDs_chi2ndof_vec.clear();
    match_PV_noDs_x_vec.clear();
    match_PV_noDs_y_vec.clear();
    match_PV_noDs_z_vec.clear();
    match_PV_noDs_xerr_vec.clear();
    match_PV_noDs_yerr_vec.clear();
    match_PV_noDs_zerr_vec.clear();

    match_Ds_PVnoDs_FDxy_vec.clear();
    match_Ds_PVnoDs_FDz_vec.clear();
    match_Ds_PVnoDs_FD_vec.clear();
    match_Ds_PVnoDs_FDxyerr_vec.clear();
    match_Ds_PVnoDs_FDxychi2_vec.clear();
    match_Ds_PVnoDs_FDzerr_vec.clear();
    match_Ds_PVnoDs_FDzchi2_vec.clear();
    match_Ds_PVnoDs_FDerr_vec.clear();
    match_Ds_PVnoDs_FDchi2_vec.clear();
    match_Ds_PVnoDs_dira_vec.clear();
    match_Ds_PVnoDs_dira_angle_vec.clear();
    match_Kp_PVnoDs_ip_vec.clear();
    match_Kp_PVnoDs_iperr_vec.clear();
    match_Kp_PVnoDs_ipchi2_vec.clear();
    match_Km_PVnoDs_ip_vec.clear();
    match_Km_PVnoDs_iperr_vec.clear();
    match_Km_PVnoDs_ipchi2_vec.clear();
    match_pi_PVnoDs_ip_vec.clear();
    match_pi_PVnoDs_iperr_vec.clear();
    match_pi_PVnoDs_ipchi2_vec.clear();
    match_phi_PVnoDs_ip_vec.clear();
    match_phi_PVnoDs_iperr_vec.clear();
    match_phi_PVnoDs_ipchi2_vec.clear();
    match_Ds_PVnoDs_ip_vec.clear();
    match_Ds_PVnoDs_iperr_vec.clear();
    match_Ds_PVnoDs_ipchi2_vec.clear();

    match_PV_withDs_IsValid_vec.clear();
    match_PV_withDs_IsFake_vec.clear();
    match_PV_withDs_chi2_vec.clear();
    match_PV_withDs_ndof_vec.clear();
    match_PV_withDs_chi2ndof_vec.clear();
    match_PV_withDs_x_vec.clear();
    match_PV_withDs_y_vec.clear();
    match_PV_withDs_z_vec.clear();
    match_PV_withDs_xerr_vec.clear();
    match_PV_withDs_yerr_vec.clear();
    match_PV_withDs_zerr_vec.clear();

    match_Ds_PVwithDs_FDxy_vec.clear();
    match_Ds_PVwithDs_FDz_vec.clear();
    match_Ds_PVwithDs_FD_vec.clear();
    match_Ds_PVwithDs_FDxyerr_vec.clear();
    match_Ds_PVwithDs_FDxychi2_vec.clear();
    match_Ds_PVwithDs_FDzerr_vec.clear();
    match_Ds_PVwithDs_FDzchi2_vec.clear();
    match_Ds_PVwithDs_FDerr_vec.clear();
    match_Ds_PVwithDs_FDchi2_vec.clear();
    match_Ds_PVwithDs_dira_vec.clear();
    match_Ds_PVwithDs_dira_angle_vec.clear();
    match_Kp_PVwithDs_ip_vec.clear();
    match_Kp_PVwithDs_iperr_vec.clear();
    match_Kp_PVwithDs_ipchi2_vec.clear();
    match_Km_PVwithDs_ip_vec.clear();
    match_Km_PVwithDs_iperr_vec.clear();
    match_Km_PVwithDs_ipchi2_vec.clear();
    match_pi_PVwithDs_ip_vec.clear();
    match_pi_PVwithDs_iperr_vec.clear();
    match_pi_PVwithDs_ipchi2_vec.clear();
    match_phi_PVwithDs_ip_vec.clear();
    match_phi_PVwithDs_iperr_vec.clear();
    match_phi_PVwithDs_ipchi2_vec.clear();
    match_Ds_PVwithDs_ip_vec.clear();
    match_Ds_PVwithDs_iperr_vec.clear();
    match_Ds_PVwithDs_ipchi2_vec.clear();
 
    match_mu_charge_vec.clear();
    match_mu_eta_vec.clear();
    match_mu_phi_vec.clear();
    match_mu_vx_vec.clear();
    match_mu_vy_vec.clear();
    match_mu_vz_vec.clear();
    match_mu_p_vec.clear();
    match_mu_pt_vec.clear();
    match_mu_px_vec.clear();
    match_mu_py_vec.clear();
    match_mu_pz_vec.clear();
    match_mu_isHighPt_vec.clear();
    match_mu_isLoose_vec.clear();
    match_mu_isMedium_vec.clear();
    match_mu_isSoft_vec.clear();
    match_mu_isTight_vec.clear();
    match_mu_isPF_vec.clear();
    match_mu_isTracker_vec.clear();
    match_mu_isGlobal_vec.clear();
    match_mu_IsoR03_sumChargedHadronPt_vec.clear();
    match_mu_IsoR03_sumChargedParticlePt_vec.clear();
    match_mu_IsoR03_sumNeutralHadronEt_vec.clear();
    match_mu_IsoR03_sumPhotonEt_vec.clear();
    match_mu_IsoR03_sumPUPt_vec.clear();
    match_mu_PFIsoR03_vec.clear();
    match_mu_IsoR04_sumChargedHadronPt_vec.clear();
    match_mu_IsoR04_sumChargedParticlePt_vec.clear();
    match_mu_IsoR04_sumNeutralHadronEt_vec.clear();
    match_mu_IsoR04_sumPhotonEt_vec.clear();
    match_mu_IsoR04_sumPUPt_vec.clear();
    match_mu_PFIsoR04_vec.clear();
    match_mu_primvtx_dxy_vec.clear();
    match_mu_primvtx_dxyerr_vec.clear();
    match_mu_primvtx_dz_vec.clear();
    match_mu_primvtx_dzerr_vec.clear();
    match_mu_primvtx_ip_vec.clear();
    match_mu_primvtx_iperr_vec.clear();
    match_mu_primvtx_ipchi2_vec.clear();
    
    num_reco_phi = 0;
    num_reco_Ds = 0;

    Kp_isIsolatedChargedHadron_vec.clear();
    Kp_charge_vec.clear();
    Kp_eta_vec.clear();
    Kp_phi_vec.clear();
    Kp_vx_vec.clear();
    Kp_vy_vec.clear();
    Kp_vz_vec.clear();
    Kp_p_vec.clear();
    Kp_pt_vec.clear();
    Kp_px_vec.clear();
    Kp_py_vec.clear();
    Kp_pz_vec.clear();

    Km_isIsolatedChargedHadron_vec.clear();
    Km_charge_vec.clear();
    Km_eta_vec.clear();
    Km_phi_vec.clear();
    Km_vx_vec.clear();
    Km_vy_vec.clear();
    Km_vz_vec.clear();
    Km_p_vec.clear();
    Km_pt_vec.clear();
    Km_px_vec.clear();
    Km_py_vec.clear();
    Km_pz_vec.clear();

    pi_isIsolatedChargedHadron_vec.clear();
    pi_charge_vec.clear();
    pi_eta_vec.clear();
    pi_phi_vec.clear();
    pi_vx_vec.clear();
    pi_vy_vec.clear();
    pi_vz_vec.clear();
    pi_p_vec.clear();
    pi_pt_vec.clear();
    pi_px_vec.clear();
    pi_py_vec.clear();
    pi_pz_vec.clear();

    phi_eta_vec.clear();
    phi_phi_vec.clear();
    phi_p_vec.clear();
    phi_pt_vec.clear();
    phi_px_vec.clear();
    phi_py_vec.clear();
    phi_pz_vec.clear();
    phi_invm_vec.clear();

    Ds_eta_vec.clear();
    Ds_phi_vec.clear();
    Ds_p_vec.clear();
    Ds_pt_vec.clear();
    Ds_px_vec.clear();
    Ds_py_vec.clear();
    Ds_pz_vec.clear();
    Ds_invm_vec.clear();

    dR_Kp_Km_vec.clear();
    dR_Kp_phi_vec.clear();
    dR_Km_phi_vec.clear();
    dR_Kp_pi_vec.clear();
    dR_Km_pi_vec.clear();
    dR_pi_phi_vec.clear();
    dR_Kp_Ds_vec.clear();
    dR_Km_Ds_vec.clear();
    dR_phi_Ds_vec.clear();
    dR_pi_Ds_vec.clear();

    dxy_phi_Ds_vec.clear();
    dz_phi_Ds_vec.clear();

    // Fit on phi
    phiFit_chi2_vec.clear();
    phiFit_ndof_vec.clear();
    phiFit_chi2ndof_vec.clear();
    phiFit_vx_vec.clear();
    phiFit_vy_vec.clear();
    phiFit_vz_vec.clear();
    phiFit_vxerr_vec.clear();
    phiFit_vyerr_vec.clear();
    phiFit_vzerr_vec.clear();

    phiFit_Kp_eta_vec.clear();
    phiFit_Kp_phi_vec.clear();
    phiFit_Kp_p_vec.clear();
    phiFit_Kp_pt_vec.clear();
    phiFit_Kp_px_vec.clear();
    phiFit_Kp_py_vec.clear();
    phiFit_Kp_pz_vec.clear();

    phiFit_Km_eta_vec.clear();
    phiFit_Km_phi_vec.clear();
    phiFit_Km_p_vec.clear();
    phiFit_Km_pt_vec.clear();
    phiFit_Km_px_vec.clear();
    phiFit_Km_py_vec.clear();
    phiFit_Km_pz_vec.clear();

    phiFit_pi_eta_vec.clear();
    phiFit_pi_phi_vec.clear();
    phiFit_pi_p_vec.clear();
    phiFit_pi_pt_vec.clear();
    phiFit_pi_px_vec.clear();
    phiFit_pi_py_vec.clear();
    phiFit_pi_pz_vec.clear();

    phiFit_phi_eta_vec.clear();
    phiFit_phi_phi_vec.clear();
    phiFit_phi_p_vec.clear();
    phiFit_phi_pt_vec.clear();
    phiFit_phi_px_vec.clear();
    phiFit_phi_py_vec.clear();
    phiFit_phi_pz_vec.clear();
    phiFit_phi_invm_vec.clear();

    phiFit_Ds_eta_vec.clear();
    phiFit_Ds_phi_vec.clear();
    phiFit_Ds_p_vec.clear();
    phiFit_Ds_pt_vec.clear();
    phiFit_Ds_px_vec.clear();
    phiFit_Ds_py_vec.clear();
    phiFit_Ds_pz_vec.clear();
    phiFit_Ds_invm_vec.clear();

    phiFit_dR_Kp_Km_vec.clear();
    phiFit_dR_Kp_phi_vec.clear();
    phiFit_dR_Km_phi_vec.clear();
    phiFit_dR_Kp_pi_vec.clear();
    phiFit_dR_Km_pi_vec.clear();
    phiFit_dR_pi_phi_vec.clear();
    phiFit_dR_Kp_Ds_vec.clear();
    phiFit_dR_Km_Ds_vec.clear();
    phiFit_dR_phi_Ds_vec.clear();
    phiFit_dR_pi_Ds_vec.clear();

    // Fit on Ds 
    DsFit_chi2_vec.clear();
    DsFit_ndof_vec.clear();
    DsFit_chi2ndof_vec.clear();
    DsFit_vx_vec.clear();
    DsFit_vy_vec.clear();
    DsFit_vz_vec.clear();
    DsFit_vxerr_vec.clear();
    DsFit_vyerr_vec.clear();
    DsFit_vzerr_vec.clear();

    DsFit_Kp_eta_vec.clear();
    DsFit_Kp_phi_vec.clear();
    DsFit_Kp_p_vec.clear();
    DsFit_Kp_pt_vec.clear();
    DsFit_Kp_px_vec.clear();
    DsFit_Kp_py_vec.clear();
    DsFit_Kp_pz_vec.clear();

    DsFit_Km_eta_vec.clear();
    DsFit_Km_phi_vec.clear();
    DsFit_Km_p_vec.clear();
    DsFit_Km_pt_vec.clear();
    DsFit_Km_px_vec.clear();
    DsFit_Km_py_vec.clear();
    DsFit_Km_pz_vec.clear();

    DsFit_pi_eta_vec.clear();
    DsFit_pi_phi_vec.clear();
    DsFit_pi_p_vec.clear();
    DsFit_pi_pt_vec.clear();
    DsFit_pi_px_vec.clear();
    DsFit_pi_py_vec.clear();
    DsFit_pi_pz_vec.clear();

    DsFit_phi_eta_vec.clear();
    DsFit_phi_phi_vec.clear();
    DsFit_phi_p_vec.clear();
    DsFit_phi_pt_vec.clear();
    DsFit_phi_px_vec.clear();
    DsFit_phi_py_vec.clear();
    DsFit_phi_pz_vec.clear();
    DsFit_phi_invm_vec.clear();

    DsFit_Ds_eta_vec.clear();
    DsFit_Ds_phi_vec.clear();
    DsFit_Ds_p_vec.clear();
    DsFit_Ds_pt_vec.clear();
    DsFit_Ds_px_vec.clear();
    DsFit_Ds_py_vec.clear();
    DsFit_Ds_pz_vec.clear();
    DsFit_Ds_invm_vec.clear();

    DsFit_dR_Kp_Km_vec.clear();
    DsFit_dR_Kp_phi_vec.clear();
    DsFit_dR_Km_phi_vec.clear();
    DsFit_dR_Kp_pi_vec.clear();
    DsFit_dR_Km_pi_vec.clear();
    DsFit_dR_pi_phi_vec.clear();
    DsFit_dR_Kp_Ds_vec.clear();
    DsFit_dR_Km_Ds_vec.clear();
    DsFit_dR_phi_Ds_vec.clear();
    DsFit_dR_pi_Ds_vec.clear();

    DsFit_Mconstraint_Ds_invm_vec.clear();

    Ds_primvtx_FDxy_vec.clear();
    Ds_primvtx_FDz_vec.clear();
    Ds_primvtx_FD_vec.clear();
    Ds_primvtx_FDxyerr_vec.clear();
    Ds_primvtx_FDxychi2_vec.clear();
    Ds_primvtx_FDzerr_vec.clear();
    Ds_primvtx_FDzchi2_vec.clear();
    Ds_primvtx_FDerr_vec.clear();
    Ds_primvtx_FDchi2_vec.clear();
    Ds_primvtx_dira_vec.clear();
    Ds_primvtx_dira_angle_vec.clear();
    Kp_primvtx_ip_vec.clear();
    Kp_primvtx_iperr_vec.clear();
    Kp_primvtx_ipchi2_vec.clear();
    Km_primvtx_ip_vec.clear();
    Km_primvtx_iperr_vec.clear();
    Km_primvtx_ipchi2_vec.clear();
    pi_primvtx_ip_vec.clear();
    pi_primvtx_iperr_vec.clear();
    pi_primvtx_ipchi2_vec.clear();
    phi_primvtx_ip_vec.clear();
    phi_primvtx_iperr_vec.clear();
    phi_primvtx_ipchi2_vec.clear();
    Ds_primvtx_ip_vec.clear();
    Ds_primvtx_iperr_vec.clear();
    Ds_primvtx_ipchi2_vec.clear();

    Ds_IsoR03_sumChargedHadronPt_vec.clear();
    Ds_IsoR03_sumNeutralHadronPt_vec.clear();
    Ds_IsoR03_sumPhotonPt_vec.clear();
    Ds_IsoR03_sumPUPt_vec.clear();
    Ds_PFIsoR03_vec.clear();
    Ds_IsoR04_sumChargedHadronPt_vec.clear();
    Ds_IsoR04_sumNeutralHadronPt_vec.clear();
    Ds_IsoR04_sumPhotonPt_vec.clear();
    Ds_IsoR04_sumPUPt_vec.clear();
    Ds_PFIsoR04_vec.clear();

    PV_noDs_IsValid_vec.clear();
    PV_noDs_IsFake_vec.clear();
    PV_noDs_chi2_vec.clear();
    PV_noDs_ndof_vec.clear();
    PV_noDs_chi2ndof_vec.clear();
    PV_noDs_x_vec.clear();
    PV_noDs_y_vec.clear();
    PV_noDs_z_vec.clear();
    PV_noDs_xerr_vec.clear();
    PV_noDs_yerr_vec.clear();
    PV_noDs_zerr_vec.clear();

    Ds_PVnoDs_FDxy_vec.clear();
    Ds_PVnoDs_FDz_vec.clear();
    Ds_PVnoDs_FD_vec.clear();
    Ds_PVnoDs_FDxyerr_vec.clear();
    Ds_PVnoDs_FDxychi2_vec.clear();
    Ds_PVnoDs_FDzerr_vec.clear();
    Ds_PVnoDs_FDzchi2_vec.clear();
    Ds_PVnoDs_FDerr_vec.clear();
    Ds_PVnoDs_FDchi2_vec.clear();
    Ds_PVnoDs_dira_vec.clear();
    Ds_PVnoDs_dira_angle_vec.clear();
    Kp_PVnoDs_ip_vec.clear();
    Kp_PVnoDs_iperr_vec.clear();
    Kp_PVnoDs_ipchi2_vec.clear();
    Km_PVnoDs_ip_vec.clear();
    Km_PVnoDs_iperr_vec.clear();
    Km_PVnoDs_ipchi2_vec.clear();
    pi_PVnoDs_ip_vec.clear();
    pi_PVnoDs_iperr_vec.clear();
    pi_PVnoDs_ipchi2_vec.clear();
    phi_PVnoDs_ip_vec.clear();
    phi_PVnoDs_iperr_vec.clear();
    phi_PVnoDs_ipchi2_vec.clear();
    Ds_PVnoDs_ip_vec.clear();
    Ds_PVnoDs_iperr_vec.clear();
    Ds_PVnoDs_ipchi2_vec.clear();

    PV_withDs_IsValid_vec.clear();
    PV_withDs_IsFake_vec.clear();
    PV_withDs_chi2_vec.clear();
    PV_withDs_ndof_vec.clear();
    PV_withDs_chi2ndof_vec.clear();
    PV_withDs_x_vec.clear();
    PV_withDs_y_vec.clear();
    PV_withDs_z_vec.clear();
    PV_withDs_xerr_vec.clear();
    PV_withDs_yerr_vec.clear();
    PV_withDs_zerr_vec.clear();

    Ds_PVwithDs_FDxy_vec.clear();
    Ds_PVwithDs_FDz_vec.clear();
    Ds_PVwithDs_FD_vec.clear();
    Ds_PVwithDs_FDxyerr_vec.clear();
    Ds_PVwithDs_FDxychi2_vec.clear();
    Ds_PVwithDs_FDzerr_vec.clear();
    Ds_PVwithDs_FDzchi2_vec.clear();
    Ds_PVwithDs_FDerr_vec.clear();
    Ds_PVwithDs_FDchi2_vec.clear();
    Ds_PVwithDs_dira_vec.clear();
    Ds_PVwithDs_dira_angle_vec.clear();
    Kp_PVwithDs_ip_vec.clear();
    Kp_PVwithDs_iperr_vec.clear();
    Kp_PVwithDs_ipchi2_vec.clear();
    Km_PVwithDs_ip_vec.clear();
    Km_PVwithDs_iperr_vec.clear();
    Km_PVwithDs_ipchi2_vec.clear();
    pi_PVwithDs_ip_vec.clear();
    pi_PVwithDs_iperr_vec.clear();
    pi_PVwithDs_ipchi2_vec.clear();
    phi_PVwithDs_ip_vec.clear();
    phi_PVwithDs_iperr_vec.clear();
    phi_PVwithDs_ipchi2_vec.clear();
    Ds_PVwithDs_ip_vec.clear();
    Ds_PVwithDs_iperr_vec.clear();
    Ds_PVwithDs_ipchi2_vec.clear();

    Kp_match_vec.clear();
    Km_match_vec.clear();
    pi_match_vec.clear();
    match_entry_vec.clear();
    non_match_entry_vec.clear();

    best_Kp_isIsolatedChargedHadron_vec.clear();
    best_Kp_charge_vec.clear();
    best_Kp_eta_vec.clear();
    best_Kp_phi_vec.clear();
    best_Kp_vx_vec.clear();
    best_Kp_vy_vec.clear();
    best_Kp_vz_vec.clear();
    best_Kp_p_vec.clear();
    best_Kp_pt_vec.clear();
    best_Kp_px_vec.clear();
    best_Kp_py_vec.clear();
    best_Kp_pz_vec.clear();

    best_Km_isIsolatedChargedHadron_vec.clear();
    best_Km_charge_vec.clear();
    best_Km_eta_vec.clear();
    best_Km_phi_vec.clear();
    best_Km_vx_vec.clear();
    best_Km_vy_vec.clear();
    best_Km_vz_vec.clear();
    best_Km_p_vec.clear();
    best_Km_pt_vec.clear();
    best_Km_px_vec.clear();
    best_Km_py_vec.clear();
    best_Km_pz_vec.clear();

    best_pi_isIsolatedChargedHadron_vec.clear();
    best_pi_charge_vec.clear();
    best_pi_eta_vec.clear();
    best_pi_phi_vec.clear();
    best_pi_vx_vec.clear();
    best_pi_vy_vec.clear();
    best_pi_vz_vec.clear();
    best_pi_p_vec.clear();
    best_pi_pt_vec.clear();
    best_pi_px_vec.clear();
    best_pi_py_vec.clear();
    best_pi_pz_vec.clear();

    best_phi_eta_vec.clear();
    best_phi_phi_vec.clear();
    best_phi_p_vec.clear();
    best_phi_pt_vec.clear();
    best_phi_px_vec.clear();
    best_phi_py_vec.clear();
    best_phi_pz_vec.clear();
    best_phi_invm_vec.clear();

    best_Ds_eta_vec.clear();
    best_Ds_phi_vec.clear();
    best_Ds_p_vec.clear();
    best_Ds_pt_vec.clear();
    best_Ds_px_vec.clear();
    best_Ds_py_vec.clear();
    best_Ds_pz_vec.clear();
    best_Ds_invm_vec.clear();

    best_dR_Kp_Km_vec.clear();
    best_dR_Kp_phi_vec.clear();
    best_dR_Km_phi_vec.clear();
    best_dR_Kp_pi_vec.clear();
    best_dR_Km_pi_vec.clear();
    best_dR_pi_phi_vec.clear();
    best_dR_Kp_Ds_vec.clear();
    best_dR_Km_Ds_vec.clear();
    best_dR_phi_Ds_vec.clear();
    best_dR_pi_Ds_vec.clear();

    best_dxy_phi_Ds_vec.clear();
    best_dz_phi_Ds_vec.clear();

    // Fit on phi
    best_phiFit_chi2_vec.clear();
    best_phiFit_ndof_vec.clear();
    best_phiFit_chi2ndof_vec.clear();
    best_phiFit_vx_vec.clear();
    best_phiFit_vy_vec.clear();
    best_phiFit_vz_vec.clear();
    best_phiFit_vxerr_vec.clear();
    best_phiFit_vyerr_vec.clear();
    best_phiFit_vzerr_vec.clear();

    best_phiFit_Kp_eta_vec.clear();
    best_phiFit_Kp_phi_vec.clear();
    best_phiFit_Kp_p_vec.clear();
    best_phiFit_Kp_pt_vec.clear();
    best_phiFit_Kp_px_vec.clear();
    best_phiFit_Kp_py_vec.clear();
    best_phiFit_Kp_pz_vec.clear();

    best_phiFit_Km_eta_vec.clear();
    best_phiFit_Km_phi_vec.clear();
    best_phiFit_Km_p_vec.clear();
    best_phiFit_Km_pt_vec.clear();
    best_phiFit_Km_px_vec.clear();
    best_phiFit_Km_py_vec.clear();
    best_phiFit_Km_pz_vec.clear();

    best_phiFit_pi_eta_vec.clear();
    best_phiFit_pi_phi_vec.clear();
    best_phiFit_pi_p_vec.clear();
    best_phiFit_pi_pt_vec.clear();
    best_phiFit_pi_px_vec.clear();
    best_phiFit_pi_py_vec.clear();
    best_phiFit_pi_pz_vec.clear();

    best_phiFit_phi_eta_vec.clear();
    best_phiFit_phi_phi_vec.clear();
    best_phiFit_phi_p_vec.clear();
    best_phiFit_phi_pt_vec.clear();
    best_phiFit_phi_px_vec.clear();
    best_phiFit_phi_py_vec.clear();
    best_phiFit_phi_pz_vec.clear();
    best_phiFit_phi_invm_vec.clear();

    best_phiFit_Ds_eta_vec.clear();
    best_phiFit_Ds_phi_vec.clear();
    best_phiFit_Ds_p_vec.clear();
    best_phiFit_Ds_pt_vec.clear();
    best_phiFit_Ds_px_vec.clear();
    best_phiFit_Ds_py_vec.clear();
    best_phiFit_Ds_pz_vec.clear();
    best_phiFit_Ds_invm_vec.clear();

    best_phiFit_dR_Kp_Km_vec.clear();
    best_phiFit_dR_Kp_phi_vec.clear();
    best_phiFit_dR_Km_phi_vec.clear();
    best_phiFit_dR_Kp_pi_vec.clear();
    best_phiFit_dR_Km_pi_vec.clear();
    best_phiFit_dR_pi_phi_vec.clear();
    best_phiFit_dR_Kp_Ds_vec.clear();
    best_phiFit_dR_Km_Ds_vec.clear();
    best_phiFit_dR_phi_Ds_vec.clear();
    best_phiFit_dR_pi_Ds_vec.clear();

    // Fit on Ds 
    best_DsFit_chi2_vec.clear();
    best_DsFit_ndof_vec.clear();
    best_DsFit_chi2ndof_vec.clear();
    best_DsFit_vx_vec.clear();
    best_DsFit_vy_vec.clear();
    best_DsFit_vz_vec.clear();
    best_DsFit_vxerr_vec.clear();
    best_DsFit_vyerr_vec.clear();
    best_DsFit_vzerr_vec.clear();

    best_DsFit_Kp_eta_vec.clear();
    best_DsFit_Kp_phi_vec.clear();
    best_DsFit_Kp_p_vec.clear();
    best_DsFit_Kp_pt_vec.clear();
    best_DsFit_Kp_px_vec.clear();
    best_DsFit_Kp_py_vec.clear();
    best_DsFit_Kp_pz_vec.clear();

    best_DsFit_Km_eta_vec.clear();
    best_DsFit_Km_phi_vec.clear();
    best_DsFit_Km_p_vec.clear();
    best_DsFit_Km_pt_vec.clear();
    best_DsFit_Km_px_vec.clear();
    best_DsFit_Km_py_vec.clear();
    best_DsFit_Km_pz_vec.clear();

    best_DsFit_pi_eta_vec.clear();
    best_DsFit_pi_phi_vec.clear();
    best_DsFit_pi_p_vec.clear();
    best_DsFit_pi_pt_vec.clear();
    best_DsFit_pi_px_vec.clear();
    best_DsFit_pi_py_vec.clear();
    best_DsFit_pi_pz_vec.clear();

    best_DsFit_phi_eta_vec.clear();
    best_DsFit_phi_phi_vec.clear();
    best_DsFit_phi_p_vec.clear();
    best_DsFit_phi_pt_vec.clear();
    best_DsFit_phi_px_vec.clear();
    best_DsFit_phi_py_vec.clear();
    best_DsFit_phi_pz_vec.clear();
    best_DsFit_phi_invm_vec.clear();

    best_DsFit_Ds_eta_vec.clear();
    best_DsFit_Ds_phi_vec.clear();
    best_DsFit_Ds_p_vec.clear();
    best_DsFit_Ds_pt_vec.clear();
    best_DsFit_Ds_px_vec.clear();
    best_DsFit_Ds_py_vec.clear();
    best_DsFit_Ds_pz_vec.clear();
    best_DsFit_Ds_invm_vec.clear();

    best_DsFit_dR_Kp_Km_vec.clear();
    best_DsFit_dR_Kp_phi_vec.clear();
    best_DsFit_dR_Km_phi_vec.clear();
    best_DsFit_dR_Kp_pi_vec.clear();
    best_DsFit_dR_Km_pi_vec.clear();
    best_DsFit_dR_pi_phi_vec.clear();
    best_DsFit_dR_Kp_Ds_vec.clear();
    best_DsFit_dR_Km_Ds_vec.clear();
    best_DsFit_dR_phi_Ds_vec.clear();
    best_DsFit_dR_pi_Ds_vec.clear();

    best_DsFit_Mconstraint_Ds_invm_vec.clear();

    best_Ds_primvtx_FDxy_vec.clear();
    best_Ds_primvtx_FDz_vec.clear();
    best_Ds_primvtx_FD_vec.clear();
    best_Ds_primvtx_FDxyerr_vec.clear();
    best_Ds_primvtx_FDxychi2_vec.clear();
    best_Ds_primvtx_FDzerr_vec.clear();
    best_Ds_primvtx_FDzchi2_vec.clear();
    best_Ds_primvtx_FDerr_vec.clear();
    best_Ds_primvtx_FDchi2_vec.clear();
    best_Ds_primvtx_dira_vec.clear();
    best_Ds_primvtx_dira_angle_vec.clear();
    best_Kp_primvtx_ip_vec.clear();
    best_Kp_primvtx_iperr_vec.clear();
    best_Kp_primvtx_ipchi2_vec.clear();
    best_Km_primvtx_ip_vec.clear();
    best_Km_primvtx_iperr_vec.clear();
    best_Km_primvtx_ipchi2_vec.clear();
    best_pi_primvtx_ip_vec.clear();
    best_pi_primvtx_iperr_vec.clear();
    best_pi_primvtx_ipchi2_vec.clear();
    best_phi_primvtx_ip_vec.clear();
    best_phi_primvtx_iperr_vec.clear();
    best_phi_primvtx_ipchi2_vec.clear();
    best_Ds_primvtx_ip_vec.clear();
    best_Ds_primvtx_iperr_vec.clear();
    best_Ds_primvtx_ipchi2_vec.clear();

    best_Ds_IsoR03_sumChargedHadronPt_vec.clear();
    best_Ds_IsoR03_sumNeutralHadronPt_vec.clear();
    best_Ds_IsoR03_sumPhotonPt_vec.clear();
    best_Ds_IsoR03_sumPUPt_vec.clear();
    best_Ds_PFIsoR03_vec.clear();
    best_Ds_IsoR04_sumChargedHadronPt_vec.clear();
    best_Ds_IsoR04_sumNeutralHadronPt_vec.clear();
    best_Ds_IsoR04_sumPhotonPt_vec.clear();
    best_Ds_IsoR04_sumPUPt_vec.clear();
    best_Ds_PFIsoR04_vec.clear();

    best_PV_noDs_IsValid_vec.clear();
    best_PV_noDs_IsFake_vec.clear();
    best_PV_noDs_chi2_vec.clear();
    best_PV_noDs_ndof_vec.clear();
    best_PV_noDs_chi2ndof_vec.clear();
    best_PV_noDs_x_vec.clear();
    best_PV_noDs_y_vec.clear();
    best_PV_noDs_z_vec.clear();
    best_PV_noDs_xerr_vec.clear();
    best_PV_noDs_yerr_vec.clear();
    best_PV_noDs_zerr_vec.clear();

    best_Ds_PVnoDs_FDxy_vec.clear();
    best_Ds_PVnoDs_FDz_vec.clear();
    best_Ds_PVnoDs_FD_vec.clear();
    best_Ds_PVnoDs_FDxyerr_vec.clear();
    best_Ds_PVnoDs_FDxychi2_vec.clear();
    best_Ds_PVnoDs_FDzerr_vec.clear();
    best_Ds_PVnoDs_FDzchi2_vec.clear();
    best_Ds_PVnoDs_FDerr_vec.clear();
    best_Ds_PVnoDs_FDchi2_vec.clear();
    best_Ds_PVnoDs_dira_vec.clear();
    best_Ds_PVnoDs_dira_angle_vec.clear();
    best_Kp_PVnoDs_ip_vec.clear();
    best_Kp_PVnoDs_iperr_vec.clear();
    best_Kp_PVnoDs_ipchi2_vec.clear();
    best_Km_PVnoDs_ip_vec.clear();
    best_Km_PVnoDs_iperr_vec.clear();
    best_Km_PVnoDs_ipchi2_vec.clear();
    best_pi_PVnoDs_ip_vec.clear();
    best_pi_PVnoDs_iperr_vec.clear();
    best_pi_PVnoDs_ipchi2_vec.clear();
    best_phi_PVnoDs_ip_vec.clear();
    best_phi_PVnoDs_iperr_vec.clear();
    best_phi_PVnoDs_ipchi2_vec.clear();
    best_Ds_PVnoDs_ip_vec.clear();
    best_Ds_PVnoDs_iperr_vec.clear();
    best_Ds_PVnoDs_ipchi2_vec.clear();

    best_PV_withDs_IsValid_vec.clear();
    best_PV_withDs_IsFake_vec.clear();
    best_PV_withDs_chi2_vec.clear();
    best_PV_withDs_ndof_vec.clear();
    best_PV_withDs_chi2ndof_vec.clear();
    best_PV_withDs_x_vec.clear();
    best_PV_withDs_y_vec.clear();
    best_PV_withDs_z_vec.clear();
    best_PV_withDs_xerr_vec.clear();
    best_PV_withDs_yerr_vec.clear();
    best_PV_withDs_zerr_vec.clear();

    best_Ds_PVwithDs_FDxy_vec.clear();
    best_Ds_PVwithDs_FDz_vec.clear();
    best_Ds_PVwithDs_FD_vec.clear();
    best_Ds_PVwithDs_FDxyerr_vec.clear();
    best_Ds_PVwithDs_FDxychi2_vec.clear();
    best_Ds_PVwithDs_FDzerr_vec.clear();
    best_Ds_PVwithDs_FDzchi2_vec.clear();
    best_Ds_PVwithDs_FDerr_vec.clear();
    best_Ds_PVwithDs_FDchi2_vec.clear();
    best_Ds_PVwithDs_dira_vec.clear();
    best_Ds_PVwithDs_dira_angle_vec.clear();
    best_Kp_PVwithDs_ip_vec.clear();
    best_Kp_PVwithDs_iperr_vec.clear();
    best_Kp_PVwithDs_ipchi2_vec.clear();
    best_Km_PVwithDs_ip_vec.clear();
    best_Km_PVwithDs_iperr_vec.clear();
    best_Km_PVwithDs_ipchi2_vec.clear();
    best_pi_PVwithDs_ip_vec.clear();
    best_pi_PVwithDs_iperr_vec.clear();
    best_pi_PVwithDs_ipchi2_vec.clear();
    best_phi_PVwithDs_ip_vec.clear();
    best_phi_PVwithDs_iperr_vec.clear();
    best_phi_PVwithDs_ipchi2_vec.clear();
    best_Ds_PVwithDs_ip_vec.clear();
    best_Ds_PVwithDs_iperr_vec.clear();
    best_Ds_PVwithDs_ipchi2_vec.clear();

    best_match_entry_vec.clear();
}

void PVTree::CreateBranches()
{
    tree->Branch("num_Gen_Kp", &num_Gen_Kp);
    tree->Branch("num_Gen_Km", &num_Gen_Km);
    tree->Branch("num_Gen_pi", &num_Gen_pi);
    tree->Branch("num_Gen_phi", &num_Gen_phi);
    tree->Branch("num_Gen_Ds", &num_Gen_Ds);
    tree->Branch("num_Gen_mu", &num_Gen_mu);
    tree->Branch("num_Gen_nu", &num_Gen_nu);
    tree->Branch("num_Gen_W", &num_Gen_W);
    tree->Branch("num_Gen_H", &num_Gen_H);

    tree->Branch("Gen_Kp_eta", &Gen_Kp_eta);
    tree->Branch("Gen_Kp_phi", &Gen_Kp_phi);
    tree->Branch("Gen_Kp_vx", &Gen_Kp_vx);
    tree->Branch("Gen_Kp_vy", &Gen_Kp_vy);
    tree->Branch("Gen_Kp_vz", &Gen_Kp_vz);
    tree->Branch("Gen_Kp_p", &Gen_Kp_p);
    tree->Branch("Gen_Kp_pt", &Gen_Kp_pt);
    tree->Branch("Gen_Kp_px", &Gen_Kp_px);
    tree->Branch("Gen_Kp_py", &Gen_Kp_py);
    tree->Branch("Gen_Kp_pz", &Gen_Kp_pz);
    tree->Branch("Gen_Kp_pp", &Gen_Kp_pp);
    tree->Branch("Gen_Kp_pl", &Gen_Kp_pl);

    tree->Branch("Gen_Km_eta", &Gen_Km_eta);
    tree->Branch("Gen_Km_phi", &Gen_Km_phi);
    tree->Branch("Gen_Km_vx", &Gen_Km_vx);
    tree->Branch("Gen_Km_vy", &Gen_Km_vy);
    tree->Branch("Gen_Km_vz", &Gen_Km_vz);
    tree->Branch("Gen_Km_p", &Gen_Km_p);
    tree->Branch("Gen_Km_pt", &Gen_Km_pt);
    tree->Branch("Gen_Km_px", &Gen_Km_px);
    tree->Branch("Gen_Km_py", &Gen_Km_py);
    tree->Branch("Gen_Km_pz", &Gen_Km_pz);
    tree->Branch("Gen_Km_pp", &Gen_Km_pp);
    tree->Branch("Gen_Km_pl", &Gen_Km_pl);

    tree->Branch("Gen_pi_eta", &Gen_pi_eta);
    tree->Branch("Gen_pi_phi", &Gen_pi_phi);
    tree->Branch("Gen_pi_vx", &Gen_pi_vx);
    tree->Branch("Gen_pi_vy", &Gen_pi_vy);
    tree->Branch("Gen_pi_vz", &Gen_pi_vz);
    tree->Branch("Gen_pi_p", &Gen_pi_p);
    tree->Branch("Gen_pi_pt", &Gen_pi_pt);
    tree->Branch("Gen_pi_px", &Gen_pi_px);
    tree->Branch("Gen_pi_py", &Gen_pi_py);
    tree->Branch("Gen_pi_pz", &Gen_pi_pz);
    tree->Branch("Gen_pi_pp", &Gen_pi_pp);
    tree->Branch("Gen_pi_pl", &Gen_pi_pl);

    tree->Branch("Gen_phi_eta", &Gen_phi_eta);
    tree->Branch("Gen_phi_phi", &Gen_phi_phi);
    tree->Branch("Gen_phi_vx", &Gen_phi_vx);
    tree->Branch("Gen_phi_vy", &Gen_phi_vy);
    tree->Branch("Gen_phi_vz", &Gen_phi_vz);
    tree->Branch("Gen_phi_p", &Gen_phi_p);
    tree->Branch("Gen_phi_pt", &Gen_phi_pt);
    tree->Branch("Gen_phi_px", &Gen_phi_px);
    tree->Branch("Gen_phi_py", &Gen_phi_py);
    tree->Branch("Gen_phi_pz", &Gen_phi_pz);
    tree->Branch("Gen_phi_pp", &Gen_phi_pp);
    tree->Branch("Gen_phi_pl", &Gen_phi_pl);

    tree->Branch("Gen_Ds_eta", &Gen_Ds_eta);
    tree->Branch("Gen_Ds_phi", &Gen_Ds_phi);
    tree->Branch("Gen_Ds_vx", &Gen_Ds_vx);
    tree->Branch("Gen_Ds_vy", &Gen_Ds_vy);
    tree->Branch("Gen_Ds_vz", &Gen_Ds_vz);
    tree->Branch("Gen_Ds_p", &Gen_Ds_p);
    tree->Branch("Gen_Ds_pt", &Gen_Ds_pt);
    tree->Branch("Gen_Ds_px", &Gen_Ds_px);
    tree->Branch("Gen_Ds_py", &Gen_Ds_py);
    tree->Branch("Gen_Ds_pz", &Gen_Ds_pz);

    tree->Branch("Gen_mu_eta", &Gen_mu_eta);
    tree->Branch("Gen_mu_phi", &Gen_mu_phi);
    tree->Branch("Gen_mu_vx", &Gen_mu_vx);
    tree->Branch("Gen_mu_vy", &Gen_mu_vy);
    tree->Branch("Gen_mu_vz", &Gen_mu_vz);
    tree->Branch("Gen_mu_p", &Gen_mu_p);
    tree->Branch("Gen_mu_pt", &Gen_mu_pt);
    tree->Branch("Gen_mu_px", &Gen_mu_px);
    tree->Branch("Gen_mu_py", &Gen_mu_py);
    tree->Branch("Gen_mu_pz", &Gen_mu_pz);

    tree->Branch("Gen_nu_eta", &Gen_nu_eta);
    tree->Branch("Gen_nu_phi", &Gen_nu_phi);
    tree->Branch("Gen_nu_vx", &Gen_nu_vx);
    tree->Branch("Gen_nu_vy", &Gen_nu_vy);
    tree->Branch("Gen_nu_vz", &Gen_nu_vz);
    tree->Branch("Gen_nu_p", &Gen_nu_p);
    tree->Branch("Gen_nu_pt", &Gen_nu_pt);
    tree->Branch("Gen_nu_px", &Gen_nu_px);
    tree->Branch("Gen_nu_py", &Gen_nu_py);
    tree->Branch("Gen_nu_pz", &Gen_nu_pz);

    tree->Branch("Gen_W_eta", &Gen_W_eta);
    tree->Branch("Gen_W_phi", &Gen_W_phi);
    tree->Branch("Gen_W_vx", &Gen_W_vx);
    tree->Branch("Gen_W_vy", &Gen_W_vy);
    tree->Branch("Gen_W_vz", &Gen_W_vz);
    tree->Branch("Gen_W_p", &Gen_W_p);
    tree->Branch("Gen_W_pt", &Gen_W_pt);
    tree->Branch("Gen_W_px", &Gen_W_px);
    tree->Branch("Gen_W_py", &Gen_W_py);
    tree->Branch("Gen_W_pz", &Gen_W_pz);

    tree->Branch("Gen_H_eta", &Gen_H_eta);
    tree->Branch("Gen_H_phi", &Gen_H_phi);
    tree->Branch("Gen_H_vx", &Gen_H_vx);
    tree->Branch("Gen_H_vy", &Gen_H_vy);
    tree->Branch("Gen_H_vz", &Gen_H_vz);
    tree->Branch("Gen_H_p", &Gen_H_p);
    tree->Branch("Gen_H_pt", &Gen_H_pt);
    tree->Branch("Gen_H_px", &Gen_H_px);
    tree->Branch("Gen_H_py", &Gen_H_py);
    tree->Branch("Gen_H_pz", &Gen_H_pz);

    tree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
    tree->Branch("Gen_dR_Kp_phi", &Gen_dR_Kp_phi);
    tree->Branch("Gen_dR_Km_phi", &Gen_dR_Km_phi);
    tree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
    tree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);
    tree->Branch("Gen_dR_pi_phi", &Gen_dR_pi_phi);
    tree->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds);
    tree->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds);
    tree->Branch("Gen_dR_phi_Ds", &Gen_dR_phi_Ds);
    tree->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds);
    tree->Branch("Gen_dR_Kp_mu", &Gen_dR_Kp_mu);
    tree->Branch("Gen_dR_Km_mu", &Gen_dR_Km_mu);
    tree->Branch("Gen_dR_phi_mu", &Gen_dR_phi_mu);
    tree->Branch("Gen_dR_pi_mu", &Gen_dR_pi_mu);
    tree->Branch("Gen_dR_Ds_mu", &Gen_dR_Ds_mu);

    tree->Branch("Gen_Ds_dx", &Gen_Ds_dx);
    tree->Branch("Gen_Ds_dy", &Gen_Ds_dy);
    tree->Branch("Gen_Ds_dz", &Gen_Ds_dz);
    tree->Branch("Gen_Ds_FDxy", &Gen_Ds_FDxy);
    tree->Branch("Gen_Ds_FD", &Gen_Ds_FD);

    tree->Branch("BS_type", &BS_type);
    tree->Branch("BS_x0", &BS_x0);
    tree->Branch("BS_y0", &BS_y0);
    tree->Branch("BS_z0", &BS_z0);
    tree->Branch("BS_sigmaZ", &BS_sigmaZ);
    tree->Branch("BS_dxdz", &BS_dxdz);
    tree->Branch("BS_dydz", &BS_dydz);
    tree->Branch("BS_BWX", &BS_BWX);
    tree->Branch("BS_BWY", &BS_BWY);
    tree->Branch("BS_x0err", &BS_x0err);
    tree->Branch("BS_y0err", &BS_y0err);
    tree->Branch("BS_z0err", &BS_z0err);
    tree->Branch("BS_sigmaZ0err", &BS_sigmaZ0err);
    tree->Branch("BS_dxdzerr", &BS_dxdzerr);
    tree->Branch("BS_dydzerr", &BS_dydzerr);
    tree->Branch("BS_BWXerr", &BS_BWXerr);
    tree->Branch("BS_BWYerr", &BS_BWYerr);
    tree->Branch("BS_emitX", &BS_emitX);
    tree->Branch("BS_emitY", &BS_emitY);
    tree->Branch("BS_betaStar", &BS_betaStar);
   
    tree->Branch("PV_vx", &PV_vx); 
    tree->Branch("PV_vy", &PV_vy); 
    tree->Branch("PV_vz", &PV_vz); 
    tree->Branch("PV_vxerr", &PV_vxerr); 
    tree->Branch("PV_vyerr", &PV_vyerr); 
    tree->Branch("PV_vzerr", &PV_vzerr); 
    
    tree->Branch("num_match_Kp", &num_match_Kp);
    tree->Branch("num_match_Km", &num_match_Km);
    tree->Branch("num_match_pi", &num_match_pi);
    tree->Branch("num_match_mu", &num_match_mu);
    tree->Branch("num_tight_match_Kp", &num_tight_match_Kp);
    tree->Branch("num_tight_match_Km", &num_tight_match_Km);
    tree->Branch("num_tight_match_pi", &num_tight_match_pi);
    tree->Branch("match_dR_Kp", &match_dR_Kp);
    tree->Branch("match_dR_Km", &match_dR_Km);
    tree->Branch("match_dR_pi", &match_dR_pi);
    tree->Branch("match_dR_mu", &match_dR_mu);

    tree->Branch("match_Kp_isIsolatedChargedHadron", &match_Kp_isIsolatedChargedHadron_vec);
    tree->Branch("match_Kp_charge", &match_Kp_charge_vec);
    tree->Branch("match_Kp_eta", &match_Kp_eta_vec);
    tree->Branch("match_Kp_phi", &match_Kp_phi_vec);
    tree->Branch("match_Kp_vx", &match_Kp_vx_vec);
    tree->Branch("match_Kp_vy", &match_Kp_vy_vec);
    tree->Branch("match_Kp_vz", &match_Kp_vz_vec);
    tree->Branch("match_Kp_p", &match_Kp_p_vec);
    tree->Branch("match_Kp_pt", &match_Kp_pt_vec);
    tree->Branch("match_Kp_px", &match_Kp_px_vec);
    tree->Branch("match_Kp_py", &match_Kp_py_vec);
    tree->Branch("match_Kp_pz", &match_Kp_pz_vec);

    tree->Branch("match_Km_isIsolatedChargedHadron", &match_Km_isIsolatedChargedHadron_vec);
    tree->Branch("match_Km_charge", &match_Km_charge_vec);
    tree->Branch("match_Km_eta", &match_Km_eta_vec);
    tree->Branch("match_Km_phi", &match_Km_phi_vec);
    tree->Branch("match_Km_vx", &match_Km_vx_vec);
    tree->Branch("match_Km_vy", &match_Km_vy_vec);
    tree->Branch("match_Km_vz", &match_Km_vz_vec);
    tree->Branch("match_Km_p", &match_Km_p_vec);
    tree->Branch("match_Km_pt", &match_Km_pt_vec);
    tree->Branch("match_Km_px", &match_Km_px_vec);
    tree->Branch("match_Km_py", &match_Km_py_vec);
    tree->Branch("match_Km_pz", &match_Km_pz_vec);

    tree->Branch("match_pi_isIsolatedChargedHadron", &match_pi_isIsolatedChargedHadron_vec);
    tree->Branch("match_pi_charge", &match_pi_charge_vec);
    tree->Branch("match_pi_eta", &match_pi_eta_vec);
    tree->Branch("match_pi_phi", &match_pi_phi_vec);
    tree->Branch("match_pi_vx", &match_pi_vx_vec);
    tree->Branch("match_pi_vy", &match_pi_vy_vec);
    tree->Branch("match_pi_vz", &match_pi_vz_vec);
    tree->Branch("match_pi_p", &match_pi_p_vec);
    tree->Branch("match_pi_pt", &match_pi_pt_vec);
    tree->Branch("match_pi_px", &match_pi_px_vec);
    tree->Branch("match_pi_py", &match_pi_py_vec);
    tree->Branch("match_pi_pz", &match_pi_pz_vec);

    tree->Branch("match_phi_eta", &match_phi_eta_vec);
    tree->Branch("match_phi_phi", &match_phi_phi_vec);
    tree->Branch("match_phi_p", &match_phi_p_vec);
    tree->Branch("match_phi_pt", &match_phi_pt_vec);
    tree->Branch("match_phi_px", &match_phi_px_vec);
    tree->Branch("match_phi_py", &match_phi_py_vec);
    tree->Branch("match_phi_pz", &match_phi_pz_vec);
    tree->Branch("match_phi_invm", &match_phi_invm_vec);

    tree->Branch("match_Ds_eta", &match_Ds_eta_vec);
    tree->Branch("match_Ds_phi", &match_Ds_phi_vec);
    tree->Branch("match_Ds_p", &match_Ds_p_vec);
    tree->Branch("match_Ds_pt", &match_Ds_pt_vec);
    tree->Branch("match_Ds_px", &match_Ds_px_vec);
    tree->Branch("match_Ds_py", &match_Ds_py_vec);
    tree->Branch("match_Ds_pz", &match_Ds_pz_vec);
    tree->Branch("match_Ds_invm", &match_Ds_invm_vec);

    tree->Branch("match_dR_Kp_Km", &match_dR_Kp_Km_vec);
    tree->Branch("match_dR_Kp_phi", &match_dR_Kp_phi_vec);
    tree->Branch("match_dR_Km_phi", &match_dR_Km_phi_vec);
    tree->Branch("match_dR_Kp_pi", &match_dR_Kp_pi_vec);
    tree->Branch("match_dR_Km_pi", &match_dR_Km_pi_vec);
    tree->Branch("match_dR_pi_phi", &match_dR_pi_phi_vec);
    tree->Branch("match_dR_Kp_Ds", &match_dR_Kp_Ds_vec);
    tree->Branch("match_dR_Km_Ds", &match_dR_Km_Ds_vec);
    tree->Branch("match_dR_phi_Ds", &match_dR_phi_Ds_vec);
    tree->Branch("match_dR_pi_Ds", &match_dR_pi_Ds_vec);

    tree->Branch("match_dxy_phi_Ds", &match_dxy_phi_Ds_vec);
    tree->Branch("match_dz_phi_Ds", &match_dz_phi_Ds_vec);

    // Fit on phi
    tree->Branch("match_phiFit_chi2", &match_phiFit_chi2_vec);
    tree->Branch("match_phiFit_ndof", &match_phiFit_ndof_vec);
    tree->Branch("match_phiFit_chi2ndof", &match_phiFit_chi2ndof_vec);
    tree->Branch("match_phiFit_vx", &match_phiFit_vx_vec);
    tree->Branch("match_phiFit_vy", &match_phiFit_vy_vec);
    tree->Branch("match_phiFit_vz", &match_phiFit_vz_vec);
    tree->Branch("match_phiFit_vxerr", &match_phiFit_vxerr_vec);
    tree->Branch("match_phiFit_vyerr", &match_phiFit_vyerr_vec);
    tree->Branch("match_phiFit_vzerr", &match_phiFit_vzerr_vec);

    tree->Branch("match_phiFit_Kp_eta", &match_phiFit_Kp_eta_vec);
    tree->Branch("match_phiFit_Kp_phi", &match_phiFit_Kp_phi_vec);
    tree->Branch("match_phiFit_Kp_p", &match_phiFit_Kp_p_vec);
    tree->Branch("match_phiFit_Kp_pt", &match_phiFit_Kp_pt_vec);
    tree->Branch("match_phiFit_Kp_px", &match_phiFit_Kp_px_vec);
    tree->Branch("match_phiFit_Kp_py", &match_phiFit_Kp_py_vec);
    tree->Branch("match_phiFit_Kp_pz", &match_phiFit_Kp_pz_vec);

    tree->Branch("match_phiFit_Km_eta", &match_phiFit_Km_eta_vec);
    tree->Branch("match_phiFit_Km_phi", &match_phiFit_Km_phi_vec);
    tree->Branch("match_phiFit_Km_p", &match_phiFit_Km_p_vec);
    tree->Branch("match_phiFit_Km_pt", &match_phiFit_Km_pt_vec);
    tree->Branch("match_phiFit_Km_px", &match_phiFit_Km_px_vec);
    tree->Branch("match_phiFit_Km_py", &match_phiFit_Km_py_vec);
    tree->Branch("match_phiFit_Km_pz", &match_phiFit_Km_pz_vec);

    tree->Branch("match_phiFit_pi_eta", &match_phiFit_pi_eta_vec);
    tree->Branch("match_phiFit_pi_phi", &match_phiFit_pi_phi_vec);
    tree->Branch("match_phiFit_pi_p", &match_phiFit_pi_p_vec);
    tree->Branch("match_phiFit_pi_pt", &match_phiFit_pi_pt_vec);
    tree->Branch("match_phiFit_pi_px", &match_phiFit_pi_px_vec);
    tree->Branch("match_phiFit_pi_py", &match_phiFit_pi_py_vec);
    tree->Branch("match_phiFit_pi_pz", &match_phiFit_pi_pz_vec);

    tree->Branch("match_phiFit_phi_eta", &match_phiFit_phi_eta_vec);
    tree->Branch("match_phiFit_phi_phi", &match_phiFit_phi_phi_vec);
    tree->Branch("match_phiFit_phi_p", &match_phiFit_phi_p_vec);
    tree->Branch("match_phiFit_phi_pt", &match_phiFit_phi_pt_vec);
    tree->Branch("match_phiFit_phi_px", &match_phiFit_phi_px_vec);
    tree->Branch("match_phiFit_phi_py", &match_phiFit_phi_py_vec);
    tree->Branch("match_phiFit_phi_pz", &match_phiFit_phi_pz_vec);
    tree->Branch("match_phiFit_phi_invm", &match_phiFit_phi_invm_vec);

    tree->Branch("match_phiFit_Ds_eta", &match_phiFit_Ds_eta_vec);
    tree->Branch("match_phiFit_Ds_phi", &match_phiFit_Ds_phi_vec);
    tree->Branch("match_phiFit_Ds_p", &match_phiFit_Ds_p_vec);
    tree->Branch("match_phiFit_Ds_pt", &match_phiFit_Ds_pt_vec);
    tree->Branch("match_phiFit_Ds_px", &match_phiFit_Ds_px_vec);
    tree->Branch("match_phiFit_Ds_py", &match_phiFit_Ds_py_vec);
    tree->Branch("match_phiFit_Ds_pz", &match_phiFit_Ds_pz_vec);
    tree->Branch("match_phiFit_Ds_invm", &match_phiFit_Ds_invm_vec);

    tree->Branch("match_phiFit_dR_Kp_Km", &match_phiFit_dR_Kp_Km_vec);
    tree->Branch("match_phiFit_dR_Kp_phi", &match_phiFit_dR_Kp_phi_vec);
    tree->Branch("match_phiFit_dR_Km_phi", &match_phiFit_dR_Km_phi_vec);
    tree->Branch("match_phiFit_dR_Kp_pi", &match_phiFit_dR_Kp_pi_vec);
    tree->Branch("match_phiFit_dR_Km_pi", &match_phiFit_dR_Km_pi_vec);
    tree->Branch("match_phiFit_dR_pi_phi", &match_phiFit_dR_pi_phi_vec);
    tree->Branch("match_phiFit_dR_Kp_Ds", &match_phiFit_dR_Kp_Ds_vec);
    tree->Branch("match_phiFit_dR_Km_Ds", &match_phiFit_dR_Km_Ds_vec);
    tree->Branch("match_phiFit_dR_phi_Ds", &match_phiFit_dR_phi_Ds_vec);
    tree->Branch("match_phiFit_dR_pi_Ds", &match_phiFit_dR_pi_Ds_vec);

    // Fit on Ds 
    tree->Branch("match_DsFit_chi2", &match_DsFit_chi2_vec);
    tree->Branch("match_DsFit_ndof", &match_DsFit_ndof_vec);
    tree->Branch("match_DsFit_chi2ndof", &match_DsFit_chi2ndof_vec);
    tree->Branch("match_DsFit_vx", &match_DsFit_vx_vec);
    tree->Branch("match_DsFit_vy", &match_DsFit_vy_vec);
    tree->Branch("match_DsFit_vz", &match_DsFit_vz_vec);
    tree->Branch("match_DsFit_vxerr", &match_DsFit_vxerr_vec);
    tree->Branch("match_DsFit_vyerr", &match_DsFit_vyerr_vec);
    tree->Branch("match_DsFit_vzerr", &match_DsFit_vzerr_vec);

    tree->Branch("match_DsFit_Kp_eta", &match_DsFit_Kp_eta_vec);
    tree->Branch("match_DsFit_Kp_phi", &match_DsFit_Kp_phi_vec);
    tree->Branch("match_DsFit_Kp_p", &match_DsFit_Kp_p_vec);
    tree->Branch("match_DsFit_Kp_pt", &match_DsFit_Kp_pt_vec);
    tree->Branch("match_DsFit_Kp_px", &match_DsFit_Kp_px_vec);
    tree->Branch("match_DsFit_Kp_py", &match_DsFit_Kp_py_vec);
    tree->Branch("match_DsFit_Kp_pz", &match_DsFit_Kp_pz_vec);

    tree->Branch("match_DsFit_Km_eta", &match_DsFit_Km_eta_vec);
    tree->Branch("match_DsFit_Km_phi", &match_DsFit_Km_phi_vec);
    tree->Branch("match_DsFit_Km_p", &match_DsFit_Km_p_vec);
    tree->Branch("match_DsFit_Km_pt", &match_DsFit_Km_pt_vec);
    tree->Branch("match_DsFit_Km_px", &match_DsFit_Km_px_vec);
    tree->Branch("match_DsFit_Km_py", &match_DsFit_Km_py_vec);
    tree->Branch("match_DsFit_Km_pz", &match_DsFit_Km_pz_vec);

    tree->Branch("match_DsFit_pi_eta", &match_DsFit_pi_eta_vec);
    tree->Branch("match_DsFit_pi_phi", &match_DsFit_pi_phi_vec);
    tree->Branch("match_DsFit_pi_p", &match_DsFit_pi_p_vec);
    tree->Branch("match_DsFit_pi_pt", &match_DsFit_pi_pt_vec);
    tree->Branch("match_DsFit_pi_px", &match_DsFit_pi_px_vec);
    tree->Branch("match_DsFit_pi_py", &match_DsFit_pi_py_vec);
    tree->Branch("match_DsFit_pi_pz", &match_DsFit_pi_pz_vec);

    tree->Branch("match_DsFit_phi_eta", &match_DsFit_phi_eta_vec);
    tree->Branch("match_DsFit_phi_phi", &match_DsFit_phi_phi_vec);
    tree->Branch("match_DsFit_phi_p", &match_DsFit_phi_p_vec);
    tree->Branch("match_DsFit_phi_pt", &match_DsFit_phi_pt_vec);
    tree->Branch("match_DsFit_phi_px", &match_DsFit_phi_px_vec);
    tree->Branch("match_DsFit_phi_py", &match_DsFit_phi_py_vec);
    tree->Branch("match_DsFit_phi_pz", &match_DsFit_phi_pz_vec);
    tree->Branch("match_DsFit_phi_invm", &match_DsFit_phi_invm_vec);

    tree->Branch("match_DsFit_Ds_eta", &match_DsFit_Ds_eta_vec);
    tree->Branch("match_DsFit_Ds_phi", &match_DsFit_Ds_phi_vec);
    tree->Branch("match_DsFit_Ds_p", &match_DsFit_Ds_p_vec);
    tree->Branch("match_DsFit_Ds_pt", &match_DsFit_Ds_pt_vec);
    tree->Branch("match_DsFit_Ds_px", &match_DsFit_Ds_px_vec);
    tree->Branch("match_DsFit_Ds_py", &match_DsFit_Ds_py_vec);
    tree->Branch("match_DsFit_Ds_pz", &match_DsFit_Ds_pz_vec);
    tree->Branch("match_DsFit_Ds_invm", &match_DsFit_Ds_invm_vec);

    tree->Branch("match_DsFit_dR_Kp_Km", &match_DsFit_dR_Kp_Km_vec);
    tree->Branch("match_DsFit_dR_Kp_phi", &match_DsFit_dR_Kp_phi_vec);
    tree->Branch("match_DsFit_dR_Km_phi", &match_DsFit_dR_Km_phi_vec);
    tree->Branch("match_DsFit_dR_Kp_pi", &match_DsFit_dR_Kp_pi_vec);
    tree->Branch("match_DsFit_dR_Km_pi", &match_DsFit_dR_Km_pi_vec);
    tree->Branch("match_DsFit_dR_pi_phi", &match_DsFit_dR_pi_phi_vec);
    tree->Branch("match_DsFit_dR_Kp_Ds", &match_DsFit_dR_Kp_Ds_vec);
    tree->Branch("match_DsFit_dR_Km_Ds", &match_DsFit_dR_Km_Ds_vec);
    tree->Branch("match_DsFit_dR_phi_Ds", &match_DsFit_dR_phi_Ds_vec);
    tree->Branch("match_DsFit_dR_pi_Ds", &match_DsFit_dR_pi_Ds_vec);

    tree->Branch("match_DsFit_Mconstraint_Ds_invm", &match_DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("match_Ds_primvtx_FDxy", &match_Ds_primvtx_FDxy_vec);
    tree->Branch("match_Ds_primvtx_FDz", &match_Ds_primvtx_FDz_vec);
    tree->Branch("match_Ds_primvtx_FD", &match_Ds_primvtx_FD_vec);
    tree->Branch("match_Ds_primvtx_FDxyerr", &match_Ds_primvtx_FDxyerr_vec);
    tree->Branch("match_Ds_primvtx_FDxychi2", &match_Ds_primvtx_FDxychi2_vec);
    tree->Branch("match_Ds_primvtx_FDzerr", &match_Ds_primvtx_FDzerr_vec);
    tree->Branch("match_Ds_primvtx_FDzchi2", &match_Ds_primvtx_FDzchi2_vec);
    tree->Branch("match_Ds_primvtx_FDerr", &match_Ds_primvtx_FDerr_vec);
    tree->Branch("match_Ds_primvtx_FDchi2", &match_Ds_primvtx_FDchi2_vec);
    tree->Branch("match_Ds_primvtx_dira", &match_Ds_primvtx_dira_vec);
    tree->Branch("match_Ds_primvtx_dira_angle", &match_Ds_primvtx_dira_angle_vec);
    tree->Branch("match_Kp_primvtx_ip", &match_Kp_primvtx_ip_vec);
    tree->Branch("match_Kp_primvtx_iperr", &match_Kp_primvtx_iperr_vec);
    tree->Branch("match_Kp_primvtx_ipchi2", &match_Kp_primvtx_ipchi2_vec);
    tree->Branch("match_Km_primvtx_ip", &match_Km_primvtx_ip_vec);
    tree->Branch("match_Km_primvtx_iperr", &match_Km_primvtx_iperr_vec);
    tree->Branch("match_Km_primvtx_ipchi2", &match_Km_primvtx_ipchi2_vec);
    tree->Branch("match_pi_primvtx_ip", &match_pi_primvtx_ip_vec);
    tree->Branch("match_pi_primvtx_iperr", &match_pi_primvtx_iperr_vec);
    tree->Branch("match_pi_primvtx_ipchi2", &match_pi_primvtx_ipchi2_vec);
    tree->Branch("match_phi_primvtx_ip", &match_phi_primvtx_ip_vec);
    tree->Branch("match_phi_primvtx_iperr", &match_phi_primvtx_iperr_vec);
    tree->Branch("match_phi_primvtx_ipchi2", &match_phi_primvtx_ipchi2_vec);
    tree->Branch("match_Ds_primvtx_ip", &match_Ds_primvtx_ip_vec);
    tree->Branch("match_Ds_primvtx_iperr", &match_Ds_primvtx_iperr_vec);
    tree->Branch("match_Ds_primvtx_ipchi2", &match_Ds_primvtx_ipchi2_vec);

    tree->Branch("match_Ds_IsoR03_sumChargedHadronPt", &match_Ds_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("match_Ds_IsoR03_sumNeutralHadronPt", &match_Ds_IsoR03_sumNeutralHadronPt_vec);
    tree->Branch("match_Ds_IsoR03_sumPhotonPt", &match_Ds_IsoR03_sumPhotonPt_vec);
    tree->Branch("match_Ds_IsoR03_sumPUPt", &match_Ds_IsoR03_sumPUPt_vec);
    tree->Branch("match_Ds_PFIsoR03", &match_Ds_PFIsoR03_vec);
    tree->Branch("match_Ds_IsoR04_sumChargedHadronPt", &match_Ds_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("match_Ds_IsoR04_sumNeutralHadronPt", &match_Ds_IsoR04_sumNeutralHadronPt_vec);
    tree->Branch("match_Ds_IsoR04_sumPhotonPt", &match_Ds_IsoR04_sumPhotonPt_vec);
    tree->Branch("match_Ds_IsoR04_sumPUPt", &match_Ds_IsoR04_sumPUPt_vec);
    tree->Branch("match_Ds_PFIsoR04", &match_Ds_PFIsoR04_vec);

    tree->Branch("match_PV_noDs_IsValid", &match_PV_noDs_IsValid_vec);
    tree->Branch("match_PV_noDs_IsFake", &match_PV_noDs_IsFake_vec);
    tree->Branch("match_PV_noDs_chi2", &match_PV_noDs_chi2_vec);
    tree->Branch("match_PV_noDs_ndof", &match_PV_noDs_ndof_vec);
    tree->Branch("match_PV_noDs_chi2ndof", &match_PV_noDs_chi2ndof_vec);
    tree->Branch("match_PV_noDs_x", &match_PV_noDs_x_vec);
    tree->Branch("match_PV_noDs_y", &match_PV_noDs_y_vec);
    tree->Branch("match_PV_noDs_z", &match_PV_noDs_z_vec);
    tree->Branch("match_PV_noDs_xerr", &match_PV_noDs_xerr_vec);
    tree->Branch("match_PV_noDs_yerr", &match_PV_noDs_yerr_vec);
    tree->Branch("match_PV_noDs_zerr", &match_PV_noDs_zerr_vec);

    tree->Branch("match_Ds_PVnoDs_FDxy", &match_Ds_PVnoDs_FDxy_vec);
    tree->Branch("match_Ds_PVnoDs_FDz", &match_Ds_PVnoDs_FDz_vec);
    tree->Branch("match_Ds_PVnoDs_FD", &match_Ds_PVnoDs_FD_vec);
    tree->Branch("match_Ds_PVnoDs_FDxyerr", &match_Ds_PVnoDs_FDxyerr_vec);
    tree->Branch("match_Ds_PVnoDs_FDxychi2", &match_Ds_PVnoDs_FDxychi2_vec);
    tree->Branch("match_Ds_PVnoDs_FDzerr", &match_Ds_PVnoDs_FDzerr_vec);
    tree->Branch("match_Ds_PVnoDs_FDzchi2", &match_Ds_PVnoDs_FDzchi2_vec);
    tree->Branch("match_Ds_PVnoDs_FDerr", &match_Ds_PVnoDs_FDerr_vec);
    tree->Branch("match_Ds_PVnoDs_FDchi2", &match_Ds_PVnoDs_FDchi2_vec);
    tree->Branch("match_Ds_PVnoDs_dira", &match_Ds_PVnoDs_dira_vec);
    tree->Branch("match_Ds_PVnoDs_dira_angle", &match_Ds_PVnoDs_dira_angle_vec);
    tree->Branch("match_Kp_PVnoDs_ip", &match_Kp_PVnoDs_ip_vec);
    tree->Branch("match_Kp_PVnoDs_iperr", &match_Kp_PVnoDs_iperr_vec);
    tree->Branch("match_Kp_PVnoDs_ipchi2", &match_Kp_PVnoDs_ipchi2_vec);
    tree->Branch("match_Km_PVnoDs_ip", &match_Km_PVnoDs_ip_vec);
    tree->Branch("match_Km_PVnoDs_iperr", &match_Km_PVnoDs_iperr_vec);
    tree->Branch("match_Km_PVnoDs_ipchi2", &match_Km_PVnoDs_ipchi2_vec);
    tree->Branch("match_pi_PVnoDs_ip", &match_pi_PVnoDs_ip_vec);
    tree->Branch("match_pi_PVnoDs_iperr", &match_pi_PVnoDs_iperr_vec);
    tree->Branch("match_pi_PVnoDs_ipchi2", &match_pi_PVnoDs_ipchi2_vec);
    tree->Branch("match_phi_PVnoDs_ip", &match_phi_PVnoDs_ip_vec);
    tree->Branch("match_phi_PVnoDs_iperr", &match_phi_PVnoDs_iperr_vec);
    tree->Branch("match_phi_PVnoDs_ipchi2", &match_phi_PVnoDs_ipchi2_vec);
    tree->Branch("match_Ds_PVnoDs_ip", &match_Ds_PVnoDs_ip_vec);
    tree->Branch("match_Ds_PVnoDs_iperr", &match_Ds_PVnoDs_iperr_vec);
    tree->Branch("match_Ds_PVnoDs_ipchi2", &match_Ds_PVnoDs_ipchi2_vec);

    tree->Branch("match_PV_withDs_IsValid", &match_PV_withDs_IsValid_vec);
    tree->Branch("match_PV_withDs_IsFake", &match_PV_withDs_IsFake_vec);
    tree->Branch("match_PV_withDs_chi2", &match_PV_withDs_chi2_vec);
    tree->Branch("match_PV_withDs_ndof", &match_PV_withDs_ndof_vec);
    tree->Branch("match_PV_withDs_chi2ndof", &match_PV_withDs_chi2ndof_vec);
    tree->Branch("match_PV_withDs_x", &match_PV_withDs_x_vec);
    tree->Branch("match_PV_withDs_y", &match_PV_withDs_y_vec);
    tree->Branch("match_PV_withDs_z", &match_PV_withDs_z_vec);
    tree->Branch("match_PV_withDs_xerr", &match_PV_withDs_xerr_vec);
    tree->Branch("match_PV_withDs_yerr", &match_PV_withDs_yerr_vec);
    tree->Branch("match_PV_withDs_zerr", &match_PV_withDs_zerr_vec);

    tree->Branch("match_Ds_PVwithDs_FDxy", &match_Ds_PVwithDs_FDxy_vec);
    tree->Branch("match_Ds_PVwithDs_FDz", &match_Ds_PVwithDs_FDz_vec);
    tree->Branch("match_Ds_PVwithDs_FD", &match_Ds_PVwithDs_FD_vec);
    tree->Branch("match_Ds_PVwithDs_FDxyerr", &match_Ds_PVwithDs_FDxyerr_vec);
    tree->Branch("match_Ds_PVwithDs_FDxychi2", &match_Ds_PVwithDs_FDxychi2_vec);
    tree->Branch("match_Ds_PVwithDs_FDzerr", &match_Ds_PVwithDs_FDzerr_vec);
    tree->Branch("match_Ds_PVwithDs_FDzchi2", &match_Ds_PVwithDs_FDzchi2_vec);
    tree->Branch("match_Ds_PVwithDs_FDerr", &match_Ds_PVwithDs_FDerr_vec);
    tree->Branch("match_Ds_PVwithDs_FDchi2", &match_Ds_PVwithDs_FDchi2_vec);
    tree->Branch("match_Ds_PVwithDs_dira", &match_Ds_PVwithDs_dira_vec);
    tree->Branch("match_Ds_PVwithDs_dira_angle", &match_Ds_PVwithDs_dira_angle_vec);
    tree->Branch("match_Kp_PVwithDs_ip", &match_Kp_PVwithDs_ip_vec);
    tree->Branch("match_Kp_PVwithDs_iperr", &match_Kp_PVwithDs_iperr_vec);
    tree->Branch("match_Kp_PVwithDs_ipchi2", &match_Kp_PVwithDs_ipchi2_vec);
    tree->Branch("match_Km_PVwithDs_ip", &match_Km_PVwithDs_ip_vec);
    tree->Branch("match_Km_PVwithDs_iperr", &match_Km_PVwithDs_iperr_vec);
    tree->Branch("match_Km_PVwithDs_ipchi2", &match_Km_PVwithDs_ipchi2_vec);
    tree->Branch("match_pi_PVwithDs_ip", &match_pi_PVwithDs_ip_vec);
    tree->Branch("match_pi_PVwithDs_iperr", &match_pi_PVwithDs_iperr_vec);
    tree->Branch("match_pi_PVwithDs_ipchi2", &match_pi_PVwithDs_ipchi2_vec);
    tree->Branch("match_phi_PVwithDs_ip", &match_phi_PVwithDs_ip_vec);
    tree->Branch("match_phi_PVwithDs_iperr", &match_phi_PVwithDs_iperr_vec);
    tree->Branch("match_phi_PVwithDs_ipchi2", &match_phi_PVwithDs_ipchi2_vec);
    tree->Branch("match_Ds_PVwithDs_ip", &match_Ds_PVwithDs_ip_vec);
    tree->Branch("match_Ds_PVwithDs_iperr", &match_Ds_PVwithDs_iperr_vec);
    tree->Branch("match_Ds_PVwithDs_ipchi2", &match_Ds_PVwithDs_ipchi2_vec);

    tree->Branch("match_mu_charge", &match_mu_charge_vec);
    tree->Branch("match_mu_eta", &match_mu_eta_vec);
    tree->Branch("match_mu_phi", &match_mu_phi_vec);
    tree->Branch("match_mu_vx", &match_mu_vx_vec);
    tree->Branch("match_mu_vy", &match_mu_vy_vec);
    tree->Branch("match_mu_vz", &match_mu_vz_vec);
    tree->Branch("match_mu_p", &match_mu_p_vec);
    tree->Branch("match_mu_pt", &match_mu_pt_vec);
    tree->Branch("match_mu_px", &match_mu_px_vec);
    tree->Branch("match_mu_py", &match_mu_py_vec);
    tree->Branch("match_mu_pz", &match_mu_pz_vec);
    tree->Branch("match_mu_isHighPt", &match_mu_isHighPt_vec);
    tree->Branch("match_mu_isLoose", &match_mu_isLoose_vec);
    tree->Branch("match_mu_isMedium", &match_mu_isMedium_vec);
    tree->Branch("match_mu_isSoft", &match_mu_isSoft_vec);
    tree->Branch("match_mu_isTight", &match_mu_isTight_vec);
    tree->Branch("match_mu_isPF", &match_mu_isPF_vec);
    tree->Branch("match_mu_isTracker", &match_mu_isTracker_vec);
    tree->Branch("match_mu_isGlobal", &match_mu_isGlobal_vec);
    tree->Branch("match_mu_IsoR03_sumChargedHadronPt", &match_mu_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR03_sumChargedParticlePt", &match_mu_IsoR03_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR03_sumNeutralHadronEt", &match_mu_IsoR03_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR03_sumPhotonEt", &match_mu_IsoR03_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR03_sumPUPt", &match_mu_IsoR03_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR03", &match_mu_PFIsoR03_vec);
    tree->Branch("match_mu_IsoR04_sumChargedHadronPt", &match_mu_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR04_sumChargedParticlePt", &match_mu_IsoR04_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR04_sumNeutralHadronEt", &match_mu_IsoR04_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR04_sumPhotonEt", &match_mu_IsoR04_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR04_sumPUPt", &match_mu_IsoR04_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR04", &match_mu_PFIsoR04_vec);
    tree->Branch("match_mu_primvtx_dxy", &match_mu_primvtx_dxy_vec);
    tree->Branch("match_mu_primvtx_dxyerr", &match_mu_primvtx_dxyerr_vec);
    tree->Branch("match_mu_primvtx_dz", &match_mu_primvtx_dz_vec);
    tree->Branch("match_mu_primvtx_dzerr", &match_mu_primvtx_dzerr_vec);
    tree->Branch("match_mu_primvtx_ip", &match_mu_primvtx_ip_vec);
    tree->Branch("match_mu_primvtx_iperr", &match_mu_primvtx_iperr_vec);
    tree->Branch("match_mu_primvtx_ipchi2", &match_mu_primvtx_ipchi2_vec);
    
    tree->Branch("match_mu_charge", &match_mu_charge_vec);
    tree->Branch("match_mu_eta", &match_mu_eta_vec);
    tree->Branch("match_mu_phi", &match_mu_phi_vec);
    tree->Branch("match_mu_vx", &match_mu_vx_vec);
    tree->Branch("match_mu_vy", &match_mu_vy_vec);
    tree->Branch("match_mu_vz", &match_mu_vz_vec);
    tree->Branch("match_mu_p", &match_mu_p_vec);
    tree->Branch("match_mu_pt", &match_mu_pt_vec);
    tree->Branch("match_mu_px", &match_mu_px_vec);
    tree->Branch("match_mu_py", &match_mu_py_vec);
    tree->Branch("match_mu_pz", &match_mu_pz_vec);
    tree->Branch("match_mu_isHighPt", &match_mu_isHighPt_vec);
    tree->Branch("match_mu_isLoose", &match_mu_isLoose_vec);
    tree->Branch("match_mu_isMedium", &match_mu_isMedium_vec);
    tree->Branch("match_mu_isSoft", &match_mu_isSoft_vec);
    tree->Branch("match_mu_isTight", &match_mu_isTight_vec);
    tree->Branch("match_mu_isPF", &match_mu_isPF_vec);
    tree->Branch("match_mu_isTracker", &match_mu_isTracker_vec);
    tree->Branch("match_mu_isGlobal", &match_mu_isGlobal_vec);
    tree->Branch("match_mu_IsoR03_sumChargedHadronPt", &match_mu_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR03_sumChargedParticlePt", &match_mu_IsoR03_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR03_sumNeutralHadronEt", &match_mu_IsoR03_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR03_sumPhotonEt", &match_mu_IsoR03_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR03_sumPUPt", &match_mu_IsoR03_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR03", &match_mu_PFIsoR03_vec);
    tree->Branch("match_mu_IsoR04_sumChargedHadronPt", &match_mu_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR04_sumChargedParticlePt", &match_mu_IsoR04_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR04_sumNeutralHadronEt", &match_mu_IsoR04_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR04_sumPhotonEt", &match_mu_IsoR04_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR04_sumPUPt", &match_mu_IsoR04_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR04", &match_mu_PFIsoR04_vec);
    tree->Branch("match_mu_primvtx_dxy", &match_mu_primvtx_dxy_vec);
    tree->Branch("match_mu_primvtx_dxyerr", &match_mu_primvtx_dxyerr_vec);
    tree->Branch("match_mu_primvtx_dz", &match_mu_primvtx_dz_vec);
    tree->Branch("match_mu_primvtx_dzerr", &match_mu_primvtx_dzerr_vec);
    tree->Branch("match_mu_primvtx_ip", &match_mu_primvtx_ip_vec);
    tree->Branch("match_mu_primvtx_iperr", &match_mu_primvtx_iperr_vec);
    tree->Branch("match_mu_primvtx_ipchi2", &match_mu_primvtx_ipchi2_vec);

    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

    tree->Branch("Kp_isIsolatedChargedHadron", &Kp_isIsolatedChargedHadron_vec);
    tree->Branch("Kp_charge", &Kp_charge_vec);
    tree->Branch("Kp_eta", &Kp_eta_vec);
    tree->Branch("Kp_phi", &Kp_phi_vec);
    tree->Branch("Kp_vx", &Kp_vx_vec);
    tree->Branch("Kp_vy", &Kp_vy_vec);
    tree->Branch("Kp_vz", &Kp_vz_vec);
    tree->Branch("Kp_p", &Kp_p_vec);
    tree->Branch("Kp_pt", &Kp_pt_vec);
    tree->Branch("Kp_px", &Kp_px_vec);
    tree->Branch("Kp_py", &Kp_py_vec);
    tree->Branch("Kp_pz", &Kp_pz_vec);

    tree->Branch("Km_isIsolatedChargedHadron", &Km_isIsolatedChargedHadron_vec);
    tree->Branch("Km_charge", &Km_charge_vec);
    tree->Branch("Km_eta", &Km_eta_vec);
    tree->Branch("Km_phi", &Km_phi_vec);
    tree->Branch("Km_vx", &Km_vx_vec);
    tree->Branch("Km_vy", &Km_vy_vec);
    tree->Branch("Km_vz", &Km_vz_vec);
    tree->Branch("Km_p", &Km_p_vec);
    tree->Branch("Km_pt", &Km_pt_vec);
    tree->Branch("Km_px", &Km_px_vec);
    tree->Branch("Km_py", &Km_py_vec);
    tree->Branch("Km_pz", &Km_pz_vec);

    tree->Branch("pi_isIsolatedChargedHadron", &pi_isIsolatedChargedHadron_vec);
    tree->Branch("pi_charge", &pi_charge_vec);
    tree->Branch("pi_eta", &pi_eta_vec);
    tree->Branch("pi_phi", &pi_phi_vec);
    tree->Branch("pi_vx", &pi_vx_vec);
    tree->Branch("pi_vy", &pi_vy_vec);
    tree->Branch("pi_vz", &pi_vz_vec);
    tree->Branch("pi_p", &pi_p_vec);
    tree->Branch("pi_pt", &pi_pt_vec);
    tree->Branch("pi_px", &pi_px_vec);
    tree->Branch("pi_py", &pi_py_vec);
    tree->Branch("pi_pz", &pi_pz_vec);

    tree->Branch("phi_eta", &phi_eta_vec);
    tree->Branch("phi_phi", &phi_phi_vec);
    tree->Branch("phi_p", &phi_p_vec);
    tree->Branch("phi_pt", &phi_pt_vec);
    tree->Branch("phi_px", &phi_px_vec);
    tree->Branch("phi_py", &phi_py_vec);
    tree->Branch("phi_pz", &phi_pz_vec);
    tree->Branch("phi_invm", &phi_invm_vec);

    tree->Branch("Ds_eta", &Ds_eta_vec);
    tree->Branch("Ds_phi", &Ds_phi_vec);
    tree->Branch("Ds_p", &Ds_p_vec);
    tree->Branch("Ds_pt", &Ds_pt_vec);
    tree->Branch("Ds_px", &Ds_px_vec);
    tree->Branch("Ds_py", &Ds_py_vec);
    tree->Branch("Ds_pz", &Ds_pz_vec);
    tree->Branch("Ds_invm", &Ds_invm_vec);

    tree->Branch("dR_Kp_Km", &dR_Kp_Km_vec);
    tree->Branch("dR_Kp_phi", &dR_Kp_phi_vec);
    tree->Branch("dR_Km_phi", &dR_Km_phi_vec);
    tree->Branch("dR_Kp_pi", &dR_Kp_pi_vec);
    tree->Branch("dR_Km_pi", &dR_Km_pi_vec);
    tree->Branch("dR_pi_phi", &dR_pi_phi_vec);
    tree->Branch("dR_Kp_Ds", &dR_Kp_Ds_vec);
    tree->Branch("dR_Km_Ds", &dR_Km_Ds_vec);
    tree->Branch("dR_phi_Ds", &dR_phi_Ds_vec);
    tree->Branch("dR_pi_Ds", &dR_pi_Ds_vec);

    tree->Branch("dxy_phi_Ds", &dxy_phi_Ds_vec);
    tree->Branch("dz_phi_Ds", &dz_phi_Ds_vec);

    // Fit on phi
    tree->Branch("phiFit_chi2", &phiFit_chi2_vec);
    tree->Branch("phiFit_ndof", &phiFit_ndof_vec);
    tree->Branch("phiFit_chi2ndof", &phiFit_chi2ndof_vec);
    tree->Branch("phiFit_vx", &phiFit_vx_vec);
    tree->Branch("phiFit_vy", &phiFit_vy_vec);
    tree->Branch("phiFit_vz", &phiFit_vz_vec);
    tree->Branch("phiFit_vxerr", &phiFit_vxerr_vec);
    tree->Branch("phiFit_vyerr", &phiFit_vyerr_vec);
    tree->Branch("phiFit_vzerr", &phiFit_vzerr_vec);

    tree->Branch("phiFit_Kp_eta", &phiFit_Kp_eta_vec);
    tree->Branch("phiFit_Kp_phi", &phiFit_Kp_phi_vec);
    tree->Branch("phiFit_Kp_p", &phiFit_Kp_p_vec);
    tree->Branch("phiFit_Kp_pt", &phiFit_Kp_pt_vec);
    tree->Branch("phiFit_Kp_px", &phiFit_Kp_px_vec);
    tree->Branch("phiFit_Kp_py", &phiFit_Kp_py_vec);
    tree->Branch("phiFit_Kp_pz", &phiFit_Kp_pz_vec);

    tree->Branch("phiFit_Km_eta", &phiFit_Km_eta_vec);
    tree->Branch("phiFit_Km_phi", &phiFit_Km_phi_vec);
    tree->Branch("phiFit_Km_p", &phiFit_Km_p_vec);
    tree->Branch("phiFit_Km_pt", &phiFit_Km_pt_vec);
    tree->Branch("phiFit_Km_px", &phiFit_Km_px_vec);
    tree->Branch("phiFit_Km_py", &phiFit_Km_py_vec);
    tree->Branch("phiFit_Km_pz", &phiFit_Km_pz_vec);

    tree->Branch("phiFit_pi_eta", &phiFit_pi_eta_vec);
    tree->Branch("phiFit_pi_phi", &phiFit_pi_phi_vec);
    tree->Branch("phiFit_pi_p", &phiFit_pi_p_vec);
    tree->Branch("phiFit_pi_pt", &phiFit_pi_pt_vec);
    tree->Branch("phiFit_pi_px", &phiFit_pi_px_vec);
    tree->Branch("phiFit_pi_py", &phiFit_pi_py_vec);
    tree->Branch("phiFit_pi_pz", &phiFit_pi_pz_vec);

    tree->Branch("phiFit_phi_eta", &phiFit_phi_eta_vec);
    tree->Branch("phiFit_phi_phi", &phiFit_phi_phi_vec);
    tree->Branch("phiFit_phi_p", &phiFit_phi_p_vec);
    tree->Branch("phiFit_phi_pt", &phiFit_phi_pt_vec);
    tree->Branch("phiFit_phi_px", &phiFit_phi_px_vec);
    tree->Branch("phiFit_phi_py", &phiFit_phi_py_vec);
    tree->Branch("phiFit_phi_pz", &phiFit_phi_pz_vec);
    tree->Branch("phiFit_phi_invm", &phiFit_phi_invm_vec);

    tree->Branch("phiFit_Ds_eta", &phiFit_Ds_eta_vec);
    tree->Branch("phiFit_Ds_phi", &phiFit_Ds_phi_vec);
    tree->Branch("phiFit_Ds_p", &phiFit_Ds_p_vec);
    tree->Branch("phiFit_Ds_pt", &phiFit_Ds_pt_vec);
    tree->Branch("phiFit_Ds_px", &phiFit_Ds_px_vec);
    tree->Branch("phiFit_Ds_py", &phiFit_Ds_py_vec);
    tree->Branch("phiFit_Ds_pz", &phiFit_Ds_pz_vec);
    tree->Branch("phiFit_Ds_invm", &phiFit_Ds_invm_vec);

    tree->Branch("phiFit_dR_Kp_Km", &phiFit_dR_Kp_Km_vec);
    tree->Branch("phiFit_dR_Kp_phi", &phiFit_dR_Kp_phi_vec);
    tree->Branch("phiFit_dR_Km_phi", &phiFit_dR_Km_phi_vec);
    tree->Branch("phiFit_dR_Kp_pi", &phiFit_dR_Kp_pi_vec);
    tree->Branch("phiFit_dR_Km_pi", &phiFit_dR_Km_pi_vec);
    tree->Branch("phiFit_dR_pi_phi", &phiFit_dR_pi_phi_vec);
    tree->Branch("phiFit_dR_Kp_Ds", &phiFit_dR_Kp_Ds_vec);
    tree->Branch("phiFit_dR_Km_Ds", &phiFit_dR_Km_Ds_vec);
    tree->Branch("phiFit_dR_phi_Ds", &phiFit_dR_phi_Ds_vec);
    tree->Branch("phiFit_dR_pi_Ds", &phiFit_dR_pi_Ds_vec);

    // Fit on Ds 
    tree->Branch("DsFit_chi2", &DsFit_chi2_vec);
    tree->Branch("DsFit_ndof", &DsFit_ndof_vec);
    tree->Branch("DsFit_chi2ndof", &DsFit_chi2ndof_vec);
    tree->Branch("DsFit_vx", &DsFit_vx_vec);
    tree->Branch("DsFit_vy", &DsFit_vy_vec);
    tree->Branch("DsFit_vz", &DsFit_vz_vec);
    tree->Branch("DsFit_vxerr", &DsFit_vxerr_vec);
    tree->Branch("DsFit_vyerr", &DsFit_vyerr_vec);
    tree->Branch("DsFit_vzerr", &DsFit_vzerr_vec);

    tree->Branch("DsFit_Kp_eta", &DsFit_Kp_eta_vec);
    tree->Branch("DsFit_Kp_phi", &DsFit_Kp_phi_vec);
    tree->Branch("DsFit_Kp_p", &DsFit_Kp_p_vec);
    tree->Branch("DsFit_Kp_pt", &DsFit_Kp_pt_vec);
    tree->Branch("DsFit_Kp_px", &DsFit_Kp_px_vec);
    tree->Branch("DsFit_Kp_py", &DsFit_Kp_py_vec);
    tree->Branch("DsFit_Kp_pz", &DsFit_Kp_pz_vec);

    tree->Branch("DsFit_Km_eta", &DsFit_Km_eta_vec);
    tree->Branch("DsFit_Km_phi", &DsFit_Km_phi_vec);
    tree->Branch("DsFit_Km_p", &DsFit_Km_p_vec);
    tree->Branch("DsFit_Km_pt", &DsFit_Km_pt_vec);
    tree->Branch("DsFit_Km_px", &DsFit_Km_px_vec);
    tree->Branch("DsFit_Km_py", &DsFit_Km_py_vec);
    tree->Branch("DsFit_Km_pz", &DsFit_Km_pz_vec);

    tree->Branch("DsFit_pi_eta", &DsFit_pi_eta_vec);
    tree->Branch("DsFit_pi_phi", &DsFit_pi_phi_vec);
    tree->Branch("DsFit_pi_p", &DsFit_pi_p_vec);
    tree->Branch("DsFit_pi_pt", &DsFit_pi_pt_vec);
    tree->Branch("DsFit_pi_px", &DsFit_pi_px_vec);
    tree->Branch("DsFit_pi_py", &DsFit_pi_py_vec);
    tree->Branch("DsFit_pi_pz", &DsFit_pi_pz_vec);

    tree->Branch("DsFit_phi_eta", &DsFit_phi_eta_vec);
    tree->Branch("DsFit_phi_phi", &DsFit_phi_phi_vec);
    tree->Branch("DsFit_phi_p", &DsFit_phi_p_vec);
    tree->Branch("DsFit_phi_pt", &DsFit_phi_pt_vec);
    tree->Branch("DsFit_phi_px", &DsFit_phi_px_vec);
    tree->Branch("DsFit_phi_py", &DsFit_phi_py_vec);
    tree->Branch("DsFit_phi_pz", &DsFit_phi_pz_vec);
    tree->Branch("DsFit_phi_invm", &DsFit_phi_invm_vec);

    tree->Branch("DsFit_Ds_eta", &DsFit_Ds_eta_vec);
    tree->Branch("DsFit_Ds_phi", &DsFit_Ds_phi_vec);
    tree->Branch("DsFit_Ds_p", &DsFit_Ds_p_vec);
    tree->Branch("DsFit_Ds_pt", &DsFit_Ds_pt_vec);
    tree->Branch("DsFit_Ds_px", &DsFit_Ds_px_vec);
    tree->Branch("DsFit_Ds_py", &DsFit_Ds_py_vec);
    tree->Branch("DsFit_Ds_pz", &DsFit_Ds_pz_vec);
    tree->Branch("DsFit_Ds_invm", &DsFit_Ds_invm_vec);

    tree->Branch("DsFit_dR_Kp_Km", &DsFit_dR_Kp_Km_vec);
    tree->Branch("DsFit_dR_Kp_phi", &DsFit_dR_Kp_phi_vec);
    tree->Branch("DsFit_dR_Km_phi", &DsFit_dR_Km_phi_vec);
    tree->Branch("DsFit_dR_Kp_pi", &DsFit_dR_Kp_pi_vec);
    tree->Branch("DsFit_dR_Km_pi", &DsFit_dR_Km_pi_vec);
    tree->Branch("DsFit_dR_pi_phi", &DsFit_dR_pi_phi_vec);
    tree->Branch("DsFit_dR_Kp_Ds", &DsFit_dR_Kp_Ds_vec);
    tree->Branch("DsFit_dR_Km_Ds", &DsFit_dR_Km_Ds_vec);
    tree->Branch("DsFit_dR_phi_Ds", &DsFit_dR_phi_Ds_vec);
    tree->Branch("DsFit_dR_pi_Ds", &DsFit_dR_pi_Ds_vec);

    tree->Branch("DsFit_Mconstraint_Ds_invm", &DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("Ds_primvtx_FDxy", &Ds_primvtx_FDxy_vec);
    tree->Branch("Ds_primvtx_FDz", &Ds_primvtx_FDz_vec);
    tree->Branch("Ds_primvtx_FD", &Ds_primvtx_FD_vec);
    tree->Branch("Ds_primvtx_FDxyerr", &Ds_primvtx_FDxyerr_vec);
    tree->Branch("Ds_primvtx_FDxychi2", &Ds_primvtx_FDxychi2_vec);
    tree->Branch("Ds_primvtx_FDzerr", &Ds_primvtx_FDzerr_vec);
    tree->Branch("Ds_primvtx_FDzchi2", &Ds_primvtx_FDzchi2_vec);
    tree->Branch("Ds_primvtx_FDerr", &Ds_primvtx_FDerr_vec);
    tree->Branch("Ds_primvtx_FDchi2", &Ds_primvtx_FDchi2_vec);
    tree->Branch("Ds_primvtx_dira", &Ds_primvtx_dira_vec);
    tree->Branch("Ds_primvtx_dira_angle", &Ds_primvtx_dira_angle_vec);
    tree->Branch("Kp_primvtx_ip", &Kp_primvtx_ip_vec);
    tree->Branch("Kp_primvtx_iperr", &Kp_primvtx_iperr_vec);
    tree->Branch("Kp_primvtx_ipchi2", &Kp_primvtx_ipchi2_vec);
    tree->Branch("Km_primvtx_ip", &Km_primvtx_ip_vec);
    tree->Branch("Km_primvtx_iperr", &Km_primvtx_iperr_vec);
    tree->Branch("Km_primvtx_ipchi2", &Km_primvtx_ipchi2_vec);
    tree->Branch("pi_primvtx_ip", &pi_primvtx_ip_vec);
    tree->Branch("pi_primvtx_iperr", &pi_primvtx_iperr_vec);
    tree->Branch("pi_primvtx_ipchi2", &pi_primvtx_ipchi2_vec);
    tree->Branch("phi_primvtx_ip", &phi_primvtx_ip_vec);
    tree->Branch("phi_primvtx_iperr", &phi_primvtx_iperr_vec);
    tree->Branch("phi_primvtx_ipchi2", &phi_primvtx_ipchi2_vec);
    tree->Branch("Ds_primvtx_ip", &Ds_primvtx_ip_vec);
    tree->Branch("Ds_primvtx_iperr", &Ds_primvtx_iperr_vec);
    tree->Branch("Ds_primvtx_ipchi2", &Ds_primvtx_ipchi2_vec);

    tree->Branch("Ds_IsoR03_sumChargedHadronPt", &Ds_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("Ds_IsoR03_sumNeutralHadronPt", &Ds_IsoR03_sumNeutralHadronPt_vec);
    tree->Branch("Ds_IsoR03_sumPhotonPt", &Ds_IsoR03_sumPhotonPt_vec);
    tree->Branch("Ds_IsoR03_sumPUPt", &Ds_IsoR03_sumPUPt_vec);
    tree->Branch("Ds_PFIsoR03", &Ds_PFIsoR03_vec);
    tree->Branch("Ds_IsoR04_sumChargedHadronPt", &Ds_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("Ds_IsoR04_sumNeutralHadronPt", &Ds_IsoR04_sumNeutralHadronPt_vec);
    tree->Branch("Ds_IsoR04_sumPhotonPt", &Ds_IsoR04_sumPhotonPt_vec);
    tree->Branch("Ds_IsoR04_sumPUPt", &Ds_IsoR04_sumPUPt_vec);
    tree->Branch("Ds_PFIsoR04", &Ds_PFIsoR04_vec);

    tree->Branch("PV_noDs_IsValid", &PV_noDs_IsValid_vec);
    tree->Branch("PV_noDs_IsFake", &PV_noDs_IsFake_vec);
    tree->Branch("PV_noDs_chi2", &PV_noDs_chi2_vec);
    tree->Branch("PV_noDs_ndof", &PV_noDs_ndof_vec);
    tree->Branch("PV_noDs_chi2ndof", &PV_noDs_chi2ndof_vec);
    tree->Branch("PV_noDs_x", &PV_noDs_x_vec);
    tree->Branch("PV_noDs_y", &PV_noDs_y_vec);
    tree->Branch("PV_noDs_z", &PV_noDs_z_vec);
    tree->Branch("PV_noDs_xerr", &PV_noDs_xerr_vec);
    tree->Branch("PV_noDs_yerr", &PV_noDs_yerr_vec);
    tree->Branch("PV_noDs_zerr", &PV_noDs_zerr_vec);

    tree->Branch("Ds_PVnoDs_FDxy", &Ds_PVnoDs_FDxy_vec);
    tree->Branch("Ds_PVnoDs_FDz", &Ds_PVnoDs_FDz_vec);
    tree->Branch("Ds_PVnoDs_FD", &Ds_PVnoDs_FD_vec);
    tree->Branch("Ds_PVnoDs_FDxyerr", &Ds_PVnoDs_FDxyerr_vec);
    tree->Branch("Ds_PVnoDs_FDxychi2", &Ds_PVnoDs_FDxychi2_vec);
    tree->Branch("Ds_PVnoDs_FDzerr", &Ds_PVnoDs_FDzerr_vec);
    tree->Branch("Ds_PVnoDs_FDzchi2", &Ds_PVnoDs_FDzchi2_vec);
    tree->Branch("Ds_PVnoDs_FDerr", &Ds_PVnoDs_FDerr_vec);
    tree->Branch("Ds_PVnoDs_FDchi2", &Ds_PVnoDs_FDchi2_vec);
    tree->Branch("Ds_PVnoDs_dira", &Ds_PVnoDs_dira_vec);
    tree->Branch("Ds_PVnoDs_dira_angle", &Ds_PVnoDs_dira_angle_vec);
    tree->Branch("Kp_PVnoDs_ip", &Kp_PVnoDs_ip_vec);
    tree->Branch("Kp_PVnoDs_iperr", &Kp_PVnoDs_iperr_vec);
    tree->Branch("Kp_PVnoDs_ipchi2", &Kp_PVnoDs_ipchi2_vec);
    tree->Branch("Km_PVnoDs_ip", &Km_PVnoDs_ip_vec);
    tree->Branch("Km_PVnoDs_iperr", &Km_PVnoDs_iperr_vec);
    tree->Branch("Km_PVnoDs_ipchi2", &Km_PVnoDs_ipchi2_vec);
    tree->Branch("pi_PVnoDs_ip", &pi_PVnoDs_ip_vec);
    tree->Branch("pi_PVnoDs_iperr", &pi_PVnoDs_iperr_vec);
    tree->Branch("pi_PVnoDs_ipchi2", &pi_PVnoDs_ipchi2_vec);
    tree->Branch("phi_PVnoDs_ip", &phi_PVnoDs_ip_vec);
    tree->Branch("phi_PVnoDs_iperr", &phi_PVnoDs_iperr_vec);
    tree->Branch("phi_PVnoDs_ipchi2", &phi_PVnoDs_ipchi2_vec);
    tree->Branch("Ds_PVnoDs_ip", &Ds_PVnoDs_ip_vec);
    tree->Branch("Ds_PVnoDs_iperr", &Ds_PVnoDs_iperr_vec);
    tree->Branch("Ds_PVnoDs_ipchi2", &Ds_PVnoDs_ipchi2_vec);

    tree->Branch("PV_withDs_IsValid", &PV_withDs_IsValid_vec);
    tree->Branch("PV_withDs_IsFake", &PV_withDs_IsFake_vec);
    tree->Branch("PV_withDs_chi2", &PV_withDs_chi2_vec);
    tree->Branch("PV_withDs_ndof", &PV_withDs_ndof_vec);
    tree->Branch("PV_withDs_chi2ndof", &PV_withDs_chi2ndof_vec);
    tree->Branch("PV_withDs_x", &PV_withDs_x_vec);
    tree->Branch("PV_withDs_y", &PV_withDs_y_vec);
    tree->Branch("PV_withDs_z", &PV_withDs_z_vec);
    tree->Branch("PV_withDs_xerr", &PV_withDs_xerr_vec);
    tree->Branch("PV_withDs_yerr", &PV_withDs_yerr_vec);
    tree->Branch("PV_withDs_zerr", &PV_withDs_zerr_vec);

    tree->Branch("Ds_PVwithDs_FDxy", &Ds_PVwithDs_FDxy_vec);
    tree->Branch("Ds_PVwithDs_FDz", &Ds_PVwithDs_FDz_vec);
    tree->Branch("Ds_PVwithDs_FD", &Ds_PVwithDs_FD_vec);
    tree->Branch("Ds_PVwithDs_FDxyerr", &Ds_PVwithDs_FDxyerr_vec);
    tree->Branch("Ds_PVwithDs_FDxychi2", &Ds_PVwithDs_FDxychi2_vec);
    tree->Branch("Ds_PVwithDs_FDzerr", &Ds_PVwithDs_FDzerr_vec);
    tree->Branch("Ds_PVwithDs_FDzchi2", &Ds_PVwithDs_FDzchi2_vec);
    tree->Branch("Ds_PVwithDs_FDerr", &Ds_PVwithDs_FDerr_vec);
    tree->Branch("Ds_PVwithDs_FDchi2", &Ds_PVwithDs_FDchi2_vec);
    tree->Branch("Ds_PVwithDs_dira", &Ds_PVwithDs_dira_vec);
    tree->Branch("Ds_PVwithDs_dira_angle", &Ds_PVwithDs_dira_angle_vec);
    tree->Branch("Kp_PVwithDs_ip", &Kp_PVwithDs_ip_vec);
    tree->Branch("Kp_PVwithDs_iperr", &Kp_PVwithDs_iperr_vec);
    tree->Branch("Kp_PVwithDs_ipchi2", &Kp_PVwithDs_ipchi2_vec);
    tree->Branch("Km_PVwithDs_ip", &Km_PVwithDs_ip_vec);
    tree->Branch("Km_PVwithDs_iperr", &Km_PVwithDs_iperr_vec);
    tree->Branch("Km_PVwithDs_ipchi2", &Km_PVwithDs_ipchi2_vec);
    tree->Branch("pi_PVwithDs_ip", &pi_PVwithDs_ip_vec);
    tree->Branch("pi_PVwithDs_iperr", &pi_PVwithDs_iperr_vec);
    tree->Branch("pi_PVwithDs_ipchi2", &pi_PVwithDs_ipchi2_vec);
    tree->Branch("phi_PVwithDs_ip", &phi_PVwithDs_ip_vec);
    tree->Branch("phi_PVwithDs_iperr", &phi_PVwithDs_iperr_vec);
    tree->Branch("phi_PVwithDs_ipchi2", &phi_PVwithDs_ipchi2_vec);
    tree->Branch("Ds_PVwithDs_ip", &Ds_PVwithDs_ip_vec);
    tree->Branch("Ds_PVwithDs_iperr", &Ds_PVwithDs_iperr_vec);
    tree->Branch("Ds_PVwithDs_ipchi2", &Ds_PVwithDs_ipchi2_vec);

    tree->Branch("Kp_match", &Kp_match_vec);
    tree->Branch("Km_match", &Km_match_vec);
    tree->Branch("pi_match", &pi_match_vec);
    tree->Branch("match_entry", &match_entry_vec);
    tree->Branch("non_match_entry", &non_match_entry_vec);

    tree->Branch("best_Kp_isIsolatedChargedHadron", &best_Kp_isIsolatedChargedHadron_vec);
    tree->Branch("best_Kp_charge", &best_Kp_charge_vec);
    tree->Branch("best_Kp_eta", &best_Kp_eta_vec);
    tree->Branch("best_Kp_phi", &best_Kp_phi_vec);
    tree->Branch("best_Kp_vx", &best_Kp_vx_vec);
    tree->Branch("best_Kp_vy", &best_Kp_vy_vec);
    tree->Branch("best_Kp_vz", &best_Kp_vz_vec);
    tree->Branch("best_Kp_p", &best_Kp_p_vec);
    tree->Branch("best_Kp_pt", &best_Kp_pt_vec);
    tree->Branch("best_Kp_px", &best_Kp_px_vec);
    tree->Branch("best_Kp_py", &best_Kp_py_vec);
    tree->Branch("best_Kp_pz", &best_Kp_pz_vec);

    tree->Branch("best_Km_isIsolatedChargedHadron", &best_Km_isIsolatedChargedHadron_vec);
    tree->Branch("best_Km_charge", &best_Km_charge_vec);
    tree->Branch("best_Km_eta", &best_Km_eta_vec);
    tree->Branch("best_Km_phi", &best_Km_phi_vec);
    tree->Branch("best_Km_vx", &best_Km_vx_vec);
    tree->Branch("best_Km_vy", &best_Km_vy_vec);
    tree->Branch("best_Km_vz", &best_Km_vz_vec);
    tree->Branch("best_Km_p", &best_Km_p_vec);
    tree->Branch("best_Km_pt", &best_Km_pt_vec);
    tree->Branch("best_Km_px", &best_Km_px_vec);
    tree->Branch("best_Km_py", &best_Km_py_vec);
    tree->Branch("best_Km_pz", &best_Km_pz_vec);

    tree->Branch("best_pi_isIsolatedChargedHadron", &best_pi_isIsolatedChargedHadron_vec);
    tree->Branch("best_pi_charge", &best_pi_charge_vec);
    tree->Branch("best_pi_eta", &best_pi_eta_vec);
    tree->Branch("best_pi_phi", &best_pi_phi_vec);
    tree->Branch("best_pi_vx", &best_pi_vx_vec);
    tree->Branch("best_pi_vy", &best_pi_vy_vec);
    tree->Branch("best_pi_vz", &best_pi_vz_vec);
    tree->Branch("best_pi_p", &best_pi_p_vec);
    tree->Branch("best_pi_pt", &best_pi_pt_vec);
    tree->Branch("best_pi_px", &best_pi_px_vec);
    tree->Branch("best_pi_py", &best_pi_py_vec);
    tree->Branch("best_pi_pz", &best_pi_pz_vec);

    tree->Branch("best_phi_eta", &best_phi_eta_vec);
    tree->Branch("best_phi_phi", &best_phi_phi_vec);
    tree->Branch("best_phi_p", &best_phi_p_vec);
    tree->Branch("best_phi_pt", &best_phi_pt_vec);
    tree->Branch("best_phi_px", &best_phi_px_vec);
    tree->Branch("best_phi_py", &best_phi_py_vec);
    tree->Branch("best_phi_pz", &best_phi_pz_vec);
    tree->Branch("best_phi_invm", &best_phi_invm_vec);

    tree->Branch("best_Ds_eta", &best_Ds_eta_vec);
    tree->Branch("best_Ds_phi", &best_Ds_phi_vec);
    tree->Branch("best_Ds_p", &best_Ds_p_vec);
    tree->Branch("best_Ds_pt", &best_Ds_pt_vec);
    tree->Branch("best_Ds_px", &best_Ds_px_vec);
    tree->Branch("best_Ds_py", &best_Ds_py_vec);
    tree->Branch("best_Ds_pz", &best_Ds_pz_vec);
    tree->Branch("best_Ds_invm", &best_Ds_invm_vec);

    tree->Branch("best_dR_Kp_Km", &best_dR_Kp_Km_vec);
    tree->Branch("best_dR_Kp_phi", &best_dR_Kp_phi_vec);
    tree->Branch("best_dR_Km_phi", &best_dR_Km_phi_vec);
    tree->Branch("best_dR_Kp_pi", &best_dR_Kp_pi_vec);
    tree->Branch("best_dR_Km_pi", &best_dR_Km_pi_vec);
    tree->Branch("best_dR_pi_phi", &best_dR_pi_phi_vec);
    tree->Branch("best_dR_Kp_Ds", &best_dR_Kp_Ds_vec);
    tree->Branch("best_dR_Km_Ds", &best_dR_Km_Ds_vec);
    tree->Branch("best_dR_phi_Ds", &best_dR_phi_Ds_vec);
    tree->Branch("best_dR_pi_Ds", &best_dR_pi_Ds_vec);

    tree->Branch("best_dxy_phi_Ds", &best_dxy_phi_Ds_vec);
    tree->Branch("best_dz_phi_Ds", &best_dz_phi_Ds_vec);

    // Fit on phi
    tree->Branch("best_phiFit_chi2", &best_phiFit_chi2_vec);
    tree->Branch("best_phiFit_ndof", &best_phiFit_ndof_vec);
    tree->Branch("best_phiFit_chi2ndof", &best_phiFit_chi2ndof_vec);
    tree->Branch("best_phiFit_vx", &best_phiFit_vx_vec);
    tree->Branch("best_phiFit_vy", &best_phiFit_vy_vec);
    tree->Branch("best_phiFit_vz", &best_phiFit_vz_vec);
    tree->Branch("best_phiFit_vxerr", &best_phiFit_vxerr_vec);
    tree->Branch("best_phiFit_vyerr", &best_phiFit_vyerr_vec);
    tree->Branch("best_phiFit_vzerr", &best_phiFit_vzerr_vec);

    tree->Branch("best_phiFit_Kp_eta", &best_phiFit_Kp_eta_vec);
    tree->Branch("best_phiFit_Kp_phi", &best_phiFit_Kp_phi_vec);
    tree->Branch("best_phiFit_Kp_p", &best_phiFit_Kp_p_vec);
    tree->Branch("best_phiFit_Kp_pt", &best_phiFit_Kp_pt_vec);
    tree->Branch("best_phiFit_Kp_px", &best_phiFit_Kp_px_vec);
    tree->Branch("best_phiFit_Kp_py", &best_phiFit_Kp_py_vec);
    tree->Branch("best_phiFit_Kp_pz", &best_phiFit_Kp_pz_vec);

    tree->Branch("best_phiFit_Km_eta", &best_phiFit_Km_eta_vec);
    tree->Branch("best_phiFit_Km_phi", &best_phiFit_Km_phi_vec);
    tree->Branch("best_phiFit_Km_p", &best_phiFit_Km_p_vec);
    tree->Branch("best_phiFit_Km_pt", &best_phiFit_Km_pt_vec);
    tree->Branch("best_phiFit_Km_px", &best_phiFit_Km_px_vec);
    tree->Branch("best_phiFit_Km_py", &best_phiFit_Km_py_vec);
    tree->Branch("best_phiFit_Km_pz", &best_phiFit_Km_pz_vec);

    tree->Branch("best_phiFit_pi_eta", &best_phiFit_pi_eta_vec);
    tree->Branch("best_phiFit_pi_phi", &best_phiFit_pi_phi_vec);
    tree->Branch("best_phiFit_pi_p", &best_phiFit_pi_p_vec);
    tree->Branch("best_phiFit_pi_pt", &best_phiFit_pi_pt_vec);
    tree->Branch("best_phiFit_pi_px", &best_phiFit_pi_px_vec);
    tree->Branch("best_phiFit_pi_py", &best_phiFit_pi_py_vec);
    tree->Branch("best_phiFit_pi_pz", &best_phiFit_pi_pz_vec);

    tree->Branch("best_phiFit_phi_eta", &best_phiFit_phi_eta_vec);
    tree->Branch("best_phiFit_phi_phi", &best_phiFit_phi_phi_vec);
    tree->Branch("best_phiFit_phi_p", &best_phiFit_phi_p_vec);
    tree->Branch("best_phiFit_phi_pt", &best_phiFit_phi_pt_vec);
    tree->Branch("best_phiFit_phi_px", &best_phiFit_phi_px_vec);
    tree->Branch("best_phiFit_phi_py", &best_phiFit_phi_py_vec);
    tree->Branch("best_phiFit_phi_pz", &best_phiFit_phi_pz_vec);
    tree->Branch("best_phiFit_phi_invm", &best_phiFit_phi_invm_vec);

    tree->Branch("best_phiFit_Ds_eta", &best_phiFit_Ds_eta_vec);
    tree->Branch("best_phiFit_Ds_phi", &best_phiFit_Ds_phi_vec);
    tree->Branch("best_phiFit_Ds_p", &best_phiFit_Ds_p_vec);
    tree->Branch("best_phiFit_Ds_pt", &best_phiFit_Ds_pt_vec);
    tree->Branch("best_phiFit_Ds_px", &best_phiFit_Ds_px_vec);
    tree->Branch("best_phiFit_Ds_py", &best_phiFit_Ds_py_vec);
    tree->Branch("best_phiFit_Ds_pz", &best_phiFit_Ds_pz_vec);
    tree->Branch("best_phiFit_Ds_invm", &best_phiFit_Ds_invm_vec);

    tree->Branch("best_phiFit_dR_Kp_Km", &best_phiFit_dR_Kp_Km_vec);
    tree->Branch("best_phiFit_dR_Kp_phi", &best_phiFit_dR_Kp_phi_vec);
    tree->Branch("best_phiFit_dR_Km_phi", &best_phiFit_dR_Km_phi_vec);
    tree->Branch("best_phiFit_dR_Kp_pi", &best_phiFit_dR_Kp_pi_vec);
    tree->Branch("best_phiFit_dR_Km_pi", &best_phiFit_dR_Km_pi_vec);
    tree->Branch("best_phiFit_dR_pi_phi", &best_phiFit_dR_pi_phi_vec);
    tree->Branch("best_phiFit_dR_Kp_Ds", &best_phiFit_dR_Kp_Ds_vec);
    tree->Branch("best_phiFit_dR_Km_Ds", &best_phiFit_dR_Km_Ds_vec);
    tree->Branch("best_phiFit_dR_phi_Ds", &best_phiFit_dR_phi_Ds_vec);
    tree->Branch("best_phiFit_dR_pi_Ds", &best_phiFit_dR_pi_Ds_vec);

    // Fit on Ds 
    tree->Branch("best_DsFit_chi2", &best_DsFit_chi2_vec);
    tree->Branch("best_DsFit_ndof", &best_DsFit_ndof_vec);
    tree->Branch("best_DsFit_chi2ndof", &best_DsFit_chi2ndof_vec);
    tree->Branch("best_DsFit_vx", &best_DsFit_vx_vec);
    tree->Branch("best_DsFit_vy", &best_DsFit_vy_vec);
    tree->Branch("best_DsFit_vz", &best_DsFit_vz_vec);
    tree->Branch("best_DsFit_vxerr", &best_DsFit_vxerr_vec);
    tree->Branch("best_DsFit_vyerr", &best_DsFit_vyerr_vec);
    tree->Branch("best_DsFit_vzerr", &best_DsFit_vzerr_vec);

    tree->Branch("best_DsFit_Kp_eta", &best_DsFit_Kp_eta_vec);
    tree->Branch("best_DsFit_Kp_phi", &best_DsFit_Kp_phi_vec);
    tree->Branch("best_DsFit_Kp_p", &best_DsFit_Kp_p_vec);
    tree->Branch("best_DsFit_Kp_pt", &best_DsFit_Kp_pt_vec);
    tree->Branch("best_DsFit_Kp_px", &best_DsFit_Kp_px_vec);
    tree->Branch("best_DsFit_Kp_py", &best_DsFit_Kp_py_vec);
    tree->Branch("best_DsFit_Kp_pz", &best_DsFit_Kp_pz_vec);

    tree->Branch("best_DsFit_Km_eta", &best_DsFit_Km_eta_vec);
    tree->Branch("best_DsFit_Km_phi", &best_DsFit_Km_phi_vec);
    tree->Branch("best_DsFit_Km_p", &best_DsFit_Km_p_vec);
    tree->Branch("best_DsFit_Km_pt", &best_DsFit_Km_pt_vec);
    tree->Branch("best_DsFit_Km_px", &best_DsFit_Km_px_vec);
    tree->Branch("best_DsFit_Km_py", &best_DsFit_Km_py_vec);
    tree->Branch("best_DsFit_Km_pz", &best_DsFit_Km_pz_vec);

    tree->Branch("best_DsFit_pi_eta", &best_DsFit_pi_eta_vec);
    tree->Branch("best_DsFit_pi_phi", &best_DsFit_pi_phi_vec);
    tree->Branch("best_DsFit_pi_p", &best_DsFit_pi_p_vec);
    tree->Branch("best_DsFit_pi_pt", &best_DsFit_pi_pt_vec);
    tree->Branch("best_DsFit_pi_px", &best_DsFit_pi_px_vec);
    tree->Branch("best_DsFit_pi_py", &best_DsFit_pi_py_vec);
    tree->Branch("best_DsFit_pi_pz", &best_DsFit_pi_pz_vec);

    tree->Branch("best_DsFit_phi_eta", &best_DsFit_phi_eta_vec);
    tree->Branch("best_DsFit_phi_phi", &best_DsFit_phi_phi_vec);
    tree->Branch("best_DsFit_phi_p", &best_DsFit_phi_p_vec);
    tree->Branch("best_DsFit_phi_pt", &best_DsFit_phi_pt_vec);
    tree->Branch("best_DsFit_phi_px", &best_DsFit_phi_px_vec);
    tree->Branch("best_DsFit_phi_py", &best_DsFit_phi_py_vec);
    tree->Branch("best_DsFit_phi_pz", &best_DsFit_phi_pz_vec);
    tree->Branch("best_DsFit_phi_invm", &best_DsFit_phi_invm_vec);

    tree->Branch("best_DsFit_Ds_eta", &best_DsFit_Ds_eta_vec);
    tree->Branch("best_DsFit_Ds_phi", &best_DsFit_Ds_phi_vec);
    tree->Branch("best_DsFit_Ds_p", &best_DsFit_Ds_p_vec);
    tree->Branch("best_DsFit_Ds_pt", &best_DsFit_Ds_pt_vec);
    tree->Branch("best_DsFit_Ds_px", &best_DsFit_Ds_px_vec);
    tree->Branch("best_DsFit_Ds_py", &best_DsFit_Ds_py_vec);
    tree->Branch("best_DsFit_Ds_pz", &best_DsFit_Ds_pz_vec);
    tree->Branch("best_DsFit_Ds_invm", &best_DsFit_Ds_invm_vec);

    tree->Branch("best_DsFit_dR_Kp_Km", &best_DsFit_dR_Kp_Km_vec);
    tree->Branch("best_DsFit_dR_Kp_phi", &best_DsFit_dR_Kp_phi_vec);
    tree->Branch("best_DsFit_dR_Km_phi", &best_DsFit_dR_Km_phi_vec);
    tree->Branch("best_DsFit_dR_Kp_pi", &best_DsFit_dR_Kp_pi_vec);
    tree->Branch("best_DsFit_dR_Km_pi", &best_DsFit_dR_Km_pi_vec);
    tree->Branch("best_DsFit_dR_pi_phi", &best_DsFit_dR_pi_phi_vec);
    tree->Branch("best_DsFit_dR_Kp_Ds", &best_DsFit_dR_Kp_Ds_vec);
    tree->Branch("best_DsFit_dR_Km_Ds", &best_DsFit_dR_Km_Ds_vec);
    tree->Branch("best_DsFit_dR_phi_Ds", &best_DsFit_dR_phi_Ds_vec);
    tree->Branch("best_DsFit_dR_pi_Ds", &best_DsFit_dR_pi_Ds_vec);

    tree->Branch("best_DsFit_Mconstraint_Ds_invm", &best_DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("best_Ds_primvtx_FDxy", &best_Ds_primvtx_FDxy_vec);
    tree->Branch("best_Ds_primvtx_FDz", &best_Ds_primvtx_FDz_vec);
    tree->Branch("best_Ds_primvtx_FD", &best_Ds_primvtx_FD_vec);
    tree->Branch("best_Ds_primvtx_FDxyerr", &best_Ds_primvtx_FDxyerr_vec);
    tree->Branch("best_Ds_primvtx_FDxychi2", &best_Ds_primvtx_FDxychi2_vec);
    tree->Branch("best_Ds_primvtx_FDzerr", &best_Ds_primvtx_FDzerr_vec);
    tree->Branch("best_Ds_primvtx_FDzchi2", &best_Ds_primvtx_FDzchi2_vec);
    tree->Branch("best_Ds_primvtx_FDerr", &best_Ds_primvtx_FDerr_vec);
    tree->Branch("best_Ds_primvtx_FDchi2", &best_Ds_primvtx_FDchi2_vec);
    tree->Branch("best_Ds_primvtx_dira", &best_Ds_primvtx_dira_vec);
    tree->Branch("best_Ds_primvtx_dira_angle", &best_Ds_primvtx_dira_angle_vec);
    tree->Branch("best_Kp_primvtx_ip", &best_Kp_primvtx_ip_vec);
    tree->Branch("best_Kp_primvtx_iperr", &best_Kp_primvtx_iperr_vec);
    tree->Branch("best_Kp_primvtx_ipchi2", &best_Kp_primvtx_ipchi2_vec);
    tree->Branch("best_Km_primvtx_ip", &best_Km_primvtx_ip_vec);
    tree->Branch("best_Km_primvtx_iperr", &best_Km_primvtx_iperr_vec);
    tree->Branch("best_Km_primvtx_ipchi2", &best_Km_primvtx_ipchi2_vec);
    tree->Branch("best_pi_primvtx_ip", &best_pi_primvtx_ip_vec);
    tree->Branch("best_pi_primvtx_iperr", &best_pi_primvtx_iperr_vec);
    tree->Branch("best_pi_primvtx_ipchi2", &best_pi_primvtx_ipchi2_vec);
    tree->Branch("best_phi_primvtx_ip", &best_phi_primvtx_ip_vec);
    tree->Branch("best_phi_primvtx_iperr", &best_phi_primvtx_iperr_vec);
    tree->Branch("best_phi_primvtx_ipchi2", &best_phi_primvtx_ipchi2_vec);
    tree->Branch("best_Ds_primvtx_ip", &best_Ds_primvtx_ip_vec);
    tree->Branch("best_Ds_primvtx_iperr", &best_Ds_primvtx_iperr_vec);
    tree->Branch("best_Ds_primvtx_ipchi2", &best_Ds_primvtx_ipchi2_vec);

    tree->Branch("best_Ds_IsoR03_sumChargedHadronPt", &best_Ds_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("best_Ds_IsoR03_sumNeutralHadronPt", &best_Ds_IsoR03_sumNeutralHadronPt_vec);
    tree->Branch("best_Ds_IsoR03_sumPhotonPt", &best_Ds_IsoR03_sumPhotonPt_vec);
    tree->Branch("best_Ds_IsoR03_sumPUPt", &best_Ds_IsoR03_sumPUPt_vec);
    tree->Branch("best_Ds_PFIsoR03", &best_Ds_PFIsoR03_vec);
    tree->Branch("best_Ds_IsoR04_sumChargedHadronPt", &best_Ds_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("best_Ds_IsoR04_sumNeutralHadronPt", &best_Ds_IsoR04_sumNeutralHadronPt_vec);
    tree->Branch("best_Ds_IsoR04_sumPhotonPt", &best_Ds_IsoR04_sumPhotonPt_vec);
    tree->Branch("best_Ds_IsoR04_sumPUPt", &best_Ds_IsoR04_sumPUPt_vec);
    tree->Branch("best_Ds_PFIsoR04", &best_Ds_PFIsoR04_vec);

    tree->Branch("best_PV_noDs_IsValid", &best_PV_noDs_IsValid_vec);
    tree->Branch("best_PV_noDs_IsFake", &best_PV_noDs_IsFake_vec);
    tree->Branch("best_PV_noDs_chi2", &best_PV_noDs_chi2_vec);
    tree->Branch("best_PV_noDs_ndof", &best_PV_noDs_ndof_vec);
    tree->Branch("best_PV_noDs_chi2ndof", &best_PV_noDs_chi2ndof_vec);
    tree->Branch("best_PV_noDs_x", &best_PV_noDs_x_vec);
    tree->Branch("best_PV_noDs_y", &best_PV_noDs_y_vec);
    tree->Branch("best_PV_noDs_z", &best_PV_noDs_z_vec);
    tree->Branch("best_PV_noDs_xerr", &best_PV_noDs_xerr_vec);
    tree->Branch("best_PV_noDs_yerr", &best_PV_noDs_yerr_vec);
    tree->Branch("best_PV_noDs_zerr", &best_PV_noDs_zerr_vec);

    tree->Branch("best_Ds_PVnoDs_FDxy", &best_Ds_PVnoDs_FDxy_vec);
    tree->Branch("best_Ds_PVnoDs_FDz", &best_Ds_PVnoDs_FDz_vec);
    tree->Branch("best_Ds_PVnoDs_FD", &best_Ds_PVnoDs_FD_vec);
    tree->Branch("best_Ds_PVnoDs_FDxyerr", &best_Ds_PVnoDs_FDxyerr_vec);
    tree->Branch("best_Ds_PVnoDs_FDxychi2", &best_Ds_PVnoDs_FDxychi2_vec);
    tree->Branch("best_Ds_PVnoDs_FDzerr", &best_Ds_PVnoDs_FDzerr_vec);
    tree->Branch("best_Ds_PVnoDs_FDzchi2", &best_Ds_PVnoDs_FDzchi2_vec);
    tree->Branch("best_Ds_PVnoDs_FDerr", &best_Ds_PVnoDs_FDerr_vec);
    tree->Branch("best_Ds_PVnoDs_FDchi2", &best_Ds_PVnoDs_FDchi2_vec);
    tree->Branch("best_Ds_PVnoDs_dira", &best_Ds_PVnoDs_dira_vec);
    tree->Branch("best_Ds_PVnoDs_dira_angle", &best_Ds_PVnoDs_dira_angle_vec);
    tree->Branch("best_Kp_PVnoDs_ip", &best_Kp_PVnoDs_ip_vec);
    tree->Branch("best_Kp_PVnoDs_iperr", &best_Kp_PVnoDs_iperr_vec);
    tree->Branch("best_Kp_PVnoDs_ipchi2", &best_Kp_PVnoDs_ipchi2_vec);
    tree->Branch("best_Km_PVnoDs_ip", &best_Km_PVnoDs_ip_vec);
    tree->Branch("best_Km_PVnoDs_iperr", &best_Km_PVnoDs_iperr_vec);
    tree->Branch("best_Km_PVnoDs_ipchi2", &best_Km_PVnoDs_ipchi2_vec);
    tree->Branch("best_pi_PVnoDs_ip", &best_pi_PVnoDs_ip_vec);
    tree->Branch("best_pi_PVnoDs_iperr", &best_pi_PVnoDs_iperr_vec);
    tree->Branch("best_pi_PVnoDs_ipchi2", &best_pi_PVnoDs_ipchi2_vec);
    tree->Branch("best_phi_PVnoDs_ip", &best_phi_PVnoDs_ip_vec);
    tree->Branch("best_phi_PVnoDs_iperr", &best_phi_PVnoDs_iperr_vec);
    tree->Branch("best_phi_PVnoDs_ipchi2", &best_phi_PVnoDs_ipchi2_vec);
    tree->Branch("best_Ds_PVnoDs_ip", &best_Ds_PVnoDs_ip_vec);
    tree->Branch("best_Ds_PVnoDs_iperr", &best_Ds_PVnoDs_iperr_vec);
    tree->Branch("best_Ds_PVnoDs_ipchi2", &best_Ds_PVnoDs_ipchi2_vec);

    tree->Branch("best_PV_withDs_IsValid", &best_PV_withDs_IsValid_vec);
    tree->Branch("best_PV_withDs_IsFake", &best_PV_withDs_IsFake_vec);
    tree->Branch("best_PV_withDs_chi2", &best_PV_withDs_chi2_vec);
    tree->Branch("best_PV_withDs_ndof", &best_PV_withDs_ndof_vec);
    tree->Branch("best_PV_withDs_chi2ndof", &best_PV_withDs_chi2ndof_vec);
    tree->Branch("best_PV_withDs_x", &best_PV_withDs_x_vec);
    tree->Branch("best_PV_withDs_y", &best_PV_withDs_y_vec);
    tree->Branch("best_PV_withDs_z", &best_PV_withDs_z_vec);
    tree->Branch("best_PV_withDs_xerr", &best_PV_withDs_xerr_vec);
    tree->Branch("best_PV_withDs_yerr", &best_PV_withDs_yerr_vec);
    tree->Branch("best_PV_withDs_zerr", &best_PV_withDs_zerr_vec);

    tree->Branch("best_Ds_PVwithDs_FDxy", &best_Ds_PVwithDs_FDxy_vec);
    tree->Branch("best_Ds_PVwithDs_FDz", &best_Ds_PVwithDs_FDz_vec);
    tree->Branch("best_Ds_PVwithDs_FD", &best_Ds_PVwithDs_FD_vec);
    tree->Branch("best_Ds_PVwithDs_FDxyerr", &best_Ds_PVwithDs_FDxyerr_vec);
    tree->Branch("best_Ds_PVwithDs_FDxychi2", &best_Ds_PVwithDs_FDxychi2_vec);
    tree->Branch("best_Ds_PVwithDs_FDzerr", &best_Ds_PVwithDs_FDzerr_vec);
    tree->Branch("best_Ds_PVwithDs_FDzchi2", &best_Ds_PVwithDs_FDzchi2_vec);
    tree->Branch("best_Ds_PVwithDs_FDerr", &best_Ds_PVwithDs_FDerr_vec);
    tree->Branch("best_Ds_PVwithDs_FDchi2", &best_Ds_PVwithDs_FDchi2_vec);
    tree->Branch("best_Ds_PVwithDs_dira", &best_Ds_PVwithDs_dira_vec);
    tree->Branch("best_Ds_PVwithDs_dira_angle", &best_Ds_PVwithDs_dira_angle_vec);
    tree->Branch("best_Kp_PVwithDs_ip", &best_Kp_PVwithDs_ip_vec);
    tree->Branch("best_Kp_PVwithDs_iperr", &best_Kp_PVwithDs_iperr_vec);
    tree->Branch("best_Kp_PVwithDs_ipchi2", &best_Kp_PVwithDs_ipchi2_vec);
    tree->Branch("best_Km_PVwithDs_ip", &best_Km_PVwithDs_ip_vec);
    tree->Branch("best_Km_PVwithDs_iperr", &best_Km_PVwithDs_iperr_vec);
    tree->Branch("best_Km_PVwithDs_ipchi2", &best_Km_PVwithDs_ipchi2_vec);
    tree->Branch("best_pi_PVwithDs_ip", &best_pi_PVwithDs_ip_vec);
    tree->Branch("best_pi_PVwithDs_iperr", &best_pi_PVwithDs_iperr_vec);
    tree->Branch("best_pi_PVwithDs_ipchi2", &best_pi_PVwithDs_ipchi2_vec);
    tree->Branch("best_phi_PVwithDs_ip", &best_phi_PVwithDs_ip_vec);
    tree->Branch("best_phi_PVwithDs_iperr", &best_phi_PVwithDs_iperr_vec);
    tree->Branch("best_phi_PVwithDs_ipchi2", &best_phi_PVwithDs_ipchi2_vec);
    tree->Branch("best_Ds_PVwithDs_ip", &best_Ds_PVwithDs_ip_vec);
    tree->Branch("best_Ds_PVwithDs_iperr", &best_Ds_PVwithDs_iperr_vec);
    tree->Branch("best_Ds_PVwithDs_ipchi2", &best_Ds_PVwithDs_ipchi2_vec);

    tree->Branch("best_match_entry", &best_match_entry_vec);

    tree->Branch("best_mu_charge", &best_mu_charge_vec);
    tree->Branch("best_mu_eta", &best_mu_eta_vec);
    tree->Branch("best_mu_phi", &best_mu_phi_vec);
    tree->Branch("best_mu_vx", &best_mu_vx_vec);
    tree->Branch("best_mu_vy", &best_mu_vy_vec);
    tree->Branch("best_mu_vz", &best_mu_vz_vec);
    tree->Branch("best_mu_p", &best_mu_p_vec);
    tree->Branch("best_mu_pt", &best_mu_pt_vec);
    tree->Branch("best_mu_px", &best_mu_px_vec);
    tree->Branch("best_mu_py", &best_mu_py_vec);
    tree->Branch("best_mu_pz", &best_mu_pz_vec);
    tree->Branch("best_mu_isHighPt", &best_mu_isHighPt_vec);
    tree->Branch("best_mu_isLoose", &best_mu_isLoose_vec);
    tree->Branch("best_mu_isMedium", &best_mu_isMedium_vec);
    tree->Branch("best_mu_isSoft", &best_mu_isSoft_vec);
    tree->Branch("best_mu_isTight", &best_mu_isTight_vec);
    tree->Branch("best_mu_isPF", &best_mu_isPF_vec);
    tree->Branch("best_mu_isTracker", &best_mu_isTracker_vec);
    tree->Branch("best_mu_isGlobal", &best_mu_isGlobal_vec);
    tree->Branch("best_mu_IsoR03_sumChargedHadronPt", &best_mu_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("best_mu_IsoR03_sumChargedParticlePt", &best_mu_IsoR03_sumChargedParticlePt_vec);
    tree->Branch("best_mu_IsoR03_sumNeutralHadronEt", &best_mu_IsoR03_sumNeutralHadronEt_vec);
    tree->Branch("best_mu_IsoR03_sumPhotonEt", &best_mu_IsoR03_sumPhotonEt_vec);
    tree->Branch("best_mu_IsoR03_sumPUPt", &best_mu_IsoR03_sumPUPt_vec);
    tree->Branch("best_mu_PFIsoR03", &best_mu_PFIsoR03_vec);
    tree->Branch("best_mu_IsoR04_sumChargedHadronPt", &best_mu_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("best_mu_IsoR04_sumChargedParticlePt", &best_mu_IsoR04_sumChargedParticlePt_vec);
    tree->Branch("best_mu_IsoR04_sumNeutralHadronEt", &best_mu_IsoR04_sumNeutralHadronEt_vec);
    tree->Branch("best_mu_IsoR04_sumPhotonEt", &best_mu_IsoR04_sumPhotonEt_vec);
    tree->Branch("best_mu_IsoR04_sumPUPt", &best_mu_IsoR04_sumPUPt_vec);
    tree->Branch("best_mu_PFIsoR04", &best_mu_PFIsoR04_vec);
    tree->Branch("best_mu_primvtx_dxy", &best_mu_primvtx_dxy_vec);
    tree->Branch("best_mu_primvtx_dxyerr", &best_mu_primvtx_dxyerr_vec);
    tree->Branch("best_mu_primvtx_dz", &best_mu_primvtx_dz_vec);
    tree->Branch("best_mu_primvtx_dzerr", &best_mu_primvtx_dzerr_vec);
    tree->Branch("best_mu_primvtx_ip", &best_mu_primvtx_ip_vec);
    tree->Branch("best_mu_primvtx_iperr", &best_mu_primvtx_iperr_vec);
    tree->Branch("best_mu_primvtx_ipchi2", &best_mu_primvtx_ipchi2_vec);
    tree->Branch("best_mu_match", &best_mu_match_vec);
}

void PVTree::Gen_Reset()
{
    Gen_Kp_eta = null;
    Gen_Kp_phi = null;
    Gen_Kp_vx = null;
    Gen_Kp_vy = null;
    Gen_Kp_vz = null;
    Gen_Kp_p = null;
    Gen_Kp_pt = null;
    Gen_Kp_px = null;
    Gen_Kp_py = null;
    Gen_Kp_pz = null;
    Gen_Kp_pp = null;
    Gen_Kp_pl = null;

    Gen_Km_eta = null;
    Gen_Km_phi = null;
    Gen_Km_vx = null;
    Gen_Km_vy = null;
    Gen_Km_vz = null;
    Gen_Km_p = null;
    Gen_Km_pt = null;
    Gen_Km_px = null;
    Gen_Km_py = null;
    Gen_Km_pz = null;
    Gen_Km_pp = null;
    Gen_Km_pl = null;

    Gen_pi_eta = null;
    Gen_pi_phi = null;
    Gen_pi_vx = null;
    Gen_pi_vy = null;
    Gen_pi_vz = null;
    Gen_pi_p = null;
    Gen_pi_pt = null;
    Gen_pi_px = null;
    Gen_pi_py = null;
    Gen_pi_pz = null;
    Gen_pi_pp = null;
    Gen_pi_pl = null;

    Gen_phi_eta = null;
    Gen_phi_phi = null;
    Gen_phi_vx = null;
    Gen_phi_vy = null;
    Gen_phi_vz = null;
    Gen_phi_p = null;
    Gen_phi_pt = null;
    Gen_phi_px = null;
    Gen_phi_py = null;
    Gen_phi_pz = null;
    Gen_phi_pp = null;
    Gen_phi_pl = null;

    Gen_Ds_eta = null;
    Gen_Ds_phi = null;
    Gen_Ds_vx = null;
    Gen_Ds_vy = null;
    Gen_Ds_vz = null;
    Gen_Ds_p = null;
    Gen_Ds_pt = null;
    Gen_Ds_px = null;
    Gen_Ds_py = null;
    Gen_Ds_pz = null;

    Gen_mu_eta = null;
    Gen_mu_phi = null;
    Gen_mu_vx = null;
    Gen_mu_vy = null;
    Gen_mu_vz = null;
    Gen_mu_p = null;
    Gen_mu_pt = null;
    Gen_mu_px = null;
    Gen_mu_py = null;
    Gen_mu_pz = null;

    Gen_nu_eta = null;
    Gen_nu_phi = null;
    Gen_nu_vx = null;
    Gen_nu_vy = null;
    Gen_nu_vz = null;
    Gen_nu_p = null;
    Gen_nu_pt = null;
    Gen_nu_px = null;
    Gen_nu_py = null;
    Gen_nu_pz = null;

    Gen_W_eta = null;
    Gen_W_phi = null;
    Gen_W_vx = null;
    Gen_W_vy = null;
    Gen_W_vz = null;
    Gen_W_p = null;
    Gen_W_pt = null;
    Gen_W_px = null;
    Gen_W_py = null;
    Gen_W_pz = null;

    Gen_H_eta = null;
    Gen_H_phi = null;
    Gen_H_vx = null;
    Gen_H_vy = null;
    Gen_H_vz = null;
    Gen_H_p = null;
    Gen_H_pt = null;
    Gen_H_px = null;
    Gen_H_py = null;
    Gen_H_pz = null;

    Gen_dR_Kp_Km = null;
    Gen_dR_Kp_phi = null;
    Gen_dR_Km_phi = null;
    Gen_dR_Kp_pi = null;
    Gen_dR_Km_pi = null;
    Gen_dR_pi_phi = null;
    Gen_dR_Kp_Ds = null;
    Gen_dR_Km_Ds = null;
    Gen_dR_phi_Ds = null;
    Gen_dR_pi_Ds = null;
    Gen_dR_Kp_mu = null;
    Gen_dR_Km_mu = null;
    Gen_dR_phi_mu = null;
    Gen_dR_pi_mu = null;
    Gen_dR_Ds_mu = null;

    Gen_Ds_dx = null;
    Gen_Ds_dy = null;
    Gen_Ds_dz = null;
    Gen_Ds_FDxy = null;
    Gen_Ds_FD = null;
}

void PVTree::BS_Reset()
{
    BS_type = null; 
    BS_x0 = null;
    BS_y0 = null;
    BS_z0 = null;
    BS_sigmaZ = null;
    BS_dxdz = null;
    BS_dydz = null;
    BS_BWX = null;
    BS_BWY = null;
    BS_x0err = null;
    BS_y0err = null;
    BS_z0err = null;
    BS_sigmaZ0err = null;
    BS_dxdzerr = null;
    BS_dydzerr = null;
    BS_BWXerr = null;
    BS_BWYerr = null;
    BS_emitX = null;
    BS_emitY = null;
    BS_betaStar = null;
}

void PVTree::PV_Reset()
{
    PV_vx = null;
    PV_vy = null;
    PV_vz = null;
    PV_vxerr = null;
    PV_vyerr = null;
    PV_vzerr = null;
}

void PVTree::Match_Reset()
{
    // Original info
    match_Kp_isIsolatedChargedHadron = false;
    match_Kp_charge = null;
    match_Kp_eta = null;
    match_Kp_phi = null;
    match_Kp_vx = null;
    match_Kp_vy = null;
    match_Kp_vz = null;
    match_Kp_p = null;
    match_Kp_pt = null;
    match_Kp_px = null;
    match_Kp_py = null;
    match_Kp_pz = null;

    match_Km_isIsolatedChargedHadron = false;
    match_Km_charge = null;
    match_Km_eta = null;
    match_Km_phi = null;
    match_Km_vx = null;
    match_Km_vy = null;
    match_Km_vz = null;
    match_Km_p = null;
    match_Km_pt = null;
    match_Km_px = null;
    match_Km_py = null;
    match_Km_pz = null;

    match_pi_isIsolatedChargedHadron = false;
    match_pi_charge = null;
    match_pi_eta = null;
    match_pi_phi = null;
    match_pi_vx = null;
    match_pi_vy = null;
    match_pi_vz = null;
    match_pi_p = null;
    match_pi_pt = null;
    match_pi_px = null;
    match_pi_py = null;
    match_pi_pz = null;

    match_phi_eta = null;
    match_phi_phi = null;
    match_phi_p = null;
    match_phi_pt = null;
    match_phi_px = null;
    match_phi_py = null;
    match_phi_pz = null;
    match_phi_invm = null;

    match_Ds_eta = null;
    match_Ds_phi = null;
    match_Ds_p = null;
    match_Ds_pt = null;
    match_Ds_px = null;
    match_Ds_py = null;
    match_Ds_pz = null;
    match_Ds_invm = null;

    match_dR_Kp_Km = null;
    match_dR_Kp_phi = null;
    match_dR_Km_phi = null;
    match_dR_Kp_pi = null;
    match_dR_Km_pi = null;
    match_dR_pi_phi = null;
    match_dR_Kp_Ds = null;
    match_dR_Km_Ds = null;
    match_dR_phi_Ds = null;
    match_dR_pi_Ds = null;

    match_dxy_phi_Ds = null;
    match_dz_phi_Ds = null;

    // Fit on phi
    match_phiFit_chi2 = null;
    match_phiFit_ndof = null;
    match_phiFit_chi2ndof = null;
    match_phiFit_vx = null;
    match_phiFit_vy = null;
    match_phiFit_vz = null;
    match_phiFit_vxerr = null;
    match_phiFit_vyerr = null;
    match_phiFit_vzerr = null;

    match_phiFit_Kp_eta = null;
    match_phiFit_Kp_phi = null;
    match_phiFit_Kp_p = null;
    match_phiFit_Kp_pt = null;
    match_phiFit_Kp_px = null;
    match_phiFit_Kp_py = null;
    match_phiFit_Kp_pz = null;

    match_phiFit_Km_eta = null;
    match_phiFit_Km_phi = null;
    match_phiFit_Km_p = null;
    match_phiFit_Km_pt = null;
    match_phiFit_Km_px = null;
    match_phiFit_Km_py = null;
    match_phiFit_Km_pz = null;

    match_phiFit_pi_eta = null;
    match_phiFit_pi_phi = null;
    match_phiFit_pi_p = null;
    match_phiFit_pi_pt = null;
    match_phiFit_pi_px = null;
    match_phiFit_pi_py = null;
    match_phiFit_pi_pz = null;

    match_phiFit_phi_eta = null;
    match_phiFit_phi_phi = null;
    match_phiFit_phi_p = null;
    match_phiFit_phi_pt = null;
    match_phiFit_phi_px = null;
    match_phiFit_phi_py = null;
    match_phiFit_phi_pz = null;
    match_phiFit_phi_invm = null;

    match_phiFit_Ds_eta = null;
    match_phiFit_Ds_phi = null;
    match_phiFit_Ds_p = null;
    match_phiFit_Ds_pt = null;
    match_phiFit_Ds_px = null;
    match_phiFit_Ds_py = null;
    match_phiFit_Ds_pz = null;
    match_phiFit_Ds_invm = null;

    match_phiFit_dR_Kp_Km = null;
    match_phiFit_dR_Kp_phi = null;
    match_phiFit_dR_Km_phi = null;
    match_phiFit_dR_Kp_pi = null;
    match_phiFit_dR_Km_pi = null;
    match_phiFit_dR_pi_phi = null;
    match_phiFit_dR_Kp_Ds = null;
    match_phiFit_dR_Km_Ds = null;
    match_phiFit_dR_phi_Ds = null;
    match_phiFit_dR_pi_Ds = null;

    // Fit on Ds 
    match_DsFit_chi2 = null;
    match_DsFit_ndof = null;
    match_DsFit_chi2ndof = null;
    match_DsFit_vx = null;
    match_DsFit_vy = null;
    match_DsFit_vz = null;
    match_DsFit_vxerr = null;
    match_DsFit_vyerr = null;
    match_DsFit_vzerr = null;

    match_DsFit_Kp_eta = null;
    match_DsFit_Kp_phi = null;
    match_DsFit_Kp_p = null;
    match_DsFit_Kp_pt = null;
    match_DsFit_Kp_px = null;
    match_DsFit_Kp_py = null;
    match_DsFit_Kp_pz = null;

    match_DsFit_Km_eta = null;
    match_DsFit_Km_phi = null;
    match_DsFit_Km_p = null;
    match_DsFit_Km_pt = null;
    match_DsFit_Km_px = null;
    match_DsFit_Km_py = null;
    match_DsFit_Km_pz = null;

    match_DsFit_pi_eta = null;
    match_DsFit_pi_phi = null;
    match_DsFit_pi_p = null;
    match_DsFit_pi_pt = null;
    match_DsFit_pi_px = null;
    match_DsFit_pi_py = null;
    match_DsFit_pi_pz = null;

    match_DsFit_phi_eta = null;
    match_DsFit_phi_phi = null;
    match_DsFit_phi_p = null;
    match_DsFit_phi_pt = null;
    match_DsFit_phi_px = null;
    match_DsFit_phi_py = null;
    match_DsFit_phi_pz = null;
    match_DsFit_phi_invm = null;

    match_DsFit_Ds_eta = null;
    match_DsFit_Ds_phi = null;
    match_DsFit_Ds_p = null;
    match_DsFit_Ds_pt = null;
    match_DsFit_Ds_px = null;
    match_DsFit_Ds_py = null;
    match_DsFit_Ds_pz = null;
    match_DsFit_Ds_invm = null;

    match_DsFit_dR_Kp_Km = null;
    match_DsFit_dR_Kp_phi = null;
    match_DsFit_dR_Km_phi = null;
    match_DsFit_dR_Kp_pi = null;
    match_DsFit_dR_Km_pi = null;
    match_DsFit_dR_pi_phi = null;
    match_DsFit_dR_Kp_Ds = null;
    match_DsFit_dR_Km_Ds = null;
    match_DsFit_dR_phi_Ds = null;
    match_DsFit_dR_pi_Ds = null;

    match_DsFit_Mconstraint_Ds_invm = null;

    match_Ds_primvtx_FDxy = null;
    match_Ds_primvtx_FDz = null;
    match_Ds_primvtx_FD = null;
    match_Ds_primvtx_FDxyerr = null;
    match_Ds_primvtx_FDxychi2 = null;
    match_Ds_primvtx_FDzerr = null;
    match_Ds_primvtx_FDzchi2 = null;
    match_Ds_primvtx_FDerr = null;
    match_Ds_primvtx_FDchi2 = null;
    match_Ds_primvtx_dira = null;
    match_Ds_primvtx_dira_angle = null;
    match_Kp_primvtx_ip = null;
    match_Kp_primvtx_iperr = null;
    match_Kp_primvtx_ipchi2 = null;
    match_Km_primvtx_ip = null;
    match_Km_primvtx_iperr = null;
    match_Km_primvtx_ipchi2 = null;
    match_pi_primvtx_ip = null;
    match_pi_primvtx_iperr = null;
    match_pi_primvtx_ipchi2 = null;
    match_phi_primvtx_ip = null;
    match_phi_primvtx_iperr = null;
    match_phi_primvtx_ipchi2 = null;
    match_Ds_primvtx_ip = null;
    match_Ds_primvtx_iperr = null;
    match_Ds_primvtx_ipchi2 = null;

    match_Ds_IsoR03_sumChargedHadronPt = 0.;
    match_Ds_IsoR03_sumNeutralHadronPt = 0.;
    match_Ds_IsoR03_sumPhotonPt = 0.;
    match_Ds_IsoR03_sumPUPt = 0.;
    match_Ds_PFIsoR03 = 0.;
    match_Ds_IsoR04_sumChargedHadronPt = 0.;
    match_Ds_IsoR04_sumNeutralHadronPt = 0.;
    match_Ds_IsoR04_sumPhotonPt = 0.;
    match_Ds_IsoR04_sumPUPt = 0.;
    match_Ds_PFIsoR04 = 0.;

    match_PV_noDs_IsValid = false;
    match_PV_noDs_IsFake = true;
    match_PV_noDs_chi2 = null;
    match_PV_noDs_ndof = null;
    match_PV_noDs_chi2ndof = null;
    match_PV_noDs_x = null;
    match_PV_noDs_y = null;
    match_PV_noDs_z = null;
    match_PV_noDs_xerr = null;
    match_PV_noDs_yerr = null;
    match_PV_noDs_zerr = null;

    match_Ds_PVnoDs_FDxy = null;
    match_Ds_PVnoDs_FDz = null;
    match_Ds_PVnoDs_FD = null;
    match_Ds_PVnoDs_FDxyerr = null;
    match_Ds_PVnoDs_FDxychi2 = null;
    match_Ds_PVnoDs_FDzerr = null;
    match_Ds_PVnoDs_FDzchi2 = null;
    match_Ds_PVnoDs_FDerr = null;
    match_Ds_PVnoDs_FDchi2 = null;
    match_Ds_PVnoDs_dira = null;
    match_Ds_PVnoDs_dira_angle = null;
    match_Kp_PVnoDs_ip = null;
    match_Kp_PVnoDs_iperr = null;
    match_Kp_PVnoDs_ipchi2 = null;
    match_Km_PVnoDs_ip = null;
    match_Km_PVnoDs_iperr = null;
    match_Km_PVnoDs_ipchi2 = null;
    match_pi_PVnoDs_ip = null;
    match_pi_PVnoDs_iperr = null;
    match_pi_PVnoDs_ipchi2 = null;
    match_phi_PVnoDs_ip = null;
    match_phi_PVnoDs_iperr = null;
    match_phi_PVnoDs_ipchi2 = null;
    match_Ds_PVnoDs_ip = null;
    match_Ds_PVnoDs_iperr = null;
    match_Ds_PVnoDs_ipchi2 = null;

    match_PV_withDs_IsValid = false;
    match_PV_withDs_IsFake = true;
    match_PV_withDs_chi2 = null;
    match_PV_withDs_ndof = null;
    match_PV_withDs_chi2ndof = null;
    match_PV_withDs_x = null;
    match_PV_withDs_y = null;
    match_PV_withDs_z = null;
    match_PV_withDs_xerr = null;
    match_PV_withDs_yerr = null;
    match_PV_withDs_zerr = null;

    match_Ds_PVwithDs_FDxy = null;
    match_Ds_PVwithDs_FDz = null;
    match_Ds_PVwithDs_FD = null;
    match_Ds_PVwithDs_FDxyerr = null;
    match_Ds_PVwithDs_FDxychi2 = null;
    match_Ds_PVwithDs_FDzerr = null;
    match_Ds_PVwithDs_FDzchi2 = null;
    match_Ds_PVwithDs_FDerr = null;
    match_Ds_PVwithDs_FDchi2 = null;
    match_Ds_PVwithDs_dira = null;
    match_Ds_PVwithDs_dira_angle = null;
    match_Kp_PVwithDs_ip = null;
    match_Kp_PVwithDs_iperr = null;
    match_Kp_PVwithDs_ipchi2 = null;
    match_Km_PVwithDs_ip = null;
    match_Km_PVwithDs_iperr = null;
    match_Km_PVwithDs_ipchi2 = null;
    match_pi_PVwithDs_ip = null;
    match_pi_PVwithDs_iperr = null;
    match_pi_PVwithDs_ipchi2 = null;
    match_phi_PVwithDs_ip = null;
    match_phi_PVwithDs_iperr = null;
    match_phi_PVwithDs_ipchi2 = null;
    match_Ds_PVwithDs_ip = null;
    match_Ds_PVwithDs_iperr = null;
    match_Ds_PVwithDs_ipchi2 = null;

    match_mu_charge = null;
    match_mu_eta = null;
    match_mu_phi = null;
    match_mu_vx = null;
    match_mu_vy = null;
    match_mu_vz = null;
    match_mu_p = null;
    match_mu_pt = null;
    match_mu_px = null;
    match_mu_py = null;
    match_mu_pz = null;
    match_mu_isHighPt = false;
    match_mu_isLoose = false;
    match_mu_isMedium = false;
    match_mu_isSoft = false;
    match_mu_isTight = false;
    match_mu_isPF = false;
    match_mu_isTracker = false;
    match_mu_isGlobal = false;
    match_mu_IsoR03_sumChargedHadronPt = null;
    match_mu_IsoR03_sumChargedParticlePt = null;
    match_mu_IsoR03_sumNeutralHadronEt = null;
    match_mu_IsoR03_sumPhotonEt = null;
    match_mu_IsoR03_sumPUPt = null;
    match_mu_PFIsoR03 = null;
    match_mu_IsoR04_sumChargedHadronPt = null;
    match_mu_IsoR04_sumChargedParticlePt = null;
    match_mu_IsoR04_sumNeutralHadronEt = null;
    match_mu_IsoR04_sumPhotonEt = null;
    match_mu_IsoR04_sumPUPt = null;
    match_mu_PFIsoR04 = null;
    match_mu_primvtx_dxy = null;
    match_mu_primvtx_dxyerr = null;
    match_mu_primvtx_dz = null;
    match_mu_primvtx_dzerr = null;
    match_mu_primvtx_ip = null;
    match_mu_primvtx_iperr = null;
    match_mu_primvtx_ipchi2 = null;
}

void PVTree::Match_Ds_Fill_Vector()
{
    match_Kp_isIsolatedChargedHadron_vec.push_back(match_Kp_isIsolatedChargedHadron);
    match_Kp_charge_vec.push_back(match_Kp_charge);
    match_Kp_eta_vec.push_back(match_Kp_eta);
    match_Kp_phi_vec.push_back(match_Kp_phi);
    match_Kp_vx_vec.push_back(match_Kp_vx);
    match_Kp_vy_vec.push_back(match_Kp_vy);
    match_Kp_vz_vec.push_back(match_Kp_vz);
    match_Kp_p_vec.push_back(match_Kp_p);
    match_Kp_pt_vec.push_back(match_Kp_pt);
    match_Kp_px_vec.push_back(match_Kp_px);
    match_Kp_py_vec.push_back(match_Kp_py);
    match_Kp_pz_vec.push_back(match_Kp_pz);

    match_Km_isIsolatedChargedHadron_vec.push_back(match_Km_isIsolatedChargedHadron);
    match_Km_charge_vec.push_back(match_Km_charge);
    match_Km_eta_vec.push_back(match_Km_eta);
    match_Km_phi_vec.push_back(match_Km_phi);
    match_Km_vx_vec.push_back(match_Km_vx);
    match_Km_vy_vec.push_back(match_Km_vy);
    match_Km_vz_vec.push_back(match_Km_vz);
    match_Km_p_vec.push_back(match_Km_p);
    match_Km_pt_vec.push_back(match_Km_pt);
    match_Km_px_vec.push_back(match_Km_px);
    match_Km_py_vec.push_back(match_Km_py);
    match_Km_pz_vec.push_back(match_Km_pz);

    match_pi_isIsolatedChargedHadron_vec.push_back(match_pi_isIsolatedChargedHadron);
    match_pi_charge_vec.push_back(match_pi_charge);
    match_pi_eta_vec.push_back(match_pi_eta);
    match_pi_phi_vec.push_back(match_pi_phi);
    match_pi_vx_vec.push_back(match_pi_vx);
    match_pi_vy_vec.push_back(match_pi_vy);
    match_pi_vz_vec.push_back(match_pi_vz);
    match_pi_p_vec.push_back(match_pi_p);
    match_pi_pt_vec.push_back(match_pi_pt);
    match_pi_px_vec.push_back(match_pi_px);
    match_pi_py_vec.push_back(match_pi_py);
    match_pi_pz_vec.push_back(match_pi_pz);

    match_phi_eta_vec.push_back(match_phi_eta);
    match_phi_phi_vec.push_back(match_phi_phi);
    match_phi_p_vec.push_back(match_phi_p);
    match_phi_pt_vec.push_back(match_phi_pt);
    match_phi_px_vec.push_back(match_phi_px);
    match_phi_py_vec.push_back(match_phi_py);
    match_phi_pz_vec.push_back(match_phi_pz);
    match_phi_invm_vec.push_back(match_phi_invm);

    match_Ds_eta_vec.push_back(match_Ds_eta);
    match_Ds_phi_vec.push_back(match_Ds_phi);
    match_Ds_p_vec.push_back(match_Ds_p);
    match_Ds_pt_vec.push_back(match_Ds_pt);
    match_Ds_px_vec.push_back(match_Ds_px);
    match_Ds_py_vec.push_back(match_Ds_py);
    match_Ds_pz_vec.push_back(match_Ds_pz);
    match_Ds_invm_vec.push_back(match_Ds_invm);

    match_dR_Kp_Km_vec.push_back(match_dR_Kp_Km);
    match_dR_Kp_phi_vec.push_back(match_dR_Kp_phi);
    match_dR_Km_phi_vec.push_back(match_dR_Km_phi);
    match_dR_Kp_pi_vec.push_back(match_dR_Kp_pi);
    match_dR_Km_pi_vec.push_back(match_dR_Km_pi);
    match_dR_pi_phi_vec.push_back(match_dR_pi_phi);
    match_dR_Kp_Ds_vec.push_back(match_dR_Kp_Ds);
    match_dR_Km_Ds_vec.push_back(match_dR_Km_Ds);
    match_dR_phi_Ds_vec.push_back(match_dR_phi_Ds);
    match_dR_pi_Ds_vec.push_back(match_dR_pi_Ds);

    match_dxy_phi_Ds_vec.push_back(match_dxy_phi_Ds);
    match_dz_phi_Ds_vec.push_back(match_dz_phi_Ds);

    // Fit on phi
    match_phiFit_chi2_vec.push_back(match_phiFit_chi2);
    match_phiFit_ndof_vec.push_back(match_phiFit_ndof);
    match_phiFit_chi2ndof_vec.push_back(match_phiFit_chi2ndof);
    match_phiFit_vx_vec.push_back(match_phiFit_vx);
    match_phiFit_vy_vec.push_back(match_phiFit_vy);
    match_phiFit_vz_vec.push_back(match_phiFit_vz);
    match_phiFit_vxerr_vec.push_back(match_phiFit_vxerr);
    match_phiFit_vyerr_vec.push_back(match_phiFit_vyerr);
    match_phiFit_vzerr_vec.push_back(match_phiFit_vzerr);

    match_phiFit_Kp_eta_vec.push_back(match_phiFit_Kp_eta);
    match_phiFit_Kp_phi_vec.push_back(match_phiFit_Kp_phi);
    match_phiFit_Kp_p_vec.push_back(match_phiFit_Kp_p);
    match_phiFit_Kp_pt_vec.push_back(match_phiFit_Kp_pt);
    match_phiFit_Kp_px_vec.push_back(match_phiFit_Kp_px);
    match_phiFit_Kp_py_vec.push_back(match_phiFit_Kp_py);
    match_phiFit_Kp_pz_vec.push_back(match_phiFit_Kp_pz);

    match_phiFit_Km_eta_vec.push_back(match_phiFit_Km_eta);
    match_phiFit_Km_phi_vec.push_back(match_phiFit_Km_phi);
    match_phiFit_Km_p_vec.push_back(match_phiFit_Km_p);
    match_phiFit_Km_pt_vec.push_back(match_phiFit_Km_pt);
    match_phiFit_Km_px_vec.push_back(match_phiFit_Km_px);
    match_phiFit_Km_py_vec.push_back(match_phiFit_Km_py);
    match_phiFit_Km_pz_vec.push_back(match_phiFit_Km_pz);

    match_phiFit_pi_eta_vec.push_back(match_phiFit_pi_eta);
    match_phiFit_pi_phi_vec.push_back(match_phiFit_pi_phi);
    match_phiFit_pi_p_vec.push_back(match_phiFit_pi_p);
    match_phiFit_pi_pt_vec.push_back(match_phiFit_pi_pt);
    match_phiFit_pi_px_vec.push_back(match_phiFit_pi_px);
    match_phiFit_pi_py_vec.push_back(match_phiFit_pi_py);
    match_phiFit_pi_pz_vec.push_back(match_phiFit_pi_pz);

    match_phiFit_phi_eta_vec.push_back(match_phiFit_phi_eta);
    match_phiFit_phi_phi_vec.push_back(match_phiFit_phi_phi);
    match_phiFit_phi_p_vec.push_back(match_phiFit_phi_p);
    match_phiFit_phi_pt_vec.push_back(match_phiFit_phi_pt);
    match_phiFit_phi_px_vec.push_back(match_phiFit_phi_px);
    match_phiFit_phi_py_vec.push_back(match_phiFit_phi_py);
    match_phiFit_phi_pz_vec.push_back(match_phiFit_phi_pz);
    match_phiFit_phi_invm_vec.push_back(match_phiFit_phi_invm);

    match_phiFit_Ds_eta_vec.push_back(match_phiFit_Ds_eta);
    match_phiFit_Ds_phi_vec.push_back(match_phiFit_Ds_phi);
    match_phiFit_Ds_p_vec.push_back(match_phiFit_Ds_p);
    match_phiFit_Ds_pt_vec.push_back(match_phiFit_Ds_pt);
    match_phiFit_Ds_px_vec.push_back(match_phiFit_Ds_px);
    match_phiFit_Ds_py_vec.push_back(match_phiFit_Ds_py);
    match_phiFit_Ds_pz_vec.push_back(match_phiFit_Ds_pz);
    match_phiFit_Ds_invm_vec.push_back(match_phiFit_Ds_invm);

    match_phiFit_dR_Kp_Km_vec.push_back(match_phiFit_dR_Kp_Km);
    match_phiFit_dR_Kp_phi_vec.push_back(match_phiFit_dR_Kp_phi);
    match_phiFit_dR_Km_phi_vec.push_back(match_phiFit_dR_Km_phi);
    match_phiFit_dR_Kp_pi_vec.push_back(match_phiFit_dR_Kp_pi);
    match_phiFit_dR_Km_pi_vec.push_back(match_phiFit_dR_Km_pi);
    match_phiFit_dR_pi_phi_vec.push_back(match_phiFit_dR_pi_phi);
    match_phiFit_dR_Kp_Ds_vec.push_back(match_phiFit_dR_Kp_Ds);
    match_phiFit_dR_Km_Ds_vec.push_back(match_phiFit_dR_Km_Ds);
    match_phiFit_dR_phi_Ds_vec.push_back(match_phiFit_dR_phi_Ds);
    match_phiFit_dR_pi_Ds_vec.push_back(match_phiFit_dR_pi_Ds);

    // Fit on Ds 
    match_DsFit_chi2_vec.push_back(match_DsFit_chi2);
    match_DsFit_ndof_vec.push_back(match_DsFit_ndof);
    match_DsFit_chi2ndof_vec.push_back(match_DsFit_chi2ndof);
    match_DsFit_vx_vec.push_back(match_DsFit_vx);
    match_DsFit_vy_vec.push_back(match_DsFit_vy);
    match_DsFit_vz_vec.push_back(match_DsFit_vz);
    match_DsFit_vxerr_vec.push_back(match_DsFit_vxerr);
    match_DsFit_vyerr_vec.push_back(match_DsFit_vyerr);
    match_DsFit_vzerr_vec.push_back(match_DsFit_vzerr);

    match_DsFit_Kp_eta_vec.push_back(match_DsFit_Kp_eta);
    match_DsFit_Kp_phi_vec.push_back(match_DsFit_Kp_phi);
    match_DsFit_Kp_p_vec.push_back(match_DsFit_Kp_p);
    match_DsFit_Kp_pt_vec.push_back(match_DsFit_Kp_pt);
    match_DsFit_Kp_px_vec.push_back(match_DsFit_Kp_px);
    match_DsFit_Kp_py_vec.push_back(match_DsFit_Kp_py);
    match_DsFit_Kp_pz_vec.push_back(match_DsFit_Kp_pz);

    match_DsFit_Km_eta_vec.push_back(match_DsFit_Km_eta);
    match_DsFit_Km_phi_vec.push_back(match_DsFit_Km_phi);
    match_DsFit_Km_p_vec.push_back(match_DsFit_Km_p);
    match_DsFit_Km_pt_vec.push_back(match_DsFit_Km_pt);
    match_DsFit_Km_px_vec.push_back(match_DsFit_Km_px);
    match_DsFit_Km_py_vec.push_back(match_DsFit_Km_py);
    match_DsFit_Km_pz_vec.push_back(match_DsFit_Km_pz);

    match_DsFit_pi_eta_vec.push_back(match_DsFit_pi_eta);
    match_DsFit_pi_phi_vec.push_back(match_DsFit_pi_phi);
    match_DsFit_pi_p_vec.push_back(match_DsFit_pi_p);
    match_DsFit_pi_pt_vec.push_back(match_DsFit_pi_pt);
    match_DsFit_pi_px_vec.push_back(match_DsFit_pi_px);
    match_DsFit_pi_py_vec.push_back(match_DsFit_pi_py);
    match_DsFit_pi_pz_vec.push_back(match_DsFit_pi_pz);

    match_DsFit_phi_eta_vec.push_back(match_DsFit_phi_eta);
    match_DsFit_phi_phi_vec.push_back(match_DsFit_phi_phi);
    match_DsFit_phi_p_vec.push_back(match_DsFit_phi_p);
    match_DsFit_phi_pt_vec.push_back(match_DsFit_phi_pt);
    match_DsFit_phi_px_vec.push_back(match_DsFit_phi_px);
    match_DsFit_phi_py_vec.push_back(match_DsFit_phi_py);
    match_DsFit_phi_pz_vec.push_back(match_DsFit_phi_pz);
    match_DsFit_phi_invm_vec.push_back(match_DsFit_phi_invm);

    match_DsFit_Ds_eta_vec.push_back(match_DsFit_Ds_eta);
    match_DsFit_Ds_phi_vec.push_back(match_DsFit_Ds_phi);
    match_DsFit_Ds_p_vec.push_back(match_DsFit_Ds_p);
    match_DsFit_Ds_pt_vec.push_back(match_DsFit_Ds_pt);
    match_DsFit_Ds_px_vec.push_back(match_DsFit_Ds_px);
    match_DsFit_Ds_py_vec.push_back(match_DsFit_Ds_py);
    match_DsFit_Ds_pz_vec.push_back(match_DsFit_Ds_pz);
    match_DsFit_Ds_invm_vec.push_back(match_DsFit_Ds_invm);

    match_DsFit_dR_Kp_Km_vec.push_back(match_DsFit_dR_Kp_Km);
    match_DsFit_dR_Kp_phi_vec.push_back(match_DsFit_dR_Kp_phi);
    match_DsFit_dR_Km_phi_vec.push_back(match_DsFit_dR_Km_phi);
    match_DsFit_dR_Kp_pi_vec.push_back(match_DsFit_dR_Kp_pi);
    match_DsFit_dR_Km_pi_vec.push_back(match_DsFit_dR_Km_pi);
    match_DsFit_dR_pi_phi_vec.push_back(match_DsFit_dR_pi_phi);
    match_DsFit_dR_Kp_Ds_vec.push_back(match_DsFit_dR_Kp_Ds);
    match_DsFit_dR_Km_Ds_vec.push_back(match_DsFit_dR_Km_Ds);
    match_DsFit_dR_phi_Ds_vec.push_back(match_DsFit_dR_phi_Ds);
    match_DsFit_dR_pi_Ds_vec.push_back(match_DsFit_dR_pi_Ds);

    match_DsFit_Mconstraint_Ds_invm_vec.push_back(match_DsFit_Mconstraint_Ds_invm);

    match_Ds_primvtx_FDxy_vec.push_back(match_Ds_primvtx_FDxy);
    match_Ds_primvtx_FDz_vec.push_back(match_Ds_primvtx_FDz);
    match_Ds_primvtx_FD_vec.push_back(match_Ds_primvtx_FD);
    match_Ds_primvtx_FDxyerr_vec.push_back(match_Ds_primvtx_FDxyerr);
    match_Ds_primvtx_FDxychi2_vec.push_back(match_Ds_primvtx_FDxychi2);
    match_Ds_primvtx_FDzerr_vec.push_back(match_Ds_primvtx_FDzerr);
    match_Ds_primvtx_FDzchi2_vec.push_back(match_Ds_primvtx_FDzchi2);
    match_Ds_primvtx_FDerr_vec.push_back(match_Ds_primvtx_FDerr);
    match_Ds_primvtx_FDchi2_vec.push_back(match_Ds_primvtx_FDchi2);
    match_Ds_primvtx_dira_vec.push_back(match_Ds_primvtx_dira);
    match_Ds_primvtx_dira_angle_vec.push_back(match_Ds_primvtx_dira_angle);
    match_Kp_primvtx_ip_vec.push_back(match_Kp_primvtx_ip);
    match_Kp_primvtx_iperr_vec.push_back(match_Kp_primvtx_iperr);
    match_Kp_primvtx_ipchi2_vec.push_back(match_Kp_primvtx_ipchi2);
    match_Km_primvtx_ip_vec.push_back(match_Km_primvtx_ip);
    match_Km_primvtx_iperr_vec.push_back(match_Km_primvtx_iperr);
    match_Km_primvtx_ipchi2_vec.push_back(match_Km_primvtx_ipchi2);
    match_pi_primvtx_ip_vec.push_back(match_pi_primvtx_ip);
    match_pi_primvtx_iperr_vec.push_back(match_pi_primvtx_iperr);
    match_pi_primvtx_ipchi2_vec.push_back(match_pi_primvtx_ipchi2);
    match_phi_primvtx_ip_vec.push_back(match_phi_primvtx_ip);
    match_phi_primvtx_iperr_vec.push_back(match_phi_primvtx_iperr);
    match_phi_primvtx_ipchi2_vec.push_back(match_phi_primvtx_ipchi2);
    match_Ds_primvtx_ip_vec.push_back(match_Ds_primvtx_ip);
    match_Ds_primvtx_iperr_vec.push_back(match_Ds_primvtx_iperr);
    match_Ds_primvtx_ipchi2_vec.push_back(match_Ds_primvtx_ipchi2);

    match_Ds_IsoR03_sumChargedHadronPt_vec.push_back(match_Ds_IsoR03_sumChargedHadronPt);
    match_Ds_IsoR03_sumNeutralHadronPt_vec.push_back(match_Ds_IsoR03_sumNeutralHadronPt);
    match_Ds_IsoR03_sumPhotonPt_vec.push_back(match_Ds_IsoR03_sumPhotonPt);
    match_Ds_IsoR03_sumPUPt_vec.push_back(match_Ds_IsoR03_sumPUPt);
    match_Ds_PFIsoR03_vec.push_back(match_Ds_PFIsoR03);
    match_Ds_IsoR04_sumChargedHadronPt_vec.push_back(match_Ds_IsoR04_sumChargedHadronPt);
    match_Ds_IsoR04_sumNeutralHadronPt_vec.push_back(match_Ds_IsoR04_sumNeutralHadronPt);
    match_Ds_IsoR04_sumPhotonPt_vec.push_back(match_Ds_IsoR04_sumPhotonPt);
    match_Ds_IsoR04_sumPUPt_vec.push_back(match_Ds_IsoR04_sumPUPt);
    match_Ds_PFIsoR04_vec.push_back(match_Ds_PFIsoR04);

    match_PV_noDs_IsValid_vec.push_back(match_PV_noDs_IsValid);
    match_PV_noDs_IsFake_vec.push_back(match_PV_noDs_IsFake);
    match_PV_noDs_chi2_vec.push_back(match_PV_noDs_chi2);
    match_PV_noDs_ndof_vec.push_back(match_PV_noDs_ndof);
    match_PV_noDs_chi2ndof_vec.push_back(match_PV_noDs_chi2ndof);
    match_PV_noDs_x_vec.push_back(match_PV_noDs_x);
    match_PV_noDs_y_vec.push_back(match_PV_noDs_y);
    match_PV_noDs_z_vec.push_back(match_PV_noDs_z);
    match_PV_noDs_xerr_vec.push_back(match_PV_noDs_xerr);
    match_PV_noDs_yerr_vec.push_back(match_PV_noDs_yerr);
    match_PV_noDs_zerr_vec.push_back(match_PV_noDs_zerr);

    match_Ds_PVnoDs_FDxy_vec.push_back(match_Ds_PVnoDs_FDxy);
    match_Ds_PVnoDs_FDz_vec.push_back(match_Ds_PVnoDs_FDz);
    match_Ds_PVnoDs_FD_vec.push_back(match_Ds_PVnoDs_FD);
    match_Ds_PVnoDs_FDxyerr_vec.push_back(match_Ds_PVnoDs_FDxyerr);
    match_Ds_PVnoDs_FDxychi2_vec.push_back(match_Ds_PVnoDs_FDxychi2);
    match_Ds_PVnoDs_FDzerr_vec.push_back(match_Ds_PVnoDs_FDzerr);
    match_Ds_PVnoDs_FDzchi2_vec.push_back(match_Ds_PVnoDs_FDzchi2);
    match_Ds_PVnoDs_FDerr_vec.push_back(match_Ds_PVnoDs_FDerr);
    match_Ds_PVnoDs_FDchi2_vec.push_back(match_Ds_PVnoDs_FDchi2);
    match_Ds_PVnoDs_dira_vec.push_back(match_Ds_PVnoDs_dira);
    match_Ds_PVnoDs_dira_angle_vec.push_back(match_Ds_PVnoDs_dira_angle);
    match_Kp_PVnoDs_ip_vec.push_back(match_Kp_PVnoDs_ip);
    match_Kp_PVnoDs_iperr_vec.push_back(match_Kp_PVnoDs_iperr);
    match_Kp_PVnoDs_ipchi2_vec.push_back(match_Kp_PVnoDs_ipchi2);
    match_Km_PVnoDs_ip_vec.push_back(match_Km_PVnoDs_ip);
    match_Km_PVnoDs_iperr_vec.push_back(match_Km_PVnoDs_iperr);
    match_Km_PVnoDs_ipchi2_vec.push_back(match_Km_PVnoDs_ipchi2);
    match_pi_PVnoDs_ip_vec.push_back(match_pi_PVnoDs_ip);
    match_pi_PVnoDs_iperr_vec.push_back(match_pi_PVnoDs_iperr);
    match_pi_PVnoDs_ipchi2_vec.push_back(match_pi_PVnoDs_ipchi2);
    match_phi_PVnoDs_ip_vec.push_back(match_phi_PVnoDs_ip);
    match_phi_PVnoDs_iperr_vec.push_back(match_phi_PVnoDs_iperr);
    match_phi_PVnoDs_ipchi2_vec.push_back(match_phi_PVnoDs_ipchi2);
    match_Ds_PVnoDs_ip_vec.push_back(match_Ds_PVnoDs_ip);
    match_Ds_PVnoDs_iperr_vec.push_back(match_Ds_PVnoDs_iperr);
    match_Ds_PVnoDs_ipchi2_vec.push_back(match_Ds_PVnoDs_ipchi2);

    match_PV_withDs_IsValid_vec.push_back(match_PV_withDs_IsValid);
    match_PV_withDs_IsFake_vec.push_back(match_PV_withDs_IsFake);
    match_PV_withDs_chi2_vec.push_back(match_PV_withDs_chi2);
    match_PV_withDs_ndof_vec.push_back(match_PV_withDs_ndof);
    match_PV_withDs_chi2ndof_vec.push_back(match_PV_withDs_chi2ndof);
    match_PV_withDs_x_vec.push_back(match_PV_withDs_x);
    match_PV_withDs_y_vec.push_back(match_PV_withDs_y);
    match_PV_withDs_z_vec.push_back(match_PV_withDs_z);
    match_PV_withDs_xerr_vec.push_back(match_PV_withDs_xerr);
    match_PV_withDs_yerr_vec.push_back(match_PV_withDs_yerr);
    match_PV_withDs_zerr_vec.push_back(match_PV_withDs_zerr);

    match_Ds_PVwithDs_FDxy_vec.push_back(match_Ds_PVwithDs_FDxy);
    match_Ds_PVwithDs_FDz_vec.push_back(match_Ds_PVwithDs_FDz);
    match_Ds_PVwithDs_FD_vec.push_back(match_Ds_PVwithDs_FD);
    match_Ds_PVwithDs_FDxyerr_vec.push_back(match_Ds_PVwithDs_FDxyerr);
    match_Ds_PVwithDs_FDxychi2_vec.push_back(match_Ds_PVwithDs_FDxychi2);
    match_Ds_PVwithDs_FDzerr_vec.push_back(match_Ds_PVwithDs_FDzerr);
    match_Ds_PVwithDs_FDzchi2_vec.push_back(match_Ds_PVwithDs_FDzchi2);
    match_Ds_PVwithDs_FDerr_vec.push_back(match_Ds_PVwithDs_FDerr);
    match_Ds_PVwithDs_FDchi2_vec.push_back(match_Ds_PVwithDs_FDchi2);
    match_Ds_PVwithDs_dira_vec.push_back(match_Ds_PVwithDs_dira);
    match_Ds_PVwithDs_dira_angle_vec.push_back(match_Ds_PVwithDs_dira_angle);
    match_Kp_PVwithDs_ip_vec.push_back(match_Kp_PVwithDs_ip);
    match_Kp_PVwithDs_iperr_vec.push_back(match_Kp_PVwithDs_iperr);
    match_Kp_PVwithDs_ipchi2_vec.push_back(match_Kp_PVwithDs_ipchi2);
    match_Km_PVwithDs_ip_vec.push_back(match_Km_PVwithDs_ip);
    match_Km_PVwithDs_iperr_vec.push_back(match_Km_PVwithDs_iperr);
    match_Km_PVwithDs_ipchi2_vec.push_back(match_Km_PVwithDs_ipchi2);
    match_pi_PVwithDs_ip_vec.push_back(match_pi_PVwithDs_ip);
    match_pi_PVwithDs_iperr_vec.push_back(match_pi_PVwithDs_iperr);
    match_pi_PVwithDs_ipchi2_vec.push_back(match_pi_PVwithDs_ipchi2);
    match_phi_PVwithDs_ip_vec.push_back(match_phi_PVwithDs_ip);
    match_phi_PVwithDs_iperr_vec.push_back(match_phi_PVwithDs_iperr);
    match_phi_PVwithDs_ipchi2_vec.push_back(match_phi_PVwithDs_ipchi2);
    match_Ds_PVwithDs_ip_vec.push_back(match_Ds_PVwithDs_ip);
    match_Ds_PVwithDs_iperr_vec.push_back(match_Ds_PVwithDs_iperr);
    match_Ds_PVwithDs_ipchi2_vec.push_back(match_Ds_PVwithDs_ipchi2);
}

void PVTree::Match_mu_Fill_Vector(){
    match_mu_charge_vec.push_back(match_mu_charge);
    match_mu_eta_vec.push_back(match_mu_eta);
    match_mu_phi_vec.push_back(match_mu_phi);
    match_mu_vx_vec.push_back(match_mu_vx);
    match_mu_vy_vec.push_back(match_mu_vy);
    match_mu_vz_vec.push_back(match_mu_vz);
    match_mu_p_vec.push_back(match_mu_p);
    match_mu_pt_vec.push_back(match_mu_pt);
    match_mu_px_vec.push_back(match_mu_px);
    match_mu_py_vec.push_back(match_mu_py);
    match_mu_pz_vec.push_back(match_mu_pz);
    match_mu_isHighPt_vec.push_back(match_mu_isHighPt);
    match_mu_isLoose_vec.push_back(match_mu_isLoose);
    match_mu_isMedium_vec.push_back(match_mu_isMedium);
    match_mu_isSoft_vec.push_back(match_mu_isSoft);
    match_mu_isTight_vec.push_back(match_mu_isTight);
    match_mu_isPF_vec.push_back(match_mu_isPF);
    match_mu_isTracker_vec.push_back(match_mu_isTracker);
    match_mu_isGlobal_vec.push_back(match_mu_isGlobal);
    match_mu_IsoR03_sumChargedHadronPt_vec.push_back(match_mu_IsoR03_sumChargedHadronPt);
    match_mu_IsoR03_sumChargedParticlePt_vec.push_back(match_mu_IsoR03_sumChargedParticlePt);
    match_mu_IsoR03_sumNeutralHadronEt_vec.push_back(match_mu_IsoR03_sumNeutralHadronEt);
    match_mu_IsoR03_sumPhotonEt_vec.push_back(match_mu_IsoR03_sumPhotonEt);
    match_mu_IsoR03_sumPUPt_vec.push_back(match_mu_IsoR03_sumPUPt);
    match_mu_PFIsoR03_vec.push_back(match_mu_PFIsoR03);
    match_mu_IsoR04_sumChargedHadronPt_vec.push_back(match_mu_IsoR04_sumChargedHadronPt);
    match_mu_IsoR04_sumChargedParticlePt_vec.push_back(match_mu_IsoR04_sumChargedParticlePt);
    match_mu_IsoR04_sumNeutralHadronEt_vec.push_back(match_mu_IsoR04_sumNeutralHadronEt);
    match_mu_IsoR04_sumPhotonEt_vec.push_back(match_mu_IsoR04_sumPhotonEt);
    match_mu_IsoR04_sumPUPt_vec.push_back(match_mu_IsoR04_sumPUPt);
    match_mu_PFIsoR04_vec.push_back(match_mu_PFIsoR04);
    match_mu_primvtx_dxy_vec.push_back(match_mu_primvtx_dxy);
    match_mu_primvtx_dxyerr_vec.push_back(match_mu_primvtx_dxyerr);
    match_mu_primvtx_dz_vec.push_back(match_mu_primvtx_dz);
    match_mu_primvtx_dzerr_vec.push_back(match_mu_primvtx_dzerr);
    match_mu_primvtx_ip_vec.push_back(match_mu_primvtx_ip);
    match_mu_primvtx_iperr_vec.push_back(match_mu_primvtx_iperr);
    match_mu_primvtx_ipchi2_vec.push_back(match_mu_primvtx_ipchi2);
}

void PVTree::Kp_Reset()
{
    Kp_isIsolatedChargedHadron = false;
    Kp_charge = null;
    Kp_eta = null;
    Kp_phi = null;
    Kp_vx = null;
    Kp_vy = null;
    Kp_vz = null;
    Kp_p = null;
    Kp_pt = null;
    Kp_px = null;
    Kp_py = null;
    Kp_pz = null;
}

void PVTree::Km_Reset()
{
    Km_isIsolatedChargedHadron = false;
    Km_charge = null;
    Km_eta = null;
    Km_phi = null;
    Km_vx = null;
    Km_vy = null;
    Km_vz = null;
    Km_p = null;
    Km_pt = null;
    Km_px = null;
    Km_py = null;
    Km_pz = null;

    phi_eta = null;
    phi_phi = null;
    phi_p = null;
    phi_pt = null;
    phi_px = null;
    phi_py = null;
    phi_pz = null;
    phi_invm = null;

    dR_Kp_Km = null;
    dR_Kp_phi = null;
    dR_Km_phi = null;
}

void PVTree::phi_Reset()
{
    phiFit_chi2 = null;
    phiFit_ndof = null;
    phiFit_chi2ndof = null;
    phiFit_vx = null;
    phiFit_vy = null;
    phiFit_vz = null;
    phiFit_vxerr = null;
    phiFit_vyerr = null;
    phiFit_vzerr = null;

    phiFit_Kp_eta = null;
    phiFit_Kp_phi = null;
    phiFit_Kp_p = null;
    phiFit_Kp_pt = null;
    phiFit_Kp_px = null;
    phiFit_Kp_py = null;
    phiFit_Kp_pz = null;

    phiFit_Km_eta = null;
    phiFit_Km_phi = null;
    phiFit_Km_p = null;
    phiFit_Km_pt = null;
    phiFit_Km_px = null;
    phiFit_Km_py = null;
    phiFit_Km_pz = null;

    phiFit_phi_eta = null;
    phiFit_phi_phi = null;
    phiFit_phi_p = null;
    phiFit_phi_pt = null;
    phiFit_phi_px = null;
    phiFit_phi_py = null;
    phiFit_phi_pz = null;
    phiFit_phi_invm = null;

    phiFit_dR_Kp_Km = null;
    phiFit_dR_Kp_phi = null;
    phiFit_dR_Km_phi = null;
}

void PVTree::pi_Reset()
{
    pi_isIsolatedChargedHadron = false;
    pi_charge = null;
    pi_eta = null;
    pi_phi = null;
    pi_vx = null;
    pi_vy = null;
    pi_vz = null;
    pi_p = null;
    pi_pt = null;
    pi_px = null;
    pi_py = null;
    pi_pz = null;

    Ds_eta = null;
    Ds_phi = null;
    Ds_p = null;
    Ds_pt = null;
    Ds_px = null;
    Ds_py = null;
    Ds_pz = null;
    Ds_invm = null;

    dR_Kp_pi = null;
    dR_Km_pi = null;
    dR_pi_phi = null;
    dR_Kp_Ds = null;
    dR_Km_Ds = null;
    dR_pi_Ds = null;
    dR_phi_Ds = null;

    phiFit_pi_eta = null;
    phiFit_pi_phi = null;
    phiFit_pi_p = null;
    phiFit_pi_pt = null;
    phiFit_pi_px = null;
    phiFit_pi_py = null;
    phiFit_pi_pz = null;

    phiFit_Ds_eta = null;
    phiFit_Ds_phi = null;
    phiFit_Ds_p = null;
    phiFit_Ds_pt = null;
    phiFit_Ds_px = null;
    phiFit_Ds_py = null;
    phiFit_Ds_pz = null;
    phiFit_Ds_invm = null;

    phiFit_dR_Kp_pi = null;
    phiFit_dR_Km_pi = null;
    phiFit_dR_pi_phi = null;
    phiFit_dR_Kp_Ds = null;
    phiFit_dR_Km_Ds = null;
    phiFit_dR_pi_Ds = null;
    phiFit_dR_phi_Ds = null;
}

void PVTree::Ds_Reset()
{
    DsFit_chi2 = null;
    DsFit_ndof = null;
    DsFit_chi2ndof = null;
    DsFit_vx = null;
    DsFit_vy = null;
    DsFit_vz = null;
    DsFit_vxerr = null;
    DsFit_vyerr = null;
    DsFit_vzerr = null;

    DsFit_Kp_eta = null;
    DsFit_Kp_phi = null;
    DsFit_Kp_p = null;
    DsFit_Kp_pt = null;
    DsFit_Kp_px = null;
    DsFit_Kp_py = null;
    DsFit_Kp_pz = null;

    DsFit_Km_eta = null;
    DsFit_Km_phi = null;
    DsFit_Km_p = null;
    DsFit_Km_pt = null;
    DsFit_Km_px = null;
    DsFit_Km_py = null;
    DsFit_Km_pz = null;

    DsFit_pi_eta = null;
    DsFit_pi_phi = null;
    DsFit_pi_p = null;
    DsFit_pi_pt = null;
    DsFit_pi_px = null;
    DsFit_pi_py = null;
    DsFit_pi_pz = null;

    DsFit_phi_eta = null;
    DsFit_phi_phi = null;
    DsFit_phi_p = null;
    DsFit_phi_pt = null;
    DsFit_phi_px = null;
    DsFit_phi_py = null;
    DsFit_phi_pz = null;
    DsFit_phi_invm = null;

    DsFit_Ds_eta = null;
    DsFit_Ds_phi = null;
    DsFit_Ds_p = null;
    DsFit_Ds_pt = null;
    DsFit_Ds_px = null;
    DsFit_Ds_py = null;
    DsFit_Ds_pz = null;
    DsFit_Ds_invm = null;

    DsFit_Mconstraint_Ds_invm = null;

    DsFit_dR_Kp_Km = null;
    DsFit_dR_Kp_pi = null;
    DsFit_dR_Km_pi = null;
    DsFit_dR_Kp_phi = null;
    DsFit_dR_Km_phi = null;
    DsFit_dR_pi_phi = null;
    DsFit_dR_Kp_Ds = null;
    DsFit_dR_Km_Ds = null;
    DsFit_dR_pi_Ds = null;
    DsFit_dR_phi_Ds = null;

    dxy_phi_Ds = null;
    dz_phi_Ds = null;

    Ds_primvtx_FDxy = null;
    Ds_primvtx_FDz = null;
    Ds_primvtx_FD = null;
    Ds_primvtx_FDxyerr = null;
    Ds_primvtx_FDxychi2 = null;
    Ds_primvtx_FDzerr = null;
    Ds_primvtx_FDzchi2 = null;
    Ds_primvtx_FDerr = null;
    Ds_primvtx_FDchi2 = null;
    Ds_primvtx_dira = null;
    Ds_primvtx_dira_angle = null;
    Kp_primvtx_ip = null;
    Kp_primvtx_iperr = null;
    Kp_primvtx_ipchi2 = null;
    Km_primvtx_ip = null;
    Km_primvtx_iperr = null;
    Km_primvtx_ipchi2 = null;
    pi_primvtx_ip = null;
    pi_primvtx_iperr = null;
    pi_primvtx_ipchi2 = null;
    phi_primvtx_ip = null;
    phi_primvtx_iperr = null;
    phi_primvtx_ipchi2 = null;
    Ds_primvtx_ip = null;
    Ds_primvtx_iperr = null;
    Ds_primvtx_ipchi2 = null;

    Ds_IsoR03_sumChargedHadronPt = 0.;
    Ds_IsoR03_sumNeutralHadronPt = 0.;
    Ds_IsoR03_sumPhotonPt = 0.;
    Ds_IsoR03_sumPUPt = 0.;
    Ds_PFIsoR03 = 0.;
    Ds_IsoR04_sumChargedHadronPt = 0.;
    Ds_IsoR04_sumNeutralHadronPt = 0.;
    Ds_IsoR04_sumPhotonPt = 0.;
    Ds_IsoR04_sumPUPt = 0.;
    Ds_PFIsoR04 = 0.;

    Kp_match = false;
    Km_match = false;
    pi_match = false;
    match_entry = false;
    non_match_entry = false; 
}

void PVTree::PV_noDs_Reset()
{
    PV_noDs_IsValid = false;
    PV_noDs_IsFake = true;
    PV_noDs_chi2 = null;
    PV_noDs_ndof = null;
    PV_noDs_chi2ndof = null;
    PV_noDs_x = null;
    PV_noDs_y = null;
    PV_noDs_z = null;
    PV_noDs_xerr = null;
    PV_noDs_yerr = null;
    PV_noDs_zerr = null;

    Ds_PVnoDs_FDxy = null;
    Ds_PVnoDs_FDz = null;
    Ds_PVnoDs_FD = null;
    Ds_PVnoDs_FDxyerr = null;
    Ds_PVnoDs_FDxychi2 = null;
    Ds_PVnoDs_FDzerr = null;
    Ds_PVnoDs_FDzchi2 = null;
    Ds_PVnoDs_FDerr = null;
    Ds_PVnoDs_FDchi2 = null;
    Ds_PVnoDs_dira = null;
    Ds_PVnoDs_dira_angle = null;
    Kp_PVnoDs_ip = null;
    Kp_PVnoDs_iperr = null;
    Kp_PVnoDs_ipchi2 = null;
    Km_PVnoDs_ip = null;
    Km_PVnoDs_iperr = null;
    Km_PVnoDs_ipchi2 = null;
    pi_PVnoDs_ip = null;
    pi_PVnoDs_iperr = null;
    pi_PVnoDs_ipchi2 = null;
    phi_PVnoDs_ip = null;
    phi_PVnoDs_iperr = null;
    phi_PVnoDs_ipchi2 = null;
    Ds_PVnoDs_ip = null;
    Ds_PVnoDs_iperr = null;
    Ds_PVnoDs_ipchi2 = null;

    PV_withDs_IsValid = false;
    PV_withDs_IsFake = true;
    PV_withDs_chi2 = null;
    PV_withDs_ndof = null;
    PV_withDs_chi2ndof = null;
    PV_withDs_x = null;
    PV_withDs_y = null;
    PV_withDs_z = null;
    PV_withDs_xerr = null;
    PV_withDs_yerr = null;
    PV_withDs_zerr = null;

    Ds_PVwithDs_FDxy = null;
    Ds_PVwithDs_FDz = null;
    Ds_PVwithDs_FD = null;
    Ds_PVwithDs_FDxyerr = null;
    Ds_PVwithDs_FDxychi2 = null;
    Ds_PVwithDs_FDzerr = null;
    Ds_PVwithDs_FDzchi2 = null;
    Ds_PVwithDs_FDerr = null;
    Ds_PVwithDs_FDchi2 = null;
    Ds_PVwithDs_dira = null;
    Ds_PVwithDs_dira_angle = null;
    Kp_PVwithDs_ip = null;
    Kp_PVwithDs_iperr = null;
    Kp_PVwithDs_ipchi2 = null;
    Km_PVwithDs_ip = null;
    Km_PVwithDs_iperr = null;
    Km_PVwithDs_ipchi2 = null;
    pi_PVwithDs_ip = null;
    pi_PVwithDs_iperr = null;
    pi_PVwithDs_ipchi2 = null;
    phi_PVwithDs_ip = null;
    phi_PVwithDs_iperr = null;
    phi_PVwithDs_ipchi2 = null;
    Ds_PVwithDs_ip = null;
    Ds_PVwithDs_iperr = null;
    Ds_PVwithDs_ipchi2 = null;
}

void PVTree::Ds_Fill_Vector()
{
    Kp_isIsolatedChargedHadron_vec.push_back(Kp_isIsolatedChargedHadron);
    Kp_charge_vec.push_back(Kp_charge);
    Kp_eta_vec.push_back(Kp_eta);
    Kp_phi_vec.push_back(Kp_phi);
    Kp_vx_vec.push_back(Kp_vx);
    Kp_vy_vec.push_back(Kp_vy);
    Kp_vz_vec.push_back(Kp_vz);
    Kp_p_vec.push_back(Kp_p);
    Kp_pt_vec.push_back(Kp_pt);
    Kp_px_vec.push_back(Kp_px);
    Kp_py_vec.push_back(Kp_py);
    Kp_pz_vec.push_back(Kp_pz);

    Km_isIsolatedChargedHadron_vec.push_back(Km_isIsolatedChargedHadron);
    Km_charge_vec.push_back(Km_charge);
    Km_eta_vec.push_back(Km_eta);
    Km_phi_vec.push_back(Km_phi);
    Km_vx_vec.push_back(Km_vx);
    Km_vy_vec.push_back(Km_vy);
    Km_vz_vec.push_back(Km_vz);
    Km_p_vec.push_back(Km_p);
    Km_pt_vec.push_back(Km_pt);
    Km_px_vec.push_back(Km_px);
    Km_py_vec.push_back(Km_py);
    Km_pz_vec.push_back(Km_pz);

    pi_isIsolatedChargedHadron_vec.push_back(pi_isIsolatedChargedHadron);
    pi_charge_vec.push_back(pi_charge);
    pi_eta_vec.push_back(pi_eta);
    pi_phi_vec.push_back(pi_phi);
    pi_vx_vec.push_back(pi_vx);
    pi_vy_vec.push_back(pi_vy);
    pi_vz_vec.push_back(pi_vz);
    pi_p_vec.push_back(pi_p);
    pi_pt_vec.push_back(pi_pt);
    pi_px_vec.push_back(pi_px);
    pi_py_vec.push_back(pi_py);
    pi_pz_vec.push_back(pi_pz);

    phi_eta_vec.push_back(phi_eta);
    phi_phi_vec.push_back(phi_phi);
    phi_p_vec.push_back(phi_p);
    phi_pt_vec.push_back(phi_pt);
    phi_px_vec.push_back(phi_px);
    phi_py_vec.push_back(phi_py);
    phi_pz_vec.push_back(phi_pz);
    phi_invm_vec.push_back(phi_invm);

    Ds_eta_vec.push_back(Ds_eta);
    Ds_phi_vec.push_back(Ds_phi);
    Ds_p_vec.push_back(Ds_p);
    Ds_pt_vec.push_back(Ds_pt);
    Ds_px_vec.push_back(Ds_px);
    Ds_py_vec.push_back(Ds_py);
    Ds_pz_vec.push_back(Ds_pz);
    Ds_invm_vec.push_back(Ds_invm);

    dR_Kp_Km_vec.push_back(dR_Kp_Km);
    dR_Kp_phi_vec.push_back(dR_Kp_phi);
    dR_Km_phi_vec.push_back(dR_Km_phi);
    dR_Kp_pi_vec.push_back(dR_Kp_pi);
    dR_Km_pi_vec.push_back(dR_Km_pi);
    dR_pi_phi_vec.push_back(dR_pi_phi);
    dR_Kp_Ds_vec.push_back(dR_Kp_Ds);
    dR_Km_Ds_vec.push_back(dR_Km_Ds);
    dR_phi_Ds_vec.push_back(dR_phi_Ds);
    dR_pi_Ds_vec.push_back(dR_pi_Ds);

    dxy_phi_Ds_vec.push_back(dxy_phi_Ds);
    dz_phi_Ds_vec.push_back(dz_phi_Ds);

    // Fit on phi
    phiFit_chi2_vec.push_back(phiFit_chi2);
    phiFit_ndof_vec.push_back(phiFit_ndof);
    phiFit_chi2ndof_vec.push_back(phiFit_chi2ndof);
    phiFit_vx_vec.push_back(phiFit_vx);
    phiFit_vy_vec.push_back(phiFit_vy);
    phiFit_vz_vec.push_back(phiFit_vz);
    phiFit_vxerr_vec.push_back(phiFit_vxerr);
    phiFit_vyerr_vec.push_back(phiFit_vyerr);
    phiFit_vzerr_vec.push_back(phiFit_vzerr);

    phiFit_Kp_eta_vec.push_back(phiFit_Kp_eta);
    phiFit_Kp_phi_vec.push_back(phiFit_Kp_phi);
    phiFit_Kp_p_vec.push_back(phiFit_Kp_p);
    phiFit_Kp_pt_vec.push_back(phiFit_Kp_pt);
    phiFit_Kp_px_vec.push_back(phiFit_Kp_px);
    phiFit_Kp_py_vec.push_back(phiFit_Kp_py);
    phiFit_Kp_pz_vec.push_back(phiFit_Kp_pz);

    phiFit_Km_eta_vec.push_back(phiFit_Km_eta);
    phiFit_Km_phi_vec.push_back(phiFit_Km_phi);
    phiFit_Km_p_vec.push_back(phiFit_Km_p);
    phiFit_Km_pt_vec.push_back(phiFit_Km_pt);
    phiFit_Km_px_vec.push_back(phiFit_Km_px);
    phiFit_Km_py_vec.push_back(phiFit_Km_py);
    phiFit_Km_pz_vec.push_back(phiFit_Km_pz);

    phiFit_pi_eta_vec.push_back(phiFit_pi_eta);
    phiFit_pi_phi_vec.push_back(phiFit_pi_phi);
    phiFit_pi_p_vec.push_back(phiFit_pi_p);
    phiFit_pi_pt_vec.push_back(phiFit_pi_pt);
    phiFit_pi_px_vec.push_back(phiFit_pi_px);
    phiFit_pi_py_vec.push_back(phiFit_pi_py);
    phiFit_pi_pz_vec.push_back(phiFit_pi_pz);

    phiFit_phi_eta_vec.push_back(phiFit_phi_eta);
    phiFit_phi_phi_vec.push_back(phiFit_phi_phi);
    phiFit_phi_p_vec.push_back(phiFit_phi_p);
    phiFit_phi_pt_vec.push_back(phiFit_phi_pt);
    phiFit_phi_px_vec.push_back(phiFit_phi_px);
    phiFit_phi_py_vec.push_back(phiFit_phi_py);
    phiFit_phi_pz_vec.push_back(phiFit_phi_pz);
    phiFit_phi_invm_vec.push_back(phiFit_phi_invm);

    phiFit_Ds_eta_vec.push_back(phiFit_Ds_eta);
    phiFit_Ds_phi_vec.push_back(phiFit_Ds_phi);
    phiFit_Ds_p_vec.push_back(phiFit_Ds_p);
    phiFit_Ds_pt_vec.push_back(phiFit_Ds_pt);
    phiFit_Ds_px_vec.push_back(phiFit_Ds_px);
    phiFit_Ds_py_vec.push_back(phiFit_Ds_py);
    phiFit_Ds_pz_vec.push_back(phiFit_Ds_pz);
    phiFit_Ds_invm_vec.push_back(phiFit_Ds_invm);

    phiFit_dR_Kp_Km_vec.push_back(phiFit_dR_Kp_Km);
    phiFit_dR_Kp_phi_vec.push_back(phiFit_dR_Kp_phi);
    phiFit_dR_Km_phi_vec.push_back(phiFit_dR_Km_phi);
    phiFit_dR_Kp_pi_vec.push_back(phiFit_dR_Kp_pi);
    phiFit_dR_Km_pi_vec.push_back(phiFit_dR_Km_pi);
    phiFit_dR_pi_phi_vec.push_back(phiFit_dR_pi_phi);
    phiFit_dR_Kp_Ds_vec.push_back(phiFit_dR_Kp_Ds);
    phiFit_dR_Km_Ds_vec.push_back(phiFit_dR_Km_Ds);
    phiFit_dR_phi_Ds_vec.push_back(phiFit_dR_phi_Ds);
    phiFit_dR_pi_Ds_vec.push_back(phiFit_dR_pi_Ds);

    // Fit on Ds 
    DsFit_chi2_vec.push_back(DsFit_chi2);
    DsFit_ndof_vec.push_back(DsFit_ndof);
    DsFit_chi2ndof_vec.push_back(DsFit_chi2ndof);
    DsFit_vx_vec.push_back(DsFit_vx);
    DsFit_vy_vec.push_back(DsFit_vy);
    DsFit_vz_vec.push_back(DsFit_vz);
    DsFit_vxerr_vec.push_back(DsFit_vxerr);
    DsFit_vyerr_vec.push_back(DsFit_vyerr);
    DsFit_vzerr_vec.push_back(DsFit_vzerr);

    DsFit_Kp_eta_vec.push_back(DsFit_Kp_eta);
    DsFit_Kp_phi_vec.push_back(DsFit_Kp_phi);
    DsFit_Kp_p_vec.push_back(DsFit_Kp_p);
    DsFit_Kp_pt_vec.push_back(DsFit_Kp_pt);
    DsFit_Kp_px_vec.push_back(DsFit_Kp_px);
    DsFit_Kp_py_vec.push_back(DsFit_Kp_py);
    DsFit_Kp_pz_vec.push_back(DsFit_Kp_pz);

    DsFit_Km_eta_vec.push_back(DsFit_Km_eta);
    DsFit_Km_phi_vec.push_back(DsFit_Km_phi);
    DsFit_Km_p_vec.push_back(DsFit_Km_p);
    DsFit_Km_pt_vec.push_back(DsFit_Km_pt);
    DsFit_Km_px_vec.push_back(DsFit_Km_px);
    DsFit_Km_py_vec.push_back(DsFit_Km_py);
    DsFit_Km_pz_vec.push_back(DsFit_Km_pz);

    DsFit_pi_eta_vec.push_back(DsFit_pi_eta);
    DsFit_pi_phi_vec.push_back(DsFit_pi_phi);
    DsFit_pi_p_vec.push_back(DsFit_pi_p);
    DsFit_pi_pt_vec.push_back(DsFit_pi_pt);
    DsFit_pi_px_vec.push_back(DsFit_pi_px);
    DsFit_pi_py_vec.push_back(DsFit_pi_py);
    DsFit_pi_pz_vec.push_back(DsFit_pi_pz);

    DsFit_phi_eta_vec.push_back(DsFit_phi_eta);
    DsFit_phi_phi_vec.push_back(DsFit_phi_phi);
    DsFit_phi_p_vec.push_back(DsFit_phi_p);
    DsFit_phi_pt_vec.push_back(DsFit_phi_pt);
    DsFit_phi_px_vec.push_back(DsFit_phi_px);
    DsFit_phi_py_vec.push_back(DsFit_phi_py);
    DsFit_phi_pz_vec.push_back(DsFit_phi_pz);
    DsFit_phi_invm_vec.push_back(DsFit_phi_invm);

    DsFit_Ds_eta_vec.push_back(DsFit_Ds_eta);
    DsFit_Ds_phi_vec.push_back(DsFit_Ds_phi);
    DsFit_Ds_p_vec.push_back(DsFit_Ds_p);
    DsFit_Ds_pt_vec.push_back(DsFit_Ds_pt);
    DsFit_Ds_px_vec.push_back(DsFit_Ds_px);
    DsFit_Ds_py_vec.push_back(DsFit_Ds_py);
    DsFit_Ds_pz_vec.push_back(DsFit_Ds_pz);
    DsFit_Ds_invm_vec.push_back(DsFit_Ds_invm);

    DsFit_dR_Kp_Km_vec.push_back(DsFit_dR_Kp_Km);
    DsFit_dR_Kp_phi_vec.push_back(DsFit_dR_Kp_phi);
    DsFit_dR_Km_phi_vec.push_back(DsFit_dR_Km_phi);
    DsFit_dR_Kp_pi_vec.push_back(DsFit_dR_Kp_pi);
    DsFit_dR_Km_pi_vec.push_back(DsFit_dR_Km_pi);
    DsFit_dR_pi_phi_vec.push_back(DsFit_dR_pi_phi);
    DsFit_dR_Kp_Ds_vec.push_back(DsFit_dR_Kp_Ds);
    DsFit_dR_Km_Ds_vec.push_back(DsFit_dR_Km_Ds);
    DsFit_dR_phi_Ds_vec.push_back(DsFit_dR_phi_Ds);
    DsFit_dR_pi_Ds_vec.push_back(DsFit_dR_pi_Ds);

    DsFit_Mconstraint_Ds_invm_vec.push_back(DsFit_Mconstraint_Ds_invm);

    Ds_primvtx_FDxy_vec.push_back(Ds_primvtx_FDxy);
    Ds_primvtx_FDz_vec.push_back(Ds_primvtx_FDz);
    Ds_primvtx_FD_vec.push_back(Ds_primvtx_FD);
    Ds_primvtx_FDxyerr_vec.push_back(Ds_primvtx_FDxyerr);
    Ds_primvtx_FDxychi2_vec.push_back(Ds_primvtx_FDxychi2);
    Ds_primvtx_FDzerr_vec.push_back(Ds_primvtx_FDzerr);
    Ds_primvtx_FDzchi2_vec.push_back(Ds_primvtx_FDzchi2);
    Ds_primvtx_FDerr_vec.push_back(Ds_primvtx_FDerr);
    Ds_primvtx_FDchi2_vec.push_back(Ds_primvtx_FDchi2);
    Ds_primvtx_dira_vec.push_back(Ds_primvtx_dira);
    Ds_primvtx_dira_angle_vec.push_back(Ds_primvtx_dira_angle);
    Kp_primvtx_ip_vec.push_back(Kp_primvtx_ip);
    Kp_primvtx_iperr_vec.push_back(Kp_primvtx_iperr);
    Kp_primvtx_ipchi2_vec.push_back(Kp_primvtx_ipchi2);
    Km_primvtx_ip_vec.push_back(Km_primvtx_ip);
    Km_primvtx_iperr_vec.push_back(Km_primvtx_iperr);
    Km_primvtx_ipchi2_vec.push_back(Km_primvtx_ipchi2);
    pi_primvtx_ip_vec.push_back(pi_primvtx_ip);
    pi_primvtx_iperr_vec.push_back(pi_primvtx_iperr);
    pi_primvtx_ipchi2_vec.push_back(pi_primvtx_ipchi2);
    phi_primvtx_ip_vec.push_back(phi_primvtx_ip);
    phi_primvtx_iperr_vec.push_back(phi_primvtx_iperr);
    phi_primvtx_ipchi2_vec.push_back(phi_primvtx_ipchi2);
    Ds_primvtx_ip_vec.push_back(Ds_primvtx_ip);
    Ds_primvtx_iperr_vec.push_back(Ds_primvtx_iperr);
    Ds_primvtx_ipchi2_vec.push_back(Ds_primvtx_ipchi2);

    Ds_IsoR03_sumChargedHadronPt_vec.push_back(Ds_IsoR03_sumChargedHadronPt);
    Ds_IsoR03_sumNeutralHadronPt_vec.push_back(Ds_IsoR03_sumNeutralHadronPt);
    Ds_IsoR03_sumPhotonPt_vec.push_back(Ds_IsoR03_sumPhotonPt);
    Ds_IsoR03_sumPUPt_vec.push_back(Ds_IsoR03_sumPUPt);
    Ds_PFIsoR03_vec.push_back(Ds_PFIsoR03);
    Ds_IsoR04_sumChargedHadronPt_vec.push_back(Ds_IsoR04_sumChargedHadronPt);
    Ds_IsoR04_sumNeutralHadronPt_vec.push_back(Ds_IsoR04_sumNeutralHadronPt);
    Ds_IsoR04_sumPhotonPt_vec.push_back(Ds_IsoR04_sumPhotonPt);
    Ds_IsoR04_sumPUPt_vec.push_back(Ds_IsoR04_sumPUPt);
    Ds_PFIsoR04_vec.push_back(Ds_PFIsoR04);

    PV_noDs_IsValid_vec.push_back(PV_noDs_IsValid);
    PV_noDs_IsFake_vec.push_back(PV_noDs_IsFake);
    PV_noDs_chi2_vec.push_back(PV_noDs_chi2);
    PV_noDs_ndof_vec.push_back(PV_noDs_ndof);
    PV_noDs_chi2ndof_vec.push_back(PV_noDs_chi2ndof);
    PV_noDs_x_vec.push_back(PV_noDs_x);
    PV_noDs_y_vec.push_back(PV_noDs_y);
    PV_noDs_z_vec.push_back(PV_noDs_z);
    PV_noDs_xerr_vec.push_back(PV_noDs_xerr);
    PV_noDs_yerr_vec.push_back(PV_noDs_yerr);
    PV_noDs_zerr_vec.push_back(PV_noDs_zerr);

    Ds_PVnoDs_FDxy_vec.push_back(Ds_PVnoDs_FDxy);
    Ds_PVnoDs_FDz_vec.push_back(Ds_PVnoDs_FDz);
    Ds_PVnoDs_FD_vec.push_back(Ds_PVnoDs_FD);
    Ds_PVnoDs_FDxyerr_vec.push_back(Ds_PVnoDs_FDxyerr);
    Ds_PVnoDs_FDxychi2_vec.push_back(Ds_PVnoDs_FDxychi2);
    Ds_PVnoDs_FDzerr_vec.push_back(Ds_PVnoDs_FDzerr);
    Ds_PVnoDs_FDzchi2_vec.push_back(Ds_PVnoDs_FDzchi2);
    Ds_PVnoDs_FDerr_vec.push_back(Ds_PVnoDs_FDerr);
    Ds_PVnoDs_FDchi2_vec.push_back(Ds_PVnoDs_FDchi2);
    Ds_PVnoDs_dira_vec.push_back(Ds_PVnoDs_dira);
    Ds_PVnoDs_dira_angle_vec.push_back(Ds_PVnoDs_dira_angle);
    Kp_PVnoDs_ip_vec.push_back(Kp_PVnoDs_ip);
    Kp_PVnoDs_iperr_vec.push_back(Kp_PVnoDs_iperr);
    Kp_PVnoDs_ipchi2_vec.push_back(Kp_PVnoDs_ipchi2);
    Km_PVnoDs_ip_vec.push_back(Km_PVnoDs_ip);
    Km_PVnoDs_iperr_vec.push_back(Km_PVnoDs_iperr);
    Km_PVnoDs_ipchi2_vec.push_back(Km_PVnoDs_ipchi2);
    pi_PVnoDs_ip_vec.push_back(pi_PVnoDs_ip);
    pi_PVnoDs_iperr_vec.push_back(pi_PVnoDs_iperr);
    pi_PVnoDs_ipchi2_vec.push_back(pi_PVnoDs_ipchi2);
    phi_PVnoDs_ip_vec.push_back(phi_PVnoDs_ip);
    phi_PVnoDs_iperr_vec.push_back(phi_PVnoDs_iperr);
    phi_PVnoDs_ipchi2_vec.push_back(phi_PVnoDs_ipchi2);
    Ds_PVnoDs_ip_vec.push_back(Ds_PVnoDs_ip);
    Ds_PVnoDs_iperr_vec.push_back(Ds_PVnoDs_iperr);
    Ds_PVnoDs_ipchi2_vec.push_back(Ds_PVnoDs_ipchi2);

    PV_withDs_IsValid_vec.push_back(PV_withDs_IsValid);
    PV_withDs_IsFake_vec.push_back(PV_withDs_IsFake);
    PV_withDs_chi2_vec.push_back(PV_withDs_chi2);
    PV_withDs_ndof_vec.push_back(PV_withDs_ndof);
    PV_withDs_chi2ndof_vec.push_back(PV_withDs_chi2ndof);
    PV_withDs_x_vec.push_back(PV_withDs_x);
    PV_withDs_y_vec.push_back(PV_withDs_y);
    PV_withDs_z_vec.push_back(PV_withDs_z);
    PV_withDs_xerr_vec.push_back(PV_withDs_xerr);
    PV_withDs_yerr_vec.push_back(PV_withDs_yerr);
    PV_withDs_zerr_vec.push_back(PV_withDs_zerr);

    Ds_PVwithDs_FDxy_vec.push_back(Ds_PVwithDs_FDxy);
    Ds_PVwithDs_FDz_vec.push_back(Ds_PVwithDs_FDz);
    Ds_PVwithDs_FD_vec.push_back(Ds_PVwithDs_FD);
    Ds_PVwithDs_FDxyerr_vec.push_back(Ds_PVwithDs_FDxyerr);
    Ds_PVwithDs_FDxychi2_vec.push_back(Ds_PVwithDs_FDxychi2);
    Ds_PVwithDs_FDzerr_vec.push_back(Ds_PVwithDs_FDzerr);
    Ds_PVwithDs_FDzchi2_vec.push_back(Ds_PVwithDs_FDzchi2);
    Ds_PVwithDs_FDerr_vec.push_back(Ds_PVwithDs_FDerr);
    Ds_PVwithDs_FDchi2_vec.push_back(Ds_PVwithDs_FDchi2);
    Ds_PVwithDs_dira_vec.push_back(Ds_PVwithDs_dira);
    Ds_PVwithDs_dira_angle_vec.push_back(Ds_PVwithDs_dira_angle);
    Kp_PVwithDs_ip_vec.push_back(Kp_PVwithDs_ip);
    Kp_PVwithDs_iperr_vec.push_back(Kp_PVwithDs_iperr);
    Kp_PVwithDs_ipchi2_vec.push_back(Kp_PVwithDs_ipchi2);
    Km_PVwithDs_ip_vec.push_back(Km_PVwithDs_ip);
    Km_PVwithDs_iperr_vec.push_back(Km_PVwithDs_iperr);
    Km_PVwithDs_ipchi2_vec.push_back(Km_PVwithDs_ipchi2);
    pi_PVwithDs_ip_vec.push_back(pi_PVwithDs_ip);
    pi_PVwithDs_iperr_vec.push_back(pi_PVwithDs_iperr);
    pi_PVwithDs_ipchi2_vec.push_back(pi_PVwithDs_ipchi2);
    phi_PVwithDs_ip_vec.push_back(phi_PVwithDs_ip);
    phi_PVwithDs_iperr_vec.push_back(phi_PVwithDs_iperr);
    phi_PVwithDs_ipchi2_vec.push_back(phi_PVwithDs_ipchi2);
    Ds_PVwithDs_ip_vec.push_back(Ds_PVwithDs_ip);
    Ds_PVwithDs_iperr_vec.push_back(Ds_PVwithDs_iperr);
    Ds_PVwithDs_ipchi2_vec.push_back(Ds_PVwithDs_ipchi2);

    Kp_match_vec.push_back(Kp_match);
    Km_match_vec.push_back(Km_match);
    pi_match_vec.push_back(pi_match);
    match_entry_vec.push_back(match_entry);
    non_match_entry_vec.push_back(non_match_entry);
}

void PVTree::Best_Ds_Fill_Vector(int idxmax)
{
    best_Kp_isIsolatedChargedHadron_vec.push_back(Kp_isIsolatedChargedHadron_vec[idxmax]);
    best_Kp_charge_vec.push_back(Kp_charge_vec[idxmax]);
    best_Kp_eta_vec.push_back(Kp_eta_vec[idxmax]);
    best_Kp_phi_vec.push_back(Kp_phi_vec[idxmax]);
    best_Kp_vx_vec.push_back(Kp_vx_vec[idxmax]);
    best_Kp_vy_vec.push_back(Kp_vy_vec[idxmax]);
    best_Kp_vz_vec.push_back(Kp_vz_vec[idxmax]);
    best_Kp_p_vec.push_back(Kp_p_vec[idxmax]);
    best_Kp_pt_vec.push_back(Kp_pt_vec[idxmax]);
    best_Kp_px_vec.push_back(Kp_px_vec[idxmax]);
    best_Kp_py_vec.push_back(Kp_py_vec[idxmax]);
    best_Kp_pz_vec.push_back(Kp_pz_vec[idxmax]);

    best_Km_isIsolatedChargedHadron_vec.push_back(Km_isIsolatedChargedHadron_vec[idxmax]);
    best_Km_charge_vec.push_back(Km_charge_vec[idxmax]);
    best_Km_eta_vec.push_back(Km_eta_vec[idxmax]);
    best_Km_phi_vec.push_back(Km_phi_vec[idxmax]);
    best_Km_vx_vec.push_back(Km_vx_vec[idxmax]);
    best_Km_vy_vec.push_back(Km_vy_vec[idxmax]);
    best_Km_vz_vec.push_back(Km_vz_vec[idxmax]);
    best_Km_p_vec.push_back(Km_p_vec[idxmax]);
    best_Km_pt_vec.push_back(Km_pt_vec[idxmax]);
    best_Km_px_vec.push_back(Km_px_vec[idxmax]);
    best_Km_py_vec.push_back(Km_py_vec[idxmax]);
    best_Km_pz_vec.push_back(Km_pz_vec[idxmax]);

    best_pi_isIsolatedChargedHadron_vec.push_back(pi_isIsolatedChargedHadron_vec[idxmax]);
    best_pi_charge_vec.push_back(pi_charge_vec[idxmax]);
    best_pi_eta_vec.push_back(pi_eta_vec[idxmax]);
    best_pi_phi_vec.push_back(pi_phi_vec[idxmax]);
    best_pi_vx_vec.push_back(pi_vx_vec[idxmax]);
    best_pi_vy_vec.push_back(pi_vy_vec[idxmax]);
    best_pi_vz_vec.push_back(pi_vz_vec[idxmax]);
    best_pi_p_vec.push_back(pi_p_vec[idxmax]);
    best_pi_pt_vec.push_back(pi_pt_vec[idxmax]);
    best_pi_px_vec.push_back(pi_px_vec[idxmax]);
    best_pi_py_vec.push_back(pi_py_vec[idxmax]);
    best_pi_pz_vec.push_back(pi_pz_vec[idxmax]);

    best_phi_eta_vec.push_back(phi_eta_vec[idxmax]);
    best_phi_phi_vec.push_back(phi_phi_vec[idxmax]);
    best_phi_p_vec.push_back(phi_p_vec[idxmax]);
    best_phi_pt_vec.push_back(phi_pt_vec[idxmax]);
    best_phi_px_vec.push_back(phi_px_vec[idxmax]);
    best_phi_py_vec.push_back(phi_py_vec[idxmax]);
    best_phi_pz_vec.push_back(phi_pz_vec[idxmax]);
    best_phi_invm_vec.push_back(phi_invm_vec[idxmax]);

    best_Ds_eta_vec.push_back(Ds_eta_vec[idxmax]);
    best_Ds_phi_vec.push_back(Ds_phi_vec[idxmax]);
    best_Ds_p_vec.push_back(Ds_p_vec[idxmax]);
    best_Ds_pt_vec.push_back(Ds_pt_vec[idxmax]);
    best_Ds_px_vec.push_back(Ds_px_vec[idxmax]);
    best_Ds_py_vec.push_back(Ds_py_vec[idxmax]);
    best_Ds_pz_vec.push_back(Ds_pz_vec[idxmax]);
    best_Ds_invm_vec.push_back(Ds_invm_vec[idxmax]);

    best_dR_Kp_Km_vec.push_back(dR_Kp_Km_vec[idxmax]);
    best_dR_Kp_phi_vec.push_back(dR_Kp_phi_vec[idxmax]);
    best_dR_Km_phi_vec.push_back(dR_Km_phi_vec[idxmax]);
    best_dR_Kp_pi_vec.push_back(dR_Kp_pi_vec[idxmax]);
    best_dR_Km_pi_vec.push_back(dR_Km_pi_vec[idxmax]);
    best_dR_pi_phi_vec.push_back(dR_pi_phi_vec[idxmax]);
    best_dR_Kp_Ds_vec.push_back(dR_Kp_Ds_vec[idxmax]);
    best_dR_Km_Ds_vec.push_back(dR_Km_Ds_vec[idxmax]);
    best_dR_phi_Ds_vec.push_back(dR_phi_Ds_vec[idxmax]);
    best_dR_pi_Ds_vec.push_back(dR_pi_Ds_vec[idxmax]);

    best_dxy_phi_Ds_vec.push_back(dxy_phi_Ds_vec[idxmax]);
    best_dz_phi_Ds_vec.push_back(dz_phi_Ds_vec[idxmax]);

    // Fit on phi
    best_phiFit_chi2_vec.push_back(phiFit_chi2_vec[idxmax]);
    best_phiFit_ndof_vec.push_back(phiFit_ndof_vec[idxmax]);
    best_phiFit_chi2ndof_vec.push_back(phiFit_chi2ndof_vec[idxmax]);
    best_phiFit_vx_vec.push_back(phiFit_vx_vec[idxmax]);
    best_phiFit_vy_vec.push_back(phiFit_vy_vec[idxmax]);
    best_phiFit_vz_vec.push_back(phiFit_vz_vec[idxmax]);
    best_phiFit_vxerr_vec.push_back(phiFit_vxerr_vec[idxmax]);
    best_phiFit_vyerr_vec.push_back(phiFit_vyerr_vec[idxmax]);
    best_phiFit_vzerr_vec.push_back(phiFit_vzerr_vec[idxmax]);

    best_phiFit_Kp_eta_vec.push_back(phiFit_Kp_eta_vec[idxmax]);
    best_phiFit_Kp_phi_vec.push_back(phiFit_Kp_phi_vec[idxmax]);
    best_phiFit_Kp_p_vec.push_back(phiFit_Kp_p_vec[idxmax]);
    best_phiFit_Kp_pt_vec.push_back(phiFit_Kp_pt_vec[idxmax]);
    best_phiFit_Kp_px_vec.push_back(phiFit_Kp_px_vec[idxmax]);
    best_phiFit_Kp_py_vec.push_back(phiFit_Kp_py_vec[idxmax]);
    best_phiFit_Kp_pz_vec.push_back(phiFit_Kp_pz_vec[idxmax]);

    best_phiFit_Km_eta_vec.push_back(phiFit_Km_eta_vec[idxmax]);
    best_phiFit_Km_phi_vec.push_back(phiFit_Km_phi_vec[idxmax]);
    best_phiFit_Km_p_vec.push_back(phiFit_Km_p_vec[idxmax]);
    best_phiFit_Km_pt_vec.push_back(phiFit_Km_pt_vec[idxmax]);
    best_phiFit_Km_px_vec.push_back(phiFit_Km_px_vec[idxmax]);
    best_phiFit_Km_py_vec.push_back(phiFit_Km_py_vec[idxmax]);
    best_phiFit_Km_pz_vec.push_back(phiFit_Km_pz_vec[idxmax]);

    best_phiFit_pi_eta_vec.push_back(phiFit_pi_eta_vec[idxmax]);
    best_phiFit_pi_phi_vec.push_back(phiFit_pi_phi_vec[idxmax]);
    best_phiFit_pi_p_vec.push_back(phiFit_pi_p_vec[idxmax]);
    best_phiFit_pi_pt_vec.push_back(phiFit_pi_pt_vec[idxmax]);
    best_phiFit_pi_px_vec.push_back(phiFit_pi_px_vec[idxmax]);
    best_phiFit_pi_py_vec.push_back(phiFit_pi_py_vec[idxmax]);
    best_phiFit_pi_pz_vec.push_back(phiFit_pi_pz_vec[idxmax]);

    best_phiFit_phi_eta_vec.push_back(phiFit_phi_eta_vec[idxmax]);
    best_phiFit_phi_phi_vec.push_back(phiFit_phi_phi_vec[idxmax]);
    best_phiFit_phi_p_vec.push_back(phiFit_phi_p_vec[idxmax]);
    best_phiFit_phi_pt_vec.push_back(phiFit_phi_pt_vec[idxmax]);
    best_phiFit_phi_px_vec.push_back(phiFit_phi_px_vec[idxmax]);
    best_phiFit_phi_py_vec.push_back(phiFit_phi_py_vec[idxmax]);
    best_phiFit_phi_pz_vec.push_back(phiFit_phi_pz_vec[idxmax]);
    best_phiFit_phi_invm_vec.push_back(phiFit_phi_invm_vec[idxmax]);

    best_phiFit_Ds_eta_vec.push_back(phiFit_Ds_eta_vec[idxmax]);
    best_phiFit_Ds_phi_vec.push_back(phiFit_Ds_phi_vec[idxmax]);
    best_phiFit_Ds_p_vec.push_back(phiFit_Ds_p_vec[idxmax]);
    best_phiFit_Ds_pt_vec.push_back(phiFit_Ds_pt_vec[idxmax]);
    best_phiFit_Ds_px_vec.push_back(phiFit_Ds_px_vec[idxmax]);
    best_phiFit_Ds_py_vec.push_back(phiFit_Ds_py_vec[idxmax]);
    best_phiFit_Ds_pz_vec.push_back(phiFit_Ds_pz_vec[idxmax]);
    best_phiFit_Ds_invm_vec.push_back(phiFit_Ds_invm_vec[idxmax]);

    best_phiFit_dR_Kp_Km_vec.push_back(phiFit_dR_Kp_Km_vec[idxmax]);
    best_phiFit_dR_Kp_phi_vec.push_back(phiFit_dR_Kp_phi_vec[idxmax]);
    best_phiFit_dR_Km_phi_vec.push_back(phiFit_dR_Km_phi_vec[idxmax]);
    best_phiFit_dR_Kp_pi_vec.push_back(phiFit_dR_Kp_pi_vec[idxmax]);
    best_phiFit_dR_Km_pi_vec.push_back(phiFit_dR_Km_pi_vec[idxmax]);
    best_phiFit_dR_pi_phi_vec.push_back(phiFit_dR_pi_phi_vec[idxmax]);
    best_phiFit_dR_Kp_Ds_vec.push_back(phiFit_dR_Kp_Ds_vec[idxmax]);
    best_phiFit_dR_Km_Ds_vec.push_back(phiFit_dR_Km_Ds_vec[idxmax]);
    best_phiFit_dR_phi_Ds_vec.push_back(phiFit_dR_phi_Ds_vec[idxmax]);
    best_phiFit_dR_pi_Ds_vec.push_back(phiFit_dR_pi_Ds_vec[idxmax]);

    // Fit on Ds 
    best_DsFit_chi2_vec.push_back(DsFit_chi2_vec[idxmax]);
    best_DsFit_ndof_vec.push_back(DsFit_ndof_vec[idxmax]);
    best_DsFit_chi2ndof_vec.push_back(DsFit_chi2ndof_vec[idxmax]);
    best_DsFit_vx_vec.push_back(DsFit_vx_vec[idxmax]);
    best_DsFit_vy_vec.push_back(DsFit_vy_vec[idxmax]);
    best_DsFit_vz_vec.push_back(DsFit_vz_vec[idxmax]);
    best_DsFit_vxerr_vec.push_back(DsFit_vxerr_vec[idxmax]);
    best_DsFit_vyerr_vec.push_back(DsFit_vyerr_vec[idxmax]);
    best_DsFit_vzerr_vec.push_back(DsFit_vzerr_vec[idxmax]);

    best_DsFit_Kp_eta_vec.push_back(DsFit_Kp_eta_vec[idxmax]);
    best_DsFit_Kp_phi_vec.push_back(DsFit_Kp_phi_vec[idxmax]);
    best_DsFit_Kp_p_vec.push_back(DsFit_Kp_p_vec[idxmax]);
    best_DsFit_Kp_pt_vec.push_back(DsFit_Kp_pt_vec[idxmax]);
    best_DsFit_Kp_px_vec.push_back(DsFit_Kp_px_vec[idxmax]);
    best_DsFit_Kp_py_vec.push_back(DsFit_Kp_py_vec[idxmax]);
    best_DsFit_Kp_pz_vec.push_back(DsFit_Kp_pz_vec[idxmax]);

    best_DsFit_Km_eta_vec.push_back(DsFit_Km_eta_vec[idxmax]);
    best_DsFit_Km_phi_vec.push_back(DsFit_Km_phi_vec[idxmax]);
    best_DsFit_Km_p_vec.push_back(DsFit_Km_p_vec[idxmax]);
    best_DsFit_Km_pt_vec.push_back(DsFit_Km_pt_vec[idxmax]);
    best_DsFit_Km_px_vec.push_back(DsFit_Km_px_vec[idxmax]);
    best_DsFit_Km_py_vec.push_back(DsFit_Km_py_vec[idxmax]);
    best_DsFit_Km_pz_vec.push_back(DsFit_Km_pz_vec[idxmax]);

    best_DsFit_pi_eta_vec.push_back(DsFit_pi_eta_vec[idxmax]);
    best_DsFit_pi_phi_vec.push_back(DsFit_pi_phi_vec[idxmax]);
    best_DsFit_pi_p_vec.push_back(DsFit_pi_p_vec[idxmax]);
    best_DsFit_pi_pt_vec.push_back(DsFit_pi_pt_vec[idxmax]);
    best_DsFit_pi_px_vec.push_back(DsFit_pi_px_vec[idxmax]);
    best_DsFit_pi_py_vec.push_back(DsFit_pi_py_vec[idxmax]);
    best_DsFit_pi_pz_vec.push_back(DsFit_pi_pz_vec[idxmax]);

    best_DsFit_phi_eta_vec.push_back(DsFit_phi_eta_vec[idxmax]);
    best_DsFit_phi_phi_vec.push_back(DsFit_phi_phi_vec[idxmax]);
    best_DsFit_phi_p_vec.push_back(DsFit_phi_p_vec[idxmax]);
    best_DsFit_phi_pt_vec.push_back(DsFit_phi_pt_vec[idxmax]);
    best_DsFit_phi_px_vec.push_back(DsFit_phi_px_vec[idxmax]);
    best_DsFit_phi_py_vec.push_back(DsFit_phi_py_vec[idxmax]);
    best_DsFit_phi_pz_vec.push_back(DsFit_phi_pz_vec[idxmax]);
    best_DsFit_phi_invm_vec.push_back(DsFit_phi_invm_vec[idxmax]);

    best_DsFit_Ds_eta_vec.push_back(DsFit_Ds_eta_vec[idxmax]);
    best_DsFit_Ds_phi_vec.push_back(DsFit_Ds_phi_vec[idxmax]);
    best_DsFit_Ds_p_vec.push_back(DsFit_Ds_p_vec[idxmax]);
    best_DsFit_Ds_pt_vec.push_back(DsFit_Ds_pt_vec[idxmax]);
    best_DsFit_Ds_px_vec.push_back(DsFit_Ds_px_vec[idxmax]);
    best_DsFit_Ds_py_vec.push_back(DsFit_Ds_py_vec[idxmax]);
    best_DsFit_Ds_pz_vec.push_back(DsFit_Ds_pz_vec[idxmax]);
    best_DsFit_Ds_invm_vec.push_back(DsFit_Ds_invm_vec[idxmax]);

    best_DsFit_dR_Kp_Km_vec.push_back(DsFit_dR_Kp_Km_vec[idxmax]);
    best_DsFit_dR_Kp_phi_vec.push_back(DsFit_dR_Kp_phi_vec[idxmax]);
    best_DsFit_dR_Km_phi_vec.push_back(DsFit_dR_Km_phi_vec[idxmax]);
    best_DsFit_dR_Kp_pi_vec.push_back(DsFit_dR_Kp_pi_vec[idxmax]);
    best_DsFit_dR_Km_pi_vec.push_back(DsFit_dR_Km_pi_vec[idxmax]);
    best_DsFit_dR_pi_phi_vec.push_back(DsFit_dR_pi_phi_vec[idxmax]);
    best_DsFit_dR_Kp_Ds_vec.push_back(DsFit_dR_Kp_Ds_vec[idxmax]);
    best_DsFit_dR_Km_Ds_vec.push_back(DsFit_dR_Km_Ds_vec[idxmax]);
    best_DsFit_dR_phi_Ds_vec.push_back(DsFit_dR_phi_Ds_vec[idxmax]);
    best_DsFit_dR_pi_Ds_vec.push_back(DsFit_dR_pi_Ds_vec[idxmax]);

    best_DsFit_Mconstraint_Ds_invm_vec.push_back(DsFit_Mconstraint_Ds_invm_vec[idxmax]);

    best_Ds_primvtx_FDxy_vec.push_back(Ds_primvtx_FDxy_vec[idxmax]);
    best_Ds_primvtx_FDz_vec.push_back(Ds_primvtx_FDz_vec[idxmax]);
    best_Ds_primvtx_FD_vec.push_back(Ds_primvtx_FD_vec[idxmax]);
    best_Ds_primvtx_FDxyerr_vec.push_back(Ds_primvtx_FDxyerr_vec[idxmax]);
    best_Ds_primvtx_FDxychi2_vec.push_back(Ds_primvtx_FDxychi2_vec[idxmax]);
    best_Ds_primvtx_FDzerr_vec.push_back(Ds_primvtx_FDzerr_vec[idxmax]);
    best_Ds_primvtx_FDzchi2_vec.push_back(Ds_primvtx_FDzchi2_vec[idxmax]);
    best_Ds_primvtx_FDerr_vec.push_back(Ds_primvtx_FDerr_vec[idxmax]);
    best_Ds_primvtx_FDchi2_vec.push_back(Ds_primvtx_FDchi2_vec[idxmax]);
    best_Ds_primvtx_dira_vec.push_back(Ds_primvtx_dira_vec[idxmax]);
    best_Ds_primvtx_dira_angle_vec.push_back(Ds_primvtx_dira_angle_vec[idxmax]);
    best_Kp_primvtx_ip_vec.push_back(Kp_primvtx_ip_vec[idxmax]);
    best_Kp_primvtx_iperr_vec.push_back(Kp_primvtx_iperr_vec[idxmax]);
    best_Kp_primvtx_ipchi2_vec.push_back(Kp_primvtx_ipchi2_vec[idxmax]);
    best_Km_primvtx_ip_vec.push_back(Km_primvtx_ip_vec[idxmax]);
    best_Km_primvtx_iperr_vec.push_back(Km_primvtx_iperr_vec[idxmax]);
    best_Km_primvtx_ipchi2_vec.push_back(Km_primvtx_ipchi2_vec[idxmax]);
    best_pi_primvtx_ip_vec.push_back(pi_primvtx_ip_vec[idxmax]);
    best_pi_primvtx_iperr_vec.push_back(pi_primvtx_iperr_vec[idxmax]);
    best_pi_primvtx_ipchi2_vec.push_back(pi_primvtx_ipchi2_vec[idxmax]);
    best_phi_primvtx_ip_vec.push_back(phi_primvtx_ip_vec[idxmax]);
    best_phi_primvtx_iperr_vec.push_back(phi_primvtx_iperr_vec[idxmax]);
    best_phi_primvtx_ipchi2_vec.push_back(phi_primvtx_ipchi2_vec[idxmax]);
    best_Ds_primvtx_ip_vec.push_back(Ds_primvtx_ip_vec[idxmax]);
    best_Ds_primvtx_iperr_vec.push_back(Ds_primvtx_iperr_vec[idxmax]);
    best_Ds_primvtx_ipchi2_vec.push_back(Ds_primvtx_ipchi2_vec[idxmax]);

    best_Ds_IsoR03_sumChargedHadronPt_vec.push_back(Ds_IsoR03_sumChargedHadronPt_vec[idxmax]);
    best_Ds_IsoR03_sumNeutralHadronPt_vec.push_back(Ds_IsoR03_sumNeutralHadronPt_vec[idxmax]);
    best_Ds_IsoR03_sumPhotonPt_vec.push_back(Ds_IsoR03_sumPhotonPt_vec[idxmax]);
    best_Ds_IsoR03_sumPUPt_vec.push_back(Ds_IsoR03_sumPUPt_vec[idxmax]);
    best_Ds_PFIsoR03_vec.push_back(Ds_PFIsoR03_vec[idxmax]);
    best_Ds_IsoR04_sumChargedHadronPt_vec.push_back(Ds_IsoR04_sumChargedHadronPt_vec[idxmax]);
    best_Ds_IsoR04_sumNeutralHadronPt_vec.push_back(Ds_IsoR04_sumNeutralHadronPt_vec[idxmax]);
    best_Ds_IsoR04_sumPhotonPt_vec.push_back(Ds_IsoR04_sumPhotonPt_vec[idxmax]);
    best_Ds_IsoR04_sumPUPt_vec.push_back(Ds_IsoR04_sumPUPt_vec[idxmax]);
    best_Ds_PFIsoR04_vec.push_back(Ds_PFIsoR04_vec[idxmax]);

    best_PV_noDs_IsValid_vec.push_back(PV_noDs_IsValid_vec[idxmax]);
    best_PV_noDs_IsFake_vec.push_back(PV_noDs_IsFake_vec[idxmax]);
    best_PV_noDs_chi2_vec.push_back(PV_noDs_chi2_vec[idxmax]);
    best_PV_noDs_ndof_vec.push_back(PV_noDs_ndof_vec[idxmax]);
    best_PV_noDs_chi2ndof_vec.push_back(PV_noDs_chi2ndof_vec[idxmax]);
    best_PV_noDs_x_vec.push_back(PV_noDs_x_vec[idxmax]);
    best_PV_noDs_y_vec.push_back(PV_noDs_y_vec[idxmax]);
    best_PV_noDs_z_vec.push_back(PV_noDs_z_vec[idxmax]);
    best_PV_noDs_xerr_vec.push_back(PV_noDs_xerr_vec[idxmax]);
    best_PV_noDs_yerr_vec.push_back(PV_noDs_yerr_vec[idxmax]);
    best_PV_noDs_zerr_vec.push_back(PV_noDs_zerr_vec[idxmax]);

    best_Ds_PVnoDs_FDxy_vec.push_back(Ds_PVnoDs_FDxy_vec[idxmax]);
    best_Ds_PVnoDs_FDz_vec.push_back(Ds_PVnoDs_FDz_vec[idxmax]);
    best_Ds_PVnoDs_FD_vec.push_back(Ds_PVnoDs_FD_vec[idxmax]);
    best_Ds_PVnoDs_FDxyerr_vec.push_back(Ds_PVnoDs_FDxyerr_vec[idxmax]);
    best_Ds_PVnoDs_FDxychi2_vec.push_back(Ds_PVnoDs_FDxychi2_vec[idxmax]);
    best_Ds_PVnoDs_FDzerr_vec.push_back(Ds_PVnoDs_FDzerr_vec[idxmax]);
    best_Ds_PVnoDs_FDzchi2_vec.push_back(Ds_PVnoDs_FDzchi2_vec[idxmax]);
    best_Ds_PVnoDs_FDerr_vec.push_back(Ds_PVnoDs_FDerr_vec[idxmax]);
    best_Ds_PVnoDs_FDchi2_vec.push_back(Ds_PVnoDs_FDchi2_vec[idxmax]);
    best_Ds_PVnoDs_dira_vec.push_back(Ds_PVnoDs_dira_vec[idxmax]);
    best_Ds_PVnoDs_dira_angle_vec.push_back(Ds_PVnoDs_dira_angle_vec[idxmax]);
    best_Kp_PVnoDs_ip_vec.push_back(Kp_PVnoDs_ip_vec[idxmax]);
    best_Kp_PVnoDs_iperr_vec.push_back(Kp_PVnoDs_iperr_vec[idxmax]);
    best_Kp_PVnoDs_ipchi2_vec.push_back(Kp_PVnoDs_ipchi2_vec[idxmax]);
    best_Km_PVnoDs_ip_vec.push_back(Km_PVnoDs_ip_vec[idxmax]);
    best_Km_PVnoDs_iperr_vec.push_back(Km_PVnoDs_iperr_vec[idxmax]);
    best_Km_PVnoDs_ipchi2_vec.push_back(Km_PVnoDs_ipchi2_vec[idxmax]);
    best_pi_PVnoDs_ip_vec.push_back(pi_PVnoDs_ip_vec[idxmax]);
    best_pi_PVnoDs_iperr_vec.push_back(pi_PVnoDs_iperr_vec[idxmax]);
    best_pi_PVnoDs_ipchi2_vec.push_back(pi_PVnoDs_ipchi2_vec[idxmax]);
    best_phi_PVnoDs_ip_vec.push_back(phi_PVnoDs_ip_vec[idxmax]);
    best_phi_PVnoDs_iperr_vec.push_back(phi_PVnoDs_iperr_vec[idxmax]);
    best_phi_PVnoDs_ipchi2_vec.push_back(phi_PVnoDs_ipchi2_vec[idxmax]);
    best_Ds_PVnoDs_ip_vec.push_back(Ds_PVnoDs_ip_vec[idxmax]);
    best_Ds_PVnoDs_iperr_vec.push_back(Ds_PVnoDs_iperr_vec[idxmax]);
    best_Ds_PVnoDs_ipchi2_vec.push_back(Ds_PVnoDs_ipchi2_vec[idxmax]);

    best_PV_withDs_IsValid_vec.push_back(PV_withDs_IsValid_vec[idxmax]);
    best_PV_withDs_IsFake_vec.push_back(PV_withDs_IsFake_vec[idxmax]);
    best_PV_withDs_chi2_vec.push_back(PV_withDs_chi2_vec[idxmax]);
    best_PV_withDs_ndof_vec.push_back(PV_withDs_ndof_vec[idxmax]);
    best_PV_withDs_chi2ndof_vec.push_back(PV_withDs_chi2ndof_vec[idxmax]);
    best_PV_withDs_x_vec.push_back(PV_withDs_x_vec[idxmax]);
    best_PV_withDs_y_vec.push_back(PV_withDs_y_vec[idxmax]);
    best_PV_withDs_z_vec.push_back(PV_withDs_z_vec[idxmax]);
    best_PV_withDs_xerr_vec.push_back(PV_withDs_xerr_vec[idxmax]);
    best_PV_withDs_yerr_vec.push_back(PV_withDs_yerr_vec[idxmax]);
    best_PV_withDs_zerr_vec.push_back(PV_withDs_zerr_vec[idxmax]);

    best_Ds_PVwithDs_FDxy_vec.push_back(Ds_PVwithDs_FDxy_vec[idxmax]);
    best_Ds_PVwithDs_FDz_vec.push_back(Ds_PVwithDs_FDz_vec[idxmax]);
    best_Ds_PVwithDs_FD_vec.push_back(Ds_PVwithDs_FD_vec[idxmax]);
    best_Ds_PVwithDs_FDxyerr_vec.push_back(Ds_PVwithDs_FDxyerr_vec[idxmax]);
    best_Ds_PVwithDs_FDxychi2_vec.push_back(Ds_PVwithDs_FDxychi2_vec[idxmax]);
    best_Ds_PVwithDs_FDzerr_vec.push_back(Ds_PVwithDs_FDzerr_vec[idxmax]);
    best_Ds_PVwithDs_FDzchi2_vec.push_back(Ds_PVwithDs_FDzchi2_vec[idxmax]);
    best_Ds_PVwithDs_FDerr_vec.push_back(Ds_PVwithDs_FDerr_vec[idxmax]);
    best_Ds_PVwithDs_FDchi2_vec.push_back(Ds_PVwithDs_FDchi2_vec[idxmax]);
    best_Ds_PVwithDs_dira_vec.push_back(Ds_PVwithDs_dira_vec[idxmax]);
    best_Ds_PVwithDs_dira_angle_vec.push_back(Ds_PVwithDs_dira_angle_vec[idxmax]);
    best_Kp_PVwithDs_ip_vec.push_back(Kp_PVwithDs_ip_vec[idxmax]);
    best_Kp_PVwithDs_iperr_vec.push_back(Kp_PVwithDs_iperr_vec[idxmax]);
    best_Kp_PVwithDs_ipchi2_vec.push_back(Kp_PVwithDs_ipchi2_vec[idxmax]);
    best_Km_PVwithDs_ip_vec.push_back(Km_PVwithDs_ip_vec[idxmax]);
    best_Km_PVwithDs_iperr_vec.push_back(Km_PVwithDs_iperr_vec[idxmax]);
    best_Km_PVwithDs_ipchi2_vec.push_back(Km_PVwithDs_ipchi2_vec[idxmax]);
    best_pi_PVwithDs_ip_vec.push_back(pi_PVwithDs_ip_vec[idxmax]);
    best_pi_PVwithDs_iperr_vec.push_back(pi_PVwithDs_iperr_vec[idxmax]);
    best_pi_PVwithDs_ipchi2_vec.push_back(pi_PVwithDs_ipchi2_vec[idxmax]);
    best_phi_PVwithDs_ip_vec.push_back(phi_PVwithDs_ip_vec[idxmax]);
    best_phi_PVwithDs_iperr_vec.push_back(phi_PVwithDs_iperr_vec[idxmax]);
    best_phi_PVwithDs_ipchi2_vec.push_back(phi_PVwithDs_ipchi2_vec[idxmax]);
    best_Ds_PVwithDs_ip_vec.push_back(Ds_PVwithDs_ip_vec[idxmax]);
    best_Ds_PVwithDs_iperr_vec.push_back(Ds_PVwithDs_iperr_vec[idxmax]);
    best_Ds_PVwithDs_ipchi2_vec.push_back(Ds_PVwithDs_ipchi2_vec[idxmax]);

    best_match_entry_vec.push_back(match_entry_vec[idxmax]);
}

void PVTree::mu_Reset()
{
    best_mu_charge = null;
    best_mu_eta = null;
    best_mu_phi = null;
    best_mu_vx = null;
    best_mu_vy = null;
    best_mu_vz = null;
    best_mu_p = null;
    best_mu_pt = null;
    best_mu_px = null;
    best_mu_py = null;
    best_mu_pz = null;
    best_mu_isHighPt = false;
    best_mu_isLoose = false;
    best_mu_isMedium = false;
    best_mu_isSoft = false;
    best_mu_isTight = false;
    best_mu_isPF = false;
    best_mu_isTracker = false;
    best_mu_isGlobal = false;
    best_mu_IsoR03_sumChargedHadronPt = null;
    best_mu_IsoR03_sumChargedParticlePt = null;
    best_mu_IsoR03_sumNeutralHadronEt = null;
    best_mu_IsoR03_sumPhotonEt = null;
    best_mu_IsoR03_sumPUPt = null;
    best_mu_PFIsoR03 = null;
    best_mu_IsoR04_sumChargedHadronPt = null;
    best_mu_IsoR04_sumChargedParticlePt = null;
    best_mu_IsoR04_sumNeutralHadronEt = null;
    best_mu_IsoR04_sumPhotonEt = null;
    best_mu_IsoR04_sumPUPt = null;
    best_mu_PFIsoR04 = null;
    best_mu_primvtx_dxy = null;
    best_mu_primvtx_dxyerr = null;
    best_mu_primvtx_dz = null;
    best_mu_primvtx_dzerr = null;
    best_mu_primvtx_ip = null;
    best_mu_primvtx_iperr = null;
    best_mu_primvtx_ipchi2 = null;
    best_mu_match = false;
}


void PVTree::Best_mu_Fill_Vector()
{
    best_mu_charge_vec.push_back(best_mu_charge);
    best_mu_eta_vec.push_back(best_mu_eta);
    best_mu_phi_vec.push_back(best_mu_phi);
    best_mu_vx_vec.push_back(best_mu_vx);
    best_mu_vy_vec.push_back(best_mu_vy);
    best_mu_vz_vec.push_back(best_mu_vz);
    best_mu_p_vec.push_back(best_mu_p);
    best_mu_pt_vec.push_back(best_mu_pt);
    best_mu_px_vec.push_back(best_mu_px);
    best_mu_py_vec.push_back(best_mu_py);
    best_mu_pz_vec.push_back(best_mu_pz);
    best_mu_isHighPt_vec.push_back(best_mu_isHighPt);
    best_mu_isLoose_vec.push_back(best_mu_isLoose);
    best_mu_isMedium_vec.push_back(best_mu_isMedium);
    best_mu_isSoft_vec.push_back(best_mu_isSoft);
    best_mu_isTight_vec.push_back(best_mu_isTight);
    best_mu_isPF_vec.push_back(best_mu_isPF);
    best_mu_isTracker_vec.push_back(best_mu_isTracker);
    best_mu_isGlobal_vec.push_back(best_mu_isGlobal);
    best_mu_IsoR03_sumChargedHadronPt_vec.push_back(best_mu_IsoR03_sumChargedHadronPt);
    best_mu_IsoR03_sumChargedParticlePt_vec.push_back(best_mu_IsoR03_sumChargedParticlePt);
    best_mu_IsoR03_sumNeutralHadronEt_vec.push_back(best_mu_IsoR03_sumNeutralHadronEt);
    best_mu_IsoR03_sumPhotonEt_vec.push_back(best_mu_IsoR03_sumPhotonEt);
    best_mu_IsoR03_sumPUPt_vec.push_back(best_mu_IsoR03_sumPUPt);
    best_mu_PFIsoR03_vec.push_back(best_mu_PFIsoR03);
    best_mu_IsoR04_sumChargedHadronPt_vec.push_back(best_mu_IsoR04_sumChargedHadronPt);
    best_mu_IsoR04_sumChargedParticlePt_vec.push_back(best_mu_IsoR04_sumChargedParticlePt);
    best_mu_IsoR04_sumNeutralHadronEt_vec.push_back(best_mu_IsoR04_sumNeutralHadronEt);
    best_mu_IsoR04_sumPhotonEt_vec.push_back(best_mu_IsoR04_sumPhotonEt);
    best_mu_IsoR04_sumPUPt_vec.push_back(best_mu_IsoR04_sumPUPt);
    best_mu_PFIsoR04_vec.push_back(best_mu_PFIsoR04);
    best_mu_primvtx_dxy_vec.push_back(best_mu_primvtx_dxy);
    best_mu_primvtx_dxyerr_vec.push_back(best_mu_primvtx_dxyerr);
    best_mu_primvtx_dz_vec.push_back(best_mu_primvtx_dz);
    best_mu_primvtx_dzerr_vec.push_back(best_mu_primvtx_dzerr);
    best_mu_primvtx_ip_vec.push_back(best_mu_primvtx_ip);
    best_mu_primvtx_iperr_vec.push_back(best_mu_primvtx_iperr);
    best_mu_primvtx_ipchi2_vec.push_back(best_mu_primvtx_ipchi2);
    best_mu_match_vec.push_back(best_mu_match);
}

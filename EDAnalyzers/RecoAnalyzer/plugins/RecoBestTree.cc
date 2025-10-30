#include "EDAnalyzers/RecoAnalyzer/interface/RecoBestTree.h"
#include <iostream>

RecoBestTree::RecoBestTree(TTree *tree_)
{
    tree = tree_;
}
void RecoBestTree::Init()
{
    /* PV_withBS_IsValid_vec.clear(); */
    /* PV_withBS_IsFake_vec.clear(); */
    /* PV_withBS_chi2_vec.clear(); */
    /* PV_withBS_ndof_vec.clear(); */
    /* PV_withBS_chi2ndof_vec.clear(); */
    /* PV_withBS_X_vec.clear(); */
    /* PV_withBS_Y_vec.clear(); */
    /* PV_withBS_Z_vec.clear(); */
    /* PV_withBS_Xerr_vec.clear(); */
    /* PV_withBS_Yerr_vec.clear(); */
    /* PV_withBS_Zerr_vec.clear(); */
    /* PV_noBS_IsValid_vec.clear(); */
    /* PV_noBS_IsFake_vec.clear(); */
    /* PV_noBS_chi2_vec.clear(); */
    /* PV_noBS_ndof_vec.clear(); */
    /* PV_noBS_chi2ndof_vec.clear(); */
    /* PV_noBS_X_vec.clear(); */
    /* PV_noBS_Y_vec.clear(); */
    /* PV_noBS_Z_vec.clear(); */
    /* PV_noBS_Xerr_vec.clear(); */
    /* PV_noBS_Yerr_vec.clear(); */
    /* PV_noBS_Zerr_vec.clear(); */

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

    Kp_pp_vec.clear();
    Kp_pl_vec.clear();
    Km_pp_vec.clear();
    Km_pl_vec.clear();

    phi_pp_vec.clear();
    phi_pl_vec.clear();
    pi_pp_vec.clear();
    pi_pl_vec.clear();

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

    dxy_Kp_Km_vec.clear();
    dxy_Kp_phi_vec.clear();
    dxy_Km_phi_vec.clear();
    dxy_Kp_pi_vec.clear();
    dxy_Km_pi_vec.clear();
    dxy_pi_phi_vec.clear();
    dxy_Kp_Ds_vec.clear();
    dxy_Km_Ds_vec.clear();
    dxy_phi_Ds_vec.clear();
    dxy_pi_Ds_vec.clear();

    dz_Kp_Km_vec.clear();
    dz_Kp_phi_vec.clear();
    dz_Km_phi_vec.clear();
    dz_Kp_pi_vec.clear();
    dz_Km_pi_vec.clear();
    dz_pi_phi_vec.clear();
    dz_Kp_Ds_vec.clear();
    dz_Km_Ds_vec.clear();
    dz_phi_Ds_vec.clear();
    dz_pi_Ds_vec.clear();

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

    phiFit_Kp_pp_vec.clear();
    phiFit_Kp_pl_vec.clear();
    phiFit_Km_pp_vec.clear();
    phiFit_Km_pl_vec.clear();

    phiFit_phi_pp_vec.clear();
    phiFit_phi_pl_vec.clear();
    phiFit_pi_pp_vec.clear();
    phiFit_pi_pl_vec.clear();

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

    DsFit_Kp_pp_vec.clear();
    DsFit_Kp_pl_vec.clear();
    DsFit_Km_pp_vec.clear();
    DsFit_Km_pl_vec.clear();

    DsFit_phi_pp_vec.clear();
    DsFit_phi_pl_vec.clear();
    DsFit_pi_pp_vec.clear();
    DsFit_pi_pl_vec.clear();

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

    /*     Ds_FDxy_vec.clear(); */
    /*     Ds_FDxy_Err_vec.clear(); */
    /*     Ds_FDxy_Chi2_vec.clear(); */
    /*     Ds_FDz_vec.clear(); */
    /*     Ds_FDz_Err_vec.clear(); */
    /*     Ds_FDz_Chi2_vec.clear(); */
    /*     Ds_FD_vec.clear(); */
    /*     Ds_FD_Err_vec.clear(); */
    /*     Ds_FD_Chi2_vec.clear(); */
    /*     Ds_DIRA_angle_vec.clear(); */
    /*     Ds_DIRA_vec.clear(); */
    /*     Kp_IP_vec.clear(); */
    /*     Kp_IP_Err_vec.clear(); */
    /*     Kp_IP_Chi2_vec.clear(); */
    /*     Km_IP_vec.clear(); */
    /*     Km_IP_Err_vec.clear(); */
    /*     Km_IP_Chi2_vec.clear(); */
    /*     pi_IP_vec.clear(); */
    /*     pi_IP_Err_vec.clear(); */
    /*     pi_IP_Chi2_vec.clear(); */

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

    best_Kp_pp_vec.clear();
    best_Kp_pl_vec.clear();
    best_Km_pp_vec.clear();
    best_Km_pl_vec.clear();

    best_phi_pp_vec.clear();
    best_phi_pl_vec.clear();
    best_pi_pp_vec.clear();
    best_pi_pl_vec.clear();

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

    best_dxy_Kp_Km_vec.clear();
    best_dxy_Kp_phi_vec.clear();
    best_dxy_Km_phi_vec.clear();
    best_dxy_Kp_pi_vec.clear();
    best_dxy_Km_pi_vec.clear();
    best_dxy_pi_phi_vec.clear();
    best_dxy_Kp_Ds_vec.clear();
    best_dxy_Km_Ds_vec.clear();
    best_dxy_phi_Ds_vec.clear();
    best_dxy_pi_Ds_vec.clear();

    best_dz_Kp_Km_vec.clear();
    best_dz_Kp_phi_vec.clear();
    best_dz_Km_phi_vec.clear();
    best_dz_Kp_pi_vec.clear();
    best_dz_Km_pi_vec.clear();
    best_dz_pi_phi_vec.clear();
    best_dz_Kp_Ds_vec.clear();
    best_dz_Km_Ds_vec.clear();
    best_dz_phi_Ds_vec.clear();
    best_dz_pi_Ds_vec.clear();

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

    best_phiFit_Kp_pp_vec.clear();
    best_phiFit_Kp_pl_vec.clear();
    best_phiFit_Km_pp_vec.clear();
    best_phiFit_Km_pl_vec.clear();

    best_phiFit_phi_pp_vec.clear();
    best_phiFit_phi_pl_vec.clear();
    best_phiFit_pi_pp_vec.clear();
    best_phiFit_pi_pl_vec.clear();

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

    best_DsFit_Kp_pp_vec.clear();
    best_DsFit_Kp_pl_vec.clear();
    best_DsFit_Km_pp_vec.clear();
    best_DsFit_Km_pl_vec.clear();

    best_DsFit_phi_pp_vec.clear();
    best_DsFit_phi_pl_vec.clear();
    best_DsFit_pi_pp_vec.clear();
    best_DsFit_pi_pl_vec.clear();

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

    /* best_Ds_FDxy_vec.clear(); */
    /* best_Ds_FDxy_Err_vec.clear(); */
    /* best_Ds_FDxy_Chi2_vec.clear(); */
    /* best_Ds_FDz_vec.clear(); */
    /* best_Ds_FDz_Err_vec.clear(); */
    /* best_Ds_FDz_Chi2_vec.clear(); */
    /* best_Ds_FD_vec.clear(); */
    /* best_Ds_FD_Err_vec.clear(); */
    /* best_Ds_FD_Chi2_vec.clear(); */
    /* best_Ds_DIRA_angle_vec.clear(); */
    /* best_Ds_DIRA_vec.clear(); */
    /* best_Kp_IP_vec.clear(); */
    /* best_Kp_IP_Err_vec.clear(); */
    /* best_Kp_IP_Chi2_vec.clear(); */
    /* best_Km_IP_vec.clear(); */
    /* best_Km_IP_Err_vec.clear(); */
    /* best_Km_IP_Chi2_vec.clear(); */
    /* best_pi_IP_vec.clear(); */
    /* best_pi_IP_Err_vec.clear(); */
    /* best_pi_IP_Chi2_vec.clear(); */

}

void RecoBestTree::Best_Fill_Vector(int idxmax)
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

    best_Kp_pp_vec.push_back(Kp_pp_vec[idxmax]);
    best_Kp_pl_vec.push_back(Kp_pl_vec[idxmax]);
    best_Km_pp_vec.push_back(Km_pp_vec[idxmax]);
    best_Km_pl_vec.push_back(Km_pl_vec[idxmax]);

    best_phi_pp_vec.push_back(phi_pp_vec[idxmax]);
    best_phi_pl_vec.push_back(phi_pl_vec[idxmax]);
    best_pi_pp_vec.push_back(pi_pp_vec[idxmax]);
    best_pi_pl_vec.push_back(pi_pl_vec[idxmax]);

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

    best_dxy_Kp_Km_vec.push_back(dxy_Kp_Km_vec[idxmax]);
    best_dxy_Kp_phi_vec.push_back(dxy_Kp_phi_vec[idxmax]);
    best_dxy_Km_phi_vec.push_back(dxy_Km_phi_vec[idxmax]);
    best_dxy_Kp_pi_vec.push_back(dxy_Kp_pi_vec[idxmax]);
    best_dxy_Km_pi_vec.push_back(dxy_Km_pi_vec[idxmax]);
    best_dxy_pi_phi_vec.push_back(dxy_pi_phi_vec[idxmax]);
    best_dxy_Kp_Ds_vec.push_back(dxy_Kp_Ds_vec[idxmax]);
    best_dxy_Km_Ds_vec.push_back(dxy_Km_Ds_vec[idxmax]);
    best_dxy_phi_Ds_vec.push_back(dxy_phi_Ds_vec[idxmax]);
    best_dxy_pi_Ds_vec.push_back(dxy_pi_Ds_vec[idxmax]);

    best_dz_Kp_Km_vec.push_back(dz_Kp_Km_vec[idxmax]);
    best_dz_Kp_phi_vec.push_back(dz_Kp_phi_vec[idxmax]);
    best_dz_Km_phi_vec.push_back(dz_Km_phi_vec[idxmax]);
    best_dz_Kp_pi_vec.push_back(dz_Kp_pi_vec[idxmax]);
    best_dz_Km_pi_vec.push_back(dz_Km_pi_vec[idxmax]);
    best_dz_pi_phi_vec.push_back(dz_pi_phi_vec[idxmax]);
    best_dz_Kp_Ds_vec.push_back(dz_Kp_Ds_vec[idxmax]);
    best_dz_Km_Ds_vec.push_back(dz_Km_Ds_vec[idxmax]);
    best_dz_phi_Ds_vec.push_back(dz_phi_Ds_vec[idxmax]);
    best_dz_pi_Ds_vec.push_back(dz_pi_Ds_vec[idxmax]);


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

    best_phiFit_Kp_pp_vec.push_back(phiFit_Kp_pp_vec[idxmax]);
    best_phiFit_Kp_pl_vec.push_back(phiFit_Kp_pl_vec[idxmax]);
    best_phiFit_Km_pp_vec.push_back(phiFit_Km_pp_vec[idxmax]);
    best_phiFit_Km_pl_vec.push_back(phiFit_Km_pl_vec[idxmax]);

    best_phiFit_phi_pp_vec.push_back(phiFit_phi_pp_vec[idxmax]);
    best_phiFit_phi_pl_vec.push_back(phiFit_phi_pl_vec[idxmax]);
    best_phiFit_pi_pp_vec.push_back(phiFit_pi_pp_vec[idxmax]);
    best_phiFit_pi_pl_vec.push_back(phiFit_pi_pl_vec[idxmax]);

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

    best_DsFit_Kp_pp_vec.push_back(DsFit_Kp_pp_vec[idxmax]);
    best_DsFit_Kp_pl_vec.push_back(DsFit_Kp_pl_vec[idxmax]);
    best_DsFit_Km_pp_vec.push_back(DsFit_Km_pp_vec[idxmax]);
    best_DsFit_Km_pl_vec.push_back(DsFit_Km_pl_vec[idxmax]);

    best_DsFit_phi_pp_vec.push_back(DsFit_phi_pp_vec[idxmax]);
    best_DsFit_phi_pl_vec.push_back(DsFit_phi_pl_vec[idxmax]);
    best_DsFit_pi_pp_vec.push_back(DsFit_pi_pp_vec[idxmax]);
    best_DsFit_pi_pl_vec.push_back(DsFit_pi_pl_vec[idxmax]);

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

    /* best_Ds_FDxy_vec.push_back(Ds_FDxy_vec[idxmax]); */
    /* best_Ds_FDxy_Err_vec.push_back(Ds_FDxy_Err_vec[idxmax]); */
    /* best_Ds_FDxy_Chi2_vec.push_back(Ds_FDxy_Chi2_vec[idxmax]); */
    /* best_Ds_FDz_vec.push_back(Ds_FDz_vec[idxmax]); */
    /* best_Ds_FDz_Err_vec.push_back(Ds_FDz_Err_vec[idxmax]); */
    /* best_Ds_FDz_Chi2_vec.push_back(Ds_FDz_Chi2_vec[idxmax]); */
    /* best_Ds_FD_vec.push_back(Ds_FD_vec[idxmax]); */
    /* best_Ds_FD_Err_vec.push_back(Ds_FD_Err_vec[idxmax]); */
    /* best_Ds_FD_Chi2_vec.push_back(Ds_FD_Chi2_vec[idxmax]); */
    /* best_Ds_DIRA_angle_vec.push_back(Ds_DIRA_angle_vec[idxmax]); */
    /* best_Ds_DIRA_vec.push_back(Ds_DIRA_vec[idxmax]); */
    /* best_Kp_IP_vec.push_back(Kp_IP_vec[idxmax]); */
    /* best_Kp_IP_Err_vec.push_back(Kp_IP_Err_vec[idxmax]); */
    /* best_Kp_IP_Chi2_vec.push_back(Kp_IP_Chi2_vec[idxmax]); */
    /* best_Km_IP_vec.push_back(Km_IP_vec[idxmax]); */
    /* best_Km_IP_Err_vec.push_back(Km_IP_Err_vec[idxmax]); */
    /* best_Km_IP_Chi2_vec.push_back(Km_IP_Chi2_vec[idxmax]); */
    /* best_pi_IP_vec.push_back(pi_IP_vec[idxmax]); */
    /* best_pi_IP_Err_vec.push_back(pi_IP_Err_vec[idxmax]); */
    /* best_pi_IP_Chi2_vec.push_back(pi_IP_Chi2_vec[idxmax]); */
}

void RecoBestTree::Fill_Vector()
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

    Kp_pp_vec.push_back(Kp_pp);
    Kp_pl_vec.push_back(Kp_pl);
    Km_pp_vec.push_back(Km_pp);
    Km_pl_vec.push_back(Km_pl);

    phi_pp_vec.push_back(phi_pp);
    phi_pl_vec.push_back(phi_pl);
    pi_pp_vec.push_back(pi_pp);
    pi_pl_vec.push_back(pi_pl);

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

    dxy_Kp_Km_vec.push_back(dxy_Kp_Km);
    dxy_Kp_phi_vec.push_back(dxy_Kp_phi);
    dxy_Km_phi_vec.push_back(dxy_Km_phi);
    dxy_Kp_pi_vec.push_back(dxy_Kp_pi);
    dxy_Km_pi_vec.push_back(dxy_Km_pi);
    dxy_pi_phi_vec.push_back(dxy_pi_phi);
    dxy_Kp_Ds_vec.push_back(dxy_Kp_Ds);
    dxy_Km_Ds_vec.push_back(dxy_Km_Ds);
    dxy_phi_Ds_vec.push_back(dxy_phi_Ds);
    dxy_pi_Ds_vec.push_back(dxy_pi_Ds);

    dz_Kp_Km_vec.push_back(dz_Kp_Km);
    dz_Kp_phi_vec.push_back(dz_Kp_phi);
    dz_Km_phi_vec.push_back(dz_Km_phi);
    dz_Kp_pi_vec.push_back(dz_Kp_pi);
    dz_Km_pi_vec.push_back(dz_Km_pi);
    dz_pi_phi_vec.push_back(dz_pi_phi);
    dz_Kp_Ds_vec.push_back(dz_Kp_Ds);
    dz_Km_Ds_vec.push_back(dz_Km_Ds);
    dz_phi_Ds_vec.push_back(dz_phi_Ds);
    dz_pi_Ds_vec.push_back(dz_pi_Ds);


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

    phiFit_Kp_pp_vec.push_back(phiFit_Kp_pp);
    phiFit_Kp_pl_vec.push_back(phiFit_Kp_pl);
    phiFit_Km_pp_vec.push_back(phiFit_Km_pp);
    phiFit_Km_pl_vec.push_back(phiFit_Km_pl);

    phiFit_phi_pp_vec.push_back(phiFit_phi_pp);
    phiFit_phi_pl_vec.push_back(phiFit_phi_pl);
    phiFit_pi_pp_vec.push_back(phiFit_pi_pp);
    phiFit_pi_pl_vec.push_back(phiFit_pi_pl);

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

    DsFit_Kp_pp_vec.push_back(DsFit_Kp_pp);
    DsFit_Kp_pl_vec.push_back(DsFit_Kp_pl);
    DsFit_Km_pp_vec.push_back(DsFit_Km_pp);
    DsFit_Km_pl_vec.push_back(DsFit_Km_pl);

    DsFit_phi_pp_vec.push_back(DsFit_phi_pp);
    DsFit_phi_pl_vec.push_back(DsFit_phi_pl);
    DsFit_pi_pp_vec.push_back(DsFit_pi_pp);
    DsFit_pi_pl_vec.push_back(DsFit_pi_pl);

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

    /* Ds_FDxy_vec.push_back(Ds_FDxy); */
    /* Ds_FDxy_Err_vec.push_back(Ds_FDxy_Err); */
    /* Ds_FDxy_Chi2_vec.push_back(Ds_FDxy_Chi2); */
    /* Ds_FDz_vec.push_back(Ds_FDz); */
    /* Ds_FDz_Err_vec.push_back(Ds_FDz_Err); */
    /* Ds_FDz_Chi2_vec.push_back(Ds_FDz_Chi2); */
    /* Ds_FD_vec.push_back(Ds_FD); */
    /* Ds_FD_Err_vec.push_back(Ds_FD_Err); */
    /* Ds_FD_Chi2_vec.push_back(Ds_FD_Chi2); */
    /* Ds_DIRA_angle_vec.push_back(Ds_DIRA_angle); */
    /* Ds_DIRA_vec.push_back(Ds_DIRA); */
    /* Kp_IP_vec.push_back(Kp_IP); */
    /* Kp_IP_Err_vec.push_back(Kp_IP_Err); */
    /* Kp_IP_Chi2_vec.push_back(Kp_IP_Chi2); */
    /* Km_IP_vec.push_back(Km_IP); */
    /* Km_IP_Err_vec.push_back(Km_IP_Err); */
    /* Km_IP_Chi2_vec.push_back(Km_IP_Chi2); */
    /* pi_IP_vec.push_back(pi_IP); */
    /* pi_IP_Err_vec.push_back(pi_IP_Err); */
    /* pi_IP_Chi2_vec.push_back(pi_IP_Chi2); */
}

/* void RecoBestTree::PV_withBS_Fill_Vector() */
/* { */
/*     PV_withBS_IsValid_vec.push_back(PV_withBS_IsValid); */
/*     PV_withBS_IsFake_vec.push_back(PV_withBS_IsFake); */
/*     PV_withBS_chi2_vec.push_back(PV_withBS_chi2); */
/*     PV_withBS_ndof_vec.push_back(PV_withBS_ndof); */
/*     PV_withBS_chi2ndof_vec.push_back(PV_withBS_chi2ndof); */
/*     PV_withBS_X_vec.push_back(PV_withBS_X); */
/*     PV_withBS_Y_vec.push_back(PV_withBS_Y); */
/*     PV_withBS_Z_vec.push_back(PV_withBS_Z); */
/*     PV_withBS_Xerr_vec.push_back(PV_withBS_Xerr); */
/*     PV_withBS_Yerr_vec.push_back(PV_withBS_Yerr); */
/*     PV_withBS_Zerr_vec.push_back(PV_withBS_Zerr); */
/* } */

/* void RecoBestTree::PV_noBS_Fill_Vector() */
/* { */
/*     PV_noBS_IsValid_vec.push_back(PV_noBS_IsValid); */
/*     PV_noBS_IsFake_vec.push_back(PV_noBS_IsFake); */
/*     PV_noBS_chi2_vec.push_back(PV_noBS_chi2); */
/*     PV_noBS_ndof_vec.push_back(PV_noBS_ndof); */
/*     PV_noBS_chi2ndof_vec.push_back(PV_noBS_chi2ndof); */
/*     PV_noBS_X_vec.push_back(PV_noBS_X); */
/*     PV_noBS_Y_vec.push_back(PV_noBS_Y); */
/*     PV_noBS_Z_vec.push_back(PV_noBS_Z); */
/*     PV_noBS_Xerr_vec.push_back(PV_noBS_Xerr); */
/*     PV_noBS_Yerr_vec.push_back(PV_noBS_Yerr); */
/*     PV_noBS_Zerr_vec.push_back(PV_noBS_Zerr); */
/* } */

void RecoBestTree::CreateBranches()
{
    /* tree->Branch("BS_type", &BS_type); */ 
    /* tree->Branch("BS_X0", &BS_X0); */
    /* tree->Branch("BS_Y0", &BS_Y0); */
    /* tree->Branch("BS_Z0", &BS_Z0); */
    /* tree->Branch("BS_SigmaZ", &BS_SigmaZ); */
    /* tree->Branch("BS_dXdZ", &BS_dXdZ); */
    /* tree->Branch("BS_dYdZ", &BS_dYdZ); */
    /* tree->Branch("BS_BWX", &BS_BWX); */
    /* tree->Branch("BS_BWY", &BS_BWY); */
    /* tree->Branch("BS_X0err", &BS_X0err); */
    /* tree->Branch("BS_Y0err", &BS_Y0err); */
    /* tree->Branch("BS_Z0err", &BS_Z0err); */
    /* tree->Branch("BS_SigmaZ0err", &BS_SigmaZ0err); */
    /* tree->Branch("BS_dXdZerr", &BS_dXdZerr); */
    /* tree->Branch("BS_dYdZerr", &BS_dYdZerr); */
    /* tree->Branch("BS_BWXerr", &BS_BWXerr); */
    /* tree->Branch("BS_BWYerr", &BS_BWYerr); */
    /* tree->Branch("BS_EmitX", &BS_EmitX); */
    /* tree->Branch("BS_EmitY", &BS_EmitY); */
    /* tree->Branch("BS_BetaStar", &BS_BetaStar); */

    /* tree->Branch("PV_withBS_IsValid", &PV_withBS_IsValid_vec); */
    /* tree->Branch("PV_withBS_IsFake", &PV_withBS_IsFake_vec); */
    /* tree->Branch("PV_withBS_chi2", &PV_withBS_chi2_vec); */
    /* tree->Branch("PV_withBS_ndof", &PV_withBS_ndof_vec); */
    /* tree->Branch("PV_withBS_chi2ndof", &PV_withBS_chi2ndof_vec); */
    /* tree->Branch("PV_withBS_X", &PV_withBS_X_vec); */
    /* tree->Branch("PV_withBS_Y", &PV_withBS_Y_vec); */
    /* tree->Branch("PV_withBS_Z", &PV_withBS_Z_vec); */
    /* tree->Branch("PV_withBS_Xerr", &PV_withBS_Xerr_vec); */
    /* tree->Branch("PV_withBS_Yerr", &PV_withBS_Yerr_vec); */
    /* tree->Branch("PV_withBS_Zerr", &PV_withBS_Zerr_vec); */
    /* tree->Branch("PV_noBS_IsValid", &PV_noBS_IsValid_vec); */
    /* tree->Branch("PV_noBS_IsFake", &PV_noBS_IsFake_vec); */
    /* tree->Branch("PV_noBS_chi2", &PV_noBS_chi2_vec); */
    /* tree->Branch("PV_noBS_ndof", &PV_noBS_ndof_vec); */
    /* tree->Branch("PV_noBS_chi2ndof", &PV_noBS_chi2ndof_vec); */
    /* tree->Branch("PV_noBS_X", &PV_noBS_X_vec); */
    /* tree->Branch("PV_noBS_Y", &PV_noBS_Y_vec); */
    /* tree->Branch("PV_noBS_Z", &PV_noBS_Z_vec); */
    /* tree->Branch("PV_noBS_Xerr", &PV_noBS_Xerr_vec); */
    /* tree->Branch("PV_noBS_Yerr", &PV_noBS_Yerr_vec); */
    /* tree->Branch("PV_noBS_Zerr", &PV_noBS_Zerr_vec); */

    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

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

    tree->Branch("best_Kp_pp", &best_Kp_pp_vec);
    tree->Branch("best_Kp_pl", &best_Kp_pl_vec);
    tree->Branch("best_Km_pp", &best_Km_pp_vec);
    tree->Branch("best_Km_pl", &best_Km_pl_vec);

    tree->Branch("best_phi_pp", &best_phi_pp_vec);
    tree->Branch("best_phi_pl", &best_phi_pl_vec);
    tree->Branch("best_pi_pp", &best_pi_pp_vec);
    tree->Branch("best_pi_pl", &best_pi_pl_vec);

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

    tree->Branch("best_dxy_Kp_Km", &best_dxy_Kp_Km_vec);
    tree->Branch("best_dxy_Kp_phi", &best_dxy_Kp_phi_vec);
    tree->Branch("best_dxy_Km_phi", &best_dxy_Km_phi_vec);
    tree->Branch("best_dxy_Kp_pi", &best_dxy_Kp_pi_vec);
    tree->Branch("best_dxy_Km_pi", &best_dxy_Km_pi_vec);
    tree->Branch("best_dxy_pi_phi", &best_dxy_pi_phi_vec);
    tree->Branch("best_dxy_Kp_Ds", &best_dxy_Kp_Ds_vec);
    tree->Branch("best_dxy_Km_Ds", &best_dxy_Km_Ds_vec);
    tree->Branch("best_dxy_phi_Ds", &best_dxy_phi_Ds_vec);
    tree->Branch("best_dxy_pi_Ds", &best_dxy_pi_Ds_vec);

    tree->Branch("best_dz_Kp_Km", &best_dz_Kp_Km_vec);
    tree->Branch("best_dz_Kp_phi", &best_dz_Kp_phi_vec);
    tree->Branch("best_dz_Km_phi", &best_dz_Km_phi_vec);
    tree->Branch("best_dz_Kp_pi", &best_dz_Kp_pi_vec);
    tree->Branch("best_dz_Km_pi", &best_dz_Km_pi_vec);
    tree->Branch("best_dz_pi_phi", &best_dz_pi_phi_vec);
    tree->Branch("best_dz_Kp_Ds", &best_dz_Kp_Ds_vec);
    tree->Branch("best_dz_Km_Ds", &best_dz_Km_Ds_vec);
    tree->Branch("best_dz_phi_Ds", &best_dz_phi_Ds_vec);
    tree->Branch("best_dz_pi_Ds", &best_dz_pi_Ds_vec);

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

    tree->Branch("best_phiFit_Kp_pp", &best_phiFit_Kp_pp_vec);
    tree->Branch("best_phiFit_Kp_pl", &best_phiFit_Kp_pl_vec);
    tree->Branch("best_phiFit_Km_pp", &best_phiFit_Km_pp_vec);
    tree->Branch("best_phiFit_Km_pl", &best_phiFit_Km_pl_vec);

    tree->Branch("best_phiFit_phi_pp", &best_phiFit_phi_pp_vec);
    tree->Branch("best_phiFit_phi_pl", &best_phiFit_phi_pl_vec);
    tree->Branch("best_phiFit_pi_pp", &best_phiFit_pi_pp_vec);
    tree->Branch("best_phiFit_pi_pl", &best_phiFit_pi_pl_vec);

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

    tree->Branch("best_DsFit_Kp_pp", &best_DsFit_Kp_pp_vec);
    tree->Branch("best_DsFit_Kp_pl", &best_DsFit_Kp_pl_vec);
    tree->Branch("best_DsFit_Km_pp", &best_DsFit_Km_pp_vec);
    tree->Branch("best_DsFit_Km_pl", &best_DsFit_Km_pl_vec);

    tree->Branch("best_DsFit_phi_pp", &best_DsFit_phi_pp_vec);
    tree->Branch("best_DsFit_phi_pl", &best_DsFit_phi_pl_vec);
    tree->Branch("best_DsFit_pi_pp", &best_DsFit_pi_pp_vec);
    tree->Branch("best_DsFit_pi_pl", &best_DsFit_pi_pl_vec);

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

    /* tree->Branch("best_Ds_FDxy", &best_Ds_FDxy_vec); */
    /* tree->Branch("best_Ds_FDxy_Err", &best_Ds_FDxy_Err_vec); */
    /* tree->Branch("best_Ds_FDxy_Chi2", &best_Ds_FDxy_Chi2_vec); */
    /* tree->Branch("best_Ds_FDz", &best_Ds_FDz_vec); */
    /* tree->Branch("best_Ds_FDz_Err", &best_Ds_FDz_Err_vec); */
    /* tree->Branch("best_Ds_FDz_Chi2", &best_Ds_FDz_Chi2_vec); */
    /* tree->Branch("best_Ds_FD", &best_Ds_FD_vec); */
    /* tree->Branch("best_Ds_FD_Err", &best_Ds_FD_Err_vec); */
    /* tree->Branch("best_Ds_FD_Chi2", &best_Ds_FD_Chi2_vec); */
    /* tree->Branch("best_Ds_DIRA_angle", &best_Ds_DIRA_angle_vec); */
    /* tree->Branch("best_Ds_DIRA", &best_Ds_DIRA_vec); */
    /* tree->Branch("best_Kp_IP", &best_Kp_IP_vec); */
    /* tree->Branch("best_Kp_IP_Err", &best_Kp_IP_Err_vec); */
    /* tree->Branch("best_Kp_IP_Chi2", &best_Kp_IP_Chi2_vec); */
    /* tree->Branch("best_Km_IP", &best_Km_IP_vec); */
    /* tree->Branch("best_Km_IP_Err", &best_Km_IP_Err_vec); */
    /* tree->Branch("best_Km_IP_Chi2", &best_Km_IP_Chi2_vec); */
    /* tree->Branch("best_pi_IP", &best_pi_IP_vec); */
    /* tree->Branch("best_pi_IP_Err", &best_pi_IP_Err_vec); */
    /* tree->Branch("best_pi_IP_Chi2", &best_pi_IP_Chi2_vec); */
}

void RecoBestTree::Kp_Reset()
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

void RecoBestTree::Km_Reset()
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

    dxy_Kp_Km = null;
    dz_Kp_Km = null;

    Kp_pp = null;
    Kp_pl = null;
    Km_pp = null;
    Km_pl = null;
}

void RecoBestTree::phi_Reset()
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

    phiFit_Kp_pp = null;
    phiFit_Kp_pl = null;
    phiFit_Km_pp = null;
    phiFit_Km_pl = null;

    phiFit_dR_Kp_Km = null;
    phiFit_dR_Kp_phi = null;
    phiFit_dR_Km_phi = null;

    dxy_Kp_phi = null;
    dxy_Km_phi = null;
    dz_Kp_phi = null;
    dz_Km_phi = null;
}

void RecoBestTree::pi_Reset()
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

    phi_pp = null;
    phi_pl = null;
    pi_pp = null;
    pi_pl = null;

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

    phiFit_phi_pp = null;
    phiFit_phi_pl = null;
    phiFit_pi_pp = null;
    phiFit_pi_pl = null;

    phiFit_dR_Kp_pi = null;
    phiFit_dR_Km_pi = null;
    phiFit_dR_pi_phi = null;
    phiFit_dR_Kp_Ds = null;
    phiFit_dR_Km_Ds = null;
    phiFit_dR_pi_Ds = null;
    phiFit_dR_phi_Ds = null;

    dxy_Kp_pi = null;
    dxy_Km_pi = null;
    dxy_pi_phi = null;
    dz_Kp_pi = null;
    dz_Km_pi = null;
    dz_pi_phi = null;
}

void RecoBestTree::Ds_Reset()
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

    DsFit_Kp_pp = null;
    DsFit_Kp_pl = null;
    DsFit_Km_pp = null;
    DsFit_Km_pl = null;

    DsFit_phi_pp = null;
    DsFit_phi_pl = null;
    DsFit_pi_pp = null;
    DsFit_pi_pl = null;

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

    dxy_Kp_Ds = null;
    dxy_Km_Ds = null;
    dxy_pi_Ds = null;
    dxy_phi_Ds = null;
    dz_Kp_Ds = null;
    dz_Km_Ds = null;
    dz_pi_Ds = null;
    dz_phi_Ds = null;

    /* Ds_FDxy = null; */
    /* Ds_FDxy_Err = null; */
    /* Ds_FDxy_Chi2 = null; */
    /* Ds_FDz = null; */
    /* Ds_FDz_Err = null; */
    /* Ds_FDz_Chi2 = null; */
    /* Ds_FD = null; */
    /* Ds_FD_Err = null; */
    /* Ds_FD_Chi2 = null; */
    /* Ds_DIRA_angle = null; */
    /* Ds_DIRA = null; */
    /* Kp_IP = null; */
    /* Kp_IP_Err = null; */
    /* Kp_IP_Chi2 = null; */
    /* Km_IP = null; */
    /* Km_IP_Err = null; */
    /* Km_IP_Chi2 = null; */
    /* pi_IP = null; */
    /* pi_IP_Err = null; */
    /* pi_IP_Chi2 = null; */
}

/* void RecoBestTree::BS_Reset() */
/* { */
/*     BS_type = null; */ 
/*     BS_X0 = null; */
/*     BS_Y0 = null; */
/*     BS_Z0 = null; */
/*     BS_SigmaZ = null; */
/*     BS_dXdZ = null; */
/*     BS_dYdZ = null; */
/*     BS_BWX = null; */
/*     BS_BWY = null; */
/*     BS_X0err = null; */
/*     BS_Y0err = null; */
/*     BS_Z0err = null; */
/*     BS_SigmaZ0err = null; */
/*     BS_dXdZerr = null; */
/*     BS_dYdZerr = null; */
/*     BS_BWXerr = null; */
/*     BS_BWYerr = null; */
/*     BS_EmitX = null; */
/*     BS_EmitY = null; */
/*     BS_BetaStar = null; */
/* } */

/* void RecoBestTree::PV_Reset() */
/* { */
/*     PV_withBS_IsValid = false; */
/*     PV_withBS_IsFake = false; */
/*     PV_withBS_chi2 = null; */
/*     PV_withBS_ndof = null; */
/*     PV_withBS_chi2ndof = null; */
/*     PV_withBS_X = null; */
/*     PV_withBS_Y = null; */
/*     PV_withBS_Z = null; */
/*     PV_withBS_Xerr = null; */
/*     PV_withBS_Yerr = null; */
/*     PV_withBS_Zerr = null; */
/*     PV_noBS_IsValid = false; */
/*     PV_noBS_IsFake = false; */
/*     PV_noBS_chi2 = null; */
/*     PV_noBS_ndof = null; */
/*     PV_noBS_chi2ndof = null; */
/*     PV_noBS_X = null; */
/*     PV_noBS_Y = null; */
/*     PV_noBS_Z = null; */
/*     PV_noBS_Xerr = null; */
/*     PV_noBS_Yerr = null; */
/*     PV_noBS_Zerr = null; */
/* } */

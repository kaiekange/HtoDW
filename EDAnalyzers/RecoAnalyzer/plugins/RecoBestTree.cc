#include "EDAnalyzers/RecoAnalyzer/interface/RecoBestTree.h"
#include <iostream>

RecoBestTree::RecoBestTree(TTree *tree_)
{
    tree = tree_;
}
void RecoBestTree::Init()
{
    PV_withBS_IsValid_vec.clear();
    PV_withBS_IsFake_vec.clear();
    PV_withBS_CHI2_vec.clear();
    PV_withBS_NDOF_vec.clear();
    PV_withBS_CHI2NDOF_vec.clear();
    PV_withBS_X_vec.clear();
    PV_withBS_Y_vec.clear();
    PV_withBS_Z_vec.clear();
    PV_withBS_XERR_vec.clear();
    PV_withBS_YERR_vec.clear();
    PV_withBS_ZERR_vec.clear();
    PV_noBS_IsValid_vec.clear();
    PV_noBS_IsFake_vec.clear();
    PV_noBS_CHI2_vec.clear();
    PV_noBS_NDOF_vec.clear();
    PV_noBS_CHI2NDOF_vec.clear();
    PV_noBS_X_vec.clear();
    PV_noBS_Y_vec.clear();
    PV_noBS_Z_vec.clear();
    PV_noBS_XERR_vec.clear();
    PV_noBS_YERR_vec.clear();
    PV_noBS_ZERR_vec.clear();
    
    num_reco_phi = 0;
    num_reco_Ds = 0;

    Kp_ETA_vec.clear();
    Kp_PHI_vec.clear();
    Kp_ORIVX_X_vec.clear();
    Kp_ORIVX_Y_vec.clear();
    Kp_ORIVX_Z_vec.clear();
    Kp_P_vec.clear();
    Kp_PT_vec.clear();
    Kp_PX_vec.clear();
    Kp_PY_vec.clear();
    Kp_PZ_vec.clear();

    Km_ETA_vec.clear();
    Km_PHI_vec.clear();
    Km_ORIVX_X_vec.clear();
    Km_ORIVX_Y_vec.clear();
    Km_ORIVX_Z_vec.clear();
    Km_P_vec.clear();
    Km_PT_vec.clear();
    Km_PX_vec.clear();
    Km_PY_vec.clear();
    Km_PZ_vec.clear();

    pi_ETA_vec.clear();
    pi_PHI_vec.clear();
    pi_ORIVX_X_vec.clear();
    pi_ORIVX_Y_vec.clear();
    pi_ORIVX_Z_vec.clear();
    pi_P_vec.clear();
    pi_PT_vec.clear();
    pi_PX_vec.clear();
    pi_PY_vec.clear();
    pi_PZ_vec.clear();

    phi_ETA_vec.clear();
    phi_PHI_vec.clear();
    phi_P_vec.clear();
    phi_PT_vec.clear();
    phi_PX_vec.clear();
    phi_PY_vec.clear();
    phi_PZ_vec.clear();
    phi_M_vec.clear();

    Ds_ETA_vec.clear();
    Ds_PHI_vec.clear();
    Ds_P_vec.clear();
    Ds_PT_vec.clear();
    Ds_PX_vec.clear();
    Ds_PY_vec.clear();
    Ds_PZ_vec.clear();
    Ds_M_vec.clear();

    Kp_PP_vec.clear();
    Kp_PL_vec.clear();
    Km_PP_vec.clear();
    Km_PL_vec.clear();

    phi_PP_vec.clear();
    phi_PL_vec.clear();
    pi_PP_vec.clear();
    pi_PL_vec.clear();

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
    phiFit_CHI2_vec.clear();
    phiFit_NDOF_vec.clear();
    phiFit_CHI2NDOF_vec.clear();
    phiFit_ENDVX_X_vec.clear();
    phiFit_ENDVX_Y_vec.clear();
    phiFit_ENDVX_Z_vec.clear();
    phiFit_ENDVX_XERR_vec.clear();
    phiFit_ENDVX_YERR_vec.clear();
    phiFit_ENDVX_ZERR_vec.clear();

    phiFit_Kp_ETA_vec.clear();
    phiFit_Kp_PHI_vec.clear();
    phiFit_Kp_P_vec.clear();
    phiFit_Kp_PT_vec.clear();
    phiFit_Kp_PX_vec.clear();
    phiFit_Kp_PY_vec.clear();
    phiFit_Kp_PZ_vec.clear();

    phiFit_Km_ETA_vec.clear();
    phiFit_Km_PHI_vec.clear();
    phiFit_Km_P_vec.clear();
    phiFit_Km_PT_vec.clear();
    phiFit_Km_PX_vec.clear();
    phiFit_Km_PY_vec.clear();
    phiFit_Km_PZ_vec.clear();

    phiFit_pi_ETA_vec.clear();
    phiFit_pi_PHI_vec.clear();
    phiFit_pi_P_vec.clear();
    phiFit_pi_PT_vec.clear();
    phiFit_pi_PX_vec.clear();
    phiFit_pi_PY_vec.clear();
    phiFit_pi_PZ_vec.clear();

    phiFit_phi_ETA_vec.clear();
    phiFit_phi_PHI_vec.clear();
    phiFit_phi_P_vec.clear();
    phiFit_phi_PT_vec.clear();
    phiFit_phi_PX_vec.clear();
    phiFit_phi_PY_vec.clear();
    phiFit_phi_PZ_vec.clear();
    phiFit_phi_M_vec.clear();

    phiFit_Ds_ETA_vec.clear();
    phiFit_Ds_PHI_vec.clear();
    phiFit_Ds_P_vec.clear();
    phiFit_Ds_PT_vec.clear();
    phiFit_Ds_PX_vec.clear();
    phiFit_Ds_PY_vec.clear();
    phiFit_Ds_PZ_vec.clear();
    phiFit_Ds_M_vec.clear();

    phiFit_Kp_PP_vec.clear();
    phiFit_Kp_PL_vec.clear();
    phiFit_Km_PP_vec.clear();
    phiFit_Km_PL_vec.clear();

    phiFit_phi_PP_vec.clear();
    phiFit_phi_PL_vec.clear();
    phiFit_pi_PP_vec.clear();
    phiFit_pi_PL_vec.clear();

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
    DsFit_CHI2_vec.clear();
    DsFit_NDOF_vec.clear();
    DsFit_CHI2NDOF_vec.clear();
    DsFit_ENDVX_X_vec.clear();
    DsFit_ENDVX_Y_vec.clear();
    DsFit_ENDVX_Z_vec.clear();
    DsFit_ENDVX_XERR_vec.clear();
    DsFit_ENDVX_YERR_vec.clear();
    DsFit_ENDVX_ZERR_vec.clear();

    DsFit_Kp_ETA_vec.clear();
    DsFit_Kp_PHI_vec.clear();
    DsFit_Kp_P_vec.clear();
    DsFit_Kp_PT_vec.clear();
    DsFit_Kp_PX_vec.clear();
    DsFit_Kp_PY_vec.clear();
    DsFit_Kp_PZ_vec.clear();

    DsFit_Km_ETA_vec.clear();
    DsFit_Km_PHI_vec.clear();
    DsFit_Km_P_vec.clear();
    DsFit_Km_PT_vec.clear();
    DsFit_Km_PX_vec.clear();
    DsFit_Km_PY_vec.clear();
    DsFit_Km_PZ_vec.clear();

    DsFit_pi_ETA_vec.clear();
    DsFit_pi_PHI_vec.clear();
    DsFit_pi_P_vec.clear();
    DsFit_pi_PT_vec.clear();
    DsFit_pi_PX_vec.clear();
    DsFit_pi_PY_vec.clear();
    DsFit_pi_PZ_vec.clear();

    DsFit_phi_ETA_vec.clear();
    DsFit_phi_PHI_vec.clear();
    DsFit_phi_P_vec.clear();
    DsFit_phi_PT_vec.clear();
    DsFit_phi_PX_vec.clear();
    DsFit_phi_PY_vec.clear();
    DsFit_phi_PZ_vec.clear();
    DsFit_phi_M_vec.clear();

    DsFit_Ds_ETA_vec.clear();
    DsFit_Ds_PHI_vec.clear();
    DsFit_Ds_P_vec.clear();
    DsFit_Ds_PT_vec.clear();
    DsFit_Ds_PX_vec.clear();
    DsFit_Ds_PY_vec.clear();
    DsFit_Ds_PZ_vec.clear();
    DsFit_Ds_M_vec.clear();

    DsFit_Kp_PP_vec.clear();
    DsFit_Kp_PL_vec.clear();
    DsFit_Km_PP_vec.clear();
    DsFit_Km_PL_vec.clear();

    DsFit_phi_PP_vec.clear();
    DsFit_phi_PL_vec.clear();
    DsFit_pi_PP_vec.clear();
    DsFit_pi_PL_vec.clear();

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

    DsFit_Mconstraint_Ds_M_vec.clear();

    best_Kp_ETA_vec.clear();
    best_Kp_PHI_vec.clear();
    best_Kp_ORIVX_X_vec.clear();
    best_Kp_ORIVX_Y_vec.clear();
    best_Kp_ORIVX_Z_vec.clear();
    best_Kp_P_vec.clear();
    best_Kp_PT_vec.clear();
    best_Kp_PX_vec.clear();
    best_Kp_PY_vec.clear();
    best_Kp_PZ_vec.clear();

    best_Km_ETA_vec.clear();
    best_Km_PHI_vec.clear();
    best_Km_ORIVX_X_vec.clear();
    best_Km_ORIVX_Y_vec.clear();
    best_Km_ORIVX_Z_vec.clear();
    best_Km_P_vec.clear();
    best_Km_PT_vec.clear();
    best_Km_PX_vec.clear();
    best_Km_PY_vec.clear();
    best_Km_PZ_vec.clear();

    best_pi_ETA_vec.clear();
    best_pi_PHI_vec.clear();
    best_pi_ORIVX_X_vec.clear();
    best_pi_ORIVX_Y_vec.clear();
    best_pi_ORIVX_Z_vec.clear();
    best_pi_P_vec.clear();
    best_pi_PT_vec.clear();
    best_pi_PX_vec.clear();
    best_pi_PY_vec.clear();
    best_pi_PZ_vec.clear();

    best_phi_ETA_vec.clear();
    best_phi_PHI_vec.clear();
    best_phi_P_vec.clear();
    best_phi_PT_vec.clear();
    best_phi_PX_vec.clear();
    best_phi_PY_vec.clear();
    best_phi_PZ_vec.clear();
    best_phi_M_vec.clear();

    best_Ds_ETA_vec.clear();
    best_Ds_PHI_vec.clear();
    best_Ds_P_vec.clear();
    best_Ds_PT_vec.clear();
    best_Ds_PX_vec.clear();
    best_Ds_PY_vec.clear();
    best_Ds_PZ_vec.clear();
    best_Ds_M_vec.clear();

    best_Kp_PP_vec.clear();
    best_Kp_PL_vec.clear();
    best_Km_PP_vec.clear();
    best_Km_PL_vec.clear();

    best_phi_PP_vec.clear();
    best_phi_PL_vec.clear();
    best_pi_PP_vec.clear();
    best_pi_PL_vec.clear();

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
    best_phiFit_CHI2_vec.clear();
    best_phiFit_NDOF_vec.clear();
    best_phiFit_CHI2NDOF_vec.clear();
    best_phiFit_ENDVX_X_vec.clear();
    best_phiFit_ENDVX_Y_vec.clear();
    best_phiFit_ENDVX_Z_vec.clear();
    best_phiFit_ENDVX_XERR_vec.clear();
    best_phiFit_ENDVX_YERR_vec.clear();
    best_phiFit_ENDVX_ZERR_vec.clear();

    best_phiFit_Kp_ETA_vec.clear();
    best_phiFit_Kp_PHI_vec.clear();
    best_phiFit_Kp_P_vec.clear();
    best_phiFit_Kp_PT_vec.clear();
    best_phiFit_Kp_PX_vec.clear();
    best_phiFit_Kp_PY_vec.clear();
    best_phiFit_Kp_PZ_vec.clear();

    best_phiFit_Km_ETA_vec.clear();
    best_phiFit_Km_PHI_vec.clear();
    best_phiFit_Km_P_vec.clear();
    best_phiFit_Km_PT_vec.clear();
    best_phiFit_Km_PX_vec.clear();
    best_phiFit_Km_PY_vec.clear();
    best_phiFit_Km_PZ_vec.clear();

    best_phiFit_pi_ETA_vec.clear();
    best_phiFit_pi_PHI_vec.clear();
    best_phiFit_pi_P_vec.clear();
    best_phiFit_pi_PT_vec.clear();
    best_phiFit_pi_PX_vec.clear();
    best_phiFit_pi_PY_vec.clear();
    best_phiFit_pi_PZ_vec.clear();

    best_phiFit_phi_ETA_vec.clear();
    best_phiFit_phi_PHI_vec.clear();
    best_phiFit_phi_P_vec.clear();
    best_phiFit_phi_PT_vec.clear();
    best_phiFit_phi_PX_vec.clear();
    best_phiFit_phi_PY_vec.clear();
    best_phiFit_phi_PZ_vec.clear();
    best_phiFit_phi_M_vec.clear();

    best_phiFit_Ds_ETA_vec.clear();
    best_phiFit_Ds_PHI_vec.clear();
    best_phiFit_Ds_P_vec.clear();
    best_phiFit_Ds_PT_vec.clear();
    best_phiFit_Ds_PX_vec.clear();
    best_phiFit_Ds_PY_vec.clear();
    best_phiFit_Ds_PZ_vec.clear();
    best_phiFit_Ds_M_vec.clear();

    best_phiFit_Kp_PP_vec.clear();
    best_phiFit_Kp_PL_vec.clear();
    best_phiFit_Km_PP_vec.clear();
    best_phiFit_Km_PL_vec.clear();

    best_phiFit_phi_PP_vec.clear();
    best_phiFit_phi_PL_vec.clear();
    best_phiFit_pi_PP_vec.clear();
    best_phiFit_pi_PL_vec.clear();

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
    best_DsFit_CHI2_vec.clear();
    best_DsFit_NDOF_vec.clear();
    best_DsFit_CHI2NDOF_vec.clear();
    best_DsFit_ENDVX_X_vec.clear();
    best_DsFit_ENDVX_Y_vec.clear();
    best_DsFit_ENDVX_Z_vec.clear();
    best_DsFit_ENDVX_XERR_vec.clear();
    best_DsFit_ENDVX_YERR_vec.clear();
    best_DsFit_ENDVX_ZERR_vec.clear();

    best_DsFit_Kp_ETA_vec.clear();
    best_DsFit_Kp_PHI_vec.clear();
    best_DsFit_Kp_P_vec.clear();
    best_DsFit_Kp_PT_vec.clear();
    best_DsFit_Kp_PX_vec.clear();
    best_DsFit_Kp_PY_vec.clear();
    best_DsFit_Kp_PZ_vec.clear();

    best_DsFit_Km_ETA_vec.clear();
    best_DsFit_Km_PHI_vec.clear();
    best_DsFit_Km_P_vec.clear();
    best_DsFit_Km_PT_vec.clear();
    best_DsFit_Km_PX_vec.clear();
    best_DsFit_Km_PY_vec.clear();
    best_DsFit_Km_PZ_vec.clear();

    best_DsFit_pi_ETA_vec.clear();
    best_DsFit_pi_PHI_vec.clear();
    best_DsFit_pi_P_vec.clear();
    best_DsFit_pi_PT_vec.clear();
    best_DsFit_pi_PX_vec.clear();
    best_DsFit_pi_PY_vec.clear();
    best_DsFit_pi_PZ_vec.clear();

    best_DsFit_phi_ETA_vec.clear();
    best_DsFit_phi_PHI_vec.clear();
    best_DsFit_phi_P_vec.clear();
    best_DsFit_phi_PT_vec.clear();
    best_DsFit_phi_PX_vec.clear();
    best_DsFit_phi_PY_vec.clear();
    best_DsFit_phi_PZ_vec.clear();
    best_DsFit_phi_M_vec.clear();

    best_DsFit_Ds_ETA_vec.clear();
    best_DsFit_Ds_PHI_vec.clear();
    best_DsFit_Ds_P_vec.clear();
    best_DsFit_Ds_PT_vec.clear();
    best_DsFit_Ds_PX_vec.clear();
    best_DsFit_Ds_PY_vec.clear();
    best_DsFit_Ds_PZ_vec.clear();
    best_DsFit_Ds_M_vec.clear();

    best_DsFit_Kp_PP_vec.clear();
    best_DsFit_Kp_PL_vec.clear();
    best_DsFit_Km_PP_vec.clear();
    best_DsFit_Km_PL_vec.clear();

    best_DsFit_phi_PP_vec.clear();
    best_DsFit_phi_PL_vec.clear();
    best_DsFit_pi_PP_vec.clear();
    best_DsFit_pi_PL_vec.clear();

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

    best_DsFit_Mconstraint_Ds_M_vec.clear();
}

void RecoBestTree::Best_Fill_Vector(int idxmax)
{
    best_Kp_ETA_vec.push_back(Kp_ETA_vec[idxmax]);
    best_Kp_PHI_vec.push_back(Kp_PHI_vec[idxmax]);
    best_Kp_ORIVX_X_vec.push_back(Kp_ORIVX_X_vec[idxmax]);
    best_Kp_ORIVX_Y_vec.push_back(Kp_ORIVX_Y_vec[idxmax]);
    best_Kp_ORIVX_Z_vec.push_back(Kp_ORIVX_Z_vec[idxmax]);
    best_Kp_P_vec.push_back(Kp_P_vec[idxmax]);
    best_Kp_PT_vec.push_back(Kp_PT_vec[idxmax]);
    best_Kp_PX_vec.push_back(Kp_PX_vec[idxmax]);
    best_Kp_PY_vec.push_back(Kp_PY_vec[idxmax]);
    best_Kp_PZ_vec.push_back(Kp_PZ_vec[idxmax]);

    best_Km_ETA_vec.push_back(Km_ETA_vec[idxmax]);
    best_Km_PHI_vec.push_back(Km_PHI_vec[idxmax]);
    best_Km_ORIVX_X_vec.push_back(Km_ORIVX_X_vec[idxmax]);
    best_Km_ORIVX_Y_vec.push_back(Km_ORIVX_Y_vec[idxmax]);
    best_Km_ORIVX_Z_vec.push_back(Km_ORIVX_Z_vec[idxmax]);
    best_Km_P_vec.push_back(Km_P_vec[idxmax]);
    best_Km_PT_vec.push_back(Km_PT_vec[idxmax]);
    best_Km_PX_vec.push_back(Km_PX_vec[idxmax]);
    best_Km_PY_vec.push_back(Km_PY_vec[idxmax]);
    best_Km_PZ_vec.push_back(Km_PZ_vec[idxmax]);

    best_pi_ETA_vec.push_back(pi_ETA_vec[idxmax]);
    best_pi_PHI_vec.push_back(pi_PHI_vec[idxmax]);
    best_pi_ORIVX_X_vec.push_back(pi_ORIVX_X_vec[idxmax]);
    best_pi_ORIVX_Y_vec.push_back(pi_ORIVX_Y_vec[idxmax]);
    best_pi_ORIVX_Z_vec.push_back(pi_ORIVX_Z_vec[idxmax]);
    best_pi_P_vec.push_back(pi_P_vec[idxmax]);
    best_pi_PT_vec.push_back(pi_PT_vec[idxmax]);
    best_pi_PX_vec.push_back(pi_PX_vec[idxmax]);
    best_pi_PY_vec.push_back(pi_PY_vec[idxmax]);
    best_pi_PZ_vec.push_back(pi_PZ_vec[idxmax]);

    best_phi_ETA_vec.push_back(phi_ETA_vec[idxmax]);
    best_phi_PHI_vec.push_back(phi_PHI_vec[idxmax]);
    best_phi_P_vec.push_back(phi_P_vec[idxmax]);
    best_phi_PT_vec.push_back(phi_PT_vec[idxmax]);
    best_phi_PX_vec.push_back(phi_PX_vec[idxmax]);
    best_phi_PY_vec.push_back(phi_PY_vec[idxmax]);
    best_phi_PZ_vec.push_back(phi_PZ_vec[idxmax]);
    best_phi_M_vec.push_back(phi_M_vec[idxmax]);

    best_Ds_ETA_vec.push_back(Ds_ETA_vec[idxmax]);
    best_Ds_PHI_vec.push_back(Ds_PHI_vec[idxmax]);
    best_Ds_P_vec.push_back(Ds_P_vec[idxmax]);
    best_Ds_PT_vec.push_back(Ds_PT_vec[idxmax]);
    best_Ds_PX_vec.push_back(Ds_PX_vec[idxmax]);
    best_Ds_PY_vec.push_back(Ds_PY_vec[idxmax]);
    best_Ds_PZ_vec.push_back(Ds_PZ_vec[idxmax]);
    best_Ds_M_vec.push_back(Ds_M_vec[idxmax]);

    best_Kp_PP_vec.push_back(Kp_PP_vec[idxmax]);
    best_Kp_PL_vec.push_back(Kp_PL_vec[idxmax]);
    best_Km_PP_vec.push_back(Km_PP_vec[idxmax]);
    best_Km_PL_vec.push_back(Km_PL_vec[idxmax]);

    best_phi_PP_vec.push_back(phi_PP_vec[idxmax]);
    best_phi_PL_vec.push_back(phi_PL_vec[idxmax]);
    best_pi_PP_vec.push_back(pi_PP_vec[idxmax]);
    best_pi_PL_vec.push_back(pi_PL_vec[idxmax]);

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
    best_phiFit_CHI2_vec.push_back(phiFit_CHI2_vec[idxmax]);
    best_phiFit_NDOF_vec.push_back(phiFit_NDOF_vec[idxmax]);
    best_phiFit_CHI2NDOF_vec.push_back(phiFit_CHI2NDOF_vec[idxmax]);
    best_phiFit_ENDVX_X_vec.push_back(phiFit_ENDVX_X_vec[idxmax]);
    best_phiFit_ENDVX_Y_vec.push_back(phiFit_ENDVX_Y_vec[idxmax]);
    best_phiFit_ENDVX_Z_vec.push_back(phiFit_ENDVX_Z_vec[idxmax]);
    best_phiFit_ENDVX_XERR_vec.push_back(phiFit_ENDVX_XERR_vec[idxmax]);
    best_phiFit_ENDVX_YERR_vec.push_back(phiFit_ENDVX_YERR_vec[idxmax]);
    best_phiFit_ENDVX_ZERR_vec.push_back(phiFit_ENDVX_ZERR_vec[idxmax]);

    best_phiFit_Kp_ETA_vec.push_back(phiFit_Kp_ETA_vec[idxmax]);
    best_phiFit_Kp_PHI_vec.push_back(phiFit_Kp_PHI_vec[idxmax]);
    best_phiFit_Kp_P_vec.push_back(phiFit_Kp_P_vec[idxmax]);
    best_phiFit_Kp_PT_vec.push_back(phiFit_Kp_PT_vec[idxmax]);
    best_phiFit_Kp_PX_vec.push_back(phiFit_Kp_PX_vec[idxmax]);
    best_phiFit_Kp_PY_vec.push_back(phiFit_Kp_PY_vec[idxmax]);
    best_phiFit_Kp_PZ_vec.push_back(phiFit_Kp_PZ_vec[idxmax]);

    best_phiFit_Km_ETA_vec.push_back(phiFit_Km_ETA_vec[idxmax]);
    best_phiFit_Km_PHI_vec.push_back(phiFit_Km_PHI_vec[idxmax]);
    best_phiFit_Km_P_vec.push_back(phiFit_Km_P_vec[idxmax]);
    best_phiFit_Km_PT_vec.push_back(phiFit_Km_PT_vec[idxmax]);
    best_phiFit_Km_PX_vec.push_back(phiFit_Km_PX_vec[idxmax]);
    best_phiFit_Km_PY_vec.push_back(phiFit_Km_PY_vec[idxmax]);
    best_phiFit_Km_PZ_vec.push_back(phiFit_Km_PZ_vec[idxmax]);

    best_phiFit_pi_ETA_vec.push_back(phiFit_pi_ETA_vec[idxmax]);
    best_phiFit_pi_PHI_vec.push_back(phiFit_pi_PHI_vec[idxmax]);
    best_phiFit_pi_P_vec.push_back(phiFit_pi_P_vec[idxmax]);
    best_phiFit_pi_PT_vec.push_back(phiFit_pi_PT_vec[idxmax]);
    best_phiFit_pi_PX_vec.push_back(phiFit_pi_PX_vec[idxmax]);
    best_phiFit_pi_PY_vec.push_back(phiFit_pi_PY_vec[idxmax]);
    best_phiFit_pi_PZ_vec.push_back(phiFit_pi_PZ_vec[idxmax]);

    best_phiFit_phi_ETA_vec.push_back(phiFit_phi_ETA_vec[idxmax]);
    best_phiFit_phi_PHI_vec.push_back(phiFit_phi_PHI_vec[idxmax]);
    best_phiFit_phi_P_vec.push_back(phiFit_phi_P_vec[idxmax]);
    best_phiFit_phi_PT_vec.push_back(phiFit_phi_PT_vec[idxmax]);
    best_phiFit_phi_PX_vec.push_back(phiFit_phi_PX_vec[idxmax]);
    best_phiFit_phi_PY_vec.push_back(phiFit_phi_PY_vec[idxmax]);
    best_phiFit_phi_PZ_vec.push_back(phiFit_phi_PZ_vec[idxmax]);
    best_phiFit_phi_M_vec.push_back(phiFit_phi_M_vec[idxmax]);

    best_phiFit_Ds_ETA_vec.push_back(phiFit_Ds_ETA_vec[idxmax]);
    best_phiFit_Ds_PHI_vec.push_back(phiFit_Ds_PHI_vec[idxmax]);
    best_phiFit_Ds_P_vec.push_back(phiFit_Ds_P_vec[idxmax]);
    best_phiFit_Ds_PT_vec.push_back(phiFit_Ds_PT_vec[idxmax]);
    best_phiFit_Ds_PX_vec.push_back(phiFit_Ds_PX_vec[idxmax]);
    best_phiFit_Ds_PY_vec.push_back(phiFit_Ds_PY_vec[idxmax]);
    best_phiFit_Ds_PZ_vec.push_back(phiFit_Ds_PZ_vec[idxmax]);
    best_phiFit_Ds_M_vec.push_back(phiFit_Ds_M_vec[idxmax]);

    best_phiFit_Kp_PP_vec.push_back(phiFit_Kp_PP_vec[idxmax]);
    best_phiFit_Kp_PL_vec.push_back(phiFit_Kp_PL_vec[idxmax]);
    best_phiFit_Km_PP_vec.push_back(phiFit_Km_PP_vec[idxmax]);
    best_phiFit_Km_PL_vec.push_back(phiFit_Km_PL_vec[idxmax]);

    best_phiFit_phi_PP_vec.push_back(phiFit_phi_PP_vec[idxmax]);
    best_phiFit_phi_PL_vec.push_back(phiFit_phi_PL_vec[idxmax]);
    best_phiFit_pi_PP_vec.push_back(phiFit_pi_PP_vec[idxmax]);
    best_phiFit_pi_PL_vec.push_back(phiFit_pi_PL_vec[idxmax]);

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
    best_DsFit_CHI2_vec.push_back(DsFit_CHI2_vec[idxmax]);
    best_DsFit_NDOF_vec.push_back(DsFit_NDOF_vec[idxmax]);
    best_DsFit_CHI2NDOF_vec.push_back(DsFit_CHI2NDOF_vec[idxmax]);
    best_DsFit_ENDVX_X_vec.push_back(DsFit_ENDVX_X_vec[idxmax]);
    best_DsFit_ENDVX_Y_vec.push_back(DsFit_ENDVX_Y_vec[idxmax]);
    best_DsFit_ENDVX_Z_vec.push_back(DsFit_ENDVX_Z_vec[idxmax]);
    best_DsFit_ENDVX_XERR_vec.push_back(DsFit_ENDVX_XERR_vec[idxmax]);
    best_DsFit_ENDVX_YERR_vec.push_back(DsFit_ENDVX_YERR_vec[idxmax]);
    best_DsFit_ENDVX_ZERR_vec.push_back(DsFit_ENDVX_ZERR_vec[idxmax]);

    best_DsFit_Kp_ETA_vec.push_back(DsFit_Kp_ETA_vec[idxmax]);
    best_DsFit_Kp_PHI_vec.push_back(DsFit_Kp_PHI_vec[idxmax]);
    best_DsFit_Kp_P_vec.push_back(DsFit_Kp_P_vec[idxmax]);
    best_DsFit_Kp_PT_vec.push_back(DsFit_Kp_PT_vec[idxmax]);
    best_DsFit_Kp_PX_vec.push_back(DsFit_Kp_PX_vec[idxmax]);
    best_DsFit_Kp_PY_vec.push_back(DsFit_Kp_PY_vec[idxmax]);
    best_DsFit_Kp_PZ_vec.push_back(DsFit_Kp_PZ_vec[idxmax]);

    best_DsFit_Km_ETA_vec.push_back(DsFit_Km_ETA_vec[idxmax]);
    best_DsFit_Km_PHI_vec.push_back(DsFit_Km_PHI_vec[idxmax]);
    best_DsFit_Km_P_vec.push_back(DsFit_Km_P_vec[idxmax]);
    best_DsFit_Km_PT_vec.push_back(DsFit_Km_PT_vec[idxmax]);
    best_DsFit_Km_PX_vec.push_back(DsFit_Km_PX_vec[idxmax]);
    best_DsFit_Km_PY_vec.push_back(DsFit_Km_PY_vec[idxmax]);
    best_DsFit_Km_PZ_vec.push_back(DsFit_Km_PZ_vec[idxmax]);

    best_DsFit_pi_ETA_vec.push_back(DsFit_pi_ETA_vec[idxmax]);
    best_DsFit_pi_PHI_vec.push_back(DsFit_pi_PHI_vec[idxmax]);
    best_DsFit_pi_P_vec.push_back(DsFit_pi_P_vec[idxmax]);
    best_DsFit_pi_PT_vec.push_back(DsFit_pi_PT_vec[idxmax]);
    best_DsFit_pi_PX_vec.push_back(DsFit_pi_PX_vec[idxmax]);
    best_DsFit_pi_PY_vec.push_back(DsFit_pi_PY_vec[idxmax]);
    best_DsFit_pi_PZ_vec.push_back(DsFit_pi_PZ_vec[idxmax]);

    best_DsFit_phi_ETA_vec.push_back(DsFit_phi_ETA_vec[idxmax]);
    best_DsFit_phi_PHI_vec.push_back(DsFit_phi_PHI_vec[idxmax]);
    best_DsFit_phi_P_vec.push_back(DsFit_phi_P_vec[idxmax]);
    best_DsFit_phi_PT_vec.push_back(DsFit_phi_PT_vec[idxmax]);
    best_DsFit_phi_PX_vec.push_back(DsFit_phi_PX_vec[idxmax]);
    best_DsFit_phi_PY_vec.push_back(DsFit_phi_PY_vec[idxmax]);
    best_DsFit_phi_PZ_vec.push_back(DsFit_phi_PZ_vec[idxmax]);
    best_DsFit_phi_M_vec.push_back(DsFit_phi_M_vec[idxmax]);

    best_DsFit_Ds_ETA_vec.push_back(DsFit_Ds_ETA_vec[idxmax]);
    best_DsFit_Ds_PHI_vec.push_back(DsFit_Ds_PHI_vec[idxmax]);
    best_DsFit_Ds_P_vec.push_back(DsFit_Ds_P_vec[idxmax]);
    best_DsFit_Ds_PT_vec.push_back(DsFit_Ds_PT_vec[idxmax]);
    best_DsFit_Ds_PX_vec.push_back(DsFit_Ds_PX_vec[idxmax]);
    best_DsFit_Ds_PY_vec.push_back(DsFit_Ds_PY_vec[idxmax]);
    best_DsFit_Ds_PZ_vec.push_back(DsFit_Ds_PZ_vec[idxmax]);
    best_DsFit_Ds_M_vec.push_back(DsFit_Ds_M_vec[idxmax]);

    best_DsFit_Kp_PP_vec.push_back(DsFit_Kp_PP_vec[idxmax]);
    best_DsFit_Kp_PL_vec.push_back(DsFit_Kp_PL_vec[idxmax]);
    best_DsFit_Km_PP_vec.push_back(DsFit_Km_PP_vec[idxmax]);
    best_DsFit_Km_PL_vec.push_back(DsFit_Km_PL_vec[idxmax]);

    best_DsFit_phi_PP_vec.push_back(DsFit_phi_PP_vec[idxmax]);
    best_DsFit_phi_PL_vec.push_back(DsFit_phi_PL_vec[idxmax]);
    best_DsFit_pi_PP_vec.push_back(DsFit_pi_PP_vec[idxmax]);
    best_DsFit_pi_PL_vec.push_back(DsFit_pi_PL_vec[idxmax]);

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

    best_DsFit_Mconstraint_Ds_M_vec.push_back(DsFit_Mconstraint_Ds_M_vec[idxmax]);
}

void RecoBestTree::Fill_Vector()
{
    Kp_ETA_vec.push_back(Kp_ETA);
    Kp_PHI_vec.push_back(Kp_PHI);
    Kp_ORIVX_X_vec.push_back(Kp_ORIVX_X);
    Kp_ORIVX_Y_vec.push_back(Kp_ORIVX_Y);
    Kp_ORIVX_Z_vec.push_back(Kp_ORIVX_Z);
    Kp_P_vec.push_back(Kp_P);
    Kp_PT_vec.push_back(Kp_PT);
    Kp_PX_vec.push_back(Kp_PX);
    Kp_PY_vec.push_back(Kp_PY);
    Kp_PZ_vec.push_back(Kp_PZ);

    Km_ETA_vec.push_back(Km_ETA);
    Km_PHI_vec.push_back(Km_PHI);
    Km_ORIVX_X_vec.push_back(Km_ORIVX_X);
    Km_ORIVX_Y_vec.push_back(Km_ORIVX_Y);
    Km_ORIVX_Z_vec.push_back(Km_ORIVX_Z);
    Km_P_vec.push_back(Km_P);
    Km_PT_vec.push_back(Km_PT);
    Km_PX_vec.push_back(Km_PX);
    Km_PY_vec.push_back(Km_PY);
    Km_PZ_vec.push_back(Km_PZ);

    pi_ETA_vec.push_back(pi_ETA);
    pi_PHI_vec.push_back(pi_PHI);
    pi_ORIVX_X_vec.push_back(pi_ORIVX_X);
    pi_ORIVX_Y_vec.push_back(pi_ORIVX_Y);
    pi_ORIVX_Z_vec.push_back(pi_ORIVX_Z);
    pi_P_vec.push_back(pi_P);
    pi_PT_vec.push_back(pi_PT);
    pi_PX_vec.push_back(pi_PX);
    pi_PY_vec.push_back(pi_PY);
    pi_PZ_vec.push_back(pi_PZ);

    phi_ETA_vec.push_back(phi_ETA);
    phi_PHI_vec.push_back(phi_PHI);
    phi_P_vec.push_back(phi_P);
    phi_PT_vec.push_back(phi_PT);
    phi_PX_vec.push_back(phi_PX);
    phi_PY_vec.push_back(phi_PY);
    phi_PZ_vec.push_back(phi_PZ);
    phi_M_vec.push_back(phi_M);

    Ds_ETA_vec.push_back(Ds_ETA);
    Ds_PHI_vec.push_back(Ds_PHI);
    Ds_P_vec.push_back(Ds_P);
    Ds_PT_vec.push_back(Ds_PT);
    Ds_PX_vec.push_back(Ds_PX);
    Ds_PY_vec.push_back(Ds_PY);
    Ds_PZ_vec.push_back(Ds_PZ);
    Ds_M_vec.push_back(Ds_M);

    Kp_PP_vec.push_back(Kp_PP);
    Kp_PL_vec.push_back(Kp_PL);
    Km_PP_vec.push_back(Km_PP);
    Km_PL_vec.push_back(Km_PL);

    phi_PP_vec.push_back(phi_PP);
    phi_PL_vec.push_back(phi_PL);
    pi_PP_vec.push_back(pi_PP);
    pi_PL_vec.push_back(pi_PL);

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
    phiFit_CHI2_vec.push_back(phiFit_CHI2);
    phiFit_NDOF_vec.push_back(phiFit_NDOF);
    phiFit_CHI2NDOF_vec.push_back(phiFit_CHI2NDOF);
    phiFit_ENDVX_X_vec.push_back(phiFit_ENDVX_X);
    phiFit_ENDVX_Y_vec.push_back(phiFit_ENDVX_Y);
    phiFit_ENDVX_Z_vec.push_back(phiFit_ENDVX_Z);
    phiFit_ENDVX_XERR_vec.push_back(phiFit_ENDVX_XERR);
    phiFit_ENDVX_YERR_vec.push_back(phiFit_ENDVX_YERR);
    phiFit_ENDVX_ZERR_vec.push_back(phiFit_ENDVX_ZERR);

    phiFit_Kp_ETA_vec.push_back(phiFit_Kp_ETA);
    phiFit_Kp_PHI_vec.push_back(phiFit_Kp_PHI);
    phiFit_Kp_P_vec.push_back(phiFit_Kp_P);
    phiFit_Kp_PT_vec.push_back(phiFit_Kp_PT);
    phiFit_Kp_PX_vec.push_back(phiFit_Kp_PX);
    phiFit_Kp_PY_vec.push_back(phiFit_Kp_PY);
    phiFit_Kp_PZ_vec.push_back(phiFit_Kp_PZ);

    phiFit_Km_ETA_vec.push_back(phiFit_Km_ETA);
    phiFit_Km_PHI_vec.push_back(phiFit_Km_PHI);
    phiFit_Km_P_vec.push_back(phiFit_Km_P);
    phiFit_Km_PT_vec.push_back(phiFit_Km_PT);
    phiFit_Km_PX_vec.push_back(phiFit_Km_PX);
    phiFit_Km_PY_vec.push_back(phiFit_Km_PY);
    phiFit_Km_PZ_vec.push_back(phiFit_Km_PZ);

    phiFit_pi_ETA_vec.push_back(phiFit_pi_ETA);
    phiFit_pi_PHI_vec.push_back(phiFit_pi_PHI);
    phiFit_pi_P_vec.push_back(phiFit_pi_P);
    phiFit_pi_PT_vec.push_back(phiFit_pi_PT);
    phiFit_pi_PX_vec.push_back(phiFit_pi_PX);
    phiFit_pi_PY_vec.push_back(phiFit_pi_PY);
    phiFit_pi_PZ_vec.push_back(phiFit_pi_PZ);

    phiFit_phi_ETA_vec.push_back(phiFit_phi_ETA);
    phiFit_phi_PHI_vec.push_back(phiFit_phi_PHI);
    phiFit_phi_P_vec.push_back(phiFit_phi_P);
    phiFit_phi_PT_vec.push_back(phiFit_phi_PT);
    phiFit_phi_PX_vec.push_back(phiFit_phi_PX);
    phiFit_phi_PY_vec.push_back(phiFit_phi_PY);
    phiFit_phi_PZ_vec.push_back(phiFit_phi_PZ);
    phiFit_phi_M_vec.push_back(phiFit_phi_M);

    phiFit_Ds_ETA_vec.push_back(phiFit_Ds_ETA);
    phiFit_Ds_PHI_vec.push_back(phiFit_Ds_PHI);
    phiFit_Ds_P_vec.push_back(phiFit_Ds_P);
    phiFit_Ds_PT_vec.push_back(phiFit_Ds_PT);
    phiFit_Ds_PX_vec.push_back(phiFit_Ds_PX);
    phiFit_Ds_PY_vec.push_back(phiFit_Ds_PY);
    phiFit_Ds_PZ_vec.push_back(phiFit_Ds_PZ);
    phiFit_Ds_M_vec.push_back(phiFit_Ds_M);

    phiFit_Kp_PP_vec.push_back(phiFit_Kp_PP);
    phiFit_Kp_PL_vec.push_back(phiFit_Kp_PL);
    phiFit_Km_PP_vec.push_back(phiFit_Km_PP);
    phiFit_Km_PL_vec.push_back(phiFit_Km_PL);

    phiFit_phi_PP_vec.push_back(phiFit_phi_PP);
    phiFit_phi_PL_vec.push_back(phiFit_phi_PL);
    phiFit_pi_PP_vec.push_back(phiFit_pi_PP);
    phiFit_pi_PL_vec.push_back(phiFit_pi_PL);

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
    DsFit_CHI2_vec.push_back(DsFit_CHI2);
    DsFit_NDOF_vec.push_back(DsFit_NDOF);
    DsFit_CHI2NDOF_vec.push_back(DsFit_CHI2NDOF);
    DsFit_ENDVX_X_vec.push_back(DsFit_ENDVX_X);
    DsFit_ENDVX_Y_vec.push_back(DsFit_ENDVX_Y);
    DsFit_ENDVX_Z_vec.push_back(DsFit_ENDVX_Z);
    DsFit_ENDVX_XERR_vec.push_back(DsFit_ENDVX_XERR);
    DsFit_ENDVX_YERR_vec.push_back(DsFit_ENDVX_YERR);
    DsFit_ENDVX_ZERR_vec.push_back(DsFit_ENDVX_ZERR);

    DsFit_Kp_ETA_vec.push_back(DsFit_Kp_ETA);
    DsFit_Kp_PHI_vec.push_back(DsFit_Kp_PHI);
    DsFit_Kp_P_vec.push_back(DsFit_Kp_P);
    DsFit_Kp_PT_vec.push_back(DsFit_Kp_PT);
    DsFit_Kp_PX_vec.push_back(DsFit_Kp_PX);
    DsFit_Kp_PY_vec.push_back(DsFit_Kp_PY);
    DsFit_Kp_PZ_vec.push_back(DsFit_Kp_PZ);

    DsFit_Km_ETA_vec.push_back(DsFit_Km_ETA);
    DsFit_Km_PHI_vec.push_back(DsFit_Km_PHI);
    DsFit_Km_P_vec.push_back(DsFit_Km_P);
    DsFit_Km_PT_vec.push_back(DsFit_Km_PT);
    DsFit_Km_PX_vec.push_back(DsFit_Km_PX);
    DsFit_Km_PY_vec.push_back(DsFit_Km_PY);
    DsFit_Km_PZ_vec.push_back(DsFit_Km_PZ);

    DsFit_pi_ETA_vec.push_back(DsFit_pi_ETA);
    DsFit_pi_PHI_vec.push_back(DsFit_pi_PHI);
    DsFit_pi_P_vec.push_back(DsFit_pi_P);
    DsFit_pi_PT_vec.push_back(DsFit_pi_PT);
    DsFit_pi_PX_vec.push_back(DsFit_pi_PX);
    DsFit_pi_PY_vec.push_back(DsFit_pi_PY);
    DsFit_pi_PZ_vec.push_back(DsFit_pi_PZ);

    DsFit_phi_ETA_vec.push_back(DsFit_phi_ETA);
    DsFit_phi_PHI_vec.push_back(DsFit_phi_PHI);
    DsFit_phi_P_vec.push_back(DsFit_phi_P);
    DsFit_phi_PT_vec.push_back(DsFit_phi_PT);
    DsFit_phi_PX_vec.push_back(DsFit_phi_PX);
    DsFit_phi_PY_vec.push_back(DsFit_phi_PY);
    DsFit_phi_PZ_vec.push_back(DsFit_phi_PZ);
    DsFit_phi_M_vec.push_back(DsFit_phi_M);

    DsFit_Ds_ETA_vec.push_back(DsFit_Ds_ETA);
    DsFit_Ds_PHI_vec.push_back(DsFit_Ds_PHI);
    DsFit_Ds_P_vec.push_back(DsFit_Ds_P);
    DsFit_Ds_PT_vec.push_back(DsFit_Ds_PT);
    DsFit_Ds_PX_vec.push_back(DsFit_Ds_PX);
    DsFit_Ds_PY_vec.push_back(DsFit_Ds_PY);
    DsFit_Ds_PZ_vec.push_back(DsFit_Ds_PZ);
    DsFit_Ds_M_vec.push_back(DsFit_Ds_M);

    DsFit_Kp_PP_vec.push_back(DsFit_Kp_PP);
    DsFit_Kp_PL_vec.push_back(DsFit_Kp_PL);
    DsFit_Km_PP_vec.push_back(DsFit_Km_PP);
    DsFit_Km_PL_vec.push_back(DsFit_Km_PL);

    DsFit_phi_PP_vec.push_back(DsFit_phi_PP);
    DsFit_phi_PL_vec.push_back(DsFit_phi_PL);
    DsFit_pi_PP_vec.push_back(DsFit_pi_PP);
    DsFit_pi_PL_vec.push_back(DsFit_pi_PL);

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

    DsFit_Mconstraint_Ds_M_vec.push_back(DsFit_Mconstraint_Ds_M);
}

void RecoBestTree::PV_withBS_Fill_Vector()
{
    PV_withBS_IsValid_vec.push_back(PV_withBS_IsValid);
    PV_withBS_IsFake_vec.push_back(PV_withBS_IsFake);
    PV_withBS_CHI2_vec.push_back(PV_withBS_CHI2);
    PV_withBS_NDOF_vec.push_back(PV_withBS_NDOF);
    PV_withBS_CHI2NDOF_vec.push_back(PV_withBS_CHI2NDOF);
    PV_withBS_X_vec.push_back(PV_withBS_X);
    PV_withBS_Y_vec.push_back(PV_withBS_Y);
    PV_withBS_Z_vec.push_back(PV_withBS_Z);
    PV_withBS_XERR_vec.push_back(PV_withBS_XERR);
    PV_withBS_YERR_vec.push_back(PV_withBS_YERR);
    PV_withBS_ZERR_vec.push_back(PV_withBS_ZERR);
}

void RecoBestTree::PV_noBS_Fill_Vector()
{
    PV_noBS_IsValid_vec.push_back(PV_noBS_IsValid);
    PV_noBS_IsFake_vec.push_back(PV_noBS_IsFake);
    PV_noBS_CHI2_vec.push_back(PV_noBS_CHI2);
    PV_noBS_NDOF_vec.push_back(PV_noBS_NDOF);
    PV_noBS_CHI2NDOF_vec.push_back(PV_noBS_CHI2NDOF);
    PV_noBS_X_vec.push_back(PV_noBS_X);
    PV_noBS_Y_vec.push_back(PV_noBS_Y);
    PV_noBS_Z_vec.push_back(PV_noBS_Z);
    PV_noBS_XERR_vec.push_back(PV_noBS_XERR);
    PV_noBS_YERR_vec.push_back(PV_noBS_YERR);
    PV_noBS_ZERR_vec.push_back(PV_noBS_ZERR);
}

void RecoBestTree::CreateBranches()
{
    tree->Branch("BS_type", &BS_type); 
    tree->Branch("BS_X0", &BS_X0);
    tree->Branch("BS_Y0", &BS_Y0);
    tree->Branch("BS_Z0", &BS_Z0);
    tree->Branch("BS_SigmaZ", &BS_SigmaZ);
    tree->Branch("BS_dXdZ", &BS_dXdZ);
    tree->Branch("BS_dYdZ", &BS_dYdZ);
    tree->Branch("BS_BWX", &BS_BWX);
    tree->Branch("BS_BWY", &BS_BWY);
    tree->Branch("BS_X0ERR", &BS_X0ERR);
    tree->Branch("BS_Y0ERR", &BS_Y0ERR);
    tree->Branch("BS_Z0ERR", &BS_Z0ERR);
    tree->Branch("BS_SigmaZ0ERR", &BS_SigmaZ0ERR);
    tree->Branch("BS_dXdZERR", &BS_dXdZERR);
    tree->Branch("BS_dYdZERR", &BS_dYdZERR);
    tree->Branch("BS_BWXERR", &BS_BWXERR);
    tree->Branch("BS_BWYERR", &BS_BWYERR);
    tree->Branch("BS_EmitX", &BS_EmitX);
    tree->Branch("BS_EmitY", &BS_EmitY);
    tree->Branch("BS_BetaStar", &BS_BetaStar);
    
    tree->Branch("PV_withBS_IsValid", &PV_withBS_IsValid_vec);
    tree->Branch("PV_withBS_IsFake", &PV_withBS_IsFake_vec);
    tree->Branch("PV_withBS_CHI2", &PV_withBS_CHI2_vec);
    tree->Branch("PV_withBS_NDOF", &PV_withBS_NDOF_vec);
    tree->Branch("PV_withBS_CHI2NDOF", &PV_withBS_CHI2NDOF_vec);
    tree->Branch("PV_withBS_X", &PV_withBS_X_vec);
    tree->Branch("PV_withBS_Y", &PV_withBS_Y_vec);
    tree->Branch("PV_withBS_Z", &PV_withBS_Z_vec);
    tree->Branch("PV_withBS_XERR", &PV_withBS_XERR_vec);
    tree->Branch("PV_withBS_YERR", &PV_withBS_YERR_vec);
    tree->Branch("PV_withBS_ZERR", &PV_withBS_ZERR_vec);
    tree->Branch("PV_noBS_IsValid", &PV_noBS_IsValid_vec);
    tree->Branch("PV_noBS_IsFake", &PV_noBS_IsFake_vec);
    tree->Branch("PV_noBS_CHI2", &PV_noBS_CHI2_vec);
    tree->Branch("PV_noBS_NDOF", &PV_noBS_NDOF_vec);
    tree->Branch("PV_noBS_CHI2NDOF", &PV_noBS_CHI2NDOF_vec);
    tree->Branch("PV_noBS_X", &PV_noBS_X_vec);
    tree->Branch("PV_noBS_Y", &PV_noBS_Y_vec);
    tree->Branch("PV_noBS_Z", &PV_noBS_Z_vec);
    tree->Branch("PV_noBS_XERR", &PV_noBS_XERR_vec);
    tree->Branch("PV_noBS_YERR", &PV_noBS_YERR_vec);
    tree->Branch("PV_noBS_ZERR", &PV_noBS_ZERR_vec);

    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

    tree->Branch("best_Kp_ETA", &best_Kp_ETA_vec);
    tree->Branch("best_Kp_PHI", &best_Kp_PHI_vec);
    tree->Branch("best_Kp_ORIVX_X", &best_Kp_ORIVX_X_vec);
    tree->Branch("best_Kp_ORIVX_Y", &best_Kp_ORIVX_Y_vec);
    tree->Branch("best_Kp_ORIVX_Z", &best_Kp_ORIVX_Z_vec);
    tree->Branch("best_Kp_P", &best_Kp_P_vec);
    tree->Branch("best_Kp_PT", &best_Kp_PT_vec);
    tree->Branch("best_Kp_PX", &best_Kp_PX_vec);
    tree->Branch("best_Kp_PY", &best_Kp_PY_vec);
    tree->Branch("best_Kp_PZ", &best_Kp_PZ_vec);

    tree->Branch("best_Km_ETA", &best_Km_ETA_vec);
    tree->Branch("best_Km_PHI", &best_Km_PHI_vec);
    tree->Branch("best_Km_ORIVX_X", &best_Km_ORIVX_X_vec);
    tree->Branch("best_Km_ORIVX_Y", &best_Km_ORIVX_Y_vec);
    tree->Branch("best_Km_ORIVX_Z", &best_Km_ORIVX_Z_vec);
    tree->Branch("best_Km_P", &best_Km_P_vec);
    tree->Branch("best_Km_PT", &best_Km_PT_vec);
    tree->Branch("best_Km_PX", &best_Km_PX_vec);
    tree->Branch("best_Km_PY", &best_Km_PY_vec);
    tree->Branch("best_Km_PZ", &best_Km_PZ_vec);

    tree->Branch("best_pi_ETA", &best_pi_ETA_vec);
    tree->Branch("best_pi_PHI", &best_pi_PHI_vec);
    tree->Branch("best_pi_ORIVX_X", &best_pi_ORIVX_X_vec);
    tree->Branch("best_pi_ORIVX_Y", &best_pi_ORIVX_Y_vec);
    tree->Branch("best_pi_ORIVX_Z", &best_pi_ORIVX_Z_vec);
    tree->Branch("best_pi_P", &best_pi_P_vec);
    tree->Branch("best_pi_PT", &best_pi_PT_vec);
    tree->Branch("best_pi_PX", &best_pi_PX_vec);
    tree->Branch("best_pi_PY", &best_pi_PY_vec);
    tree->Branch("best_pi_PZ", &best_pi_PZ_vec);

    tree->Branch("best_phi_ETA", &best_phi_ETA_vec);
    tree->Branch("best_phi_PHI", &best_phi_PHI_vec);
    tree->Branch("best_phi_P", &best_phi_P_vec);
    tree->Branch("best_phi_PT", &best_phi_PT_vec);
    tree->Branch("best_phi_PX", &best_phi_PX_vec);
    tree->Branch("best_phi_PY", &best_phi_PY_vec);
    tree->Branch("best_phi_PZ", &best_phi_PZ_vec);
    tree->Branch("best_phi_M", &best_phi_M_vec);

    tree->Branch("best_Ds_ETA", &best_Ds_ETA_vec);
    tree->Branch("best_Ds_PHI", &best_Ds_PHI_vec);
    tree->Branch("best_Ds_P", &best_Ds_P_vec);
    tree->Branch("best_Ds_PT", &best_Ds_PT_vec);
    tree->Branch("best_Ds_PX", &best_Ds_PX_vec);
    tree->Branch("best_Ds_PY", &best_Ds_PY_vec);
    tree->Branch("best_Ds_PZ", &best_Ds_PZ_vec);
    tree->Branch("best_Ds_M", &best_Ds_M_vec);

    tree->Branch("best_Kp_PP", &best_Kp_PP_vec);
    tree->Branch("best_Kp_PL", &best_Kp_PL_vec);
    tree->Branch("best_Km_PP", &best_Km_PP_vec);
    tree->Branch("best_Km_PL", &best_Km_PL_vec);

    tree->Branch("best_phi_PP", &best_phi_PP_vec);
    tree->Branch("best_phi_PL", &best_phi_PL_vec);
    tree->Branch("best_pi_PP", &best_pi_PP_vec);
    tree->Branch("best_pi_PL", &best_pi_PL_vec);

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

    tree->Branch("best_phiFit_CHI2", &best_phiFit_CHI2_vec);
    tree->Branch("best_phiFit_NDOF", &best_phiFit_NDOF_vec);
    tree->Branch("best_phiFit_CHI2NDOF", &best_phiFit_CHI2NDOF_vec);
    tree->Branch("best_phiFit_ENDVX_X", &best_phiFit_ENDVX_X_vec);
    tree->Branch("best_phiFit_ENDVX_Y", &best_phiFit_ENDVX_Y_vec);
    tree->Branch("best_phiFit_ENDVX_Z", &best_phiFit_ENDVX_Z_vec);
    tree->Branch("best_phiFit_ENDVX_XERR", &best_phiFit_ENDVX_XERR_vec);
    tree->Branch("best_phiFit_ENDVX_YERR", &best_phiFit_ENDVX_YERR_vec);
    tree->Branch("best_phiFit_ENDVX_ZERR", &best_phiFit_ENDVX_ZERR_vec);

    tree->Branch("best_phiFit_Kp_ETA", &best_phiFit_Kp_ETA_vec);
    tree->Branch("best_phiFit_Kp_PHI", &best_phiFit_Kp_PHI_vec);
    tree->Branch("best_phiFit_Kp_P", &best_phiFit_Kp_P_vec);
    tree->Branch("best_phiFit_Kp_PT", &best_phiFit_Kp_PT_vec);
    tree->Branch("best_phiFit_Kp_PX", &best_phiFit_Kp_PX_vec);
    tree->Branch("best_phiFit_Kp_PY", &best_phiFit_Kp_PY_vec);
    tree->Branch("best_phiFit_Kp_PZ", &best_phiFit_Kp_PZ_vec);

    tree->Branch("best_phiFit_Km_ETA", &best_phiFit_Km_ETA_vec);
    tree->Branch("best_phiFit_Km_PHI", &best_phiFit_Km_PHI_vec);
    tree->Branch("best_phiFit_Km_P", &best_phiFit_Km_P_vec);
    tree->Branch("best_phiFit_Km_PT", &best_phiFit_Km_PT_vec);
    tree->Branch("best_phiFit_Km_PX", &best_phiFit_Km_PX_vec);
    tree->Branch("best_phiFit_Km_PY", &best_phiFit_Km_PY_vec);
    tree->Branch("best_phiFit_Km_PZ", &best_phiFit_Km_PZ_vec);

    tree->Branch("best_phiFit_pi_ETA", &best_phiFit_pi_ETA_vec);
    tree->Branch("best_phiFit_pi_PHI", &best_phiFit_pi_PHI_vec);
    tree->Branch("best_phiFit_pi_P", &best_phiFit_pi_P_vec);
    tree->Branch("best_phiFit_pi_PT", &best_phiFit_pi_PT_vec);
    tree->Branch("best_phiFit_pi_PX", &best_phiFit_pi_PX_vec);
    tree->Branch("best_phiFit_pi_PY", &best_phiFit_pi_PY_vec);
    tree->Branch("best_phiFit_pi_PZ", &best_phiFit_pi_PZ_vec);

    tree->Branch("best_phiFit_phi_ETA", &best_phiFit_phi_ETA_vec);
    tree->Branch("best_phiFit_phi_PHI", &best_phiFit_phi_PHI_vec);
    tree->Branch("best_phiFit_phi_P", &best_phiFit_phi_P_vec);
    tree->Branch("best_phiFit_phi_PT", &best_phiFit_phi_PT_vec);
    tree->Branch("best_phiFit_phi_PX", &best_phiFit_phi_PX_vec);
    tree->Branch("best_phiFit_phi_PY", &best_phiFit_phi_PY_vec);
    tree->Branch("best_phiFit_phi_PZ", &best_phiFit_phi_PZ_vec);
    tree->Branch("best_phiFit_phi_M", &best_phiFit_phi_M_vec);

    tree->Branch("best_phiFit_Ds_ETA", &best_phiFit_Ds_ETA_vec);
    tree->Branch("best_phiFit_Ds_PHI", &best_phiFit_Ds_PHI_vec);
    tree->Branch("best_phiFit_Ds_P", &best_phiFit_Ds_P_vec);
    tree->Branch("best_phiFit_Ds_PT", &best_phiFit_Ds_PT_vec);
    tree->Branch("best_phiFit_Ds_PX", &best_phiFit_Ds_PX_vec);
    tree->Branch("best_phiFit_Ds_PY", &best_phiFit_Ds_PY_vec);
    tree->Branch("best_phiFit_Ds_PZ", &best_phiFit_Ds_PZ_vec);
    tree->Branch("best_phiFit_Ds_M", &best_phiFit_Ds_M_vec);

    tree->Branch("best_phiFit_Kp_PP", &best_phiFit_Kp_PP_vec);
    tree->Branch("best_phiFit_Kp_PL", &best_phiFit_Kp_PL_vec);
    tree->Branch("best_phiFit_Km_PP", &best_phiFit_Km_PP_vec);
    tree->Branch("best_phiFit_Km_PL", &best_phiFit_Km_PL_vec);

    tree->Branch("best_phiFit_phi_PP", &best_phiFit_phi_PP_vec);
    tree->Branch("best_phiFit_phi_PL", &best_phiFit_phi_PL_vec);
    tree->Branch("best_phiFit_pi_PP", &best_phiFit_pi_PP_vec);
    tree->Branch("best_phiFit_pi_PL", &best_phiFit_pi_PL_vec);

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
    tree->Branch("best_DsFit_CHI2", &best_DsFit_CHI2_vec);
    tree->Branch("best_DsFit_NDOF", &best_DsFit_NDOF_vec);
    tree->Branch("best_DsFit_CHI2NDOF", &best_DsFit_CHI2NDOF_vec);
    tree->Branch("best_DsFit_ENDVX_X", &best_DsFit_ENDVX_X_vec);
    tree->Branch("best_DsFit_ENDVX_Y", &best_DsFit_ENDVX_Y_vec);
    tree->Branch("best_DsFit_ENDVX_Z", &best_DsFit_ENDVX_Z_vec);
    tree->Branch("best_DsFit_ENDVX_XERR", &best_DsFit_ENDVX_XERR_vec);
    tree->Branch("best_DsFit_ENDVX_YERR", &best_DsFit_ENDVX_YERR_vec);
    tree->Branch("best_DsFit_ENDVX_ZERR", &best_DsFit_ENDVX_ZERR_vec);

    tree->Branch("best_DsFit_Kp_ETA", &best_DsFit_Kp_ETA_vec);
    tree->Branch("best_DsFit_Kp_PHI", &best_DsFit_Kp_PHI_vec);
    tree->Branch("best_DsFit_Kp_P", &best_DsFit_Kp_P_vec);
    tree->Branch("best_DsFit_Kp_PT", &best_DsFit_Kp_PT_vec);
    tree->Branch("best_DsFit_Kp_PX", &best_DsFit_Kp_PX_vec);
    tree->Branch("best_DsFit_Kp_PY", &best_DsFit_Kp_PY_vec);
    tree->Branch("best_DsFit_Kp_PZ", &best_DsFit_Kp_PZ_vec);

    tree->Branch("best_DsFit_Km_ETA", &best_DsFit_Km_ETA_vec);
    tree->Branch("best_DsFit_Km_PHI", &best_DsFit_Km_PHI_vec);
    tree->Branch("best_DsFit_Km_P", &best_DsFit_Km_P_vec);
    tree->Branch("best_DsFit_Km_PT", &best_DsFit_Km_PT_vec);
    tree->Branch("best_DsFit_Km_PX", &best_DsFit_Km_PX_vec);
    tree->Branch("best_DsFit_Km_PY", &best_DsFit_Km_PY_vec);
    tree->Branch("best_DsFit_Km_PZ", &best_DsFit_Km_PZ_vec);

    tree->Branch("best_DsFit_pi_ETA", &best_DsFit_pi_ETA_vec);
    tree->Branch("best_DsFit_pi_PHI", &best_DsFit_pi_PHI_vec);
    tree->Branch("best_DsFit_pi_P", &best_DsFit_pi_P_vec);
    tree->Branch("best_DsFit_pi_PT", &best_DsFit_pi_PT_vec);
    tree->Branch("best_DsFit_pi_PX", &best_DsFit_pi_PX_vec);
    tree->Branch("best_DsFit_pi_PY", &best_DsFit_pi_PY_vec);
    tree->Branch("best_DsFit_pi_PZ", &best_DsFit_pi_PZ_vec);

    tree->Branch("best_DsFit_phi_ETA", &best_DsFit_phi_ETA_vec);
    tree->Branch("best_DsFit_phi_PHI", &best_DsFit_phi_PHI_vec);
    tree->Branch("best_DsFit_phi_P", &best_DsFit_phi_P_vec);
    tree->Branch("best_DsFit_phi_PT", &best_DsFit_phi_PT_vec);
    tree->Branch("best_DsFit_phi_PX", &best_DsFit_phi_PX_vec);
    tree->Branch("best_DsFit_phi_PY", &best_DsFit_phi_PY_vec);
    tree->Branch("best_DsFit_phi_PZ", &best_DsFit_phi_PZ_vec);
    tree->Branch("best_DsFit_phi_M", &best_DsFit_phi_M_vec);

    tree->Branch("best_DsFit_Ds_ETA", &best_DsFit_Ds_ETA_vec);
    tree->Branch("best_DsFit_Ds_PHI", &best_DsFit_Ds_PHI_vec);
    tree->Branch("best_DsFit_Ds_P", &best_DsFit_Ds_P_vec);
    tree->Branch("best_DsFit_Ds_PT", &best_DsFit_Ds_PT_vec);
    tree->Branch("best_DsFit_Ds_PX", &best_DsFit_Ds_PX_vec);
    tree->Branch("best_DsFit_Ds_PY", &best_DsFit_Ds_PY_vec);
    tree->Branch("best_DsFit_Ds_PZ", &best_DsFit_Ds_PZ_vec);
    tree->Branch("best_DsFit_Ds_M", &best_DsFit_Ds_M_vec);

    tree->Branch("best_DsFit_Kp_PP", &best_DsFit_Kp_PP_vec);
    tree->Branch("best_DsFit_Kp_PL", &best_DsFit_Kp_PL_vec);
    tree->Branch("best_DsFit_Km_PP", &best_DsFit_Km_PP_vec);
    tree->Branch("best_DsFit_Km_PL", &best_DsFit_Km_PL_vec);

    tree->Branch("best_DsFit_phi_PP", &best_DsFit_phi_PP_vec);
    tree->Branch("best_DsFit_phi_PL", &best_DsFit_phi_PL_vec);
    tree->Branch("best_DsFit_pi_PP", &best_DsFit_pi_PP_vec);
    tree->Branch("best_DsFit_pi_PL", &best_DsFit_pi_PL_vec);

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

    tree->Branch("best_DsFit_Mconstraint_Ds_M", &best_DsFit_Mconstraint_Ds_M_vec);
}

void RecoBestTree::Kp_Reset()
{
    Kp_ETA = null;
    Kp_PHI = null;
    Kp_ORIVX_X = null;
    Kp_ORIVX_Y = null;
    Kp_ORIVX_Z = null;
    Kp_P = null;
    Kp_PT = null;
    Kp_PX = null;
    Kp_PY = null;
    Kp_PZ = null;
}

void RecoBestTree::Km_Reset()
{
    Km_ETA = null;
    Km_PHI = null;
    Km_ORIVX_X = null;
    Km_ORIVX_Y = null;
    Km_ORIVX_Z = null;
    Km_P = null;
    Km_PT = null;
    Km_PX = null;
    Km_PY = null;
    Km_PZ = null;

    phi_ETA = null;
    phi_PHI = null;
    phi_P = null;
    phi_PT = null;
    phi_PX = null;
    phi_PY = null;
    phi_PZ = null;
    phi_M = null;

    Kp_PP = null;
    Kp_PL = null;
    Km_PP = null;
    Km_PL = null;

    dR_Kp_Km = null;
    dR_Kp_phi = null;
    dR_Km_phi = null;

    phiFit_CHI2 = null;
    phiFit_NDOF = null;
    phiFit_CHI2NDOF = null;
    phiFit_ENDVX_X = null;
    phiFit_ENDVX_Y = null;
    phiFit_ENDVX_Z = null;
    phiFit_ENDVX_XERR = null;
    phiFit_ENDVX_YERR = null;
    phiFit_ENDVX_ZERR = null;

    phiFit_Kp_ETA = null;
    phiFit_Kp_PHI = null;
    phiFit_Kp_P = null;
    phiFit_Kp_PT = null;
    phiFit_Kp_PX = null;
    phiFit_Kp_PY = null;
    phiFit_Kp_PZ = null;

    phiFit_Km_ETA = null;
    phiFit_Km_PHI = null;
    phiFit_Km_P = null;
    phiFit_Km_PT = null;
    phiFit_Km_PX = null;
    phiFit_Km_PY = null;
    phiFit_Km_PZ = null;

    phiFit_phi_ETA = null;
    phiFit_phi_PHI = null;
    phiFit_phi_P = null;
    phiFit_phi_PT = null;
    phiFit_phi_PX = null;
    phiFit_phi_PY = null;
    phiFit_phi_PZ = null;
    phiFit_phi_M = null;

    phiFit_Kp_PP = null;
    phiFit_Kp_PL = null;
    phiFit_Km_PP = null;
    phiFit_Km_PL = null;

    phiFit_dR_Kp_Km = null;
    phiFit_dR_Kp_phi = null;
    phiFit_dR_Km_phi = null;

    dxy_Kp_Km = null;
    dxy_Kp_phi = null;
    dxy_Km_phi = null;

    dz_Kp_Km = null;
    dz_Kp_phi = null;
    dz_Km_phi = null;
}

void RecoBestTree::pi_Reset()
{
    pi_ETA = null;
    pi_PHI = null;
    pi_ORIVX_X = null;
    pi_ORIVX_Y = null;
    pi_ORIVX_Z = null;
    pi_P = null;
    pi_PT = null;
    pi_PX = null;
    pi_PY = null;
    pi_PZ = null;

    Ds_ETA = null;
    Ds_PHI = null;
    Ds_P = null;
    Ds_PT = null;
    Ds_PX = null;
    Ds_PY = null;
    Ds_PZ = null;
    Ds_M = null;

    phi_PP = null;
    phi_PL = null;
    pi_PP = null;
    pi_PL = null;

    dR_Kp_pi = null;
    dR_Km_pi = null;
    dR_pi_phi = null;
    dR_Kp_Ds = null;
    dR_Km_Ds = null;
    dR_pi_Ds = null;
    dR_phi_Ds = null;

    phiFit_pi_ETA = null;
    phiFit_pi_PHI = null;
    phiFit_pi_P = null;
    phiFit_pi_PT = null;
    phiFit_pi_PX = null;
    phiFit_pi_PY = null;
    phiFit_pi_PZ = null;

    phiFit_Ds_ETA = null;
    phiFit_Ds_PHI = null;
    phiFit_Ds_P = null;
    phiFit_Ds_PT = null;
    phiFit_Ds_PX = null;
    phiFit_Ds_PY = null;
    phiFit_Ds_PZ = null;
    phiFit_Ds_M = null;

    phiFit_phi_PP = null;
    phiFit_phi_PL = null;
    phiFit_pi_PP = null;
    phiFit_pi_PL = null;

    phiFit_dR_Kp_pi = null;
    phiFit_dR_Km_pi = null;
    phiFit_dR_pi_phi = null;
    phiFit_dR_Kp_Ds = null;
    phiFit_dR_Km_Ds = null;
    phiFit_dR_pi_Ds = null;
    phiFit_dR_phi_Ds = null;

    DsFit_CHI2 = null;
    DsFit_NDOF = null;
    DsFit_CHI2NDOF = null;
    DsFit_ENDVX_X = null;
    DsFit_ENDVX_Y = null;
    DsFit_ENDVX_Z = null;
    DsFit_ENDVX_XERR = null;
    DsFit_ENDVX_YERR = null;
    DsFit_ENDVX_ZERR = null;

    DsFit_Kp_ETA = null;
    DsFit_Kp_PHI = null;
    DsFit_Kp_P = null;
    DsFit_Kp_PT = null;
    DsFit_Kp_PX = null;
    DsFit_Kp_PY = null;
    DsFit_Kp_PZ = null;

    DsFit_Km_ETA = null;
    DsFit_Km_PHI = null;
    DsFit_Km_P = null;
    DsFit_Km_PT = null;
    DsFit_Km_PX = null;
    DsFit_Km_PY = null;
    DsFit_Km_PZ = null;

    DsFit_pi_ETA = null;
    DsFit_pi_PHI = null;
    DsFit_pi_P = null;
    DsFit_pi_PT = null;
    DsFit_pi_PX = null;
    DsFit_pi_PY = null;
    DsFit_pi_PZ = null;

    DsFit_phi_ETA = null;
    DsFit_phi_PHI = null;
    DsFit_phi_P = null;
    DsFit_phi_PT = null;
    DsFit_phi_PX = null;
    DsFit_phi_PY = null;
    DsFit_phi_PZ = null;
    DsFit_phi_M = null;

    DsFit_Ds_ETA = null;
    DsFit_Ds_PHI = null;
    DsFit_Ds_P = null;
    DsFit_Ds_PT = null;
    DsFit_Ds_PX = null;
    DsFit_Ds_PY = null;
    DsFit_Ds_PZ = null;
    DsFit_Ds_M = null;

    DsFit_Mconstraint_Ds_M = null;

    DsFit_Kp_PP = null;
    DsFit_Kp_PL = null;
    DsFit_Km_PP = null;
    DsFit_Km_PL = null;

    DsFit_phi_PP = null;
    DsFit_phi_PL = null;
    DsFit_pi_PP = null;
    DsFit_pi_PL = null;

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

    dxy_Kp_pi = null;
    dxy_Km_pi = null;
    dxy_pi_phi = null;
    dxy_Kp_Ds = null;
    dxy_Km_Ds = null;
    dxy_pi_Ds = null;
    dxy_phi_Ds = null;

    dz_Kp_pi = null;
    dz_Km_pi = null;
    dz_pi_phi = null;
    dz_Kp_Ds = null;
    dz_Km_Ds = null;
    dz_pi_Ds = null;
    dz_phi_Ds = null;
}

void RecoBestTree::BS_Reset()
{
    BS_type = null; 
    BS_X0 = null;
    BS_Y0 = null;
    BS_Z0 = null;
    BS_SigmaZ = null;
    BS_dXdZ = null;
    BS_dYdZ = null;
    BS_BWX = null;
    BS_BWY = null;
    BS_X0ERR = null;
    BS_Y0ERR = null;
    BS_Z0ERR = null;
    BS_SigmaZ0ERR = null;
    BS_dXdZERR = null;
    BS_dYdZERR = null;
    BS_BWXERR = null;
    BS_BWYERR = null;
    BS_EmitX = null;
    BS_EmitY = null;
    BS_BetaStar = null;
}

void RecoBestTree::PV_Reset()
{
    PV_withBS_IsValid = false;
    PV_withBS_IsFake = false;
    PV_withBS_CHI2 = null;
    PV_withBS_NDOF = null;
    PV_withBS_CHI2NDOF = null;
    PV_withBS_X = null;
    PV_withBS_Y = null;
    PV_withBS_Z = null;
    PV_withBS_XERR = null;
    PV_withBS_YERR = null;
    PV_withBS_ZERR = null;
    PV_noBS_IsValid = false;
    PV_noBS_IsFake = false;
    PV_noBS_CHI2 = null;
    PV_noBS_NDOF = null;
    PV_noBS_CHI2NDOF = null;
    PV_noBS_X = null;
    PV_noBS_Y = null;
    PV_noBS_Z = null;
    PV_noBS_XERR = null;
    PV_noBS_YERR = null;
    PV_noBS_ZERR = null;
}

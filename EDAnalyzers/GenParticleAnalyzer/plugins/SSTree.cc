#include "EDAnalyzers/GenParticleAnalyzer/interface/SSTree.h"

SSTree::SSTree(TTree *tree_)
{
    tree = tree_;
}
void SSTree::Init()
{
    // Gen particles
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    num_Gen_phi = 0;
    num_Gen_Ds = 0;

    // Matched particles
    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;

    match_Kp_idx = null;
    match_Km_idx = null;
    match_pi_idx = null;

    Gen_Kp_ETA_vec.clear();
    Gen_Kp_PHI_vec.clear();
    Gen_Kp_ORIVX_X_vec.clear();
    Gen_Kp_ORIVX_Y_vec.clear();
    Gen_Kp_ORIVX_Z_vec.clear();
    Gen_Kp_P_vec.clear();
    Gen_Kp_PT_vec.clear();
    Gen_Kp_PX_vec.clear();
    Gen_Kp_PY_vec.clear();
    Gen_Kp_PZ_vec.clear();
    Gen_Kp_PP_vec.clear();
    Gen_Kp_PL_vec.clear();

    Gen_Km_ETA_vec.clear();
    Gen_Km_PHI_vec.clear();
    Gen_Km_ORIVX_X_vec.clear();
    Gen_Km_ORIVX_Y_vec.clear();
    Gen_Km_ORIVX_Z_vec.clear();
    Gen_Km_P_vec.clear();
    Gen_Km_PT_vec.clear();
    Gen_Km_PX_vec.clear();
    Gen_Km_PY_vec.clear();
    Gen_Km_PZ_vec.clear();
    Gen_Km_PP_vec.clear();
    Gen_Km_PL_vec.clear();

    Gen_pi_ETA_vec.clear();
    Gen_pi_PHI_vec.clear();
    Gen_pi_ORIVX_X_vec.clear();
    Gen_pi_ORIVX_Y_vec.clear();
    Gen_pi_ORIVX_Z_vec.clear();
    Gen_pi_P_vec.clear();
    Gen_pi_PT_vec.clear();
    Gen_pi_PX_vec.clear();
    Gen_pi_PY_vec.clear();
    Gen_pi_PZ_vec.clear();
    Gen_pi_PP_vec.clear();
    Gen_pi_PL_vec.clear();

    Gen_phi_ETA_vec.clear();
    Gen_phi_PHI_vec.clear();
    Gen_phi_ORIVX_X_vec.clear();
    Gen_phi_ORIVX_Y_vec.clear();
    Gen_phi_ORIVX_Z_vec.clear();
    Gen_phi_P_vec.clear();
    Gen_phi_PT_vec.clear();
    Gen_phi_PX_vec.clear();
    Gen_phi_PY_vec.clear();
    Gen_phi_PZ_vec.clear();
    Gen_phi_PP_vec.clear();
    Gen_phi_PL_vec.clear();

    Gen_Ds_ETA_vec.clear();
    Gen_Ds_PHI_vec.clear();
    Gen_Ds_ORIVX_X_vec.clear();
    Gen_Ds_ORIVX_Y_vec.clear();
    Gen_Ds_ORIVX_Z_vec.clear();
    Gen_Ds_P_vec.clear();
    Gen_Ds_PT_vec.clear();
    Gen_Ds_PX_vec.clear();
    Gen_Ds_PY_vec.clear();
    Gen_Ds_PZ_vec.clear();

    Gen_dR_Kp_Km_vec.clear();
    Gen_dR_Kp_phi_vec.clear();
    Gen_dR_Km_phi_vec.clear();
    Gen_dR_Kp_pi_vec.clear();
    Gen_dR_Km_pi_vec.clear();
    Gen_dR_pi_phi_vec.clear();
    Gen_dR_Kp_Ds_vec.clear();
    Gen_dR_Km_Ds_vec.clear();
    Gen_dR_phi_Ds_vec.clear();
    Gen_dR_pi_Ds_vec.clear();

    Gen_dxy_Kp_Km_vec.clear();
    Gen_dxy_Kp_pi_vec.clear();
    Gen_dxy_Km_pi_vec.clear();
    Gen_dz_Kp_Km_vec.clear();
    Gen_dz_Kp_pi_vec.clear();
    Gen_dz_Km_pi_vec.clear();

    // Original info
    match_Kp_ETA_vec.clear();
    match_Kp_PHI_vec.clear();
    match_Kp_ORIVX_X_vec.clear();
    match_Kp_ORIVX_Y_vec.clear();
    match_Kp_ORIVX_Z_vec.clear();
    match_Kp_P_vec.clear();
    match_Kp_PT_vec.clear();
    match_Kp_PX_vec.clear();
    match_Kp_PY_vec.clear();
    match_Kp_PZ_vec.clear();

    match_Km_ETA_vec.clear();
    match_Km_PHI_vec.clear();
    match_Km_ORIVX_X_vec.clear();
    match_Km_ORIVX_Y_vec.clear();
    match_Km_ORIVX_Z_vec.clear();
    match_Km_P_vec.clear();
    match_Km_PT_vec.clear();
    match_Km_PX_vec.clear();
    match_Km_PY_vec.clear();
    match_Km_PZ_vec.clear();

    match_pi_ETA_vec.clear();
    match_pi_PHI_vec.clear();
    match_pi_ORIVX_X_vec.clear();
    match_pi_ORIVX_Y_vec.clear();
    match_pi_ORIVX_Z_vec.clear();
    match_pi_P_vec.clear();
    match_pi_PT_vec.clear();
    match_pi_PX_vec.clear();
    match_pi_PY_vec.clear();
    match_pi_PZ_vec.clear();

    match_phi_ETA_vec.clear();
    match_phi_PHI_vec.clear();
    match_phi_P_vec.clear();
    match_phi_PT_vec.clear();
    match_phi_PX_vec.clear();
    match_phi_PY_vec.clear();
    match_phi_PZ_vec.clear();
    match_phi_M_vec.clear();

    match_Ds_ETA_vec.clear();
    match_Ds_PHI_vec.clear();
    match_Ds_P_vec.clear();
    match_Ds_PT_vec.clear();
    match_Ds_PX_vec.clear();
    match_Ds_PY_vec.clear();
    match_Ds_PZ_vec.clear();
    match_Ds_M_vec.clear();

    match_Kp_PP_vec.clear();
    match_Kp_PL_vec.clear();
    match_Km_PP_vec.clear();
    match_Km_PL_vec.clear();

    match_phi_PP_vec.clear();
    match_phi_PL_vec.clear();
    match_pi_PP_vec.clear();
    match_pi_PL_vec.clear();

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

    match_dxy_Kp_Km_vec.clear();
    match_dxy_Kp_phi_vec.clear();
    match_dxy_Km_phi_vec.clear();
    match_dxy_Kp_pi_vec.clear();
    match_dxy_Km_pi_vec.clear();
    match_dxy_pi_phi_vec.clear();
    match_dxy_Kp_Ds_vec.clear();
    match_dxy_Km_Ds_vec.clear();
    match_dxy_phi_Ds_vec.clear();
    match_dxy_pi_Ds_vec.clear();

    match_dz_Kp_Km_vec.clear();
    match_dz_Kp_phi_vec.clear();
    match_dz_Km_phi_vec.clear();
    match_dz_Kp_pi_vec.clear();
    match_dz_Km_pi_vec.clear();
    match_dz_pi_phi_vec.clear();
    match_dz_Kp_Ds_vec.clear();
    match_dz_Km_Ds_vec.clear();
    match_dz_phi_Ds_vec.clear();
    match_dz_pi_Ds_vec.clear();

    // Fit on phi
    match_phiFit_CHI2_vec.clear();
    match_phiFit_NDOF_vec.clear();
    match_phiFit_CHI2NDOF_vec.clear();
    match_phiFit_ENDVX_X_vec.clear();
    match_phiFit_ENDVX_Y_vec.clear();
    match_phiFit_ENDVX_Z_vec.clear();

    match_phiFit_Kp_ETA_vec.clear();
    match_phiFit_Kp_PHI_vec.clear();
    match_phiFit_Kp_P_vec.clear();
    match_phiFit_Kp_PT_vec.clear();
    match_phiFit_Kp_PX_vec.clear();
    match_phiFit_Kp_PY_vec.clear();
    match_phiFit_Kp_PZ_vec.clear();

    match_phiFit_Km_ETA_vec.clear();
    match_phiFit_Km_PHI_vec.clear();
    match_phiFit_Km_P_vec.clear();
    match_phiFit_Km_PT_vec.clear();
    match_phiFit_Km_PX_vec.clear();
    match_phiFit_Km_PY_vec.clear();
    match_phiFit_Km_PZ_vec.clear();

    match_phiFit_pi_ETA_vec.clear();
    match_phiFit_pi_PHI_vec.clear();
    match_phiFit_pi_P_vec.clear();
    match_phiFit_pi_PT_vec.clear();
    match_phiFit_pi_PX_vec.clear();
    match_phiFit_pi_PY_vec.clear();
    match_phiFit_pi_PZ_vec.clear();

    match_phiFit_phi_ETA_vec.clear();
    match_phiFit_phi_PHI_vec.clear();
    match_phiFit_phi_P_vec.clear();
    match_phiFit_phi_PT_vec.clear();
    match_phiFit_phi_PX_vec.clear();
    match_phiFit_phi_PY_vec.clear();
    match_phiFit_phi_PZ_vec.clear();
    match_phiFit_phi_M_vec.clear();

    match_phiFit_Ds_ETA_vec.clear();
    match_phiFit_Ds_PHI_vec.clear();
    match_phiFit_Ds_P_vec.clear();
    match_phiFit_Ds_PT_vec.clear();
    match_phiFit_Ds_PX_vec.clear();
    match_phiFit_Ds_PY_vec.clear();
    match_phiFit_Ds_PZ_vec.clear();
    match_phiFit_Ds_M_vec.clear();

    match_phiFit_Kp_PP_vec.clear();
    match_phiFit_Kp_PL_vec.clear();
    match_phiFit_Km_PP_vec.clear();
    match_phiFit_Km_PL_vec.clear();

    match_phiFit_phi_PP_vec.clear();
    match_phiFit_phi_PL_vec.clear();
    match_phiFit_pi_PP_vec.clear();
    match_phiFit_pi_PL_vec.clear();

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
    match_DsFit_CHI2_vec.clear();
    match_DsFit_NDOF_vec.clear();
    match_DsFit_CHI2NDOF_vec.clear();
    match_DsFit_ENDVX_X_vec.clear();
    match_DsFit_ENDVX_Y_vec.clear();
    match_DsFit_ENDVX_Z_vec.clear();

    match_DsFit_Kp_ETA_vec.clear();
    match_DsFit_Kp_PHI_vec.clear();
    match_DsFit_Kp_P_vec.clear();
    match_DsFit_Kp_PT_vec.clear();
    match_DsFit_Kp_PX_vec.clear();
    match_DsFit_Kp_PY_vec.clear();
    match_DsFit_Kp_PZ_vec.clear();

    match_DsFit_Km_ETA_vec.clear();
    match_DsFit_Km_PHI_vec.clear();
    match_DsFit_Km_P_vec.clear();
    match_DsFit_Km_PT_vec.clear();
    match_DsFit_Km_PX_vec.clear();
    match_DsFit_Km_PY_vec.clear();
    match_DsFit_Km_PZ_vec.clear();

    match_DsFit_pi_ETA_vec.clear();
    match_DsFit_pi_PHI_vec.clear();
    match_DsFit_pi_P_vec.clear();
    match_DsFit_pi_PT_vec.clear();
    match_DsFit_pi_PX_vec.clear();
    match_DsFit_pi_PY_vec.clear();
    match_DsFit_pi_PZ_vec.clear();

    match_DsFit_phi_ETA_vec.clear();
    match_DsFit_phi_PHI_vec.clear();
    match_DsFit_phi_P_vec.clear();
    match_DsFit_phi_PT_vec.clear();
    match_DsFit_phi_PX_vec.clear();
    match_DsFit_phi_PY_vec.clear();
    match_DsFit_phi_PZ_vec.clear();
    match_DsFit_phi_M_vec.clear();

    match_DsFit_Ds_ETA_vec.clear();
    match_DsFit_Ds_PHI_vec.clear();
    match_DsFit_Ds_P_vec.clear();
    match_DsFit_Ds_PT_vec.clear();
    match_DsFit_Ds_PX_vec.clear();
    match_DsFit_Ds_PY_vec.clear();
    match_DsFit_Ds_PZ_vec.clear();
    match_DsFit_Ds_M_vec.clear();

    match_DsFit_Kp_PP_vec.clear();
    match_DsFit_Kp_PL_vec.clear();
    match_DsFit_Km_PP_vec.clear();
    match_DsFit_Km_PL_vec.clear();

    match_DsFit_phi_PP_vec.clear();
    match_DsFit_phi_PL_vec.clear();
    match_DsFit_pi_PP_vec.clear();
    match_DsFit_pi_PL_vec.clear();

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

    match_DsFit_Mconstraint_Ds_M_vec.clear();

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

    signal_entry_vec.clear();
}

void SSTree::Gen_Fill_Vector()
{
    Gen_Kp_ETA_vec.push_back(Gen_Kp_ETA);
    Gen_Kp_PHI_vec.push_back(Gen_Kp_PHI);
    Gen_Kp_ORIVX_X_vec.push_back(Gen_Kp_ORIVX_X);
    Gen_Kp_ORIVX_Y_vec.push_back(Gen_Kp_ORIVX_Y);
    Gen_Kp_ORIVX_Z_vec.push_back(Gen_Kp_ORIVX_Z);
    Gen_Kp_P_vec.push_back(Gen_Kp_P);
    Gen_Kp_PT_vec.push_back(Gen_Kp_PT);
    Gen_Kp_PX_vec.push_back(Gen_Kp_PX);
    Gen_Kp_PY_vec.push_back(Gen_Kp_PY);
    Gen_Kp_PZ_vec.push_back(Gen_Kp_PZ);
    Gen_Kp_PP_vec.push_back(Gen_Kp_PP);
    Gen_Kp_PL_vec.push_back(Gen_Kp_PL);

    Gen_Km_ETA_vec.push_back(Gen_Km_ETA);
    Gen_Km_PHI_vec.push_back(Gen_Km_PHI);
    Gen_Km_ORIVX_X_vec.push_back(Gen_Km_ORIVX_X);
    Gen_Km_ORIVX_Y_vec.push_back(Gen_Km_ORIVX_Y);
    Gen_Km_ORIVX_Z_vec.push_back(Gen_Km_ORIVX_Z);
    Gen_Km_P_vec.push_back(Gen_Km_P);
    Gen_Km_PT_vec.push_back(Gen_Km_PT);
    Gen_Km_PX_vec.push_back(Gen_Km_PX);
    Gen_Km_PY_vec.push_back(Gen_Km_PY);
    Gen_Km_PZ_vec.push_back(Gen_Km_PZ);
    Gen_Km_PP_vec.push_back(Gen_Km_PP);
    Gen_Km_PL_vec.push_back(Gen_Km_PL);

    Gen_pi_ETA_vec.push_back(Gen_pi_ETA);
    Gen_pi_PHI_vec.push_back(Gen_pi_PHI);
    Gen_pi_ORIVX_X_vec.push_back(Gen_pi_ORIVX_X);
    Gen_pi_ORIVX_Y_vec.push_back(Gen_pi_ORIVX_Y);
    Gen_pi_ORIVX_Z_vec.push_back(Gen_pi_ORIVX_Z);
    Gen_pi_P_vec.push_back(Gen_pi_P);
    Gen_pi_PT_vec.push_back(Gen_pi_PT);
    Gen_pi_PX_vec.push_back(Gen_pi_PX);
    Gen_pi_PY_vec.push_back(Gen_pi_PY);
    Gen_pi_PZ_vec.push_back(Gen_pi_PZ);
    Gen_pi_PP_vec.push_back(Gen_pi_PP);
    Gen_pi_PL_vec.push_back(Gen_pi_PL);

    Gen_phi_ETA_vec.push_back(Gen_phi_ETA);
    Gen_phi_PHI_vec.push_back(Gen_phi_PHI);
    Gen_phi_ORIVX_X_vec.push_back(Gen_phi_ORIVX_X);
    Gen_phi_ORIVX_Y_vec.push_back(Gen_phi_ORIVX_Y);
    Gen_phi_ORIVX_Z_vec.push_back(Gen_phi_ORIVX_Z);
    Gen_phi_P_vec.push_back(Gen_phi_P);
    Gen_phi_PT_vec.push_back(Gen_phi_PT);
    Gen_phi_PX_vec.push_back(Gen_phi_PX);
    Gen_phi_PY_vec.push_back(Gen_phi_PY);
    Gen_phi_PZ_vec.push_back(Gen_phi_PZ);
    Gen_phi_PP_vec.push_back(Gen_phi_PP);
    Gen_phi_PL_vec.push_back(Gen_phi_PL);

    Gen_Ds_ETA_vec.push_back(Gen_Ds_ETA);
    Gen_Ds_PHI_vec.push_back(Gen_Ds_PHI);
    Gen_Ds_ORIVX_X_vec.push_back(Gen_Ds_ORIVX_X);
    Gen_Ds_ORIVX_Y_vec.push_back(Gen_Ds_ORIVX_Y);
    Gen_Ds_ORIVX_Z_vec.push_back(Gen_Ds_ORIVX_Z);
    Gen_Ds_P_vec.push_back(Gen_Ds_P);
    Gen_Ds_PT_vec.push_back(Gen_Ds_PT);
    Gen_Ds_PX_vec.push_back(Gen_Ds_PX);
    Gen_Ds_PY_vec.push_back(Gen_Ds_PY);
    Gen_Ds_PZ_vec.push_back(Gen_Ds_PZ);

    Gen_dR_Kp_Km_vec.push_back(Gen_dR_Kp_Km);
    Gen_dR_Kp_phi_vec.push_back(Gen_dR_Kp_phi);
    Gen_dR_Km_phi_vec.push_back(Gen_dR_Km_phi);
    Gen_dR_Kp_pi_vec.push_back(Gen_dR_Kp_pi);
    Gen_dR_Km_pi_vec.push_back(Gen_dR_Km_pi);
    Gen_dR_pi_phi_vec.push_back(Gen_dR_pi_phi);
    Gen_dR_Kp_Ds_vec.push_back(Gen_dR_Kp_Ds);
    Gen_dR_Km_Ds_vec.push_back(Gen_dR_Km_Ds);
    Gen_dR_phi_Ds_vec.push_back(Gen_dR_phi_Ds);
    Gen_dR_pi_Ds_vec.push_back(Gen_dR_pi_Ds);

    Gen_dxy_Kp_Km_vec.push_back(Gen_dxy_Kp_Km);
    Gen_dxy_Kp_pi_vec.push_back(Gen_dxy_Kp_pi);
    Gen_dxy_Km_pi_vec.push_back(Gen_dxy_Km_pi);
    Gen_dz_Kp_Km_vec.push_back(Gen_dz_Kp_Km);
    Gen_dz_Kp_pi_vec.push_back(Gen_dz_Kp_pi);
    Gen_dz_Km_pi_vec.push_back(Gen_dz_Km_pi);

}

void SSTree::Match_Fill_Vector()
{
    match_Kp_ETA_vec.push_back(match_Kp_ETA);
    match_Kp_PHI_vec.push_back(match_Kp_PHI);
    match_Kp_ORIVX_X_vec.push_back(match_Kp_ORIVX_X);
    match_Kp_ORIVX_Y_vec.push_back(match_Kp_ORIVX_Y);
    match_Kp_ORIVX_Z_vec.push_back(match_Kp_ORIVX_Z);
    match_Kp_P_vec.push_back(match_Kp_P);
    match_Kp_PT_vec.push_back(match_Kp_PT);
    match_Kp_PX_vec.push_back(match_Kp_PX);
    match_Kp_PY_vec.push_back(match_Kp_PY);
    match_Kp_PZ_vec.push_back(match_Kp_PZ);

    match_Km_ETA_vec.push_back(match_Km_ETA);
    match_Km_PHI_vec.push_back(match_Km_PHI);
    match_Km_ORIVX_X_vec.push_back(match_Km_ORIVX_X);
    match_Km_ORIVX_Y_vec.push_back(match_Km_ORIVX_Y);
    match_Km_ORIVX_Z_vec.push_back(match_Km_ORIVX_Z);
    match_Km_P_vec.push_back(match_Km_P);
    match_Km_PT_vec.push_back(match_Km_PT);
    match_Km_PX_vec.push_back(match_Km_PX);
    match_Km_PY_vec.push_back(match_Km_PY);
    match_Km_PZ_vec.push_back(match_Km_PZ);

    match_pi_ETA_vec.push_back(match_pi_ETA);
    match_pi_PHI_vec.push_back(match_pi_PHI);
    match_pi_ORIVX_X_vec.push_back(match_pi_ORIVX_X);
    match_pi_ORIVX_Y_vec.push_back(match_pi_ORIVX_Y);
    match_pi_ORIVX_Z_vec.push_back(match_pi_ORIVX_Z);
    match_pi_P_vec.push_back(match_pi_P);
    match_pi_PT_vec.push_back(match_pi_PT);
    match_pi_PX_vec.push_back(match_pi_PX);
    match_pi_PY_vec.push_back(match_pi_PY);
    match_pi_PZ_vec.push_back(match_pi_PZ);

    match_phi_ETA_vec.push_back(match_phi_ETA);
    match_phi_PHI_vec.push_back(match_phi_PHI);
    match_phi_P_vec.push_back(match_phi_P);
    match_phi_PT_vec.push_back(match_phi_PT);
    match_phi_PX_vec.push_back(match_phi_PX);
    match_phi_PY_vec.push_back(match_phi_PY);
    match_phi_PZ_vec.push_back(match_phi_PZ);
    match_phi_M_vec.push_back(match_phi_M);

    match_Ds_ETA_vec.push_back(match_Ds_ETA);
    match_Ds_PHI_vec.push_back(match_Ds_PHI);
    match_Ds_P_vec.push_back(match_Ds_P);
    match_Ds_PT_vec.push_back(match_Ds_PT);
    match_Ds_PX_vec.push_back(match_Ds_PX);
    match_Ds_PY_vec.push_back(match_Ds_PY);
    match_Ds_PZ_vec.push_back(match_Ds_PZ);
    match_Ds_M_vec.push_back(match_Ds_M);

    match_Kp_PP_vec.push_back(match_Kp_PP);
    match_Kp_PL_vec.push_back(match_Kp_PL);
    match_Km_PP_vec.push_back(match_Km_PP);
    match_Km_PL_vec.push_back(match_Km_PL);

    match_phi_PP_vec.push_back(match_phi_PP);
    match_phi_PL_vec.push_back(match_phi_PL);
    match_pi_PP_vec.push_back(match_pi_PP);
    match_pi_PL_vec.push_back(match_pi_PL);

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

    match_dxy_Kp_Km_vec.push_back(match_dxy_Kp_Km);
    match_dxy_Kp_phi_vec.push_back(match_dxy_Kp_phi);
    match_dxy_Km_phi_vec.push_back(match_dxy_Km_phi);
    match_dxy_Kp_pi_vec.push_back(match_dxy_Kp_pi);
    match_dxy_Km_pi_vec.push_back(match_dxy_Km_pi);
    match_dxy_pi_phi_vec.push_back(match_dxy_pi_phi);
    match_dxy_Kp_Ds_vec.push_back(match_dxy_Kp_Ds);
    match_dxy_Km_Ds_vec.push_back(match_dxy_Km_Ds);
    match_dxy_phi_Ds_vec.push_back(match_dxy_phi_Ds);
    match_dxy_pi_Ds_vec.push_back(match_dxy_pi_Ds);

    match_dz_Kp_Km_vec.push_back(match_dz_Kp_Km);
    match_dz_Kp_phi_vec.push_back(match_dz_Kp_phi);
    match_dz_Km_phi_vec.push_back(match_dz_Km_phi);
    match_dz_Kp_pi_vec.push_back(match_dz_Kp_pi);
    match_dz_Km_pi_vec.push_back(match_dz_Km_pi);
    match_dz_pi_phi_vec.push_back(match_dz_pi_phi);
    match_dz_Kp_Ds_vec.push_back(match_dz_Kp_Ds);
    match_dz_Km_Ds_vec.push_back(match_dz_Km_Ds);
    match_dz_phi_Ds_vec.push_back(match_dz_phi_Ds);
    match_dz_pi_Ds_vec.push_back(match_dz_pi_Ds);


    // Fit on phi
    match_phiFit_CHI2_vec.push_back(match_phiFit_CHI2);
    match_phiFit_NDOF_vec.push_back(match_phiFit_NDOF);
    match_phiFit_CHI2NDOF_vec.push_back(match_phiFit_CHI2NDOF);
    match_phiFit_ENDVX_X_vec.push_back(match_phiFit_ENDVX_X);
    match_phiFit_ENDVX_Y_vec.push_back(match_phiFit_ENDVX_Y);
    match_phiFit_ENDVX_Z_vec.push_back(match_phiFit_ENDVX_Z);

    match_phiFit_Kp_ETA_vec.push_back(match_phiFit_Kp_ETA);
    match_phiFit_Kp_PHI_vec.push_back(match_phiFit_Kp_PHI);
    match_phiFit_Kp_P_vec.push_back(match_phiFit_Kp_P);
    match_phiFit_Kp_PT_vec.push_back(match_phiFit_Kp_PT);
    match_phiFit_Kp_PX_vec.push_back(match_phiFit_Kp_PX);
    match_phiFit_Kp_PY_vec.push_back(match_phiFit_Kp_PY);
    match_phiFit_Kp_PZ_vec.push_back(match_phiFit_Kp_PZ);

    match_phiFit_Km_ETA_vec.push_back(match_phiFit_Km_ETA);
    match_phiFit_Km_PHI_vec.push_back(match_phiFit_Km_PHI);
    match_phiFit_Km_P_vec.push_back(match_phiFit_Km_P);
    match_phiFit_Km_PT_vec.push_back(match_phiFit_Km_PT);
    match_phiFit_Km_PX_vec.push_back(match_phiFit_Km_PX);
    match_phiFit_Km_PY_vec.push_back(match_phiFit_Km_PY);
    match_phiFit_Km_PZ_vec.push_back(match_phiFit_Km_PZ);

    match_phiFit_pi_ETA_vec.push_back(match_phiFit_pi_ETA);
    match_phiFit_pi_PHI_vec.push_back(match_phiFit_pi_PHI);
    match_phiFit_pi_P_vec.push_back(match_phiFit_pi_P);
    match_phiFit_pi_PT_vec.push_back(match_phiFit_pi_PT);
    match_phiFit_pi_PX_vec.push_back(match_phiFit_pi_PX);
    match_phiFit_pi_PY_vec.push_back(match_phiFit_pi_PY);
    match_phiFit_pi_PZ_vec.push_back(match_phiFit_pi_PZ);

    match_phiFit_phi_ETA_vec.push_back(match_phiFit_phi_ETA);
    match_phiFit_phi_PHI_vec.push_back(match_phiFit_phi_PHI);
    match_phiFit_phi_P_vec.push_back(match_phiFit_phi_P);
    match_phiFit_phi_PT_vec.push_back(match_phiFit_phi_PT);
    match_phiFit_phi_PX_vec.push_back(match_phiFit_phi_PX);
    match_phiFit_phi_PY_vec.push_back(match_phiFit_phi_PY);
    match_phiFit_phi_PZ_vec.push_back(match_phiFit_phi_PZ);
    match_phiFit_phi_M_vec.push_back(match_phiFit_phi_M);

    match_phiFit_Ds_ETA_vec.push_back(match_phiFit_Ds_ETA);
    match_phiFit_Ds_PHI_vec.push_back(match_phiFit_Ds_PHI);
    match_phiFit_Ds_P_vec.push_back(match_phiFit_Ds_P);
    match_phiFit_Ds_PT_vec.push_back(match_phiFit_Ds_PT);
    match_phiFit_Ds_PX_vec.push_back(match_phiFit_Ds_PX);
    match_phiFit_Ds_PY_vec.push_back(match_phiFit_Ds_PY);
    match_phiFit_Ds_PZ_vec.push_back(match_phiFit_Ds_PZ);
    match_phiFit_Ds_M_vec.push_back(match_phiFit_Ds_M);

    match_phiFit_Kp_PP_vec.push_back(match_phiFit_Kp_PP);
    match_phiFit_Kp_PL_vec.push_back(match_phiFit_Kp_PL);
    match_phiFit_Km_PP_vec.push_back(match_phiFit_Km_PP);
    match_phiFit_Km_PL_vec.push_back(match_phiFit_Km_PL);

    match_phiFit_phi_PP_vec.push_back(match_phiFit_phi_PP);
    match_phiFit_phi_PL_vec.push_back(match_phiFit_phi_PL);
    match_phiFit_pi_PP_vec.push_back(match_phiFit_pi_PP);
    match_phiFit_pi_PL_vec.push_back(match_phiFit_pi_PL);

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
    match_DsFit_CHI2_vec.push_back(match_DsFit_CHI2);
    match_DsFit_NDOF_vec.push_back(match_DsFit_NDOF);
    match_DsFit_CHI2NDOF_vec.push_back(match_DsFit_CHI2NDOF);
    match_DsFit_ENDVX_X_vec.push_back(match_DsFit_ENDVX_X);
    match_DsFit_ENDVX_Y_vec.push_back(match_DsFit_ENDVX_Y);
    match_DsFit_ENDVX_Z_vec.push_back(match_DsFit_ENDVX_Z);

    match_DsFit_Kp_ETA_vec.push_back(match_DsFit_Kp_ETA);
    match_DsFit_Kp_PHI_vec.push_back(match_DsFit_Kp_PHI);
    match_DsFit_Kp_P_vec.push_back(match_DsFit_Kp_P);
    match_DsFit_Kp_PT_vec.push_back(match_DsFit_Kp_PT);
    match_DsFit_Kp_PX_vec.push_back(match_DsFit_Kp_PX);
    match_DsFit_Kp_PY_vec.push_back(match_DsFit_Kp_PY);
    match_DsFit_Kp_PZ_vec.push_back(match_DsFit_Kp_PZ);

    match_DsFit_Km_ETA_vec.push_back(match_DsFit_Km_ETA);
    match_DsFit_Km_PHI_vec.push_back(match_DsFit_Km_PHI);
    match_DsFit_Km_P_vec.push_back(match_DsFit_Km_P);
    match_DsFit_Km_PT_vec.push_back(match_DsFit_Km_PT);
    match_DsFit_Km_PX_vec.push_back(match_DsFit_Km_PX);
    match_DsFit_Km_PY_vec.push_back(match_DsFit_Km_PY);
    match_DsFit_Km_PZ_vec.push_back(match_DsFit_Km_PZ);

    match_DsFit_pi_ETA_vec.push_back(match_DsFit_pi_ETA);
    match_DsFit_pi_PHI_vec.push_back(match_DsFit_pi_PHI);
    match_DsFit_pi_P_vec.push_back(match_DsFit_pi_P);
    match_DsFit_pi_PT_vec.push_back(match_DsFit_pi_PT);
    match_DsFit_pi_PX_vec.push_back(match_DsFit_pi_PX);
    match_DsFit_pi_PY_vec.push_back(match_DsFit_pi_PY);
    match_DsFit_pi_PZ_vec.push_back(match_DsFit_pi_PZ);

    match_DsFit_phi_ETA_vec.push_back(match_DsFit_phi_ETA);
    match_DsFit_phi_PHI_vec.push_back(match_DsFit_phi_PHI);
    match_DsFit_phi_P_vec.push_back(match_DsFit_phi_P);
    match_DsFit_phi_PT_vec.push_back(match_DsFit_phi_PT);
    match_DsFit_phi_PX_vec.push_back(match_DsFit_phi_PX);
    match_DsFit_phi_PY_vec.push_back(match_DsFit_phi_PY);
    match_DsFit_phi_PZ_vec.push_back(match_DsFit_phi_PZ);
    match_DsFit_phi_M_vec.push_back(match_DsFit_phi_M);

    match_DsFit_Ds_ETA_vec.push_back(match_DsFit_Ds_ETA);
    match_DsFit_Ds_PHI_vec.push_back(match_DsFit_Ds_PHI);
    match_DsFit_Ds_P_vec.push_back(match_DsFit_Ds_P);
    match_DsFit_Ds_PT_vec.push_back(match_DsFit_Ds_PT);
    match_DsFit_Ds_PX_vec.push_back(match_DsFit_Ds_PX);
    match_DsFit_Ds_PY_vec.push_back(match_DsFit_Ds_PY);
    match_DsFit_Ds_PZ_vec.push_back(match_DsFit_Ds_PZ);
    match_DsFit_Ds_M_vec.push_back(match_DsFit_Ds_M);

    match_DsFit_Kp_PP_vec.push_back(match_DsFit_Kp_PP);
    match_DsFit_Kp_PL_vec.push_back(match_DsFit_Kp_PL);
    match_DsFit_Km_PP_vec.push_back(match_DsFit_Km_PP);
    match_DsFit_Km_PL_vec.push_back(match_DsFit_Km_PL);

    match_DsFit_phi_PP_vec.push_back(match_DsFit_phi_PP);
    match_DsFit_phi_PL_vec.push_back(match_DsFit_phi_PL);
    match_DsFit_pi_PP_vec.push_back(match_DsFit_pi_PP);
    match_DsFit_pi_PL_vec.push_back(match_DsFit_pi_PL);

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

    match_DsFit_Mconstraint_Ds_M_vec.push_back(match_DsFit_Mconstraint_Ds_M);

}


void SSTree::Fill_Vector()
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

    signal_entry_vec.push_back(signal_entry);
}

void SSTree::CreateBranches()
{
    tree->Branch("num_Gen_Kp", &num_Gen_Kp);
    tree->Branch("num_Gen_Km", &num_Gen_Km);
    tree->Branch("num_Gen_pi", &num_Gen_pi);
    tree->Branch("num_Gen_phi", &num_Gen_phi);
    tree->Branch("num_Gen_Ds", &num_Gen_Ds);

    tree->Branch("Gen_Kp_ETA", &Gen_Kp_ETA_vec);
    tree->Branch("Gen_Kp_PHI", &Gen_Kp_PHI_vec);
    tree->Branch("Gen_Kp_ORIVX_X", &Gen_Kp_ORIVX_X_vec);
    tree->Branch("Gen_Kp_ORIVX_Y", &Gen_Kp_ORIVX_Y_vec);
    tree->Branch("Gen_Kp_ORIVX_Z", &Gen_Kp_ORIVX_Z_vec);
    tree->Branch("Gen_Kp_P", &Gen_Kp_P_vec);
    tree->Branch("Gen_Kp_PT", &Gen_Kp_PT_vec);
    tree->Branch("Gen_Kp_PX", &Gen_Kp_PX_vec);
    tree->Branch("Gen_Kp_PY", &Gen_Kp_PY_vec);
    tree->Branch("Gen_Kp_PZ", &Gen_Kp_PZ_vec);
    tree->Branch("Gen_Kp_PP", &Gen_Kp_PP_vec);
    tree->Branch("Gen_Kp_PL", &Gen_Kp_PL_vec);

    tree->Branch("Gen_Km_ETA", &Gen_Km_ETA_vec);
    tree->Branch("Gen_Km_PHI", &Gen_Km_PHI_vec);
    tree->Branch("Gen_Km_ORIVX_X", &Gen_Km_ORIVX_X_vec);
    tree->Branch("Gen_Km_ORIVX_Y", &Gen_Km_ORIVX_Y_vec);
    tree->Branch("Gen_Km_ORIVX_Z", &Gen_Km_ORIVX_Z_vec);
    tree->Branch("Gen_Km_P", &Gen_Km_P_vec);
    tree->Branch("Gen_Km_PT", &Gen_Km_PT_vec);
    tree->Branch("Gen_Km_PX", &Gen_Km_PX_vec);
    tree->Branch("Gen_Km_PY", &Gen_Km_PY_vec);
    tree->Branch("Gen_Km_PZ", &Gen_Km_PZ_vec);
    tree->Branch("Gen_Km_PP", &Gen_Km_PP_vec);
    tree->Branch("Gen_Km_PL", &Gen_Km_PL_vec);

    tree->Branch("Gen_pi_ETA", &Gen_pi_ETA_vec);
    tree->Branch("Gen_pi_PHI", &Gen_pi_PHI_vec);
    tree->Branch("Gen_pi_ORIVX_X", &Gen_pi_ORIVX_X_vec);
    tree->Branch("Gen_pi_ORIVX_Y", &Gen_pi_ORIVX_Y_vec);
    tree->Branch("Gen_pi_ORIVX_Z", &Gen_pi_ORIVX_Z_vec);
    tree->Branch("Gen_pi_P", &Gen_pi_P_vec);
    tree->Branch("Gen_pi_PT", &Gen_pi_PT_vec);
    tree->Branch("Gen_pi_PX", &Gen_pi_PX_vec);
    tree->Branch("Gen_pi_PY", &Gen_pi_PY_vec);
    tree->Branch("Gen_pi_PZ", &Gen_pi_PZ_vec);
    tree->Branch("Gen_pi_PP", &Gen_pi_PP_vec);
    tree->Branch("Gen_pi_PL", &Gen_pi_PL_vec);

    tree->Branch("Gen_phi_ETA", &Gen_phi_ETA_vec);
    tree->Branch("Gen_phi_PHI", &Gen_phi_PHI_vec);
    tree->Branch("Gen_phi_ORIVX_X", &Gen_phi_ORIVX_X_vec);
    tree->Branch("Gen_phi_ORIVX_Y", &Gen_phi_ORIVX_Y_vec);
    tree->Branch("Gen_phi_ORIVX_Z", &Gen_phi_ORIVX_Z_vec);
    tree->Branch("Gen_phi_P", &Gen_phi_P_vec);
    tree->Branch("Gen_phi_PT", &Gen_phi_PT_vec);
    tree->Branch("Gen_phi_PX", &Gen_phi_PX_vec);
    tree->Branch("Gen_phi_PY", &Gen_phi_PY_vec);
    tree->Branch("Gen_phi_PZ", &Gen_phi_PZ_vec);
    tree->Branch("Gen_phi_PP", &Gen_phi_PP_vec);
    tree->Branch("Gen_phi_PL", &Gen_phi_PL_vec);

    tree->Branch("Gen_Ds_ETA", &Gen_Ds_ETA_vec);
    tree->Branch("Gen_Ds_PHI", &Gen_Ds_PHI_vec);
    tree->Branch("Gen_Ds_ORIVX_X", &Gen_Ds_ORIVX_X_vec);
    tree->Branch("Gen_Ds_ORIVX_Y", &Gen_Ds_ORIVX_Y_vec);
    tree->Branch("Gen_Ds_ORIVX_Z", &Gen_Ds_ORIVX_Z_vec);
    tree->Branch("Gen_Ds_P", &Gen_Ds_P_vec);
    tree->Branch("Gen_Ds_PT", &Gen_Ds_PT_vec);
    tree->Branch("Gen_Ds_PX", &Gen_Ds_PX_vec);
    tree->Branch("Gen_Ds_PY", &Gen_Ds_PY_vec);
    tree->Branch("Gen_Ds_PZ", &Gen_Ds_PZ_vec);

    tree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km_vec);
    tree->Branch("Gen_dR_Kp_phi", &Gen_dR_Kp_phi_vec);
    tree->Branch("Gen_dR_Km_phi", &Gen_dR_Km_phi_vec);
    tree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi_vec);
    tree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi_vec);
    tree->Branch("Gen_dR_pi_phi", &Gen_dR_pi_phi_vec);
    tree->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds_vec);
    tree->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds_vec);
    tree->Branch("Gen_dR_phi_Ds", &Gen_dR_phi_Ds_vec);
    tree->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds_vec);

    tree->Branch("Gen_dxy_Kp_Km", &Gen_dxy_Kp_Km_vec);
    tree->Branch("Gen_dxy_Kp_pi", &Gen_dxy_Kp_pi_vec);
    tree->Branch("Gen_dxy_Km_pi", &Gen_dxy_Km_pi_vec);
    tree->Branch("Gen_dz_Kp_Km", &Gen_dz_Kp_Km_vec);
    tree->Branch("Gen_dz_Kp_pi", &Gen_dz_Kp_pi_vec);
    tree->Branch("Gen_dz_Km_pi", &Gen_dz_Km_pi_vec);

    tree->Branch("num_match_Kp", &num_match_Kp);
    tree->Branch("num_match_Km", &num_match_Km);
    tree->Branch("num_match_pi", &num_match_pi);
    tree->Branch("match_Kp_idx", &match_Kp_idx);
    tree->Branch("match_Km_idx", &match_Km_idx);
    tree->Branch("match_pi_idx", &match_pi_idx);

    tree->Branch("match_Kp_ETA", &match_Kp_ETA_vec);
    tree->Branch("match_Kp_PHI", &match_Kp_PHI_vec);
    tree->Branch("match_Kp_ORIVX_X", &match_Kp_ORIVX_X_vec);
    tree->Branch("match_Kp_ORIVX_Y", &match_Kp_ORIVX_Y_vec);
    tree->Branch("match_Kp_ORIVX_Z", &match_Kp_ORIVX_Z_vec);
    tree->Branch("match_Kp_P", &match_Kp_P_vec);
    tree->Branch("match_Kp_PT", &match_Kp_PT_vec);
    tree->Branch("match_Kp_PX", &match_Kp_PX_vec);
    tree->Branch("match_Kp_PY", &match_Kp_PY_vec);
    tree->Branch("match_Kp_PZ", &match_Kp_PZ_vec);

    tree->Branch("match_Km_ETA", &match_Km_ETA_vec);
    tree->Branch("match_Km_PHI", &match_Km_PHI_vec);
    tree->Branch("match_Km_ORIVX_X", &match_Km_ORIVX_X_vec);
    tree->Branch("match_Km_ORIVX_Y", &match_Km_ORIVX_Y_vec);
    tree->Branch("match_Km_ORIVX_Z", &match_Km_ORIVX_Z_vec);
    tree->Branch("match_Km_P", &match_Km_P_vec);
    tree->Branch("match_Km_PT", &match_Km_PT_vec);
    tree->Branch("match_Km_PX", &match_Km_PX_vec);
    tree->Branch("match_Km_PY", &match_Km_PY_vec);
    tree->Branch("match_Km_PZ", &match_Km_PZ_vec);

    tree->Branch("match_pi_ETA", &match_pi_ETA_vec);
    tree->Branch("match_pi_PHI", &match_pi_PHI_vec);
    tree->Branch("match_pi_ORIVX_X", &match_pi_ORIVX_X_vec);
    tree->Branch("match_pi_ORIVX_Y", &match_pi_ORIVX_Y_vec);
    tree->Branch("match_pi_ORIVX_Z", &match_pi_ORIVX_Z_vec);
    tree->Branch("match_pi_P", &match_pi_P_vec);
    tree->Branch("match_pi_PT", &match_pi_PT_vec);
    tree->Branch("match_pi_PX", &match_pi_PX_vec);
    tree->Branch("match_pi_PY", &match_pi_PY_vec);
    tree->Branch("match_pi_PZ", &match_pi_PZ_vec);

    tree->Branch("match_phi_ETA", &match_phi_ETA_vec);
    tree->Branch("match_phi_PHI", &match_phi_PHI_vec);
    tree->Branch("match_phi_P", &match_phi_P_vec);
    tree->Branch("match_phi_PT", &match_phi_PT_vec);
    tree->Branch("match_phi_PX", &match_phi_PX_vec);
    tree->Branch("match_phi_PY", &match_phi_PY_vec);
    tree->Branch("match_phi_PZ", &match_phi_PZ_vec);
    tree->Branch("match_phi_M", &match_phi_M_vec);

    tree->Branch("match_Ds_ETA", &match_Ds_ETA_vec);
    tree->Branch("match_Ds_PHI", &match_Ds_PHI_vec);
    tree->Branch("match_Ds_P", &match_Ds_P_vec);
    tree->Branch("match_Ds_PT", &match_Ds_PT_vec);
    tree->Branch("match_Ds_PX", &match_Ds_PX_vec);
    tree->Branch("match_Ds_PY", &match_Ds_PY_vec);
    tree->Branch("match_Ds_PZ", &match_Ds_PZ_vec);
    tree->Branch("match_Ds_M", &match_Ds_M_vec);

    tree->Branch("match_Kp_PP", &match_Kp_PP_vec);
    tree->Branch("match_Kp_PL", &match_Kp_PL_vec);
    tree->Branch("match_Km_PP", &match_Km_PP_vec);
    tree->Branch("match_Km_PL", &match_Km_PL_vec);

    tree->Branch("match_phi_PP", &match_phi_PP_vec);
    tree->Branch("match_phi_PL", &match_phi_PL_vec);
    tree->Branch("match_pi_PP", &match_pi_PP_vec);
    tree->Branch("match_pi_PL", &match_pi_PL_vec);

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

    tree->Branch("match_dxy_Kp_Km", &match_dxy_Kp_Km_vec);
    tree->Branch("match_dxy_Kp_phi", &match_dxy_Kp_phi_vec);
    tree->Branch("match_dxy_Km_phi", &match_dxy_Km_phi_vec);
    tree->Branch("match_dxy_Kp_pi", &match_dxy_Kp_pi_vec);
    tree->Branch("match_dxy_Km_pi", &match_dxy_Km_pi_vec);
    tree->Branch("match_dxy_pi_phi", &match_dxy_pi_phi_vec);
    tree->Branch("match_dxy_Kp_Ds", &match_dxy_Kp_Ds_vec);
    tree->Branch("match_dxy_Km_Ds", &match_dxy_Km_Ds_vec);
    tree->Branch("match_dxy_phi_Ds", &match_dxy_phi_Ds_vec);
    tree->Branch("match_dxy_pi_Ds", &match_dxy_pi_Ds_vec);
    
    tree->Branch("match_dz_Kp_Km", &match_dz_Kp_Km_vec);
    tree->Branch("match_dz_Kp_phi", &match_dz_Kp_phi_vec);
    tree->Branch("match_dz_Km_phi", &match_dz_Km_phi_vec);
    tree->Branch("match_dz_Kp_pi", &match_dz_Kp_pi_vec);
    tree->Branch("match_dz_Km_pi", &match_dz_Km_pi_vec);
    tree->Branch("match_dz_pi_phi", &match_dz_pi_phi_vec);
    tree->Branch("match_dz_Kp_Ds", &match_dz_Kp_Ds_vec);
    tree->Branch("match_dz_Km_Ds", &match_dz_Km_Ds_vec);
    tree->Branch("match_dz_phi_Ds", &match_dz_phi_Ds_vec);
    tree->Branch("match_dz_pi_Ds", &match_dz_pi_Ds_vec);

    tree->Branch("match_phiFit_CHI2", &match_phiFit_CHI2_vec);
    tree->Branch("match_phiFit_NDOF", &match_phiFit_NDOF_vec);
    tree->Branch("match_phiFit_CHI2NDOF", &match_phiFit_CHI2NDOF_vec);
    tree->Branch("match_phiFit_ENDVX_X", &match_phiFit_ENDVX_X_vec);
    tree->Branch("match_phiFit_ENDVX_Y", &match_phiFit_ENDVX_Y_vec);
    tree->Branch("match_phiFit_ENDVX_Z", &match_phiFit_ENDVX_Z_vec);

    tree->Branch("match_phiFit_Kp_ETA", &match_phiFit_Kp_ETA_vec);
    tree->Branch("match_phiFit_Kp_PHI", &match_phiFit_Kp_PHI_vec);
    tree->Branch("match_phiFit_Kp_P", &match_phiFit_Kp_P_vec);
    tree->Branch("match_phiFit_Kp_PT", &match_phiFit_Kp_PT_vec);
    tree->Branch("match_phiFit_Kp_PX", &match_phiFit_Kp_PX_vec);
    tree->Branch("match_phiFit_Kp_PY", &match_phiFit_Kp_PY_vec);
    tree->Branch("match_phiFit_Kp_PZ", &match_phiFit_Kp_PZ_vec);

    tree->Branch("match_phiFit_Km_ETA", &match_phiFit_Km_ETA_vec);
    tree->Branch("match_phiFit_Km_PHI", &match_phiFit_Km_PHI_vec);
    tree->Branch("match_phiFit_Km_P", &match_phiFit_Km_P_vec);
    tree->Branch("match_phiFit_Km_PT", &match_phiFit_Km_PT_vec);
    tree->Branch("match_phiFit_Km_PX", &match_phiFit_Km_PX_vec);
    tree->Branch("match_phiFit_Km_PY", &match_phiFit_Km_PY_vec);
    tree->Branch("match_phiFit_Km_PZ", &match_phiFit_Km_PZ_vec);

    tree->Branch("match_phiFit_pi_ETA", &match_phiFit_pi_ETA_vec);
    tree->Branch("match_phiFit_pi_PHI", &match_phiFit_pi_PHI_vec);
    tree->Branch("match_phiFit_pi_P", &match_phiFit_pi_P_vec);
    tree->Branch("match_phiFit_pi_PT", &match_phiFit_pi_PT_vec);
    tree->Branch("match_phiFit_pi_PX", &match_phiFit_pi_PX_vec);
    tree->Branch("match_phiFit_pi_PY", &match_phiFit_pi_PY_vec);
    tree->Branch("match_phiFit_pi_PZ", &match_phiFit_pi_PZ_vec);

    tree->Branch("match_phiFit_phi_ETA", &match_phiFit_phi_ETA_vec);
    tree->Branch("match_phiFit_phi_PHI", &match_phiFit_phi_PHI_vec);
    tree->Branch("match_phiFit_phi_P", &match_phiFit_phi_P_vec);
    tree->Branch("match_phiFit_phi_PT", &match_phiFit_phi_PT_vec);
    tree->Branch("match_phiFit_phi_PX", &match_phiFit_phi_PX_vec);
    tree->Branch("match_phiFit_phi_PY", &match_phiFit_phi_PY_vec);
    tree->Branch("match_phiFit_phi_PZ", &match_phiFit_phi_PZ_vec);
    tree->Branch("match_phiFit_phi_M", &match_phiFit_phi_M_vec);

    tree->Branch("match_phiFit_Ds_ETA", &match_phiFit_Ds_ETA_vec);
    tree->Branch("match_phiFit_Ds_PHI", &match_phiFit_Ds_PHI_vec);
    tree->Branch("match_phiFit_Ds_P", &match_phiFit_Ds_P_vec);
    tree->Branch("match_phiFit_Ds_PT", &match_phiFit_Ds_PT_vec);
    tree->Branch("match_phiFit_Ds_PX", &match_phiFit_Ds_PX_vec);
    tree->Branch("match_phiFit_Ds_PY", &match_phiFit_Ds_PY_vec);
    tree->Branch("match_phiFit_Ds_PZ", &match_phiFit_Ds_PZ_vec);
    tree->Branch("match_phiFit_Ds_M", &match_phiFit_Ds_M_vec);

    tree->Branch("match_phiFit_Kp_PP", &match_phiFit_Kp_PP_vec);
    tree->Branch("match_phiFit_Kp_PL", &match_phiFit_Kp_PL_vec);
    tree->Branch("match_phiFit_Km_PP", &match_phiFit_Km_PP_vec);
    tree->Branch("match_phiFit_Km_PL", &match_phiFit_Km_PL_vec);

    tree->Branch("match_phiFit_phi_PP", &match_phiFit_phi_PP_vec);
    tree->Branch("match_phiFit_phi_PL", &match_phiFit_phi_PL_vec);
    tree->Branch("match_phiFit_pi_PP", &match_phiFit_pi_PP_vec);
    tree->Branch("match_phiFit_pi_PL", &match_phiFit_pi_PL_vec);

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
    tree->Branch("match_DsFit_CHI2", &match_DsFit_CHI2_vec);
    tree->Branch("match_DsFit_NDOF", &match_DsFit_NDOF_vec);
    tree->Branch("match_DsFit_CHI2NDOF", &match_DsFit_CHI2NDOF_vec);
    tree->Branch("match_DsFit_ENDVX_X", &match_DsFit_ENDVX_X_vec);
    tree->Branch("match_DsFit_ENDVX_Y", &match_DsFit_ENDVX_Y_vec);
    tree->Branch("match_DsFit_ENDVX_Z", &match_DsFit_ENDVX_Z_vec);

    tree->Branch("match_DsFit_Kp_ETA", &match_DsFit_Kp_ETA_vec);
    tree->Branch("match_DsFit_Kp_PHI", &match_DsFit_Kp_PHI_vec);
    tree->Branch("match_DsFit_Kp_P", &match_DsFit_Kp_P_vec);
    tree->Branch("match_DsFit_Kp_PT", &match_DsFit_Kp_PT_vec);
    tree->Branch("match_DsFit_Kp_PX", &match_DsFit_Kp_PX_vec);
    tree->Branch("match_DsFit_Kp_PY", &match_DsFit_Kp_PY_vec);
    tree->Branch("match_DsFit_Kp_PZ", &match_DsFit_Kp_PZ_vec);

    tree->Branch("match_DsFit_Km_ETA", &match_DsFit_Km_ETA_vec);
    tree->Branch("match_DsFit_Km_PHI", &match_DsFit_Km_PHI_vec);
    tree->Branch("match_DsFit_Km_P", &match_DsFit_Km_P_vec);
    tree->Branch("match_DsFit_Km_PT", &match_DsFit_Km_PT_vec);
    tree->Branch("match_DsFit_Km_PX", &match_DsFit_Km_PX_vec);
    tree->Branch("match_DsFit_Km_PY", &match_DsFit_Km_PY_vec);
    tree->Branch("match_DsFit_Km_PZ", &match_DsFit_Km_PZ_vec);

    tree->Branch("match_DsFit_pi_ETA", &match_DsFit_pi_ETA_vec);
    tree->Branch("match_DsFit_pi_PHI", &match_DsFit_pi_PHI_vec);
    tree->Branch("match_DsFit_pi_P", &match_DsFit_pi_P_vec);
    tree->Branch("match_DsFit_pi_PT", &match_DsFit_pi_PT_vec);
    tree->Branch("match_DsFit_pi_PX", &match_DsFit_pi_PX_vec);
    tree->Branch("match_DsFit_pi_PY", &match_DsFit_pi_PY_vec);
    tree->Branch("match_DsFit_pi_PZ", &match_DsFit_pi_PZ_vec);

    tree->Branch("match_DsFit_phi_ETA", &match_DsFit_phi_ETA_vec);
    tree->Branch("match_DsFit_phi_PHI", &match_DsFit_phi_PHI_vec);
    tree->Branch("match_DsFit_phi_P", &match_DsFit_phi_P_vec);
    tree->Branch("match_DsFit_phi_PT", &match_DsFit_phi_PT_vec);
    tree->Branch("match_DsFit_phi_PX", &match_DsFit_phi_PX_vec);
    tree->Branch("match_DsFit_phi_PY", &match_DsFit_phi_PY_vec);
    tree->Branch("match_DsFit_phi_PZ", &match_DsFit_phi_PZ_vec);
    tree->Branch("match_DsFit_phi_M", &match_DsFit_phi_M_vec);

    tree->Branch("match_DsFit_Ds_ETA", &match_DsFit_Ds_ETA_vec);
    tree->Branch("match_DsFit_Ds_PHI", &match_DsFit_Ds_PHI_vec);
    tree->Branch("match_DsFit_Ds_P", &match_DsFit_Ds_P_vec);
    tree->Branch("match_DsFit_Ds_PT", &match_DsFit_Ds_PT_vec);
    tree->Branch("match_DsFit_Ds_PX", &match_DsFit_Ds_PX_vec);
    tree->Branch("match_DsFit_Ds_PY", &match_DsFit_Ds_PY_vec);
    tree->Branch("match_DsFit_Ds_PZ", &match_DsFit_Ds_PZ_vec);
    tree->Branch("match_DsFit_Ds_M", &match_DsFit_Ds_M_vec);

    tree->Branch("match_DsFit_Kp_PP", &match_DsFit_Kp_PP_vec);
    tree->Branch("match_DsFit_Kp_PL", &match_DsFit_Kp_PL_vec);
    tree->Branch("match_DsFit_Km_PP", &match_DsFit_Km_PP_vec);
    tree->Branch("match_DsFit_Km_PL", &match_DsFit_Km_PL_vec);

    tree->Branch("match_DsFit_phi_PP", &match_DsFit_phi_PP_vec);
    tree->Branch("match_DsFit_phi_PL", &match_DsFit_phi_PL_vec);
    tree->Branch("match_DsFit_pi_PP", &match_DsFit_pi_PP_vec);
    tree->Branch("match_DsFit_pi_PL", &match_DsFit_pi_PL_vec);

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

    tree->Branch("match_DsFit_Mconstraint_Ds_M", &match_DsFit_Mconstraint_Ds_M_vec);

    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

    tree->Branch("Kp_ETA", &Kp_ETA_vec);
    tree->Branch("Kp_PHI", &Kp_PHI_vec);
    tree->Branch("Kp_ORIVX_X", &Kp_ORIVX_X_vec);
    tree->Branch("Kp_ORIVX_Y", &Kp_ORIVX_Y_vec);
    tree->Branch("Kp_ORIVX_Z", &Kp_ORIVX_Z_vec);
    tree->Branch("Kp_P", &Kp_P_vec);
    tree->Branch("Kp_PT", &Kp_PT_vec);
    tree->Branch("Kp_PX", &Kp_PX_vec);
    tree->Branch("Kp_PY", &Kp_PY_vec);
    tree->Branch("Kp_PZ", &Kp_PZ_vec);

    tree->Branch("Km_ETA", &Km_ETA_vec);
    tree->Branch("Km_PHI", &Km_PHI_vec);
    tree->Branch("Km_ORIVX_X", &Km_ORIVX_X_vec);
    tree->Branch("Km_ORIVX_Y", &Km_ORIVX_Y_vec);
    tree->Branch("Km_ORIVX_Z", &Km_ORIVX_Z_vec);
    tree->Branch("Km_P", &Km_P_vec);
    tree->Branch("Km_PT", &Km_PT_vec);
    tree->Branch("Km_PX", &Km_PX_vec);
    tree->Branch("Km_PY", &Km_PY_vec);
    tree->Branch("Km_PZ", &Km_PZ_vec);

    tree->Branch("pi_ETA", &pi_ETA_vec);
    tree->Branch("pi_PHI", &pi_PHI_vec);
    tree->Branch("pi_ORIVX_X", &pi_ORIVX_X_vec);
    tree->Branch("pi_ORIVX_Y", &pi_ORIVX_Y_vec);
    tree->Branch("pi_ORIVX_Z", &pi_ORIVX_Z_vec);
    tree->Branch("pi_P", &pi_P_vec);
    tree->Branch("pi_PT", &pi_PT_vec);
    tree->Branch("pi_PX", &pi_PX_vec);
    tree->Branch("pi_PY", &pi_PY_vec);
    tree->Branch("pi_PZ", &pi_PZ_vec);

    tree->Branch("phi_ETA", &phi_ETA_vec);
    tree->Branch("phi_PHI", &phi_PHI_vec);
    tree->Branch("phi_P", &phi_P_vec);
    tree->Branch("phi_PT", &phi_PT_vec);
    tree->Branch("phi_PX", &phi_PX_vec);
    tree->Branch("phi_PY", &phi_PY_vec);
    tree->Branch("phi_PZ", &phi_PZ_vec);
    tree->Branch("phi_M", &phi_M_vec);

    tree->Branch("Ds_ETA", &Ds_ETA_vec);
    tree->Branch("Ds_PHI", &Ds_PHI_vec);
    tree->Branch("Ds_P", &Ds_P_vec);
    tree->Branch("Ds_PT", &Ds_PT_vec);
    tree->Branch("Ds_PX", &Ds_PX_vec);
    tree->Branch("Ds_PY", &Ds_PY_vec);
    tree->Branch("Ds_PZ", &Ds_PZ_vec);
    tree->Branch("Ds_M", &Ds_M_vec);

    tree->Branch("Kp_PP", &Kp_PP_vec);
    tree->Branch("Kp_PL", &Kp_PL_vec);
    tree->Branch("Km_PP", &Km_PP_vec);
    tree->Branch("Km_PL", &Km_PL_vec);

    tree->Branch("phi_PP", &phi_PP_vec);
    tree->Branch("phi_PL", &phi_PL_vec);
    tree->Branch("pi_PP", &pi_PP_vec);
    tree->Branch("pi_PL", &pi_PL_vec);

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

    tree->Branch("dxy_Kp_Km", &dxy_Kp_Km_vec);
    tree->Branch("dxy_Kp_phi", &dxy_Kp_phi_vec);
    tree->Branch("dxy_Km_phi", &dxy_Km_phi_vec);
    tree->Branch("dxy_Kp_pi", &dxy_Kp_pi_vec);
    tree->Branch("dxy_Km_pi", &dxy_Km_pi_vec);
    tree->Branch("dxy_pi_phi", &dxy_pi_phi_vec);
    tree->Branch("dxy_Kp_Ds", &dxy_Kp_Ds_vec);
    tree->Branch("dxy_Km_Ds", &dxy_Km_Ds_vec);
    tree->Branch("dxy_phi_Ds", &dxy_phi_Ds_vec);
    tree->Branch("dxy_pi_Ds", &dxy_pi_Ds_vec);
    
    tree->Branch("dz_Kp_Km", &dz_Kp_Km_vec);
    tree->Branch("dz_Kp_phi", &dz_Kp_phi_vec);
    tree->Branch("dz_Km_phi", &dz_Km_phi_vec);
    tree->Branch("dz_Kp_pi", &dz_Kp_pi_vec);
    tree->Branch("dz_Km_pi", &dz_Km_pi_vec);
    tree->Branch("dz_pi_phi", &dz_pi_phi_vec);
    tree->Branch("dz_Kp_Ds", &dz_Kp_Ds_vec);
    tree->Branch("dz_Km_Ds", &dz_Km_Ds_vec);
    tree->Branch("dz_phi_Ds", &dz_phi_Ds_vec);
    tree->Branch("dz_pi_Ds", &dz_pi_Ds_vec);

    tree->Branch("phiFit_CHI2", &phiFit_CHI2_vec);
    tree->Branch("phiFit_NDOF", &phiFit_NDOF_vec);
    tree->Branch("phiFit_CHI2NDOF", &phiFit_CHI2NDOF_vec);
    tree->Branch("phiFit_ENDVX_X", &phiFit_ENDVX_X_vec);
    tree->Branch("phiFit_ENDVX_Y", &phiFit_ENDVX_Y_vec);
    tree->Branch("phiFit_ENDVX_Z", &phiFit_ENDVX_Z_vec);

    tree->Branch("phiFit_Kp_ETA", &phiFit_Kp_ETA_vec);
    tree->Branch("phiFit_Kp_PHI", &phiFit_Kp_PHI_vec);
    tree->Branch("phiFit_Kp_P", &phiFit_Kp_P_vec);
    tree->Branch("phiFit_Kp_PT", &phiFit_Kp_PT_vec);
    tree->Branch("phiFit_Kp_PX", &phiFit_Kp_PX_vec);
    tree->Branch("phiFit_Kp_PY", &phiFit_Kp_PY_vec);
    tree->Branch("phiFit_Kp_PZ", &phiFit_Kp_PZ_vec);

    tree->Branch("phiFit_Km_ETA", &phiFit_Km_ETA_vec);
    tree->Branch("phiFit_Km_PHI", &phiFit_Km_PHI_vec);
    tree->Branch("phiFit_Km_P", &phiFit_Km_P_vec);
    tree->Branch("phiFit_Km_PT", &phiFit_Km_PT_vec);
    tree->Branch("phiFit_Km_PX", &phiFit_Km_PX_vec);
    tree->Branch("phiFit_Km_PY", &phiFit_Km_PY_vec);
    tree->Branch("phiFit_Km_PZ", &phiFit_Km_PZ_vec);

    tree->Branch("phiFit_pi_ETA", &phiFit_pi_ETA_vec);
    tree->Branch("phiFit_pi_PHI", &phiFit_pi_PHI_vec);
    tree->Branch("phiFit_pi_P", &phiFit_pi_P_vec);
    tree->Branch("phiFit_pi_PT", &phiFit_pi_PT_vec);
    tree->Branch("phiFit_pi_PX", &phiFit_pi_PX_vec);
    tree->Branch("phiFit_pi_PY", &phiFit_pi_PY_vec);
    tree->Branch("phiFit_pi_PZ", &phiFit_pi_PZ_vec);

    tree->Branch("phiFit_phi_ETA", &phiFit_phi_ETA_vec);
    tree->Branch("phiFit_phi_PHI", &phiFit_phi_PHI_vec);
    tree->Branch("phiFit_phi_P", &phiFit_phi_P_vec);
    tree->Branch("phiFit_phi_PT", &phiFit_phi_PT_vec);
    tree->Branch("phiFit_phi_PX", &phiFit_phi_PX_vec);
    tree->Branch("phiFit_phi_PY", &phiFit_phi_PY_vec);
    tree->Branch("phiFit_phi_PZ", &phiFit_phi_PZ_vec);
    tree->Branch("phiFit_phi_M", &phiFit_phi_M_vec);

    tree->Branch("phiFit_Ds_ETA", &phiFit_Ds_ETA_vec);
    tree->Branch("phiFit_Ds_PHI", &phiFit_Ds_PHI_vec);
    tree->Branch("phiFit_Ds_P", &phiFit_Ds_P_vec);
    tree->Branch("phiFit_Ds_PT", &phiFit_Ds_PT_vec);
    tree->Branch("phiFit_Ds_PX", &phiFit_Ds_PX_vec);
    tree->Branch("phiFit_Ds_PY", &phiFit_Ds_PY_vec);
    tree->Branch("phiFit_Ds_PZ", &phiFit_Ds_PZ_vec);
    tree->Branch("phiFit_Ds_M", &phiFit_Ds_M_vec);

    tree->Branch("phiFit_Kp_PP", &phiFit_Kp_PP_vec);
    tree->Branch("phiFit_Kp_PL", &phiFit_Kp_PL_vec);
    tree->Branch("phiFit_Km_PP", &phiFit_Km_PP_vec);
    tree->Branch("phiFit_Km_PL", &phiFit_Km_PL_vec);

    tree->Branch("phiFit_phi_PP", &phiFit_phi_PP_vec);
    tree->Branch("phiFit_phi_PL", &phiFit_phi_PL_vec);
    tree->Branch("phiFit_pi_PP", &phiFit_pi_PP_vec);
    tree->Branch("phiFit_pi_PL", &phiFit_pi_PL_vec);

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
    tree->Branch("DsFit_CHI2", &DsFit_CHI2_vec);
    tree->Branch("DsFit_NDOF", &DsFit_NDOF_vec);
    tree->Branch("DsFit_CHI2NDOF", &DsFit_CHI2NDOF_vec);
    tree->Branch("DsFit_ENDVX_X", &DsFit_ENDVX_X_vec);
    tree->Branch("DsFit_ENDVX_Y", &DsFit_ENDVX_Y_vec);
    tree->Branch("DsFit_ENDVX_Z", &DsFit_ENDVX_Z_vec);

    tree->Branch("DsFit_Kp_ETA", &DsFit_Kp_ETA_vec);
    tree->Branch("DsFit_Kp_PHI", &DsFit_Kp_PHI_vec);
    tree->Branch("DsFit_Kp_P", &DsFit_Kp_P_vec);
    tree->Branch("DsFit_Kp_PT", &DsFit_Kp_PT_vec);
    tree->Branch("DsFit_Kp_PX", &DsFit_Kp_PX_vec);
    tree->Branch("DsFit_Kp_PY", &DsFit_Kp_PY_vec);
    tree->Branch("DsFit_Kp_PZ", &DsFit_Kp_PZ_vec);

    tree->Branch("DsFit_Km_ETA", &DsFit_Km_ETA_vec);
    tree->Branch("DsFit_Km_PHI", &DsFit_Km_PHI_vec);
    tree->Branch("DsFit_Km_P", &DsFit_Km_P_vec);
    tree->Branch("DsFit_Km_PT", &DsFit_Km_PT_vec);
    tree->Branch("DsFit_Km_PX", &DsFit_Km_PX_vec);
    tree->Branch("DsFit_Km_PY", &DsFit_Km_PY_vec);
    tree->Branch("DsFit_Km_PZ", &DsFit_Km_PZ_vec);

    tree->Branch("DsFit_pi_ETA", &DsFit_pi_ETA_vec);
    tree->Branch("DsFit_pi_PHI", &DsFit_pi_PHI_vec);
    tree->Branch("DsFit_pi_P", &DsFit_pi_P_vec);
    tree->Branch("DsFit_pi_PT", &DsFit_pi_PT_vec);
    tree->Branch("DsFit_pi_PX", &DsFit_pi_PX_vec);
    tree->Branch("DsFit_pi_PY", &DsFit_pi_PY_vec);
    tree->Branch("DsFit_pi_PZ", &DsFit_pi_PZ_vec);

    tree->Branch("DsFit_phi_ETA", &DsFit_phi_ETA_vec);
    tree->Branch("DsFit_phi_PHI", &DsFit_phi_PHI_vec);
    tree->Branch("DsFit_phi_P", &DsFit_phi_P_vec);
    tree->Branch("DsFit_phi_PT", &DsFit_phi_PT_vec);
    tree->Branch("DsFit_phi_PX", &DsFit_phi_PX_vec);
    tree->Branch("DsFit_phi_PY", &DsFit_phi_PY_vec);
    tree->Branch("DsFit_phi_PZ", &DsFit_phi_PZ_vec);
    tree->Branch("DsFit_phi_M", &DsFit_phi_M_vec);

    tree->Branch("DsFit_Ds_ETA", &DsFit_Ds_ETA_vec);
    tree->Branch("DsFit_Ds_PHI", &DsFit_Ds_PHI_vec);
    tree->Branch("DsFit_Ds_P", &DsFit_Ds_P_vec);
    tree->Branch("DsFit_Ds_PT", &DsFit_Ds_PT_vec);
    tree->Branch("DsFit_Ds_PX", &DsFit_Ds_PX_vec);
    tree->Branch("DsFit_Ds_PY", &DsFit_Ds_PY_vec);
    tree->Branch("DsFit_Ds_PZ", &DsFit_Ds_PZ_vec);
    tree->Branch("DsFit_Ds_M", &DsFit_Ds_M_vec);

    tree->Branch("DsFit_Kp_PP", &DsFit_Kp_PP_vec);
    tree->Branch("DsFit_Kp_PL", &DsFit_Kp_PL_vec);
    tree->Branch("DsFit_Km_PP", &DsFit_Km_PP_vec);
    tree->Branch("DsFit_Km_PL", &DsFit_Km_PL_vec);

    tree->Branch("DsFit_phi_PP", &DsFit_phi_PP_vec);
    tree->Branch("DsFit_phi_PL", &DsFit_phi_PL_vec);
    tree->Branch("DsFit_pi_PP", &DsFit_pi_PP_vec);
    tree->Branch("DsFit_pi_PL", &DsFit_pi_PL_vec);

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

    tree->Branch("DsFit_Mconstraint_Ds_M", &DsFit_Mconstraint_Ds_M_vec);
    
    tree->Branch("signal_entry", &signal_entry_vec);
}

void SSTree::Gen_Reset()
{
    Gen_Kp_ETA = null;
    Gen_Kp_PHI = null;
    Gen_Kp_ORIVX_X = null;
    Gen_Kp_ORIVX_Y = null;
    Gen_Kp_ORIVX_Z = null;
    Gen_Kp_P = null;
    Gen_Kp_PT = null;
    Gen_Kp_PX = null;
    Gen_Kp_PY = null;
    Gen_Kp_PZ = null;
    Gen_Kp_PP = null;
    Gen_Kp_PL = null;

    Gen_Km_ETA = null;
    Gen_Km_PHI = null;
    Gen_Km_ORIVX_X = null;
    Gen_Km_ORIVX_Y = null;
    Gen_Km_ORIVX_Z = null;
    Gen_Km_P = null;
    Gen_Km_PT = null;
    Gen_Km_PX = null;
    Gen_Km_PY = null;
    Gen_Km_PZ = null;
    Gen_Km_PP = null;
    Gen_Km_PL = null;

    Gen_pi_ETA = null;
    Gen_pi_PHI = null;
    Gen_pi_ORIVX_X = null;
    Gen_pi_ORIVX_Y = null;
    Gen_pi_ORIVX_Z = null;
    Gen_pi_P = null;
    Gen_pi_PT = null;
    Gen_pi_PX = null;
    Gen_pi_PY = null;
    Gen_pi_PZ = null;
    Gen_pi_PP = null;
    Gen_pi_PL = null;

    Gen_phi_ETA = null;
    Gen_phi_PHI = null;
    Gen_phi_ORIVX_X = null;
    Gen_phi_ORIVX_Y = null;
    Gen_phi_ORIVX_Z = null;
    Gen_phi_P = null;
    Gen_phi_PT = null;
    Gen_phi_PX = null;
    Gen_phi_PY = null;
    Gen_phi_PZ = null;
    Gen_phi_PP = null;
    Gen_phi_PL = null;

    Gen_Ds_ETA = null;
    Gen_Ds_PHI = null;
    Gen_Ds_ORIVX_X = null;
    Gen_Ds_ORIVX_Y = null;
    Gen_Ds_ORIVX_Z = null;
    Gen_Ds_P = null;
    Gen_Ds_PT = null;
    Gen_Ds_PX = null;
    Gen_Ds_PY = null;
    Gen_Ds_PZ = null;

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

    Gen_dxy_Kp_Km = null;
    Gen_dxy_Kp_pi = null;
    Gen_dxy_Km_pi = null;
    Gen_dz_Kp_Km = null;
    Gen_dz_Kp_pi = null;
    Gen_dz_Km_pi = null;
}

void SSTree::Match_Reset()
{
    // Original info
    match_Kp_ETA = null;
    match_Kp_PHI = null;
    match_Kp_ORIVX_X = null;
    match_Kp_ORIVX_Y = null;
    match_Kp_ORIVX_Z = null;
    match_Kp_P = null;
    match_Kp_PT = null;
    match_Kp_PX = null;
    match_Kp_PY = null;
    match_Kp_PZ = null;

    match_Km_ETA = null;
    match_Km_PHI = null;
    match_Km_ORIVX_X = null;
    match_Km_ORIVX_Y = null;
    match_Km_ORIVX_Z = null;
    match_Km_P = null;
    match_Km_PT = null;
    match_Km_PX = null;
    match_Km_PY = null;
    match_Km_PZ = null;

    match_pi_ETA = null;
    match_pi_PHI = null;
    match_pi_ORIVX_X = null;
    match_pi_ORIVX_Y = null;
    match_pi_ORIVX_Z = null;
    match_pi_P = null;
    match_pi_PT = null;
    match_pi_PX = null;
    match_pi_PY = null;
    match_pi_PZ = null;

    match_phi_ETA = null;
    match_phi_PHI = null;
    match_phi_P = null;
    match_phi_PT = null;
    match_phi_PX = null;
    match_phi_PY = null;
    match_phi_PZ = null;
    match_phi_M = null;

    match_Ds_ETA = null;
    match_Ds_PHI = null;
    match_Ds_P = null;
    match_Ds_PT = null;
    match_Ds_PX = null;
    match_Ds_PY = null;
    match_Ds_PZ = null;
    match_Ds_M = null;

    match_Kp_PP = null;
    match_Kp_PL = null;
    match_Km_PP = null;
    match_Km_PL = null;

    match_phi_PP = null;
    match_phi_PL = null;
    match_pi_PP = null;
    match_pi_PL = null;

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

    match_dxy_Kp_Km = null;
    match_dxy_Kp_phi = null;
    match_dxy_Km_phi = null;
    match_dxy_Kp_pi = null;
    match_dxy_Km_pi = null;
    match_dxy_pi_phi = null;
    match_dxy_Kp_Ds = null;
    match_dxy_Km_Ds = null;
    match_dxy_phi_Ds = null;
    match_dxy_pi_Ds = null;

    match_dz_Kp_Km = null;
    match_dz_Kp_phi = null;
    match_dz_Km_phi = null;
    match_dz_Kp_pi = null;
    match_dz_Km_pi = null;
    match_dz_pi_phi = null;
    match_dz_Kp_Ds = null;
    match_dz_Km_Ds = null;
    match_dz_phi_Ds = null;
    match_dz_pi_Ds = null;

    // Fit on phi
    match_phiFit_CHI2 = null;
    match_phiFit_NDOF = null;
    match_phiFit_CHI2NDOF = null;
    match_phiFit_ENDVX_X = null;
    match_phiFit_ENDVX_Y = null;
    match_phiFit_ENDVX_Z = null;

    match_phiFit_Kp_ETA = null;
    match_phiFit_Kp_PHI = null;
    match_phiFit_Kp_P = null;
    match_phiFit_Kp_PT = null;
    match_phiFit_Kp_PX = null;
    match_phiFit_Kp_PY = null;
    match_phiFit_Kp_PZ = null;

    match_phiFit_Km_ETA = null;
    match_phiFit_Km_PHI = null;
    match_phiFit_Km_P = null;
    match_phiFit_Km_PT = null;
    match_phiFit_Km_PX = null;
    match_phiFit_Km_PY = null;
    match_phiFit_Km_PZ = null;

    match_phiFit_pi_ETA = null;
    match_phiFit_pi_PHI = null;
    match_phiFit_pi_P = null;
    match_phiFit_pi_PT = null;
    match_phiFit_pi_PX = null;
    match_phiFit_pi_PY = null;
    match_phiFit_pi_PZ = null;

    match_phiFit_phi_ETA = null;
    match_phiFit_phi_PHI = null;
    match_phiFit_phi_P = null;
    match_phiFit_phi_PT = null;
    match_phiFit_phi_PX = null;
    match_phiFit_phi_PY = null;
    match_phiFit_phi_PZ = null;
    match_phiFit_phi_M = null;

    match_phiFit_Ds_ETA = null;
    match_phiFit_Ds_PHI = null;
    match_phiFit_Ds_P = null;
    match_phiFit_Ds_PT = null;
    match_phiFit_Ds_PX = null;
    match_phiFit_Ds_PY = null;
    match_phiFit_Ds_PZ = null;
    match_phiFit_Ds_M = null;

    match_phiFit_Kp_PP = null;
    match_phiFit_Kp_PL = null;
    match_phiFit_Km_PP = null;
    match_phiFit_Km_PL = null;

    match_phiFit_phi_PP = null;
    match_phiFit_phi_PL = null;
    match_phiFit_pi_PP = null;
    match_phiFit_pi_PL = null;

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
    match_DsFit_CHI2 = null;
    match_DsFit_NDOF = null;
    match_DsFit_CHI2NDOF = null;
    match_DsFit_ENDVX_X = null;
    match_DsFit_ENDVX_Y = null;
    match_DsFit_ENDVX_Z = null;

    match_DsFit_Kp_ETA = null;
    match_DsFit_Kp_PHI = null;
    match_DsFit_Kp_P = null;
    match_DsFit_Kp_PT = null;
    match_DsFit_Kp_PX = null;
    match_DsFit_Kp_PY = null;
    match_DsFit_Kp_PZ = null;

    match_DsFit_Km_ETA = null;
    match_DsFit_Km_PHI = null;
    match_DsFit_Km_P = null;
    match_DsFit_Km_PT = null;
    match_DsFit_Km_PX = null;
    match_DsFit_Km_PY = null;
    match_DsFit_Km_PZ = null;

    match_DsFit_pi_ETA = null;
    match_DsFit_pi_PHI = null;
    match_DsFit_pi_P = null;
    match_DsFit_pi_PT = null;
    match_DsFit_pi_PX = null;
    match_DsFit_pi_PY = null;
    match_DsFit_pi_PZ = null;

    match_DsFit_phi_ETA = null;
    match_DsFit_phi_PHI = null;
    match_DsFit_phi_P = null;
    match_DsFit_phi_PT = null;
    match_DsFit_phi_PX = null;
    match_DsFit_phi_PY = null;
    match_DsFit_phi_PZ = null;
    match_DsFit_phi_M = null;

    match_DsFit_Ds_ETA = null;
    match_DsFit_Ds_PHI = null;
    match_DsFit_Ds_P = null;
    match_DsFit_Ds_PT = null;
    match_DsFit_Ds_PX = null;
    match_DsFit_Ds_PY = null;
    match_DsFit_Ds_PZ = null;
    match_DsFit_Ds_M = null;

    match_DsFit_Kp_PP = null;
    match_DsFit_Kp_PL = null;
    match_DsFit_Km_PP = null;
    match_DsFit_Km_PL = null;

    match_DsFit_phi_PP = null;
    match_DsFit_phi_PL = null;
    match_DsFit_pi_PP = null;
    match_DsFit_pi_PL = null;

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

    match_DsFit_Mconstraint_Ds_M = null;
}

void SSTree::Kp_Reset()
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

void SSTree::Km_Reset()
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

void SSTree::pi_Reset()
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
    
    signal_entry = false;
}


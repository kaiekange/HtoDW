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

    match_phiFit_phi_M_vec.push_back(match_phiFit_phi_M);
    match_phiFit_Ds_M_vec.push_back(match_phiFit_Ds_M);

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

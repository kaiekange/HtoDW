#ifndef SSTREE_H 
#define SSTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class SSTree 
{
    public:

        SSTree(TTree *tree_);

        TTree *tree;

        void Init();
        void Gen_Reset();
        void Match_Reset();
        /* void Kp_Reset(); */
        /* void Km_Reset(); */
        /* void phi_Reset(); */
        /* void pi_Reset(); */
        /* void Ds_Reset(); */
        void Gen_Fill_Vector();
        void Match_Fill_Vector();
        void CreateBranches();

        // Gen particles

        int num_Gen_Kp;
        int num_Gen_Km;
        int num_Gen_pi;
        int num_Gen_phi;
        int num_Gen_Ds;

        double Gen_Kp_ETA;
        double Gen_Kp_PHI;
        double Gen_Kp_ORIVX_X;
        double Gen_Kp_ORIVX_Y;
        double Gen_Kp_ORIVX_Z;
        double Gen_Kp_P;
        double Gen_Kp_PT;
        double Gen_Kp_PX;
        double Gen_Kp_PY;
        double Gen_Kp_PZ;
        double Gen_Kp_PP;
        double Gen_Kp_PL;

        double Gen_Km_ETA;
        double Gen_Km_PHI;
        double Gen_Km_ORIVX_X;
        double Gen_Km_ORIVX_Y;
        double Gen_Km_ORIVX_Z;
        double Gen_Km_P;
        double Gen_Km_PT;
        double Gen_Km_PX;
        double Gen_Km_PY;
        double Gen_Km_PZ;
        double Gen_Km_PP;
        double Gen_Km_PL;

        double Gen_pi_ETA;
        double Gen_pi_PHI;
        double Gen_pi_ORIVX_X;
        double Gen_pi_ORIVX_Y;
        double Gen_pi_ORIVX_Z;
        double Gen_pi_P;
        double Gen_pi_PT;
        double Gen_pi_PX;
        double Gen_pi_PY;
        double Gen_pi_PZ;
        double Gen_pi_PP;
        double Gen_pi_PL;

        double Gen_phi_ETA;
        double Gen_phi_PHI;
        double Gen_phi_ORIVX_X;
        double Gen_phi_ORIVX_Y;
        double Gen_phi_ORIVX_Z;
        double Gen_phi_P;
        double Gen_phi_PT;
        double Gen_phi_PX;
        double Gen_phi_PY;
        double Gen_phi_PZ;
        double Gen_phi_PP;
        double Gen_phi_PL;

        double Gen_Ds_ETA;
        double Gen_Ds_PHI;
        double Gen_Ds_ORIVX_X;
        double Gen_Ds_ORIVX_Y;
        double Gen_Ds_ORIVX_Z;
        double Gen_Ds_P;
        double Gen_Ds_PT;
        double Gen_Ds_PX;
        double Gen_Ds_PY;
        double Gen_Ds_PZ;

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

        double Gen_dxy_Kp_Km;
        double Gen_dxy_Kp_pi;
        double Gen_dxy_Km_pi;
        double Gen_dz_Kp_Km;
        double Gen_dz_Kp_pi;
        double Gen_dz_Km_pi;

        // Matched particles
        int num_match_Kp;
        int num_match_Km;
        int num_match_pi;
        int match_Kp_idx;
        int match_Km_idx;
        int match_pi_idx;

        // Original info
        double match_Kp_ETA;
        double match_Kp_PHI;
        double match_Kp_ORIVX_X;
        double match_Kp_ORIVX_Y;
        double match_Kp_ORIVX_Z;
        double match_Kp_P;
        double match_Kp_PT;
        double match_Kp_PX;
        double match_Kp_PY;
        double match_Kp_PZ;

        double match_Km_ETA;
        double match_Km_PHI;
        double match_Km_ORIVX_X;
        double match_Km_ORIVX_Y;
        double match_Km_ORIVX_Z;
        double match_Km_P;
        double match_Km_PT;
        double match_Km_PX;
        double match_Km_PY;
        double match_Km_PZ;

        double match_pi_ETA;
        double match_pi_PHI;
        double match_pi_ORIVX_X;
        double match_pi_ORIVX_Y;
        double match_pi_ORIVX_Z;
        double match_pi_P;
        double match_pi_PT;
        double match_pi_PX;
        double match_pi_PY;
        double match_pi_PZ;

        double match_phi_ETA;
        double match_phi_PHI;
        double match_phi_P;
        double match_phi_PT;
        double match_phi_PX;
        double match_phi_PY;
        double match_phi_PZ;
        double match_phi_M;

        double match_Ds_ETA;
        double match_Ds_PHI;
        double match_Ds_P;
        double match_Ds_PT;
        double match_Ds_PX;
        double match_Ds_PY;
        double match_Ds_PZ;
        double match_Ds_M;

        double match_Kp_PP;
        double match_Kp_PL;
        double match_Km_PP;
        double match_Km_PL;

        double match_phi_PP;
        double match_phi_PL;
        double match_pi_PP;
        double match_pi_PL;

        double match_dR_Kp_Km;
        double match_dR_Kp_phi;
        double match_dR_Km_phi;
        double match_dR_Kp_pi;
        double match_dR_Km_pi;
        double match_dR_pi_phi;
        double match_dR_Kp_Ds;
        double match_dR_Km_Ds;
        double match_dR_phi_Ds;
        double match_dR_pi_Ds;

        double match_dxy_Kp_Km;
        double match_dxy_Kp_pi;
        double match_dxy_Km_pi;
        double match_dxy_Kp_phi;
        double match_dxy_Km_phi;
        double match_dxy_pi_phi;
        double match_dxy_Kp_Ds;
        double match_dxy_Km_Ds;
        double match_dxy_pi_Ds;
        double match_dxy_phi_Ds;

        double match_dz_Kp_Km;
        double match_dz_Kp_pi;
        double match_dz_Km_pi;
        double match_dz_Kp_phi;
        double match_dz_Km_phi;
        double match_dz_pi_phi;
        double match_dz_Kp_Ds;
        double match_dz_Km_Ds;
        double match_dz_pi_Ds;
        double match_dz_phi_Ds;

        // Fit on phi
        double match_phiFit_CHI2;
        double match_phiFit_NDOF;
        double match_phiFit_CHI2NDOF;
        double match_phiFit_ENDVX_X;
        double match_phiFit_ENDVX_Y;
        double match_phiFit_ENDVX_Z;

        double match_phiFit_Kp_ETA;
        double match_phiFit_Kp_PHI;
        double match_phiFit_Kp_P;
        double match_phiFit_Kp_PT;
        double match_phiFit_Kp_PX;
        double match_phiFit_Kp_PY;
        double match_phiFit_Kp_PZ;

        double match_phiFit_Km_ETA;
        double match_phiFit_Km_PHI;
        double match_phiFit_Km_P;
        double match_phiFit_Km_PT;
        double match_phiFit_Km_PX;
        double match_phiFit_Km_PY;
        double match_phiFit_Km_PZ;

        double match_phiFit_pi_ETA;
        double match_phiFit_pi_PHI;
        double match_phiFit_pi_P;
        double match_phiFit_pi_PT;
        double match_phiFit_pi_PX;
        double match_phiFit_pi_PY;
        double match_phiFit_pi_PZ;

        double match_phiFit_phi_ETA;
        double match_phiFit_phi_PHI;
        double match_phiFit_phi_P;
        double match_phiFit_phi_PT;
        double match_phiFit_phi_PX;
        double match_phiFit_phi_PY;
        double match_phiFit_phi_PZ;
        double match_phiFit_phi_M;

        double match_phiFit_Ds_ETA;
        double match_phiFit_Ds_PHI;
        double match_phiFit_Ds_P;
        double match_phiFit_Ds_PT;
        double match_phiFit_Ds_PX;
        double match_phiFit_Ds_PY;
        double match_phiFit_Ds_PZ;
        double match_phiFit_Ds_M;

        double match_phiFit_Kp_PP;
        double match_phiFit_Kp_PL;
        double match_phiFit_Km_PP;
        double match_phiFit_Km_PL;

        double match_phiFit_phi_PP;
        double match_phiFit_phi_PL;
        double match_phiFit_pi_PP;
        double match_phiFit_pi_PL;

        double match_phiFit_dR_Kp_Km;
        double match_phiFit_dR_Kp_phi;
        double match_phiFit_dR_Km_phi;
        double match_phiFit_dR_Kp_pi;
        double match_phiFit_dR_Km_pi;
        double match_phiFit_dR_pi_phi;
        double match_phiFit_dR_Kp_Ds;
        double match_phiFit_dR_Km_Ds;
        double match_phiFit_dR_phi_Ds;
        double match_phiFit_dR_pi_Ds;

        // Fit on Ds 
        double match_DsFit_CHI2;
        double match_DsFit_NDOF;
        double match_DsFit_CHI2NDOF;
        double match_DsFit_ENDVX_X;
        double match_DsFit_ENDVX_Y;
        double match_DsFit_ENDVX_Z;

        double match_DsFit_Kp_ETA;
        double match_DsFit_Kp_PHI;
        double match_DsFit_Kp_P;
        double match_DsFit_Kp_PT;
        double match_DsFit_Kp_PX;
        double match_DsFit_Kp_PY;
        double match_DsFit_Kp_PZ;

        double match_DsFit_Km_ETA;
        double match_DsFit_Km_PHI;
        double match_DsFit_Km_P;
        double match_DsFit_Km_PT;
        double match_DsFit_Km_PX;
        double match_DsFit_Km_PY;
        double match_DsFit_Km_PZ;

        double match_DsFit_pi_ETA;
        double match_DsFit_pi_PHI;
        double match_DsFit_pi_P;
        double match_DsFit_pi_PT;
        double match_DsFit_pi_PX;
        double match_DsFit_pi_PY;
        double match_DsFit_pi_PZ;

        double match_DsFit_phi_ETA;
        double match_DsFit_phi_PHI;
        double match_DsFit_phi_P;
        double match_DsFit_phi_PT;
        double match_DsFit_phi_PX;
        double match_DsFit_phi_PY;
        double match_DsFit_phi_PZ;
        double match_DsFit_phi_M;

        double match_DsFit_Ds_ETA;
        double match_DsFit_Ds_PHI;
        double match_DsFit_Ds_P;
        double match_DsFit_Ds_PT;
        double match_DsFit_Ds_PX;
        double match_DsFit_Ds_PY;
        double match_DsFit_Ds_PZ;
        double match_DsFit_Ds_M;

        double match_DsFit_Kp_PP;
        double match_DsFit_Kp_PL;
        double match_DsFit_Km_PP;
        double match_DsFit_Km_PL;

        double match_DsFit_phi_PP;
        double match_DsFit_phi_PL;
        double match_DsFit_pi_PP;
        double match_DsFit_pi_PL;

        double match_DsFit_dR_Kp_Km;
        double match_DsFit_dR_Kp_phi;
        double match_DsFit_dR_Km_phi;
        double match_DsFit_dR_Kp_pi;
        double match_DsFit_dR_Km_pi;
        double match_DsFit_dR_pi_phi;
        double match_DsFit_dR_Kp_Ds;
        double match_DsFit_dR_Km_Ds;
        double match_DsFit_dR_phi_Ds;
        double match_DsFit_dR_pi_Ds;

        double match_DsFit_Mconstraint_Ds_M;

        std::vector<double> Gen_Kp_ETA_vec;
        std::vector<double> Gen_Kp_PHI_vec;
        std::vector<double> Gen_Kp_ORIVX_X_vec;
        std::vector<double> Gen_Kp_ORIVX_Y_vec;
        std::vector<double> Gen_Kp_ORIVX_Z_vec;
        std::vector<double> Gen_Kp_P_vec;
        std::vector<double> Gen_Kp_PT_vec;
        std::vector<double> Gen_Kp_PX_vec;
        std::vector<double> Gen_Kp_PY_vec;
        std::vector<double> Gen_Kp_PZ_vec;
        std::vector<double> Gen_Kp_PP_vec;
        std::vector<double> Gen_Kp_PL_vec;

        std::vector<double> Gen_Km_ETA_vec;
        std::vector<double> Gen_Km_PHI_vec;
        std::vector<double> Gen_Km_ORIVX_X_vec;
        std::vector<double> Gen_Km_ORIVX_Y_vec;
        std::vector<double> Gen_Km_ORIVX_Z_vec;
        std::vector<double> Gen_Km_P_vec;
        std::vector<double> Gen_Km_PT_vec;
        std::vector<double> Gen_Km_PX_vec;
        std::vector<double> Gen_Km_PY_vec;
        std::vector<double> Gen_Km_PZ_vec;
        std::vector<double> Gen_Km_PP_vec;
        std::vector<double> Gen_Km_PL_vec;

        std::vector<double> Gen_pi_ETA_vec;
        std::vector<double> Gen_pi_PHI_vec;
        std::vector<double> Gen_pi_ORIVX_X_vec;
        std::vector<double> Gen_pi_ORIVX_Y_vec;
        std::vector<double> Gen_pi_ORIVX_Z_vec;
        std::vector<double> Gen_pi_P_vec;
        std::vector<double> Gen_pi_PT_vec;
        std::vector<double> Gen_pi_PX_vec;
        std::vector<double> Gen_pi_PY_vec;
        std::vector<double> Gen_pi_PZ_vec;
        std::vector<double> Gen_pi_PP_vec;
        std::vector<double> Gen_pi_PL_vec;

        std::vector<double> Gen_phi_ETA_vec;
        std::vector<double> Gen_phi_PHI_vec;
        std::vector<double> Gen_phi_ORIVX_X_vec;
        std::vector<double> Gen_phi_ORIVX_Y_vec;
        std::vector<double> Gen_phi_ORIVX_Z_vec;
        std::vector<double> Gen_phi_P_vec;
        std::vector<double> Gen_phi_PT_vec;
        std::vector<double> Gen_phi_PX_vec;
        std::vector<double> Gen_phi_PY_vec;
        std::vector<double> Gen_phi_PZ_vec;
        std::vector<double> Gen_phi_PP_vec;
        std::vector<double> Gen_phi_PL_vec;

        std::vector<double> Gen_Ds_ETA_vec;
        std::vector<double> Gen_Ds_PHI_vec;
        std::vector<double> Gen_Ds_ORIVX_X_vec;
        std::vector<double> Gen_Ds_ORIVX_Y_vec;
        std::vector<double> Gen_Ds_ORIVX_Z_vec;
        std::vector<double> Gen_Ds_P_vec;
        std::vector<double> Gen_Ds_PT_vec;
        std::vector<double> Gen_Ds_PX_vec;
        std::vector<double> Gen_Ds_PY_vec;
        std::vector<double> Gen_Ds_PZ_vec;

        std::vector<double> Gen_dR_Kp_Km_vec;
        std::vector<double> Gen_dR_Kp_phi_vec;
        std::vector<double> Gen_dR_Km_phi_vec;
        std::vector<double> Gen_dR_Kp_pi_vec;
        std::vector<double> Gen_dR_Km_pi_vec;
        std::vector<double> Gen_dR_pi_phi_vec;
        std::vector<double> Gen_dR_Kp_Ds_vec;
        std::vector<double> Gen_dR_Km_Ds_vec;
        std::vector<double> Gen_dR_phi_Ds_vec;
        std::vector<double> Gen_dR_pi_Ds_vec;

        std::vector<double> Gen_dxy_Kp_Km_vec;
        std::vector<double> Gen_dxy_Kp_pi_vec;
        std::vector<double> Gen_dxy_Km_pi_vec;
        std::vector<double> Gen_dz_Kp_Km_vec;
        std::vector<double> Gen_dz_Kp_pi_vec;
        std::vector<double> Gen_dz_Km_pi_vec;

        // Original info
        std::vector<double> match_Kp_ETA_vec;
        std::vector<double> match_Kp_PHI_vec;
        std::vector<double> match_Kp_ORIVX_X_vec;
        std::vector<double> match_Kp_ORIVX_Y_vec;
        std::vector<double> match_Kp_ORIVX_Z_vec;
        std::vector<double> match_Kp_P_vec;
        std::vector<double> match_Kp_PT_vec;
        std::vector<double> match_Kp_PX_vec;
        std::vector<double> match_Kp_PY_vec;
        std::vector<double> match_Kp_PZ_vec;

        std::vector<double> match_Km_ETA_vec;
        std::vector<double> match_Km_PHI_vec;
        std::vector<double> match_Km_ORIVX_X_vec;
        std::vector<double> match_Km_ORIVX_Y_vec;
        std::vector<double> match_Km_ORIVX_Z_vec;
        std::vector<double> match_Km_P_vec;
        std::vector<double> match_Km_PT_vec;
        std::vector<double> match_Km_PX_vec;
        std::vector<double> match_Km_PY_vec;
        std::vector<double> match_Km_PZ_vec;

        std::vector<double> match_pi_ETA_vec;
        std::vector<double> match_pi_PHI_vec;
        std::vector<double> match_pi_ORIVX_X_vec;
        std::vector<double> match_pi_ORIVX_Y_vec;
        std::vector<double> match_pi_ORIVX_Z_vec;
        std::vector<double> match_pi_P_vec;
        std::vector<double> match_pi_PT_vec;
        std::vector<double> match_pi_PX_vec;
        std::vector<double> match_pi_PY_vec;
        std::vector<double> match_pi_PZ_vec;

        std::vector<double> match_phi_ETA_vec;
        std::vector<double> match_phi_PHI_vec;
        std::vector<double> match_phi_P_vec;
        std::vector<double> match_phi_PT_vec;
        std::vector<double> match_phi_PX_vec;
        std::vector<double> match_phi_PY_vec;
        std::vector<double> match_phi_PZ_vec;
        std::vector<double> match_phi_M_vec;

        std::vector<double> match_Ds_ETA_vec;
        std::vector<double> match_Ds_PHI_vec;
        std::vector<double> match_Ds_P_vec;
        std::vector<double> match_Ds_PT_vec;
        std::vector<double> match_Ds_PX_vec;
        std::vector<double> match_Ds_PY_vec;
        std::vector<double> match_Ds_PZ_vec;
        std::vector<double> match_Ds_M_vec;

        std::vector<double> match_Kp_PP_vec;
        std::vector<double> match_Kp_PL_vec;
        std::vector<double> match_Km_PP_vec;
        std::vector<double> match_Km_PL_vec;

        std::vector<double> match_phi_PP_vec;
        std::vector<double> match_phi_PL_vec;
        std::vector<double> match_pi_PP_vec;
        std::vector<double> match_pi_PL_vec;

        std::vector<double> match_dR_Kp_Km_vec;
        std::vector<double> match_dR_Kp_phi_vec;
        std::vector<double> match_dR_Km_phi_vec;
        std::vector<double> match_dR_Kp_pi_vec;
        std::vector<double> match_dR_Km_pi_vec;
        std::vector<double> match_dR_pi_phi_vec;
        std::vector<double> match_dR_Kp_Ds_vec;
        std::vector<double> match_dR_Km_Ds_vec;
        std::vector<double> match_dR_phi_Ds_vec;
        std::vector<double> match_dR_pi_Ds_vec;

        std::vector<double> match_dxy_Kp_Km_vec;
        std::vector<double> match_dxy_Kp_phi_vec;
        std::vector<double> match_dxy_Km_phi_vec;
        std::vector<double> match_dxy_Kp_pi_vec;
        std::vector<double> match_dxy_Km_pi_vec;
        std::vector<double> match_dxy_pi_phi_vec;
        std::vector<double> match_dxy_Kp_Ds_vec;
        std::vector<double> match_dxy_Km_Ds_vec;
        std::vector<double> match_dxy_phi_Ds_vec;
        std::vector<double> match_dxy_pi_Ds_vec;

        std::vector<double> match_dz_Kp_Km_vec;
        std::vector<double> match_dz_Kp_phi_vec;
        std::vector<double> match_dz_Km_phi_vec;
        std::vector<double> match_dz_Kp_pi_vec;
        std::vector<double> match_dz_Km_pi_vec;
        std::vector<double> match_dz_pi_phi_vec;
        std::vector<double> match_dz_Kp_Ds_vec;
        std::vector<double> match_dz_Km_Ds_vec;
        std::vector<double> match_dz_phi_Ds_vec;
        std::vector<double> match_dz_pi_Ds_vec;

        // Fit on phi
        std::vector<double> match_phiFit_CHI2_vec;
        std::vector<double> match_phiFit_NDOF_vec;
        std::vector<double> match_phiFit_CHI2NDOF_vec;
        std::vector<double> match_phiFit_ENDVX_X_vec;
        std::vector<double> match_phiFit_ENDVX_Y_vec;
        std::vector<double> match_phiFit_ENDVX_Z_vec;

        std::vector<double> match_phiFit_Kp_ETA_vec;
        std::vector<double> match_phiFit_Kp_PHI_vec;
        std::vector<double> match_phiFit_Kp_P_vec;
        std::vector<double> match_phiFit_Kp_PT_vec;
        std::vector<double> match_phiFit_Kp_PX_vec;
        std::vector<double> match_phiFit_Kp_PY_vec;
        std::vector<double> match_phiFit_Kp_PZ_vec;

        std::vector<double> match_phiFit_Km_ETA_vec;
        std::vector<double> match_phiFit_Km_PHI_vec;
        std::vector<double> match_phiFit_Km_P_vec;
        std::vector<double> match_phiFit_Km_PT_vec;
        std::vector<double> match_phiFit_Km_PX_vec;
        std::vector<double> match_phiFit_Km_PY_vec;
        std::vector<double> match_phiFit_Km_PZ_vec;

        std::vector<double> match_phiFit_pi_ETA_vec;
        std::vector<double> match_phiFit_pi_PHI_vec;
        std::vector<double> match_phiFit_pi_P_vec;
        std::vector<double> match_phiFit_pi_PT_vec;
        std::vector<double> match_phiFit_pi_PX_vec;
        std::vector<double> match_phiFit_pi_PY_vec;
        std::vector<double> match_phiFit_pi_PZ_vec;

        std::vector<double> match_phiFit_phi_ETA_vec;
        std::vector<double> match_phiFit_phi_PHI_vec;
        std::vector<double> match_phiFit_phi_P_vec;
        std::vector<double> match_phiFit_phi_PT_vec;
        std::vector<double> match_phiFit_phi_PX_vec;
        std::vector<double> match_phiFit_phi_PY_vec;
        std::vector<double> match_phiFit_phi_PZ_vec;
        std::vector<double> match_phiFit_phi_M_vec;

        std::vector<double> match_phiFit_Ds_ETA_vec;
        std::vector<double> match_phiFit_Ds_PHI_vec;
        std::vector<double> match_phiFit_Ds_P_vec;
        std::vector<double> match_phiFit_Ds_PT_vec;
        std::vector<double> match_phiFit_Ds_PX_vec;
        std::vector<double> match_phiFit_Ds_PY_vec;
        std::vector<double> match_phiFit_Ds_PZ_vec;
        std::vector<double> match_phiFit_Ds_M_vec;

        std::vector<double> match_phiFit_Kp_PP_vec;
        std::vector<double> match_phiFit_Kp_PL_vec;
        std::vector<double> match_phiFit_Km_PP_vec;
        std::vector<double> match_phiFit_Km_PL_vec;

        std::vector<double> match_phiFit_phi_PP_vec;
        std::vector<double> match_phiFit_phi_PL_vec;
        std::vector<double> match_phiFit_pi_PP_vec;
        std::vector<double> match_phiFit_pi_PL_vec;

        std::vector<double> match_phiFit_dR_Kp_Km_vec;
        std::vector<double> match_phiFit_dR_Kp_phi_vec;
        std::vector<double> match_phiFit_dR_Km_phi_vec;
        std::vector<double> match_phiFit_dR_Kp_pi_vec;
        std::vector<double> match_phiFit_dR_Km_pi_vec;
        std::vector<double> match_phiFit_dR_pi_phi_vec;
        std::vector<double> match_phiFit_dR_Kp_Ds_vec;
        std::vector<double> match_phiFit_dR_Km_Ds_vec;
        std::vector<double> match_phiFit_dR_phi_Ds_vec;
        std::vector<double> match_phiFit_dR_pi_Ds_vec;

        // Fit on Ds 
        std::vector<double> match_DsFit_CHI2_vec;
        std::vector<double> match_DsFit_NDOF_vec;
        std::vector<double> match_DsFit_CHI2NDOF_vec;
        std::vector<double> match_DsFit_ENDVX_X_vec;
        std::vector<double> match_DsFit_ENDVX_Y_vec;
        std::vector<double> match_DsFit_ENDVX_Z_vec;

        std::vector<double> match_DsFit_Kp_ETA_vec;
        std::vector<double> match_DsFit_Kp_PHI_vec;
        std::vector<double> match_DsFit_Kp_P_vec;
        std::vector<double> match_DsFit_Kp_PT_vec;
        std::vector<double> match_DsFit_Kp_PX_vec;
        std::vector<double> match_DsFit_Kp_PY_vec;
        std::vector<double> match_DsFit_Kp_PZ_vec;

        std::vector<double> match_DsFit_Km_ETA_vec;
        std::vector<double> match_DsFit_Km_PHI_vec;
        std::vector<double> match_DsFit_Km_P_vec;
        std::vector<double> match_DsFit_Km_PT_vec;
        std::vector<double> match_DsFit_Km_PX_vec;
        std::vector<double> match_DsFit_Km_PY_vec;
        std::vector<double> match_DsFit_Km_PZ_vec;

        std::vector<double> match_DsFit_pi_ETA_vec;
        std::vector<double> match_DsFit_pi_PHI_vec;
        std::vector<double> match_DsFit_pi_P_vec;
        std::vector<double> match_DsFit_pi_PT_vec;
        std::vector<double> match_DsFit_pi_PX_vec;
        std::vector<double> match_DsFit_pi_PY_vec;
        std::vector<double> match_DsFit_pi_PZ_vec;

        std::vector<double> match_DsFit_phi_ETA_vec;
        std::vector<double> match_DsFit_phi_PHI_vec;
        std::vector<double> match_DsFit_phi_P_vec;
        std::vector<double> match_DsFit_phi_PT_vec;
        std::vector<double> match_DsFit_phi_PX_vec;
        std::vector<double> match_DsFit_phi_PY_vec;
        std::vector<double> match_DsFit_phi_PZ_vec;
        std::vector<double> match_DsFit_phi_M_vec;

        std::vector<double> match_DsFit_Ds_ETA_vec;
        std::vector<double> match_DsFit_Ds_PHI_vec;
        std::vector<double> match_DsFit_Ds_P_vec;
        std::vector<double> match_DsFit_Ds_PT_vec;
        std::vector<double> match_DsFit_Ds_PX_vec;
        std::vector<double> match_DsFit_Ds_PY_vec;
        std::vector<double> match_DsFit_Ds_PZ_vec;
        std::vector<double> match_DsFit_Ds_M_vec;

        std::vector<double> match_DsFit_Kp_PP_vec;
        std::vector<double> match_DsFit_Kp_PL_vec;
        std::vector<double> match_DsFit_Km_PP_vec;
        std::vector<double> match_DsFit_Km_PL_vec;

        std::vector<double> match_DsFit_phi_PP_vec;
        std::vector<double> match_DsFit_phi_PL_vec;
        std::vector<double> match_DsFit_pi_PP_vec;
        std::vector<double> match_DsFit_pi_PL_vec;

        std::vector<double> match_DsFit_dR_Kp_Km_vec;
        std::vector<double> match_DsFit_dR_Kp_phi_vec;
        std::vector<double> match_DsFit_dR_Km_phi_vec;
        std::vector<double> match_DsFit_dR_Kp_pi_vec;
        std::vector<double> match_DsFit_dR_Km_pi_vec;
        std::vector<double> match_DsFit_dR_pi_phi_vec;
        std::vector<double> match_DsFit_dR_Kp_Ds_vec;
        std::vector<double> match_DsFit_dR_Km_Ds_vec;
        std::vector<double> match_DsFit_dR_phi_Ds_vec;
        std::vector<double> match_DsFit_dR_pi_Ds_vec;

        std::vector<double> match_DsFit_Mconstraint_Ds_M_vec;


};

#endif

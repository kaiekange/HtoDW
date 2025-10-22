#ifndef RECOBESTTREE_H 
#define RECOBESTTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class RecoBestTree 
{
    public:

        RecoBestTree(TTree *tree_);

        TTree *tree;

        void Init();
        void Kp_Reset();
        void Km_Reset();
        void pi_Reset();
        void BS_Reset();
        void PV_Reset();
        void Best_Fill_Vector(int idxmax);
        void Fill_Vector();
        void PV_withBS_Fill_Vector(); 
        void PV_noBS_Fill_Vector(); 
        void CreateBranches();

        // Gen particles

        int BS_type;
        double BS_X0;
        double BS_Y0;
        double BS_Z0;
        double BS_SigmaZ;
        double BS_dXdZ;
        double BS_dYdZ;
        double BS_BWX;
        double BS_BWY;
        double BS_X0ERR;
        double BS_Y0ERR;
        double BS_Z0ERR;
        double BS_SigmaZ0ERR;
        double BS_dXdZERR;
        double BS_dYdZERR;
        double BS_BWXERR;
        double BS_BWYERR;
        double BS_EmitX;
        double BS_EmitY;
        double BS_BetaStar;

        bool PV_withBS_IsValid;
        bool PV_withBS_IsFake;
        double PV_withBS_CHI2;
        double PV_withBS_NDOF;
        double PV_withBS_CHI2NDOF;
        double PV_withBS_X;
        double PV_withBS_Y;
        double PV_withBS_Z;
        double PV_withBS_XERR;
        double PV_withBS_YERR;
        double PV_withBS_ZERR;
        bool PV_noBS_IsValid;
        bool PV_noBS_IsFake;
        double PV_noBS_CHI2;
        double PV_noBS_NDOF;
        double PV_noBS_CHI2NDOF;
        double PV_noBS_X;
        double PV_noBS_Y;
        double PV_noBS_Z;
        double PV_noBS_XERR;
        double PV_noBS_YERR;
        double PV_noBS_ZERR;

        int num_reco_phi;
        int num_reco_Ds;

        double Kp_ETA;
        double Kp_PHI;
        double Kp_ORIVX_X;
        double Kp_ORIVX_Y;
        double Kp_ORIVX_Z;
        double Kp_P;
        double Kp_PT;
        double Kp_PX;
        double Kp_PY;
        double Kp_PZ;

        double Km_ETA;
        double Km_PHI;
        double Km_ORIVX_X;
        double Km_ORIVX_Y;
        double Km_ORIVX_Z;
        double Km_P;
        double Km_PT;
        double Km_PX;
        double Km_PY;
        double Km_PZ;

        double pi_ETA;
        double pi_PHI;
        double pi_ORIVX_X;
        double pi_ORIVX_Y;
        double pi_ORIVX_Z;
        double pi_P;
        double pi_PT;
        double pi_PX;
        double pi_PY;
        double pi_PZ;

        double phi_ETA;
        double phi_PHI;
        double phi_P;
        double phi_PT;
        double phi_PX;
        double phi_PY;
        double phi_PZ;
        double phi_M;

        double Ds_ETA;
        double Ds_PHI;
        double Ds_P;
        double Ds_PT;
        double Ds_PX;
        double Ds_PY;
        double Ds_PZ;
        double Ds_M;

        double Kp_PP;
        double Kp_PL;
        double Km_PP;
        double Km_PL;

        double phi_PP;
        double phi_PL;
        double pi_PP;
        double pi_PL;

        double dR_Kp_Km;
        double dR_Kp_phi;
        double dR_Km_phi;
        double dR_Kp_pi;
        double dR_Km_pi;
        double dR_pi_phi;
        double dR_Kp_Ds;
        double dR_Km_Ds;
        double dR_phi_Ds;
        double dR_pi_Ds;

        double dxy_Kp_Km;
        double dxy_Kp_pi;
        double dxy_Km_pi;
        double dxy_Kp_phi;
        double dxy_Km_phi;
        double dxy_pi_phi;
        double dxy_Kp_Ds;
        double dxy_Km_Ds;
        double dxy_pi_Ds;
        double dxy_phi_Ds;

        double dz_Kp_Km;
        double dz_Kp_pi;
        double dz_Km_pi;
        double dz_Kp_phi;
        double dz_Km_phi;
        double dz_pi_phi;
        double dz_Kp_Ds;
        double dz_Km_Ds;
        double dz_pi_Ds;
        double dz_phi_Ds;

        // Fit on phi
        double phiFit_CHI2;
        double phiFit_NDOF;
        double phiFit_CHI2NDOF;
        double phiFit_ENDVX_X;
        double phiFit_ENDVX_Y;
        double phiFit_ENDVX_Z;
        double phiFit_ENDVX_XERR;
        double phiFit_ENDVX_YERR;
        double phiFit_ENDVX_ZERR;

        double phiFit_Kp_ETA;
        double phiFit_Kp_PHI;
        double phiFit_Kp_P;
        double phiFit_Kp_PT;
        double phiFit_Kp_PX;
        double phiFit_Kp_PY;
        double phiFit_Kp_PZ;

        double phiFit_Km_ETA;
        double phiFit_Km_PHI;
        double phiFit_Km_P;
        double phiFit_Km_PT;
        double phiFit_Km_PX;
        double phiFit_Km_PY;
        double phiFit_Km_PZ;

        double phiFit_pi_ETA;
        double phiFit_pi_PHI;
        double phiFit_pi_P;
        double phiFit_pi_PT;
        double phiFit_pi_PX;
        double phiFit_pi_PY;
        double phiFit_pi_PZ;

        double phiFit_phi_ETA;
        double phiFit_phi_PHI;
        double phiFit_phi_P;
        double phiFit_phi_PT;
        double phiFit_phi_PX;
        double phiFit_phi_PY;
        double phiFit_phi_PZ;
        double phiFit_phi_M;

        double phiFit_Ds_ETA;
        double phiFit_Ds_PHI;
        double phiFit_Ds_P;
        double phiFit_Ds_PT;
        double phiFit_Ds_PX;
        double phiFit_Ds_PY;
        double phiFit_Ds_PZ;
        double phiFit_Ds_M;

        double phiFit_Kp_PP;
        double phiFit_Kp_PL;
        double phiFit_Km_PP;
        double phiFit_Km_PL;

        double phiFit_phi_PP;
        double phiFit_phi_PL;
        double phiFit_pi_PP;
        double phiFit_pi_PL;

        double phiFit_dR_Kp_Km;
        double phiFit_dR_Kp_phi;
        double phiFit_dR_Km_phi;
        double phiFit_dR_Kp_pi;
        double phiFit_dR_Km_pi;
        double phiFit_dR_pi_phi;
        double phiFit_dR_Kp_Ds;
        double phiFit_dR_Km_Ds;
        double phiFit_dR_phi_Ds;
        double phiFit_dR_pi_Ds;

        // Fit on Ds 
        double DsFit_CHI2;
        double DsFit_NDOF;
        double DsFit_CHI2NDOF;
        double DsFit_ENDVX_X;
        double DsFit_ENDVX_Y;
        double DsFit_ENDVX_Z;
        double DsFit_ENDVX_XERR;
        double DsFit_ENDVX_YERR;
        double DsFit_ENDVX_ZERR;

        double DsFit_Kp_ETA;
        double DsFit_Kp_PHI;
        double DsFit_Kp_P;
        double DsFit_Kp_PT;
        double DsFit_Kp_PX;
        double DsFit_Kp_PY;
        double DsFit_Kp_PZ;

        double DsFit_Km_ETA;
        double DsFit_Km_PHI;
        double DsFit_Km_P;
        double DsFit_Km_PT;
        double DsFit_Km_PX;
        double DsFit_Km_PY;
        double DsFit_Km_PZ;

        double DsFit_pi_ETA;
        double DsFit_pi_PHI;
        double DsFit_pi_P;
        double DsFit_pi_PT;
        double DsFit_pi_PX;
        double DsFit_pi_PY;
        double DsFit_pi_PZ;

        double DsFit_phi_ETA;
        double DsFit_phi_PHI;
        double DsFit_phi_P;
        double DsFit_phi_PT;
        double DsFit_phi_PX;
        double DsFit_phi_PY;
        double DsFit_phi_PZ;
        double DsFit_phi_M;

        double DsFit_Ds_ETA;
        double DsFit_Ds_PHI;
        double DsFit_Ds_P;
        double DsFit_Ds_PT;
        double DsFit_Ds_PX;
        double DsFit_Ds_PY;
        double DsFit_Ds_PZ;
        double DsFit_Ds_M;

        double DsFit_Kp_PP;
        double DsFit_Kp_PL;
        double DsFit_Km_PP;
        double DsFit_Km_PL;

        double DsFit_phi_PP;
        double DsFit_phi_PL;
        double DsFit_pi_PP;
        double DsFit_pi_PL;

        double DsFit_dR_Kp_Km;
        double DsFit_dR_Kp_phi;
        double DsFit_dR_Km_phi;
        double DsFit_dR_Kp_pi;
        double DsFit_dR_Km_pi;
        double DsFit_dR_pi_phi;
        double DsFit_dR_Kp_Ds;
        double DsFit_dR_Km_Ds;
        double DsFit_dR_phi_Ds;
        double DsFit_dR_pi_Ds;

        double DsFit_Mconstraint_Ds_M;

        double Ds_FDxy;
        double Ds_FDxy_Err;
        double Ds_FDxy_Chi2;
        double Ds_FDz;
        double Ds_FDz_Err;
        double Ds_FDz_Chi2;
        double Ds_FD;
        double Ds_FD_Err;
        double Ds_FD_Chi2;
        double Ds_DIRA_angle;
        double Ds_DIRA;
        double Kp_IP;
        double Kp_IP_Err;
        double Kp_IP_Chi2;
        double Km_IP;
        double Km_IP_Err;
        double Km_IP_Chi2;
        double pi_IP;
        double pi_IP_Err;
        double pi_IP_Chi2;

        std::vector<bool> PV_withBS_IsValid_vec;
        std::vector<bool> PV_withBS_IsFake_vec;
        std::vector<double> PV_withBS_CHI2_vec;
        std::vector<double> PV_withBS_NDOF_vec;
        std::vector<double> PV_withBS_CHI2NDOF_vec;
        std::vector<double> PV_withBS_X_vec;
        std::vector<double> PV_withBS_Y_vec;
        std::vector<double> PV_withBS_Z_vec;
        std::vector<double> PV_withBS_XERR_vec;
        std::vector<double> PV_withBS_YERR_vec;
        std::vector<double> PV_withBS_ZERR_vec;
        std::vector<bool> PV_noBS_IsValid_vec;
        std::vector<bool> PV_noBS_IsFake_vec;
        std::vector<double> PV_noBS_CHI2_vec;
        std::vector<double> PV_noBS_NDOF_vec;
        std::vector<double> PV_noBS_CHI2NDOF_vec;
        std::vector<double> PV_noBS_X_vec;
        std::vector<double> PV_noBS_Y_vec;
        std::vector<double> PV_noBS_Z_vec;
        std::vector<double> PV_noBS_XERR_vec;
        std::vector<double> PV_noBS_YERR_vec;
        std::vector<double> PV_noBS_ZERR_vec;

        std::vector<double> Kp_ETA_vec;
        std::vector<double> Kp_PHI_vec;
        std::vector<double> Kp_ORIVX_X_vec;
        std::vector<double> Kp_ORIVX_Y_vec;
        std::vector<double> Kp_ORIVX_Z_vec;
        std::vector<double> Kp_P_vec;
        std::vector<double> Kp_PT_vec;
        std::vector<double> Kp_PX_vec;
        std::vector<double> Kp_PY_vec;
        std::vector<double> Kp_PZ_vec;

        std::vector<double> Km_ETA_vec;
        std::vector<double> Km_PHI_vec;
        std::vector<double> Km_ORIVX_X_vec;
        std::vector<double> Km_ORIVX_Y_vec;
        std::vector<double> Km_ORIVX_Z_vec;
        std::vector<double> Km_P_vec;
        std::vector<double> Km_PT_vec;
        std::vector<double> Km_PX_vec;
        std::vector<double> Km_PY_vec;
        std::vector<double> Km_PZ_vec;

        std::vector<double> pi_ETA_vec;
        std::vector<double> pi_PHI_vec;
        std::vector<double> pi_ORIVX_X_vec;
        std::vector<double> pi_ORIVX_Y_vec;
        std::vector<double> pi_ORIVX_Z_vec;
        std::vector<double> pi_P_vec;
        std::vector<double> pi_PT_vec;
        std::vector<double> pi_PX_vec;
        std::vector<double> pi_PY_vec;
        std::vector<double> pi_PZ_vec;

        std::vector<double> phi_ETA_vec;
        std::vector<double> phi_PHI_vec;
        std::vector<double> phi_P_vec;
        std::vector<double> phi_PT_vec;
        std::vector<double> phi_PX_vec;
        std::vector<double> phi_PY_vec;
        std::vector<double> phi_PZ_vec;
        std::vector<double> phi_M_vec;

        std::vector<double> Ds_ETA_vec;
        std::vector<double> Ds_PHI_vec;
        std::vector<double> Ds_P_vec;
        std::vector<double> Ds_PT_vec;
        std::vector<double> Ds_PX_vec;
        std::vector<double> Ds_PY_vec;
        std::vector<double> Ds_PZ_vec;
        std::vector<double> Ds_M_vec;

        std::vector<double> Kp_PP_vec;
        std::vector<double> Kp_PL_vec;
        std::vector<double> Km_PP_vec;
        std::vector<double> Km_PL_vec;

        std::vector<double> phi_PP_vec;
        std::vector<double> phi_PL_vec;
        std::vector<double> pi_PP_vec;
        std::vector<double> pi_PL_vec;

        std::vector<double> dR_Kp_Km_vec;
        std::vector<double> dR_Kp_phi_vec;
        std::vector<double> dR_Km_phi_vec;
        std::vector<double> dR_Kp_pi_vec;
        std::vector<double> dR_Km_pi_vec;
        std::vector<double> dR_pi_phi_vec;
        std::vector<double> dR_Kp_Ds_vec;
        std::vector<double> dR_Km_Ds_vec;
        std::vector<double> dR_phi_Ds_vec;
        std::vector<double> dR_pi_Ds_vec;

        std::vector<double> dxy_Kp_Km_vec;
        std::vector<double> dxy_Kp_phi_vec;
        std::vector<double> dxy_Km_phi_vec;
        std::vector<double> dxy_Kp_pi_vec;
        std::vector<double> dxy_Km_pi_vec;
        std::vector<double> dxy_pi_phi_vec;
        std::vector<double> dxy_Kp_Ds_vec;
        std::vector<double> dxy_Km_Ds_vec;
        std::vector<double> dxy_phi_Ds_vec;
        std::vector<double> dxy_pi_Ds_vec;

        std::vector<double> dz_Kp_Km_vec;
        std::vector<double> dz_Kp_phi_vec;
        std::vector<double> dz_Km_phi_vec;
        std::vector<double> dz_Kp_pi_vec;
        std::vector<double> dz_Km_pi_vec;
        std::vector<double> dz_pi_phi_vec;
        std::vector<double> dz_Kp_Ds_vec;
        std::vector<double> dz_Km_Ds_vec;
        std::vector<double> dz_phi_Ds_vec;
        std::vector<double> dz_pi_Ds_vec;

        // Fit on phi
        std::vector<double> phiFit_CHI2_vec;
        std::vector<double> phiFit_NDOF_vec;
        std::vector<double> phiFit_CHI2NDOF_vec;
        std::vector<double> phiFit_ENDVX_X_vec;
        std::vector<double> phiFit_ENDVX_Y_vec;
        std::vector<double> phiFit_ENDVX_Z_vec;
        std::vector<double> phiFit_ENDVX_XERR_vec;
        std::vector<double> phiFit_ENDVX_YERR_vec;
        std::vector<double> phiFit_ENDVX_ZERR_vec;

        std::vector<double> phiFit_Kp_ETA_vec;
        std::vector<double> phiFit_Kp_PHI_vec;
        std::vector<double> phiFit_Kp_P_vec;
        std::vector<double> phiFit_Kp_PT_vec;
        std::vector<double> phiFit_Kp_PX_vec;
        std::vector<double> phiFit_Kp_PY_vec;
        std::vector<double> phiFit_Kp_PZ_vec;

        std::vector<double> phiFit_Km_ETA_vec;
        std::vector<double> phiFit_Km_PHI_vec;
        std::vector<double> phiFit_Km_P_vec;
        std::vector<double> phiFit_Km_PT_vec;
        std::vector<double> phiFit_Km_PX_vec;
        std::vector<double> phiFit_Km_PY_vec;
        std::vector<double> phiFit_Km_PZ_vec;

        std::vector<double> phiFit_pi_ETA_vec;
        std::vector<double> phiFit_pi_PHI_vec;
        std::vector<double> phiFit_pi_P_vec;
        std::vector<double> phiFit_pi_PT_vec;
        std::vector<double> phiFit_pi_PX_vec;
        std::vector<double> phiFit_pi_PY_vec;
        std::vector<double> phiFit_pi_PZ_vec;

        std::vector<double> phiFit_phi_ETA_vec;
        std::vector<double> phiFit_phi_PHI_vec;
        std::vector<double> phiFit_phi_P_vec;
        std::vector<double> phiFit_phi_PT_vec;
        std::vector<double> phiFit_phi_PX_vec;
        std::vector<double> phiFit_phi_PY_vec;
        std::vector<double> phiFit_phi_PZ_vec;
        std::vector<double> phiFit_phi_M_vec;

        std::vector<double> phiFit_Ds_ETA_vec;
        std::vector<double> phiFit_Ds_PHI_vec;
        std::vector<double> phiFit_Ds_P_vec;
        std::vector<double> phiFit_Ds_PT_vec;
        std::vector<double> phiFit_Ds_PX_vec;
        std::vector<double> phiFit_Ds_PY_vec;
        std::vector<double> phiFit_Ds_PZ_vec;
        std::vector<double> phiFit_Ds_M_vec;

        std::vector<double> phiFit_Kp_PP_vec;
        std::vector<double> phiFit_Kp_PL_vec;
        std::vector<double> phiFit_Km_PP_vec;
        std::vector<double> phiFit_Km_PL_vec;

        std::vector<double> phiFit_phi_PP_vec;
        std::vector<double> phiFit_phi_PL_vec;
        std::vector<double> phiFit_pi_PP_vec;
        std::vector<double> phiFit_pi_PL_vec;

        std::vector<double> phiFit_dR_Kp_Km_vec;
        std::vector<double> phiFit_dR_Kp_phi_vec;
        std::vector<double> phiFit_dR_Km_phi_vec;
        std::vector<double> phiFit_dR_Kp_pi_vec;
        std::vector<double> phiFit_dR_Km_pi_vec;
        std::vector<double> phiFit_dR_pi_phi_vec;
        std::vector<double> phiFit_dR_Kp_Ds_vec;
        std::vector<double> phiFit_dR_Km_Ds_vec;
        std::vector<double> phiFit_dR_phi_Ds_vec;
        std::vector<double> phiFit_dR_pi_Ds_vec;

        // Fit on Ds 
        std::vector<double> DsFit_CHI2_vec;
        std::vector<double> DsFit_NDOF_vec;
        std::vector<double> DsFit_CHI2NDOF_vec;
        std::vector<double> DsFit_ENDVX_X_vec;
        std::vector<double> DsFit_ENDVX_Y_vec;
        std::vector<double> DsFit_ENDVX_Z_vec;
        std::vector<double> DsFit_ENDVX_XERR_vec;
        std::vector<double> DsFit_ENDVX_YERR_vec;
        std::vector<double> DsFit_ENDVX_ZERR_vec;

        std::vector<double> DsFit_Kp_ETA_vec;
        std::vector<double> DsFit_Kp_PHI_vec;
        std::vector<double> DsFit_Kp_P_vec;
        std::vector<double> DsFit_Kp_PT_vec;
        std::vector<double> DsFit_Kp_PX_vec;
        std::vector<double> DsFit_Kp_PY_vec;
        std::vector<double> DsFit_Kp_PZ_vec;

        std::vector<double> DsFit_Km_ETA_vec;
        std::vector<double> DsFit_Km_PHI_vec;
        std::vector<double> DsFit_Km_P_vec;
        std::vector<double> DsFit_Km_PT_vec;
        std::vector<double> DsFit_Km_PX_vec;
        std::vector<double> DsFit_Km_PY_vec;
        std::vector<double> DsFit_Km_PZ_vec;

        std::vector<double> DsFit_pi_ETA_vec;
        std::vector<double> DsFit_pi_PHI_vec;
        std::vector<double> DsFit_pi_P_vec;
        std::vector<double> DsFit_pi_PT_vec;
        std::vector<double> DsFit_pi_PX_vec;
        std::vector<double> DsFit_pi_PY_vec;
        std::vector<double> DsFit_pi_PZ_vec;

        std::vector<double> DsFit_phi_ETA_vec;
        std::vector<double> DsFit_phi_PHI_vec;
        std::vector<double> DsFit_phi_P_vec;
        std::vector<double> DsFit_phi_PT_vec;
        std::vector<double> DsFit_phi_PX_vec;
        std::vector<double> DsFit_phi_PY_vec;
        std::vector<double> DsFit_phi_PZ_vec;
        std::vector<double> DsFit_phi_M_vec;

        std::vector<double> DsFit_Ds_ETA_vec;
        std::vector<double> DsFit_Ds_PHI_vec;
        std::vector<double> DsFit_Ds_P_vec;
        std::vector<double> DsFit_Ds_PT_vec;
        std::vector<double> DsFit_Ds_PX_vec;
        std::vector<double> DsFit_Ds_PY_vec;
        std::vector<double> DsFit_Ds_PZ_vec;
        std::vector<double> DsFit_Ds_M_vec;

        std::vector<double> DsFit_Kp_PP_vec;
        std::vector<double> DsFit_Kp_PL_vec;
        std::vector<double> DsFit_Km_PP_vec;
        std::vector<double> DsFit_Km_PL_vec;

        std::vector<double> DsFit_phi_PP_vec;
        std::vector<double> DsFit_phi_PL_vec;
        std::vector<double> DsFit_pi_PP_vec;
        std::vector<double> DsFit_pi_PL_vec;

        std::vector<double> DsFit_dR_Kp_Km_vec;
        std::vector<double> DsFit_dR_Kp_phi_vec;
        std::vector<double> DsFit_dR_Km_phi_vec;
        std::vector<double> DsFit_dR_Kp_pi_vec;
        std::vector<double> DsFit_dR_Km_pi_vec;
        std::vector<double> DsFit_dR_pi_phi_vec;
        std::vector<double> DsFit_dR_Kp_Ds_vec;
        std::vector<double> DsFit_dR_Km_Ds_vec;
        std::vector<double> DsFit_dR_phi_Ds_vec;
        std::vector<double> DsFit_dR_pi_Ds_vec;

        std::vector<double> DsFit_Mconstraint_Ds_M_vec;

        std::vector<double> Ds_FDxy_vec;
        std::vector<double> Ds_FDxy_Err_vec;
        std::vector<double> Ds_FDxy_Chi2_vec;
        std::vector<double> Ds_FDz_vec;
        std::vector<double> Ds_FDz_Err_vec;
        std::vector<double> Ds_FDz_Chi2_vec;
        std::vector<double> Ds_FD_vec;
        std::vector<double> Ds_FD_Err_vec;
        std::vector<double> Ds_FD_Chi2_vec;
        std::vector<double> Ds_DIRA_angle_vec;
        std::vector<double> Ds_DIRA_vec;
        std::vector<double> Kp_IP_vec;
        std::vector<double> Kp_IP_Err_vec;
        std::vector<double> Kp_IP_Chi2_vec;
        std::vector<double> Km_IP_vec;
        std::vector<double> Km_IP_Err_vec;
        std::vector<double> Km_IP_Chi2_vec;
        std::vector<double> pi_IP_vec;
        std::vector<double> pi_IP_Err_vec;
        std::vector<double> pi_IP_Chi2_vec;

        std::vector<double> best_Kp_ETA_vec;
        std::vector<double> best_Kp_PHI_vec;
        std::vector<double> best_Kp_ORIVX_X_vec;
        std::vector<double> best_Kp_ORIVX_Y_vec;
        std::vector<double> best_Kp_ORIVX_Z_vec;
        std::vector<double> best_Kp_P_vec;
        std::vector<double> best_Kp_PT_vec;
        std::vector<double> best_Kp_PX_vec;
        std::vector<double> best_Kp_PY_vec;
        std::vector<double> best_Kp_PZ_vec;

        std::vector<double> best_Km_ETA_vec;
        std::vector<double> best_Km_PHI_vec;
        std::vector<double> best_Km_ORIVX_X_vec;
        std::vector<double> best_Km_ORIVX_Y_vec;
        std::vector<double> best_Km_ORIVX_Z_vec;
        std::vector<double> best_Km_P_vec;
        std::vector<double> best_Km_PT_vec;
        std::vector<double> best_Km_PX_vec;
        std::vector<double> best_Km_PY_vec;
        std::vector<double> best_Km_PZ_vec;

        std::vector<double> best_pi_ETA_vec;
        std::vector<double> best_pi_PHI_vec;
        std::vector<double> best_pi_ORIVX_X_vec;
        std::vector<double> best_pi_ORIVX_Y_vec;
        std::vector<double> best_pi_ORIVX_Z_vec;
        std::vector<double> best_pi_P_vec;
        std::vector<double> best_pi_PT_vec;
        std::vector<double> best_pi_PX_vec;
        std::vector<double> best_pi_PY_vec;
        std::vector<double> best_pi_PZ_vec;

        std::vector<double> best_phi_ETA_vec;
        std::vector<double> best_phi_PHI_vec;
        std::vector<double> best_phi_P_vec;
        std::vector<double> best_phi_PT_vec;
        std::vector<double> best_phi_PX_vec;
        std::vector<double> best_phi_PY_vec;
        std::vector<double> best_phi_PZ_vec;
        std::vector<double> best_phi_M_vec;

        std::vector<double> best_Ds_ETA_vec;
        std::vector<double> best_Ds_PHI_vec;
        std::vector<double> best_Ds_P_vec;
        std::vector<double> best_Ds_PT_vec;
        std::vector<double> best_Ds_PX_vec;
        std::vector<double> best_Ds_PY_vec;
        std::vector<double> best_Ds_PZ_vec;
        std::vector<double> best_Ds_M_vec;

        std::vector<double> best_Kp_PP_vec;
        std::vector<double> best_Kp_PL_vec;
        std::vector<double> best_Km_PP_vec;
        std::vector<double> best_Km_PL_vec;

        std::vector<double> best_phi_PP_vec;
        std::vector<double> best_phi_PL_vec;
        std::vector<double> best_pi_PP_vec;
        std::vector<double> best_pi_PL_vec;

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

        std::vector<double> best_dxy_Kp_Km_vec;
        std::vector<double> best_dxy_Kp_phi_vec;
        std::vector<double> best_dxy_Km_phi_vec;
        std::vector<double> best_dxy_Kp_pi_vec;
        std::vector<double> best_dxy_Km_pi_vec;
        std::vector<double> best_dxy_pi_phi_vec;
        std::vector<double> best_dxy_Kp_Ds_vec;
        std::vector<double> best_dxy_Km_Ds_vec;
        std::vector<double> best_dxy_phi_Ds_vec;
        std::vector<double> best_dxy_pi_Ds_vec;

        std::vector<double> best_dz_Kp_Km_vec;
        std::vector<double> best_dz_Kp_phi_vec;
        std::vector<double> best_dz_Km_phi_vec;
        std::vector<double> best_dz_Kp_pi_vec;
        std::vector<double> best_dz_Km_pi_vec;
        std::vector<double> best_dz_pi_phi_vec;
        std::vector<double> best_dz_Kp_Ds_vec;
        std::vector<double> best_dz_Km_Ds_vec;
        std::vector<double> best_dz_phi_Ds_vec;
        std::vector<double> best_dz_pi_Ds_vec;

        // Fit on phi
        std::vector<double> best_phiFit_CHI2_vec;
        std::vector<double> best_phiFit_NDOF_vec;
        std::vector<double> best_phiFit_CHI2NDOF_vec;
        std::vector<double> best_phiFit_ENDVX_X_vec;
        std::vector<double> best_phiFit_ENDVX_Y_vec;
        std::vector<double> best_phiFit_ENDVX_Z_vec;
        std::vector<double> best_phiFit_ENDVX_XERR_vec;
        std::vector<double> best_phiFit_ENDVX_YERR_vec;
        std::vector<double> best_phiFit_ENDVX_ZERR_vec;

        std::vector<double> best_phiFit_Kp_ETA_vec;
        std::vector<double> best_phiFit_Kp_PHI_vec;
        std::vector<double> best_phiFit_Kp_P_vec;
        std::vector<double> best_phiFit_Kp_PT_vec;
        std::vector<double> best_phiFit_Kp_PX_vec;
        std::vector<double> best_phiFit_Kp_PY_vec;
        std::vector<double> best_phiFit_Kp_PZ_vec;

        std::vector<double> best_phiFit_Km_ETA_vec;
        std::vector<double> best_phiFit_Km_PHI_vec;
        std::vector<double> best_phiFit_Km_P_vec;
        std::vector<double> best_phiFit_Km_PT_vec;
        std::vector<double> best_phiFit_Km_PX_vec;
        std::vector<double> best_phiFit_Km_PY_vec;
        std::vector<double> best_phiFit_Km_PZ_vec;

        std::vector<double> best_phiFit_pi_ETA_vec;
        std::vector<double> best_phiFit_pi_PHI_vec;
        std::vector<double> best_phiFit_pi_P_vec;
        std::vector<double> best_phiFit_pi_PT_vec;
        std::vector<double> best_phiFit_pi_PX_vec;
        std::vector<double> best_phiFit_pi_PY_vec;
        std::vector<double> best_phiFit_pi_PZ_vec;

        std::vector<double> best_phiFit_phi_ETA_vec;
        std::vector<double> best_phiFit_phi_PHI_vec;
        std::vector<double> best_phiFit_phi_P_vec;
        std::vector<double> best_phiFit_phi_PT_vec;
        std::vector<double> best_phiFit_phi_PX_vec;
        std::vector<double> best_phiFit_phi_PY_vec;
        std::vector<double> best_phiFit_phi_PZ_vec;
        std::vector<double> best_phiFit_phi_M_vec;

        std::vector<double> best_phiFit_Ds_ETA_vec;
        std::vector<double> best_phiFit_Ds_PHI_vec;
        std::vector<double> best_phiFit_Ds_P_vec;
        std::vector<double> best_phiFit_Ds_PT_vec;
        std::vector<double> best_phiFit_Ds_PX_vec;
        std::vector<double> best_phiFit_Ds_PY_vec;
        std::vector<double> best_phiFit_Ds_PZ_vec;
        std::vector<double> best_phiFit_Ds_M_vec;

        std::vector<double> best_phiFit_Kp_PP_vec;
        std::vector<double> best_phiFit_Kp_PL_vec;
        std::vector<double> best_phiFit_Km_PP_vec;
        std::vector<double> best_phiFit_Km_PL_vec;

        std::vector<double> best_phiFit_phi_PP_vec;
        std::vector<double> best_phiFit_phi_PL_vec;
        std::vector<double> best_phiFit_pi_PP_vec;
        std::vector<double> best_phiFit_pi_PL_vec;

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

        // Fit on Ds 
        std::vector<double> best_DsFit_CHI2_vec;
        std::vector<double> best_DsFit_NDOF_vec;
        std::vector<double> best_DsFit_CHI2NDOF_vec;
        std::vector<double> best_DsFit_ENDVX_X_vec;
        std::vector<double> best_DsFit_ENDVX_Y_vec;
        std::vector<double> best_DsFit_ENDVX_Z_vec;
        std::vector<double> best_DsFit_ENDVX_XERR_vec;
        std::vector<double> best_DsFit_ENDVX_YERR_vec;
        std::vector<double> best_DsFit_ENDVX_ZERR_vec;

        std::vector<double> best_DsFit_Kp_ETA_vec;
        std::vector<double> best_DsFit_Kp_PHI_vec;
        std::vector<double> best_DsFit_Kp_P_vec;
        std::vector<double> best_DsFit_Kp_PT_vec;
        std::vector<double> best_DsFit_Kp_PX_vec;
        std::vector<double> best_DsFit_Kp_PY_vec;
        std::vector<double> best_DsFit_Kp_PZ_vec;

        std::vector<double> best_DsFit_Km_ETA_vec;
        std::vector<double> best_DsFit_Km_PHI_vec;
        std::vector<double> best_DsFit_Km_P_vec;
        std::vector<double> best_DsFit_Km_PT_vec;
        std::vector<double> best_DsFit_Km_PX_vec;
        std::vector<double> best_DsFit_Km_PY_vec;
        std::vector<double> best_DsFit_Km_PZ_vec;

        std::vector<double> best_DsFit_pi_ETA_vec;
        std::vector<double> best_DsFit_pi_PHI_vec;
        std::vector<double> best_DsFit_pi_P_vec;
        std::vector<double> best_DsFit_pi_PT_vec;
        std::vector<double> best_DsFit_pi_PX_vec;
        std::vector<double> best_DsFit_pi_PY_vec;
        std::vector<double> best_DsFit_pi_PZ_vec;

        std::vector<double> best_DsFit_phi_ETA_vec;
        std::vector<double> best_DsFit_phi_PHI_vec;
        std::vector<double> best_DsFit_phi_P_vec;
        std::vector<double> best_DsFit_phi_PT_vec;
        std::vector<double> best_DsFit_phi_PX_vec;
        std::vector<double> best_DsFit_phi_PY_vec;
        std::vector<double> best_DsFit_phi_PZ_vec;
        std::vector<double> best_DsFit_phi_M_vec;

        std::vector<double> best_DsFit_Ds_ETA_vec;
        std::vector<double> best_DsFit_Ds_PHI_vec;
        std::vector<double> best_DsFit_Ds_P_vec;
        std::vector<double> best_DsFit_Ds_PT_vec;
        std::vector<double> best_DsFit_Ds_PX_vec;
        std::vector<double> best_DsFit_Ds_PY_vec;
        std::vector<double> best_DsFit_Ds_PZ_vec;
        std::vector<double> best_DsFit_Ds_M_vec;

        std::vector<double> best_DsFit_Kp_PP_vec;
        std::vector<double> best_DsFit_Kp_PL_vec;
        std::vector<double> best_DsFit_Km_PP_vec;
        std::vector<double> best_DsFit_Km_PL_vec;

        std::vector<double> best_DsFit_phi_PP_vec;
        std::vector<double> best_DsFit_phi_PL_vec;
        std::vector<double> best_DsFit_pi_PP_vec;
        std::vector<double> best_DsFit_pi_PL_vec;

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

        std::vector<double> best_DsFit_Mconstraint_Ds_M_vec;

        std::vector<double> best_Ds_FDxy_vec;
        std::vector<double> best_Ds_FDxy_Err_vec;
        std::vector<double> best_Ds_FDxy_Chi2_vec;
        std::vector<double> best_Ds_FDz_vec;
        std::vector<double> best_Ds_FDz_Err_vec;
        std::vector<double> best_Ds_FDz_Chi2_vec;
        std::vector<double> best_Ds_FD_vec;
        std::vector<double> best_Ds_FD_Err_vec;
        std::vector<double> best_Ds_FD_Chi2_vec;
        std::vector<double> best_Ds_DIRA_angle_vec;
        std::vector<double> best_Ds_DIRA_vec;
        std::vector<double> best_Kp_IP_vec;
        std::vector<double> best_Kp_IP_Err_vec;
        std::vector<double> best_Kp_IP_Chi2_vec;
        std::vector<double> best_Km_IP_vec;
        std::vector<double> best_Km_IP_Err_vec;
        std::vector<double> best_Km_IP_Chi2_vec;
        std::vector<double> best_pi_IP_vec;
        std::vector<double> best_pi_IP_Err_vec;
        std::vector<double> best_pi_IP_Chi2_vec;

};

#endif

#ifndef PFRECOTREE_H 
#define PFRECOTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class PFRecoTree 
{
    public:

        PFRecoTree(TTree *tree_);

        TTree *tree;

        void Init();
        void Kp_Reset();
        void Km_Reset();
        void phi_Reset();
        void pi_Reset();
        void Ds_Reset();
        void Fill_Vector();
        void CreateBranches();

        // Matched PF Candidates
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
        double Kp_PP;
        double Kp_PL;

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
        double Km_PP;
        double Km_PL;

        double phi_ETA;
        double phi_PHI;
        double phi_P;
        double phi_PT;
        double phi_PX;
        double phi_PY;
        double phi_PZ;
        double phi_PP;
        double phi_PL;
        double phi_CHI2;
        double phi_M;
        double phi_ENDVX_X;
        double phi_ENDVX_Y;
        double phi_ENDVX_Z;

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
        double pi_PP;
        double pi_PL;

        double Ds_ETA;
        double Ds_PHI;
        double Ds_P;
        double Ds_PT;
        double Ds_PX;
        double Ds_PY;
        double Ds_PZ;
        double Ds_CHI2;
        double Ds_M;
        double Ds_ENDVX_X;
        double Ds_ENDVX_Y;
        double Ds_ENDVX_Z;

        double dR_Kp_Km;
        double dR_Kp_phi;
        double dR_Km_phi;
        double dR_Kp_pi;
        double dR_Km_pi;
        double dR_phi_pi;
        double dR_Kp_Ds;
        double dR_Km_Ds;
        double dR_phi_Ds;
        double dR_pi_Ds;

        double d_Kp_Km;
        double d_Kp_phi;
        double d_Km_phi;
        double d_Kp_pi;
        double d_Km_pi;
        double d_phi_pi;
        double d_Kp_Ds;
        double d_Km_Ds;
        double d_phi_Ds;
        double d_pi_Ds;

        double phiFit_Kp_ETA;
        double phiFit_Kp_PHI;
        double phiFit_Kp_P;
        double phiFit_Kp_PT;
        double phiFit_Kp_PX;
        double phiFit_Kp_PY;
        double phiFit_Kp_PZ;
        double phiFit_Kp_PP;
        double phiFit_Kp_PL;

        double phiFit_Km_ETA;
        double phiFit_Km_PHI;
        double phiFit_Km_P;
        double phiFit_Km_PT;
        double phiFit_Km_PX;
        double phiFit_Km_PY;
        double phiFit_Km_PZ;
        double phiFit_Km_PP;
        double phiFit_Km_PL;

        double DsFit_Kp_ETA;
        double DsFit_Kp_PHI;
        double DsFit_Kp_P;
        double DsFit_Kp_PT;
        double DsFit_Kp_PX;
        double DsFit_Kp_PY;
        double DsFit_Kp_PZ;
        double DsFit_Kp_PP;
        double DsFit_Kp_PL;

        double DsFit_Km_ETA;
        double DsFit_Km_PHI;
        double DsFit_Km_P;
        double DsFit_Km_PT;
        double DsFit_Km_PX;
        double DsFit_Km_PY;
        double DsFit_Km_PZ;
        double DsFit_Km_PP;
        double DsFit_Km_PL;

        double DsFit_phi_ETA;
        double DsFit_phi_PHI;
        double DsFit_phi_P;
        double DsFit_phi_PT;
        double DsFit_phi_PX;
        double DsFit_phi_PY;
        double DsFit_phi_PZ;
        double DsFit_phi_PP;
        double DsFit_phi_PL;
        double DsFit_phi_M;

        double DsFit_pi_ETA;
        double DsFit_pi_PHI;
        double DsFit_pi_P;
        double DsFit_pi_PT;
        double DsFit_pi_PX;
        double DsFit_pi_PY;
        double DsFit_pi_PZ;
        double DsFit_pi_PP;
        double DsFit_pi_PL;

        double DsFit_Mconstraint_Ds_M;
        
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
        std::vector<double> Kp_PP_vec;
        std::vector<double> Kp_PL_vec;

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
        std::vector<double> Km_PP_vec;
        std::vector<double> Km_PL_vec;

        std::vector<double> phi_ETA_vec;
        std::vector<double> phi_PHI_vec;
        std::vector<double> phi_P_vec;
        std::vector<double> phi_PT_vec;
        std::vector<double> phi_PX_vec;
        std::vector<double> phi_PY_vec;
        std::vector<double> phi_PZ_vec;
        std::vector<double> phi_PP_vec;
        std::vector<double> phi_PL_vec;
        std::vector<double> phi_CHI2_vec;
        std::vector<double> phi_M_vec;
        std::vector<double> phi_ENDVX_X_vec;
        std::vector<double> phi_ENDVX_Y_vec;
        std::vector<double> phi_ENDVX_Z_vec;

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
        std::vector<double> pi_PP_vec;
        std::vector<double> pi_PL_vec;

        std::vector<double> Ds_ETA_vec;
        std::vector<double> Ds_PHI_vec;
        std::vector<double> Ds_P_vec;
        std::vector<double> Ds_PT_vec;
        std::vector<double> Ds_PX_vec;
        std::vector<double> Ds_PY_vec;
        std::vector<double> Ds_PZ_vec;
        std::vector<double> Ds_CHI2_vec;
        std::vector<double> Ds_M_vec;
        std::vector<double> Ds_ENDVX_X_vec;
        std::vector<double> Ds_ENDVX_Y_vec;
        std::vector<double> Ds_ENDVX_Z_vec;

        std::vector<double> dR_Kp_Km_vec;
        std::vector<double> dR_Kp_phi_vec;
        std::vector<double> dR_Km_phi_vec;
        std::vector<double> dR_Kp_pi_vec;
        std::vector<double> dR_Km_pi_vec;
        std::vector<double> dR_phi_pi_vec;
        std::vector<double> dR_Kp_Ds_vec;
        std::vector<double> dR_Km_Ds_vec;
        std::vector<double> dR_phi_Ds_vec;
        std::vector<double> dR_pi_Ds_vec;

        std::vector<double> d_Kp_Km_vec;
        std::vector<double> d_Kp_phi_vec;
        std::vector<double> d_Km_phi_vec;
        std::vector<double> d_Kp_pi_vec;
        std::vector<double> d_Km_pi_vec;
        std::vector<double> d_phi_pi_vec;
        std::vector<double> d_Kp_Ds_vec;
        std::vector<double> d_Km_Ds_vec;
        std::vector<double> d_phi_Ds_vec;
        std::vector<double> d_pi_Ds_vec;

        std::vector<double> phiFit_Kp_ETA_vec;
        std::vector<double> phiFit_Kp_PHI_vec;
        std::vector<double> phiFit_Kp_P_vec;
        std::vector<double> phiFit_Kp_PT_vec;
        std::vector<double> phiFit_Kp_PX_vec;
        std::vector<double> phiFit_Kp_PY_vec;
        std::vector<double> phiFit_Kp_PZ_vec;
        std::vector<double> phiFit_Kp_PP_vec;
        std::vector<double> phiFit_Kp_PL_vec;

        std::vector<double> phiFit_Km_ETA_vec;
        std::vector<double> phiFit_Km_PHI_vec;
        std::vector<double> phiFit_Km_P_vec;
        std::vector<double> phiFit_Km_PT_vec;
        std::vector<double> phiFit_Km_PX_vec;
        std::vector<double> phiFit_Km_PY_vec;
        std::vector<double> phiFit_Km_PZ_vec;
        std::vector<double> phiFit_Km_PP_vec;
        std::vector<double> phiFit_Km_PL_vec;

        std::vector<double> DsFit_Kp_ETA_vec;
        std::vector<double> DsFit_Kp_PHI_vec;
        std::vector<double> DsFit_Kp_P_vec;
        std::vector<double> DsFit_Kp_PT_vec;
        std::vector<double> DsFit_Kp_PX_vec;
        std::vector<double> DsFit_Kp_PY_vec;
        std::vector<double> DsFit_Kp_PZ_vec;
        std::vector<double> DsFit_Kp_PP_vec;
        std::vector<double> DsFit_Kp_PL_vec;

        std::vector<double> DsFit_Km_ETA_vec;
        std::vector<double> DsFit_Km_PHI_vec;
        std::vector<double> DsFit_Km_P_vec;
        std::vector<double> DsFit_Km_PT_vec;
        std::vector<double> DsFit_Km_PX_vec;
        std::vector<double> DsFit_Km_PY_vec;
        std::vector<double> DsFit_Km_PZ_vec;
        std::vector<double> DsFit_Km_PP_vec;
        std::vector<double> DsFit_Km_PL_vec;

        std::vector<double> DsFit_phi_ETA_vec;
        std::vector<double> DsFit_phi_PHI_vec;
        std::vector<double> DsFit_phi_P_vec;
        std::vector<double> DsFit_phi_PT_vec;
        std::vector<double> DsFit_phi_PX_vec;
        std::vector<double> DsFit_phi_PY_vec;
        std::vector<double> DsFit_phi_PZ_vec;
        std::vector<double> DsFit_phi_PP_vec;
        std::vector<double> DsFit_phi_PL_vec;
        std::vector<double> DsFit_phi_M_vec;

        std::vector<double> DsFit_pi_ETA_vec;
        std::vector<double> DsFit_pi_PHI_vec;
        std::vector<double> DsFit_pi_P_vec;
        std::vector<double> DsFit_pi_PT_vec;
        std::vector<double> DsFit_pi_PX_vec;
        std::vector<double> DsFit_pi_PY_vec;
        std::vector<double> DsFit_pi_PZ_vec;
        std::vector<double> DsFit_pi_PP_vec;
        std::vector<double> DsFit_pi_PL_vec;

        std::vector<double> DsFit_Mconstraint_Ds_M_vec;
};

#endif

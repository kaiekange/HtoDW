#ifndef PVTREE_H 
#define PVTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class PVTree 
{
    public:

        PVTree(TTree *tree_);

        TTree *tree;

        void Init();
        void Gen_Reset();
        void Match_Reset();
        void Kp_Reset();
        void Km_Reset();
        void pi_Reset();
        void Gen_Fill_Vector();
        void Match_Fill_Vector();
        void Match_Fill_PV();
        void Fill_PV();
        void Best_Fill_Vector(int idxmax);
        void Fill_Vector();
        void CreateBranches();

        // Gen particles

        int num_Gen_Kp;
        int num_Gen_Km;
        int num_Gen_pi;
        int num_Gen_phi;
        int num_Gen_Ds;
        int num_Gen_mu;
        int num_Gen_nu;
        int num_Gen_W;
        int num_Gen_H;

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

        double Gen_mu_ETA;
        double Gen_mu_PHI;
        double Gen_mu_ORIVX_X;
        double Gen_mu_ORIVX_Y;
        double Gen_mu_ORIVX_Z;
        double Gen_mu_P;
        double Gen_mu_PT;
        double Gen_mu_PX;
        double Gen_mu_PY;
        double Gen_mu_PZ;

        double Gen_nu_ETA;
        double Gen_nu_PHI;
        double Gen_nu_ORIVX_X;
        double Gen_nu_ORIVX_Y;
        double Gen_nu_ORIVX_Z;
        double Gen_nu_P;
        double Gen_nu_PT;
        double Gen_nu_PX;
        double Gen_nu_PY;
        double Gen_nu_PZ;

        double Gen_W_ETA;
        double Gen_W_PHI;
        double Gen_W_ORIVX_X;
        double Gen_W_ORIVX_Y;
        double Gen_W_ORIVX_Z;
        double Gen_W_P;
        double Gen_W_PT;
        double Gen_W_PX;
        double Gen_W_PY;
        double Gen_W_PZ;
        
        double Gen_H_ETA;
        double Gen_H_PHI;
        double Gen_H_ORIVX_X;
        double Gen_H_ORIVX_Y;
        double Gen_H_ORIVX_Z;
        double Gen_H_P;
        double Gen_H_PT;
        double Gen_H_PX;
        double Gen_H_PY;
        double Gen_H_PZ;
        
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
        double match_phiFit_ENDVX_XERR;
        double match_phiFit_ENDVX_YERR;
        double match_phiFit_ENDVX_ZERR;

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
        double match_DsFit_ENDVX_XERR;
        double match_DsFit_ENDVX_YERR;
        double match_DsFit_ENDVX_ZERR;

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

        double match_PV_CHI2;
        double match_PV_NDOF;
        double match_PV_CHI2NDOF;
        double match_PV_X;
        double match_PV_Y;
        double match_PV_Z;
        double match_PV_XERR;
        double match_PV_YERR;
        double match_PV_ZERR;

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

        bool Kp_match;
        bool Km_match;
        bool pi_match;
        bool match_entry;
        bool non_match_entry;

        double PV_CHI2;
        double PV_NDOF;
        double PV_CHI2NDOF;
        double PV_X;
        double PV_Y;
        double PV_Z;
        double PV_XERR;
        double PV_YERR;
        double PV_ZERR;
        
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

        std::vector<double> Gen_mu_ETA_vec;
        std::vector<double> Gen_mu_PHI_vec;
        std::vector<double> Gen_mu_ORIVX_X_vec;
        std::vector<double> Gen_mu_ORIVX_Y_vec;
        std::vector<double> Gen_mu_ORIVX_Z_vec;
        std::vector<double> Gen_mu_P_vec;
        std::vector<double> Gen_mu_PT_vec;
        std::vector<double> Gen_mu_PX_vec;
        std::vector<double> Gen_mu_PY_vec;
        std::vector<double> Gen_mu_PZ_vec;

        std::vector<double> Gen_nu_ETA_vec;
        std::vector<double> Gen_nu_PHI_vec;
        std::vector<double> Gen_nu_ORIVX_X_vec;
        std::vector<double> Gen_nu_ORIVX_Y_vec;
        std::vector<double> Gen_nu_ORIVX_Z_vec;
        std::vector<double> Gen_nu_P_vec;
        std::vector<double> Gen_nu_PT_vec;
        std::vector<double> Gen_nu_PX_vec;
        std::vector<double> Gen_nu_PY_vec;
        std::vector<double> Gen_nu_PZ_vec;

        std::vector<double> Gen_W_ETA_vec;
        std::vector<double> Gen_W_PHI_vec;
        std::vector<double> Gen_W_ORIVX_X_vec;
        std::vector<double> Gen_W_ORIVX_Y_vec;
        std::vector<double> Gen_W_ORIVX_Z_vec;
        std::vector<double> Gen_W_P_vec;
        std::vector<double> Gen_W_PT_vec;
        std::vector<double> Gen_W_PX_vec;
        std::vector<double> Gen_W_PY_vec;
        std::vector<double> Gen_W_PZ_vec;

        std::vector<double> Gen_H_ETA_vec;
        std::vector<double> Gen_H_PHI_vec;
        std::vector<double> Gen_H_ORIVX_X_vec;
        std::vector<double> Gen_H_ORIVX_Y_vec;
        std::vector<double> Gen_H_ORIVX_Z_vec;
        std::vector<double> Gen_H_P_vec;
        std::vector<double> Gen_H_PT_vec;
        std::vector<double> Gen_H_PX_vec;
        std::vector<double> Gen_H_PY_vec;
        std::vector<double> Gen_H_PZ_vec;

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
        std::vector<double> match_phiFit_ENDVX_XERR_vec;
        std::vector<double> match_phiFit_ENDVX_YERR_vec;
        std::vector<double> match_phiFit_ENDVX_ZERR_vec;

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
        std::vector<double> match_DsFit_ENDVX_XERR_vec;
        std::vector<double> match_DsFit_ENDVX_YERR_vec;
        std::vector<double> match_DsFit_ENDVX_ZERR_vec;

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

        std::vector<double> match_PV_CHI2_vec;
        std::vector<double> match_PV_NDOF_vec;
        std::vector<double> match_PV_CHI2NDOF_vec;
        std::vector<double> match_PV_X_vec;
        std::vector<double> match_PV_Y_vec;
        std::vector<double> match_PV_Z_vec;
        std::vector<double> match_PV_XERR_vec;
        std::vector<double> match_PV_YERR_vec;
        std::vector<double> match_PV_ZERR_vec;

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

        std::vector<bool> Kp_match_vec;
        std::vector<bool> Km_match_vec;
        std::vector<bool> pi_match_vec;
        std::vector<bool> match_entry_vec;
        std::vector<bool> non_match_entry_vec;
        
        std::vector<double> PV_CHI2_vec;
        std::vector<double> PV_NDOF_vec;
        std::vector<double> PV_CHI2NDOF_vec;
        std::vector<double> PV_X_vec;
        std::vector<double> PV_Y_vec;
        std::vector<double> PV_Z_vec;
        std::vector<double> PV_XERR_vec;
        std::vector<double> PV_YERR_vec;
        std::vector<double> PV_ZERR_vec;

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

        std::vector<bool> best_match_entry_vec;
    
        std::vector<double> best_PV_CHI2_vec;
        std::vector<double> best_PV_NDOF_vec;
        std::vector<double> best_PV_CHI2NDOF_vec;
        std::vector<double> best_PV_X_vec;
        std::vector<double> best_PV_Y_vec;
        std::vector<double> best_PV_Z_vec;
        std::vector<double> best_PV_XERR_vec;
        std::vector<double> best_PV_YERR_vec;
        std::vector<double> best_PV_ZERR_vec;

};

#endif

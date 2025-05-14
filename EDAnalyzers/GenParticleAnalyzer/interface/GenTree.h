#ifndef GENTREE_H
#define GENTREE_H

#include <TTree.h>
#include <vector>

class GenTree
{
    public:

        GenTree(TTree *tree_);

        TTree *tree;

        void Init();
        void CreateBranches();

        // Gen Particles
        int num_Gen_Kp;
        int num_Gen_Km;
        int num_Gen_pi;
        int num_Gen_phi;
        int num_Gen_Ds;

        std::vector<double> Gen_Kp_ETA;
        std::vector<double> Gen_Kp_PHI;
        std::vector<double> Gen_Kp_ORIVX_X;
        std::vector<double> Gen_Kp_ORIVX_Y;
        std::vector<double> Gen_Kp_ORIVX_Z;
        std::vector<double> Gen_Kp_P;
        std::vector<double> Gen_Kp_PT;
        std::vector<double> Gen_Kp_PX;
        std::vector<double> Gen_Kp_PY;
        std::vector<double> Gen_Kp_PZ;
        std::vector<double> Gen_Kp_PP;
        std::vector<double> Gen_Kp_PL;

        std::vector<double> Gen_Km_ETA;
        std::vector<double> Gen_Km_PHI;
        std::vector<double> Gen_Km_ORIVX_X;
        std::vector<double> Gen_Km_ORIVX_Y;
        std::vector<double> Gen_Km_ORIVX_Z;
        std::vector<double> Gen_Km_P;
        std::vector<double> Gen_Km_PT;
        std::vector<double> Gen_Km_PX;
        std::vector<double> Gen_Km_PY;
        std::vector<double> Gen_Km_PZ;
        std::vector<double> Gen_Km_PP;
        std::vector<double> Gen_Km_PL;

        std::vector<double> Gen_pi_ETA;
        std::vector<double> Gen_pi_PHI;
        std::vector<double> Gen_pi_ORIVX_X;
        std::vector<double> Gen_pi_ORIVX_Y;
        std::vector<double> Gen_pi_ORIVX_Z;
        std::vector<double> Gen_pi_P;
        std::vector<double> Gen_pi_PT;
        std::vector<double> Gen_pi_PX;
        std::vector<double> Gen_pi_PY;
        std::vector<double> Gen_pi_PZ;
        std::vector<double> Gen_pi_PP;
        std::vector<double> Gen_pi_PL;

        std::vector<double> Gen_phi_ETA;
        std::vector<double> Gen_phi_PHI;
        std::vector<double> Gen_phi_ORIVX_X;
        std::vector<double> Gen_phi_ORIVX_Y;
        std::vector<double> Gen_phi_ORIVX_Z;
        std::vector<double> Gen_phi_P;
        std::vector<double> Gen_phi_PT;
        std::vector<double> Gen_phi_PX;
        std::vector<double> Gen_phi_PY;
        std::vector<double> Gen_phi_PZ;
        std::vector<double> Gen_phi_PP;
        std::vector<double> Gen_phi_PL;

        std::vector<double> Gen_Ds_ETA;
        std::vector<double> Gen_Ds_PHI;
        std::vector<double> Gen_Ds_ORIVX_X;
        std::vector<double> Gen_Ds_ORIVX_Y;
        std::vector<double> Gen_Ds_ORIVX_Z;
        std::vector<double> Gen_Ds_P;
        std::vector<double> Gen_Ds_PT;
        std::vector<double> Gen_Ds_PX;
        std::vector<double> Gen_Ds_PY;
        std::vector<double> Gen_Ds_PZ;

        std::vector<double> Gen_dR_Kp_Km;
        std::vector<double> Gen_dR_Kp_phi;
        std::vector<double> Gen_dR_Km_phi;
        std::vector<double> Gen_dR_Kp_pi;
        std::vector<double> Gen_dR_Km_pi;
        std::vector<double> Gen_dR_phi_pi;
        std::vector<double> Gen_dR_Kp_Ds;
        std::vector<double> Gen_dR_Km_Ds;
        std::vector<double> Gen_dR_phi_Ds;
        std::vector<double> Gen_dR_pi_Ds;

        std::vector<double> Gen_d_Kp_Km;
        std::vector<double> Gen_d_Kp_pi;
        std::vector<double> Gen_d_Km_pi;

        // Matched PF Candidates
        int num_match_Kp;
        int num_match_Km;
        int num_match_pi;
        
        std::vector<double> Kp_ETA;
        std::vector<double> Kp_PHI;
        std::vector<double> Kp_ORIVX_X;
        std::vector<double> Kp_ORIVX_Y;
        std::vector<double> Kp_ORIVX_Z;
        std::vector<double> Kp_P;
        std::vector<double> Kp_PT;
        std::vector<double> Kp_PX;
        std::vector<double> Kp_PY;
        std::vector<double> Kp_PZ;
        std::vector<double> Kp_PP;
        std::vector<double> Kp_PL;

        std::vector<double> Km_ETA;
        std::vector<double> Km_PHI;
        std::vector<double> Km_ORIVX_X;
        std::vector<double> Km_ORIVX_Y;
        std::vector<double> Km_ORIVX_Z;
        std::vector<double> Km_P;
        std::vector<double> Km_PT;
        std::vector<double> Km_PX;
        std::vector<double> Km_PY;
        std::vector<double> Km_PZ;
        std::vector<double> Km_PP;
        std::vector<double> Km_PL;

        std::vector<double> phi_ETA;
        std::vector<double> phi_PHI;
        std::vector<double> phi_P;
        std::vector<double> phi_PT;
        std::vector<double> phi_PX;
        std::vector<double> phi_PY;
        std::vector<double> phi_PZ;
        std::vector<double> phi_PP;
        std::vector<double> phi_PL;
        std::vector<double> phi_CHI2;
        std::vector<double> phi_M;
        std::vector<double> phi_ENDVX_X;
        std::vector<double> phi_ENDVX_Y;
        std::vector<double> phi_ENDVX_Z;

        std::vector<double> pi_ETA;
        std::vector<double> pi_PHI;
        std::vector<double> pi_ORIVX_X;
        std::vector<double> pi_ORIVX_Y;
        std::vector<double> pi_ORIVX_Z;
        std::vector<double> pi_P;
        std::vector<double> pi_PT;
        std::vector<double> pi_PX;
        std::vector<double> pi_PY;
        std::vector<double> pi_PZ;
        std::vector<double> pi_PP;
        std::vector<double> pi_PL;

        std::vector<double> Ds_ETA;
        std::vector<double> Ds_PHI;
        std::vector<double> Ds_P;
        std::vector<double> Ds_PT;
        std::vector<double> Ds_PX;
        std::vector<double> Ds_PY;
        std::vector<double> Ds_PZ;
        std::vector<double> Ds_CHI2;
        std::vector<double> Ds_M;
        std::vector<double> Ds_ENDVX_X;
        std::vector<double> Ds_ENDVX_Y;
        std::vector<double> Ds_ENDVX_Z;

        std::vector<double> dR_Kp_Km;
        std::vector<double> dR_Kp_phi;
        std::vector<double> dR_Km_phi;
        std::vector<double> dR_Kp_pi;
        std::vector<double> dR_Km_pi;
        std::vector<double> dR_phi_pi;
        std::vector<double> dR_Kp_Ds;
        std::vector<double> dR_Km_Ds;
        std::vector<double> dR_phi_Ds;
        std::vector<double> dR_pi_Ds;

        std::vector<double> d_Kp_Km;
        std::vector<double> d_Kp_phi;
        std::vector<double> d_Km_phi;
        std::vector<double> d_Kp_pi;
        std::vector<double> d_Km_pi;
        std::vector<double> d_phi_pi;
        std::vector<double> d_Kp_Ds;
        std::vector<double> d_Km_Ds;
        std::vector<double> d_phi_Ds;
        std::vector<double> d_pi_Ds;

        std::vector<double> phiFit_Kp_ETA;
        std::vector<double> phiFit_Kp_PHI;
        std::vector<double> phiFit_Kp_P;
        std::vector<double> phiFit_Kp_PT;
        std::vector<double> phiFit_Kp_PX;
        std::vector<double> phiFit_Kp_PY;
        std::vector<double> phiFit_Kp_PZ;
        std::vector<double> phiFit_Kp_PP;
        std::vector<double> phiFit_Kp_PL;

        std::vector<double> phiFit_Km_ETA;
        std::vector<double> phiFit_Km_PHI;
        std::vector<double> phiFit_Km_P;
        std::vector<double> phiFit_Km_PT;
        std::vector<double> phiFit_Km_PX;
        std::vector<double> phiFit_Km_PY;
        std::vector<double> phiFit_Km_PZ;
        std::vector<double> phiFit_Km_PP;
        std::vector<double> phiFit_Km_PL;

        std::vector<double> DsFit_Kp_ETA;
        std::vector<double> DsFit_Kp_PHI;
        std::vector<double> DsFit_Kp_P;
        std::vector<double> DsFit_Kp_PT;
        std::vector<double> DsFit_Kp_PX;
        std::vector<double> DsFit_Kp_PY;
        std::vector<double> DsFit_Kp_PZ;
        std::vector<double> DsFit_Kp_PP;
        std::vector<double> DsFit_Kp_PL;

        std::vector<double> DsFit_Km_ETA;
        std::vector<double> DsFit_Km_PHI;
        std::vector<double> DsFit_Km_P;
        std::vector<double> DsFit_Km_PT;
        std::vector<double> DsFit_Km_PX;
        std::vector<double> DsFit_Km_PY;
        std::vector<double> DsFit_Km_PZ;
        std::vector<double> DsFit_Km_PP;
        std::vector<double> DsFit_Km_PL;

        std::vector<double> DsFit_phi_ETA;
        std::vector<double> DsFit_phi_PHI;
        std::vector<double> DsFit_phi_P;
        std::vector<double> DsFit_phi_PT;
        std::vector<double> DsFit_phi_PX;
        std::vector<double> DsFit_phi_PY;
        std::vector<double> DsFit_phi_PZ;
        std::vector<double> DsFit_phi_PP;
        std::vector<double> DsFit_phi_PL;
        std::vector<double> DsFit_phi_M;

        std::vector<double> DsFit_pi_ETA;
        std::vector<double> DsFit_pi_PHI;
        std::vector<double> DsFit_pi_P;
        std::vector<double> DsFit_pi_PT;
        std::vector<double> DsFit_pi_PX;
        std::vector<double> DsFit_pi_PY;
        std::vector<double> DsFit_pi_PZ;
        std::vector<double> DsFit_pi_PP;
        std::vector<double> DsFit_pi_PL;

        std::vector<double> DsFit_Mconstraint_Ds_M;
};

#endif

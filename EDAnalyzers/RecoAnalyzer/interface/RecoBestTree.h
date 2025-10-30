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
        void phi_Reset();
        void pi_Reset();
        void Ds_Reset();
        /* void BS_Reset(); */
        /* void PV_Reset(); */
        void Best_Fill_Vector(int idxmax);
        void Fill_Vector();
        /* void PV_withBS_Fill_Vector(); */ 
        /* void PV_noBS_Fill_Vector(); */ 
        void CreateBranches();

        // Gen particles

        /* int BS_type; */
        /* double BS_X0; */
        /* double BS_Y0; */
        /* double BS_Z0; */
        /* double BS_SigmaZ; */
        /* double BS_dXdZ; */
        /* double BS_dYdZ; */
        /* double BS_BWX; */
        /* double BS_BWY; */
        /* double BS_X0err; */
        /* double BS_Y0err; */
        /* double BS_Z0err; */
        /* double BS_SigmaZ0err; */
        /* double BS_dXdZerr; */
        /* double BS_dYdZerr; */
        /* double BS_BWXerr; */
        /* double BS_BWYerr; */
        /* double BS_EmitX; */
        /* double BS_EmitY; */
        /* double BS_BetaStar; */

        /* bool PV_withBS_IsValid; */
        /* bool PV_withBS_IsFake; */
        /* double PV_withBS_chi2; */
        /* double PV_withBS_ndof; */
        /* double PV_withBS_chi2ndof; */
        /* double PV_withBS_X; */
        /* double PV_withBS_Y; */
        /* double PV_withBS_Z; */
        /* double PV_withBS_Xerr; */
        /* double PV_withBS_Yerr; */
        /* double PV_withBS_Zerr; */
        /* bool PV_noBS_IsValid; */
        /* bool PV_noBS_IsFake; */
        /* double PV_noBS_chi2; */
        /* double PV_noBS_ndof; */
        /* double PV_noBS_chi2ndof; */
        /* double PV_noBS_X; */
        /* double PV_noBS_Y; */
        /* double PV_noBS_Z; */
        /* double PV_noBS_Xerr; */
        /* double PV_noBS_Yerr; */
        /* double PV_noBS_Zerr; */

        int num_reco_phi;
        int num_reco_Ds;

        bool Kp_isIsolatedChargedHadron;
        int Kp_charge;
        double Kp_eta;
        double Kp_phi;
        double Kp_vx;
        double Kp_vy;
        double Kp_vz;
        double Kp_p;
        double Kp_pt;
        double Kp_px;
        double Kp_py;
        double Kp_pz;

        bool Km_isIsolatedChargedHadron;
        int Km_charge;
        double Km_eta;
        double Km_phi;
        double Km_vx;
        double Km_vy;
        double Km_vz;
        double Km_p;
        double Km_pt;
        double Km_px;
        double Km_py;
        double Km_pz;

        bool pi_isIsolatedChargedHadron;
        int pi_charge;
        double pi_eta;
        double pi_phi;
        double pi_vx;
        double pi_vy;
        double pi_vz;
        double pi_p;
        double pi_pt;
        double pi_px;
        double pi_py;
        double pi_pz;

        double phi_eta;
        double phi_phi;
        double phi_p;
        double phi_pt;
        double phi_px;
        double phi_py;
        double phi_pz;
        double phi_invm;

        double Ds_eta;
        double Ds_phi;
        double Ds_p;
        double Ds_pt;
        double Ds_px;
        double Ds_py;
        double Ds_pz;
        double Ds_invm;

        double Kp_pp;
        double Kp_pl;
        double Km_pp;
        double Km_pl;

        double phi_pp;
        double phi_pl;
        double pi_pp;
        double pi_pl;

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
        double phiFit_chi2;
        double phiFit_ndof;
        double phiFit_chi2ndof;
        double phiFit_vx;
        double phiFit_vy;
        double phiFit_vz;
        double phiFit_vxerr;
        double phiFit_vyerr;
        double phiFit_vzerr;

        double phiFit_Kp_eta;
        double phiFit_Kp_phi;
        double phiFit_Kp_p;
        double phiFit_Kp_pt;
        double phiFit_Kp_px;
        double phiFit_Kp_py;
        double phiFit_Kp_pz;

        double phiFit_Km_eta;
        double phiFit_Km_phi;
        double phiFit_Km_p;
        double phiFit_Km_pt;
        double phiFit_Km_px;
        double phiFit_Km_py;
        double phiFit_Km_pz;

        double phiFit_pi_eta;
        double phiFit_pi_phi;
        double phiFit_pi_p;
        double phiFit_pi_pt;
        double phiFit_pi_px;
        double phiFit_pi_py;
        double phiFit_pi_pz;

        double phiFit_phi_eta;
        double phiFit_phi_phi;
        double phiFit_phi_p;
        double phiFit_phi_pt;
        double phiFit_phi_px;
        double phiFit_phi_py;
        double phiFit_phi_pz;
        double phiFit_phi_invm;

        double phiFit_Ds_eta;
        double phiFit_Ds_phi;
        double phiFit_Ds_p;
        double phiFit_Ds_pt;
        double phiFit_Ds_px;
        double phiFit_Ds_py;
        double phiFit_Ds_pz;
        double phiFit_Ds_invm;

        double phiFit_Kp_pp;
        double phiFit_Kp_pl;
        double phiFit_Km_pp;
        double phiFit_Km_pl;

        double phiFit_phi_pp;
        double phiFit_phi_pl;
        double phiFit_pi_pp;
        double phiFit_pi_pl;

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
        double DsFit_chi2;
        double DsFit_ndof;
        double DsFit_chi2ndof;
        double DsFit_vx;
        double DsFit_vy;
        double DsFit_vz;
        double DsFit_vxerr;
        double DsFit_vyerr;
        double DsFit_vzerr;

        double DsFit_Kp_eta;
        double DsFit_Kp_phi;
        double DsFit_Kp_p;
        double DsFit_Kp_pt;
        double DsFit_Kp_px;
        double DsFit_Kp_py;
        double DsFit_Kp_pz;

        double DsFit_Km_eta;
        double DsFit_Km_phi;
        double DsFit_Km_p;
        double DsFit_Km_pt;
        double DsFit_Km_px;
        double DsFit_Km_py;
        double DsFit_Km_pz;

        double DsFit_pi_eta;
        double DsFit_pi_phi;
        double DsFit_pi_p;
        double DsFit_pi_pt;
        double DsFit_pi_px;
        double DsFit_pi_py;
        double DsFit_pi_pz;

        double DsFit_phi_eta;
        double DsFit_phi_phi;
        double DsFit_phi_p;
        double DsFit_phi_pt;
        double DsFit_phi_px;
        double DsFit_phi_py;
        double DsFit_phi_pz;
        double DsFit_phi_invm;

        double DsFit_Ds_eta;
        double DsFit_Ds_phi;
        double DsFit_Ds_p;
        double DsFit_Ds_pt;
        double DsFit_Ds_px;
        double DsFit_Ds_py;
        double DsFit_Ds_pz;
        double DsFit_Ds_invm;

        double DsFit_Kp_pp;
        double DsFit_Kp_pl;
        double DsFit_Km_pp;
        double DsFit_Km_pl;

        double DsFit_phi_pp;
        double DsFit_phi_pl;
        double DsFit_pi_pp;
        double DsFit_pi_pl;

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

        double DsFit_Mconstraint_Ds_invm;

        /* double Ds_FDxy; */
        /* double Ds_FDxy_Err; */
        /* double Ds_FDxy_Chi2; */
        /* double Ds_FDz; */
        /* double Ds_FDz_Err; */
        /* double Ds_FDz_Chi2; */
        /* double Ds_FD; */
        /* double Ds_FD_Err; */
        /* double Ds_FD_Chi2; */
        /* double Ds_DIRA_angle; */
        /* double Ds_DIRA; */
        /* double Kp_IP; */
        /* double Kp_IP_Err; */
        /* double Kp_IP_Chi2; */
        /* double Km_IP; */
        /* double Km_IP_Err; */
        /* double Km_IP_Chi2; */
        /* double pi_IP; */
        /* double pi_IP_Err; */
        /* double pi_IP_Chi2; */

        /* std::vector<bool> PV_withBS_IsValid_vec; */
        /* std::vector<bool> PV_withBS_IsFake_vec; */
        /* std::vector<double> PV_withBS_chi2_vec; */
        /* std::vector<double> PV_withBS_ndof_vec; */
        /* std::vector<double> PV_withBS_chi2ndof_vec; */
        /* std::vector<double> PV_withBS_X_vec; */
        /* std::vector<double> PV_withBS_Y_vec; */
        /* std::vector<double> PV_withBS_Z_vec; */
        /* std::vector<double> PV_withBS_Xerr_vec; */
        /* std::vector<double> PV_withBS_Yerr_vec; */
        /* std::vector<double> PV_withBS_Zerr_vec; */
        /* std::vector<bool> PV_noBS_IsValid_vec; */
        /* std::vector<bool> PV_noBS_IsFake_vec; */
        /* std::vector<double> PV_noBS_chi2_vec; */
        /* std::vector<double> PV_noBS_ndof_vec; */
        /* std::vector<double> PV_noBS_chi2ndof_vec; */
        /* std::vector<double> PV_noBS_X_vec; */
        /* std::vector<double> PV_noBS_Y_vec; */
        /* std::vector<double> PV_noBS_Z_vec; */
        /* std::vector<double> PV_noBS_Xerr_vec; */
        /* std::vector<double> PV_noBS_Yerr_vec; */
        /* std::vector<double> PV_noBS_Zerr_vec; */

        std::vector<bool> Kp_isIsolatedChargedHadron_vec;
        std::vector<int> Kp_charge_vec;
        std::vector<double> Kp_eta_vec;
        std::vector<double> Kp_phi_vec;
        std::vector<double> Kp_vx_vec;
        std::vector<double> Kp_vy_vec;
        std::vector<double> Kp_vz_vec;
        std::vector<double> Kp_p_vec;
        std::vector<double> Kp_pt_vec;
        std::vector<double> Kp_px_vec;
        std::vector<double> Kp_py_vec;
        std::vector<double> Kp_pz_vec;

        std::vector<bool> Km_isIsolatedChargedHadron_vec;
        std::vector<int> Km_charge_vec;
        std::vector<double> Km_eta_vec;
        std::vector<double> Km_phi_vec;
        std::vector<double> Km_vx_vec;
        std::vector<double> Km_vy_vec;
        std::vector<double> Km_vz_vec;
        std::vector<double> Km_p_vec;
        std::vector<double> Km_pt_vec;
        std::vector<double> Km_px_vec;
        std::vector<double> Km_py_vec;
        std::vector<double> Km_pz_vec;

        std::vector<bool> pi_isIsolatedChargedHadron_vec;
        std::vector<int> pi_charge_vec;
        std::vector<double> pi_eta_vec;
        std::vector<double> pi_phi_vec;
        std::vector<double> pi_vx_vec;
        std::vector<double> pi_vy_vec;
        std::vector<double> pi_vz_vec;
        std::vector<double> pi_p_vec;
        std::vector<double> pi_pt_vec;
        std::vector<double> pi_px_vec;
        std::vector<double> pi_py_vec;
        std::vector<double> pi_pz_vec;

        std::vector<double> phi_eta_vec;
        std::vector<double> phi_phi_vec;
        std::vector<double> phi_p_vec;
        std::vector<double> phi_pt_vec;
        std::vector<double> phi_px_vec;
        std::vector<double> phi_py_vec;
        std::vector<double> phi_pz_vec;
        std::vector<double> phi_invm_vec;

        std::vector<double> Ds_eta_vec;
        std::vector<double> Ds_phi_vec;
        std::vector<double> Ds_p_vec;
        std::vector<double> Ds_pt_vec;
        std::vector<double> Ds_px_vec;
        std::vector<double> Ds_py_vec;
        std::vector<double> Ds_pz_vec;
        std::vector<double> Ds_invm_vec;

        std::vector<double> Kp_pp_vec;
        std::vector<double> Kp_pl_vec;
        std::vector<double> Km_pp_vec;
        std::vector<double> Km_pl_vec;

        std::vector<double> phi_pp_vec;
        std::vector<double> phi_pl_vec;
        std::vector<double> pi_pp_vec;
        std::vector<double> pi_pl_vec;

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
        std::vector<double> phiFit_chi2_vec;
        std::vector<double> phiFit_ndof_vec;
        std::vector<double> phiFit_chi2ndof_vec;
        std::vector<double> phiFit_vx_vec;
        std::vector<double> phiFit_vy_vec;
        std::vector<double> phiFit_vz_vec;
        std::vector<double> phiFit_vxerr_vec;
        std::vector<double> phiFit_vyerr_vec;
        std::vector<double> phiFit_vzerr_vec;

        std::vector<double> phiFit_Kp_eta_vec;
        std::vector<double> phiFit_Kp_phi_vec;
        std::vector<double> phiFit_Kp_p_vec;
        std::vector<double> phiFit_Kp_pt_vec;
        std::vector<double> phiFit_Kp_px_vec;
        std::vector<double> phiFit_Kp_py_vec;
        std::vector<double> phiFit_Kp_pz_vec;

        std::vector<double> phiFit_Km_eta_vec;
        std::vector<double> phiFit_Km_phi_vec;
        std::vector<double> phiFit_Km_p_vec;
        std::vector<double> phiFit_Km_pt_vec;
        std::vector<double> phiFit_Km_px_vec;
        std::vector<double> phiFit_Km_py_vec;
        std::vector<double> phiFit_Km_pz_vec;

        std::vector<double> phiFit_pi_eta_vec;
        std::vector<double> phiFit_pi_phi_vec;
        std::vector<double> phiFit_pi_p_vec;
        std::vector<double> phiFit_pi_pt_vec;
        std::vector<double> phiFit_pi_px_vec;
        std::vector<double> phiFit_pi_py_vec;
        std::vector<double> phiFit_pi_pz_vec;

        std::vector<double> phiFit_phi_eta_vec;
        std::vector<double> phiFit_phi_phi_vec;
        std::vector<double> phiFit_phi_p_vec;
        std::vector<double> phiFit_phi_pt_vec;
        std::vector<double> phiFit_phi_px_vec;
        std::vector<double> phiFit_phi_py_vec;
        std::vector<double> phiFit_phi_pz_vec;
        std::vector<double> phiFit_phi_invm_vec;

        std::vector<double> phiFit_Ds_eta_vec;
        std::vector<double> phiFit_Ds_phi_vec;
        std::vector<double> phiFit_Ds_p_vec;
        std::vector<double> phiFit_Ds_pt_vec;
        std::vector<double> phiFit_Ds_px_vec;
        std::vector<double> phiFit_Ds_py_vec;
        std::vector<double> phiFit_Ds_pz_vec;
        std::vector<double> phiFit_Ds_invm_vec;

        std::vector<double> phiFit_Kp_pp_vec;
        std::vector<double> phiFit_Kp_pl_vec;
        std::vector<double> phiFit_Km_pp_vec;
        std::vector<double> phiFit_Km_pl_vec;

        std::vector<double> phiFit_phi_pp_vec;
        std::vector<double> phiFit_phi_pl_vec;
        std::vector<double> phiFit_pi_pp_vec;
        std::vector<double> phiFit_pi_pl_vec;

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
        std::vector<double> DsFit_chi2_vec;
        std::vector<double> DsFit_ndof_vec;
        std::vector<double> DsFit_chi2ndof_vec;
        std::vector<double> DsFit_vx_vec;
        std::vector<double> DsFit_vy_vec;
        std::vector<double> DsFit_vz_vec;
        std::vector<double> DsFit_vxerr_vec;
        std::vector<double> DsFit_vyerr_vec;
        std::vector<double> DsFit_vzerr_vec;

        std::vector<double> DsFit_Kp_eta_vec;
        std::vector<double> DsFit_Kp_phi_vec;
        std::vector<double> DsFit_Kp_p_vec;
        std::vector<double> DsFit_Kp_pt_vec;
        std::vector<double> DsFit_Kp_px_vec;
        std::vector<double> DsFit_Kp_py_vec;
        std::vector<double> DsFit_Kp_pz_vec;

        std::vector<double> DsFit_Km_eta_vec;
        std::vector<double> DsFit_Km_phi_vec;
        std::vector<double> DsFit_Km_p_vec;
        std::vector<double> DsFit_Km_pt_vec;
        std::vector<double> DsFit_Km_px_vec;
        std::vector<double> DsFit_Km_py_vec;
        std::vector<double> DsFit_Km_pz_vec;

        std::vector<double> DsFit_pi_eta_vec;
        std::vector<double> DsFit_pi_phi_vec;
        std::vector<double> DsFit_pi_p_vec;
        std::vector<double> DsFit_pi_pt_vec;
        std::vector<double> DsFit_pi_px_vec;
        std::vector<double> DsFit_pi_py_vec;
        std::vector<double> DsFit_pi_pz_vec;

        std::vector<double> DsFit_phi_eta_vec;
        std::vector<double> DsFit_phi_phi_vec;
        std::vector<double> DsFit_phi_p_vec;
        std::vector<double> DsFit_phi_pt_vec;
        std::vector<double> DsFit_phi_px_vec;
        std::vector<double> DsFit_phi_py_vec;
        std::vector<double> DsFit_phi_pz_vec;
        std::vector<double> DsFit_phi_invm_vec;

        std::vector<double> DsFit_Ds_eta_vec;
        std::vector<double> DsFit_Ds_phi_vec;
        std::vector<double> DsFit_Ds_p_vec;
        std::vector<double> DsFit_Ds_pt_vec;
        std::vector<double> DsFit_Ds_px_vec;
        std::vector<double> DsFit_Ds_py_vec;
        std::vector<double> DsFit_Ds_pz_vec;
        std::vector<double> DsFit_Ds_invm_vec;

        std::vector<double> DsFit_Kp_pp_vec;
        std::vector<double> DsFit_Kp_pl_vec;
        std::vector<double> DsFit_Km_pp_vec;
        std::vector<double> DsFit_Km_pl_vec;

        std::vector<double> DsFit_phi_pp_vec;
        std::vector<double> DsFit_phi_pl_vec;
        std::vector<double> DsFit_pi_pp_vec;
        std::vector<double> DsFit_pi_pl_vec;

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

        std::vector<double> DsFit_Mconstraint_Ds_invm_vec;

        /* std::vector<double> Ds_FDxy_vec; */
        /* std::vector<double> Ds_FDxy_Err_vec; */
        /* std::vector<double> Ds_FDxy_Chi2_vec; */
        /* std::vector<double> Ds_FDz_vec; */
        /* std::vector<double> Ds_FDz_Err_vec; */
        /* std::vector<double> Ds_FDz_Chi2_vec; */
        /* std::vector<double> Ds_FD_vec; */
        /* std::vector<double> Ds_FD_Err_vec; */
        /* std::vector<double> Ds_FD_Chi2_vec; */
        /* std::vector<double> Ds_DIRA_angle_vec; */
        /* std::vector<double> Ds_DIRA_vec; */
        /* std::vector<double> Kp_IP_vec; */
        /* std::vector<double> Kp_IP_Err_vec; */
        /* std::vector<double> Kp_IP_Chi2_vec; */
        /* std::vector<double> Km_IP_vec; */
        /* std::vector<double> Km_IP_Err_vec; */
        /* std::vector<double> Km_IP_Chi2_vec; */
        /* std::vector<double> pi_IP_vec; */
        /* std::vector<double> pi_IP_Err_vec; */
        /* std::vector<double> pi_IP_Chi2_vec; */

        std::vector<bool> best_Kp_isIsolatedChargedHadron_vec;
        std::vector<int> best_Kp_charge_vec;
        std::vector<double> best_Kp_eta_vec;
        std::vector<double> best_Kp_phi_vec;
        std::vector<double> best_Kp_vx_vec;
        std::vector<double> best_Kp_vy_vec;
        std::vector<double> best_Kp_vz_vec;
        std::vector<double> best_Kp_p_vec;
        std::vector<double> best_Kp_pt_vec;
        std::vector<double> best_Kp_px_vec;
        std::vector<double> best_Kp_py_vec;
        std::vector<double> best_Kp_pz_vec;

        std::vector<bool> best_Km_isIsolatedChargedHadron_vec;
        std::vector<int> best_Km_charge_vec;
        std::vector<double> best_Km_eta_vec;
        std::vector<double> best_Km_phi_vec;
        std::vector<double> best_Km_vx_vec;
        std::vector<double> best_Km_vy_vec;
        std::vector<double> best_Km_vz_vec;
        std::vector<double> best_Km_p_vec;
        std::vector<double> best_Km_pt_vec;
        std::vector<double> best_Km_px_vec;
        std::vector<double> best_Km_py_vec;
        std::vector<double> best_Km_pz_vec;

        std::vector<bool> best_pi_isIsolatedChargedHadron_vec;
        std::vector<int> best_pi_charge_vec;
        std::vector<double> best_pi_eta_vec;
        std::vector<double> best_pi_phi_vec;
        std::vector<double> best_pi_vx_vec;
        std::vector<double> best_pi_vy_vec;
        std::vector<double> best_pi_vz_vec;
        std::vector<double> best_pi_p_vec;
        std::vector<double> best_pi_pt_vec;
        std::vector<double> best_pi_px_vec;
        std::vector<double> best_pi_py_vec;
        std::vector<double> best_pi_pz_vec;

        std::vector<double> best_phi_eta_vec;
        std::vector<double> best_phi_phi_vec;
        std::vector<double> best_phi_p_vec;
        std::vector<double> best_phi_pt_vec;
        std::vector<double> best_phi_px_vec;
        std::vector<double> best_phi_py_vec;
        std::vector<double> best_phi_pz_vec;
        std::vector<double> best_phi_invm_vec;

        std::vector<double> best_Ds_eta_vec;
        std::vector<double> best_Ds_phi_vec;
        std::vector<double> best_Ds_p_vec;
        std::vector<double> best_Ds_pt_vec;
        std::vector<double> best_Ds_px_vec;
        std::vector<double> best_Ds_py_vec;
        std::vector<double> best_Ds_pz_vec;
        std::vector<double> best_Ds_invm_vec;

        std::vector<double> best_Kp_pp_vec;
        std::vector<double> best_Kp_pl_vec;
        std::vector<double> best_Km_pp_vec;
        std::vector<double> best_Km_pl_vec;

        std::vector<double> best_phi_pp_vec;
        std::vector<double> best_phi_pl_vec;
        std::vector<double> best_pi_pp_vec;
        std::vector<double> best_pi_pl_vec;

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
        std::vector<double> best_phiFit_chi2_vec;
        std::vector<double> best_phiFit_ndof_vec;
        std::vector<double> best_phiFit_chi2ndof_vec;
        std::vector<double> best_phiFit_vx_vec;
        std::vector<double> best_phiFit_vy_vec;
        std::vector<double> best_phiFit_vz_vec;
        std::vector<double> best_phiFit_vxerr_vec;
        std::vector<double> best_phiFit_vyerr_vec;
        std::vector<double> best_phiFit_vzerr_vec;

        std::vector<double> best_phiFit_Kp_eta_vec;
        std::vector<double> best_phiFit_Kp_phi_vec;
        std::vector<double> best_phiFit_Kp_p_vec;
        std::vector<double> best_phiFit_Kp_pt_vec;
        std::vector<double> best_phiFit_Kp_px_vec;
        std::vector<double> best_phiFit_Kp_py_vec;
        std::vector<double> best_phiFit_Kp_pz_vec;

        std::vector<double> best_phiFit_Km_eta_vec;
        std::vector<double> best_phiFit_Km_phi_vec;
        std::vector<double> best_phiFit_Km_p_vec;
        std::vector<double> best_phiFit_Km_pt_vec;
        std::vector<double> best_phiFit_Km_px_vec;
        std::vector<double> best_phiFit_Km_py_vec;
        std::vector<double> best_phiFit_Km_pz_vec;

        std::vector<double> best_phiFit_pi_eta_vec;
        std::vector<double> best_phiFit_pi_phi_vec;
        std::vector<double> best_phiFit_pi_p_vec;
        std::vector<double> best_phiFit_pi_pt_vec;
        std::vector<double> best_phiFit_pi_px_vec;
        std::vector<double> best_phiFit_pi_py_vec;
        std::vector<double> best_phiFit_pi_pz_vec;

        std::vector<double> best_phiFit_phi_eta_vec;
        std::vector<double> best_phiFit_phi_phi_vec;
        std::vector<double> best_phiFit_phi_p_vec;
        std::vector<double> best_phiFit_phi_pt_vec;
        std::vector<double> best_phiFit_phi_px_vec;
        std::vector<double> best_phiFit_phi_py_vec;
        std::vector<double> best_phiFit_phi_pz_vec;
        std::vector<double> best_phiFit_phi_invm_vec;

        std::vector<double> best_phiFit_Ds_eta_vec;
        std::vector<double> best_phiFit_Ds_phi_vec;
        std::vector<double> best_phiFit_Ds_p_vec;
        std::vector<double> best_phiFit_Ds_pt_vec;
        std::vector<double> best_phiFit_Ds_px_vec;
        std::vector<double> best_phiFit_Ds_py_vec;
        std::vector<double> best_phiFit_Ds_pz_vec;
        std::vector<double> best_phiFit_Ds_invm_vec;

        std::vector<double> best_phiFit_Kp_pp_vec;
        std::vector<double> best_phiFit_Kp_pl_vec;
        std::vector<double> best_phiFit_Km_pp_vec;
        std::vector<double> best_phiFit_Km_pl_vec;

        std::vector<double> best_phiFit_phi_pp_vec;
        std::vector<double> best_phiFit_phi_pl_vec;
        std::vector<double> best_phiFit_pi_pp_vec;
        std::vector<double> best_phiFit_pi_pl_vec;

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
        std::vector<double> best_DsFit_chi2_vec;
        std::vector<double> best_DsFit_ndof_vec;
        std::vector<double> best_DsFit_chi2ndof_vec;
        std::vector<double> best_DsFit_vx_vec;
        std::vector<double> best_DsFit_vy_vec;
        std::vector<double> best_DsFit_vz_vec;
        std::vector<double> best_DsFit_vxerr_vec;
        std::vector<double> best_DsFit_vyerr_vec;
        std::vector<double> best_DsFit_vzerr_vec;

        std::vector<double> best_DsFit_Kp_eta_vec;
        std::vector<double> best_DsFit_Kp_phi_vec;
        std::vector<double> best_DsFit_Kp_p_vec;
        std::vector<double> best_DsFit_Kp_pt_vec;
        std::vector<double> best_DsFit_Kp_px_vec;
        std::vector<double> best_DsFit_Kp_py_vec;
        std::vector<double> best_DsFit_Kp_pz_vec;

        std::vector<double> best_DsFit_Km_eta_vec;
        std::vector<double> best_DsFit_Km_phi_vec;
        std::vector<double> best_DsFit_Km_p_vec;
        std::vector<double> best_DsFit_Km_pt_vec;
        std::vector<double> best_DsFit_Km_px_vec;
        std::vector<double> best_DsFit_Km_py_vec;
        std::vector<double> best_DsFit_Km_pz_vec;

        std::vector<double> best_DsFit_pi_eta_vec;
        std::vector<double> best_DsFit_pi_phi_vec;
        std::vector<double> best_DsFit_pi_p_vec;
        std::vector<double> best_DsFit_pi_pt_vec;
        std::vector<double> best_DsFit_pi_px_vec;
        std::vector<double> best_DsFit_pi_py_vec;
        std::vector<double> best_DsFit_pi_pz_vec;

        std::vector<double> best_DsFit_phi_eta_vec;
        std::vector<double> best_DsFit_phi_phi_vec;
        std::vector<double> best_DsFit_phi_p_vec;
        std::vector<double> best_DsFit_phi_pt_vec;
        std::vector<double> best_DsFit_phi_px_vec;
        std::vector<double> best_DsFit_phi_py_vec;
        std::vector<double> best_DsFit_phi_pz_vec;
        std::vector<double> best_DsFit_phi_invm_vec;

        std::vector<double> best_DsFit_Ds_eta_vec;
        std::vector<double> best_DsFit_Ds_phi_vec;
        std::vector<double> best_DsFit_Ds_p_vec;
        std::vector<double> best_DsFit_Ds_pt_vec;
        std::vector<double> best_DsFit_Ds_px_vec;
        std::vector<double> best_DsFit_Ds_py_vec;
        std::vector<double> best_DsFit_Ds_pz_vec;
        std::vector<double> best_DsFit_Ds_invm_vec;

        std::vector<double> best_DsFit_Kp_pp_vec;
        std::vector<double> best_DsFit_Kp_pl_vec;
        std::vector<double> best_DsFit_Km_pp_vec;
        std::vector<double> best_DsFit_Km_pl_vec;

        std::vector<double> best_DsFit_phi_pp_vec;
        std::vector<double> best_DsFit_phi_pl_vec;
        std::vector<double> best_DsFit_pi_pp_vec;
        std::vector<double> best_DsFit_pi_pl_vec;

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

        std::vector<double> best_DsFit_Mconstraint_Ds_invm_vec;

        /* std::vector<double> best_Ds_FDxy_vec; */
        /* std::vector<double> best_Ds_FDxy_Err_vec; */
        /* std::vector<double> best_Ds_FDxy_Chi2_vec; */
        /* std::vector<double> best_Ds_FDz_vec; */
        /* std::vector<double> best_Ds_FDz_Err_vec; */
        /* std::vector<double> best_Ds_FDz_Chi2_vec; */
        /* std::vector<double> best_Ds_FD_vec; */
        /* std::vector<double> best_Ds_FD_Err_vec; */
        /* std::vector<double> best_Ds_FD_Chi2_vec; */
        /* std::vector<double> best_Ds_DIRA_angle_vec; */
        /* std::vector<double> best_Ds_DIRA_vec; */
        /* std::vector<double> best_Kp_IP_vec; */
        /* std::vector<double> best_Kp_IP_Err_vec; */
        /* std::vector<double> best_Kp_IP_Chi2_vec; */
        /* std::vector<double> best_Km_IP_vec; */
        /* std::vector<double> best_Km_IP_Err_vec; */
        /* std::vector<double> best_Km_IP_Chi2_vec; */
        /* std::vector<double> best_pi_IP_vec; */
        /* std::vector<double> best_pi_IP_Err_vec; */
        /* std::vector<double> best_pi_IP_Chi2_vec; */

};

#endif

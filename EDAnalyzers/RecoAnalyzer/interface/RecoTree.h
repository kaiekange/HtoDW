#ifndef RECOTREE_H 
#define RECOTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class RecoTree 
{
    public:

        RecoTree(TTree *tree_);

        TTree *tree;

        void Init();
        void CreateBranches();
        void BS_Reset();
        void PV_Reset();
        void Kp_Reset();
        void Km_Reset();
        void phi_Reset();
        void pi_Reset();
        void Ds_Reset();
        void PV_noDs_Reset();
        void Fill_Vector();
        void Best_Fill_Vector(int idxmax);

        int BS_type;
        double BS_x0;
        double BS_y0;
        double BS_z0;
        double BS_sigmaZ;
        double BS_dxdz;
        double BS_dydz;
        double BS_BWX;
        double BS_BWY;
        double BS_x0err;
        double BS_y0err;
        double BS_z0err;
        double BS_sigmaZ0err;
        double BS_dxdzerr;
        double BS_dydzerr;
        double BS_BWXerr;
        double BS_BWYerr;
        double BS_emitX;
        double BS_emitY;
        double BS_betaStar;

        double PV_vx;
        double PV_vy;
        double PV_vz;
        double PV_vxerr;
        double PV_vyerr;
        double PV_vzerr;

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

        double alpha_phi;
        double beta_phi;
        double APvar_phi;

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

        double alpha_Ds;
        double beta_Ds;
        double APvar_Ds;

        double Ds_primvtx_FDxy;
        double Ds_primvtx_FDz;
        double Ds_primvtx_FD;
        double Ds_primvtx_FDxyerr;
        double Ds_primvtx_FDxychi2;
        double Ds_primvtx_FDzerr;
        double Ds_primvtx_FDzchi2;
        double Ds_primvtx_FDerr;
        double Ds_primvtx_FDchi2;
        double Ds_primvtx_dira;
        double Ds_primvtx_dira_angle;
        double Kp_primvtx_ip;
        double Kp_primvtx_iperr;
        double Kp_primvtx_ipchi2;
        double Km_primvtx_ip;
        double Km_primvtx_iperr;
        double Km_primvtx_ipchi2;
        double pi_primvtx_ip;
        double pi_primvtx_iperr;
        double pi_primvtx_ipchi2;

        double Ds_Iso_R0p3;
        double Ds_Iso_R0p4;
        double Ds_IsoRel_R0p3;
        double Ds_IsoRel_R0p4;

        bool PV_noDs_withBS_IsValid;
        bool PV_noDs_withBS_IsFake;
        double PV_noDs_withBS_chi2;
        double PV_noDs_withBS_ndof;
        double PV_noDs_withBS_chi2ndof;
        double PV_noDs_withBS_x;
        double PV_noDs_withBS_y;
        double PV_noDs_withBS_z;
        double PV_noDs_withBS_xerr;
        double PV_noDs_withBS_yerr;
        double PV_noDs_withBS_zerr;

        bool PV_noDs_noBS_IsValid;
        bool PV_noDs_noBS_IsFake;
        double PV_noDs_noBS_chi2;
        double PV_noDs_noBS_ndof;
        double PV_noDs_noBS_chi2ndof;
        double PV_noDs_noBS_x;
        double PV_noDs_noBS_y;
        double PV_noDs_noBS_z;
        double PV_noDs_noBS_xerr;
        double PV_noDs_noBS_yerr;
        double PV_noDs_noBS_zerr;

        double Ds_PVnoDs_FDxy;
        double Ds_PVnoDs_FDz;
        double Ds_PVnoDs_FD;
        double Ds_PVnoDs_FDxyerr;
        double Ds_PVnoDs_FDxychi2;
        double Ds_PVnoDs_FDzerr;
        double Ds_PVnoDs_FDzchi2;
        double Ds_PVnoDs_FDerr;
        double Ds_PVnoDs_FDchi2;
        double Ds_PVnoDs_dira;
        double Ds_PVnoDs_dira_angle;
        double Kp_PVnoDs_ip;
        double Kp_PVnoDs_iperr;
        double Kp_PVnoDs_ipchi2;
        double Km_PVnoDs_ip;
        double Km_PVnoDs_iperr;
        double Km_PVnoDs_ipchi2;
        double pi_PVnoDs_ip;
        double pi_PVnoDs_iperr;
        double pi_PVnoDs_ipchi2;

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

        std::vector<double> alpha_phi_vec;
        std::vector<double> beta_phi_vec;
        std::vector<double> APvar_phi_vec;

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

        std::vector<double> alpha_Ds_vec;
        std::vector<double> beta_Ds_vec;
        std::vector<double> APvar_Ds_vec;

        std::vector<double> Ds_primvtx_FDxy_vec;
        std::vector<double> Ds_primvtx_FDz_vec;
        std::vector<double> Ds_primvtx_FD_vec;
        std::vector<double> Ds_primvtx_FDxyerr_vec;
        std::vector<double> Ds_primvtx_FDxychi2_vec;
        std::vector<double> Ds_primvtx_FDzerr_vec;
        std::vector<double> Ds_primvtx_FDzchi2_vec;
        std::vector<double> Ds_primvtx_FDerr_vec;
        std::vector<double> Ds_primvtx_FDchi2_vec;
        std::vector<double> Ds_primvtx_dira_vec;
        std::vector<double> Ds_primvtx_dira_angle_vec;
        std::vector<double> Kp_primvtx_ip_vec;
        std::vector<double> Kp_primvtx_iperr_vec;
        std::vector<double> Kp_primvtx_ipchi2_vec;
        std::vector<double> Km_primvtx_ip_vec;
        std::vector<double> Km_primvtx_iperr_vec;
        std::vector<double> Km_primvtx_ipchi2_vec;
        std::vector<double> pi_primvtx_ip_vec;
        std::vector<double> pi_primvtx_iperr_vec;
        std::vector<double> pi_primvtx_ipchi2_vec;

        std::vector<double> Ds_Iso_R0p3_vec;
        std::vector<double> Ds_Iso_R0p4_vec;
        std::vector<double> Ds_IsoRel_R0p3_vec;
        std::vector<double> Ds_IsoRel_R0p4_vec;

        std::vector<bool> PV_noDs_withBS_IsValid_vec;
        std::vector<bool> PV_noDs_withBS_IsFake_vec;
        std::vector<double> PV_noDs_withBS_chi2_vec;
        std::vector<double> PV_noDs_withBS_ndof_vec;
        std::vector<double> PV_noDs_withBS_chi2ndof_vec;
        std::vector<double> PV_noDs_withBS_x_vec;
        std::vector<double> PV_noDs_withBS_y_vec;
        std::vector<double> PV_noDs_withBS_z_vec;
        std::vector<double> PV_noDs_withBS_xerr_vec;
        std::vector<double> PV_noDs_withBS_yerr_vec;
        std::vector<double> PV_noDs_withBS_zerr_vec;

        std::vector<bool> PV_noDs_noBS_IsValid_vec;
        std::vector<bool> PV_noDs_noBS_IsFake_vec;
        std::vector<double> PV_noDs_noBS_chi2_vec;
        std::vector<double> PV_noDs_noBS_ndof_vec;
        std::vector<double> PV_noDs_noBS_chi2ndof_vec;
        std::vector<double> PV_noDs_noBS_x_vec;
        std::vector<double> PV_noDs_noBS_y_vec;
        std::vector<double> PV_noDs_noBS_z_vec;
        std::vector<double> PV_noDs_noBS_xerr_vec;
        std::vector<double> PV_noDs_noBS_yerr_vec;
        std::vector<double> PV_noDs_noBS_zerr_vec;

        std::vector<double> Ds_PVnoDs_FDxy_vec;
        std::vector<double> Ds_PVnoDs_FDz_vec;
        std::vector<double> Ds_PVnoDs_FD_vec;
        std::vector<double> Ds_PVnoDs_FDxyerr_vec;
        std::vector<double> Ds_PVnoDs_FDxychi2_vec;
        std::vector<double> Ds_PVnoDs_FDzerr_vec;
        std::vector<double> Ds_PVnoDs_FDzchi2_vec;
        std::vector<double> Ds_PVnoDs_FDerr_vec;
        std::vector<double> Ds_PVnoDs_FDchi2_vec;
        std::vector<double> Ds_PVnoDs_dira_vec;
        std::vector<double> Ds_PVnoDs_dira_angle_vec;
        std::vector<double> Kp_PVnoDs_ip_vec;
        std::vector<double> Kp_PVnoDs_iperr_vec;
        std::vector<double> Kp_PVnoDs_ipchi2_vec;
        std::vector<double> Km_PVnoDs_ip_vec;
        std::vector<double> Km_PVnoDs_iperr_vec;
        std::vector<double> Km_PVnoDs_ipchi2_vec;
        std::vector<double> pi_PVnoDs_ip_vec;
        std::vector<double> pi_PVnoDs_iperr_vec;
        std::vector<double> pi_PVnoDs_ipchi2_vec;

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

        std::vector<double> best_alpha_phi_vec;
        std::vector<double> best_beta_phi_vec;
        std::vector<double> best_APvar_phi_vec;

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

        std::vector<double> best_alpha_Ds_vec;
        std::vector<double> best_beta_Ds_vec;
        std::vector<double> best_APvar_Ds_vec;

        std::vector<double> best_Ds_primvtx_FDxy_vec;
        std::vector<double> best_Ds_primvtx_FDz_vec;
        std::vector<double> best_Ds_primvtx_FD_vec;
        std::vector<double> best_Ds_primvtx_FDxyerr_vec;
        std::vector<double> best_Ds_primvtx_FDxychi2_vec;
        std::vector<double> best_Ds_primvtx_FDzerr_vec;
        std::vector<double> best_Ds_primvtx_FDzchi2_vec;
        std::vector<double> best_Ds_primvtx_FDerr_vec;
        std::vector<double> best_Ds_primvtx_FDchi2_vec;
        std::vector<double> best_Ds_primvtx_dira_vec;
        std::vector<double> best_Ds_primvtx_dira_angle_vec;
        std::vector<double> best_Kp_primvtx_ip_vec;
        std::vector<double> best_Kp_primvtx_iperr_vec;
        std::vector<double> best_Kp_primvtx_ipchi2_vec;
        std::vector<double> best_Km_primvtx_ip_vec;
        std::vector<double> best_Km_primvtx_iperr_vec;
        std::vector<double> best_Km_primvtx_ipchi2_vec;
        std::vector<double> best_pi_primvtx_ip_vec;
        std::vector<double> best_pi_primvtx_iperr_vec;
        std::vector<double> best_pi_primvtx_ipchi2_vec;

        std::vector<double> best_Ds_Iso_R0p3_vec;
        std::vector<double> best_Ds_Iso_R0p4_vec;
        std::vector<double> best_Ds_IsoRel_R0p3_vec;
        std::vector<double> best_Ds_IsoRel_R0p4_vec;

        std::vector<bool> best_PV_noDs_withBS_IsValid_vec;
        std::vector<bool> best_PV_noDs_withBS_IsFake_vec;
        std::vector<double> best_PV_noDs_withBS_chi2_vec;
        std::vector<double> best_PV_noDs_withBS_ndof_vec;
        std::vector<double> best_PV_noDs_withBS_chi2ndof_vec;
        std::vector<double> best_PV_noDs_withBS_x_vec;
        std::vector<double> best_PV_noDs_withBS_y_vec;
        std::vector<double> best_PV_noDs_withBS_z_vec;
        std::vector<double> best_PV_noDs_withBS_xerr_vec;
        std::vector<double> best_PV_noDs_withBS_yerr_vec;
        std::vector<double> best_PV_noDs_withBS_zerr_vec;

        std::vector<bool> best_PV_noDs_noBS_IsValid_vec;
        std::vector<bool> best_PV_noDs_noBS_IsFake_vec;
        std::vector<double> best_PV_noDs_noBS_chi2_vec;
        std::vector<double> best_PV_noDs_noBS_ndof_vec;
        std::vector<double> best_PV_noDs_noBS_chi2ndof_vec;
        std::vector<double> best_PV_noDs_noBS_x_vec;
        std::vector<double> best_PV_noDs_noBS_y_vec;
        std::vector<double> best_PV_noDs_noBS_z_vec;
        std::vector<double> best_PV_noDs_noBS_xerr_vec;
        std::vector<double> best_PV_noDs_noBS_yerr_vec;
        std::vector<double> best_PV_noDs_noBS_zerr_vec;

        std::vector<double> best_Ds_PVnoDs_FDxy_vec;
        std::vector<double> best_Ds_PVnoDs_FDz_vec;
        std::vector<double> best_Ds_PVnoDs_FD_vec;
        std::vector<double> best_Ds_PVnoDs_FDxyerr_vec;
        std::vector<double> best_Ds_PVnoDs_FDxychi2_vec;
        std::vector<double> best_Ds_PVnoDs_FDzerr_vec;
        std::vector<double> best_Ds_PVnoDs_FDzchi2_vec;
        std::vector<double> best_Ds_PVnoDs_FDerr_vec;
        std::vector<double> best_Ds_PVnoDs_FDchi2_vec;
        std::vector<double> best_Ds_PVnoDs_dira_vec;
        std::vector<double> best_Ds_PVnoDs_dira_angle_vec;
        std::vector<double> best_Kp_PVnoDs_ip_vec;
        std::vector<double> best_Kp_PVnoDs_iperr_vec;
        std::vector<double> best_Kp_PVnoDs_ipchi2_vec;
        std::vector<double> best_Km_PVnoDs_ip_vec;
        std::vector<double> best_Km_PVnoDs_iperr_vec;
        std::vector<double> best_Km_PVnoDs_ipchi2_vec;
        std::vector<double> best_pi_PVnoDs_ip_vec;
        std::vector<double> best_pi_PVnoDs_iperr_vec;
        std::vector<double> best_pi_PVnoDs_ipchi2_vec;

};

#endif

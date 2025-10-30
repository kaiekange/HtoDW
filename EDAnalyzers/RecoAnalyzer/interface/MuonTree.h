#ifndef PVTREE_H 
#define PVTREE_H

#include <TTree.h>
#include <vector>

#define null -777

class MuonTree 
{
    public:

        MuonTree(TTree *tree_);

        TTree *tree;

        void Init();
        void CreateBranches();
        void Gen_Reset();
        void BS_Reset();
        void PV_Reset();
        void Match_Reset();
        void Match_Fill_Vector();
        void mu_Reset();
        void Fill_Vector();

        // Gen particles
        int num_Gen_mu;
        int num_Gen_nu;
        int num_Gen_W;

        double Gen_mu_eta;
        double Gen_mu_phi;
        double Gen_mu_vx;
        double Gen_mu_vy;
        double Gen_mu_vz;
        double Gen_mu_p;
        double Gen_mu_pt;
        double Gen_mu_px;
        double Gen_mu_py;
        double Gen_mu_pz;

        double Gen_nu_eta;
        double Gen_nu_phi;
        double Gen_nu_vx;
        double Gen_nu_vy;
        double Gen_nu_vz;
        double Gen_nu_p;
        double Gen_nu_pt;
        double Gen_nu_px;
        double Gen_nu_py;
        double Gen_nu_pz;

        double Gen_W_eta;
        double Gen_W_phi;
        double Gen_W_vx;
        double Gen_W_vy;
        double Gen_W_vz;
        double Gen_W_p;
        double Gen_W_pt;
        double Gen_W_px;
        double Gen_W_py;
        double Gen_W_pz;

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

        // Matched particles
        int num_match_mu;
        int match_mu_idx;

        // Original info
        float match_mu_dR;
        int match_mu_charge;
        double match_mu_eta;
        double match_mu_phi;
        double match_mu_vx;
        double match_mu_vy;
        double match_mu_vz;
        double match_mu_p;
        double match_mu_pt;
        double match_mu_px;
        double match_mu_py;
        double match_mu_pz;
        bool match_mu_isHighPt;
        bool match_mu_isLoose;
        bool match_mu_isMedium;
        bool match_mu_isSoft;
        bool match_mu_isTight;
        bool match_mu_isPF;
        bool match_mu_isTracker;
        bool match_mu_isGlobal;
        double match_mu_IsoR03_sumChargedHadronPt;
        double match_mu_IsoR03_sumChargedParticlePt;
        double match_mu_IsoR03_sumNeutralHadronEt;
        double match_mu_IsoR03_sumPhotonEt;
        double match_mu_IsoR03_sumPUPt;
        double match_mu_PFIsoR03;
        double match_mu_IsoR04_sumChargedHadronPt;
        double match_mu_IsoR04_sumChargedParticlePt;
        double match_mu_IsoR04_sumNeutralHadronEt;
        double match_mu_IsoR04_sumPhotonEt;
        double match_mu_IsoR04_sumPUPt;
        double match_mu_PFIsoR04;
        double match_mu_primvtx_dxy;
        double match_mu_primvtx_dxyerr;
        double match_mu_primvtx_dz;
        double match_mu_primvtx_dzerr;
        double match_mu_primvtx_ip;
        double match_mu_primvtx_iperr;
        double match_mu_primvtx_ipchi2;

        std::vector<float> match_mu_dR_vec;
        std::vector<int> match_mu_charge_vec;
        std::vector<double> match_mu_eta_vec;
        std::vector<double> match_mu_phi_vec;
        std::vector<double> match_mu_vx_vec;
        std::vector<double> match_mu_vy_vec;
        std::vector<double> match_mu_vz_vec;
        std::vector<double> match_mu_p_vec;
        std::vector<double> match_mu_pt_vec;
        std::vector<double> match_mu_px_vec;
        std::vector<double> match_mu_py_vec;
        std::vector<double> match_mu_pz_vec;
        std::vector<bool> match_mu_isHighPt_vec;
        std::vector<bool> match_mu_isLoose_vec;
        std::vector<bool> match_mu_isMedium_vec;
        std::vector<bool> match_mu_isSoft_vec;
        std::vector<bool> match_mu_isTight_vec;
        std::vector<bool> match_mu_isPF_vec;
        std::vector<bool> match_mu_isTracker_vec;
        std::vector<bool> match_mu_isGlobal_vec;
        std::vector<double> match_mu_IsoR03_sumChargedHadronPt_vec;
        std::vector<double> match_mu_IsoR03_sumChargedParticlePt_vec;
        std::vector<double> match_mu_IsoR03_sumNeutralHadronEt_vec;
        std::vector<double> match_mu_IsoR03_sumPhotonEt_vec;
        std::vector<double> match_mu_IsoR03_sumPUPt_vec;
        std::vector<double> match_mu_PFIsoR03_vec;
        std::vector<double> match_mu_IsoR04_sumChargedHadronPt_vec;
        std::vector<double> match_mu_IsoR04_sumChargedParticlePt_vec;
        std::vector<double> match_mu_IsoR04_sumNeutralHadronEt_vec;
        std::vector<double> match_mu_IsoR04_sumPhotonEt_vec;
        std::vector<double> match_mu_IsoR04_sumPUPt_vec;
        std::vector<double> match_mu_PFIsoR04_vec;
        std::vector<double> match_mu_primvtx_dxy_vec;
        std::vector<double> match_mu_primvtx_dxyerr_vec;
        std::vector<double> match_mu_primvtx_dz_vec;
        std::vector<double> match_mu_primvtx_dzerr_vec;
        std::vector<double> match_mu_primvtx_ip_vec;
        std::vector<double> match_mu_primvtx_iperr_vec;
        std::vector<double> match_mu_primvtx_ipchi2_vec;

        int num_mu;
        
        int mu_charge;
        double mu_eta;
        double mu_phi;
        double mu_vx;
        double mu_vy;
        double mu_vz;
        double mu_p;
        double mu_pt;
        double mu_px;
        double mu_py;
        double mu_pz;
        bool mu_isHighPt;
        bool mu_isLoose;
        bool mu_isMedium;
        bool mu_isSoft;
        bool mu_isTight;
        bool mu_isPF;
        bool mu_isTracker;
        bool mu_isGlobal;
        double mu_IsoR03_sumChargedHadronPt;
        double mu_IsoR03_sumChargedParticlePt;
        double mu_IsoR03_sumNeutralHadronEt;
        double mu_IsoR03_sumPhotonEt;
        double mu_IsoR03_sumPUPt;
        double mu_PFIsoR03;
        double mu_IsoR04_sumChargedHadronPt;
        double mu_IsoR04_sumChargedParticlePt;
        double mu_IsoR04_sumNeutralHadronEt;
        double mu_IsoR04_sumPhotonEt;
        double mu_IsoR04_sumPUPt;
        double mu_PFIsoR04;
        double mu_primvtx_dxy;
        double mu_primvtx_dxyerr;
        double mu_primvtx_dz;
        double mu_primvtx_dzerr;
        double mu_primvtx_ip;
        double mu_primvtx_iperr;
        double mu_primvtx_ipchi2;
        
        bool mu_match;

        std::vector<int> mu_charge_vec;
        std::vector<double> mu_eta_vec;
        std::vector<double> mu_phi_vec;
        std::vector<double> mu_vx_vec;
        std::vector<double> mu_vy_vec;
        std::vector<double> mu_vz_vec;
        std::vector<double> mu_p_vec;
        std::vector<double> mu_pt_vec;
        std::vector<double> mu_px_vec;
        std::vector<double> mu_py_vec;
        std::vector<double> mu_pz_vec;
        std::vector<bool> mu_isHighPt_vec;
        std::vector<bool> mu_isLoose_vec;
        std::vector<bool> mu_isMedium_vec;
        std::vector<bool> mu_isSoft_vec;
        std::vector<bool> mu_isTight_vec;
        std::vector<bool> mu_isPF_vec;
        std::vector<bool> mu_isTracker_vec;
        std::vector<bool> mu_isGlobal_vec;
        std::vector<double> mu_IsoR03_sumChargedHadronPt_vec;
        std::vector<double> mu_IsoR03_sumChargedParticlePt_vec;
        std::vector<double> mu_IsoR03_sumNeutralHadronEt_vec;
        std::vector<double> mu_IsoR03_sumPhotonEt_vec;
        std::vector<double> mu_IsoR03_sumPUPt_vec;
        std::vector<double> mu_PFIsoR03_vec;
        std::vector<double> mu_IsoR04_sumChargedHadronPt_vec;
        std::vector<double> mu_IsoR04_sumChargedParticlePt_vec;
        std::vector<double> mu_IsoR04_sumNeutralHadronEt_vec;
        std::vector<double> mu_IsoR04_sumPhotonEt_vec;
        std::vector<double> mu_IsoR04_sumPUPt_vec;
        std::vector<double> mu_PFIsoR04_vec;
        std::vector<double> mu_primvtx_dxy_vec;
        std::vector<double> mu_primvtx_dxyerr_vec;
        std::vector<double> mu_primvtx_dz_vec;
        std::vector<double> mu_primvtx_dzerr_vec;
        std::vector<double> mu_primvtx_ip_vec;
        std::vector<double> mu_primvtx_iperr_vec;
        std::vector<double> mu_primvtx_ipchi2_vec;

        std::vector<bool> mu_match_vec;
};

#endif

#ifndef RECO_DSSTW_KKPI_TREE
#define RECO_DSSTW_KKPI_TREE

#include <TTree.h>

#define null -777

class RecoAll_DsstW_KKpimu_Tree
{
public:
	RecoAll_DsstW_KKpimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;

	Float_t Gen_PV_vx;
	Float_t Gen_PV_vy;
	Float_t Gen_PV_vz;

	Float_t Gen_SV_vx;
	Float_t Gen_SV_vy;
	Float_t Gen_SV_vz;

	Int_t Gen_H_charge;
	Float_t Gen_H_eta;
	Float_t Gen_H_phi;
	Float_t Gen_H_px;
	Float_t Gen_H_py;
	Float_t Gen_H_pz;
	Float_t Gen_H_pt;
	Float_t Gen_H_p;
	Float_t Gen_H_mass;

	Int_t Gen_Dsst_charge;
	Float_t Gen_Dsst_eta;
	Float_t Gen_Dsst_phi;
	Float_t Gen_Dsst_px;
	Float_t Gen_Dsst_py;
	Float_t Gen_Dsst_pz;
	Float_t Gen_Dsst_pt;
	Float_t Gen_Dsst_p;
	Float_t Gen_Dsst_mass;

	Int_t Gen_W_charge;
	Float_t Gen_W_eta;
	Float_t Gen_W_phi;
	Float_t Gen_W_px;
	Float_t Gen_W_py;
	Float_t Gen_W_pz;
	Float_t Gen_W_pt;
	Float_t Gen_W_p;
	Float_t Gen_W_mass;

	Int_t Gen_Ds_charge;
	Float_t Gen_Ds_eta;
	Float_t Gen_Ds_phi;
	Float_t Gen_Ds_px;
	Float_t Gen_Ds_py;
	Float_t Gen_Ds_pz;
	Float_t Gen_Ds_pt;
	Float_t Gen_Ds_p;
	Float_t Gen_Ds_mass;

	Int_t Gen_gamma_charge;
	Float_t Gen_gamma_eta;
	Float_t Gen_gamma_phi;
	Float_t Gen_gamma_px;
	Float_t Gen_gamma_py;
	Float_t Gen_gamma_pz;
	Float_t Gen_gamma_pt;
	Float_t Gen_gamma_p;
	Float_t Gen_gamma_mass;

	Int_t Gen_Kp_charge;
	Float_t Gen_Kp_eta;
	Float_t Gen_Kp_phi;
	Float_t Gen_Kp_px;
	Float_t Gen_Kp_py;
	Float_t Gen_Kp_pz;
	Float_t Gen_Kp_pt;
	Float_t Gen_Kp_p;
	Float_t Gen_Kp_mass;

	Int_t Gen_Km_charge;
	Float_t Gen_Km_eta;
	Float_t Gen_Km_phi;
	Float_t Gen_Km_px;
	Float_t Gen_Km_py;
	Float_t Gen_Km_pz;
	Float_t Gen_Km_pt;
	Float_t Gen_Km_p;
	Float_t Gen_Km_mass;

	Int_t Gen_pi_charge;
	Float_t Gen_pi_eta;
	Float_t Gen_pi_phi;
	Float_t Gen_pi_px;
	Float_t Gen_pi_py;
	Float_t Gen_pi_pz;
	Float_t Gen_pi_pt;
	Float_t Gen_pi_p;
	Float_t Gen_pi_mass;

	Int_t Gen_mu_charge;
	Float_t Gen_mu_eta;
	Float_t Gen_mu_phi;
	Float_t Gen_mu_px;
	Float_t Gen_mu_py;
	Float_t Gen_mu_pz;
	Float_t Gen_mu_pt;
	Float_t Gen_mu_p;
	Float_t Gen_mu_mass;

	Int_t Gen_Dsmu_charge;
	Float_t Gen_Dsmu_eta;
	Float_t Gen_Dsmu_phi;
	Float_t Gen_Dsmu_px;
	Float_t Gen_Dsmu_py;
	Float_t Gen_Dsmu_pz;
	Float_t Gen_Dsmu_pt;
	Float_t Gen_Dsmu_p;
	Float_t Gen_Dsmu_mass;

	Int_t Gen_vis_charge;
	Float_t Gen_vis_eta;
	Float_t Gen_vis_phi;
	Float_t Gen_vis_px;
	Float_t Gen_vis_py;
	Float_t Gen_vis_pz;
	Float_t Gen_vis_pt;
	Float_t Gen_vis_p;
	Float_t Gen_vis_mass;

	Float_t Gen_dR_Kp_Km;
	Float_t Gen_dR_Kp_pi;
	Float_t Gen_dR_Kp_Ds;
	Float_t Gen_dR_Km_pi;
	Float_t Gen_dR_Km_Ds;
	Float_t Gen_dR_pi_Ds;

	Float_t Gen_dR_Kp_gamma;
	Float_t Gen_dR_Kp_Dsst;
	Float_t Gen_dR_Km_gamma;
	Float_t Gen_dR_Km_Dsst;
	Float_t Gen_dR_pi_gamma;
	Float_t Gen_dR_pi_Dsst;
	Float_t Gen_dR_Ds_gamma;
	Float_t Gen_dR_Ds_Dsst;
	Float_t Gen_dR_gamma_Dsst;

	Float_t Gen_dR_Kp_mu;
	Float_t Gen_dR_Km_mu;
	Float_t Gen_dR_pi_mu;
	Float_t Gen_dR_Ds_mu;
	Float_t Gen_dR_gamma_mu;
	Float_t Gen_dR_Dsst_mu;

	Float_t Gen_met_phi;
	Float_t Gen_met_px;
	Float_t Gen_met_py;
	Float_t Gen_met_pt;
	Float_t Gen_met_sumEt;

	// trigger
	Bool_t trig_HLT_IsoMu24;
	Bool_t trig_HLT_Mu50;

	// PV
	Bool_t primvtx_isValid;
	Bool_t primvtx_isFake;
	Float_t primvtx_chi2;
	Int_t primvtx_ndof;
	Float_t primvtx_chi2ndof;
	Float_t primvtx_vx;
	Float_t primvtx_vy;
	Float_t primvtx_vz;
	Float_t primvtx_vxerr;
	Float_t primvtx_vyerr;
	Float_t primvtx_vzerr;

	// MET kinematics
	Float_t met_sig;
	Float_t met_phi;
	Float_t met_px;
	Float_t met_py;
	Float_t met_pt;
	Float_t met_sumEt;

	// Jets
	Int_t num_jets;
	Int_t num_jets_pt20;
	Int_t num_jets_pt30;
	Int_t num_jets_pt40;
	Int_t num_jets_pt50;
	Int_t num_jets_pt60;
	Int_t num_jets_pt70;
	Int_t num_jets_pt80;
	Int_t num_jets_pt90;
	Int_t num_jets_pt100;

	// muons
	Int_t num_rs_mu;
	Int_t num_ws_mu;

	// muon kinematics
	Int_t mu_charge;
	Float_t mu_eta;
	Float_t mu_phi;
	Float_t mu_px;
	Float_t mu_py;
	Float_t mu_pz;
	Float_t mu_pt;
	Float_t mu_p;
	Bool_t mu_isHighPt;
	Bool_t mu_isLoose;
	Bool_t mu_isMedium;
	Bool_t mu_isSoft;
	Bool_t mu_isTight;
	Bool_t mu_isPF;
	Bool_t mu_isTracker;
	Bool_t mu_isGlobal;

	Float_t mu_ip;
	Float_t mu_iperr;
	Float_t mu_ipsig;
	Float_t mu_ipxy;
	Float_t mu_ipxyerr;
	Float_t mu_ipxysig;

	// muon isolation
	Float_t mu_IsoR03_sumChargedHadronPt;
	Float_t mu_IsoR03_sumNeutralHadronEt;
	Float_t mu_IsoR03_sumPhotonEt;
	Float_t mu_IsoR03_sumPUPt;
	Float_t mu_IsoR03_PFIso;
	Float_t mu_IsoR03_RelIso;

	Float_t mu_IsoR04_sumChargedHadronPt;
	Float_t mu_IsoR04_sumNeutralHadronEt;
	Float_t mu_IsoR04_sumPhotonEt;
	Float_t mu_IsoR04_sumPUPt;
	Float_t mu_IsoR04_PFIso;
	Float_t mu_IsoR04_RelIso;

	Bool_t mu_match;

	// refitted info from Kinematic Fit
	Int_t Kp_charge;
	Float_t Kp_eta;
	Float_t Kp_phi;
	Float_t Kp_px;
	Float_t Kp_py;
	Float_t Kp_pz;
	Float_t Kp_pt;
	Float_t Kp_p;
	Float_t Kp_mass;
	Int_t Kp_PVassoc;

	Int_t Km_charge;
	Float_t Km_eta;
	Float_t Km_phi;
	Float_t Km_px;
	Float_t Km_py;
	Float_t Km_pz;
	Float_t Km_pt;
	Float_t Km_p;
	Float_t Km_mass;
	Int_t Km_PVassoc;

	Int_t pi_charge;
	Float_t pi_eta;
	Float_t pi_phi;
	Float_t pi_px;
	Float_t pi_py;
	Float_t pi_pz;
	Float_t pi_pt;
	Float_t pi_p;
	Float_t pi_mass;
	Int_t pi_PVassoc;

	Int_t KK_charge;
	Float_t KK_eta;
	Float_t KK_phi;
	Float_t KK_px;
	Float_t KK_py;
	Float_t KK_pz;
	Float_t KK_pt;
	Float_t KK_p;
	Float_t KK_mass;

	Int_t Kppi_charge;
	Float_t Kppi_eta;
	Float_t Kppi_phi;
	Float_t Kppi_px;
	Float_t Kppi_py;
	Float_t Kppi_pz;
	Float_t Kppi_pt;
	Float_t Kppi_p;
	Float_t Kppi_mass;

	Int_t Kmpi_charge;
	Float_t Kmpi_eta;
	Float_t Kmpi_phi;
	Float_t Kmpi_px;
	Float_t Kmpi_py;
	Float_t Kmpi_pz;
	Float_t Kmpi_pt;
	Float_t Kmpi_p;
	Float_t Kmpi_mass;

	Int_t Ds_charge;
	Float_t Ds_eta;
	Float_t Ds_phi;
	Float_t Ds_px;
	Float_t Ds_py;
	Float_t Ds_pz;
	Float_t Ds_pt;
	Float_t Ds_p;
	Float_t Ds_mass;

	Int_t Dsmu_charge;
	Float_t Dsmu_eta;
	Float_t Dsmu_phi;
	Float_t Dsmu_px;
	Float_t Dsmu_py;
	Float_t Dsmu_pz;
	Float_t Dsmu_pt;
	Float_t Dsmu_p;
	Float_t Dsmu_mass;

	Int_t gamma_charge;
	Float_t gamma_eta;
	Float_t gamma_phi;
	Float_t gamma_px;
	Float_t gamma_py;
	Float_t gamma_pz;
	Float_t gamma_pt;
	Float_t gamma_p;
	Float_t gamma_mass;
	Int_t gamma_PVassoc;

	Int_t Dsst_charge;
	Float_t Dsst_eta;
	Float_t Dsst_phi;
	Float_t Dsst_px;
	Float_t Dsst_py;
	Float_t Dsst_pz;
	Float_t Dsst_pt;
	Float_t Dsst_p;
	Float_t Dsst_mass;

	Int_t vis_charge;
	Float_t vis_eta;
	Float_t vis_phi;
	Float_t vis_px;
	Float_t vis_py;
	Float_t vis_pz;
	Float_t vis_pt;
	Float_t vis_p;
	Float_t vis_mass;

	// refitted vertex
	Float_t Ds_vx;
	Float_t Ds_vy;
	Float_t Ds_vz;
	Float_t Ds_vxerr;
	Float_t Ds_vyerr;
	Float_t Ds_vzerr;
	Float_t Ds_chi2;
	Int_t Ds_ndf;
	Float_t Ds_chi2ndf;
	Float_t Ds_prob;

	// angular distance
	Float_t dR_Kp_Km;
	Float_t dR_Kp_pi;
	Float_t dR_Kp_Ds;
	Float_t dR_Km_pi;
	Float_t dR_Km_Ds;
	Float_t dR_pi_Ds;

	Float_t dR_Kp_gamma;
	Float_t dR_Kp_Dsst;
	Float_t dR_Km_gamma;
	Float_t dR_Km_Dsst;
	Float_t dR_pi_gamma;
	Float_t dR_pi_Dsst;
	Float_t dR_Ds_gamma;
	Float_t dR_Ds_Dsst;
	Float_t dR_gamma_Dsst;

	Float_t dR_Kp_mu;
	Float_t dR_Km_mu;
	Float_t dR_pi_mu;
	Float_t dR_Ds_mu;
	Float_t dR_gamma_mu;
	Float_t dR_Dsst_mu;

	Float_t Kp_ip;
	Float_t Kp_iperr;
	Float_t Kp_ipsig;
	Float_t Kp_ipxy;
	Float_t Kp_ipxyerr;
	Float_t Kp_ipxysig;
	Float_t Km_ip;
	Float_t Km_iperr;
	Float_t Km_ipsig;
	Float_t Km_ipxy;
	Float_t Km_ipxyerr;
	Float_t Km_ipxysig;
	Float_t pi_ip;
	Float_t pi_iperr;
	Float_t pi_ipsig;
	Float_t pi_ipxy;
	Float_t pi_ipxyerr;
	Float_t pi_ipxysig;
	Float_t Ds_ip;
	Float_t Ds_iperr;
	Float_t Ds_ipsig;
	Float_t Ds_ipxy;
	Float_t Ds_ipxyerr;
	Float_t Ds_ipxysig;

	// PV-SV variable
	Float_t dira;

	Float_t Ds_FDxy;
	Float_t Ds_FDxyerr;
	Float_t Ds_FDxysig;

	Float_t Ds_FD;
	Float_t Ds_FDerr;
	Float_t Ds_FDsig;

	// Ds isolation
	Float_t Ds_IsoR03_sumChargedHadronPt;
	Float_t Ds_IsoR03_sumNeutralHadronEt;
	Float_t Ds_IsoR03_sumPhotonEt;
	Float_t Ds_IsoR03_sumPUPt;
	Float_t Ds_IsoR03_PFIso;
	Float_t Ds_IsoR03_RelIso;

	Float_t Ds_IsoR04_sumChargedHadronPt;
	Float_t Ds_IsoR04_sumNeutralHadronEt;
	Float_t Ds_IsoR04_sumPhotonEt;
	Float_t Ds_IsoR04_sumPUPt;
	Float_t Ds_IsoR04_PFIso;
	Float_t Ds_IsoR04_RelIso;

	Float_t Dsst_IsoR03_sumChargedHadronPt;
	Float_t Dsst_IsoR03_sumNeutralHadronEt;
	Float_t Dsst_IsoR03_sumPhotonEt;
	Float_t Dsst_IsoR03_sumPUPt;
	Float_t Dsst_IsoR03_PFIso;
	Float_t Dsst_IsoR03_RelIso;

	Float_t Dsst_IsoR04_sumChargedHadronPt;
	Float_t Dsst_IsoR04_sumNeutralHadronEt;
	Float_t Dsst_IsoR04_sumPhotonEt;
	Float_t Dsst_IsoR04_sumPUPt;
	Float_t Dsst_IsoR04_PFIso;
	Float_t Dsst_IsoR04_RelIso;

	Bool_t Kp_match;
	Bool_t Km_match;
	Bool_t pi_match;
	Bool_t Ds_match;
	Bool_t gamma_match;
	Bool_t Dsst_match;
	Bool_t event_match;
};
#endif

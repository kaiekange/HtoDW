#ifndef H_PHIPIMU_TREE
#define H_PHIPIMU_TREE

#include <TTree.h>

#define null -777

class H_phipimu_Tree
{
public:
	H_phipimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();
	void InitGamma();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;
	Int_t cand_id;

	// Gen PV / SV
	Float_t Gen_PV_vx;
	Float_t Gen_PV_vy;
	Float_t Gen_PV_vz;
	Float_t Gen_SV_vx;
	Float_t Gen_SV_vy;
	Float_t Gen_SV_vz;

	// Gen H
	Int_t Gen_H_charge;
	Float_t Gen_H_eta;
	Float_t Gen_H_phi;
	Float_t Gen_H_px;
	Float_t Gen_H_py;
	Float_t Gen_H_pz;
	Float_t Gen_H_pt;
	Float_t Gen_H_p;
	Float_t Gen_H_mass;

	// Gen Ds* (DsstW channel only, else null)
	Int_t Gen_Dsst_charge;
	Float_t Gen_Dsst_eta;
	Float_t Gen_Dsst_phi;
	Float_t Gen_Dsst_px;
	Float_t Gen_Dsst_py;
	Float_t Gen_Dsst_pz;
	Float_t Gen_Dsst_pt;
	Float_t Gen_Dsst_p;
	Float_t Gen_Dsst_mass;

	// Gen Ds
	Int_t Gen_Ds_charge;
	Float_t Gen_Ds_eta;
	Float_t Gen_Ds_phi;
	Float_t Gen_Ds_px;
	Float_t Gen_Ds_py;
	Float_t Gen_Ds_pz;
	Float_t Gen_Ds_pt;
	Float_t Gen_Ds_p;
	Float_t Gen_Ds_mass;

	// Gen gamma (DsstW channel only, else null)
	Int_t Gen_gamma_charge;
	Float_t Gen_gamma_eta;
	Float_t Gen_gamma_phi;
	Float_t Gen_gamma_px;
	Float_t Gen_gamma_py;
	Float_t Gen_gamma_pz;
	Float_t Gen_gamma_pt;
	Float_t Gen_gamma_p;
	Float_t Gen_gamma_mass;

	// Gen W
	Int_t Gen_W_charge;
	Float_t Gen_W_eta;
	Float_t Gen_W_phi;
	Float_t Gen_W_px;
	Float_t Gen_W_py;
	Float_t Gen_W_pz;
	Float_t Gen_W_pt;
	Float_t Gen_W_p;
	Float_t Gen_W_mass;

	// Gen K+
	Int_t Gen_Kp_charge;
	Float_t Gen_Kp_eta;
	Float_t Gen_Kp_phi;
	Float_t Gen_Kp_px;
	Float_t Gen_Kp_py;
	Float_t Gen_Kp_pz;
	Float_t Gen_Kp_pt;
	Float_t Gen_Kp_p;
	Float_t Gen_Kp_mass;

	// Gen K-
	Int_t Gen_Km_charge;
	Float_t Gen_Km_eta;
	Float_t Gen_Km_phi;
	Float_t Gen_Km_px;
	Float_t Gen_Km_py;
	Float_t Gen_Km_pz;
	Float_t Gen_Km_pt;
	Float_t Gen_Km_p;
	Float_t Gen_Km_mass;

	// Gen pi
	Int_t Gen_pi_charge;
	Float_t Gen_pi_eta;
	Float_t Gen_pi_phi;
	Float_t Gen_pi_px;
	Float_t Gen_pi_py;
	Float_t Gen_pi_pz;
	Float_t Gen_pi_pt;
	Float_t Gen_pi_p;
	Float_t Gen_pi_mass;

	// Gen mu
	Int_t Gen_mu_charge;
	Float_t Gen_mu_eta;
	Float_t Gen_mu_phi;
	Float_t Gen_mu_px;
	Float_t Gen_mu_py;
	Float_t Gen_mu_pz;
	Float_t Gen_mu_pt;
	Float_t Gen_mu_p;
	Float_t Gen_mu_mass;

	// Gen Ds+mu
	Int_t Gen_Dsmu_charge;
	Float_t Gen_Dsmu_eta;
	Float_t Gen_Dsmu_phi;
	Float_t Gen_Dsmu_px;
	Float_t Gen_Dsmu_py;
	Float_t Gen_Dsmu_pz;
	Float_t Gen_Dsmu_pt;
	Float_t Gen_Dsmu_p;
	Float_t Gen_Dsmu_mass;

	// Gen Ds*+mu (DsstW channel only, else null)
	Int_t Gen_Dsstmu_charge;
	Float_t Gen_Dsstmu_eta;
	Float_t Gen_Dsstmu_phi;
	Float_t Gen_Dsstmu_px;
	Float_t Gen_Dsstmu_py;
	Float_t Gen_Dsstmu_pz;
	Float_t Gen_Dsstmu_pt;
	Float_t Gen_Dsstmu_p;
	Float_t Gen_Dsstmu_mass;

	// Gen dR (Ds daughters)
	Float_t Gen_dR_Kp_Km;
	Float_t Gen_dR_Kp_pi;
	Float_t Gen_dR_Kp_Ds;
	Float_t Gen_dR_Km_pi;
	Float_t Gen_dR_Km_Ds;
	Float_t Gen_dR_pi_Ds;

	// Gen dR (gamma/Dsst, DsstW channel only, else null)
	Float_t Gen_dR_Kp_gamma;
	Float_t Gen_dR_Kp_Dsst;
	Float_t Gen_dR_Km_gamma;
	Float_t Gen_dR_Km_Dsst;
	Float_t Gen_dR_pi_gamma;
	Float_t Gen_dR_pi_Dsst;
	Float_t Gen_dR_Ds_gamma;
	Float_t Gen_dR_Ds_Dsst;
	Float_t Gen_dR_gamma_Dsst;

	// Gen dR (mu)
	Float_t Gen_dR_Kp_mu;
	Float_t Gen_dR_Km_mu;
	Float_t Gen_dR_pi_mu;
	Float_t Gen_dR_Ds_mu;
	Float_t Gen_dR_gamma_mu;   // DsstW channel only
	Float_t Gen_dR_Dsst_mu;    // DsstW channel only

	// Gen MET
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

	// MET
	Float_t met_sig;
	Float_t met_phi;
	Float_t met_px;
	Float_t met_py;
	Float_t met_pt;
	Float_t met_sumEt;
	Float_t W_tmass;

	// jets
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

	// muon count
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

	// refitted K+
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

	// refitted K-
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

	// refitted pi
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

	// phi = K+K-
	Int_t phi_charge;
	Float_t phi_eta;
	Float_t phi_phi;
	Float_t phi_px;
	Float_t phi_py;
	Float_t phi_pz;
	Float_t phi_pt;
	Float_t phi_p;
	Float_t phi_mass;

	// Ds
	Int_t Ds_charge;
	Float_t Ds_eta;
	Float_t Ds_phi;
	Float_t Ds_px;
	Float_t Ds_py;
	Float_t Ds_pz;
	Float_t Ds_pt;
	Float_t Ds_p;
	Float_t Ds_mass;

	// Ds + mu (always filled)
	Int_t Dsmu_charge;
	Float_t Dsmu_eta;
	Float_t Dsmu_phi;
	Float_t Dsmu_px;
	Float_t Dsmu_py;
	Float_t Dsmu_pz;
	Float_t Dsmu_pt;
	Float_t Dsmu_p;
	Float_t Dsmu_mass;
	Float_t Dsmu_Ht;

	// unconstrained gamma summary per Ds candidate (no iso/mass cut, only dM window [50,250] MeV)
	Int_t   n_gamma_window;
	Float_t gamma_uc_dM;
	Float_t gamma_uc_pt;
	Float_t gamma_uc_dR_Ds;

	// gamma (null if has_gamma == false)
	Bool_t has_gamma;
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

	// Ds* = Ds + gamma (null if has_gamma == false)
	Int_t Dsst_charge;
	Float_t Dsst_eta;
	Float_t Dsst_phi;
	Float_t Dsst_px;
	Float_t Dsst_py;
	Float_t Dsst_pz;
	Float_t Dsst_pt;
	Float_t Dsst_p;
	Float_t Dsst_mass;

	// Ds* + mu (null if has_gamma == false)
	Int_t Dsstmu_charge;
	Float_t Dsstmu_eta;
	Float_t Dsstmu_phi;
	Float_t Dsstmu_px;
	Float_t Dsstmu_py;
	Float_t Dsstmu_pz;
	Float_t Dsstmu_pt;
	Float_t Dsstmu_p;
	Float_t Dsstmu_mass;
	Float_t Dsstmu_Ht;

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

	// angular distances (Ds system)
	Float_t dR_Kp_Km;
	Float_t dR_Kp_pi;
	Float_t dR_Kp_Ds;
	Float_t dR_Km_pi;
	Float_t dR_Km_Ds;
	Float_t dR_pi_Ds;
	Float_t dR_phi_pi;
	Float_t dR_phi_Ds;
	Float_t dR_Kp_mu;
	Float_t dR_Km_mu;
	Float_t dR_pi_mu;
	Float_t dR_phi_mu;
	Float_t dR_Ds_mu;

	// angular distances (gamma/Dsst, null if has_gamma == false)
	Float_t dR_Kp_gamma;
	Float_t dR_Kp_Dsst;
	Float_t dR_Km_gamma;
	Float_t dR_Km_Dsst;
	Float_t dR_pi_gamma;
	Float_t dR_pi_Dsst;
	Float_t dR_phi_gamma;
	Float_t dR_phi_Dsst;
	Float_t dR_Ds_gamma;
	Float_t dR_Ds_Dsst;
	Float_t dR_gamma_Dsst;
	Float_t dR_gamma_mu;
	Float_t dR_Dsst_mu;

	// IP
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

	// PV-SV
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

	// Ds* isolation (null if has_gamma == false)
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

	// match flags
	Bool_t Kp_match;
	Bool_t Km_match;
	Bool_t phi_match;
	Bool_t pi_match;
	Bool_t Ds_match;
	Bool_t gamma_match;   // false if has_gamma==false or channel!="DsstW"
	Bool_t Dsst_match;    // false if has_gamma==false or channel!="DsstW"
	Bool_t event_match;
};
#endif

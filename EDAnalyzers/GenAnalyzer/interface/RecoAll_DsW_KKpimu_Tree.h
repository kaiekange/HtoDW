#ifndef RECO_DSW_KKPI_TREE
#define RECO_DSW_KKPI_TREE

#include <TTree.h>

#define null -777

class RecoAll_DsW_KKpimu_Tree
{
public:
	RecoAll_DsW_KKpimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();
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

	Int_t Gen_Ds_charge;
	Float_t Gen_Ds_eta;
	Float_t Gen_Ds_phi;
	Float_t Gen_Ds_px;
	Float_t Gen_Ds_py;
	Float_t Gen_Ds_pz;
	Float_t Gen_Ds_pt;
	Float_t Gen_Ds_p;
	Float_t Gen_Ds_mass;

	Int_t Gen_W_charge;
	Float_t Gen_W_eta;
	Float_t Gen_W_phi;
	Float_t Gen_W_px;
	Float_t Gen_W_py;
	Float_t Gen_W_pz;
	Float_t Gen_W_pt;
	Float_t Gen_W_p;
	Float_t Gen_W_mass;

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
	Float_t Gen_dR_Kp_mu;
	Float_t Gen_dR_Km_mu;
	Float_t Gen_dR_pi_mu;
	Float_t Gen_dR_Ds_mu;

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
	std::vector<Int_t> Kp_charge;
	std::vector<Float_t> Kp_eta;
	std::vector<Float_t> Kp_phi;
	std::vector<Float_t> Kp_px;
	std::vector<Float_t> Kp_py;
	std::vector<Float_t> Kp_pz;
	std::vector<Float_t> Kp_pt;
	std::vector<Float_t> Kp_p;
	std::vector<Float_t> Kp_mass;
	std::vector<Int_t> Kp_PVassoc;

	std::vector<Int_t> Km_charge;
	std::vector<Float_t> Km_eta;
	std::vector<Float_t> Km_phi;
	std::vector<Float_t> Km_px;
	std::vector<Float_t> Km_py;
	std::vector<Float_t> Km_pz;
	std::vector<Float_t> Km_pt;
	std::vector<Float_t> Km_p;
	std::vector<Float_t> Km_mass;
	std::vector<Int_t> Km_PVassoc;

	std::vector<Int_t> pi_charge;
	std::vector<Float_t> pi_eta;
	std::vector<Float_t> pi_phi;
	std::vector<Float_t> pi_px;
	std::vector<Float_t> pi_py;
	std::vector<Float_t> pi_pz;
	std::vector<Float_t> pi_pt;
	std::vector<Float_t> pi_p;
	std::vector<Float_t> pi_mass;
	std::vector<Int_t> pi_PVassoc;

	std::vector<Int_t> KK_charge;
	std::vector<Float_t> KK_eta;
	std::vector<Float_t> KK_phi;
	std::vector<Float_t> KK_px;
	std::vector<Float_t> KK_py;
	std::vector<Float_t> KK_pz;
	std::vector<Float_t> KK_pt;
	std::vector<Float_t> KK_p;
	std::vector<Float_t> KK_mass;

	std::vector<Int_t> Kppi_charge;
	std::vector<Float_t> Kppi_eta;
	std::vector<Float_t> Kppi_phi;
	std::vector<Float_t> Kppi_px;
	std::vector<Float_t> Kppi_py;
	std::vector<Float_t> Kppi_pz;
	std::vector<Float_t> Kppi_pt;
	std::vector<Float_t> Kppi_p;
	std::vector<Float_t> Kppi_mass;

	std::vector<Int_t> Kmpi_charge;
	std::vector<Float_t> Kmpi_eta;
	std::vector<Float_t> Kmpi_phi;
	std::vector<Float_t> Kmpi_px;
	std::vector<Float_t> Kmpi_py;
	std::vector<Float_t> Kmpi_pz;
	std::vector<Float_t> Kmpi_pt;
	std::vector<Float_t> Kmpi_p;
	std::vector<Float_t> Kmpi_mass;

	std::vector<Int_t> Ds_charge;
	std::vector<Float_t> Ds_eta;
	std::vector<Float_t> Ds_phi;
	std::vector<Float_t> Ds_px;
	std::vector<Float_t> Ds_py;
	std::vector<Float_t> Ds_pz;
	std::vector<Float_t> Ds_pt;
	std::vector<Float_t> Ds_p;
	std::vector<Float_t> Ds_mass;

	std::vector<Int_t> vis_charge;
	std::vector<Float_t> vis_eta;
	std::vector<Float_t> vis_phi;
	std::vector<Float_t> vis_px;
	std::vector<Float_t> vis_py;
	std::vector<Float_t> vis_pz;
	std::vector<Float_t> vis_pt;
	std::vector<Float_t> vis_p;
	std::vector<Float_t> vis_mass;

	// refitted vertex
	std::vector<Float_t> Ds_vx;
	std::vector<Float_t> Ds_vy;
	std::vector<Float_t> Ds_vz;
	std::vector<Float_t> Ds_vxerr;
	std::vector<Float_t> Ds_vyerr;
	std::vector<Float_t> Ds_vzerr;
	std::vector<Float_t> Ds_chi2;
	std::vector<Int_t> Ds_ndf;
	std::vector<Float_t> Ds_chi2ndf;
	std::vector<Float_t> Ds_prob;

	// angular distance
	std::vector<Float_t> dR_Kp_Km;
	std::vector<Float_t> dR_Kp_pi;
	std::vector<Float_t> dR_Kp_Ds;
	std::vector<Float_t> dR_Km_pi;
	std::vector<Float_t> dR_Km_Ds;
	std::vector<Float_t> dR_pi_Ds;
	std::vector<Float_t> dR_Kp_mu;
	std::vector<Float_t> dR_Km_mu;
	std::vector<Float_t> dR_pi_mu;
	std::vector<Float_t> dR_Ds_mu;

	std::vector<Float_t> Kp_ip;
	std::vector<Float_t> Kp_iperr;
	std::vector<Float_t> Kp_ipsig;
	std::vector<Float_t> Kp_ipxy;
	std::vector<Float_t> Kp_ipxyerr;
	std::vector<Float_t> Kp_ipxysig;
	std::vector<Float_t> Km_ip;
	std::vector<Float_t> Km_iperr;
	std::vector<Float_t> Km_ipsig;
	std::vector<Float_t> Km_ipxy;
	std::vector<Float_t> Km_ipxyerr;
	std::vector<Float_t> Km_ipxysig;
	std::vector<Float_t> pi_ip;
	std::vector<Float_t> pi_iperr;
	std::vector<Float_t> pi_ipsig;
	std::vector<Float_t> pi_ipxy;
	std::vector<Float_t> pi_ipxyerr;
	std::vector<Float_t> pi_ipxysig;
	std::vector<Float_t> Ds_ip;
	std::vector<Float_t> Ds_iperr;
	std::vector<Float_t> Ds_ipsig;
	std::vector<Float_t> Ds_ipxy;
	std::vector<Float_t> Ds_ipxyerr;
	std::vector<Float_t> Ds_ipxysig;

	// PV-SV variable
	std::vector<Float_t> dira;

	std::vector<Float_t> Ds_FDxy;
	std::vector<Float_t> Ds_FDxyerr;
	std::vector<Float_t> Ds_FDxysig;

	std::vector<Float_t> Ds_FD;
	std::vector<Float_t> Ds_FDerr;
	std::vector<Float_t> Ds_FDsig;

	// Ds ioslation
	std::vector<Float_t> Ds_IsoR03_sumChargedHadronPt;
	std::vector<Float_t> Ds_IsoR03_sumNeutralHadronEt;
	std::vector<Float_t> Ds_IsoR03_sumPhotonEt;
	std::vector<Float_t> Ds_IsoR03_sumPUPt;
	std::vector<Float_t> Ds_IsoR03_PFIso;
	std::vector<Float_t> Ds_IsoR03_RelIso;

	std::vector<Float_t> Ds_IsoR04_sumChargedHadronPt;
	std::vector<Float_t> Ds_IsoR04_sumNeutralHadronEt;
	std::vector<Float_t> Ds_IsoR04_sumPhotonEt;
	std::vector<Float_t> Ds_IsoR04_sumPUPt;
	std::vector<Float_t> Ds_IsoR04_PFIso;
	std::vector<Float_t> Ds_IsoR04_RelIso;

	std::vector<Bool_t> Kp_match;
	std::vector<Bool_t> Km_match;
	std::vector<Bool_t> pi_match;
	std::vector<Bool_t> Ds_match;
	std::vector<Bool_t> event_match;

	Int_t nCandidates;
};
#endif

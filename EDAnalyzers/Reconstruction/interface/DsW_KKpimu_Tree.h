#ifndef DSW_KKPI_TREE
#define DSW_KKPI_TREE

#include <TTree.h>

#define null -777

class DsW_KKpimu_Tree
{
public:
	DsW_KKpimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;

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
	Float_t dR_Kp_mu;
	Float_t dR_Km_mu;
	Float_t dR_pi_mu;
	Float_t dR_Ds_mu;

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
};
#endif

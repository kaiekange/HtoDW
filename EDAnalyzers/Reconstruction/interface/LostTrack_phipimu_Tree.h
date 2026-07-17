#ifndef LOSTTRACK_PHIPIMU_TREE
#define LOSTTRACK_PHIPIMU_TREE

#include <TTree.h>

#define null -777

class LostTrack_phipimu_Tree
{
public:
	LostTrack_phipimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;

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

	// gen match to packedPFCandidates
	Bool_t Kp_pf_match;
	Bool_t Km_pf_match;
	Bool_t pi_pf_match;
	Bool_t gamma_pf_match;   // DsstW channel only

	// gen match to lostTracks. No gamma: lostTracks holds only charged tracks
	// (|pdgId| 211/13), so a photon can never be recovered from it.
	Bool_t Kp_lt_match;
	Bool_t Km_lt_match;
	Bool_t pi_lt_match;

	// gen match surviving the charged-hadron pool cuts of H_phipimu
	// (hasTrackDetails, trackHighPurity, |pdgId|=211, pt>2 for K and pt>1 for pi)
	Bool_t Kp_pf_pool;
	Bool_t Km_pf_pool;
	Bool_t pi_pf_pool;
	Bool_t Kp_lt_pool;
	Bool_t Km_lt_pool;
	Bool_t pi_lt_pool;

	// matched lostTrack properties (null if no match)
	Float_t Kp_lt_pt;
	Float_t Kp_lt_eta;
	Float_t Km_lt_pt;
	Float_t Km_lt_eta;
	Float_t pi_lt_pt;
	Float_t pi_lt_eta;

	// gen-match quality (null if no match). dR to the gen daughter and the
	// momentum residual (pt_reco - pt_gen) / pt_gen. Both are bounded by the
	// matching criteria themselves: dR < 0.05 and |dpt| < 0.3.
	Float_t Kp_pf_dR;
	Float_t Kp_pf_dpt;
	Float_t Km_pf_dR;
	Float_t Km_pf_dpt;
	Float_t pi_pf_dR;
	Float_t pi_pf_dpt;

	Float_t Kp_lt_dR;
	Float_t Kp_lt_dpt;
	Float_t Km_lt_dR;
	Float_t Km_lt_dpt;
	Float_t pi_lt_dR;
	Float_t pi_lt_dpt;
};

#endif

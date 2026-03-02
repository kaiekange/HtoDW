#ifndef FULLRECOTREE_H
#define FULLRECOTREE_H

#include <TTree.h>

#define null -777

class FullRecoTree
{
public:
	FullRecoTree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	Bool_t primvtx_isValid;
	Bool_t primvtx_isFake;
	Float_t primvtx_chi2;
	Float_t primvtx_ndof;
	Float_t primvtx_chi2ndof;
	Float_t primvtx_vx;
	Float_t primvtx_vy;
	Float_t primvtx_vz;
	Float_t primvtx_vxerr;
	Float_t primvtx_vyerr;
	Float_t primvtx_vzerr;

	Bool_t trig_HLT_IsoMu24;
	Bool_t trig_HLT_Mu50;

	Int_t Kp_charge;
	Float_t Kp_eta;
	Float_t Kp_phi;
	Float_t Kp_px;
	Float_t Kp_py;
	Float_t Kp_pz;
	Float_t Kp_pt;
	Float_t Kp_p;
	Float_t Kp_mass;

	Int_t Km_charge;
	Float_t Km_eta;
	Float_t Km_phi;
	Float_t Km_px;
	Float_t Km_py;
	Float_t Km_pz;
	Float_t Km_pt;
	Float_t Km_p;
	Float_t Km_mass;

	Int_t pi_charge;
	Float_t pi_eta;
	Float_t pi_phi;
	Float_t pi_px;
	Float_t pi_py;
	Float_t pi_pz;
	Float_t pi_pt;
	Float_t pi_p;
	Float_t pi_mass;

	Int_t phi_charge;
	Float_t phi_eta;
	Float_t phi_phi;
	Float_t phi_px;
	Float_t phi_py;
	Float_t phi_pz;
	Float_t phi_pt;
	Float_t phi_p;
	Float_t phi_mass;

	Int_t Ds_charge;
	Float_t Ds_eta;
	Float_t Ds_phi;
	Float_t Ds_px;
	Float_t Ds_py;
	Float_t Ds_pz;
	Float_t Ds_pt;
	Float_t Ds_p;
	Float_t Ds_mass;

	Float_t Ds_vx;
	Float_t Ds_vy;
	Float_t Ds_vz;
	Float_t Ds_vxerr;
	Float_t Ds_vyerr;
	Float_t Ds_vzerr;
	Float_t Ds_chi2;
	Float_t Ds_ndf;
	Float_t Ds_chi2ndf;
	Float_t Ds_prob;

	Float_t dR_Kp_Km;
	Float_t dR_Kp_pi;
	Float_t dR_Kp_phi;
	Float_t dR_Kp_Ds;
	Float_t dR_Km_pi;
	Float_t dR_Km_phi;
	Float_t dR_Km_Ds;
	Float_t dR_pi_phi;
	Float_t dR_pi_Ds;
	Float_t dR_phi_Ds;

	Float_t Kp_ip;
	Float_t Kp_iperr;
	Float_t Kp_ipsig;
	Float_t Km_ip;
	Float_t Km_iperr;
	Float_t Km_ipsig;
	Float_t pi_ip;
	Float_t pi_iperr;
	Float_t pi_ipsig;
	Float_t Ds_ip;
	Float_t Ds_iperr;
	Float_t Ds_ipsig;

	Float_t dira;

	Float_t Ds_FDxy;
	Float_t Ds_FDxyerr;
	Float_t Ds_FDxysig;

	Float_t Ds_FD;
	Float_t Ds_FDerr;
	Float_t Ds_FDsig;

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

	Int_t num_ss_mu;
	Int_t num_os_mu;

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

	Float_t dR_Kp_mu;
	Float_t dR_Km_mu;
	Float_t dR_pi_mu;
	Float_t dR_phi_mu;
	Float_t dR_Ds_mu;

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

	Float_t met_sig;
	Float_t met_phi;
	Float_t met_px;
	Float_t met_py;
	Float_t met_pt;
	Float_t met_sumEt;
	Float_t met_corphi;
	Float_t met_corpx;
	Float_t met_corpy;
	Float_t met_corpt;
	Float_t met_corsumEt;

	Float_t num_jets;
	Float_t num_jets_pt30;
	Float_t deltaR_Ds_Jet;
};
#endif
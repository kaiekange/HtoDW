#ifndef KINETREE_H
#define KINETREE_H

#include <TTree.h>

#define null -777

class KineTree
{
public:
	KineTree(TTree *tree_);
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

	Int_t Kp_charge;
	Float_t Kp_eta;
	Float_t Kp_phi;
	Float_t Kp_px;
	Float_t Kp_py;
	Float_t Kp_pz;
	Float_t Kp_pt;
	Float_t Kp_p;

	Int_t Km_charge;
	Float_t Km_eta;
	Float_t Km_phi;
	Float_t Km_px;
	Float_t Km_py;
	Float_t Km_pz;
	Float_t Km_pt;
	Float_t Km_p;

	Int_t pi_charge;
	Float_t pi_eta;
	Float_t pi_phi;
	Float_t pi_px;
	Float_t pi_py;
	Float_t pi_pz;
	Float_t pi_pt;
	Float_t pi_p;

	Int_t phiFit_Kp_charge;
	Float_t phiFit_Kp_eta;
	Float_t phiFit_Kp_phi;
	Float_t phiFit_Kp_px;
	Float_t phiFit_Kp_py;
	Float_t phiFit_Kp_pz;
	Float_t phiFit_Kp_pt;
	Float_t phiFit_Kp_p;
	Float_t phiFit_Kp_mass;

	Int_t phiFit_Km_charge;
	Float_t phiFit_Km_eta;
	Float_t phiFit_Km_phi;
	Float_t phiFit_Km_px;
	Float_t phiFit_Km_py;
	Float_t phiFit_Km_pz;
	Float_t phiFit_Km_pt;
	Float_t phiFit_Km_p;
	Float_t phiFit_Km_mass;

	Int_t phiFit_phi_charge;
	Float_t phiFit_phi_eta;
	Float_t phiFit_phi_phi;
	Float_t phiFit_phi_px;
	Float_t phiFit_phi_py;
	Float_t phiFit_phi_pz;
	Float_t phiFit_phi_pt;
	Float_t phiFit_phi_p;
	Float_t phiFit_phi_mass;

	Int_t DsFit_Kp_charge;
	Float_t DsFit_Kp_eta;
	Float_t DsFit_Kp_phi;
	Float_t DsFit_Kp_px;
	Float_t DsFit_Kp_py;
	Float_t DsFit_Kp_pz;
	Float_t DsFit_Kp_pt;
	Float_t DsFit_Kp_p;
	Float_t DsFit_Kp_mass;

	Int_t DsFit_Km_charge;
	Float_t DsFit_Km_eta;
	Float_t DsFit_Km_phi;
	Float_t DsFit_Km_px;
	Float_t DsFit_Km_py;
	Float_t DsFit_Km_pz;
	Float_t DsFit_Km_pt;
	Float_t DsFit_Km_p;
	Float_t DsFit_Km_mass;

	Int_t DsFit_pi_charge;
	Float_t DsFit_pi_eta;
	Float_t DsFit_pi_phi;
	Float_t DsFit_pi_px;
	Float_t DsFit_pi_py;
	Float_t DsFit_pi_pz;
	Float_t DsFit_pi_pt;
	Float_t DsFit_pi_p;
	Float_t DsFit_pi_mass;

	Int_t DsFit_phi_charge;
	Float_t DsFit_phi_eta;
	Float_t DsFit_phi_phi;
	Float_t DsFit_phi_px;
	Float_t DsFit_phi_py;
	Float_t DsFit_phi_pz;
	Float_t DsFit_phi_pt;
	Float_t DsFit_phi_p;
	Float_t DsFit_phi_mass;

	Int_t DsFit_Ds_charge;
	Float_t DsFit_Ds_eta;
	Float_t DsFit_Ds_phi;
	Float_t DsFit_Ds_px;
	Float_t DsFit_Ds_py;
	Float_t DsFit_Ds_pz;
	Float_t DsFit_Ds_pt;
	Float_t DsFit_Ds_p;
	Float_t DsFit_Ds_mass;

	Float_t phi_vx;
	Float_t phi_vy;
	Float_t phi_vz;
	Float_t phi_vxerr;
	Float_t phi_vyerr;
	Float_t phi_vzerr;
	Float_t phi_chi2;
	Float_t phi_ndf;
	Float_t phi_chi2ndf;
	Float_t phi_prob;

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

	Float_t dxy_phi_Ds;
	Float_t dz_phi_Ds;

	Float_t dR_Kp_Km;
	Float_t dR_Kp_pi;
	Float_t dR_Km_pi;

	Float_t phiFit_dR_Kp_Km;
	Float_t phiFit_dR_Kp_pi;
	Float_t phiFit_dR_Kp_phi;
	Float_t phiFit_dR_Km_pi;
	Float_t phiFit_dR_Km_phi;
	Float_t phiFit_dR_pi_phi;

	Float_t DsFit_dR_Kp_Km;
	Float_t DsFit_dR_Kp_pi;
	Float_t DsFit_dR_Kp_phi;
	Float_t DsFit_dR_Kp_Ds;
	Float_t DsFit_dR_Km_pi;
	Float_t DsFit_dR_Km_phi;
	Float_t DsFit_dR_Km_Ds;
	Float_t DsFit_dR_pi_phi;
	Float_t DsFit_dR_pi_Ds;
	Float_t DsFit_dR_phi_Ds;

	Float_t Kp_ip;
	Float_t Kp_iperr;
	Float_t Kp_ipsig;
	Float_t Km_ip;
	Float_t Km_iperr;
	Float_t Km_ipsig;
	Float_t pi_ip;
	Float_t pi_iperr;
	Float_t pi_ipsig;
	Float_t phi_ip;
	Float_t phi_iperr;
	Float_t phi_ipsig;
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

};
#endif
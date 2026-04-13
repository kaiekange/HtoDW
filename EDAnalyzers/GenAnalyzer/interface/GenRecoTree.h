#ifndef GENRECOTREE_H
#define GENRECOTREE_H

#include <TTree.h>

#define null -777

class GenRecoTree
{
public:
	GenRecoTree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	Int_t num_Gen_H;
	Int_t num_Gen_Ds;
	Int_t num_Gen_W;
	Int_t num_Gen_l;
	Int_t num_Gen_Kp;
	Int_t num_Gen_Km;
	Int_t num_Gen_pi;
	Int_t MU_EVENT;

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

	Int_t Gen_num_jets;
	Int_t Gen_num_jets_pt20;
	Int_t Gen_num_jets_pt30;
	Int_t Gen_num_jets_pt40;
	Int_t Gen_num_jets_pt50;
	Int_t Gen_num_jets_pt60;
	Int_t Gen_num_jets_pt70;
	Int_t Gen_num_jets_pt80;
	Int_t Gen_num_jets_pt90;
	Int_t Gen_num_jets_pt100;

	Float_t met_phi;
	Float_t met_px;
	Float_t met_py;
	Float_t met_pt;
	Float_t met_sumEt;
};
#endif
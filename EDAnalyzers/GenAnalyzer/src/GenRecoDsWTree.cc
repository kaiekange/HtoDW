#include "EDAnalyzers/GenAnalyzer/interface/GenRecoDsWTree.h"
#include <iostream>

GenRecoDsWTree::GenRecoDsWTree(TTree *tree_)
{
	tree = tree_;
}

void GenRecoDsWTree::Init()
{
	MU_EVENT = false;
	
	Gen_PV_vx = null;
	Gen_PV_vy = null;
	Gen_PV_vz = null;

	Gen_SV_vx = null;
	Gen_SV_vy = null;
	Gen_SV_vz = null;

	Gen_H_charge = null;
	Gen_H_eta = null;
	Gen_H_phi = null;
	Gen_H_px = null;
	Gen_H_py = null;
	Gen_H_pz = null;
	Gen_H_pt = null;
	Gen_H_p = null;
	Gen_H_mass = null;

	Gen_Ds_charge = null;
	Gen_Ds_eta = null;
	Gen_Ds_phi = null;
	Gen_Ds_px = null;
	Gen_Ds_py = null;
	Gen_Ds_pz = null;
	Gen_Ds_pt = null;
	Gen_Ds_p = null;
	Gen_Ds_mass = null;

	Gen_W_charge = null;
	Gen_W_eta = null;
	Gen_W_phi = null;
	Gen_W_px = null;
	Gen_W_py = null;
	Gen_W_pz = null;
	Gen_W_pt = null;
	Gen_W_p = null;
	Gen_W_mass = null;

	Gen_Kp_charge = null;
	Gen_Kp_eta = null;
	Gen_Kp_phi = null;
	Gen_Kp_px = null;
	Gen_Kp_py = null;
	Gen_Kp_pz = null;
	Gen_Kp_pt = null;
	Gen_Kp_p = null;
	Gen_Kp_mass = null;

	Gen_Km_charge = null;
	Gen_Km_eta = null;
	Gen_Km_phi = null;
	Gen_Km_px = null;
	Gen_Km_py = null;
	Gen_Km_pz = null;
	Gen_Km_pt = null;
	Gen_Km_p = null;
	Gen_Km_mass = null;

	Gen_pi_charge = null;
	Gen_pi_eta = null;
	Gen_pi_phi = null;
	Gen_pi_px = null;
	Gen_pi_py = null;
	Gen_pi_pz = null;
	Gen_pi_pt = null;
	Gen_pi_p = null;
	Gen_pi_mass = null;

	Gen_l_charge = null;
	Gen_l_eta = null;
	Gen_l_phi = null;
	Gen_l_px = null;
	Gen_l_py = null;
	Gen_l_pz = null;
	Gen_l_pt = null;
	Gen_l_p = null;
	Gen_l_mass = null;

	Gen_vis_charge = null;
	Gen_vis_eta = null;
	Gen_vis_phi = null;
	Gen_vis_px = null;
	Gen_vis_py = null;
	Gen_vis_pz = null;
	Gen_vis_pt = null;
	Gen_vis_p = null;
	Gen_vis_mass = null;

	Gen_dR_Kp_Km = null;
	Gen_dR_Kp_pi = null;
	Gen_dR_Kp_Ds = null;
	Gen_dR_Km_pi = null;
	Gen_dR_Km_Ds = null;
	Gen_dR_pi_Ds = null;
	Gen_dR_Kp_l = null;
	Gen_dR_Km_l = null;
	Gen_dR_pi_l = null;
	Gen_dR_Ds_l = null;

	Gen_num_jets = null;
	Gen_num_jets_pt20 = null;
	Gen_num_jets_pt30 = null;
	Gen_num_jets_pt40 = null;
	Gen_num_jets_pt50 = null;
	Gen_num_jets_pt60 = null;
	Gen_num_jets_pt70 = null;
	Gen_num_jets_pt80 = null;
	Gen_num_jets_pt90 = null;
	Gen_num_jets_pt100 = null;

	met_phi = null;
	met_px = null;
	met_py = null;
	met_pt = null;
	met_sumEt = null;
}

void GenRecoDsWTree::CreateBranches()
{
	tree->Branch("MU_EVENT", &MU_EVENT);

	tree->Branch("Gen_PV_vx", &Gen_PV_vx);
	tree->Branch("Gen_PV_vy", &Gen_PV_vy);
	tree->Branch("Gen_PV_vz", &Gen_PV_vz);

	tree->Branch("Gen_SV_vx", &Gen_SV_vx);
	tree->Branch("Gen_SV_vy", &Gen_SV_vy);
	tree->Branch("Gen_SV_vz", &Gen_SV_vz);

	tree->Branch("Gen_H_charge", &Gen_H_charge);
	tree->Branch("Gen_H_eta", &Gen_H_eta);
	tree->Branch("Gen_H_phi", &Gen_H_phi);
	tree->Branch("Gen_H_px", &Gen_H_px);
	tree->Branch("Gen_H_py", &Gen_H_py);
	tree->Branch("Gen_H_pz", &Gen_H_pz);
	tree->Branch("Gen_H_pt", &Gen_H_pt);
	tree->Branch("Gen_H_p", &Gen_H_p);
	tree->Branch("Gen_H_mass", &Gen_H_mass);

	tree->Branch("Gen_Ds_charge", &Gen_Ds_charge);
	tree->Branch("Gen_Ds_eta", &Gen_Ds_eta);
	tree->Branch("Gen_Ds_phi", &Gen_Ds_phi);
	tree->Branch("Gen_Ds_px", &Gen_Ds_px);
	tree->Branch("Gen_Ds_py", &Gen_Ds_py);
	tree->Branch("Gen_Ds_pz", &Gen_Ds_pz);
	tree->Branch("Gen_Ds_pt", &Gen_Ds_pt);
	tree->Branch("Gen_Ds_p", &Gen_Ds_p);
	tree->Branch("Gen_Ds_mass", &Gen_Ds_mass);

	tree->Branch("Gen_W_charge", &Gen_W_charge);
	tree->Branch("Gen_W_eta", &Gen_W_eta);
	tree->Branch("Gen_W_phi", &Gen_W_phi);
	tree->Branch("Gen_W_px", &Gen_W_px);
	tree->Branch("Gen_W_py", &Gen_W_py);
	tree->Branch("Gen_W_pz", &Gen_W_pz);
	tree->Branch("Gen_W_pt", &Gen_W_pt);
	tree->Branch("Gen_W_p", &Gen_W_p);
	tree->Branch("Gen_W_mass", &Gen_W_mass);

	tree->Branch("Gen_Kp_charge", &Gen_Kp_charge);
	tree->Branch("Gen_Kp_eta", &Gen_Kp_eta);
	tree->Branch("Gen_Kp_phi", &Gen_Kp_phi);
	tree->Branch("Gen_Kp_px", &Gen_Kp_px);
	tree->Branch("Gen_Kp_py", &Gen_Kp_py);
	tree->Branch("Gen_Kp_pz", &Gen_Kp_pz);
	tree->Branch("Gen_Kp_pt", &Gen_Kp_pt);
	tree->Branch("Gen_Kp_p", &Gen_Kp_p);
	tree->Branch("Gen_Kp_mass", &Gen_Kp_mass);

	tree->Branch("Gen_Km_charge", &Gen_Km_charge);
	tree->Branch("Gen_Km_eta", &Gen_Km_eta);
	tree->Branch("Gen_Km_phi", &Gen_Km_phi);
	tree->Branch("Gen_Km_px", &Gen_Km_px);
	tree->Branch("Gen_Km_py", &Gen_Km_py);
	tree->Branch("Gen_Km_pz", &Gen_Km_pz);
	tree->Branch("Gen_Km_pt", &Gen_Km_pt);
	tree->Branch("Gen_Km_p", &Gen_Km_p);
	tree->Branch("Gen_Km_mass", &Gen_Km_mass);

	tree->Branch("Gen_pi_charge", &Gen_pi_charge);
	tree->Branch("Gen_pi_eta", &Gen_pi_eta);
	tree->Branch("Gen_pi_phi", &Gen_pi_phi);
	tree->Branch("Gen_pi_px", &Gen_pi_px);
	tree->Branch("Gen_pi_py", &Gen_pi_py);
	tree->Branch("Gen_pi_pz", &Gen_pi_pz);
	tree->Branch("Gen_pi_pt", &Gen_pi_pt);
	tree->Branch("Gen_pi_p", &Gen_pi_p);
	tree->Branch("Gen_pi_mass", &Gen_pi_mass);

	tree->Branch("Gen_l_charge", &Gen_l_charge);
	tree->Branch("Gen_l_eta", &Gen_l_eta);
	tree->Branch("Gen_l_phi", &Gen_l_phi);
	tree->Branch("Gen_l_px", &Gen_l_px);
	tree->Branch("Gen_l_py", &Gen_l_py);
	tree->Branch("Gen_l_pz", &Gen_l_pz);
	tree->Branch("Gen_l_pt", &Gen_l_pt);
	tree->Branch("Gen_l_p", &Gen_l_p);
	tree->Branch("Gen_l_mass", &Gen_l_mass);

	tree->Branch("Gen_vis_charge", &Gen_vis_charge);
	tree->Branch("Gen_vis_eta", &Gen_vis_eta);
	tree->Branch("Gen_vis_phi", &Gen_vis_phi);
	tree->Branch("Gen_vis_px", &Gen_vis_px);
	tree->Branch("Gen_vis_py", &Gen_vis_py);
	tree->Branch("Gen_vis_pz", &Gen_vis_pz);
	tree->Branch("Gen_vis_pt", &Gen_vis_pt);
	tree->Branch("Gen_vis_p", &Gen_vis_p);
	tree->Branch("Gen_vis_mass", &Gen_vis_mass);

	tree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
	tree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
	tree->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds);
	tree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);
	tree->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds);
	tree->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds);
	tree->Branch("Gen_dR_Kp_l", &Gen_dR_Kp_l);
	tree->Branch("Gen_dR_Km_l", &Gen_dR_Km_l);
	tree->Branch("Gen_dR_pi_l", &Gen_dR_pi_l);
	tree->Branch("Gen_dR_Ds_l", &Gen_dR_Ds_l);

	tree->Branch("Gen_num_jets", &Gen_num_jets);
	tree->Branch("Gen_num_jets_pt20", &Gen_num_jets_pt20);
	tree->Branch("Gen_num_jets_pt30", &Gen_num_jets_pt30);
	tree->Branch("Gen_num_jets_pt40", &Gen_num_jets_pt40);
	tree->Branch("Gen_num_jets_pt50", &Gen_num_jets_pt50);
	tree->Branch("Gen_num_jets_pt60", &Gen_num_jets_pt60);
	tree->Branch("Gen_num_jets_pt70", &Gen_num_jets_pt70);
	tree->Branch("Gen_num_jets_pt80", &Gen_num_jets_pt80);
	tree->Branch("Gen_num_jets_pt90", &Gen_num_jets_pt90);
	tree->Branch("Gen_num_jets_pt100", &Gen_num_jets_pt100);

	tree->Branch("met_phi", &met_phi);
	tree->Branch("met_px", &met_px);
	tree->Branch("met_py", &met_py);
	tree->Branch("met_pt", &met_pt);
	tree->Branch("met_sumEt", &met_sumEt);
}

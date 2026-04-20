#include "EDAnalyzers/GenAnalyzer/interface/Reco_DsW_KKpimu_Tree.h"
#include <iostream>

Reco_DsW_KKpimu_Tree::Reco_DsW_KKpimu_Tree(TTree *tree_)
{
	tree = tree_;
}

void Reco_DsW_KKpimu_Tree::Init()
{
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

	Gen_mu_charge = null;
	Gen_mu_eta = null;
	Gen_mu_phi = null;
	Gen_mu_px = null;
	Gen_mu_py = null;
	Gen_mu_pz = null;
	Gen_mu_pt = null;
	Gen_mu_p = null;
	Gen_mu_mass = null;

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
	Gen_dR_Kp_mu = null;
	Gen_dR_Km_mu = null;
	Gen_dR_pi_mu = null;
	Gen_dR_Ds_mu = null;

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

	Gen_met_phi = null;
	Gen_met_px = null;
	Gen_met_py = null;
	Gen_met_pt = null;
	Gen_met_sumEt = null;

	trig_HLT_IsoMu24 = false;
	trig_HLT_Mu50 = false;

	primvtx_isValid = false;
	primvtx_isFake = false;
	primvtx_chi2 = null;
	primvtx_ndof = null;
	primvtx_chi2ndof = null;
	primvtx_vx = null;
	primvtx_vy = null;
	primvtx_vz = null;
	primvtx_vxerr = null;
	primvtx_vyerr = null;
	primvtx_vzerr = null;

	// refitted info from Kinematic Fit
	Kp_charge = null;
	Kp_eta = null;
	Kp_phi = null;
	Kp_px = null;
	Kp_py = null;
	Kp_pz = null;
	Kp_pt = null;
	Kp_p = null;
	Kp_mass = null;

	Km_charge = null;
	Km_eta = null;
	Km_phi = null;
	Km_px = null;
	Km_py = null;
	Km_pz = null;
	Km_pt = null;
	Km_p = null;
	Km_mass = null;

	pi_charge = null;
	pi_eta = null;
	pi_phi = null;
	pi_px = null;
	pi_py = null;
	pi_pz = null;
	pi_pt = null;
	pi_p = null;
	pi_mass = null;

	KK_charge = null;
	KK_eta = null;
	KK_phi = null;
	KK_px = null;
	KK_py = null;
	KK_pz = null;
	KK_pt = null;
	KK_p = null;
	KK_mass = null;

	Kppi_charge = null;
	Kppi_eta = null;
	Kppi_phi = null;
	Kppi_px = null;
	Kppi_py = null;
	Kppi_pz = null;
	Kppi_pt = null;
	Kppi_p = null;
	Kppi_mass = null;

	Kmpi_charge = null;
	Kmpi_eta = null;
	Kmpi_phi = null;
	Kmpi_px = null;
	Kmpi_py = null;
	Kmpi_pz = null;
	Kmpi_pt = null;
	Kmpi_p = null;
	Kmpi_mass = null;

	Ds_charge = null;
	Ds_eta = null;
	Ds_phi = null;
	Ds_px = null;
	Ds_py = null;
	Ds_pz = null;
	Ds_pt = null;
	Ds_p = null;
	Ds_mass = null;

	vis_charge = null;
	vis_eta = null;
	vis_phi = null;
	vis_px = null;
	vis_py = null;
	vis_pz = null;
	vis_pt = null;
	vis_p = null;
	vis_mass = null;

	// refitted vertex
	Ds_vx = null;
	Ds_vy = null;
	Ds_vz = null;
	Ds_vxerr = null;
	Ds_vyerr = null;
	Ds_vzerr = null;
	Ds_chi2 = null;
	Ds_ndf = null;
	Ds_chi2ndf = null;
	Ds_prob = null;

	// angular distance
	dR_Kp_Km = null;
	dR_Kp_pi = null;
	dR_Kp_Ds = null;
	dR_Km_pi = null;
	dR_Km_Ds = null;
	dR_pi_Ds = null;

	Kp_ip = null;
	Kp_iperr = null;
	Kp_ipsig = null;
	Km_ip = null;
	Km_iperr = null;
	Km_ipsig = null;
	pi_ip = null;
	pi_iperr = null;
	pi_ipsig = null;
	Ds_ip = null;
	Ds_iperr = null;
	Ds_ipsig = null;

	// PV-SV variable
	dira = null;

	Ds_FDxy = null;
	Ds_FDxyerr = null;
	Ds_FDxysig = null;

	Ds_FD = null;
	Ds_FDerr = null;
	Ds_FDsig = null;

	// Ds ioslation
	Ds_IsoR03_sumChargedHadronPt = null;
	Ds_IsoR03_sumNeutralHadronEt = null;
	Ds_IsoR03_sumPhotonEt = null;
	Ds_IsoR03_sumPUPt = null;
	Ds_IsoR03_PFIso = null;
	Ds_IsoR03_RelIso = null;

	Ds_IsoR04_sumChargedHadronPt = null;
	Ds_IsoR04_sumNeutralHadronEt = null;
	Ds_IsoR04_sumPhotonEt = null;
	Ds_IsoR04_sumPUPt = null;
	Ds_IsoR04_PFIso = null;
	Ds_IsoR04_RelIso = null;

	// muon kinematics
	mu_charge = null;
	mu_eta = null;
	mu_phi = null;
	mu_px = null;
	mu_py = null;
	mu_pz = null;
	mu_pt = null;
	mu_p = null;
	mu_isHighPt = false;
	mu_isLoose = false;
	mu_isMedium = false;
	mu_isSoft = false;
	mu_isTight = false;
	mu_isPF = false;
	mu_isTracker = false;
	mu_isGlobal = false;

	mu_ip = null;
	mu_iperr = null;
	mu_ipsig = null;

	num_rs_mu = null;
	num_ws_mu = null;

	dR_Kp_mu = null;
	dR_Km_mu = null;
	dR_pi_mu = null;
	dR_Ds_mu = null;

	// muon isolation
	mu_IsoR03_sumChargedHadronPt = null;
	mu_IsoR03_sumNeutralHadronEt = null;
	mu_IsoR03_sumPhotonEt = null;
	mu_IsoR03_sumPUPt = null;
	mu_IsoR03_PFIso = null;
	mu_IsoR03_RelIso = null;

	mu_IsoR04_sumChargedHadronPt = null;
	mu_IsoR04_sumNeutralHadronEt = null;
	mu_IsoR04_sumPhotonEt = null;
	mu_IsoR04_sumPUPt = null;
	mu_IsoR04_PFIso = null;
	mu_IsoR04_RelIso = null;

	// MET kinematics
	met_sig = null;
	met_phi = null;
	met_px = null;
	met_py = null;
	met_pt = null;
	met_sumEt = null;

	// Jets
	num_jets = null;
	num_jets_pt20 = null;
	num_jets_pt30 = null;
	num_jets_pt40 = null;
	num_jets_pt50 = null;
	num_jets_pt60 = null;
	num_jets_pt70 = null;
	num_jets_pt80 = null;
	num_jets_pt90 = null;
	num_jets_pt100 = null;

	Kp_match = false;
	Km_match = false;
	pi_match = false;
	mu_match = false;
	event_match = false;
}

void Reco_DsW_KKpimu_Tree::CreateBranches()
{
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

	tree->Branch("Gen_mu_charge", &Gen_mu_charge);
	tree->Branch("Gen_mu_eta", &Gen_mu_eta);
	tree->Branch("Gen_mu_phi", &Gen_mu_phi);
	tree->Branch("Gen_mu_px", &Gen_mu_px);
	tree->Branch("Gen_mu_py", &Gen_mu_py);
	tree->Branch("Gen_mu_pz", &Gen_mu_pz);
	tree->Branch("Gen_mu_pt", &Gen_mu_pt);
	tree->Branch("Gen_mu_p", &Gen_mu_p);
	tree->Branch("Gen_mu_mass", &Gen_mu_mass);

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
	tree->Branch("Gen_dR_Kp_mu", &Gen_dR_Kp_mu);
	tree->Branch("Gen_dR_Km_mu", &Gen_dR_Km_mu);
	tree->Branch("Gen_dR_pi_mu", &Gen_dR_pi_mu);
	tree->Branch("Gen_dR_Ds_mu", &Gen_dR_Ds_mu);

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

	tree->Branch("Gen_met_phi", &Gen_met_phi);
	tree->Branch("Gen_met_px", &Gen_met_px);
	tree->Branch("Gen_met_py", &Gen_met_py);
	tree->Branch("Gen_met_pt", &Gen_met_pt);
	tree->Branch("Gen_met_sumEt", &Gen_met_sumEt);

	tree->Branch("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	tree->Branch("trig_HLT_Mu50", &trig_HLT_Mu50);

	tree->Branch("primvtx_isValid", &primvtx_isValid);
	tree->Branch("primvtx_isFake", &primvtx_isFake);
	tree->Branch("primvtx_chi2", &primvtx_chi2);
	tree->Branch("primvtx_ndof", &primvtx_ndof);
	tree->Branch("primvtx_chi2ndof", &primvtx_chi2ndof);
	tree->Branch("primvtx_vx", &primvtx_vx);
	tree->Branch("primvtx_vy", &primvtx_vy);
	tree->Branch("primvtx_vz", &primvtx_vz);
	tree->Branch("primvtx_vxerr", &primvtx_vxerr);
	tree->Branch("primvtx_vyerr", &primvtx_vyerr);
	tree->Branch("primvtx_vzerr", &primvtx_vzerr);

	// refitted info from Kinematic Fit
	tree->Branch("Kp_charge", &Kp_charge);
	tree->Branch("Kp_eta", &Kp_eta);
	tree->Branch("Kp_phi", &Kp_phi);
	tree->Branch("Kp_px", &Kp_px);
	tree->Branch("Kp_py", &Kp_py);
	tree->Branch("Kp_pz", &Kp_pz);
	tree->Branch("Kp_pt", &Kp_pt);
	tree->Branch("Kp_p", &Kp_p);
	tree->Branch("Kp_mass", &Kp_mass);

	tree->Branch("Km_charge", &Km_charge);
	tree->Branch("Km_eta", &Km_eta);
	tree->Branch("Km_phi", &Km_phi);
	tree->Branch("Km_px", &Km_px);
	tree->Branch("Km_py", &Km_py);
	tree->Branch("Km_pz", &Km_pz);
	tree->Branch("Km_pt", &Km_pt);
	tree->Branch("Km_p", &Km_p);
	tree->Branch("Km_mass", &Km_mass);

	tree->Branch("pi_charge", &pi_charge);
	tree->Branch("pi_eta", &pi_eta);
	tree->Branch("pi_phi", &pi_phi);
	tree->Branch("pi_px", &pi_px);
	tree->Branch("pi_py", &pi_py);
	tree->Branch("pi_pz", &pi_pz);
	tree->Branch("pi_pt", &pi_pt);
	tree->Branch("pi_p", &pi_p);
	tree->Branch("pi_mass", &pi_mass);

	tree->Branch("KK_charge", &KK_charge);
	tree->Branch("KK_eta", &KK_eta);
	tree->Branch("KK_phi", &KK_phi);
	tree->Branch("KK_px", &KK_px);
	tree->Branch("KK_py", &KK_py);
	tree->Branch("KK_pz", &KK_pz);
	tree->Branch("KK_pt", &KK_pt);
	tree->Branch("KK_p", &KK_p);
	tree->Branch("KK_mass", &KK_mass);

	tree->Branch("Kppi_charge", &Kppi_charge);
	tree->Branch("Kppi_eta", &Kppi_eta);
	tree->Branch("Kppi_phi", &Kppi_phi);
	tree->Branch("Kppi_px", &Kppi_px);
	tree->Branch("Kppi_py", &Kppi_py);
	tree->Branch("Kppi_pz", &Kppi_pz);
	tree->Branch("Kppi_pt", &Kppi_pt);
	tree->Branch("Kppi_p", &Kppi_p);
	tree->Branch("Kppi_mass", &Kppi_mass);

	tree->Branch("Kmpi_charge", &Kmpi_charge);
	tree->Branch("Kmpi_eta", &Kmpi_eta);
	tree->Branch("Kmpi_phi", &Kmpi_phi);
	tree->Branch("Kmpi_px", &Kmpi_px);
	tree->Branch("Kmpi_py", &Kmpi_py);
	tree->Branch("Kmpi_pz", &Kmpi_pz);
	tree->Branch("Kmpi_pt", &Kmpi_pt);
	tree->Branch("Kmpi_p", &Kmpi_p);
	tree->Branch("Kmpi_mass", &Kmpi_mass);

	tree->Branch("Ds_charge", &Ds_charge);
	tree->Branch("Ds_eta", &Ds_eta);
	tree->Branch("Ds_phi", &Ds_phi);
	tree->Branch("Ds_px", &Ds_px);
	tree->Branch("Ds_py", &Ds_py);
	tree->Branch("Ds_pz", &Ds_pz);
	tree->Branch("Ds_pt", &Ds_pt);
	tree->Branch("Ds_p", &Ds_p);
	tree->Branch("Ds_mass", &Ds_mass);

	tree->Branch("vis_charge", &vis_charge);
	tree->Branch("vis_eta", &vis_eta);
	tree->Branch("vis_phi", &vis_phi);
	tree->Branch("vis_px", &vis_px);
	tree->Branch("vis_py", &vis_py);
	tree->Branch("vis_pz", &vis_pz);
	tree->Branch("vis_pt", &vis_pt);
	tree->Branch("vis_p", &vis_p);
	tree->Branch("vis_mass", &vis_mass);

	// refitted vertex
	tree->Branch("Ds_vx", &Ds_vx);
	tree->Branch("Ds_vy", &Ds_vy);
	tree->Branch("Ds_vz", &Ds_vz);
	tree->Branch("Ds_vxerr", &Ds_vxerr);
	tree->Branch("Ds_vyerr", &Ds_vyerr);
	tree->Branch("Ds_vzerr", &Ds_vzerr);
	tree->Branch("Ds_chi2", &Ds_chi2);
	tree->Branch("Ds_ndf", &Ds_ndf);
	tree->Branch("Ds_chi2ndf", &Ds_chi2ndf);
	tree->Branch("Ds_prob", &Ds_prob);

	// angular distance
	tree->Branch("dR_Kp_Km", &dR_Kp_Km);
	tree->Branch("dR_Kp_pi", &dR_Kp_pi);
	tree->Branch("dR_Kp_Ds", &dR_Kp_Ds);
	tree->Branch("dR_Km_pi", &dR_Km_pi);
	tree->Branch("dR_Km_Ds", &dR_Km_Ds);
	tree->Branch("dR_pi_Ds", &dR_pi_Ds);

	tree->Branch("Kp_ip", &Kp_ip);
	tree->Branch("Kp_iperr", &Kp_iperr);
	tree->Branch("Kp_ipsig", &Kp_ipsig);
	tree->Branch("Km_ip", &Km_ip);
	tree->Branch("Km_iperr", &Km_iperr);
	tree->Branch("Km_ipsig", &Km_ipsig);
	tree->Branch("pi_ip", &pi_ip);
	tree->Branch("pi_iperr", &pi_iperr);
	tree->Branch("pi_ipsig", &pi_ipsig);
	tree->Branch("Ds_ip", &Ds_ip);
	tree->Branch("Ds_iperr", &Ds_iperr);
	tree->Branch("Ds_ipsig", &Ds_ipsig);

	// PV-SV variable
	tree->Branch("dira", &dira);

	tree->Branch("Ds_FDxy", &Ds_FDxy);
	tree->Branch("Ds_FDxyerr", &Ds_FDxyerr);
	tree->Branch("Ds_FDxysig", &Ds_FDxysig);

	tree->Branch("Ds_FD", &Ds_FD);
	tree->Branch("Ds_FDerr", &Ds_FDerr);
	tree->Branch("Ds_FDsig", &Ds_FDsig);

	// Ds ioslation
	tree->Branch("Ds_IsoR03_sumChargedHadronPt", &Ds_IsoR03_sumChargedHadronPt);
	tree->Branch("Ds_IsoR03_sumNeutralHadronEt", &Ds_IsoR03_sumNeutralHadronEt);
	tree->Branch("Ds_IsoR03_sumPhotonEt", &Ds_IsoR03_sumPhotonEt);
	tree->Branch("Ds_IsoR03_sumPUPt", &Ds_IsoR03_sumPUPt);
	tree->Branch("Ds_IsoR03_PFIso", &Ds_IsoR03_PFIso);
	tree->Branch("Ds_IsoR03_RelIso", &Ds_IsoR03_RelIso);

	tree->Branch("Ds_IsoR04_sumChargedHadronPt", &Ds_IsoR04_sumChargedHadronPt);
	tree->Branch("Ds_IsoR04_sumNeutralHadronEt", &Ds_IsoR04_sumNeutralHadronEt);
	tree->Branch("Ds_IsoR04_sumPhotonEt", &Ds_IsoR04_sumPhotonEt);
	tree->Branch("Ds_IsoR04_sumPUPt", &Ds_IsoR04_sumPUPt);
	tree->Branch("Ds_IsoR04_PFIso", &Ds_IsoR04_PFIso);
	tree->Branch("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso);

	// muon kinematics
	tree->Branch("mu_charge", &mu_charge);
	tree->Branch("mu_eta", &mu_eta);
	tree->Branch("mu_phi", &mu_phi);
	tree->Branch("mu_px", &mu_px);
	tree->Branch("mu_py", &mu_py);
	tree->Branch("mu_pz", &mu_pz);
	tree->Branch("mu_pt", &mu_pt);
	tree->Branch("mu_p", &mu_p);
	tree->Branch("mu_isHighPt", &mu_isHighPt);
	tree->Branch("mu_isLoose", &mu_isLoose);
	tree->Branch("mu_isMedium", &mu_isMedium);
	tree->Branch("mu_isSoft", &mu_isSoft);
	tree->Branch("mu_isTight", &mu_isTight);
	tree->Branch("mu_isPF", &mu_isPF);
	tree->Branch("mu_isTracker", &mu_isTracker);
	tree->Branch("mu_isGlobal", &mu_isGlobal);

	tree->Branch("mu_ip", &mu_ip);
	tree->Branch("mu_iperr", &mu_iperr);
	tree->Branch("mu_ipsig", &mu_ipsig);

	tree->Branch("dR_Kp_mu", &dR_Kp_mu);
	tree->Branch("dR_Km_mu", &dR_Km_mu);
	tree->Branch("dR_pi_mu", &dR_pi_mu);
	tree->Branch("dR_Ds_mu", &dR_Ds_mu);

	tree->Branch("num_rs_mu", &num_rs_mu);
	tree->Branch("num_ws_mu", &num_ws_mu);

	// muon isolation
	tree->Branch("mu_IsoR03_sumChargedHadronPt", &mu_IsoR03_sumChargedHadronPt);
	tree->Branch("mu_IsoR03_sumNeutralHadronEt", &mu_IsoR03_sumNeutralHadronEt);
	tree->Branch("mu_IsoR03_sumPhotonEt", &mu_IsoR03_sumPhotonEt);
	tree->Branch("mu_IsoR03_sumPUPt", &mu_IsoR03_sumPUPt);
	tree->Branch("mu_IsoR03_PFIso", &mu_IsoR03_PFIso);
	tree->Branch("mu_IsoR03_RelIso", &mu_IsoR03_RelIso);

	tree->Branch("mu_IsoR04_sumChargedHadronPt", &mu_IsoR04_sumChargedHadronPt);
	tree->Branch("mu_IsoR04_sumNeutralHadronEt", &mu_IsoR04_sumNeutralHadronEt);
	tree->Branch("mu_IsoR04_sumPhotonEt", &mu_IsoR04_sumPhotonEt);
	tree->Branch("mu_IsoR04_sumPUPt", &mu_IsoR04_sumPUPt);
	tree->Branch("mu_IsoR04_PFIso", &mu_IsoR04_PFIso);
	tree->Branch("mu_IsoR04_RelIso", &mu_IsoR04_RelIso);

	// MET kinematics
	tree->Branch("met_sig", &met_sig);
	tree->Branch("met_phi", &met_phi);
	tree->Branch("met_px", &met_px);
	tree->Branch("met_py", &met_py);
	tree->Branch("met_pt", &met_pt);
	tree->Branch("met_sumEt", &met_sumEt);

	// Jets
	tree->Branch("num_jets", &num_jets);
	tree->Branch("num_jets_pt20", &num_jets_pt20);
	tree->Branch("num_jets_pt30", &num_jets_pt30);
	tree->Branch("num_jets_pt40", &num_jets_pt40);
	tree->Branch("num_jets_pt50", &num_jets_pt50);
	tree->Branch("num_jets_pt60", &num_jets_pt60);
	tree->Branch("num_jets_pt70", &num_jets_pt70);
	tree->Branch("num_jets_pt80", &num_jets_pt80);
	tree->Branch("num_jets_pt90", &num_jets_pt90);
	tree->Branch("num_jets_pt100", &num_jets_pt100);

	tree->Branch("Kp_match", &Kp_match);
	tree->Branch("Km_match", &Km_match);
	tree->Branch("pi_match", &pi_match);
	tree->Branch("mu_match", &mu_match);
	tree->Branch("event_match", &event_match);
}

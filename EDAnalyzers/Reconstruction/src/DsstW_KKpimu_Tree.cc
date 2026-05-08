#include "EDAnalyzers/Reconstruction/interface/DsstW_KKpimu_Tree.h"
#include <iostream>

DsstW_KKpimu_Tree::DsstW_KKpimu_Tree(TTree *tree_)
{
	tree = tree_;
}

void DsstW_KKpimu_Tree::Init()
{
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

	// MET kinematics
	met_sig = null;
	met_phi = null;
	met_px = null;
	met_py = null;
	met_pt = null;
	met_sumEt = null;

	// Jets
	num_jets = 0;
	num_jets_pt20 = 0;
	num_jets_pt30 = 0;
	num_jets_pt40 = 0;
	num_jets_pt50 = 0;
	num_jets_pt60 = 0;
	num_jets_pt70 = 0;
	num_jets_pt80 = 0;
	num_jets_pt90 = 0;
	num_jets_pt100 = 0;

	// muon kinematics
	num_rs_mu = 0;
	num_ws_mu = 0;

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
	mu_ipxy = null;
	mu_ipxyerr = null;
	mu_ipxysig = null;

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

	event_id = 0;
	run_id = 0;
	lumi_id = 0;

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
	Kp_PVassoc = null;

	Km_charge = null;
	Km_eta = null;
	Km_phi = null;
	Km_px = null;
	Km_py = null;
	Km_pz = null;
	Km_pt = null;
	Km_p = null;
	Km_mass = null;
	Km_PVassoc = null;

	pi_charge = null;
	pi_eta = null;
	pi_phi = null;
	pi_px = null;
	pi_py = null;
	pi_pz = null;
	pi_pt = null;
	pi_p = null;
	pi_mass = null;
	pi_PVassoc = null;

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

	Dsmu_charge = null;
	Dsmu_eta = null;
	Dsmu_phi = null;
	Dsmu_px = null;
	Dsmu_py = null;
	Dsmu_pz = null;
	Dsmu_pt = null;
	Dsmu_p = null;
	Dsmu_mass = null;

	gamma_charge = null;
	gamma_eta = null;
	gamma_phi = null;
	gamma_px = null;
	gamma_py = null;
	gamma_pz = null;
	gamma_pt = null;
	gamma_p = null;
	gamma_mass = null;
	gamma_PVassoc = null;

	Dsst_charge = null;
	Dsst_eta = null;
	Dsst_phi = null;
	Dsst_px = null;
	Dsst_py = null;
	Dsst_pz = null;
	Dsst_pt = null;
	Dsst_p = null;
	Dsst_mass = null;

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

	dR_Kp_gamma = null;
	dR_Kp_Dsst = null;
	dR_Km_gamma = null;
	dR_Km_Dsst = null;
	dR_pi_gamma = null;
	dR_pi_Dsst = null;
	dR_Ds_gamma = null;
	dR_Ds_Dsst = null;
	dR_gamma_Dsst = null;

	dR_Kp_mu = null;
	dR_Km_mu = null;
	dR_pi_mu = null;
	dR_Ds_mu = null;
	dR_gamma_mu = null;
	dR_Dsst_mu = null;

	Kp_ip = null;
	Kp_iperr = null;
	Kp_ipsig = null;
	Kp_ipxy = null;
	Kp_ipxyerr = null;
	Kp_ipxysig = null;
	Km_ip = null;
	Km_iperr = null;
	Km_ipsig = null;
	Km_ipxy = null;
	Km_ipxyerr = null;
	Km_ipxysig = null;
	pi_ip = null;
	pi_iperr = null;
	pi_ipsig = null;
	pi_ipxy = null;
	pi_ipxyerr = null;
	pi_ipxysig = null;
	Ds_ip = null;
	Ds_iperr = null;
	Ds_ipsig = null;
	Ds_ipxy = null;
	Ds_ipxyerr = null;
	Ds_ipxysig = null;

	// PV-SV variable
	dira = null;

	Ds_FDxy = null;
	Ds_FDxyerr = null;
	Ds_FDxysig = null;

	Ds_FD = null;
	Ds_FDerr = null;
	Ds_FDsig = null;

	// Ds isolation
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

	// Dsst isolation
	Dsst_IsoR03_sumChargedHadronPt = null;
	Dsst_IsoR03_sumNeutralHadronEt = null;
	Dsst_IsoR03_sumPhotonEt = null;
	Dsst_IsoR03_sumPUPt = null;
	Dsst_IsoR03_PFIso = null;
	Dsst_IsoR03_RelIso = null;

	Dsst_IsoR04_sumChargedHadronPt = null;
	Dsst_IsoR04_sumNeutralHadronEt = null;
	Dsst_IsoR04_sumPhotonEt = null;
	Dsst_IsoR04_sumPUPt = null;
	Dsst_IsoR04_PFIso = null;
	Dsst_IsoR04_RelIso = null;
}

void DsstW_KKpimu_Tree::CreateBranches()
{
	tree->Branch("event_id", &event_id);
	tree->Branch("run_id", &run_id);
	tree->Branch("lumi_id", &lumi_id);

	// trigger
	tree->Branch("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	tree->Branch("trig_HLT_Mu50", &trig_HLT_Mu50);

	// PV
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

	// muons
	tree->Branch("num_rs_mu", &num_rs_mu);
	tree->Branch("num_ws_mu", &num_ws_mu);

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
	tree->Branch("mu_ipxy", &mu_ipxy);
	tree->Branch("mu_ipxyerr", &mu_ipxyerr);
	tree->Branch("mu_ipxysig", &mu_ipxysig);

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
	tree->Branch("Kp_PVassoc", &Kp_PVassoc);

	tree->Branch("Km_charge", &Km_charge);
	tree->Branch("Km_eta", &Km_eta);
	tree->Branch("Km_phi", &Km_phi);
	tree->Branch("Km_px", &Km_px);
	tree->Branch("Km_py", &Km_py);
	tree->Branch("Km_pz", &Km_pz);
	tree->Branch("Km_pt", &Km_pt);
	tree->Branch("Km_p", &Km_p);
	tree->Branch("Km_mass", &Km_mass);
	tree->Branch("Km_PVassoc", &Km_PVassoc);

	tree->Branch("pi_charge", &pi_charge);
	tree->Branch("pi_eta", &pi_eta);
	tree->Branch("pi_phi", &pi_phi);
	tree->Branch("pi_px", &pi_px);
	tree->Branch("pi_py", &pi_py);
	tree->Branch("pi_pz", &pi_pz);
	tree->Branch("pi_pt", &pi_pt);
	tree->Branch("pi_p", &pi_p);
	tree->Branch("pi_mass", &pi_mass);
	tree->Branch("pi_PVassoc", &pi_PVassoc);

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

	tree->Branch("Dsmu_charge", &Dsmu_charge);
	tree->Branch("Dsmu_eta", &Dsmu_eta);
	tree->Branch("Dsmu_phi", &Dsmu_phi);
	tree->Branch("Dsmu_px", &Dsmu_px);
	tree->Branch("Dsmu_py", &Dsmu_py);
	tree->Branch("Dsmu_pz", &Dsmu_pz);
	tree->Branch("Dsmu_pt", &Dsmu_pt);
	tree->Branch("Dsmu_p", &Dsmu_p);
	tree->Branch("Dsmu_mass", &Dsmu_mass);

	tree->Branch("gamma_charge", &gamma_charge);
	tree->Branch("gamma_eta", &gamma_eta);
	tree->Branch("gamma_phi", &gamma_phi);
	tree->Branch("gamma_px", &gamma_px);
	tree->Branch("gamma_py", &gamma_py);
	tree->Branch("gamma_pz", &gamma_pz);
	tree->Branch("gamma_pt", &gamma_pt);
	tree->Branch("gamma_p", &gamma_p);
	tree->Branch("gamma_mass", &gamma_mass);
	tree->Branch("gamma_PVassoc", &gamma_PVassoc);

	tree->Branch("Dsst_charge", &Dsst_charge);
	tree->Branch("Dsst_eta", &Dsst_eta);
	tree->Branch("Dsst_phi", &Dsst_phi);
	tree->Branch("Dsst_px", &Dsst_px);
	tree->Branch("Dsst_py", &Dsst_py);
	tree->Branch("Dsst_pz", &Dsst_pz);
	tree->Branch("Dsst_pt", &Dsst_pt);
	tree->Branch("Dsst_p", &Dsst_p);
	tree->Branch("Dsst_mass", &Dsst_mass);

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

	tree->Branch("dR_Kp_gamma", &dR_Kp_gamma);
	tree->Branch("dR_Kp_Dsst", &dR_Kp_Dsst);
	tree->Branch("dR_Km_gamma", &dR_Km_gamma);
	tree->Branch("dR_Km_Dsst", &dR_Km_Dsst);
	tree->Branch("dR_pi_gamma", &dR_pi_gamma);
	tree->Branch("dR_pi_Dsst", &dR_pi_Dsst);
	tree->Branch("dR_Ds_gamma", &dR_Ds_gamma);
	tree->Branch("dR_Ds_Dsst", &dR_Ds_Dsst);
	tree->Branch("dR_gamma_Dsst", &dR_gamma_Dsst);

	tree->Branch("dR_Kp_mu", &dR_Kp_mu);
	tree->Branch("dR_Km_mu", &dR_Km_mu);
	tree->Branch("dR_pi_mu", &dR_pi_mu);
	tree->Branch("dR_Ds_mu", &dR_Ds_mu);
	tree->Branch("dR_gamma_mu", &dR_gamma_mu);
	tree->Branch("dR_Dsst_mu", &dR_Dsst_mu);

	tree->Branch("Kp_ip", &Kp_ip);
	tree->Branch("Kp_iperr", &Kp_iperr);
	tree->Branch("Kp_ipsig", &Kp_ipsig);
	tree->Branch("Kp_ipxy", &Kp_ipxy);
	tree->Branch("Kp_ipxyerr", &Kp_ipxyerr);
	tree->Branch("Kp_ipxysig", &Kp_ipxysig);
	tree->Branch("Km_ip", &Km_ip);
	tree->Branch("Km_iperr", &Km_iperr);
	tree->Branch("Km_ipsig", &Km_ipsig);
	tree->Branch("Km_ipxy", &Km_ipxy);
	tree->Branch("Km_ipxyerr", &Km_ipxyerr);
	tree->Branch("Km_ipxysig", &Km_ipxysig);
	tree->Branch("pi_ip", &pi_ip);
	tree->Branch("pi_iperr", &pi_iperr);
	tree->Branch("pi_ipsig", &pi_ipsig);
	tree->Branch("pi_ipxy", &pi_ipxy);
	tree->Branch("pi_ipxyerr", &pi_ipxyerr);
	tree->Branch("pi_ipxysig", &pi_ipxysig);
	tree->Branch("Ds_ip", &Ds_ip);
	tree->Branch("Ds_iperr", &Ds_iperr);
	tree->Branch("Ds_ipsig", &Ds_ipsig);
	tree->Branch("Ds_ipxy", &Ds_ipxy);
	tree->Branch("Ds_ipxyerr", &Ds_ipxyerr);
	tree->Branch("Ds_ipxysig", &Ds_ipxysig);

	// PV-SV variable
	tree->Branch("dira", &dira);

	tree->Branch("Ds_FDxy", &Ds_FDxy);
	tree->Branch("Ds_FDxyerr", &Ds_FDxyerr);
	tree->Branch("Ds_FDxysig", &Ds_FDxysig);

	tree->Branch("Ds_FD", &Ds_FD);
	tree->Branch("Ds_FDerr", &Ds_FDerr);
	tree->Branch("Ds_FDsig", &Ds_FDsig);

	// Ds isolation
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

	// Dsst isolation
	tree->Branch("Dsst_IsoR03_sumChargedHadronPt", &Dsst_IsoR03_sumChargedHadronPt);
	tree->Branch("Dsst_IsoR03_sumNeutralHadronEt", &Dsst_IsoR03_sumNeutralHadronEt);
	tree->Branch("Dsst_IsoR03_sumPhotonEt", &Dsst_IsoR03_sumPhotonEt);
	tree->Branch("Dsst_IsoR03_sumPUPt", &Dsst_IsoR03_sumPUPt);
	tree->Branch("Dsst_IsoR03_PFIso", &Dsst_IsoR03_PFIso);
	tree->Branch("Dsst_IsoR03_RelIso", &Dsst_IsoR03_RelIso);

	tree->Branch("Dsst_IsoR04_sumChargedHadronPt", &Dsst_IsoR04_sumChargedHadronPt);
	tree->Branch("Dsst_IsoR04_sumNeutralHadronEt", &Dsst_IsoR04_sumNeutralHadronEt);
	tree->Branch("Dsst_IsoR04_sumPhotonEt", &Dsst_IsoR04_sumPhotonEt);
	tree->Branch("Dsst_IsoR04_sumPUPt", &Dsst_IsoR04_sumPUPt);
	tree->Branch("Dsst_IsoR04_PFIso", &Dsst_IsoR04_PFIso);
	tree->Branch("Dsst_IsoR04_RelIso", &Dsst_IsoR04_RelIso);
}

#include "EDAnalyzers/Reconstruction/interface/FullRecoTree.h"
#include <iostream>

FullRecoTree::FullRecoTree(TTree *tree_)
{
	tree = tree_;
}

void FullRecoTree::Init()
{
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

	trig_HLT_IsoMu24 = false;
	trig_HLT_Mu50 = false;

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

	phi_charge = null;
	phi_eta = null;
	phi_phi = null;
	phi_px = null;
	phi_py = null;
	phi_pz = null;
	phi_pt = null;
	phi_p = null;
	phi_mass = null;

	Ds_charge = null;
	Ds_eta = null;
	Ds_phi = null;
	Ds_px = null;
	Ds_py = null;
	Ds_pz = null;
	Ds_pt = null;
	Ds_p = null;
	Ds_mass = null;

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

	dR_Kp_Km = null;
	dR_Kp_pi = null;
	dR_Kp_phi = null;
	dR_Kp_Ds = null;
	dR_Km_pi = null;
	dR_Km_phi = null;
	dR_Km_Ds = null;
	dR_pi_phi = null;
	dR_pi_Ds = null;
	dR_phi_Ds = null;

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

	dira = null;

	Ds_FDxy = null;
	Ds_FDxyerr = null;
	Ds_FDxysig = null;

	Ds_FD = null;
	Ds_FDerr = null;
	Ds_FDsig = null;

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

	num_ss_mu = 0;
	num_os_mu = 0;

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

	dR_Kp_mu = null;
	dR_Km_mu = null;
	dR_pi_mu = null;
	dR_phi_mu = null;
	dR_Ds_mu = null;

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

	met_sig = null;
	met_phi = null;
	met_px = null;
	met_py = null;
	met_pt = null;
	met_sumEt = null;
	met_corphi = null;
	met_corpx = null;
	met_corpy = null;
	met_corpt = null;
	met_corsumEt = null;

	num_jets = null;
	num_jets_pt30 = null;
	deltaR_Ds_Jet = null;
}

void FullRecoTree::CreateBranches()
{
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

	tree->Branch("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	tree->Branch("trig_HLT_Mu50", &trig_HLT_Mu50);

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

	tree->Branch("phi_charge", &phi_charge);
	tree->Branch("phi_eta", &phi_eta);
	tree->Branch("phi_phi", &phi_phi);
	tree->Branch("phi_px", &phi_px);
	tree->Branch("phi_py", &phi_py);
	tree->Branch("phi_pz", &phi_pz);
	tree->Branch("phi_pt", &phi_pt);
	tree->Branch("phi_p", &phi_p);
	tree->Branch("phi_mass", &phi_mass);

	tree->Branch("Ds_charge", &Ds_charge);
	tree->Branch("Ds_eta", &Ds_eta);
	tree->Branch("Ds_phi", &Ds_phi);
	tree->Branch("Ds_px", &Ds_px);
	tree->Branch("Ds_py", &Ds_py);
	tree->Branch("Ds_pz", &Ds_pz);
	tree->Branch("Ds_pt", &Ds_pt);
	tree->Branch("Ds_p", &Ds_p);
	tree->Branch("Ds_mass", &Ds_mass);

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

	tree->Branch("dR_Kp_Km", &dR_Kp_Km);
	tree->Branch("dR_Kp_pi", &dR_Kp_pi);
	tree->Branch("dR_Kp_phi", &dR_Kp_phi);
	tree->Branch("dR_Kp_Ds", &dR_Kp_Ds);
	tree->Branch("dR_Km_pi", &dR_Km_pi);
	tree->Branch("dR_Km_phi", &dR_Km_phi);
	tree->Branch("dR_Km_Ds", &dR_Km_Ds);
	tree->Branch("dR_pi_phi", &dR_pi_phi);
	tree->Branch("dR_pi_Ds", &dR_pi_Ds);
	tree->Branch("dR_phi_Ds", &dR_phi_Ds);

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

	tree->Branch("dira", &dira);

	tree->Branch("Ds_FDxy", &Ds_FDxy);
	tree->Branch("Ds_FDxyerr", &Ds_FDxyerr);
	tree->Branch("Ds_FDxysig", &Ds_FDxysig);

	tree->Branch("Ds_FD", &Ds_FD);
	tree->Branch("Ds_FDerr", &Ds_FDerr);
	tree->Branch("Ds_FDsig", &Ds_FDsig);

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

	tree->Branch("num_ss_mu", &num_ss_mu);
	tree->Branch("num_os_mu", &num_os_mu);

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
	tree->Branch("dR_phi_mu", &dR_phi_mu);
	tree->Branch("dR_Ds_mu", &dR_Ds_mu);

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

	tree->Branch("met_sig", &met_sig);
	tree->Branch("met_phi", &met_phi);
	tree->Branch("met_px", &met_px);
	tree->Branch("met_py", &met_py);
	tree->Branch("met_pt", &met_pt);
	tree->Branch("met_sumEt", &met_sumEt);
	tree->Branch("met_corphi", &met_corphi);
	tree->Branch("met_corpx", &met_corpx);
	tree->Branch("met_corpy", &met_corpy);
	tree->Branch("met_corpt", &met_corpt);
	tree->Branch("met_corsumEt", &met_corsumEt);
	
	tree->Branch("num_jets", &num_jets);
	tree->Branch("num_jets_pt30", &num_jets_pt30);
	tree->Branch("deltaR_Ds_Jet", &deltaR_Ds_Jet);
}

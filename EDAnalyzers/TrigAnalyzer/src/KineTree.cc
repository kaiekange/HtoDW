#include "EDAnalyzers/TrigAnalyzer/interface/KineTree.h"
#include <iostream>

KineTree::KineTree(TTree *tree_)
{
	tree = tree_;
}

void KineTree::Init()
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

	Kp_charge = null;
	Kp_eta = null;
	Kp_phi = null;
	Kp_px = null;
	Kp_py = null;
	Kp_pz = null;
	Kp_pt = null;
	Kp_p = null;

	Km_charge = null;
	Km_eta = null;
	Km_phi = null;
	Km_px = null;
	Km_py = null;
	Km_pz = null;
	Km_pt = null;
	Km_p = null;

	pi_charge = null;
	pi_eta = null;
	pi_phi = null;
	pi_px = null;
	pi_py = null;
	pi_pz = null;
	pi_pt = null;
	pi_p = null;

	phiFit_Kp_charge = null;
	phiFit_Kp_eta = null;
	phiFit_Kp_phi = null;
	phiFit_Kp_px = null;
	phiFit_Kp_py = null;
	phiFit_Kp_pz = null;
	phiFit_Kp_pt = null;
	phiFit_Kp_p = null;
	phiFit_Kp_mass = null;

	phiFit_Km_charge = null;
	phiFit_Km_eta = null;
	phiFit_Km_phi = null;
	phiFit_Km_px = null;
	phiFit_Km_py = null;
	phiFit_Km_pz = null;
	phiFit_Km_pt = null;
	phiFit_Km_p = null;
	phiFit_Km_mass = null;

	phiFit_phi_charge = null;
	phiFit_phi_eta = null;
	phiFit_phi_phi = null;
	phiFit_phi_px = null;
	phiFit_phi_py = null;
	phiFit_phi_pz = null;
	phiFit_phi_pt = null;
	phiFit_phi_p = null;
	phiFit_phi_mass = null;

	DsFit_Kp_charge = null;
	DsFit_Kp_eta = null;
	DsFit_Kp_phi = null;
	DsFit_Kp_px = null;
	DsFit_Kp_py = null;
	DsFit_Kp_pz = null;
	DsFit_Kp_pt = null;
	DsFit_Kp_p = null;
	DsFit_Kp_mass = null;

	DsFit_Km_charge = null;
	DsFit_Km_eta = null;
	DsFit_Km_phi = null;
	DsFit_Km_px = null;
	DsFit_Km_py = null;
	DsFit_Km_pz = null;
	DsFit_Km_pt = null;
	DsFit_Km_p = null;
	DsFit_Km_mass = null;

	DsFit_pi_charge = null;
	DsFit_pi_eta = null;
	DsFit_pi_phi = null;
	DsFit_pi_px = null;
	DsFit_pi_py = null;
	DsFit_pi_pz = null;
	DsFit_pi_pt = null;
	DsFit_pi_p = null;
	DsFit_pi_mass = null;

	DsFit_phi_charge = null;
	DsFit_phi_eta = null;
	DsFit_phi_phi = null;
	DsFit_phi_px = null;
	DsFit_phi_py = null;
	DsFit_phi_pz = null;
	DsFit_phi_pt = null;
	DsFit_phi_p = null;
	DsFit_phi_mass = null;

	DsFit_Ds_charge = null;
	DsFit_Ds_eta = null;
	DsFit_Ds_phi = null;
	DsFit_Ds_px = null;
	DsFit_Ds_py = null;
	DsFit_Ds_pz = null;
	DsFit_Ds_pt = null;
	DsFit_Ds_p = null;
	DsFit_Ds_mass = null;

	phi_vx = null;
	phi_vy = null;
	phi_vz = null;
	phi_vxerr = null;
	phi_vyerr = null;
	phi_vzerr = null;
	phi_chi2 = null;
	phi_ndf = null;
	phi_chi2ndf = null;
	phi_prob = null;

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

	dxy_phi_Ds = null;
	dz_phi_Ds = null;

	dR_Kp_Km = null;
	dR_Kp_pi = null;
	dR_Km_pi = null;

	phiFit_dR_Kp_Km = null;
	phiFit_dR_Kp_pi = null;
	phiFit_dR_Kp_phi = null;
	phiFit_dR_Km_pi = null;
	phiFit_dR_Km_phi = null;
	phiFit_dR_pi_phi = null;

	DsFit_dR_Kp_Km = null;
	DsFit_dR_Kp_pi = null;
	DsFit_dR_Kp_phi = null;
	DsFit_dR_Kp_Ds = null;
	DsFit_dR_Km_pi = null;
	DsFit_dR_Km_phi = null;
	DsFit_dR_Km_Ds = null;
	DsFit_dR_pi_phi = null;
	DsFit_dR_pi_Ds = null;
	DsFit_dR_phi_Ds = null;

	Kp_ip = null;
	Kp_iperr = null;
	Kp_ipsig = null;
	Km_ip = null;
	Km_iperr = null;
	Km_ipsig = null;
	pi_ip = null;
	pi_iperr = null;
	pi_ipsig = null;
	phi_ip = null;
	phi_iperr = null;
	phi_ipsig = null;
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
}

void KineTree::CreateBranches()
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

	tree->Branch("Kp_charge", &Kp_charge);
	tree->Branch("Kp_eta", &Kp_eta);
	tree->Branch("Kp_phi", &Kp_phi);
	tree->Branch("Kp_px", &Kp_px);
	tree->Branch("Kp_py", &Kp_py);
	tree->Branch("Kp_pz", &Kp_pz);
	tree->Branch("Kp_pt", &Kp_pt);
	tree->Branch("Kp_p", &Kp_p);

	tree->Branch("Km_charge", &Km_charge);
	tree->Branch("Km_eta", &Km_eta);
	tree->Branch("Km_phi", &Km_phi);
	tree->Branch("Km_px", &Km_px);
	tree->Branch("Km_py", &Km_py);
	tree->Branch("Km_pz", &Km_pz);
	tree->Branch("Km_pt", &Km_pt);
	tree->Branch("Km_p", &Km_p);

	tree->Branch("pi_charge", &pi_charge);
	tree->Branch("pi_eta", &pi_eta);
	tree->Branch("pi_phi", &pi_phi);
	tree->Branch("pi_px", &pi_px);
	tree->Branch("pi_py", &pi_py);
	tree->Branch("pi_pz", &pi_pz);
	tree->Branch("pi_pt", &pi_pt);
	tree->Branch("pi_p", &pi_p);

	tree->Branch("phiFit_Kp_charge", &phiFit_Kp_charge);
	tree->Branch("phiFit_Kp_eta", &phiFit_Kp_eta);
	tree->Branch("phiFit_Kp_phi", &phiFit_Kp_phi);
	tree->Branch("phiFit_Kp_px", &phiFit_Kp_px);
	tree->Branch("phiFit_Kp_py", &phiFit_Kp_py);
	tree->Branch("phiFit_Kp_pz", &phiFit_Kp_pz);
	tree->Branch("phiFit_Kp_pt", &phiFit_Kp_pt);
	tree->Branch("phiFit_Kp_p", &phiFit_Kp_p);
	tree->Branch("phiFit_Kp_mass", &phiFit_Kp_mass);

	tree->Branch("phiFit_Km_charge", &phiFit_Km_charge);
	tree->Branch("phiFit_Km_eta", &phiFit_Km_eta);
	tree->Branch("phiFit_Km_phi", &phiFit_Km_phi);
	tree->Branch("phiFit_Km_px", &phiFit_Km_px);
	tree->Branch("phiFit_Km_py", &phiFit_Km_py);
	tree->Branch("phiFit_Km_pz", &phiFit_Km_pz);
	tree->Branch("phiFit_Km_pt", &phiFit_Km_pt);
	tree->Branch("phiFit_Km_p", &phiFit_Km_p);
	tree->Branch("phiFit_Km_mass", &phiFit_Km_mass);

	tree->Branch("phiFit_phi_charge", &phiFit_phi_charge);
	tree->Branch("phiFit_phi_eta", &phiFit_phi_eta);
	tree->Branch("phiFit_phi_phi", &phiFit_phi_phi);
	tree->Branch("phiFit_phi_px", &phiFit_phi_px);
	tree->Branch("phiFit_phi_py", &phiFit_phi_py);
	tree->Branch("phiFit_phi_pz", &phiFit_phi_pz);
	tree->Branch("phiFit_phi_pt", &phiFit_phi_pt);
	tree->Branch("phiFit_phi_p", &phiFit_phi_p);
	tree->Branch("phiFit_phi_mass", &phiFit_phi_mass);

	tree->Branch("DsFit_Kp_charge", &DsFit_Kp_charge);
	tree->Branch("DsFit_Kp_eta", &DsFit_Kp_eta);
	tree->Branch("DsFit_Kp_phi", &DsFit_Kp_phi);
	tree->Branch("DsFit_Kp_px", &DsFit_Kp_px);
	tree->Branch("DsFit_Kp_py", &DsFit_Kp_py);
	tree->Branch("DsFit_Kp_pz", &DsFit_Kp_pz);
	tree->Branch("DsFit_Kp_pt", &DsFit_Kp_pt);
	tree->Branch("DsFit_Kp_p", &DsFit_Kp_p);
	tree->Branch("DsFit_Kp_mass", &DsFit_Kp_mass);

	tree->Branch("DsFit_Km_charge", &DsFit_Km_charge);
	tree->Branch("DsFit_Km_eta", &DsFit_Km_eta);
	tree->Branch("DsFit_Km_phi", &DsFit_Km_phi);
	tree->Branch("DsFit_Km_px", &DsFit_Km_px);
	tree->Branch("DsFit_Km_py", &DsFit_Km_py);
	tree->Branch("DsFit_Km_pz", &DsFit_Km_pz);
	tree->Branch("DsFit_Km_pt", &DsFit_Km_pt);
	tree->Branch("DsFit_Km_p", &DsFit_Km_p);
	tree->Branch("DsFit_Km_mass", &DsFit_Km_mass);

	tree->Branch("DsFit_pi_charge", &DsFit_pi_charge);
	tree->Branch("DsFit_pi_eta", &DsFit_pi_eta);
	tree->Branch("DsFit_pi_phi", &DsFit_pi_phi);
	tree->Branch("DsFit_pi_px", &DsFit_pi_px);
	tree->Branch("DsFit_pi_py", &DsFit_pi_py);
	tree->Branch("DsFit_pi_pz", &DsFit_pi_pz);
	tree->Branch("DsFit_pi_pt", &DsFit_pi_pt);
	tree->Branch("DsFit_pi_p", &DsFit_pi_p);
	tree->Branch("DsFit_pi_mass", &DsFit_pi_mass);

	tree->Branch("DsFit_phi_charge", &DsFit_phi_charge);
	tree->Branch("DsFit_phi_eta", &DsFit_phi_eta);
	tree->Branch("DsFit_phi_phi", &DsFit_phi_phi);
	tree->Branch("DsFit_phi_px", &DsFit_phi_px);
	tree->Branch("DsFit_phi_py", &DsFit_phi_py);
	tree->Branch("DsFit_phi_pz", &DsFit_phi_pz);
	tree->Branch("DsFit_phi_pt", &DsFit_phi_pt);
	tree->Branch("DsFit_phi_p", &DsFit_phi_p);
	tree->Branch("DsFit_phi_mass", &DsFit_phi_mass);

	tree->Branch("DsFit_Ds_charge", &DsFit_Ds_charge);
	tree->Branch("DsFit_Ds_eta", &DsFit_Ds_eta);
	tree->Branch("DsFit_Ds_phi", &DsFit_Ds_phi);
	tree->Branch("DsFit_Ds_px", &DsFit_Ds_px);
	tree->Branch("DsFit_Ds_py", &DsFit_Ds_py);
	tree->Branch("DsFit_Ds_pz", &DsFit_Ds_pz);
	tree->Branch("DsFit_Ds_pt", &DsFit_Ds_pt);
	tree->Branch("DsFit_Ds_p", &DsFit_Ds_p);
	tree->Branch("DsFit_Ds_mass", &DsFit_Ds_mass);

	tree->Branch("phi_vx", &phi_vx);
	tree->Branch("phi_vy", &phi_vy);
	tree->Branch("phi_vz", &phi_vz);
	tree->Branch("phi_vxerr", &phi_vxerr);
	tree->Branch("phi_vyerr", &phi_vyerr);
	tree->Branch("phi_vzerr", &phi_vzerr);
	tree->Branch("phi_chi2", &phi_chi2);
	tree->Branch("phi_ndf", &phi_ndf);
	tree->Branch("phi_chi2ndf", &phi_chi2ndf);
	tree->Branch("phi_prob", &phi_prob);

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

	tree->Branch("dxy_phi_Ds", &dxy_phi_Ds);
	tree->Branch("dz_phi_Ds", &dz_phi_Ds);

	tree->Branch("dR_Kp_Km", &dR_Kp_Km);
	tree->Branch("dR_Kp_pi", &dR_Kp_pi);
	tree->Branch("dR_Km_pi", &dR_Km_pi);

	tree->Branch("phiFit_dR_Kp_Km", &phiFit_dR_Kp_Km);
	tree->Branch("phiFit_dR_Kp_pi", &phiFit_dR_Kp_pi);
	tree->Branch("phiFit_dR_Kp_phi", &phiFit_dR_Kp_phi);
	tree->Branch("phiFit_dR_Km_pi", &phiFit_dR_Km_pi);
	tree->Branch("phiFit_dR_Km_phi", &phiFit_dR_Km_phi);
	tree->Branch("phiFit_dR_pi_phi", &phiFit_dR_pi_phi);

	tree->Branch("DsFit_dR_Kp_Km", &DsFit_dR_Kp_Km);
	tree->Branch("DsFit_dR_Kp_pi", &DsFit_dR_Kp_pi);
	tree->Branch("DsFit_dR_Kp_phi", &DsFit_dR_Kp_phi);
	tree->Branch("DsFit_dR_Kp_Ds", &DsFit_dR_Kp_Ds);
	tree->Branch("DsFit_dR_Km_pi", &DsFit_dR_Km_pi);
	tree->Branch("DsFit_dR_Km_phi", &DsFit_dR_Km_phi);
	tree->Branch("DsFit_dR_Km_Ds", &DsFit_dR_Km_Ds);
	tree->Branch("DsFit_dR_pi_phi", &DsFit_dR_pi_phi);
	tree->Branch("DsFit_dR_pi_Ds", &DsFit_dR_pi_Ds);
	tree->Branch("DsFit_dR_phi_Ds", &DsFit_dR_phi_Ds);

	tree->Branch("Kp_ip", &Kp_ip);
	tree->Branch("Kp_iperr", &Kp_iperr);
	tree->Branch("Kp_ipsig", &Kp_ipsig);
	tree->Branch("Km_ip", &Km_ip);
	tree->Branch("Km_iperr", &Km_iperr);
	tree->Branch("Km_ipsig", &Km_ipsig);
	tree->Branch("pi_ip", &pi_ip);
	tree->Branch("pi_iperr", &pi_iperr);
	tree->Branch("pi_ipsig", &pi_ipsig);
	tree->Branch("phi_ip", &phi_ip);
	tree->Branch("phi_iperr", &phi_iperr);
	tree->Branch("phi_ipsig", &phi_ipsig);
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
}

#include "EDAnalyzers/GenAnalyzer/interface/ThreeFitters_Tree.h"
#include <iostream>

ThreeFitters_Tree::ThreeFitters_Tree(TTree *tree_)
{
	tree = tree_;
}

void ThreeFitters_Tree::Init()
{
	Gen_PV_vx = null;
	Gen_PV_vy = null;
	Gen_PV_vz = null;

	Gen_SV_vx = null;
	Gen_SV_vy = null;
	Gen_SV_vz = null;

	Gen_Ds_eta = null;
	Gen_Ds_phi = null;
	Gen_Ds_px = null;
	Gen_Ds_py = null;
	Gen_Ds_pz = null;
	Gen_Ds_pt = null;
	Gen_Ds_p = null;
	Gen_Ds_mass = null;

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

	KVF_KK_mass = null;
	KVF_Kppi_mass = null;
	KVF_Kmpi_mass = null;

	KVF_Ds_eta = null;
	KVF_Ds_phi = null;
	KVF_Ds_px = null;
	KVF_Ds_py = null;
	KVF_Ds_pz = null;
	KVF_Ds_pt = null;
	KVF_Ds_p = null;
	KVF_Ds_mass = null;

	KVF_Ds_vx = null;
	KVF_Ds_vy = null;
	KVF_Ds_vz = null;
	KVF_Ds_vxerr = null;
	KVF_Ds_vyerr = null;
	KVF_Ds_vzerr = null;
	KVF_Ds_chi2 = null;
	KVF_Ds_ndf = null;
	KVF_Ds_chi2ndf = null;
	KVF_Ds_prob = null;

	KVF_dira = null;
	KVF_Ds_FDxy = null;
	KVF_Ds_FDxyerr = null;
	KVF_Ds_FDxysig = null;
	KVF_Ds_FD = null;
	KVF_Ds_FDerr = null;
	KVF_Ds_FDsig = null;

	KLM_KK_mass = null;
	KLM_Kppi_mass = null;
	KLM_Kmpi_mass = null;

	KLM_Ds_eta = null;
	KLM_Ds_phi = null;
	KLM_Ds_px = null;
	KLM_Ds_py = null;
	KLM_Ds_pz = null;
	KLM_Ds_pt = null;
	KLM_Ds_p = null;
	KLM_Ds_mass = null;

	KLM_Ds_vx = null;
	KLM_Ds_vy = null;
	KLM_Ds_vz = null;
	KLM_Ds_vxerr = null;
	KLM_Ds_vyerr = null;
	KLM_Ds_vzerr = null;
	KLM_Ds_chi2 = null;
	KLM_Ds_ndf = null;
	KLM_Ds_chi2ndf = null;
	KLM_Ds_prob = null;

	KLM_dira = null;
	KLM_Ds_FDxy = null;
	KLM_Ds_FDxyerr = null;
	KLM_Ds_FDxysig = null;
	KLM_Ds_FD = null;
	KLM_Ds_FDerr = null;
	KLM_Ds_FDsig = null;

	AVF_KK_mass = null;
	AVF_Kppi_mass = null;
	AVF_Kmpi_mass = null;

	AVF_Ds_eta = null;
	AVF_Ds_phi = null;
	AVF_Ds_px = null;
	AVF_Ds_py = null;
	AVF_Ds_pz = null;
	AVF_Ds_pt = null;
	AVF_Ds_p = null;
	AVF_Ds_mass = null;

	AVF_Ds_vx = null;
	AVF_Ds_vy = null;
	AVF_Ds_vz = null;
	AVF_Ds_vxerr = null;
	AVF_Ds_vyerr = null;
	AVF_Ds_vzerr = null;
	AVF_Ds_chi2 = null;
	AVF_Ds_ndf = null;
	AVF_Ds_chi2ndf = null;
	AVF_Ds_prob = null;

	AVF_dira = null;
	AVF_Ds_FDxy = null;
	AVF_Ds_FDxyerr = null;
	AVF_Ds_FDxysig = null;
	AVF_Ds_FD = null;
	AVF_Ds_FDerr = null;
	AVF_Ds_FDsig = null;
}

void ThreeFitters_Tree::CreateBranches()
{
	tree->Branch("Gen_PV_vx", &Gen_PV_vx);
	tree->Branch("Gen_PV_vy", &Gen_PV_vy);
	tree->Branch("Gen_PV_vz", &Gen_PV_vz);

	tree->Branch("Gen_SV_vx", &Gen_SV_vx);
	tree->Branch("Gen_SV_vy", &Gen_SV_vy);
	tree->Branch("Gen_SV_vz", &Gen_SV_vz);

	tree->Branch("Gen_Ds_eta", &Gen_Ds_eta);
	tree->Branch("Gen_Ds_phi", &Gen_Ds_phi);
	tree->Branch("Gen_Ds_px", &Gen_Ds_px);
	tree->Branch("Gen_Ds_py", &Gen_Ds_py);
	tree->Branch("Gen_Ds_pz", &Gen_Ds_pz);
	tree->Branch("Gen_Ds_pt", &Gen_Ds_pt);
	tree->Branch("Gen_Ds_p", &Gen_Ds_p);
	tree->Branch("Gen_Ds_mass", &Gen_Ds_mass);

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

	tree->Branch("KVF_KK_mass", &KVF_KK_mass);
	tree->Branch("KVF_Kppi_mass", &KVF_Kppi_mass);
	tree->Branch("KVF_Kmpi_mass", &KVF_Kmpi_mass);

	tree->Branch("KVF_Ds_eta", &KVF_Ds_eta);
	tree->Branch("KVF_Ds_phi", &KVF_Ds_phi);
	tree->Branch("KVF_Ds_px", &KVF_Ds_px);
	tree->Branch("KVF_Ds_py", &KVF_Ds_py);
	tree->Branch("KVF_Ds_pz", &KVF_Ds_pz);
	tree->Branch("KVF_Ds_pt", &KVF_Ds_pt);
	tree->Branch("KVF_Ds_p", &KVF_Ds_p);
	tree->Branch("KVF_Ds_mass", &KVF_Ds_mass);

	tree->Branch("KVF_Ds_vx", &KVF_Ds_vx);
	tree->Branch("KVF_Ds_vy", &KVF_Ds_vy);
	tree->Branch("KVF_Ds_vz", &KVF_Ds_vz);
	tree->Branch("KVF_Ds_vxerr", &KVF_Ds_vxerr);
	tree->Branch("KVF_Ds_vyerr", &KVF_Ds_vyerr);
	tree->Branch("KVF_Ds_vzerr", &KVF_Ds_vzerr);
	tree->Branch("KVF_Ds_chi2", &KVF_Ds_chi2);
	tree->Branch("KVF_Ds_ndf", &KVF_Ds_ndf);
	tree->Branch("KVF_Ds_chi2ndf", &KVF_Ds_chi2ndf);
	tree->Branch("KVF_Ds_prob", &KVF_Ds_prob);

	tree->Branch("KVF_dira", &KVF_dira);
	tree->Branch("KVF_Ds_FDxy", &KVF_Ds_FDxy);
	tree->Branch("KVF_Ds_FDxyerr", &KVF_Ds_FDxyerr);
	tree->Branch("KVF_Ds_FDxysig", &KVF_Ds_FDxysig);
	tree->Branch("KVF_Ds_FD", &KVF_Ds_FD);
	tree->Branch("KVF_Ds_FDerr", &KVF_Ds_FDerr);
	tree->Branch("KVF_Ds_FDsig", &KVF_Ds_FDsig);

	tree->Branch("KLM_KK_mass", &KLM_KK_mass);
	tree->Branch("KLM_Kppi_mass", &KLM_Kppi_mass);
	tree->Branch("KLM_Kmpi_mass", &KLM_Kmpi_mass);

	tree->Branch("KLM_Ds_eta", &KLM_Ds_eta);
	tree->Branch("KLM_Ds_phi", &KLM_Ds_phi);
	tree->Branch("KLM_Ds_px", &KLM_Ds_px);
	tree->Branch("KLM_Ds_py", &KLM_Ds_py);
	tree->Branch("KLM_Ds_pz", &KLM_Ds_pz);
	tree->Branch("KLM_Ds_pt", &KLM_Ds_pt);
	tree->Branch("KLM_Ds_p", &KLM_Ds_p);
	tree->Branch("KLM_Ds_mass", &KLM_Ds_mass);

	tree->Branch("KLM_Ds_vx", &KLM_Ds_vx);
	tree->Branch("KLM_Ds_vy", &KLM_Ds_vy);
	tree->Branch("KLM_Ds_vz", &KLM_Ds_vz);
	tree->Branch("KLM_Ds_vxerr", &KLM_Ds_vxerr);
	tree->Branch("KLM_Ds_vyerr", &KLM_Ds_vyerr);
	tree->Branch("KLM_Ds_vzerr", &KLM_Ds_vzerr);
	tree->Branch("KLM_Ds_chi2", &KLM_Ds_chi2);
	tree->Branch("KLM_Ds_ndf", &KLM_Ds_ndf);
	tree->Branch("KLM_Ds_chi2ndf", &KLM_Ds_chi2ndf);
	tree->Branch("KLM_Ds_prob", &KLM_Ds_prob);

	tree->Branch("KLM_dira", &KLM_dira);
	tree->Branch("KLM_Ds_FDxy", &KLM_Ds_FDxy);
	tree->Branch("KLM_Ds_FDxyerr", &KLM_Ds_FDxyerr);
	tree->Branch("KLM_Ds_FDxysig", &KLM_Ds_FDxysig);
	tree->Branch("KLM_Ds_FD", &KLM_Ds_FD);
	tree->Branch("KLM_Ds_FDerr", &KLM_Ds_FDerr);
	tree->Branch("KLM_Ds_FDsig", &KLM_Ds_FDsig);

	tree->Branch("AVF_KK_mass", &AVF_KK_mass);
	tree->Branch("AVF_Kppi_mass", &AVF_Kppi_mass);
	tree->Branch("AVF_Kmpi_mass", &AVF_Kmpi_mass);

	tree->Branch("AVF_Ds_eta", &AVF_Ds_eta);
	tree->Branch("AVF_Ds_phi", &AVF_Ds_phi);
	tree->Branch("AVF_Ds_px", &AVF_Ds_px);
	tree->Branch("AVF_Ds_py", &AVF_Ds_py);
	tree->Branch("AVF_Ds_pz", &AVF_Ds_pz);
	tree->Branch("AVF_Ds_pt", &AVF_Ds_pt);
	tree->Branch("AVF_Ds_p", &AVF_Ds_p);
	tree->Branch("AVF_Ds_mass", &AVF_Ds_mass);

	tree->Branch("AVF_Ds_vx", &AVF_Ds_vx);
	tree->Branch("AVF_Ds_vy", &AVF_Ds_vy);
	tree->Branch("AVF_Ds_vz", &AVF_Ds_vz);
	tree->Branch("AVF_Ds_vxerr", &AVF_Ds_vxerr);
	tree->Branch("AVF_Ds_vyerr", &AVF_Ds_vyerr);
	tree->Branch("AVF_Ds_vzerr", &AVF_Ds_vzerr);
	tree->Branch("AVF_Ds_chi2", &AVF_Ds_chi2);
	tree->Branch("AVF_Ds_ndf", &AVF_Ds_ndf);
	tree->Branch("AVF_Ds_chi2ndf", &AVF_Ds_chi2ndf);
	tree->Branch("AVF_Ds_prob", &AVF_Ds_prob);

	tree->Branch("AVF_dira", &AVF_dira);
	tree->Branch("AVF_Ds_FDxy", &AVF_Ds_FDxy);
	tree->Branch("AVF_Ds_FDxyerr", &AVF_Ds_FDxyerr);
	tree->Branch("AVF_Ds_FDxysig", &AVF_Ds_FDxysig);
	tree->Branch("AVF_Ds_FD", &AVF_Ds_FD);
	tree->Branch("AVF_Ds_FDerr", &AVF_Ds_FDerr);
	tree->Branch("AVF_Ds_FDsig", &AVF_Ds_FDsig);
}

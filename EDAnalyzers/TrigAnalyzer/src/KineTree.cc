#include "EDAnalyzers/TrigAnalyzer/interface/KineTree.h"
#include <iostream>

KineTree::KineTree(TTree *tree_)
{
    tree = tree_;
}

void KineTree::Init()
{
    phi_mass = null;
	phi_mass_refitted = null;
	phi_px = null;
	phi_py = null;
	phi_pz = null;
	phi_pt = null;
	phi_p = null;
	phi_eta = null;
	phi_phi = null;
	phi_vx = null;
	phi_vy = null;
	phi_vz = null;
	phi_chi2 = null;
	phi_ndf = null;
	phi_chi2ndf = null;
	phi_prob = null;

    Ds_mass = null;
	Ds_px = null;
	Ds_py = null;
	Ds_pz = null;
	Ds_pt = null;
	Ds_p = null;
	Ds_eta = null;
	Ds_phi = null;
	Ds_vx = null;
	Ds_vy = null;
	Ds_vz = null;
	Ds_chi2 = null;
	Ds_ndf = null;
	Ds_chi2ndf = null;
	Ds_prob = null;
}

void KineTree::CreateBranches()
{
    tree->Branch("phi_mass", &phi_mass);;
	tree->Branch("phi_mass_refitted", &phi_mass_refitted);;
	tree->Branch("phi_px", &phi_px);;
	tree->Branch("phi_py", &phi_py);;
	tree->Branch("phi_pz", &phi_pz);;
	tree->Branch("phi_pt", &phi_pt);;
	tree->Branch("phi_p", &phi_p);;
	tree->Branch("phi_eta", &phi_eta);;
	tree->Branch("phi_phi", &phi_phi);;
	tree->Branch("phi_vx", &phi_vx);;
	tree->Branch("phi_vy", &phi_vy);;
	tree->Branch("phi_vz", &phi_vz);;
	tree->Branch("phi_chi2", &phi_chi2);;
	tree->Branch("phi_ndf", &phi_ndf);;
	tree->Branch("phi_chi2ndf", &phi_chi2ndf);;
	tree->Branch("phi_prob", &phi_prob);;

    tree->Branch("Ds_mass", &Ds_mass);;
	tree->Branch("Ds_px", &Ds_px);;
	tree->Branch("Ds_py", &Ds_py);;
	tree->Branch("Ds_pz", &Ds_pz);;
	tree->Branch("Ds_pt", &Ds_pt);;
	tree->Branch("Ds_p", &Ds_p);;
	tree->Branch("Ds_eta", &Ds_eta);;
	tree->Branch("Ds_phi", &Ds_phi);;
	tree->Branch("Ds_vx", &Ds_vx);;
	tree->Branch("Ds_vy", &Ds_vy);;
	tree->Branch("Ds_vz", &Ds_vz);;
	tree->Branch("Ds_chi2", &Ds_chi2);;
	tree->Branch("Ds_ndf", &Ds_ndf);;
	tree->Branch("Ds_chi2ndf", &Ds_chi2ndf);;
	tree->Branch("Ds_prob", &Ds_prob);;
}

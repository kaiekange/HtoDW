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

	Float_t phi_mass;
	Float_t phi_mass_refitted;
	Float_t phi_px;
	Float_t phi_py;
	Float_t phi_pz;
	Float_t phi_pt;
	Float_t phi_p;
	Float_t phi_eta;
	Float_t phi_phi;
	Float_t phi_vx;
	Float_t phi_vy;
	Float_t phi_vz;
	Float_t phi_chi2;
	Float_t phi_ndf;
	Float_t phi_chi2ndf;
	Float_t phi_prob;

    Float_t Ds_mass;
	Float_t Ds_px;
	Float_t Ds_py;
	Float_t Ds_pz;
	Float_t Ds_pt;
	Float_t Ds_p;
	Float_t Ds_eta;
	Float_t Ds_phi;
	Float_t Ds_vx;
	Float_t Ds_vy;
	Float_t Ds_vz;
	Float_t Ds_chi2;
	Float_t Ds_ndf;
	Float_t Ds_chi2ndf;
	Float_t Ds_prob;
};
#endif
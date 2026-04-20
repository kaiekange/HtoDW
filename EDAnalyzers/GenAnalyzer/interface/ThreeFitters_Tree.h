#ifndef THREEFITTERS_TREE
#define THREEFITTERS_TREE

#include <TTree.h>

#define null -777

class ThreeFitters_Tree
{
public:
	ThreeFitters_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	Float_t Gen_PV_vx;
	Float_t Gen_PV_vy;
	Float_t Gen_PV_vz;

	Float_t Gen_SV_vx;
	Float_t Gen_SV_vy;
	Float_t Gen_SV_vz;

	Float_t Gen_Ds_eta;
	Float_t Gen_Ds_phi;
	Float_t Gen_Ds_px;
	Float_t Gen_Ds_py;
	Float_t Gen_Ds_pz;
	Float_t Gen_Ds_pt;
	Float_t Gen_Ds_p;
	Float_t Gen_Ds_mass;

	Bool_t primvtx_isValid;
	Bool_t primvtx_isFake;
	Float_t primvtx_chi2;
	Int_t primvtx_ndof;
	Float_t primvtx_chi2ndof;
	Float_t primvtx_vx;
	Float_t primvtx_vy;
	Float_t primvtx_vz;
	Float_t primvtx_vxerr;
	Float_t primvtx_vyerr;
	Float_t primvtx_vzerr;

	Float_t KVF_KK_mass;
	Float_t KVF_Kppi_mass;
	Float_t KVF_Kmpi_mass;

	Float_t KVF_Ds_eta;
	Float_t KVF_Ds_phi;
	Float_t KVF_Ds_px;
	Float_t KVF_Ds_py;
	Float_t KVF_Ds_pz;
	Float_t KVF_Ds_pt;
	Float_t KVF_Ds_p;
	Float_t KVF_Ds_mass;

	Float_t KVF_Ds_vx;
	Float_t KVF_Ds_vy;
	Float_t KVF_Ds_vz;
	Float_t KVF_Ds_vxerr;
	Float_t KVF_Ds_vyerr;
	Float_t KVF_Ds_vzerr;
	Float_t KVF_Ds_chi2;
	Int_t KVF_Ds_ndf;
	Float_t KVF_Ds_chi2ndf;
	Float_t KVF_Ds_prob;

	Float_t KVF_dira;
	Float_t KVF_Ds_FDxy;
	Float_t KVF_Ds_FDxyerr;
	Float_t KVF_Ds_FDxysig;
	Float_t KVF_Ds_FD;
	Float_t KVF_Ds_FDerr;
	Float_t KVF_Ds_FDsig;

	Float_t KLM_KK_mass;
	Float_t KLM_Kppi_mass;
	Float_t KLM_Kmpi_mass;

	Float_t KLM_Ds_eta;
	Float_t KLM_Ds_phi;
	Float_t KLM_Ds_px;
	Float_t KLM_Ds_py;
	Float_t KLM_Ds_pz;
	Float_t KLM_Ds_pt;
	Float_t KLM_Ds_p;
	Float_t KLM_Ds_mass;

	Float_t KLM_Ds_vx;
	Float_t KLM_Ds_vy;
	Float_t KLM_Ds_vz;
	Float_t KLM_Ds_vxerr;
	Float_t KLM_Ds_vyerr;
	Float_t KLM_Ds_vzerr;
	Float_t KLM_Ds_chi2;
	Int_t KLM_Ds_ndf;
	Float_t KLM_Ds_chi2ndf;
	Float_t KLM_Ds_prob;

	Float_t KLM_dira;
	Float_t KLM_Ds_FDxy;
	Float_t KLM_Ds_FDxyerr;
	Float_t KLM_Ds_FDxysig;
	Float_t KLM_Ds_FD;
	Float_t KLM_Ds_FDerr;
	Float_t KLM_Ds_FDsig;

	Float_t AVF_KK_mass;
	Float_t AVF_Kppi_mass;
	Float_t AVF_Kmpi_mass;

	Float_t AVF_Ds_eta;
	Float_t AVF_Ds_phi;
	Float_t AVF_Ds_px;
	Float_t AVF_Ds_py;
	Float_t AVF_Ds_pz;
	Float_t AVF_Ds_pt;
	Float_t AVF_Ds_p;
	Float_t AVF_Ds_mass;

	Float_t AVF_Ds_vx;
	Float_t AVF_Ds_vy;
	Float_t AVF_Ds_vz;
	Float_t AVF_Ds_vxerr;
	Float_t AVF_Ds_vyerr;
	Float_t AVF_Ds_vzerr;
	Float_t AVF_Ds_chi2;
	Int_t AVF_Ds_ndf;
	Float_t AVF_Ds_chi2ndf;
	Float_t AVF_Ds_prob;

	Float_t AVF_dira;
	Float_t AVF_Ds_FDxy;
	Float_t AVF_Ds_FDxyerr;
	Float_t AVF_Ds_FDxysig;
	Float_t AVF_Ds_FD;
	Float_t AVF_Ds_FDerr;
	Float_t AVF_Ds_FDsig;
};
#endif

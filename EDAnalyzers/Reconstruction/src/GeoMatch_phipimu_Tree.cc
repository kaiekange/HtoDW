#include "EDAnalyzers/Reconstruction/interface/GeoMatch_phipimu_Tree.h"

GeoMatch_phipimu_Tree::GeoMatch_phipimu_Tree(TTree *tree_)
{
	tree = tree_;
}

void GeoMatch_phipimu_Tree::Init()
{
	Gen_PV_vx = null; Gen_PV_vy = null; Gen_PV_vz = null;
	Gen_SV_vx = null; Gen_SV_vy = null; Gen_SV_vz = null;

	Gen_H_charge = null;
	Gen_H_eta = null; Gen_H_phi = null; Gen_H_px = null; Gen_H_py = null; Gen_H_pz = null;
	Gen_H_pt = null; Gen_H_p = null; Gen_H_mass = null;

	Gen_Dsst_charge = null;
	Gen_Dsst_eta = null; Gen_Dsst_phi = null; Gen_Dsst_px = null; Gen_Dsst_py = null; Gen_Dsst_pz = null;
	Gen_Dsst_pt = null; Gen_Dsst_p = null; Gen_Dsst_mass = null;

	Gen_Ds_charge = null;
	Gen_Ds_eta = null; Gen_Ds_phi = null; Gen_Ds_px = null; Gen_Ds_py = null; Gen_Ds_pz = null;
	Gen_Ds_pt = null; Gen_Ds_p = null; Gen_Ds_mass = null;

	Gen_gamma_charge = null;
	Gen_gamma_eta = null; Gen_gamma_phi = null; Gen_gamma_px = null; Gen_gamma_py = null; Gen_gamma_pz = null;
	Gen_gamma_pt = null; Gen_gamma_p = null; Gen_gamma_mass = null;

	Gen_W_charge = null;
	Gen_W_eta = null; Gen_W_phi = null; Gen_W_px = null; Gen_W_py = null; Gen_W_pz = null;
	Gen_W_pt = null; Gen_W_p = null; Gen_W_mass = null;

	Gen_Kp_charge = null;
	Gen_Kp_eta = null; Gen_Kp_phi = null; Gen_Kp_px = null; Gen_Kp_py = null; Gen_Kp_pz = null;
	Gen_Kp_pt = null; Gen_Kp_p = null; Gen_Kp_mass = null;

	Gen_Km_charge = null;
	Gen_Km_eta = null; Gen_Km_phi = null; Gen_Km_px = null; Gen_Km_py = null; Gen_Km_pz = null;
	Gen_Km_pt = null; Gen_Km_p = null; Gen_Km_mass = null;

	Gen_pi_charge = null;
	Gen_pi_eta = null; Gen_pi_phi = null; Gen_pi_px = null; Gen_pi_py = null; Gen_pi_pz = null;
	Gen_pi_pt = null; Gen_pi_p = null; Gen_pi_mass = null;

	Gen_mu_charge = null;
	Gen_mu_eta = null; Gen_mu_phi = null; Gen_mu_px = null; Gen_mu_py = null; Gen_mu_pz = null;
	Gen_mu_pt = null; Gen_mu_p = null; Gen_mu_mass = null;

	Gen_Dsmu_charge = null;
	Gen_Dsmu_eta = null; Gen_Dsmu_phi = null; Gen_Dsmu_px = null; Gen_Dsmu_py = null; Gen_Dsmu_pz = null;
	Gen_Dsmu_pt = null; Gen_Dsmu_p = null; Gen_Dsmu_mass = null;

	Gen_Dsstmu_charge = null;
	Gen_Dsstmu_eta = null; Gen_Dsstmu_phi = null; Gen_Dsstmu_px = null; Gen_Dsstmu_py = null; Gen_Dsstmu_pz = null;
	Gen_Dsstmu_pt = null; Gen_Dsstmu_p = null; Gen_Dsstmu_mass = null;

	Gen_dR_Kp_Km = null; Gen_dR_Kp_pi = null; Gen_dR_Kp_Ds = null;
	Gen_dR_Km_pi = null; Gen_dR_Km_Ds = null; Gen_dR_pi_Ds = null;
	Gen_dR_Kp_gamma = null; Gen_dR_Kp_Dsst = null;
	Gen_dR_Km_gamma = null; Gen_dR_Km_Dsst = null;
	Gen_dR_pi_gamma = null; Gen_dR_pi_Dsst = null;
	Gen_dR_Ds_gamma = null; Gen_dR_Ds_Dsst = null;
	Gen_dR_gamma_Dsst = null;
	Gen_dR_Kp_mu = null; Gen_dR_Km_mu = null; Gen_dR_pi_mu = null; Gen_dR_Ds_mu = null;
	Gen_dR_gamma_mu = null; Gen_dR_Dsst_mu = null;

	Gen_met_phi = null; Gen_met_px = null; Gen_met_py = null; Gen_met_pt = null; Gen_met_sumEt = null;

	trig_HLT_IsoMu24 = false;
	trig_HLT_Mu50 = false;

	primvtx_isValid = false;
	primvtx_isFake = false;
	primvtx_chi2 = null; primvtx_ndof = null; primvtx_chi2ndof = null;
	primvtx_vx = null; primvtx_vy = null; primvtx_vz = null;
	primvtx_vxerr = null; primvtx_vyerr = null; primvtx_vzerr = null;

	num_rs_mu = 0; num_ws_mu = 0;

	mu_charge = null;
	mu_eta = null; mu_phi = null; mu_px = null; mu_py = null; mu_pz = null; mu_pt = null; mu_p = null;
	mu_isHighPt = false; mu_isLoose = false; mu_isMedium = false; mu_isSoft = false; mu_isTight = false;
	mu_isPF = false; mu_isTracker = false; mu_isGlobal = false;
	mu_ip = null; mu_iperr = null; mu_ipsig = null;
	mu_ipxy = null; mu_ipxyerr = null; mu_ipxysig = null;
	mu_IsoR03_sumChargedHadronPt = null; mu_IsoR03_sumNeutralHadronEt = null;
	mu_IsoR03_sumPhotonEt = null; mu_IsoR03_sumPUPt = null;
	mu_IsoR03_PFIso = null; mu_IsoR03_RelIso = null;
	mu_IsoR04_sumChargedHadronPt = null; mu_IsoR04_sumNeutralHadronEt = null;
	mu_IsoR04_sumPhotonEt = null; mu_IsoR04_sumPUPt = null;
	mu_IsoR04_PFIso = null; mu_IsoR04_RelIso = null;
	mu_match = false;

	Kp_pf_match = false; Km_pf_match = false; pi_pf_match = false; gamma_pf_match = false;
	Kp_geo_match = false; Km_geo_match = false; pi_geo_match = false;

	Kp_pf_idx = null; Km_pf_idx = null; pi_pf_idx = null;
	Kp_pf_key = null; Km_pf_key = null; pi_pf_key = null;
	Kp_geo_idx = null; Km_geo_idx = null; pi_geo_idx = null;

	Kp_pf_dR = null; Kp_pf_dpt = null;
	Km_pf_dR = null; Km_pf_dpt = null;
	pi_pf_dR = null; pi_pf_dpt = null;
	Kp_geo_dR = null; Kp_geo_dpt = null;
	Km_geo_dR = null; Km_geo_dpt = null;
	pi_geo_dR = null; pi_geo_dpt = null;

	Kp_pf_pool = false; Km_pf_pool = false; pi_pf_pool = false;
	Kp_geo_pool = false; Km_geo_pool = false; pi_geo_pool = false;

	Kp_geo_charge = null;
	Kp_geo_px = null; Kp_geo_py = null; Kp_geo_pz = null; Kp_geo_pt = null; Kp_geo_eta = null;
	Km_geo_charge = null;
	Km_geo_px = null; Km_geo_py = null; Km_geo_pz = null; Km_geo_pt = null; Km_geo_eta = null;
	pi_geo_charge = null;
	pi_geo_px = null; pi_geo_py = null; pi_geo_pz = null; pi_geo_pt = null; pi_geo_eta = null;

	Kp_geo_assoc_key = null; Km_geo_assoc_key = null; pi_geo_assoc_key = null;
	Kp_geo_assoc_pdgId = null; Km_geo_assoc_pdgId = null; pi_geo_assoc_pdgId = null;

	Gen_Kp_idx = null; Gen_Km_idx = null; Gen_pi_idx = null;

	event_id = 0; run_id = 0; lumi_id = 0;
}

void GeoMatch_phipimu_Tree::CreateBranches()
{
	tree->Branch("event_id", &event_id);
	tree->Branch("run_id", &run_id);
	tree->Branch("lumi_id", &lumi_id);

	tree->Branch("Gen_PV_vx", &Gen_PV_vx);
	tree->Branch("Gen_PV_vy", &Gen_PV_vy);
	tree->Branch("Gen_PV_vz", &Gen_PV_vz);
	tree->Branch("Gen_SV_vx", &Gen_SV_vx);
	tree->Branch("Gen_SV_vy", &Gen_SV_vy);
	tree->Branch("Gen_SV_vz", &Gen_SV_vz);

	tree->Branch("Gen_H_charge", &Gen_H_charge);
	tree->Branch("Gen_H_eta", &Gen_H_eta); tree->Branch("Gen_H_phi", &Gen_H_phi);
	tree->Branch("Gen_H_px", &Gen_H_px); tree->Branch("Gen_H_py", &Gen_H_py); tree->Branch("Gen_H_pz", &Gen_H_pz);
	tree->Branch("Gen_H_pt", &Gen_H_pt); tree->Branch("Gen_H_p", &Gen_H_p); tree->Branch("Gen_H_mass", &Gen_H_mass);

	tree->Branch("Gen_Dsst_charge", &Gen_Dsst_charge);
	tree->Branch("Gen_Dsst_eta", &Gen_Dsst_eta); tree->Branch("Gen_Dsst_phi", &Gen_Dsst_phi);
	tree->Branch("Gen_Dsst_px", &Gen_Dsst_px); tree->Branch("Gen_Dsst_py", &Gen_Dsst_py); tree->Branch("Gen_Dsst_pz", &Gen_Dsst_pz);
	tree->Branch("Gen_Dsst_pt", &Gen_Dsst_pt); tree->Branch("Gen_Dsst_p", &Gen_Dsst_p); tree->Branch("Gen_Dsst_mass", &Gen_Dsst_mass);

	tree->Branch("Gen_Ds_charge", &Gen_Ds_charge);
	tree->Branch("Gen_Ds_eta", &Gen_Ds_eta); tree->Branch("Gen_Ds_phi", &Gen_Ds_phi);
	tree->Branch("Gen_Ds_px", &Gen_Ds_px); tree->Branch("Gen_Ds_py", &Gen_Ds_py); tree->Branch("Gen_Ds_pz", &Gen_Ds_pz);
	tree->Branch("Gen_Ds_pt", &Gen_Ds_pt); tree->Branch("Gen_Ds_p", &Gen_Ds_p); tree->Branch("Gen_Ds_mass", &Gen_Ds_mass);

	tree->Branch("Gen_gamma_charge", &Gen_gamma_charge);
	tree->Branch("Gen_gamma_eta", &Gen_gamma_eta); tree->Branch("Gen_gamma_phi", &Gen_gamma_phi);
	tree->Branch("Gen_gamma_px", &Gen_gamma_px); tree->Branch("Gen_gamma_py", &Gen_gamma_py); tree->Branch("Gen_gamma_pz", &Gen_gamma_pz);
	tree->Branch("Gen_gamma_pt", &Gen_gamma_pt); tree->Branch("Gen_gamma_p", &Gen_gamma_p); tree->Branch("Gen_gamma_mass", &Gen_gamma_mass);

	tree->Branch("Gen_W_charge", &Gen_W_charge);
	tree->Branch("Gen_W_eta", &Gen_W_eta); tree->Branch("Gen_W_phi", &Gen_W_phi);
	tree->Branch("Gen_W_px", &Gen_W_px); tree->Branch("Gen_W_py", &Gen_W_py); tree->Branch("Gen_W_pz", &Gen_W_pz);
	tree->Branch("Gen_W_pt", &Gen_W_pt); tree->Branch("Gen_W_p", &Gen_W_p); tree->Branch("Gen_W_mass", &Gen_W_mass);

	tree->Branch("Gen_Kp_charge", &Gen_Kp_charge);
	tree->Branch("Gen_Kp_eta", &Gen_Kp_eta); tree->Branch("Gen_Kp_phi", &Gen_Kp_phi);
	tree->Branch("Gen_Kp_px", &Gen_Kp_px); tree->Branch("Gen_Kp_py", &Gen_Kp_py); tree->Branch("Gen_Kp_pz", &Gen_Kp_pz);
	tree->Branch("Gen_Kp_pt", &Gen_Kp_pt); tree->Branch("Gen_Kp_p", &Gen_Kp_p); tree->Branch("Gen_Kp_mass", &Gen_Kp_mass);

	tree->Branch("Gen_Km_charge", &Gen_Km_charge);
	tree->Branch("Gen_Km_eta", &Gen_Km_eta); tree->Branch("Gen_Km_phi", &Gen_Km_phi);
	tree->Branch("Gen_Km_px", &Gen_Km_px); tree->Branch("Gen_Km_py", &Gen_Km_py); tree->Branch("Gen_Km_pz", &Gen_Km_pz);
	tree->Branch("Gen_Km_pt", &Gen_Km_pt); tree->Branch("Gen_Km_p", &Gen_Km_p); tree->Branch("Gen_Km_mass", &Gen_Km_mass);

	tree->Branch("Gen_pi_charge", &Gen_pi_charge);
	tree->Branch("Gen_pi_eta", &Gen_pi_eta); tree->Branch("Gen_pi_phi", &Gen_pi_phi);
	tree->Branch("Gen_pi_px", &Gen_pi_px); tree->Branch("Gen_pi_py", &Gen_pi_py); tree->Branch("Gen_pi_pz", &Gen_pi_pz);
	tree->Branch("Gen_pi_pt", &Gen_pi_pt); tree->Branch("Gen_pi_p", &Gen_pi_p); tree->Branch("Gen_pi_mass", &Gen_pi_mass);

	tree->Branch("Gen_mu_charge", &Gen_mu_charge);
	tree->Branch("Gen_mu_eta", &Gen_mu_eta); tree->Branch("Gen_mu_phi", &Gen_mu_phi);
	tree->Branch("Gen_mu_px", &Gen_mu_px); tree->Branch("Gen_mu_py", &Gen_mu_py); tree->Branch("Gen_mu_pz", &Gen_mu_pz);
	tree->Branch("Gen_mu_pt", &Gen_mu_pt); tree->Branch("Gen_mu_p", &Gen_mu_p); tree->Branch("Gen_mu_mass", &Gen_mu_mass);

	tree->Branch("Gen_Dsmu_charge", &Gen_Dsmu_charge);
	tree->Branch("Gen_Dsmu_eta", &Gen_Dsmu_eta); tree->Branch("Gen_Dsmu_phi", &Gen_Dsmu_phi);
	tree->Branch("Gen_Dsmu_px", &Gen_Dsmu_px); tree->Branch("Gen_Dsmu_py", &Gen_Dsmu_py); tree->Branch("Gen_Dsmu_pz", &Gen_Dsmu_pz);
	tree->Branch("Gen_Dsmu_pt", &Gen_Dsmu_pt); tree->Branch("Gen_Dsmu_p", &Gen_Dsmu_p); tree->Branch("Gen_Dsmu_mass", &Gen_Dsmu_mass);

	tree->Branch("Gen_Dsstmu_charge", &Gen_Dsstmu_charge);
	tree->Branch("Gen_Dsstmu_eta", &Gen_Dsstmu_eta); tree->Branch("Gen_Dsstmu_phi", &Gen_Dsstmu_phi);
	tree->Branch("Gen_Dsstmu_px", &Gen_Dsstmu_px); tree->Branch("Gen_Dsstmu_py", &Gen_Dsstmu_py); tree->Branch("Gen_Dsstmu_pz", &Gen_Dsstmu_pz);
	tree->Branch("Gen_Dsstmu_pt", &Gen_Dsstmu_pt); tree->Branch("Gen_Dsstmu_p", &Gen_Dsstmu_p); tree->Branch("Gen_Dsstmu_mass", &Gen_Dsstmu_mass);

	tree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km); tree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi); tree->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds);
	tree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi); tree->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds); tree->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds);
	tree->Branch("Gen_dR_Kp_gamma", &Gen_dR_Kp_gamma); tree->Branch("Gen_dR_Kp_Dsst", &Gen_dR_Kp_Dsst);
	tree->Branch("Gen_dR_Km_gamma", &Gen_dR_Km_gamma); tree->Branch("Gen_dR_Km_Dsst", &Gen_dR_Km_Dsst);
	tree->Branch("Gen_dR_pi_gamma", &Gen_dR_pi_gamma); tree->Branch("Gen_dR_pi_Dsst", &Gen_dR_pi_Dsst);
	tree->Branch("Gen_dR_Ds_gamma", &Gen_dR_Ds_gamma); tree->Branch("Gen_dR_Ds_Dsst", &Gen_dR_Ds_Dsst);
	tree->Branch("Gen_dR_gamma_Dsst", &Gen_dR_gamma_Dsst);
	tree->Branch("Gen_dR_Kp_mu", &Gen_dR_Kp_mu); tree->Branch("Gen_dR_Km_mu", &Gen_dR_Km_mu);
	tree->Branch("Gen_dR_pi_mu", &Gen_dR_pi_mu); tree->Branch("Gen_dR_Ds_mu", &Gen_dR_Ds_mu);
	tree->Branch("Gen_dR_gamma_mu", &Gen_dR_gamma_mu); tree->Branch("Gen_dR_Dsst_mu", &Gen_dR_Dsst_mu);

	tree->Branch("Gen_met_phi", &Gen_met_phi); tree->Branch("Gen_met_px", &Gen_met_px); tree->Branch("Gen_met_py", &Gen_met_py);
	tree->Branch("Gen_met_pt", &Gen_met_pt); tree->Branch("Gen_met_sumEt", &Gen_met_sumEt);

	tree->Branch("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	tree->Branch("trig_HLT_Mu50", &trig_HLT_Mu50);

	tree->Branch("primvtx_isValid", &primvtx_isValid); tree->Branch("primvtx_isFake", &primvtx_isFake);
	tree->Branch("primvtx_chi2", &primvtx_chi2); tree->Branch("primvtx_ndof", &primvtx_ndof); tree->Branch("primvtx_chi2ndof", &primvtx_chi2ndof);
	tree->Branch("primvtx_vx", &primvtx_vx); tree->Branch("primvtx_vy", &primvtx_vy); tree->Branch("primvtx_vz", &primvtx_vz);
	tree->Branch("primvtx_vxerr", &primvtx_vxerr); tree->Branch("primvtx_vyerr", &primvtx_vyerr); tree->Branch("primvtx_vzerr", &primvtx_vzerr);

	tree->Branch("num_rs_mu", &num_rs_mu); tree->Branch("num_ws_mu", &num_ws_mu);

	tree->Branch("mu_charge", &mu_charge);
	tree->Branch("mu_eta", &mu_eta); tree->Branch("mu_phi", &mu_phi);
	tree->Branch("mu_px", &mu_px); tree->Branch("mu_py", &mu_py); tree->Branch("mu_pz", &mu_pz);
	tree->Branch("mu_pt", &mu_pt); tree->Branch("mu_p", &mu_p);
	tree->Branch("mu_isHighPt", &mu_isHighPt); tree->Branch("mu_isLoose", &mu_isLoose);
	tree->Branch("mu_isMedium", &mu_isMedium); tree->Branch("mu_isSoft", &mu_isSoft);
	tree->Branch("mu_isTight", &mu_isTight); tree->Branch("mu_isPF", &mu_isPF);
	tree->Branch("mu_isTracker", &mu_isTracker); tree->Branch("mu_isGlobal", &mu_isGlobal);
	tree->Branch("mu_ip", &mu_ip); tree->Branch("mu_iperr", &mu_iperr); tree->Branch("mu_ipsig", &mu_ipsig);
	tree->Branch("mu_ipxy", &mu_ipxy); tree->Branch("mu_ipxyerr", &mu_ipxyerr); tree->Branch("mu_ipxysig", &mu_ipxysig);
	tree->Branch("mu_IsoR03_sumChargedHadronPt", &mu_IsoR03_sumChargedHadronPt);
	tree->Branch("mu_IsoR03_sumNeutralHadronEt", &mu_IsoR03_sumNeutralHadronEt);
	tree->Branch("mu_IsoR03_sumPhotonEt", &mu_IsoR03_sumPhotonEt);
	tree->Branch("mu_IsoR03_sumPUPt", &mu_IsoR03_sumPUPt);
	tree->Branch("mu_IsoR03_PFIso", &mu_IsoR03_PFIso); tree->Branch("mu_IsoR03_RelIso", &mu_IsoR03_RelIso);
	tree->Branch("mu_IsoR04_sumChargedHadronPt", &mu_IsoR04_sumChargedHadronPt);
	tree->Branch("mu_IsoR04_sumNeutralHadronEt", &mu_IsoR04_sumNeutralHadronEt);
	tree->Branch("mu_IsoR04_sumPhotonEt", &mu_IsoR04_sumPhotonEt);
	tree->Branch("mu_IsoR04_sumPUPt", &mu_IsoR04_sumPUPt);
	tree->Branch("mu_IsoR04_PFIso", &mu_IsoR04_PFIso); tree->Branch("mu_IsoR04_RelIso", &mu_IsoR04_RelIso);
	tree->Branch("mu_match", &mu_match);

	tree->Branch("Kp_pf_match", &Kp_pf_match); tree->Branch("Km_pf_match", &Km_pf_match); tree->Branch("pi_pf_match", &pi_pf_match);
	tree->Branch("gamma_pf_match", &gamma_pf_match);
	tree->Branch("Kp_geo_match", &Kp_geo_match); tree->Branch("Km_geo_match", &Km_geo_match); tree->Branch("pi_geo_match", &pi_geo_match);

	tree->Branch("Kp_pf_idx", &Kp_pf_idx); tree->Branch("Km_pf_idx", &Km_pf_idx); tree->Branch("pi_pf_idx", &pi_pf_idx);
	tree->Branch("Kp_pf_key", &Kp_pf_key); tree->Branch("Km_pf_key", &Km_pf_key); tree->Branch("pi_pf_key", &pi_pf_key);
	tree->Branch("Kp_geo_idx", &Kp_geo_idx); tree->Branch("Km_geo_idx", &Km_geo_idx); tree->Branch("pi_geo_idx", &pi_geo_idx);

	tree->Branch("Kp_pf_dR", &Kp_pf_dR); tree->Branch("Kp_pf_dpt", &Kp_pf_dpt);
	tree->Branch("Km_pf_dR", &Km_pf_dR); tree->Branch("Km_pf_dpt", &Km_pf_dpt);
	tree->Branch("pi_pf_dR", &pi_pf_dR); tree->Branch("pi_pf_dpt", &pi_pf_dpt);
	tree->Branch("Kp_geo_dR", &Kp_geo_dR); tree->Branch("Kp_geo_dpt", &Kp_geo_dpt);
	tree->Branch("Km_geo_dR", &Km_geo_dR); tree->Branch("Km_geo_dpt", &Km_geo_dpt);
	tree->Branch("pi_geo_dR", &pi_geo_dR); tree->Branch("pi_geo_dpt", &pi_geo_dpt);

	tree->Branch("Kp_pf_pool", &Kp_pf_pool); tree->Branch("Km_pf_pool", &Km_pf_pool); tree->Branch("pi_pf_pool", &pi_pf_pool);
	tree->Branch("Kp_geo_pool", &Kp_geo_pool); tree->Branch("Km_geo_pool", &Km_geo_pool); tree->Branch("pi_geo_pool", &pi_geo_pool);

	tree->Branch("Kp_geo_charge", &Kp_geo_charge);
	tree->Branch("Kp_geo_px", &Kp_geo_px); tree->Branch("Kp_geo_py", &Kp_geo_py); tree->Branch("Kp_geo_pz", &Kp_geo_pz);
	tree->Branch("Kp_geo_pt", &Kp_geo_pt); tree->Branch("Kp_geo_eta", &Kp_geo_eta);
	tree->Branch("Km_geo_charge", &Km_geo_charge);
	tree->Branch("Km_geo_px", &Km_geo_px); tree->Branch("Km_geo_py", &Km_geo_py); tree->Branch("Km_geo_pz", &Km_geo_pz);
	tree->Branch("Km_geo_pt", &Km_geo_pt); tree->Branch("Km_geo_eta", &Km_geo_eta);
	tree->Branch("pi_geo_charge", &pi_geo_charge);
	tree->Branch("pi_geo_px", &pi_geo_px); tree->Branch("pi_geo_py", &pi_geo_py); tree->Branch("pi_geo_pz", &pi_geo_pz);
	tree->Branch("pi_geo_pt", &pi_geo_pt); tree->Branch("pi_geo_eta", &pi_geo_eta);

	tree->Branch("Kp_geo_assoc_key", &Kp_geo_assoc_key); tree->Branch("Km_geo_assoc_key", &Km_geo_assoc_key); tree->Branch("pi_geo_assoc_key", &pi_geo_assoc_key);
	tree->Branch("Kp_geo_assoc_pdgId", &Kp_geo_assoc_pdgId); tree->Branch("Km_geo_assoc_pdgId", &Km_geo_assoc_pdgId); tree->Branch("pi_geo_assoc_pdgId", &pi_geo_assoc_pdgId);

	tree->Branch("Gen_Kp_idx", &Gen_Kp_idx); tree->Branch("Gen_Km_idx", &Gen_Km_idx); tree->Branch("Gen_pi_idx", &Gen_pi_idx);
}

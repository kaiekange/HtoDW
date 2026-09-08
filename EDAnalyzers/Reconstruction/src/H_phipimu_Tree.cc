#include "EDAnalyzers/Reconstruction/interface/H_phipimu_Tree.h"

H_phipimu_Tree::H_phipimu_Tree(TTree *tree_)
{
	tree = tree_;
}

void H_phipimu_Tree::Init()
{
	Gen_W_lep_pdgId = null;
	Gen_tau_decay_mode = null;
	Gen_Dsst_mode = null;
	Gen_pi0gamma1_match = false; Gen_pi0gamma2_match = false;

	trig_HLT_IsoMu24 = false;
	trig_HLT_Mu50 = false;

	primvtx_isValid = false;
	primvtx_isFake = false;
	primvtx_chi2 = null; primvtx_ndof = null; primvtx_chi2ndof = null;
	primvtx_vx = null; primvtx_vy = null; primvtx_vz = null;
	primvtx_vxerr = null; primvtx_vyerr = null; primvtx_vzerr = null;

	met_sig = null; met_phi = null; met_px = null; met_py = null; met_pt = null; met_sumEt = null; W_tmass = null;

	num_jets = 0;
	num_jets_pt20 = 0; num_jets_pt30 = 0; num_jets_pt40 = 0; num_jets_pt50 = 0;
	num_jets_pt60 = 0; num_jets_pt70 = 0; num_jets_pt80 = 0; num_jets_pt90 = 0; num_jets_pt100 = 0;
	jet_ht = null; jet_ht_pt20 = null; jet_leadpt = null; jet_maxbscore = null;
	jet_vetomap_pass = false;
	gen_weight = 1.0;
	nTrueInt = null;
	nPUInt = null;

	num_os_mu_loose = 0; num_ss_mu_loose = 0;
	num_os_mu_medium = 0; num_ss_mu_medium = 0;
	num_os_mu_tight = 0; num_ss_mu_tight = 0;
	num_os_e_cutloose = 0; num_ss_e_cutloose = 0;
	num_os_e_cutmedium = 0; num_ss_e_cutmedium = 0;
	num_os_e_cuttight = 0; num_ss_e_cuttight = 0;
	num_os_e_mvawp80 = 0; num_ss_e_mvawp80 = 0;
	num_os_e_mvawp90 = 0; num_ss_e_mvawp90 = 0;

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

	event_id = 0; run_id = 0; lumi_id = 0; cand_id = 0;

	Kp_charge = null; Kp_eta = null; Kp_phi = null; Kp_px = null; Kp_py = null; Kp_pz = null;
	Kp_pt = null; Kp_p = null; Kp_mass = null; Kp_PVassoc = null;

	Km_charge = null; Km_eta = null; Km_phi = null; Km_px = null; Km_py = null; Km_pz = null;
	Km_pt = null; Km_p = null; Km_mass = null; Km_PVassoc = null;

	pi_charge = null; pi_eta = null; pi_phi = null; pi_px = null; pi_py = null; pi_pz = null;
	pi_pt = null; pi_p = null; pi_mass = null; pi_PVassoc = null;

	phi_charge = null; phi_eta = null; phi_phi = null; phi_px = null; phi_py = null; phi_pz = null;
	phi_pt = null; phi_p = null; phi_mass = null;

	Ds_charge = null; Ds_eta = null; Ds_phi = null; Ds_px = null; Ds_py = null; Ds_pz = null;
	Ds_pt = null; Ds_p = null; Ds_mass = null;

	Dsmu_charge = null; Dsmu_eta = null; Dsmu_phi = null; Dsmu_px = null; Dsmu_py = null; Dsmu_pz = null;
	Dsmu_pt = null; Dsmu_p = null; Dsmu_mass = null; Dsmu_Ht = null;

	n_gamma_window = 0; gamma_uc_dM = null; gamma_uc_pt = null; gamma_uc_dR_Ds = null;

	has_gamma = false;
	gamma_charge = null; gamma_eta = null; gamma_phi = null; gamma_px = null; gamma_py = null; gamma_pz = null;
	gamma_pt = null; gamma_p = null; gamma_mass = null; gamma_PVassoc = null;

	Dsst_charge = null; Dsst_eta = null; Dsst_phi = null; Dsst_px = null; Dsst_py = null; Dsst_pz = null;
	Dsst_pt = null; Dsst_p = null; Dsst_mass = null;

	Dsstmu_charge = null; Dsstmu_eta = null; Dsstmu_phi = null; Dsstmu_px = null; Dsstmu_py = null; Dsstmu_pz = null;
	Dsstmu_pt = null; Dsstmu_p = null; Dsstmu_mass = null; Dsstmu_Ht = null;

	Ds_vx = null; Ds_vy = null; Ds_vz = null;
	Ds_vxerr = null; Ds_vyerr = null; Ds_vzerr = null;
	Ds_chi2 = null; Ds_ndf = null; Ds_chi2ndf = null; Ds_prob = null;

	dR_Kp_Km = null; dR_Kp_pi = null; dR_Kp_Ds = null;
	dR_Km_pi = null; dR_Km_Ds = null; dR_pi_Ds = null;
	dR_phi_pi = null; dR_phi_Ds = null;
	dR_Kp_mu = null; dR_Km_mu = null; dR_pi_mu = null; dR_phi_mu = null; dR_Ds_mu = null;
	dR_Kp_gamma = null; dR_Kp_Dsst = null;
	dR_Km_gamma = null; dR_Km_Dsst = null;
	dR_pi_gamma = null; dR_pi_Dsst = null;
	dR_phi_gamma = null; dR_phi_Dsst = null;
	dR_Ds_gamma = null; dR_Ds_Dsst = null;
	dR_gamma_Dsst = null; dR_gamma_mu = null; dR_Dsst_mu = null;

	Kp_ip = null; Kp_iperr = null; Kp_ipsig = null;
	Kp_ipxy = null; Kp_ipxyerr = null; Kp_ipxysig = null;
	Km_ip = null; Km_iperr = null; Km_ipsig = null;
	Km_ipxy = null; Km_ipxyerr = null; Km_ipxysig = null;
	pi_ip = null; pi_iperr = null; pi_ipsig = null;
	pi_ipxy = null; pi_ipxyerr = null; pi_ipxysig = null;
	Ds_ip = null; Ds_iperr = null; Ds_ipsig = null;
	Ds_ipxy = null; Ds_ipxyerr = null; Ds_ipxysig = null;

	dira = null;
	Ds_FDxy = null; Ds_FDxyerr = null; Ds_FDxysig = null;
	Ds_FD = null; Ds_FDerr = null; Ds_FDsig = null;

	Ds_IsoR03_sumChargedHadronPt = null; Ds_IsoR03_sumNeutralHadronEt = null;
	Ds_IsoR03_sumPhotonEt = null; Ds_IsoR03_sumPUPt = null;
	Ds_IsoR03_PFIso = null; Ds_IsoR03_RelIso = null;
	Ds_IsoR04_sumChargedHadronPt = null; Ds_IsoR04_sumNeutralHadronEt = null;
	Ds_IsoR04_sumPhotonEt = null; Ds_IsoR04_sumPUPt = null;
	Ds_IsoR04_PFIso = null; Ds_IsoR04_RelIso = null;
	Ds_nearestjet_dR = null; Ds_nearestjet_bscore = null;
	n_b_jets_loose_dR03 = 0; n_b_jets_loose_dR04 = 0;
	n_b_jets_medium_dR03 = 0; n_b_jets_medium_dR04 = 0;
	n_b_jets_tight_dR03 = 0; n_b_jets_tight_dR04 = 0;
	n_b_jets_loose_pt20_dR03 = 0; n_b_jets_loose_pt20_dR04 = 0;
	n_b_jets_medium_pt20_dR03 = 0; n_b_jets_medium_pt20_dR04 = 0;
	n_b_jets_tight_pt20_dR03 = 0; n_b_jets_tight_pt20_dR04 = 0;

	Dsst_IsoR03_sumChargedHadronPt = null; Dsst_IsoR03_sumNeutralHadronEt = null;
	Dsst_IsoR03_sumPhotonEt = null; Dsst_IsoR03_sumPUPt = null;
	Dsst_IsoR03_PFIso = null; Dsst_IsoR03_RelIso = null;
	Dsst_IsoR04_sumChargedHadronPt = null; Dsst_IsoR04_sumNeutralHadronEt = null;
	Dsst_IsoR04_sumPhotonEt = null; Dsst_IsoR04_sumPUPt = null;
	Dsst_IsoR04_PFIso = null; Dsst_IsoR04_RelIso = null;

	Kp_match = false; Km_match = false; phi_match = false; pi_match = false; Ds_match = false;
	gamma_match = false; Dsst_match = false; mu_match = false; event_match = false;
}

void H_phipimu_Tree::InitGamma()
{
	has_gamma = false;
	gamma_charge = null; gamma_eta = null; gamma_phi = null; gamma_px = null; gamma_py = null; gamma_pz = null;
	gamma_pt = null; gamma_p = null; gamma_mass = null; gamma_PVassoc = null;

	Dsst_charge = null; Dsst_eta = null; Dsst_phi = null; Dsst_px = null; Dsst_py = null; Dsst_pz = null;
	Dsst_pt = null; Dsst_p = null; Dsst_mass = null;

	Dsstmu_charge = null; Dsstmu_eta = null; Dsstmu_phi = null; Dsstmu_px = null; Dsstmu_py = null; Dsstmu_pz = null;
	Dsstmu_pt = null; Dsstmu_p = null; Dsstmu_mass = null; Dsstmu_Ht = null;

	dR_Kp_gamma = null; dR_Kp_Dsst = null;
	dR_Km_gamma = null; dR_Km_Dsst = null;
	dR_pi_gamma = null; dR_pi_Dsst = null;
	dR_phi_gamma = null; dR_phi_Dsst = null;
	dR_Ds_gamma = null; dR_Ds_Dsst = null;
	dR_gamma_Dsst = null; dR_gamma_mu = null; dR_Dsst_mu = null;

	Dsst_IsoR03_sumChargedHadronPt = null; Dsst_IsoR03_sumNeutralHadronEt = null;
	Dsst_IsoR03_sumPhotonEt = null; Dsst_IsoR03_sumPUPt = null;
	Dsst_IsoR03_PFIso = null; Dsst_IsoR03_RelIso = null;
	Dsst_IsoR04_sumChargedHadronPt = null; Dsst_IsoR04_sumNeutralHadronEt = null;
	Dsst_IsoR04_sumPhotonEt = null; Dsst_IsoR04_sumPUPt = null;
	Dsst_IsoR04_PFIso = null; Dsst_IsoR04_RelIso = null;

	gamma_match = false; Dsst_match = false; event_match = false;
}

void H_phipimu_Tree::CreateBranches()
{
	tree->Branch("event_id", &event_id);
	tree->Branch("run_id", &run_id);
	tree->Branch("lumi_id", &lumi_id);
	tree->Branch("cand_id", &cand_id);

	tree->Branch("Gen_W_lep_pdgId", &Gen_W_lep_pdgId);
	tree->Branch("Gen_tau_decay_mode", &Gen_tau_decay_mode);
	tree->Branch("Gen_Dsst_mode", &Gen_Dsst_mode);
	tree->Branch("Gen_pi0gamma1_match", &Gen_pi0gamma1_match); tree->Branch("Gen_pi0gamma2_match", &Gen_pi0gamma2_match);

	tree->Branch("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	tree->Branch("trig_HLT_Mu50", &trig_HLT_Mu50);

	tree->Branch("primvtx_isValid", &primvtx_isValid); tree->Branch("primvtx_isFake", &primvtx_isFake);
	tree->Branch("primvtx_chi2", &primvtx_chi2); tree->Branch("primvtx_ndof", &primvtx_ndof); tree->Branch("primvtx_chi2ndof", &primvtx_chi2ndof);
	tree->Branch("primvtx_vx", &primvtx_vx); tree->Branch("primvtx_vy", &primvtx_vy); tree->Branch("primvtx_vz", &primvtx_vz);
	tree->Branch("primvtx_vxerr", &primvtx_vxerr); tree->Branch("primvtx_vyerr", &primvtx_vyerr); tree->Branch("primvtx_vzerr", &primvtx_vzerr);

	tree->Branch("met_sig", &met_sig); tree->Branch("met_phi", &met_phi);
	tree->Branch("met_px", &met_px); tree->Branch("met_py", &met_py); tree->Branch("met_pt", &met_pt); tree->Branch("met_sumEt", &met_sumEt);
	tree->Branch("W_tmass", &W_tmass);

	tree->Branch("num_jets", &num_jets);
	tree->Branch("num_jets_pt20", &num_jets_pt20); tree->Branch("num_jets_pt30", &num_jets_pt30);
	tree->Branch("num_jets_pt40", &num_jets_pt40); tree->Branch("num_jets_pt50", &num_jets_pt50);
	tree->Branch("num_jets_pt60", &num_jets_pt60); tree->Branch("num_jets_pt70", &num_jets_pt70);
	tree->Branch("num_jets_pt80", &num_jets_pt80); tree->Branch("num_jets_pt90", &num_jets_pt90);
	tree->Branch("num_jets_pt100", &num_jets_pt100);
	tree->Branch("jet_ht", &jet_ht); tree->Branch("jet_ht_pt20", &jet_ht_pt20); tree->Branch("jet_leadpt", &jet_leadpt); tree->Branch("jet_maxbscore", &jet_maxbscore);
	tree->Branch("jet_vetomap_pass", &jet_vetomap_pass);
	tree->Branch("gen_weight", &gen_weight);
	tree->Branch("nTrueInt", &nTrueInt);
	tree->Branch("nPUInt", &nPUInt);

	tree->Branch("num_os_mu_loose", &num_os_mu_loose); tree->Branch("num_ss_mu_loose", &num_ss_mu_loose);
	tree->Branch("num_os_mu_medium", &num_os_mu_medium); tree->Branch("num_ss_mu_medium", &num_ss_mu_medium);
	tree->Branch("num_os_mu_tight", &num_os_mu_tight); tree->Branch("num_ss_mu_tight", &num_ss_mu_tight);
	tree->Branch("num_os_e_cutloose", &num_os_e_cutloose); tree->Branch("num_ss_e_cutloose", &num_ss_e_cutloose);
	tree->Branch("num_os_e_cutmedium", &num_os_e_cutmedium); tree->Branch("num_ss_e_cutmedium", &num_ss_e_cutmedium);
	tree->Branch("num_os_e_cuttight", &num_os_e_cuttight); tree->Branch("num_ss_e_cuttight", &num_ss_e_cuttight);
	tree->Branch("num_os_e_mvawp80", &num_os_e_mvawp80); tree->Branch("num_ss_e_mvawp80", &num_ss_e_mvawp80);
	tree->Branch("num_os_e_mvawp90", &num_os_e_mvawp90); tree->Branch("num_ss_e_mvawp90", &num_ss_e_mvawp90);

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

	tree->Branch("Kp_charge", &Kp_charge);
	tree->Branch("Kp_eta", &Kp_eta); tree->Branch("Kp_phi", &Kp_phi);
	tree->Branch("Kp_px", &Kp_px); tree->Branch("Kp_py", &Kp_py); tree->Branch("Kp_pz", &Kp_pz);
	tree->Branch("Kp_pt", &Kp_pt); tree->Branch("Kp_p", &Kp_p); tree->Branch("Kp_mass", &Kp_mass);
	tree->Branch("Kp_PVassoc", &Kp_PVassoc);

	tree->Branch("Km_charge", &Km_charge);
	tree->Branch("Km_eta", &Km_eta); tree->Branch("Km_phi", &Km_phi);
	tree->Branch("Km_px", &Km_px); tree->Branch("Km_py", &Km_py); tree->Branch("Km_pz", &Km_pz);
	tree->Branch("Km_pt", &Km_pt); tree->Branch("Km_p", &Km_p); tree->Branch("Km_mass", &Km_mass);
	tree->Branch("Km_PVassoc", &Km_PVassoc);

	tree->Branch("pi_charge", &pi_charge);
	tree->Branch("pi_eta", &pi_eta); tree->Branch("pi_phi", &pi_phi);
	tree->Branch("pi_px", &pi_px); tree->Branch("pi_py", &pi_py); tree->Branch("pi_pz", &pi_pz);
	tree->Branch("pi_pt", &pi_pt); tree->Branch("pi_p", &pi_p); tree->Branch("pi_mass", &pi_mass);
	tree->Branch("pi_PVassoc", &pi_PVassoc);

	tree->Branch("phi_charge", &phi_charge);
	tree->Branch("phi_eta", &phi_eta); tree->Branch("phi_phi", &phi_phi);
	tree->Branch("phi_px", &phi_px); tree->Branch("phi_py", &phi_py); tree->Branch("phi_pz", &phi_pz);
	tree->Branch("phi_pt", &phi_pt); tree->Branch("phi_p", &phi_p); tree->Branch("phi_mass", &phi_mass);

	tree->Branch("Ds_charge", &Ds_charge);
	tree->Branch("Ds_eta", &Ds_eta); tree->Branch("Ds_phi", &Ds_phi);
	tree->Branch("Ds_px", &Ds_px); tree->Branch("Ds_py", &Ds_py); tree->Branch("Ds_pz", &Ds_pz);
	tree->Branch("Ds_pt", &Ds_pt); tree->Branch("Ds_p", &Ds_p); tree->Branch("Ds_mass", &Ds_mass);

	tree->Branch("Dsmu_charge", &Dsmu_charge);
	tree->Branch("Dsmu_eta", &Dsmu_eta); tree->Branch("Dsmu_phi", &Dsmu_phi);
	tree->Branch("Dsmu_px", &Dsmu_px); tree->Branch("Dsmu_py", &Dsmu_py); tree->Branch("Dsmu_pz", &Dsmu_pz);
	tree->Branch("Dsmu_pt", &Dsmu_pt); tree->Branch("Dsmu_p", &Dsmu_p);
	tree->Branch("Dsmu_mass", &Dsmu_mass); tree->Branch("Dsmu_Ht", &Dsmu_Ht);

	tree->Branch("n_gamma_window", &n_gamma_window);
	tree->Branch("gamma_uc_dM", &gamma_uc_dM); tree->Branch("gamma_uc_pt", &gamma_uc_pt); tree->Branch("gamma_uc_dR_Ds", &gamma_uc_dR_Ds);

	tree->Branch("has_gamma", &has_gamma);
	tree->Branch("gamma_charge", &gamma_charge);
	tree->Branch("gamma_eta", &gamma_eta); tree->Branch("gamma_phi", &gamma_phi);
	tree->Branch("gamma_px", &gamma_px); tree->Branch("gamma_py", &gamma_py); tree->Branch("gamma_pz", &gamma_pz);
	tree->Branch("gamma_pt", &gamma_pt); tree->Branch("gamma_p", &gamma_p); tree->Branch("gamma_mass", &gamma_mass);
	tree->Branch("gamma_PVassoc", &gamma_PVassoc);

	tree->Branch("Dsst_charge", &Dsst_charge);
	tree->Branch("Dsst_eta", &Dsst_eta); tree->Branch("Dsst_phi", &Dsst_phi);
	tree->Branch("Dsst_px", &Dsst_px); tree->Branch("Dsst_py", &Dsst_py); tree->Branch("Dsst_pz", &Dsst_pz);
	tree->Branch("Dsst_pt", &Dsst_pt); tree->Branch("Dsst_p", &Dsst_p); tree->Branch("Dsst_mass", &Dsst_mass);

	tree->Branch("Dsstmu_charge", &Dsstmu_charge);
	tree->Branch("Dsstmu_eta", &Dsstmu_eta); tree->Branch("Dsstmu_phi", &Dsstmu_phi);
	tree->Branch("Dsstmu_px", &Dsstmu_px); tree->Branch("Dsstmu_py", &Dsstmu_py); tree->Branch("Dsstmu_pz", &Dsstmu_pz);
	tree->Branch("Dsstmu_pt", &Dsstmu_pt); tree->Branch("Dsstmu_p", &Dsstmu_p);
	tree->Branch("Dsstmu_mass", &Dsstmu_mass); tree->Branch("Dsstmu_Ht", &Dsstmu_Ht);

	tree->Branch("Ds_vx", &Ds_vx); tree->Branch("Ds_vy", &Ds_vy); tree->Branch("Ds_vz", &Ds_vz);
	tree->Branch("Ds_vxerr", &Ds_vxerr); tree->Branch("Ds_vyerr", &Ds_vyerr); tree->Branch("Ds_vzerr", &Ds_vzerr);
	tree->Branch("Ds_chi2", &Ds_chi2); tree->Branch("Ds_ndf", &Ds_ndf);
	tree->Branch("Ds_chi2ndf", &Ds_chi2ndf); tree->Branch("Ds_prob", &Ds_prob);

	tree->Branch("dR_Kp_Km", &dR_Kp_Km); tree->Branch("dR_Kp_pi", &dR_Kp_pi); tree->Branch("dR_Kp_Ds", &dR_Kp_Ds);
	tree->Branch("dR_Km_pi", &dR_Km_pi); tree->Branch("dR_Km_Ds", &dR_Km_Ds); tree->Branch("dR_pi_Ds", &dR_pi_Ds);
	tree->Branch("dR_phi_pi", &dR_phi_pi); tree->Branch("dR_phi_Ds", &dR_phi_Ds);
	tree->Branch("dR_Kp_mu", &dR_Kp_mu); tree->Branch("dR_Km_mu", &dR_Km_mu);
	tree->Branch("dR_pi_mu", &dR_pi_mu); tree->Branch("dR_phi_mu", &dR_phi_mu); tree->Branch("dR_Ds_mu", &dR_Ds_mu);
	tree->Branch("dR_Kp_gamma", &dR_Kp_gamma); tree->Branch("dR_Kp_Dsst", &dR_Kp_Dsst);
	tree->Branch("dR_Km_gamma", &dR_Km_gamma); tree->Branch("dR_Km_Dsst", &dR_Km_Dsst);
	tree->Branch("dR_pi_gamma", &dR_pi_gamma); tree->Branch("dR_pi_Dsst", &dR_pi_Dsst);
	tree->Branch("dR_phi_gamma", &dR_phi_gamma); tree->Branch("dR_phi_Dsst", &dR_phi_Dsst);
	tree->Branch("dR_Ds_gamma", &dR_Ds_gamma); tree->Branch("dR_Ds_Dsst", &dR_Ds_Dsst);
	tree->Branch("dR_gamma_Dsst", &dR_gamma_Dsst); tree->Branch("dR_gamma_mu", &dR_gamma_mu); tree->Branch("dR_Dsst_mu", &dR_Dsst_mu);

	tree->Branch("Kp_ip", &Kp_ip); tree->Branch("Kp_iperr", &Kp_iperr); tree->Branch("Kp_ipsig", &Kp_ipsig);
	tree->Branch("Kp_ipxy", &Kp_ipxy); tree->Branch("Kp_ipxyerr", &Kp_ipxyerr); tree->Branch("Kp_ipxysig", &Kp_ipxysig);
	tree->Branch("Km_ip", &Km_ip); tree->Branch("Km_iperr", &Km_iperr); tree->Branch("Km_ipsig", &Km_ipsig);
	tree->Branch("Km_ipxy", &Km_ipxy); tree->Branch("Km_ipxyerr", &Km_ipxyerr); tree->Branch("Km_ipxysig", &Km_ipxysig);
	tree->Branch("pi_ip", &pi_ip); tree->Branch("pi_iperr", &pi_iperr); tree->Branch("pi_ipsig", &pi_ipsig);
	tree->Branch("pi_ipxy", &pi_ipxy); tree->Branch("pi_ipxyerr", &pi_ipxyerr); tree->Branch("pi_ipxysig", &pi_ipxysig);
	tree->Branch("Ds_ip", &Ds_ip); tree->Branch("Ds_iperr", &Ds_iperr); tree->Branch("Ds_ipsig", &Ds_ipsig);
	tree->Branch("Ds_ipxy", &Ds_ipxy); tree->Branch("Ds_ipxyerr", &Ds_ipxyerr); tree->Branch("Ds_ipxysig", &Ds_ipxysig);

	tree->Branch("dira", &dira);
	tree->Branch("Ds_FDxy", &Ds_FDxy); tree->Branch("Ds_FDxyerr", &Ds_FDxyerr); tree->Branch("Ds_FDxysig", &Ds_FDxysig);
	tree->Branch("Ds_FD", &Ds_FD); tree->Branch("Ds_FDerr", &Ds_FDerr); tree->Branch("Ds_FDsig", &Ds_FDsig);

	tree->Branch("Ds_IsoR03_sumChargedHadronPt", &Ds_IsoR03_sumChargedHadronPt);
	tree->Branch("Ds_IsoR03_sumNeutralHadronEt", &Ds_IsoR03_sumNeutralHadronEt);
	tree->Branch("Ds_IsoR03_sumPhotonEt", &Ds_IsoR03_sumPhotonEt);
	tree->Branch("Ds_IsoR03_sumPUPt", &Ds_IsoR03_sumPUPt);
	tree->Branch("Ds_IsoR03_PFIso", &Ds_IsoR03_PFIso); tree->Branch("Ds_IsoR03_RelIso", &Ds_IsoR03_RelIso);
	tree->Branch("Ds_IsoR04_sumChargedHadronPt", &Ds_IsoR04_sumChargedHadronPt);
	tree->Branch("Ds_IsoR04_sumNeutralHadronEt", &Ds_IsoR04_sumNeutralHadronEt);
	tree->Branch("Ds_IsoR04_sumPhotonEt", &Ds_IsoR04_sumPhotonEt);
	tree->Branch("Ds_IsoR04_sumPUPt", &Ds_IsoR04_sumPUPt);
	tree->Branch("Ds_IsoR04_PFIso", &Ds_IsoR04_PFIso); tree->Branch("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso);
	tree->Branch("Ds_nearestjet_dR", &Ds_nearestjet_dR); tree->Branch("Ds_nearestjet_bscore", &Ds_nearestjet_bscore);
	tree->Branch("n_b_jets_loose_dR03", &n_b_jets_loose_dR03); tree->Branch("n_b_jets_loose_dR04", &n_b_jets_loose_dR04);
	tree->Branch("n_b_jets_medium_dR03", &n_b_jets_medium_dR03); tree->Branch("n_b_jets_medium_dR04", &n_b_jets_medium_dR04);
	tree->Branch("n_b_jets_tight_dR03", &n_b_jets_tight_dR03); tree->Branch("n_b_jets_tight_dR04", &n_b_jets_tight_dR04);
	tree->Branch("n_b_jets_loose_pt20_dR03", &n_b_jets_loose_pt20_dR03); tree->Branch("n_b_jets_loose_pt20_dR04", &n_b_jets_loose_pt20_dR04);
	tree->Branch("n_b_jets_medium_pt20_dR03", &n_b_jets_medium_pt20_dR03); tree->Branch("n_b_jets_medium_pt20_dR04", &n_b_jets_medium_pt20_dR04);
	tree->Branch("n_b_jets_tight_pt20_dR03", &n_b_jets_tight_pt20_dR03); tree->Branch("n_b_jets_tight_pt20_dR04", &n_b_jets_tight_pt20_dR04);

	tree->Branch("Dsst_IsoR03_sumChargedHadronPt", &Dsst_IsoR03_sumChargedHadronPt);
	tree->Branch("Dsst_IsoR03_sumNeutralHadronEt", &Dsst_IsoR03_sumNeutralHadronEt);
	tree->Branch("Dsst_IsoR03_sumPhotonEt", &Dsst_IsoR03_sumPhotonEt);
	tree->Branch("Dsst_IsoR03_sumPUPt", &Dsst_IsoR03_sumPUPt);
	tree->Branch("Dsst_IsoR03_PFIso", &Dsst_IsoR03_PFIso); tree->Branch("Dsst_IsoR03_RelIso", &Dsst_IsoR03_RelIso);
	tree->Branch("Dsst_IsoR04_sumChargedHadronPt", &Dsst_IsoR04_sumChargedHadronPt);
	tree->Branch("Dsst_IsoR04_sumNeutralHadronEt", &Dsst_IsoR04_sumNeutralHadronEt);
	tree->Branch("Dsst_IsoR04_sumPhotonEt", &Dsst_IsoR04_sumPhotonEt);
	tree->Branch("Dsst_IsoR04_sumPUPt", &Dsst_IsoR04_sumPUPt);
	tree->Branch("Dsst_IsoR04_PFIso", &Dsst_IsoR04_PFIso); tree->Branch("Dsst_IsoR04_RelIso", &Dsst_IsoR04_RelIso);

	tree->Branch("Kp_match", &Kp_match); tree->Branch("Km_match", &Km_match);
	tree->Branch("phi_match", &phi_match); tree->Branch("pi_match", &pi_match); tree->Branch("Ds_match", &Ds_match);
	tree->Branch("gamma_match", &gamma_match); tree->Branch("Dsst_match", &Dsst_match);
	tree->Branch("mu_match", &mu_match); tree->Branch("event_match", &event_match);
}

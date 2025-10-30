#include "EDAnalyzers/RecoAnalyzer/interface/MuonTree.h"
#include <iostream>

MuonTree::MuonTree(TTree *tree_)
{
    tree = tree_;
}
void MuonTree::Init()
{
    // Gen particles
    num_Gen_mu = 0;
    num_Gen_nu = 0;
    num_Gen_W = 0;

    // Matched particles
    num_match_mu = 0;

    match_mu_idx = null;

    match_mu_charge_vec.clear();
    match_mu_eta_vec.clear();
    match_mu_phi_vec.clear();
    match_mu_vx_vec.clear();
    match_mu_vy_vec.clear();
    match_mu_vz_vec.clear();
    match_mu_p_vec.clear();
    match_mu_pt_vec.clear();
    match_mu_px_vec.clear();
    match_mu_py_vec.clear();
    match_mu_pz_vec.clear();
    match_mu_isHighPt_vec.clear();
    match_mu_isLoose_vec.clear();
    match_mu_isMedium_vec.clear();
    match_mu_isSoft_vec.clear();
    match_mu_isTight_vec.clear();
    match_mu_isPF_vec.clear();
    match_mu_isTracker_vec.clear();
    match_mu_isGlobal_vec.clear();
    match_mu_IsoR03_sumChargedHadronPt_vec.clear();
    match_mu_IsoR03_sumChargedParticlePt_vec.clear();
    match_mu_IsoR03_sumNeutralHadronEt_vec.clear();
    match_mu_IsoR03_sumPhotonEt_vec.clear();
    match_mu_IsoR03_sumPUPt_vec.clear();
    match_mu_PFIsoR03_vec.clear();
    match_mu_IsoR04_sumChargedHadronPt_vec.clear();
    match_mu_IsoR04_sumChargedParticlePt_vec.clear();
    match_mu_IsoR04_sumNeutralHadronEt_vec.clear();
    match_mu_IsoR04_sumPhotonEt_vec.clear();
    match_mu_IsoR04_sumPUPt_vec.clear();
    match_mu_PFIsoR04_vec.clear();
    match_mu_primvtx_dxy_vec.clear();
    match_mu_primvtx_dxyerr_vec.clear();
    match_mu_primvtx_dz_vec.clear();
    match_mu_primvtx_dzerr_vec.clear();
    match_mu_primvtx_ip_vec.clear();
    match_mu_primvtx_iperr_vec.clear();
    match_mu_primvtx_ipchi2_vec.clear();

    num_mu = 0;

    mu_charge_vec.clear();
    mu_eta_vec.clear();
    mu_phi_vec.clear();
    mu_vx_vec.clear();
    mu_vy_vec.clear();
    mu_vz_vec.clear();
    mu_p_vec.clear();
    mu_pt_vec.clear();
    mu_px_vec.clear();
    mu_py_vec.clear();
    mu_pz_vec.clear();
    mu_isHighPt_vec.clear();
    mu_isLoose_vec.clear();
    mu_isMedium_vec.clear();
    mu_isSoft_vec.clear();
    mu_isTight_vec.clear();
    mu_isPF_vec.clear();
    mu_isTracker_vec.clear();
    mu_isGlobal_vec.clear();
    mu_IsoR03_sumChargedHadronPt_vec.clear();
    mu_IsoR03_sumChargedParticlePt_vec.clear();
    mu_IsoR03_sumNeutralHadronEt_vec.clear();
    mu_IsoR03_sumPhotonEt_vec.clear();
    mu_IsoR03_sumPUPt_vec.clear();
    mu_PFIsoR03_vec.clear();
    mu_IsoR04_sumChargedHadronPt_vec.clear();
    mu_IsoR04_sumChargedParticlePt_vec.clear();
    mu_IsoR04_sumNeutralHadronEt_vec.clear();
    mu_IsoR04_sumPhotonEt_vec.clear();
    mu_IsoR04_sumPUPt_vec.clear();
    mu_PFIsoR04_vec.clear();
    mu_primvtx_dxy_vec.clear();
    mu_primvtx_dxyerr_vec.clear();
    mu_primvtx_dz_vec.clear();
    mu_primvtx_dzerr_vec.clear();
    mu_primvtx_ip_vec.clear();
    mu_primvtx_iperr_vec.clear();
    mu_primvtx_ipchi2_vec.clear();

    mu_match_vec.clear();
}

void MuonTree::CreateBranches()
{
    tree->Branch("num_Gen_mu", &num_Gen_mu);
    tree->Branch("num_Gen_nu", &num_Gen_nu);
    tree->Branch("num_Gen_W", &num_Gen_W);

    tree->Branch("Gen_mu_eta", &Gen_mu_eta);
    tree->Branch("Gen_mu_phi", &Gen_mu_phi);
    tree->Branch("Gen_mu_vx", &Gen_mu_vx);
    tree->Branch("Gen_mu_vy", &Gen_mu_vy);
    tree->Branch("Gen_mu_vz", &Gen_mu_vz);
    tree->Branch("Gen_mu_p", &Gen_mu_p);
    tree->Branch("Gen_mu_pt", &Gen_mu_pt);
    tree->Branch("Gen_mu_px", &Gen_mu_px);
    tree->Branch("Gen_mu_py", &Gen_mu_py);
    tree->Branch("Gen_mu_pz", &Gen_mu_pz);

    tree->Branch("Gen_nu_eta", &Gen_nu_eta);
    tree->Branch("Gen_nu_phi", &Gen_nu_phi);
    tree->Branch("Gen_nu_vx", &Gen_nu_vx);
    tree->Branch("Gen_nu_vy", &Gen_nu_vy);
    tree->Branch("Gen_nu_vz", &Gen_nu_vz);
    tree->Branch("Gen_nu_p", &Gen_nu_p);
    tree->Branch("Gen_nu_pt", &Gen_nu_pt);
    tree->Branch("Gen_nu_px", &Gen_nu_px);
    tree->Branch("Gen_nu_py", &Gen_nu_py);
    tree->Branch("Gen_nu_pz", &Gen_nu_pz);

    tree->Branch("Gen_W_eta", &Gen_W_eta);
    tree->Branch("Gen_W_phi", &Gen_W_phi);
    tree->Branch("Gen_W_vx", &Gen_W_vx);
    tree->Branch("Gen_W_vy", &Gen_W_vy);
    tree->Branch("Gen_W_vz", &Gen_W_vz);
    tree->Branch("Gen_W_p", &Gen_W_p);
    tree->Branch("Gen_W_pt", &Gen_W_pt);
    tree->Branch("Gen_W_px", &Gen_W_px);
    tree->Branch("Gen_W_py", &Gen_W_py);
    tree->Branch("Gen_W_pz", &Gen_W_pz);

    tree->Branch("BS_type", &BS_type);
    tree->Branch("BS_x0", &BS_x0);
    tree->Branch("BS_y0", &BS_y0);
    tree->Branch("BS_z0", &BS_z0);
    tree->Branch("BS_sigmaZ", &BS_sigmaZ);
    tree->Branch("BS_dxdz", &BS_dxdz);
    tree->Branch("BS_dydz", &BS_dydz);
    tree->Branch("BS_BWX", &BS_BWX);
    tree->Branch("BS_BWY", &BS_BWY);
    tree->Branch("BS_x0err", &BS_x0err);
    tree->Branch("BS_y0err", &BS_y0err);
    tree->Branch("BS_z0err", &BS_z0err);
    tree->Branch("BS_sigmaZ0err", &BS_sigmaZ0err);
    tree->Branch("BS_dxdzerr", &BS_dxdzerr);
    tree->Branch("BS_dydzerr", &BS_dydzerr);
    tree->Branch("BS_BWXerr", &BS_BWXerr);
    tree->Branch("BS_BWYerr", &BS_BWYerr);
    tree->Branch("BS_emitX", &BS_emitX);
    tree->Branch("BS_emitY", &BS_emitY);
    tree->Branch("BS_betaStar", &BS_betaStar);
   
    tree->Branch("PV_vx", &PV_vx); 
    tree->Branch("PV_vy", &PV_vy); 
    tree->Branch("PV_vz", &PV_vz); 
    tree->Branch("PV_vxerr", &PV_vxerr); 
    tree->Branch("PV_vyerr", &PV_vyerr); 
    tree->Branch("PV_vzerr", &PV_vzerr); 
    
    tree->Branch("num_match_mu", &num_match_mu);

    tree->Branch("match_mu_idx", &match_mu_idx);

    tree->Branch("match_mu_dR", &match_mu_dR_vec);
    tree->Branch("match_mu_charge", &match_mu_charge_vec);
    tree->Branch("match_mu_eta", &match_mu_eta_vec);
    tree->Branch("match_mu_phi", &match_mu_phi_vec);
    tree->Branch("match_mu_vx", &match_mu_vx_vec);
    tree->Branch("match_mu_vy", &match_mu_vy_vec);
    tree->Branch("match_mu_vz", &match_mu_vz_vec);
    tree->Branch("match_mu_p", &match_mu_p_vec);
    tree->Branch("match_mu_pt", &match_mu_pt_vec);
    tree->Branch("match_mu_px", &match_mu_px_vec);
    tree->Branch("match_mu_py", &match_mu_py_vec);
    tree->Branch("match_mu_pz", &match_mu_pz_vec);
    tree->Branch("match_mu_isHighPt", &match_mu_isHighPt_vec);
    tree->Branch("match_mu_isLoose", &match_mu_isLoose_vec);
    tree->Branch("match_mu_isMedium", &match_mu_isMedium_vec);
    tree->Branch("match_mu_isSoft", &match_mu_isSoft_vec);
    tree->Branch("match_mu_isTight", &match_mu_isTight_vec);
    tree->Branch("match_mu_isPF", &match_mu_isPF_vec);
    tree->Branch("match_mu_isTracker", &match_mu_isTracker_vec);
    tree->Branch("match_mu_isGlobal", &match_mu_isGlobal_vec);
    tree->Branch("match_mu_IsoR03_sumChargedHadronPt", &match_mu_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR03_sumChargedParticlePt", &match_mu_IsoR03_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR03_sumNeutralHadronEt", &match_mu_IsoR03_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR03_sumPhotonEt", &match_mu_IsoR03_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR03_sumPUPt", &match_mu_IsoR03_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR03", &match_mu_PFIsoR03_vec);
    tree->Branch("match_mu_IsoR04_sumChargedHadronPt", &match_mu_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("match_mu_IsoR04_sumChargedParticlePt", &match_mu_IsoR04_sumChargedParticlePt_vec);
    tree->Branch("match_mu_IsoR04_sumNeutralHadronEt", &match_mu_IsoR04_sumNeutralHadronEt_vec);
    tree->Branch("match_mu_IsoR04_sumPhotonEt", &match_mu_IsoR04_sumPhotonEt_vec);
    tree->Branch("match_mu_IsoR04_sumPUPt", &match_mu_IsoR04_sumPUPt_vec);
    tree->Branch("match_mu_PFIsoR04", &match_mu_PFIsoR03_vec);
    tree->Branch("match_mu_primvtx_dxy", &match_mu_primvtx_dxy_vec);
    tree->Branch("match_mu_primvtx_dxyerr", &match_mu_primvtx_dxyerr_vec);
    tree->Branch("match_mu_primvtx_dz", &match_mu_primvtx_dz_vec);
    tree->Branch("match_mu_primvtx_dzerr", &match_mu_primvtx_dzerr_vec);
    tree->Branch("match_mu_primvtx_ip", &match_mu_primvtx_ip_vec);
    tree->Branch("match_mu_primvtx_iperr", &match_mu_primvtx_iperr_vec);
    tree->Branch("match_mu_primvtx_ipchi2", &match_mu_primvtx_ipchi2_vec);
    
    tree->Branch("num_mu", &num_mu);

    tree->Branch("mu_charge", &mu_charge_vec);
    tree->Branch("mu_eta", &mu_eta_vec);
    tree->Branch("mu_phi", &mu_phi_vec);
    tree->Branch("mu_vx", &mu_vx_vec);
    tree->Branch("mu_vy", &mu_vy_vec);
    tree->Branch("mu_vz", &mu_vz_vec);
    tree->Branch("mu_p", &mu_p_vec);
    tree->Branch("mu_pt", &mu_pt_vec);
    tree->Branch("mu_px", &mu_px_vec);
    tree->Branch("mu_py", &mu_py_vec);
    tree->Branch("mu_pz", &mu_pz_vec);
    tree->Branch("mu_isHighPt", &mu_isHighPt_vec);
    tree->Branch("mu_isLoose", &mu_isLoose_vec);
    tree->Branch("mu_isMedium", &mu_isMedium_vec);
    tree->Branch("mu_isSoft", &mu_isSoft_vec);
    tree->Branch("mu_isTight", &mu_isTight_vec);
    tree->Branch("mu_isPF", &mu_isPF_vec);
    tree->Branch("mu_isTracker", &mu_isTracker_vec);
    tree->Branch("mu_isGlobal", &mu_isGlobal_vec);
    tree->Branch("mu_IsoR03_sumChargedHadronPt", &mu_IsoR03_sumChargedHadronPt_vec);
    tree->Branch("mu_IsoR03_sumChargedParticlePt", &mu_IsoR03_sumChargedParticlePt_vec);
    tree->Branch("mu_IsoR03_sumNeutralHadronEt", &mu_IsoR03_sumNeutralHadronEt_vec);
    tree->Branch("mu_IsoR03_sumPhotonEt", &mu_IsoR03_sumPhotonEt_vec);
    tree->Branch("mu_IsoR03_sumPUPt", &mu_IsoR03_sumPUPt_vec);
    tree->Branch("mu_PFIsoR03", &mu_PFIsoR03_vec);
    tree->Branch("mu_IsoR04_sumChargedHadronPt", &mu_IsoR04_sumChargedHadronPt_vec);
    tree->Branch("mu_IsoR04_sumChargedParticlePt", &mu_IsoR04_sumChargedParticlePt_vec);
    tree->Branch("mu_IsoR04_sumNeutralHadronEt", &mu_IsoR04_sumNeutralHadronEt_vec);
    tree->Branch("mu_IsoR04_sumPhotonEt", &mu_IsoR04_sumPhotonEt_vec);
    tree->Branch("mu_IsoR04_sumPUPt", &mu_IsoR04_sumPUPt_vec);
    tree->Branch("mu_PFIsoR04", &mu_PFIsoR03_vec);
    tree->Branch("mu_primvtx_dxy", &mu_primvtx_dxy_vec);
    tree->Branch("mu_primvtx_dxyerr", &mu_primvtx_dxyerr_vec);
    tree->Branch("mu_primvtx_dz", &mu_primvtx_dz_vec);
    tree->Branch("mu_primvtx_dzerr", &mu_primvtx_dzerr_vec);
    tree->Branch("mu_primvtx_ip", &mu_primvtx_ip_vec);
    tree->Branch("mu_primvtx_iperr", &mu_primvtx_iperr_vec);
    tree->Branch("mu_primvtx_ipchi2", &mu_primvtx_ipchi2_vec);

    tree->Branch("mu_match", &mu_match_vec);
}

void MuonTree::Gen_Reset()
{
    Gen_mu_eta = null;
    Gen_mu_phi = null;
    Gen_mu_vx = null;
    Gen_mu_vy = null;
    Gen_mu_vz = null;
    Gen_mu_p = null;
    Gen_mu_pt = null;
    Gen_mu_px = null;
    Gen_mu_py = null;
    Gen_mu_pz = null;

    Gen_nu_eta = null;
    Gen_nu_phi = null;
    Gen_nu_vx = null;
    Gen_nu_vy = null;
    Gen_nu_vz = null;
    Gen_nu_p = null;
    Gen_nu_pt = null;
    Gen_nu_px = null;
    Gen_nu_py = null;
    Gen_nu_pz = null;

    Gen_W_eta = null;
    Gen_W_phi = null;
    Gen_W_vx = null;
    Gen_W_vy = null;
    Gen_W_vz = null;
    Gen_W_p = null;
    Gen_W_pt = null;
    Gen_W_px = null;
    Gen_W_py = null;
    Gen_W_pz = null;
}

void MuonTree::BS_Reset()
{
    BS_type = null; 
    BS_x0 = null;
    BS_y0 = null;
    BS_z0 = null;
    BS_sigmaZ = null;
    BS_dxdz = null;
    BS_dydz = null;
    BS_BWX = null;
    BS_BWY = null;
    BS_x0err = null;
    BS_y0err = null;
    BS_z0err = null;
    BS_sigmaZ0err = null;
    BS_dxdzerr = null;
    BS_dydzerr = null;
    BS_BWXerr = null;
    BS_BWYerr = null;
    BS_emitX = null;
    BS_emitY = null;
    BS_betaStar = null;
}

void MuonTree::PV_Reset()
{
    PV_vx = null;
    PV_vy = null;
    PV_vz = null;
    PV_vxerr = null;
    PV_vyerr = null;
    PV_vzerr = null;
}

void MuonTree::Match_Reset()
{
    match_mu_dR = null;
    match_mu_charge = null;
    match_mu_eta = null;
    match_mu_phi = null;
    match_mu_vx = null;
    match_mu_vy = null;
    match_mu_vz = null;
    match_mu_p = null;
    match_mu_pt = null;
    match_mu_px = null;
    match_mu_py = null;
    match_mu_pz = null;
    match_mu_isHighPt = false;
    match_mu_isLoose = false;
    match_mu_isMedium = false;
    match_mu_isSoft = false;
    match_mu_isTight = false;
    match_mu_isPF = false;
    match_mu_isTracker = false;
    match_mu_isGlobal = false;
    match_mu_IsoR03_sumChargedHadronPt = null;
    match_mu_IsoR03_sumChargedParticlePt = null;
    match_mu_IsoR03_sumNeutralHadronEt = null;
    match_mu_IsoR03_sumPhotonEt = null;
    match_mu_IsoR03_sumPUPt = null;
    match_mu_PFIsoR03 = null;
    match_mu_IsoR04_sumChargedHadronPt = null;
    match_mu_IsoR04_sumChargedParticlePt = null;
    match_mu_IsoR04_sumNeutralHadronEt = null;
    match_mu_IsoR04_sumPhotonEt = null;
    match_mu_IsoR04_sumPUPt = null;
    match_mu_PFIsoR04 = null;
    match_mu_primvtx_dxy = null;
    match_mu_primvtx_dxyerr = null;
    match_mu_primvtx_dz = null;
    match_mu_primvtx_dzerr = null;
    match_mu_primvtx_ip = null;
    match_mu_primvtx_iperr = null;
    match_mu_primvtx_ipchi2 = null;
}

void MuonTree::Match_Fill_Vector()
{
    match_mu_dR_vec.push_back(match_mu_dR);
    match_mu_charge_vec.push_back(match_mu_charge);
    match_mu_eta_vec.push_back(match_mu_eta);
    match_mu_phi_vec.push_back(match_mu_phi);
    match_mu_vx_vec.push_back(match_mu_vx);
    match_mu_vy_vec.push_back(match_mu_vy);
    match_mu_vz_vec.push_back(match_mu_vz);
    match_mu_p_vec.push_back(match_mu_p);
    match_mu_pt_vec.push_back(match_mu_pt);
    match_mu_px_vec.push_back(match_mu_px);
    match_mu_py_vec.push_back(match_mu_py);
    match_mu_pz_vec.push_back(match_mu_pz);
    match_mu_isHighPt_vec.push_back(match_mu_isHighPt);
    match_mu_isLoose_vec.push_back(match_mu_isLoose);
    match_mu_isMedium_vec.push_back(match_mu_isMedium);
    match_mu_isSoft_vec.push_back(match_mu_isSoft);
    match_mu_isTight_vec.push_back(match_mu_isTight);
    match_mu_isPF_vec.push_back(match_mu_isPF);
    match_mu_isTracker_vec.push_back(match_mu_isTracker);
    match_mu_isGlobal_vec.push_back(match_mu_isGlobal);
    match_mu_IsoR03_sumChargedHadronPt_vec.push_back(match_mu_IsoR03_sumChargedHadronPt);
    match_mu_IsoR03_sumChargedParticlePt_vec.push_back(match_mu_IsoR03_sumChargedParticlePt);
    match_mu_IsoR03_sumNeutralHadronEt_vec.push_back(match_mu_IsoR03_sumNeutralHadronEt);
    match_mu_IsoR03_sumPhotonEt_vec.push_back(match_mu_IsoR03_sumPhotonEt);
    match_mu_IsoR03_sumPUPt_vec.push_back(match_mu_IsoR03_sumPUPt);
    match_mu_PFIsoR03_vec.push_back(match_mu_PFIsoR03);
    match_mu_IsoR04_sumChargedHadronPt_vec.push_back(match_mu_IsoR04_sumChargedHadronPt);
    match_mu_IsoR04_sumChargedParticlePt_vec.push_back(match_mu_IsoR04_sumChargedParticlePt);
    match_mu_IsoR04_sumNeutralHadronEt_vec.push_back(match_mu_IsoR04_sumNeutralHadronEt);
    match_mu_IsoR04_sumPhotonEt_vec.push_back(match_mu_IsoR04_sumPhotonEt);
    match_mu_IsoR04_sumPUPt_vec.push_back(match_mu_IsoR04_sumPUPt);
    match_mu_PFIsoR04_vec.push_back(match_mu_PFIsoR04);
    match_mu_primvtx_dxy_vec.push_back(match_mu_primvtx_dxy);
    match_mu_primvtx_dxyerr_vec.push_back(match_mu_primvtx_dxyerr);
    match_mu_primvtx_dz_vec.push_back(match_mu_primvtx_dz);
    match_mu_primvtx_dzerr_vec.push_back(match_mu_primvtx_dzerr);
    match_mu_primvtx_ip_vec.push_back(match_mu_primvtx_ip);
    match_mu_primvtx_iperr_vec.push_back(match_mu_primvtx_iperr);
    match_mu_primvtx_ipchi2_vec.push_back(match_mu_primvtx_ipchi2);
}

void MuonTree::mu_Reset()
{
    mu_charge = null;
    mu_eta = null;
    mu_phi = null;
    mu_vx = null;
    mu_vy = null;
    mu_vz = null;
    mu_p = null;
    mu_pt = null;
    mu_px = null;
    mu_py = null;
    mu_pz = null;
    mu_isHighPt = false;
    mu_isLoose = false;
    mu_isMedium = false;
    mu_isSoft = false;
    mu_isTight = false;
    mu_isPF = false;
    mu_isTracker = false;
    mu_isGlobal = false;
    mu_IsoR03_sumChargedHadronPt = null;
    mu_IsoR03_sumChargedParticlePt = null;
    mu_IsoR03_sumNeutralHadronEt = null;
    mu_IsoR03_sumPhotonEt = null;
    mu_IsoR03_sumPUPt = null;
    mu_PFIsoR03 = null;
    mu_IsoR04_sumChargedHadronPt = null;
    mu_IsoR04_sumChargedParticlePt = null;
    mu_IsoR04_sumNeutralHadronEt = null;
    mu_IsoR04_sumPhotonEt = null;
    mu_IsoR04_sumPUPt = null;
    mu_PFIsoR04 = null;
    mu_primvtx_dxy = null;
    mu_primvtx_dxyerr = null;
    mu_primvtx_dz = null;
    mu_primvtx_dzerr = null;
    mu_primvtx_ip = null;
    mu_primvtx_iperr = null;
    mu_primvtx_ipchi2 = null;
    mu_match = false;
}

void MuonTree::Fill_Vector()
{
    mu_charge_vec.push_back(mu_charge);
    mu_eta_vec.push_back(mu_eta);
    mu_phi_vec.push_back(mu_phi);
    mu_vx_vec.push_back(mu_vx);
    mu_vy_vec.push_back(mu_vy);
    mu_vz_vec.push_back(mu_vz);
    mu_p_vec.push_back(mu_p);
    mu_pt_vec.push_back(mu_pt);
    mu_px_vec.push_back(mu_px);
    mu_py_vec.push_back(mu_py);
    mu_pz_vec.push_back(mu_pz);
    mu_isHighPt_vec.push_back(mu_isHighPt);
    mu_isLoose_vec.push_back(mu_isLoose);
    mu_isMedium_vec.push_back(mu_isMedium);
    mu_isSoft_vec.push_back(mu_isSoft);
    mu_isTight_vec.push_back(mu_isTight);
    mu_isPF_vec.push_back(mu_isPF);
    mu_isTracker_vec.push_back(mu_isTracker);
    mu_isGlobal_vec.push_back(mu_isGlobal);
    mu_IsoR03_sumChargedHadronPt_vec.push_back(mu_IsoR03_sumChargedHadronPt);
    mu_IsoR03_sumChargedParticlePt_vec.push_back(mu_IsoR03_sumChargedParticlePt);
    mu_IsoR03_sumNeutralHadronEt_vec.push_back(mu_IsoR03_sumNeutralHadronEt);
    mu_IsoR03_sumPhotonEt_vec.push_back(mu_IsoR03_sumPhotonEt);
    mu_IsoR03_sumPUPt_vec.push_back(mu_IsoR03_sumPUPt);
    mu_PFIsoR03_vec.push_back(mu_PFIsoR03);
    mu_IsoR04_sumChargedHadronPt_vec.push_back(mu_IsoR04_sumChargedHadronPt);
    mu_IsoR04_sumChargedParticlePt_vec.push_back(mu_IsoR04_sumChargedParticlePt);
    mu_IsoR04_sumNeutralHadronEt_vec.push_back(mu_IsoR04_sumNeutralHadronEt);
    mu_IsoR04_sumPhotonEt_vec.push_back(mu_IsoR04_sumPhotonEt);
    mu_IsoR04_sumPUPt_vec.push_back(mu_IsoR04_sumPUPt);
    mu_PFIsoR04_vec.push_back(mu_PFIsoR04);
    mu_primvtx_dxy_vec.push_back(mu_primvtx_dxy);
    mu_primvtx_dxyerr_vec.push_back(mu_primvtx_dxyerr);
    mu_primvtx_dz_vec.push_back(mu_primvtx_dz);
    mu_primvtx_dzerr_vec.push_back(mu_primvtx_dzerr);
    mu_primvtx_ip_vec.push_back(mu_primvtx_ip);
    mu_primvtx_iperr_vec.push_back(mu_primvtx_iperr);
    mu_primvtx_ipchi2_vec.push_back(mu_primvtx_ipchi2);
    mu_match_vec.push_back(mu_match);
}

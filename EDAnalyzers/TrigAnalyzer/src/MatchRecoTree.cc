#include "EDAnalyzers/TrigAnalyzer/interface/MatchRecoTree.h"
#include <iostream>

MatchRecoTree::MatchRecoTree(TTree *tree_)
{
    tree = tree_;
}

void MatchRecoTree::Init()
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

    Km_charge = null;
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

    phi_charge = null;
    phi_eta = null;
    phi_phi = null;
    phi_px = null;
    phi_py = null;
    phi_pz = null;
    phi_pt = null;
    phi_p = null;
    phi_m = null;

    Ds_charge = null;
    Ds_eta = null;
    Ds_phi = null;
    Ds_px = null;
    Ds_py = null;
    Ds_pz = null;
    Ds_pt = null;
    Ds_p = null;
    Ds_m = null;

    phiFit_isValid = false;
    phiFit_isFake = false;
    phiFit_chi2 = null;
    phiFit_ndof = null;
    phiFit_chi2ndof = null;
    phiFit_vx = null;
    phiFit_vy = null;
    phiFit_vz = null;
    phiFit_vxerr = null;
    phiFit_vyerr = null;
    phiFit_vzerr = null;

    DsFit_isValid = false;
    DsFit_isFake = false;
    DsFit_chi2 = null;
    DsFit_ndof = null;
    DsFit_chi2ndof = null;
    DsFit_vx = null;
    DsFit_vy = null;
    DsFit_vz = null;
    DsFit_vxerr = null;
    DsFit_vyerr = null;
    DsFit_vzerr = null;

    orig_dR_Kp_Km = null;
    orig_dR_Kp_pi = null;
    orig_dR_Km_pi = null;

    phiFit_dR_Kp_Km = null;
    phiFit_dR_Kp_pi = null;
    phiFit_dR_Kp_phi = null;
    phiFit_dR_Km_pi = null;
    phiFit_dR_Km_phi = null;
    phiFit_dR_pi_phi = null;

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
    dxy_phi_Ds = null;
    dz_phi_Ds = null;

    Kp_ip = null;
    Kp_iperr = null;
    Kp_ipchi2 = null;
    Km_ip = null;
    Km_iperr = null;
    Km_ipchi2 = null;
    pi_ip = null;
    pi_iperr = null;
    pi_ipchi2 = null;
    Ds_FDxy = null;
    Ds_FDxyerr = null;
    Ds_FDxysig = null;
    Ds_FD = null;
    Ds_FDerr = null;
    Ds_FDsig = null;
    dira = null;
    dira_angle = null;

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

    dR_Kp_mu = null;
    dR_Km_mu = null;
    dR_pi_mu = null;
    dR_phi_mu = null;
    dR_Ds_mu = null;

    mu_ip = null;
    mu_iperr = null;
    mu_ipchi2 = null;

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
}

void MatchRecoTree::CreateBranches()
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

    tree->Branch("phi_charge", &phi_charge);
    tree->Branch("phi_eta", &phi_eta);
    tree->Branch("phi_phi", &phi_phi);
    tree->Branch("phi_px", &phi_px);
    tree->Branch("phi_py", &phi_py);
    tree->Branch("phi_pz", &phi_pz);
    tree->Branch("phi_pt", &phi_pt);
    tree->Branch("phi_p", &phi_p);
    tree->Branch("phi_m", &phi_m);

    tree->Branch("Ds_charge", &Ds_charge);
    tree->Branch("Ds_eta", &Ds_eta);
    tree->Branch("Ds_phi", &Ds_phi);
    tree->Branch("Ds_px", &Ds_px);
    tree->Branch("Ds_py", &Ds_py);
    tree->Branch("Ds_pz", &Ds_pz);
    tree->Branch("Ds_pt", &Ds_pt);
    tree->Branch("Ds_p", &Ds_p);
    tree->Branch("Ds_m", &Ds_m);

    tree->Branch("phiFit_isValid", &phiFit_isValid);
    tree->Branch("phiFit_isFake", &phiFit_isFake);
    tree->Branch("phiFit_chi2", &phiFit_chi2);
    tree->Branch("phiFit_ndof", &phiFit_ndof);
    tree->Branch("phiFit_chi2ndof", &phiFit_chi2ndof);
    tree->Branch("phiFit_vx", &phiFit_vx);
    tree->Branch("phiFit_vy", &phiFit_vy);
    tree->Branch("phiFit_vz", &phiFit_vz);
    tree->Branch("phiFit_vxerr", &phiFit_vxerr);
    tree->Branch("phiFit_vyerr", &phiFit_vyerr);
    tree->Branch("phiFit_vzerr", &phiFit_vzerr);

    tree->Branch("DsFit_isValid", &DsFit_isValid);
    tree->Branch("DsFit_isFake", &DsFit_isFake);
    tree->Branch("DsFit_chi2", &DsFit_chi2);
    tree->Branch("DsFit_ndof", &DsFit_ndof);
    tree->Branch("DsFit_chi2ndof", &DsFit_chi2ndof);
    tree->Branch("DsFit_vx", &DsFit_vx);
    tree->Branch("DsFit_vy", &DsFit_vy);
    tree->Branch("DsFit_vz", &DsFit_vz);
    tree->Branch("DsFit_vxerr", &DsFit_vxerr);
    tree->Branch("DsFit_vyerr", &DsFit_vyerr);
    tree->Branch("DsFit_vzerr", &DsFit_vzerr);

    tree->Branch("orig_dR_Kp_Km", &orig_dR_Kp_Km);
    tree->Branch("orig_dR_Kp_pi", &orig_dR_Kp_pi);
    tree->Branch("orig_dR_Km_pi", &orig_dR_Km_pi);

    tree->Branch("phiFit_dR_Kp_Km", &phiFit_dR_Kp_Km);
    tree->Branch("phiFit_dR_Kp_pi", &phiFit_dR_Kp_pi);
    tree->Branch("phiFit_dR_Kp_phi", &phiFit_dR_Kp_phi);
    tree->Branch("phiFit_dR_Km_pi", &phiFit_dR_Km_pi);
    tree->Branch("phiFit_dR_Km_phi", &phiFit_dR_Km_phi);
    tree->Branch("phiFit_dR_pi_phi", &phiFit_dR_pi_phi);

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
    tree->Branch("dxy_phi_Ds", &dxy_phi_Ds);
    tree->Branch("dz_phi_Ds", &dz_phi_Ds);

    tree->Branch("Kp_ip", &Kp_ip);
    tree->Branch("Kp_iperr", &Kp_iperr);
    tree->Branch("Kp_ipchi2", &Kp_ipchi2);
    tree->Branch("Km_ip", &Km_ip);
    tree->Branch("Km_iperr", &Km_iperr);
    tree->Branch("Km_ipchi2", &Km_ipchi2);
    tree->Branch("pi_ip", &pi_ip);
    tree->Branch("pi_iperr", &pi_iperr);
    tree->Branch("pi_ipchi2", &pi_ipchi2);
    tree->Branch("Ds_FDxy", &Ds_FDxy);
    tree->Branch("Ds_FDxyerr", &Ds_FDxyerr);
    tree->Branch("Ds_FDxysig", &Ds_FDxysig);
    tree->Branch("Ds_FD", &Ds_FD);
    tree->Branch("Ds_FDerr", &Ds_FDerr);
    tree->Branch("Ds_FDsig", &Ds_FDsig);
    tree->Branch("dira", &dira);
    tree->Branch("dira_angle", &dira_angle);

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

    tree->Branch("dR_Kp_mu", &dR_Kp_mu);
    tree->Branch("dR_Km_mu", &dR_Km_mu);
    tree->Branch("dR_pi_mu", &dR_pi_mu);
    tree->Branch("dR_phi_mu", &dR_phi_mu);
    tree->Branch("dR_Ds_mu", &dR_Ds_mu);

    tree->Branch("mu_ip", &mu_ip);
    tree->Branch("mu_iperr", &mu_iperr);
    tree->Branch("mu_ipchi2", &mu_ipchi2);

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
}

#include "EDAnalyzers/Reconstruction/interface/RecoTree.h"
#include <iostream>

RecoTree::RecoTree(TTree *tree_)
{
    tree = tree_;
}

void RecoTree::Init()
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

    K1_charge.clear();
    K1_eta.clear();
    K1_phi.clear();
    K1_px.clear();
    K1_py.clear();
    K1_pz.clear();
    K1_pt.clear();
    K1_p.clear();

    K2_charge.clear();
    K2_eta.clear();
    K2_phi.clear();
    K2_px.clear();
    K2_py.clear();
    K2_pz.clear();
    K2_pt.clear();
    K2_p.clear();

    pi_charge.clear();
    pi_eta.clear();
    pi_phi.clear();
    pi_px.clear();
    pi_py.clear();
    pi_pz.clear();
    pi_pt.clear();
    pi_p.clear();

    phi_charge.clear();
    phi_eta.clear();
    phi_phi.clear();
    phi_px.clear();
    phi_py.clear();
    phi_pz.clear();
    phi_pt.clear();
    phi_p.clear();
    phi_m.clear();

    Ds_charge.clear();
    Ds_eta.clear();
    Ds_phi.clear();
    Ds_px.clear();
    Ds_py.clear();
    Ds_pz.clear();
    Ds_pt.clear();
    Ds_p.clear();
    Ds_m.clear();
    Ds_m_phiM.clear();

    mu_charge.clear();
    mu_eta.clear();
    mu_phi.clear();
    mu_px.clear();
    mu_py.clear();
    mu_pz.clear();
    mu_pt.clear();
    mu_p.clear();
    mu_isHighPt.clear();
    mu_isLoose.clear();
    mu_isMedium.clear();
    mu_isSoft.clear();
    mu_isTight.clear();
    mu_isPF.clear();
    mu_isTracker.clear();
    mu_isGlobal.clear();

    dR_K1_K2.clear();
    dR_K1_pi.clear();
    dR_K1_phi.clear();
    dR_K1_Ds.clear();
    dR_K1_mu.clear();
    dR_K2_pi.clear();
    dR_K2_phi.clear();
    dR_K2_Ds.clear();
    dR_K2_mu.clear();
    dR_pi_phi.clear();
    dR_pi_Ds.clear();
    dR_pi_mu.clear();
    dR_phi_Ds.clear();
    dR_phi_mu.clear();
    dR_Ds_mu.clear();

    phiFit_isValid.clear();
    phiFit_isFake.clear();
    phiFit_chi2.clear();
    phiFit_ndof.clear();
    phiFit_chi2ndof.clear();
    phiFit_vx.clear();
    phiFit_vy.clear();
    phiFit_vz.clear();
    phiFit_vxerr.clear();
    phiFit_vyerr.clear();
    phiFit_vzerr.clear();

    DsFit_isValid.clear();
    DsFit_isFake.clear();
    DsFit_chi2.clear();
    DsFit_ndof.clear();
    DsFit_chi2ndof.clear();
    DsFit_vx.clear();
    DsFit_vy.clear();
    DsFit_vz.clear();
    DsFit_vxerr.clear();
    DsFit_vyerr.clear();
    DsFit_vzerr.clear();

    dxy_phi_Ds.clear();
    dz_phi_Ds.clear();

    K1_ip.clear();
    K1_iperr.clear();
    K1_ipchi2.clear();
    K2_ip.clear();
    K2_iperr.clear();
    K2_ipchi2.clear();
    pi_ip.clear();
    pi_iperr.clear();
    pi_ipchi2.clear();
    mu_ip.clear();
    mu_iperr.clear();
    mu_ipchi2.clear();

    dira.clear();
    dira_angle.clear();

    Ds_FDxy.clear();
    Ds_FDxyerr.clear();
    Ds_FDxysig.clear();

    Ds_FD.clear();
    Ds_FDerr.clear();
    Ds_FDsig.clear();

    Ds_IsoR03_sumChargedHadronPt.clear();
    Ds_IsoR03_sumNeutralHadronEt.clear();
    Ds_IsoR03_sumPhotonEt.clear();
    Ds_IsoR03_sumPUPt.clear();
    Ds_IsoR03_PFIso.clear();
    Ds_IsoR03_RelIso.clear();

    Ds_IsoR04_sumChargedHadronPt.clear();
    Ds_IsoR04_sumNeutralHadronEt.clear();
    Ds_IsoR04_sumPhotonEt.clear();
    Ds_IsoR04_sumPUPt.clear();
    Ds_IsoR04_PFIso.clear();
    Ds_IsoR04_RelIso.clear();

    mu_IsoR03_sumChargedHadronPt.clear();
    mu_IsoR03_sumNeutralHadronEt.clear();
    mu_IsoR03_sumPhotonEt.clear();
    mu_IsoR03_sumPUPt.clear();
    mu_IsoR03_PFIso.clear();
    mu_IsoR03_RelIso.clear();

    mu_IsoR04_sumChargedHadronPt.clear();
    mu_IsoR04_sumNeutralHadronEt.clear();
    mu_IsoR04_sumPhotonEt.clear();
    mu_IsoR04_sumPUPt.clear();
    mu_IsoR04_PFIso.clear();
    mu_IsoR04_RelIso.clear();

    Kp_match.clear();
    Km_match.clear();
    pi_match.clear();
    mu_match.clear();
}

void RecoTree::CreateBranches()
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

    tree->Branch("K1_charge", &K1_charge);
    tree->Branch("K1_eta", &K1_eta);
    tree->Branch("K1_phi", &K1_phi);
    tree->Branch("K1_px", &K1_px);
    tree->Branch("K1_py", &K1_py);
    tree->Branch("K1_pz", &K1_pz);
    tree->Branch("K1_pt", &K1_pt);
    tree->Branch("K1_p", &K1_p);

    tree->Branch("K2_charge", &K2_charge);
    tree->Branch("K2_eta", &K2_eta);
    tree->Branch("K2_phi", &K2_phi);
    tree->Branch("K2_px", &K2_px);
    tree->Branch("K2_py", &K2_py);
    tree->Branch("K2_pz", &K2_pz);
    tree->Branch("K2_pt", &K2_pt);
    tree->Branch("K2_p", &K2_p);

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
    tree->Branch("Ds_m_phiM", &Ds_m_phiM);

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

    tree->Branch("dR_K1_K2", &dR_K1_K2);
    tree->Branch("dR_K1_pi", &dR_K1_pi);
    tree->Branch("dR_K1_phi", &dR_K1_phi);
    tree->Branch("dR_K1_Ds", &dR_K1_Ds);
    tree->Branch("dR_K1_mu", &dR_K1_mu);
    tree->Branch("dR_K2_pi", &dR_K2_pi);
    tree->Branch("dR_K2_phi", &dR_K2_phi);
    tree->Branch("dR_K2_Ds", &dR_K2_Ds);
    tree->Branch("dR_K2_mu", &dR_K2_mu);
    tree->Branch("dR_pi_phi", &dR_pi_phi);
    tree->Branch("dR_pi_Ds", &dR_pi_Ds);
    tree->Branch("dR_pi_mu", &dR_pi_mu);
    tree->Branch("dR_phi_Ds", &dR_phi_Ds);
    tree->Branch("dR_phi_mu", &dR_phi_mu);
    tree->Branch("dR_Ds_mu", &dR_Ds_mu);

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

    tree->Branch("dxy_phi_Ds", &dxy_phi_Ds);
    tree->Branch("dz_phi_Ds", &dz_phi_Ds);

    tree->Branch("K1_ip", &K1_ip);
    tree->Branch("K1_iperr", &K1_iperr);
    tree->Branch("K1_ipchi2", &K1_ipchi2);
    tree->Branch("K2_ip", &K2_ip);
    tree->Branch("K2_iperr", &K2_iperr);
    tree->Branch("K2_ipchi2", &K2_ipchi2);
    tree->Branch("pi_ip", &pi_ip);
    tree->Branch("pi_iperr", &pi_iperr);
    tree->Branch("pi_ipchi2", &pi_ipchi2);
    tree->Branch("mu_ip", &mu_ip);
    tree->Branch("mu_iperr", &mu_iperr);
    tree->Branch("mu_ipchi2", &mu_ipchi2);

    tree->Branch("dira", &dira);
    tree->Branch("dira_angle", &dira_angle);

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

    tree->Branch("Kp_match", &Kp_match);
    tree->Branch("Km_match", &Km_match);
    tree->Branch("pi_match", &pi_match);
    tree->Branch("mu_match", &mu_match);
}

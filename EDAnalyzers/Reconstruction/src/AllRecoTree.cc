#include "EDAnalyzers/Reconstruction/interface/AllRecoTree.h"
#include <iostream>

AllRecoTree::AllRecoTree(TTree *tree_)
{
    tree = tree_;
}

void AllRecoTree::Init()
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

    Kp_charge.clear();
    Kp_eta.clear();
    Kp_phi.clear();
    Kp_px.clear();
    Kp_py.clear();
    Kp_pz.clear();
    Kp_pt.clear();
    Kp_p.clear();

    Km_charge.clear();
    Km_eta.clear();
    Km_phi.clear();
    Km_px.clear();
    Km_py.clear();
    Km_pz.clear();
    Km_pt.clear();
    Km_p.clear();

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


    orig_dR_Kp_Km.clear();
    orig_dR_Kp_pi.clear();
    orig_dR_Km_pi.clear();

    phiFit_dR_Kp_Km.clear();
    phiFit_dR_Kp_pi.clear();
    phiFit_dR_Kp_phi.clear();
    phiFit_dR_Km_pi.clear();
    phiFit_dR_Km_phi.clear();
    phiFit_dR_pi_phi.clear();

    dR_Kp_Km.clear();
    dR_Kp_pi.clear();
    dR_Kp_phi.clear();
    dR_Kp_Ds.clear();
    dR_Km_pi.clear();
    dR_Km_phi.clear();
    dR_Km_Ds.clear();
    dR_pi_phi.clear();
    dR_pi_Ds.clear();
    dR_phi_Ds.clear();

    dxy_phi_Ds.clear();
    dz_phi_Ds.clear();

    Kp_ip.clear();
    Kp_iperr.clear();
    Kp_ipchi2.clear();
    Km_ip.clear();
    Km_iperr.clear();
    Km_ipchi2.clear();
    pi_ip.clear();
    pi_iperr.clear();
    pi_ipchi2.clear();

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

    Kp_match.clear();
    Km_match.clear();
    pi_match.clear();

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

    mu_ip.clear();
    mu_iperr.clear();
    mu_ipchi2.clear();

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

    mu_match.clear();
}

void AllRecoTree::CreateBranches()
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
    tree->Branch("Ds_m_phiM", &Ds_m_phiM);

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

    tree->Branch("Kp_match", &Kp_match);
    tree->Branch("Km_match", &Km_match);
    tree->Branch("pi_match", &pi_match);

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

    tree->Branch("mu_match", &mu_match);
}

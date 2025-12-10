#include "EDAnalyzers/RecoAnalyzer/interface/RecoTree.h"
#include <iostream>

RecoTree::RecoTree(TTree *tree_)
{
    tree = tree_;
}

void RecoTree::BSInfo::fillAll( const reco::BeamSpot& beamspot )
{
    type       = beamspot.type();
    x0         = beamspot.x0();
    y0         = beamspot.y0();
    z0         = beamspot.z0();
    sigmaZ     = beamspot.sigmaZ();
    dxdz       = beamspot.dxdz();
    dydz       = beamspot.dydz();
    BWX        = beamspot.BeamWidthX();
    BWY        = beamspot.BeamWidthY();
    x0err      = beamspot.x0Error();
    y0err      = beamspot.y0Error();
    z0err      = beamspot.z0Error();
    sigmaZ0err = beamspot.sigmaZ0Error();
    dxdzerr    = beamspot.dxdzError();
    dydzerr    = beamspot.dydzError();
    BWXerr     = beamspot.BeamWidthXError();
    BWYerr     = beamspot.BeamWidthYError();
    emitX      = beamspot.emittanceX();
    emitY      = beamspot.emittanceY();
    betaStar   = beamspot.betaStar();
}

void RecoTree::PFInfo::fillAll( const pat::PackedCandidate& pfCands )
{
    isIsolatedChargedHadron = pfCands.isIsolatedChargedHadron();
    charge                  = pfCands.charge();
    eta                     = pfCands.eta();
    phi                     = pfCands.phi();
    vx                      = pfCands.vx();
    vy                      = pfCands.vy();
    vz                      = pfCands.vz();
    p                       = pfCands.p();
    pt                      = pfCands.pt();
    px                      = pfCands.px();
    py                      = pfCands.py();
    pz                      = pfCands.pz();
}

void RecoTree::TLVInfo::fillAll( const TLorentzVector& TLV )
{
    eta  = TLV.Eta();
    phi  = TLV.Phi();
    p    = TLV.P();
    pt   = TLV.Pt();
    px   = TLV.Px();
    py   = TLV.Py();
    pz   = TLV.Pz();
    invm = TLV.M();
}

void RecoTree::VtxInfo::fillAll( const reco::Vertex& Vtx )
{
    isValid  = Vtx.isValid();
    isFake   = Vtx.isFake();
    chi2     = Vtx.chi2();
    ndof     = Vtx.ndof();
    chi2ndof = Vtx.chi2()/Vtx.ndof();
    vx       = Vtx.x();
    vy       = Vtx.y();
    vz       = Vtx.z();
    vxerr    = Vtx.xError();
    vyerr    = Vtx.yError();
    vzerr    = Vtx.zError();
}

void RecoTree::FDInfo::fillAll( const GlobalVector& Pos, const GlobalError& Cov, const GlobalVector& PDirec )
{
    AlgebraicSymMatrix33 cov_mtrx = Cov.matrix();
    
    FDxy = std::hypot(Pos.x(), Pos.y());
    FDxyerr = std::sqrt( Pos.x()*Pos.x()*cov_mtrx(0,0) + Pos.y()*Pos.y()*cov_mtrx(1,1) + 2*Pos.x()*Pos.y()*cov_mtrx(0,1) ) / FDxy;
    Measurement1D FDxy_withErr(FDxy, FDxyerr);
    FDxychi2 = FDxy_withErr.significance();

    FDzerr = std::sqrt( cov_mtrx(2,2) );
    FDz = std::abs(Pos.z());
    Measurement1D FDz_withErr(FDz, FDzerr);
    FDzchi2 = FDz_withErr.significance();

    FD = std::hypot(FDxy, FDz);
    FDerr = std::sqrt( Pos.x()*Pos.x()*cov_mtrx(0,0) + Pos.y()*Pos.y()*cov_mtrx(1,1) + Pos.z()*Pos.z()*cov_mtrx(2,2) + 2*Pos.x()*Pos.y()*cov_mtrx(0,1) + 2*Pos.x()*Pos.z()*cov_mtrx(0,2) + 2*Pos.y()*Pos.z()*cov_mtrx(1,2)) / FD;
    Measurement1D FD_withErr(FD, FDerr);
    FDchi2 = FD_withErr.significance();

    dira = Pos.dot(PDirec) / ( Pos.mag() * PDirec.mag() );
    dira_angle = std::acos(dira);
}

void RecoTree::IPInfo::fillAll( const reco::TransientTrack& ttrack, const reco::Vertex& vtx )
{
    std::pair<bool, Measurement1D> result = IPTools::absoluteImpactParameter3D(ttrack, vtx);
    if( result.first ){
        ip     = result.second.value();
        iperr  = result.second.error();
        ipchi2 = result.second.significance();
    }
}

void RecoTree::MuonInfo::fillAll( const pat::Muon& muon, const reco::Vertex & Vtx )
{
    charge    = muon.charge();
    eta       = muon.eta();
    phi       = muon.phi();
    vx        = muon.vx();
    vy        = muon.vy();
    vz        = muon.vz();
    p         = muon.p();
    pt        = muon.pt();
    px        = muon.px();
    py        = muon.py();
    pz        = muon.pz();
    dxy       = muon.muonBestTrack()->dxy(Vtx.position());
    dxyerr    = muon.muonBestTrack()->dxyError();
    dz        = muon.muonBestTrack()->dz(Vtx.position());
    dzerr     = muon.muonBestTrack()->dzError();
    isHighPt  = muon.isHighPtMuon(Vtx);
    isLoose   = muon.isLooseMuon();
    isMedium  = muon.isMediumMuon();
    isSoft    = muon.isSoftMuon(Vtx);
    isTight   = muon.isTightMuon(Vtx);
    isPF      = muon.isPFMuon();
    isTracker = muon.isTrackerMuon();
    isGlobal  = muon.isGlobalMuon();
}


void RecoTree::Init()
{
    num_reco_phi = 0;
    num_reco_Ds = 0;

    Kp_vec.clearAll();
    Km_vec.clearAll();
    pi_vec.clearAll();
    phi_vec.clearAll();
    Ds_vec.clearAll();
    phiFit_Kp_vec.clearAll();
    phiFit_Km_vec.clearAll();
    phiFit_pi_vec.clearAll();
    phiFit_phi_vec.clearAll();
    phiFit_Ds_vec.clearAll();
    DsFit_Kp_vec.clearAll();
    DsFit_Km_vec.clearAll();
    DsFit_pi_vec.clearAll();
    DsFit_phi_vec.clearAll();
    DsFit_Ds_vec.clearAll();
   
    dR_Kp_Km_vec.clear();
    dR_Kp_phi_vec.clear();
    dR_Km_phi_vec.clear();
    dR_Kp_pi_vec.clear();
    dR_Km_pi_vec.clear();
    dR_pi_phi_vec.clear();
    dR_Kp_Ds_vec.clear();
    dR_Km_Ds_vec.clear();
    dR_phi_Ds_vec.clear();
    dR_pi_Ds_vec.clear();

    phiFit_dR_Kp_Km_vec.clear();
    phiFit_dR_Kp_phi_vec.clear();
    phiFit_dR_Km_phi_vec.clear();
    phiFit_dR_Kp_pi_vec.clear();
    phiFit_dR_Km_pi_vec.clear();
    phiFit_dR_pi_phi_vec.clear();
    phiFit_dR_Kp_Ds_vec.clear();
    phiFit_dR_Km_Ds_vec.clear();
    phiFit_dR_phi_Ds_vec.clear();
    phiFit_dR_pi_Ds_vec.clear();

    DsFit_dR_Kp_Km_vec.clear();
    DsFit_dR_Kp_phi_vec.clear();
    DsFit_dR_Km_phi_vec.clear();
    DsFit_dR_Kp_pi_vec.clear();
    DsFit_dR_Km_pi_vec.clear();
    DsFit_dR_pi_phi_vec.clear();
    DsFit_dR_Kp_Ds_vec.clear();
    DsFit_dR_Km_Ds_vec.clear();
    DsFit_dR_phi_Ds_vec.clear();
    DsFit_dR_pi_Ds_vec.clear();

    dxy_phi_Ds_vec.clear();
    dz_phi_Ds_vec.clear();

    DsFit_Mconstraint_Ds_invm_vec.clear();

    phiFit_vec.clearAll();
    DsFit_vec.clearAll();
    PVnoDs_vec.clearAll();
    PVwithDs_vec.clearAll();

    Ds_primvtx_FD_vec.clearAll();
    Kp_primvtx_ip_vec.clearAll();
    Km_primvtx_ip_vec.clearAll();
    pi_primvtx_ip_vec.clearAll();
    phi_primvtx_ip_vec.clearAll();
    Ds_primvtx_ip_vec.clearAll();

    Ds_PVnoDs_FD_vec.clearAll();
    Kp_PVnoDs_ip_vec.clearAll();
    Km_PVnoDs_ip_vec.clearAll();
    pi_PVnoDs_ip_vec.clearAll();
    phi_PVnoDs_ip_vec.clearAll();
    Ds_PVnoDs_ip_vec.clearAll();

    Ds_PVwithDs_FD_vec.clearAll();
    Kp_PVwithDs_ip_vec.clearAll();
    Km_PVwithDs_ip_vec.clearAll();
    pi_PVwithDs_ip_vec.clearAll();
    phi_PVwithDs_ip_vec.clearAll();
    Ds_PVwithDs_ip_vec.clearAll();

    Ds_IsoR03_vec.clearAll();
    Ds_IsoR04_vec.clearAll();

    best_Kp_vec.clearAll();
    best_Km_vec.clearAll();
    best_pi_vec.clearAll();
    best_phi_vec.clearAll();
    best_Ds_vec.clearAll();
    best_phiFit_Kp_vec.clearAll();
    best_phiFit_Km_vec.clearAll();
    best_phiFit_pi_vec.clearAll();
    best_phiFit_phi_vec.clearAll();
    best_phiFit_Ds_vec.clearAll();
    best_DsFit_Kp_vec.clearAll();
    best_DsFit_Km_vec.clearAll();
    best_DsFit_pi_vec.clearAll();
    best_DsFit_phi_vec.clearAll();
    best_DsFit_Ds_vec.clearAll();

    best_dR_Kp_Km_vec.clear();
    best_dR_Kp_phi_vec.clear();
    best_dR_Km_phi_vec.clear();
    best_dR_Kp_pi_vec.clear();
    best_dR_Km_pi_vec.clear();
    best_dR_pi_phi_vec.clear();
    best_dR_Kp_Ds_vec.clear();
    best_dR_Km_Ds_vec.clear();
    best_dR_phi_Ds_vec.clear();
    best_dR_pi_Ds_vec.clear();

    best_phiFit_dR_Kp_Km_vec.clear();
    best_phiFit_dR_Kp_phi_vec.clear();
    best_phiFit_dR_Km_phi_vec.clear();
    best_phiFit_dR_Kp_pi_vec.clear();
    best_phiFit_dR_Km_pi_vec.clear();
    best_phiFit_dR_pi_phi_vec.clear();
    best_phiFit_dR_Kp_Ds_vec.clear();
    best_phiFit_dR_Km_Ds_vec.clear();
    best_phiFit_dR_phi_Ds_vec.clear();
    best_phiFit_dR_pi_Ds_vec.clear();

    best_DsFit_dR_Kp_Km_vec.clear();
    best_DsFit_dR_Kp_phi_vec.clear();
    best_DsFit_dR_Km_phi_vec.clear();
    best_DsFit_dR_Kp_pi_vec.clear();
    best_DsFit_dR_Km_pi_vec.clear();
    best_DsFit_dR_pi_phi_vec.clear();
    best_DsFit_dR_Kp_Ds_vec.clear();
    best_DsFit_dR_Km_Ds_vec.clear();
    best_DsFit_dR_phi_Ds_vec.clear();
    best_DsFit_dR_pi_Ds_vec.clear();

    best_dxy_phi_Ds_vec.clear();
    best_dz_phi_Ds_vec.clear();

    best_DsFit_Mconstraint_Ds_invm_vec.clear();

    best_phiFit_vec.clearAll();
    best_DsFit_vec.clearAll();
    best_PVnoDs_vec.clearAll();
    best_PVwithDs_vec.clearAll();

    best_Ds_primvtx_FD_vec.clearAll();
    best_Kp_primvtx_ip_vec.clearAll();
    best_Km_primvtx_ip_vec.clearAll();
    best_pi_primvtx_ip_vec.clearAll();
    best_phi_primvtx_ip_vec.clearAll();
    best_Ds_primvtx_ip_vec.clearAll();

    best_Ds_PVnoDs_FD_vec.clearAll();
    best_Kp_PVnoDs_ip_vec.clearAll();
    best_Km_PVnoDs_ip_vec.clearAll();
    best_pi_PVnoDs_ip_vec.clearAll();
    best_phi_PVnoDs_ip_vec.clearAll();
    best_Ds_PVnoDs_ip_vec.clearAll();

    best_Ds_PVwithDs_FD_vec.clearAll();
    best_Kp_PVwithDs_ip_vec.clearAll();
    best_Km_PVwithDs_ip_vec.clearAll();
    best_pi_PVwithDs_ip_vec.clearAll();
    best_phi_PVwithDs_ip_vec.clearAll();
    best_Ds_PVwithDs_ip_vec.clearAll();

    best_Ds_IsoR03_vec.clearAll();
    best_Ds_IsoR04_vec.clearAll();
 
    mu_vec.clearAll();
    mu_IsoR03_vec.clearAll();
    mu_IsoR04_vec.clearAll();
    mu_primvtx_ip_vec.clearAll();
}

void RecoTree::CreateBranches()
{
    tree->Branch("beamspot_type",       &beamspotInfo.type,       "beamspot_type/I");
    tree->Branch("beamspot_x0",         &beamspotInfo.x0,         "beamspot_x0/D");
    tree->Branch("beamspot_y0",         &beamspotInfo.y0,         "beamspot_y0/D");
    tree->Branch("beamspot_z0",         &beamspotInfo.z0,         "beamspot_z0/D");
    tree->Branch("beamspot_sigmaZ",     &beamspotInfo.sigmaZ,     "beamspot_sigmaZ/D");
    tree->Branch("beamspot_dxdz",       &beamspotInfo.dxdz,       "beamspot_dxdz/D");
    tree->Branch("beamspot_dydz",       &beamspotInfo.dydz,       "beamspot_dydz/D");
    tree->Branch("beamspot_BWX",        &beamspotInfo.BWX,        "beamspot_BWX/D");
    tree->Branch("beamspot_BWY",        &beamspotInfo.BWY,        "beamspot_BWY/D");
    tree->Branch("beamspot_x0err",      &beamspotInfo.x0err,      "beamspot_x0err/D");
    tree->Branch("beamspot_y0err",      &beamspotInfo.y0err,      "beamspot_y0err/D");
    tree->Branch("beamspot_z0err",      &beamspotInfo.z0err,      "beamspot_z0err/D");
    tree->Branch("beamspot_sigmaZ0err", &beamspotInfo.sigmaZ0err, "beamspot_sigmaZ0err/D");
    tree->Branch("beamspot_dxdzerr",    &beamspotInfo.dxdzerr,    "beamspot_dxdzerr/D");
    tree->Branch("beamspot_dydzerr",    &beamspotInfo.dydzerr,    "beamspot_dydzerr/D");
    tree->Branch("beamspot_BWXerr",     &beamspotInfo.BWXerr,     "beamspot_BWXerr/D");
    tree->Branch("beamspot_BWYerr",     &beamspotInfo.BWYerr,     "beamspot_BWYerr/D");
    tree->Branch("beamspot_emitX",      &beamspotInfo.emitX,      "beamspot_emitX/D");
    tree->Branch("beamspot_emitY",      &beamspotInfo.emitY,      "beamspot_emitY/D");
    tree->Branch("beamspot_betaStar",   &beamspotInfo.betaStar,   "beamspot_betaStar/D");
    
    tree->Branch("primvtx_isValid",  &primvtxInfo.isValid,  "primvtx_isValid/O");
    tree->Branch("primvtx_isFake",   &primvtxInfo.isFake,   "primvtx_isFake/O");
    tree->Branch("primvtx_chi2",     &primvtxInfo.chi2,     "primvtx_chi2/D");
    tree->Branch("primvtx_ndof",     &primvtxInfo.ndof,     "primvtx_ndof/D");
    tree->Branch("primvtx_chi2ndof", &primvtxInfo.chi2ndof, "primvtx_chi2ndof/D");
    tree->Branch("primvtx_vx",       &primvtxInfo.vx,       "primvtx_vx/D");
    tree->Branch("primvtx_vy",       &primvtxInfo.vy,       "primvtx_vy/D");
    tree->Branch("primvtx_vz",       &primvtxInfo.vz,       "primvtx_vz/D");
    tree->Branch("primvtx_vxerr",    &primvtxInfo.vxerr,    "primvtx_vxerr/D");
    tree->Branch("primvtx_vyerr",    &primvtxInfo.vyerr,    "primvtx_vyerr/D");
    tree->Branch("primvtx_vzerr",    &primvtxInfo.vzerr,    "primvtx_vzerr/D");

    tree->Branch("best_Kp_isIsolatedChargedHadron", &best_Kp_vec.isIsolatedChargedHadron);
    tree->Branch("best_Kp_charge",                  &best_Kp_vec.charge);
    tree->Branch("best_Kp_eta",                     &best_Kp_vec.eta);
    tree->Branch("best_Kp_phi",                     &best_Kp_vec.phi);
    tree->Branch("best_Kp_vx",                      &best_Kp_vec.vx);
    tree->Branch("best_Kp_vy",                      &best_Kp_vec.vy);
    tree->Branch("best_Kp_vz",                      &best_Kp_vec.vz);
    tree->Branch("best_Kp_p",                       &best_Kp_vec.p);
    tree->Branch("best_Kp_pt",                      &best_Kp_vec.pt);
    tree->Branch("best_Kp_px",                      &best_Kp_vec.px);
    tree->Branch("best_Kp_py",                      &best_Kp_vec.py);
    tree->Branch("best_Kp_pz",                      &best_Kp_vec.pz);

    tree->Branch("best_Km_isIsolatedChargedHadron", &best_Km_vec.isIsolatedChargedHadron);
    tree->Branch("best_Km_charge",                  &best_Km_vec.charge);
    tree->Branch("best_Km_eta",                     &best_Km_vec.eta);
    tree->Branch("best_Km_phi",                     &best_Km_vec.phi);
    tree->Branch("best_Km_vx",                      &best_Km_vec.vx);
    tree->Branch("best_Km_vy",                      &best_Km_vec.vy);
    tree->Branch("best_Km_vz",                      &best_Km_vec.vz);
    tree->Branch("best_Km_p",                       &best_Km_vec.p);
    tree->Branch("best_Km_pt",                      &best_Km_vec.pt);
    tree->Branch("best_Km_px",                      &best_Km_vec.px);
    tree->Branch("best_Km_py",                      &best_Km_vec.py);
    tree->Branch("best_Km_pz",                      &best_Km_vec.pz);

    tree->Branch("best_pi_isIsolatedChargedHadron", &best_pi_vec.isIsolatedChargedHadron);
    tree->Branch("best_pi_charge",                  &best_pi_vec.charge);
    tree->Branch("best_pi_eta",                     &best_pi_vec.eta);
    tree->Branch("best_pi_phi",                     &best_pi_vec.phi);
    tree->Branch("best_pi_vx",                      &best_pi_vec.vx);
    tree->Branch("best_pi_vy",                      &best_pi_vec.vy);
    tree->Branch("best_pi_vz",                      &best_pi_vec.vz);
    tree->Branch("best_pi_p",                       &best_pi_vec.p);
    tree->Branch("best_pi_pt",                      &best_pi_vec.pt);
    tree->Branch("best_pi_px",                      &best_pi_vec.px);
    tree->Branch("best_pi_py",                      &best_pi_vec.py);
    tree->Branch("best_pi_pz",                      &best_pi_vec.pz);

    tree->Branch("best_phi_eta",  &best_phi_vec.eta);
    tree->Branch("best_phi_phi",  &best_phi_vec.phi);
    tree->Branch("best_phi_p",    &best_phi_vec.p);
    tree->Branch("best_phi_pt",   &best_phi_vec.pt);
    tree->Branch("best_phi_px",   &best_phi_vec.px);
    tree->Branch("best_phi_py",   &best_phi_vec.py);
    tree->Branch("best_phi_pz",   &best_phi_vec.pz);
    tree->Branch("best_phi_invm", &best_phi_vec.invm);

    tree->Branch("best_Ds_eta",  &best_Ds_vec.eta);
    tree->Branch("best_Ds_phi",  &best_Ds_vec.phi);
    tree->Branch("best_Ds_p",    &best_Ds_vec.p);
    tree->Branch("best_Ds_pt",   &best_Ds_vec.pt);
    tree->Branch("best_Ds_px",   &best_Ds_vec.px);
    tree->Branch("best_Ds_py",   &best_Ds_vec.py);
    tree->Branch("best_Ds_pz",   &best_Ds_vec.pz);
    tree->Branch("best_Ds_invm", &best_Ds_vec.invm);

    tree->Branch("best_phiFit_Kp_eta",  &best_phiFit_Kp_vec.eta);
    tree->Branch("best_phiFit_Kp_phi",  &best_phiFit_Kp_vec.phi);
    tree->Branch("best_phiFit_Kp_p",    &best_phiFit_Kp_vec.p);
    tree->Branch("best_phiFit_Kp_pt",   &best_phiFit_Kp_vec.pt);
    tree->Branch("best_phiFit_Kp_px",   &best_phiFit_Kp_vec.px);
    tree->Branch("best_phiFit_Kp_py",   &best_phiFit_Kp_vec.py);
    tree->Branch("best_phiFit_Kp_pz",   &best_phiFit_Kp_vec.pz);
    tree->Branch("best_phiFit_Kp_invm", &best_phiFit_Kp_vec.invm);

    tree->Branch("best_phiFit_Km_eta",  &best_phiFit_Km_vec.eta);
    tree->Branch("best_phiFit_Km_phi",  &best_phiFit_Km_vec.phi);
    tree->Branch("best_phiFit_Km_p",    &best_phiFit_Km_vec.p);
    tree->Branch("best_phiFit_Km_pt",   &best_phiFit_Km_vec.pt);
    tree->Branch("best_phiFit_Km_px",   &best_phiFit_Km_vec.px);
    tree->Branch("best_phiFit_Km_py",   &best_phiFit_Km_vec.py);
    tree->Branch("best_phiFit_Km_pz",   &best_phiFit_Km_vec.pz);
    tree->Branch("best_phiFit_Km_invm", &best_phiFit_Km_vec.invm);

    tree->Branch("best_phiFit_pi_eta",  &best_phiFit_pi_vec.eta);
    tree->Branch("best_phiFit_pi_phi",  &best_phiFit_pi_vec.phi);
    tree->Branch("best_phiFit_pi_p",    &best_phiFit_pi_vec.p);
    tree->Branch("best_phiFit_pi_pt",   &best_phiFit_pi_vec.pt);
    tree->Branch("best_phiFit_pi_px",   &best_phiFit_pi_vec.px);
    tree->Branch("best_phiFit_pi_py",   &best_phiFit_pi_vec.py);
    tree->Branch("best_phiFit_pi_pz",   &best_phiFit_pi_vec.pz);
    tree->Branch("best_phiFit_pi_invm", &best_phiFit_pi_vec.invm);

    tree->Branch("best_phiFit_phi_eta",  &best_phiFit_phi_vec.eta);
    tree->Branch("best_phiFit_phi_phi",  &best_phiFit_phi_vec.phi);
    tree->Branch("best_phiFit_phi_p",    &best_phiFit_phi_vec.p);
    tree->Branch("best_phiFit_phi_pt",   &best_phiFit_phi_vec.pt);
    tree->Branch("best_phiFit_phi_px",   &best_phiFit_phi_vec.px);
    tree->Branch("best_phiFit_phi_py",   &best_phiFit_phi_vec.py);
    tree->Branch("best_phiFit_phi_pz",   &best_phiFit_phi_vec.pz);
    tree->Branch("best_phiFit_phi_invm", &best_phiFit_phi_vec.invm);

    tree->Branch("best_phiFit_Ds_eta",  &best_phiFit_Ds_vec.eta);
    tree->Branch("best_phiFit_Ds_phi",  &best_phiFit_Ds_vec.phi);
    tree->Branch("best_phiFit_Ds_p",    &best_phiFit_Ds_vec.p);
    tree->Branch("best_phiFit_Ds_pt",   &best_phiFit_Ds_vec.pt);
    tree->Branch("best_phiFit_Ds_px",   &best_phiFit_Ds_vec.px);
    tree->Branch("best_phiFit_Ds_py",   &best_phiFit_Ds_vec.py);
    tree->Branch("best_phiFit_Ds_pz",   &best_phiFit_Ds_vec.pz);
    tree->Branch("best_phiFit_Ds_invm", &best_phiFit_Ds_vec.invm);

    tree->Branch("best_DsFit_Kp_eta",  &best_DsFit_Kp_vec.eta);
    tree->Branch("best_DsFit_Kp_phi",  &best_DsFit_Kp_vec.phi);
    tree->Branch("best_DsFit_Kp_p",    &best_DsFit_Kp_vec.p);
    tree->Branch("best_DsFit_Kp_pt",   &best_DsFit_Kp_vec.pt);
    tree->Branch("best_DsFit_Kp_px",   &best_DsFit_Kp_vec.px);
    tree->Branch("best_DsFit_Kp_py",   &best_DsFit_Kp_vec.py);
    tree->Branch("best_DsFit_Kp_pz",   &best_DsFit_Kp_vec.pz);
    tree->Branch("best_DsFit_Kp_invm", &best_DsFit_Kp_vec.invm);

    tree->Branch("best_DsFit_Km_eta",  &best_DsFit_Km_vec.eta);
    tree->Branch("best_DsFit_Km_phi",  &best_DsFit_Km_vec.phi);
    tree->Branch("best_DsFit_Km_p",    &best_DsFit_Km_vec.p);
    tree->Branch("best_DsFit_Km_pt",   &best_DsFit_Km_vec.pt);
    tree->Branch("best_DsFit_Km_px",   &best_DsFit_Km_vec.px);
    tree->Branch("best_DsFit_Km_py",   &best_DsFit_Km_vec.py);
    tree->Branch("best_DsFit_Km_pz",   &best_DsFit_Km_vec.pz);
    tree->Branch("best_DsFit_Km_invm", &best_DsFit_Km_vec.invm);

    tree->Branch("best_DsFit_pi_eta",  &best_DsFit_pi_vec.eta);
    tree->Branch("best_DsFit_pi_phi",  &best_DsFit_pi_vec.phi);
    tree->Branch("best_DsFit_pi_p",    &best_DsFit_pi_vec.p);
    tree->Branch("best_DsFit_pi_pt",   &best_DsFit_pi_vec.pt);
    tree->Branch("best_DsFit_pi_px",   &best_DsFit_pi_vec.px);
    tree->Branch("best_DsFit_pi_py",   &best_DsFit_pi_vec.py);
    tree->Branch("best_DsFit_pi_pz",   &best_DsFit_pi_vec.pz);
    tree->Branch("best_DsFit_pi_invm", &best_DsFit_pi_vec.invm);

    tree->Branch("best_DsFit_phi_eta",  &best_DsFit_phi_vec.eta);
    tree->Branch("best_DsFit_phi_phi",  &best_DsFit_phi_vec.phi);
    tree->Branch("best_DsFit_phi_p",    &best_DsFit_phi_vec.p);
    tree->Branch("best_DsFit_phi_pt",   &best_DsFit_phi_vec.pt);
    tree->Branch("best_DsFit_phi_px",   &best_DsFit_phi_vec.px);
    tree->Branch("best_DsFit_phi_py",   &best_DsFit_phi_vec.py);
    tree->Branch("best_DsFit_phi_pz",   &best_DsFit_phi_vec.pz);
    tree->Branch("best_DsFit_phi_invm", &best_DsFit_phi_vec.invm);

    tree->Branch("best_DsFit_Ds_eta",  &best_DsFit_Ds_vec.eta);
    tree->Branch("best_DsFit_Ds_phi",  &best_DsFit_Ds_vec.phi);
    tree->Branch("best_DsFit_Ds_p",    &best_DsFit_Ds_vec.p);
    tree->Branch("best_DsFit_Ds_pt",   &best_DsFit_Ds_vec.pt);
    tree->Branch("best_DsFit_Ds_px",   &best_DsFit_Ds_vec.px);
    tree->Branch("best_DsFit_Ds_py",   &best_DsFit_Ds_vec.py);
    tree->Branch("best_DsFit_Ds_pz",   &best_DsFit_Ds_vec.pz);
    tree->Branch("best_DsFit_Ds_invm", &best_DsFit_Ds_vec.invm);

    tree->Branch("best_dR_Kp_Km",  &best_dR_Kp_Km_vec);
    tree->Branch("best_dR_Kp_phi", &best_dR_Kp_phi_vec);
    tree->Branch("best_dR_Km_phi", &best_dR_Km_phi_vec);
    tree->Branch("best_dR_Kp_pi",  &best_dR_Kp_pi_vec);
    tree->Branch("best_dR_Km_pi",  &best_dR_Km_pi_vec);
    tree->Branch("best_dR_pi_phi", &best_dR_pi_phi_vec);
    tree->Branch("best_dR_Kp_Ds",  &best_dR_Kp_Ds_vec);
    tree->Branch("best_dR_Km_Ds",  &best_dR_Km_Ds_vec);
    tree->Branch("best_dR_phi_Ds", &best_dR_phi_Ds_vec);
    tree->Branch("best_dR_pi_Ds",  &best_dR_pi_Ds_vec);
    
    tree->Branch("best_phiFit_dR_Kp_Km",  &best_phiFit_dR_Kp_Km_vec);
    tree->Branch("best_phiFit_dR_Kp_phi", &best_phiFit_dR_Kp_phi_vec);
    tree->Branch("best_phiFit_dR_Km_phi", &best_phiFit_dR_Km_phi_vec);
    tree->Branch("best_phiFit_dR_Kp_pi",  &best_phiFit_dR_Kp_pi_vec);
    tree->Branch("best_phiFit_dR_Km_pi",  &best_phiFit_dR_Km_pi_vec);
    tree->Branch("best_phiFit_dR_pi_phi", &best_phiFit_dR_pi_phi_vec);
    tree->Branch("best_phiFit_dR_Kp_Ds",  &best_phiFit_dR_Kp_Ds_vec);
    tree->Branch("best_phiFit_dR_Km_Ds",  &best_phiFit_dR_Km_Ds_vec);
    tree->Branch("best_phiFit_dR_phi_Ds", &best_phiFit_dR_phi_Ds_vec);
    tree->Branch("best_phiFit_dR_pi_Ds",  &best_phiFit_dR_pi_Ds_vec);


    tree->Branch("best_DsFit_dR_Kp_Km",  &best_DsFit_dR_Kp_Km_vec);
    tree->Branch("best_DsFit_dR_Kp_phi", &best_DsFit_dR_Kp_phi_vec);
    tree->Branch("best_DsFit_dR_Km_phi", &best_DsFit_dR_Km_phi_vec);
    tree->Branch("best_DsFit_dR_Kp_pi",  &best_DsFit_dR_Kp_pi_vec);
    tree->Branch("best_DsFit_dR_Km_pi",  &best_DsFit_dR_Km_pi_vec);
    tree->Branch("best_DsFit_dR_pi_phi", &best_DsFit_dR_pi_phi_vec);
    tree->Branch("best_DsFit_dR_Kp_Ds",  &best_DsFit_dR_Kp_Ds_vec);
    tree->Branch("best_DsFit_dR_Km_Ds",  &best_DsFit_dR_Km_Ds_vec);
    tree->Branch("best_DsFit_dR_phi_Ds", &best_DsFit_dR_phi_Ds_vec);
    tree->Branch("best_DsFit_dR_pi_Ds",  &best_DsFit_dR_pi_Ds_vec);

    tree->Branch("best_dxy_phi_Ds", &best_dxy_phi_Ds_vec);
    tree->Branch("best_dz_phi_Ds",  &best_dz_phi_Ds_vec);

    tree->Branch("best_DsFit_Mconstraint_Ds_invm", &best_DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("best_phiFit_isValid",  &best_phiFit_vec.isValid);
    tree->Branch("best_phiFit_isFake",   &best_phiFit_vec.isFake);
    tree->Branch("best_phiFit_chi2",     &best_phiFit_vec.chi2);
    tree->Branch("best_phiFit_ndof",     &best_phiFit_vec.ndof);
    tree->Branch("best_phiFit_chi2ndof", &best_phiFit_vec.chi2ndof);
    tree->Branch("best_phiFit_vx",       &best_phiFit_vec.vx);
    tree->Branch("best_phiFit_vy",       &best_phiFit_vec.vy);
    tree->Branch("best_phiFit_vz",       &best_phiFit_vec.vz);
    tree->Branch("best_phiFit_vxerr",    &best_phiFit_vec.vxerr);
    tree->Branch("best_phiFit_vyerr",    &best_phiFit_vec.vyerr);
    tree->Branch("best_phiFit_vzerr",    &best_phiFit_vec.vzerr);

    tree->Branch("best_DsFit_isValid",  &best_DsFit_vec.isValid);
    tree->Branch("best_DsFit_isFake",   &best_DsFit_vec.isFake);
    tree->Branch("best_DsFit_chi2",     &best_DsFit_vec.chi2);
    tree->Branch("best_DsFit_ndof",     &best_DsFit_vec.ndof);
    tree->Branch("best_DsFit_chi2ndof", &best_DsFit_vec.chi2ndof);
    tree->Branch("best_DsFit_vx",       &best_DsFit_vec.vx);
    tree->Branch("best_DsFit_vy",       &best_DsFit_vec.vy);
    tree->Branch("best_DsFit_vz",       &best_DsFit_vec.vz);
    tree->Branch("best_DsFit_vxerr",    &best_DsFit_vec.vxerr);
    tree->Branch("best_DsFit_vyerr",    &best_DsFit_vec.vyerr);
    tree->Branch("best_DsFit_vzerr",    &best_DsFit_vec.vzerr);

    tree->Branch("best_PVnoDs_isValid",  &best_PVnoDs_vec.isValid);
    tree->Branch("best_PVnoDs_isFake",   &best_PVnoDs_vec.isFake);
    tree->Branch("best_PVnoDs_chi2",     &best_PVnoDs_vec.chi2);
    tree->Branch("best_PVnoDs_ndof",     &best_PVnoDs_vec.ndof);
    tree->Branch("best_PVnoDs_chi2ndof", &best_PVnoDs_vec.chi2ndof);
    tree->Branch("best_PVnoDs_vx",       &best_PVnoDs_vec.vx);
    tree->Branch("best_PVnoDs_vy",       &best_PVnoDs_vec.vy);
    tree->Branch("best_PVnoDs_vz",       &best_PVnoDs_vec.vz);
    tree->Branch("best_PVnoDs_vxerr",    &best_PVnoDs_vec.vxerr);
    tree->Branch("best_PVnoDs_vyerr",    &best_PVnoDs_vec.vyerr);
    tree->Branch("best_PVnoDs_vzerr",    &best_PVnoDs_vec.vzerr);

    tree->Branch("best_PVwithDs_isValid",  &best_PVwithDs_vec.isValid);
    tree->Branch("best_PVwithDs_isFake",   &best_PVwithDs_vec.isFake);
    tree->Branch("best_PVwithDs_chi2",     &best_PVwithDs_vec.chi2);
    tree->Branch("best_PVwithDs_ndof",     &best_PVwithDs_vec.ndof);
    tree->Branch("best_PVwithDs_chi2ndof", &best_PVwithDs_vec.chi2ndof);
    tree->Branch("best_PVwithDs_vx",       &best_PVwithDs_vec.vx);
    tree->Branch("best_PVwithDs_vy",       &best_PVwithDs_vec.vy);
    tree->Branch("best_PVwithDs_vz",       &best_PVwithDs_vec.vz);
    tree->Branch("best_PVwithDs_vxerr",    &best_PVwithDs_vec.vxerr);
    tree->Branch("best_PVwithDs_vyerr",    &best_PVwithDs_vec.vyerr);
    tree->Branch("best_PVwithDs_vzerr",    &best_PVwithDs_vec.vzerr);

    tree->Branch("best_Ds_primvtx_FDxy",       &best_Ds_primvtx_FD_vec.FDxy);
    tree->Branch("best_Ds_primvtx_FDz",        &best_Ds_primvtx_FD_vec.FDz);
    tree->Branch("best_Ds_primvtx_FD",         &best_Ds_primvtx_FD_vec.FD);
    tree->Branch("best_Ds_primvtx_FDxyerr",    &best_Ds_primvtx_FD_vec.FDxyerr);
    tree->Branch("best_Ds_primvtx_FDzerr",     &best_Ds_primvtx_FD_vec.FDzerr);
    tree->Branch("best_Ds_primvtx_FDerr",      &best_Ds_primvtx_FD_vec.FDerr);
    tree->Branch("best_Ds_primvtx_FDxychi2",   &best_Ds_primvtx_FD_vec.FDxychi2);
    tree->Branch("best_Ds_primvtx_FDzchi2",    &best_Ds_primvtx_FD_vec.FDzchi2);
    tree->Branch("best_Ds_primvtx_FDchi2",     &best_Ds_primvtx_FD_vec.FDchi2);
    tree->Branch("best_Ds_primvtx_dira",       &best_Ds_primvtx_FD_vec.dira);
    tree->Branch("best_Ds_primvtx_dira_angle", &best_Ds_primvtx_FD_vec.dira_angle);
    tree->Branch("best_Kp_primvtx_ip",         &best_Kp_primvtx_ip_vec.ip);
    tree->Branch("best_Kp_primvtx_iperr",      &best_Kp_primvtx_ip_vec.iperr);
    tree->Branch("best_Kp_primvtx_ipchi2",     &best_Kp_primvtx_ip_vec.ipchi2);
    tree->Branch("best_Km_primvtx_ip",         &best_Km_primvtx_ip_vec.ip);
    tree->Branch("best_Km_primvtx_iperr",      &best_Km_primvtx_ip_vec.iperr);
    tree->Branch("best_Km_primvtx_ipchi2",     &best_Km_primvtx_ip_vec.ipchi2);
    tree->Branch("best_pi_primvtx_ip",         &best_pi_primvtx_ip_vec.ip);
    tree->Branch("best_pi_primvtx_iperr",      &best_pi_primvtx_ip_vec.iperr);
    tree->Branch("best_pi_primvtx_ipchi2",     &best_pi_primvtx_ip_vec.ipchi2);
    tree->Branch("best_phi_primvtx_ip",        &best_phi_primvtx_ip_vec.ip);
    tree->Branch("best_phi_primvtx_iperr",     &best_phi_primvtx_ip_vec.iperr);
    tree->Branch("best_phi_primvtx_ipchi2",    &best_phi_primvtx_ip_vec.ipchi2);
    tree->Branch("best_Ds_primvtx_ip",         &best_Ds_primvtx_ip_vec.ip);
    tree->Branch("best_Ds_primvtx_iperr",      &best_Ds_primvtx_ip_vec.iperr);
    tree->Branch("best_Ds_primvtx_ipchi2",     &best_Ds_primvtx_ip_vec.ipchi2);

    tree->Branch("best_Ds_PVnoDs_FDxy",       &best_Ds_PVnoDs_FD_vec.FDxy);
    tree->Branch("best_Ds_PVnoDs_FDz",        &best_Ds_PVnoDs_FD_vec.FDz);
    tree->Branch("best_Ds_PVnoDs_FD",         &best_Ds_PVnoDs_FD_vec.FD);
    tree->Branch("best_Ds_PVnoDs_FDxyerr",    &best_Ds_PVnoDs_FD_vec.FDxyerr);
    tree->Branch("best_Ds_PVnoDs_FDzerr",     &best_Ds_PVnoDs_FD_vec.FDzerr);
    tree->Branch("best_Ds_PVnoDs_FDerr",      &best_Ds_PVnoDs_FD_vec.FDerr);
    tree->Branch("best_Ds_PVnoDs_FDxychi2",   &best_Ds_PVnoDs_FD_vec.FDxychi2);
    tree->Branch("best_Ds_PVnoDs_FDzchi2",    &best_Ds_PVnoDs_FD_vec.FDzchi2);
    tree->Branch("best_Ds_PVnoDs_FDchi2",     &best_Ds_PVnoDs_FD_vec.FDchi2);
    tree->Branch("best_Ds_PVnoDs_dira",       &best_Ds_PVnoDs_FD_vec.dira);
    tree->Branch("best_Ds_PVnoDs_dira_angle", &best_Ds_PVnoDs_FD_vec.dira_angle);
    tree->Branch("best_Kp_PVnoDs_ip",         &best_Kp_PVnoDs_ip_vec.ip);
    tree->Branch("best_Kp_PVnoDs_iperr",      &best_Kp_PVnoDs_ip_vec.iperr);
    tree->Branch("best_Kp_PVnoDs_ipchi2",     &best_Kp_PVnoDs_ip_vec.ipchi2);
    tree->Branch("best_Km_PVnoDs_ip",         &best_Km_PVnoDs_ip_vec.ip);
    tree->Branch("best_Km_PVnoDs_iperr",      &best_Km_PVnoDs_ip_vec.iperr);
    tree->Branch("best_Km_PVnoDs_ipchi2",     &best_Km_PVnoDs_ip_vec.ipchi2);
    tree->Branch("best_pi_PVnoDs_ip",         &best_pi_PVnoDs_ip_vec.ip);
    tree->Branch("best_pi_PVnoDs_iperr",      &best_pi_PVnoDs_ip_vec.iperr);
    tree->Branch("best_pi_PVnoDs_ipchi2",     &best_pi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("best_phi_PVnoDs_ip",        &best_phi_PVnoDs_ip_vec.ip);
    tree->Branch("best_phi_PVnoDs_iperr",     &best_phi_PVnoDs_ip_vec.iperr);
    tree->Branch("best_phi_PVnoDs_ipchi2",    &best_phi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("best_Ds_PVnoDs_ip",         &best_Ds_PVnoDs_ip_vec.ip);
    tree->Branch("best_Ds_PVnoDs_iperr",      &best_Ds_PVnoDs_ip_vec.iperr);
    tree->Branch("best_Ds_PVnoDs_ipchi2",     &best_Ds_PVnoDs_ip_vec.ipchi2);

    tree->Branch("best_Ds_PVwithDs_FDxy",       &best_Ds_PVwithDs_FD_vec.FDxy);
    tree->Branch("best_Ds_PVwithDs_FDz",        &best_Ds_PVwithDs_FD_vec.FDz);
    tree->Branch("best_Ds_PVwithDs_FD",         &best_Ds_PVwithDs_FD_vec.FD);
    tree->Branch("best_Ds_PVwithDs_FDxyerr",    &best_Ds_PVwithDs_FD_vec.FDxyerr);
    tree->Branch("best_Ds_PVwithDs_FDzerr",     &best_Ds_PVwithDs_FD_vec.FDzerr);
    tree->Branch("best_Ds_PVwithDs_FDerr",      &best_Ds_PVwithDs_FD_vec.FDerr);
    tree->Branch("best_Ds_PVwithDs_FDxychi2",   &best_Ds_PVwithDs_FD_vec.FDxychi2);
    tree->Branch("best_Ds_PVwithDs_FDzchi2",    &best_Ds_PVwithDs_FD_vec.FDzchi2);
    tree->Branch("best_Ds_PVwithDs_FDchi2",     &best_Ds_PVwithDs_FD_vec.FDchi2);
    tree->Branch("best_Ds_PVwithDs_dira",       &best_Ds_PVwithDs_FD_vec.dira);
    tree->Branch("best_Ds_PVwithDs_dira_angle", &best_Ds_PVwithDs_FD_vec.dira_angle);
    tree->Branch("best_Kp_PVwithDs_ip",         &best_Kp_PVwithDs_ip_vec.ip);
    tree->Branch("best_Kp_PVwithDs_iperr",      &best_Kp_PVwithDs_ip_vec.iperr);
    tree->Branch("best_Kp_PVwithDs_ipchi2",     &best_Kp_PVwithDs_ip_vec.ipchi2);
    tree->Branch("best_Km_PVwithDs_ip",         &best_Km_PVwithDs_ip_vec.ip);
    tree->Branch("best_Km_PVwithDs_iperr",      &best_Km_PVwithDs_ip_vec.iperr);
    tree->Branch("best_Km_PVwithDs_ipchi2",     &best_Km_PVwithDs_ip_vec.ipchi2);
    tree->Branch("best_pi_PVwithDs_ip",         &best_pi_PVwithDs_ip_vec.ip);
    tree->Branch("best_pi_PVwithDs_iperr",      &best_pi_PVwithDs_ip_vec.iperr);
    tree->Branch("best_pi_PVwithDs_ipchi2",     &best_pi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("best_phi_PVwithDs_ip",        &best_phi_PVwithDs_ip_vec.ip);
    tree->Branch("best_phi_PVwithDs_iperr",     &best_phi_PVwithDs_ip_vec.iperr);
    tree->Branch("best_phi_PVwithDs_ipchi2",    &best_phi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("best_Ds_PVwithDs_ip",         &best_Ds_PVwithDs_ip_vec.ip);
    tree->Branch("best_Ds_PVwithDs_iperr",      &best_Ds_PVwithDs_ip_vec.iperr);
    tree->Branch("best_Ds_PVwithDs_ipchi2",     &best_Ds_PVwithDs_ip_vec.ipchi2);

    tree->Branch("best_Ds_IsoR03_sumChargedHadronPt", &best_Ds_IsoR03_vec.sumChargedHadronPt);
    tree->Branch("best_Ds_IsoR03_sumNeutralHadronEt", &best_Ds_IsoR03_vec.sumNeutralHadronEt);
    tree->Branch("best_Ds_IsoR03_sumPhotonEt",        &best_Ds_IsoR03_vec.sumPhotonEt);
    tree->Branch("best_Ds_IsoR03_sumPUPt",            &best_Ds_IsoR03_vec.sumPUPt);
    tree->Branch("best_Ds_IsoR03_PFIso",              &best_Ds_IsoR03_vec.PFIso);

    tree->Branch("best_Ds_IsoR04_sumChargedHadronPt", &best_Ds_IsoR04_vec.sumChargedHadronPt);
    tree->Branch("best_Ds_IsoR04_sumNeutralHadronEt", &best_Ds_IsoR04_vec.sumNeutralHadronEt);
    tree->Branch("best_Ds_IsoR04_sumPhotonEt",        &best_Ds_IsoR04_vec.sumPhotonEt);
    tree->Branch("best_Ds_IsoR04_sumPUPt",            &best_Ds_IsoR04_vec.sumPUPt);
    tree->Branch("best_Ds_IsoR04_PFIso",              &best_Ds_IsoR04_vec.PFIso);

    tree->Branch("best_mu_charge",    &mu_vec.charge);
    tree->Branch("best_mu_eta",       &mu_vec.eta);
    tree->Branch("best_mu_phi",       &mu_vec.phi);
    tree->Branch("best_mu_vx",        &mu_vec.vx);
    tree->Branch("best_mu_vy",        &mu_vec.vy);
    tree->Branch("best_mu_vz",        &mu_vec.vz);
    tree->Branch("best_mu_p",         &mu_vec.p);
    tree->Branch("best_mu_pt",        &mu_vec.pt);
    tree->Branch("best_mu_px",        &mu_vec.px);
    tree->Branch("best_mu_py",        &mu_vec.py);
    tree->Branch("best_mu_pz",        &mu_vec.pz);
    tree->Branch("best_mu_dxy",       &mu_vec.dxy);
    tree->Branch("best_mu_dxyerr",    &mu_vec.dxyerr);
    tree->Branch("best_mu_dz",        &mu_vec.dz);
    tree->Branch("best_mu_dzerr",     &mu_vec.dzerr);
    tree->Branch("best_mu_isHighPt",  &mu_vec.isHighPt);
    tree->Branch("best_mu_isLoose",   &mu_vec.isLoose);
    tree->Branch("best_mu_isMedium",  &mu_vec.isMedium);
    tree->Branch("best_mu_isSoft",    &mu_vec.isSoft);
    tree->Branch("best_mu_isTight",   &mu_vec.isTight);
    tree->Branch("best_mu_isPF",      &mu_vec.isPF);
    tree->Branch("best_mu_isTracker", &mu_vec.isTracker);
    tree->Branch("best_mu_isGlobal",  &mu_vec.isGlobal);

    tree->Branch("best_mu_IsoR03_sumChargedHadronPt", &mu_IsoR03_vec.sumChargedHadronPt);
    tree->Branch("best_mu_IsoR03_sumNeutralHadronEt", &mu_IsoR03_vec.sumNeutralHadronEt);
    tree->Branch("best_mu_IsoR03_sumPhotonEt",        &mu_IsoR03_vec.sumPhotonEt);
    tree->Branch("best_mu_IsoR03_sumPUPt",            &mu_IsoR03_vec.sumPUPt);
    tree->Branch("best_mu_IsoR03_PFIso",              &mu_IsoR03_vec.PFIso);

    tree->Branch("best_mu_IsoR04_sumChargedHadronPt", &mu_IsoR04_vec.sumChargedHadronPt);
    tree->Branch("best_mu_IsoR04_sumNeutralHadronEt", &mu_IsoR04_vec.sumNeutralHadronEt);
    tree->Branch("best_mu_IsoR04_sumPhotonEt",        &mu_IsoR04_vec.sumPhotonEt);
    tree->Branch("best_mu_IsoR04_sumPUPt",            &mu_IsoR04_vec.sumPUPt);
    tree->Branch("best_mu_IsoR04_PFIso",              &mu_IsoR04_vec.PFIso);

    tree->Branch("best_mu_primvtx_ip",     &mu_primvtx_ip_vec.ip);
    tree->Branch("best_mu_primvtx_iperr",  &mu_primvtx_ip_vec.iperr);
    tree->Branch("best_mu_primvtx_ipchi2", &mu_primvtx_ip_vec.ipchi2);
}

void RecoTree::FillDs()
{
    Kp_vec.fillAll( Kp );
    Km_vec.fillAll( Km );
    pi_vec.fillAll( pi );
    phi_vec.fillAll( phi );
    Ds_vec.fillAll( Ds );
    phiFit_Kp_vec.fillAll( phiFit_Kp );
    phiFit_Km_vec.fillAll( phiFit_Km );
    phiFit_pi_vec.fillAll( phiFit_pi );
    phiFit_phi_vec.fillAll( phiFit_phi );
    phiFit_Ds_vec.fillAll( phiFit_Ds );
    DsFit_Kp_vec.fillAll( DsFit_Kp );
    DsFit_Km_vec.fillAll( DsFit_Km );
    DsFit_pi_vec.fillAll( DsFit_pi );
    DsFit_phi_vec.fillAll( DsFit_phi );
    DsFit_Ds_vec.fillAll( DsFit_Ds );

    dR_Kp_Km_vec.push_back(dR_Kp_Km);
    dR_Kp_phi_vec.push_back(dR_Kp_phi);
    dR_Km_phi_vec.push_back(dR_Km_phi);
    dR_Kp_pi_vec.push_back(dR_Kp_pi);
    dR_Km_pi_vec.push_back(dR_Km_pi);
    dR_pi_phi_vec.push_back(dR_pi_phi);
    dR_Kp_Ds_vec.push_back(dR_Kp_Ds);
    dR_Km_Ds_vec.push_back(dR_Km_Ds);
    dR_phi_Ds_vec.push_back(dR_phi_Ds);
    dR_pi_Ds_vec.push_back(dR_pi_Ds);

    dxy_phi_Ds_vec.push_back(dxy_phi_Ds);
    dz_phi_Ds_vec.push_back(dz_phi_Ds);

    phiFit_dR_Kp_Km_vec.push_back(phiFit_dR_Kp_Km);
    phiFit_dR_Kp_phi_vec.push_back(phiFit_dR_Kp_phi);
    phiFit_dR_Km_phi_vec.push_back(phiFit_dR_Km_phi);
    phiFit_dR_Kp_pi_vec.push_back(phiFit_dR_Kp_pi);
    phiFit_dR_Km_pi_vec.push_back(phiFit_dR_Km_pi);
    phiFit_dR_pi_phi_vec.push_back(phiFit_dR_pi_phi);
    phiFit_dR_Kp_Ds_vec.push_back(phiFit_dR_Kp_Ds);
    phiFit_dR_Km_Ds_vec.push_back(phiFit_dR_Km_Ds);
    phiFit_dR_phi_Ds_vec.push_back(phiFit_dR_phi_Ds);
    phiFit_dR_pi_Ds_vec.push_back(phiFit_dR_pi_Ds);

    DsFit_dR_Kp_Km_vec.push_back(DsFit_dR_Kp_Km);
    DsFit_dR_Kp_phi_vec.push_back(DsFit_dR_Kp_phi);
    DsFit_dR_Km_phi_vec.push_back(DsFit_dR_Km_phi);
    DsFit_dR_Kp_pi_vec.push_back(DsFit_dR_Kp_pi);
    DsFit_dR_Km_pi_vec.push_back(DsFit_dR_Km_pi);
    DsFit_dR_pi_phi_vec.push_back(DsFit_dR_pi_phi);
    DsFit_dR_Kp_Ds_vec.push_back(DsFit_dR_Kp_Ds);
    DsFit_dR_Km_Ds_vec.push_back(DsFit_dR_Km_Ds);
    DsFit_dR_phi_Ds_vec.push_back(DsFit_dR_phi_Ds);
    DsFit_dR_pi_Ds_vec.push_back(DsFit_dR_pi_Ds);

    DsFit_Mconstraint_Ds_invm_vec.push_back(DsFit_Mconstraint_Ds_invm);

    phiFit_vec.fillAll( phiFit );
    DsFit_vec.fillAll( DsFit );
    PVnoDs_vec.fillAll( PVnoDs );
    PVwithDs_vec.fillAll( PVwithDs );

    Ds_primvtx_FD_vec.fillAll( Ds_primvtx_FD );
    Kp_primvtx_ip_vec.fillAll( Kp_primvtx_ip );
    Km_primvtx_ip_vec.fillAll( Km_primvtx_ip );
    pi_primvtx_ip_vec.fillAll( pi_primvtx_ip );
    phi_primvtx_ip_vec.fillAll( phi_primvtx_ip );
    Ds_primvtx_ip_vec.fillAll( Ds_primvtx_ip );

    Ds_PVnoDs_FD_vec.fillAll( Ds_PVnoDs_FD );
    Kp_PVnoDs_ip_vec.fillAll( Kp_PVnoDs_ip );
    Km_PVnoDs_ip_vec.fillAll( Km_PVnoDs_ip );
    pi_PVnoDs_ip_vec.fillAll( pi_PVnoDs_ip );
    phi_PVnoDs_ip_vec.fillAll( phi_PVnoDs_ip );
    Ds_PVnoDs_ip_vec.fillAll( Ds_PVnoDs_ip );

    Ds_PVwithDs_FD_vec.fillAll( Ds_PVwithDs_FD );
    Kp_PVwithDs_ip_vec.fillAll( Kp_PVwithDs_ip );
    Km_PVwithDs_ip_vec.fillAll( Km_PVwithDs_ip );
    pi_PVwithDs_ip_vec.fillAll( pi_PVwithDs_ip );
    phi_PVwithDs_ip_vec.fillAll( phi_PVwithDs_ip );
    Ds_PVwithDs_ip_vec.fillAll( Ds_PVwithDs_ip );

    Ds_IsoR03_vec.fillAll( Ds_IsoR03 );
    Ds_IsoR04_vec.fillAll( Ds_IsoR04 );
}

void RecoTree::FillBestDs(int idxmax)
{
    best_Kp_vec.fillAll( Kp_vec.findIdx(idxmax) );
    best_Km_vec.fillAll( Km_vec.findIdx(idxmax) );
    best_pi_vec.fillAll( pi_vec.findIdx(idxmax) );
    best_phi_vec.fillAll( phi_vec.findIdx(idxmax) );
    best_Ds_vec.fillAll( Ds_vec.findIdx(idxmax) );
    best_phiFit_Kp_vec.fillAll( phiFit_Kp_vec.findIdx(idxmax) );
    best_phiFit_Km_vec.fillAll( phiFit_Km_vec.findIdx(idxmax) );
    best_phiFit_pi_vec.fillAll( phiFit_pi_vec.findIdx(idxmax) );
    best_phiFit_phi_vec.fillAll( phiFit_phi_vec.findIdx(idxmax) );
    best_phiFit_Ds_vec.fillAll( phiFit_Ds_vec.findIdx(idxmax) );
    best_DsFit_Kp_vec.fillAll( DsFit_Kp_vec.findIdx(idxmax) );
    best_DsFit_Km_vec.fillAll( DsFit_Km_vec.findIdx(idxmax) );
    best_DsFit_pi_vec.fillAll( DsFit_pi_vec.findIdx(idxmax) );
    best_DsFit_phi_vec.fillAll( DsFit_phi_vec.findIdx(idxmax) );
    best_DsFit_Ds_vec.fillAll( DsFit_Ds_vec.findIdx(idxmax) );

    best_dR_Kp_Km_vec.push_back( dR_Kp_Km_vec[idxmax] );
    best_dR_Kp_phi_vec.push_back( dR_Kp_phi_vec[idxmax] );
    best_dR_Km_phi_vec.push_back( dR_Km_phi_vec[idxmax] );
    best_dR_Kp_pi_vec.push_back( dR_Kp_pi_vec[idxmax] );
    best_dR_Km_pi_vec.push_back( dR_Km_pi_vec[idxmax] );
    best_dR_pi_phi_vec.push_back( dR_pi_phi_vec[idxmax] );
    best_dR_Kp_Ds_vec.push_back( dR_Kp_Ds_vec[idxmax] );
    best_dR_Km_Ds_vec.push_back( dR_Km_Ds_vec[idxmax] );
    best_dR_phi_Ds_vec.push_back( dR_phi_Ds_vec[idxmax] );
    best_dR_pi_Ds_vec.push_back( dR_pi_Ds_vec[idxmax] );

    best_phiFit_dR_Kp_Km_vec.push_back( phiFit_dR_Kp_Km_vec[idxmax] );
    best_phiFit_dR_Kp_phi_vec.push_back( phiFit_dR_Kp_phi_vec[idxmax] );
    best_phiFit_dR_Km_phi_vec.push_back( phiFit_dR_Km_phi_vec[idxmax] );
    best_phiFit_dR_Kp_pi_vec.push_back( phiFit_dR_Kp_pi_vec[idxmax] );
    best_phiFit_dR_Km_pi_vec.push_back( phiFit_dR_Km_pi_vec[idxmax] );
    best_phiFit_dR_pi_phi_vec.push_back( phiFit_dR_pi_phi_vec[idxmax] );
    best_phiFit_dR_Kp_Ds_vec.push_back( phiFit_dR_Kp_Ds_vec[idxmax] );
    best_phiFit_dR_Km_Ds_vec.push_back( phiFit_dR_Km_Ds_vec[idxmax] );
    best_phiFit_dR_phi_Ds_vec.push_back( phiFit_dR_phi_Ds_vec[idxmax] );
    best_phiFit_dR_pi_Ds_vec.push_back( phiFit_dR_pi_Ds_vec[idxmax] );

    best_DsFit_dR_Kp_Km_vec.push_back( DsFit_dR_Kp_Km_vec[idxmax] );
    best_DsFit_dR_Kp_phi_vec.push_back( DsFit_dR_Kp_phi_vec[idxmax] );
    best_DsFit_dR_Km_phi_vec.push_back( DsFit_dR_Km_phi_vec[idxmax] );
    best_DsFit_dR_Kp_pi_vec.push_back( DsFit_dR_Kp_pi_vec[idxmax] );
    best_DsFit_dR_Km_pi_vec.push_back( DsFit_dR_Km_pi_vec[idxmax] );
    best_DsFit_dR_pi_phi_vec.push_back( DsFit_dR_pi_phi_vec[idxmax] );
    best_DsFit_dR_Kp_Ds_vec.push_back( DsFit_dR_Kp_Ds_vec[idxmax] );
    best_DsFit_dR_Km_Ds_vec.push_back( DsFit_dR_Km_Ds_vec[idxmax] );
    best_DsFit_dR_phi_Ds_vec.push_back( DsFit_dR_phi_Ds_vec[idxmax] );
    best_DsFit_dR_pi_Ds_vec.push_back( DsFit_dR_pi_Ds_vec[idxmax] );

    best_dxy_phi_Ds_vec.push_back( dxy_phi_Ds_vec[idxmax] );
    best_dz_phi_Ds_vec.push_back( dz_phi_Ds_vec[idxmax] );

    best_DsFit_Mconstraint_Ds_invm_vec.push_back( DsFit_Mconstraint_Ds_invm_vec[idxmax] );

    best_phiFit_vec.fillAll( phiFit_vec.findIdx(idxmax) );
    best_DsFit_vec.fillAll( DsFit_vec.findIdx(idxmax) );
    best_PVnoDs_vec.fillAll( PVnoDs_vec.findIdx(idxmax) );
    best_PVwithDs_vec.fillAll( PVwithDs_vec.findIdx(idxmax) );

    best_Ds_primvtx_FD_vec.fillAll( Ds_primvtx_FD_vec.findIdx(idxmax) );
    best_Kp_primvtx_ip_vec.fillAll( Kp_primvtx_ip_vec.findIdx(idxmax) );
    best_Km_primvtx_ip_vec.fillAll( Km_primvtx_ip_vec.findIdx(idxmax) );
    best_pi_primvtx_ip_vec.fillAll( pi_primvtx_ip_vec.findIdx(idxmax) );
    best_phi_primvtx_ip_vec.fillAll( phi_primvtx_ip_vec.findIdx(idxmax) );
    best_Ds_primvtx_ip_vec.fillAll( Ds_primvtx_ip_vec.findIdx(idxmax) );

    best_Ds_PVnoDs_FD_vec.fillAll( Ds_PVnoDs_FD_vec.findIdx(idxmax) );
    best_Kp_PVnoDs_ip_vec.fillAll( Kp_PVnoDs_ip_vec.findIdx(idxmax) );
    best_Km_PVnoDs_ip_vec.fillAll( Km_PVnoDs_ip_vec.findIdx(idxmax) );
    best_pi_PVnoDs_ip_vec.fillAll( pi_PVnoDs_ip_vec.findIdx(idxmax) );
    best_phi_PVnoDs_ip_vec.fillAll( phi_PVnoDs_ip_vec.findIdx(idxmax) );
    best_Ds_PVnoDs_ip_vec.fillAll( Ds_PVnoDs_ip_vec.findIdx(idxmax) );

    best_Ds_PVwithDs_FD_vec.fillAll( Ds_PVwithDs_FD_vec.findIdx(idxmax) );
    best_Kp_PVwithDs_ip_vec.fillAll( Kp_PVwithDs_ip_vec.findIdx(idxmax) );
    best_Km_PVwithDs_ip_vec.fillAll( Km_PVwithDs_ip_vec.findIdx(idxmax) );
    best_pi_PVwithDs_ip_vec.fillAll( pi_PVwithDs_ip_vec.findIdx(idxmax) );
    best_phi_PVwithDs_ip_vec.fillAll( phi_PVwithDs_ip_vec.findIdx(idxmax) );
    best_Ds_PVwithDs_ip_vec.fillAll( Ds_PVwithDs_ip_vec.findIdx(idxmax) );

    best_Ds_IsoR03_vec.fillAll( Ds_IsoR03_vec.findIdx(idxmax) );
    best_Ds_IsoR04_vec.fillAll( Ds_IsoR04_vec.findIdx(idxmax) );
}

void RecoTree::FillBestmu()
{
    mu_vec.fillAll( mu );
    mu_IsoR03_vec.fillAll( mu_IsoR03 );
    mu_IsoR04_vec.fillAll( mu_IsoR04 );
    mu_primvtx_ip_vec.fillAll( mu_primvtx_ip );
}

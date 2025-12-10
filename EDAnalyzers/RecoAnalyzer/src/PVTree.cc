#include "EDAnalyzers/RecoAnalyzer/interface/PVTree.h"
#include <iostream>

PVTree::PVTree(TTree *tree_)
{
    tree = tree_;
}

void PVTree::GenInfo::fillAll( const reco::GenParticle& pruned )
{
    eta = pruned.eta();
    phi = pruned.phi();
    vx  = pruned.vx();
    vy  = pruned.vy();
    vz  = pruned.vz();
    p   = pruned.p();
    pt  = pruned.pt();
    px  = pruned.px();
    py  = pruned.py();
    pz  = pruned.pz();
}

void PVTree::BSInfo::fillAll( const reco::BeamSpot& beamspot )
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

void PVTree::PFInfo::fillAll( const pat::PackedCandidate& pfCands )
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

void PVTree::TLVInfo::fillAll( const TLorentzVector& TLV )
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

void PVTree::VtxInfo::fillAll( const reco::Vertex& Vtx )
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

void PVTree::FDInfo::fillAll( const GlobalVector& Pos, const GlobalError& Cov, const GlobalVector& PDirec )
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

void PVTree::IPInfo::fillAll( const reco::TransientTrack& ttrack, const reco::Vertex& vtx )
{
    std::pair<bool, Measurement1D> result = IPTools::absoluteImpactParameter3D(ttrack, vtx);
    if( result.first ){
        ip     = result.second.value();
        iperr  = result.second.error();
        ipchi2 = result.second.significance();
    }
}

void PVTree::MuonInfo::fillAll( const pat::Muon& muon, const reco::Vertex & Vtx )
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


void PVTree::Init()
{
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    num_Gen_phi = 0;
    num_Gen_Ds = 0;
    num_Gen_mu = 0;
    num_Gen_nu = 0;
    num_Gen_W = 0;
    num_Gen_H = 0;

    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;
    num_match_mu = 0;
    num_tight_match_Kp = 0;
    num_tight_match_Km = 0;
    num_tight_match_pi = 0;

    match_dR_Kp = null;
    match_dR_Km = null;
    match_dR_pi = null;
    match_dR_mu = null;

    match_Kp_vec.clearAll();
    match_Km_vec.clearAll();
    match_pi_vec.clearAll();
    match_phi_vec.clearAll();
    match_Ds_vec.clearAll();
    match_phiFit_Kp_vec.clearAll();
    match_phiFit_Km_vec.clearAll();
    match_phiFit_pi_vec.clearAll();
    match_phiFit_phi_vec.clearAll();
    match_phiFit_Ds_vec.clearAll();
    match_DsFit_Kp_vec.clearAll();
    match_DsFit_Km_vec.clearAll();
    match_DsFit_pi_vec.clearAll();
    match_DsFit_phi_vec.clearAll();
    match_DsFit_Ds_vec.clearAll();

    match_dR_Kp_Km_vec.clear();
    match_dR_Kp_phi_vec.clear();
    match_dR_Km_phi_vec.clear();
    match_dR_Kp_pi_vec.clear();
    match_dR_Km_pi_vec.clear();
    match_dR_pi_phi_vec.clear();
    match_dR_Kp_Ds_vec.clear();
    match_dR_Km_Ds_vec.clear();
    match_dR_phi_Ds_vec.clear();
    match_dR_pi_Ds_vec.clear();

    match_phiFit_dR_Kp_Km_vec.clear();
    match_phiFit_dR_Kp_phi_vec.clear();
    match_phiFit_dR_Km_phi_vec.clear();
    match_phiFit_dR_Kp_pi_vec.clear();
    match_phiFit_dR_Km_pi_vec.clear();
    match_phiFit_dR_pi_phi_vec.clear();
    match_phiFit_dR_Kp_Ds_vec.clear();
    match_phiFit_dR_Km_Ds_vec.clear();
    match_phiFit_dR_phi_Ds_vec.clear();
    match_phiFit_dR_pi_Ds_vec.clear();

    match_DsFit_dR_Kp_Km_vec.clear();
    match_DsFit_dR_Kp_phi_vec.clear();
    match_DsFit_dR_Km_phi_vec.clear();
    match_DsFit_dR_Kp_pi_vec.clear();
    match_DsFit_dR_Km_pi_vec.clear();
    match_DsFit_dR_pi_phi_vec.clear();
    match_DsFit_dR_Kp_Ds_vec.clear();
    match_DsFit_dR_Km_Ds_vec.clear();
    match_DsFit_dR_phi_Ds_vec.clear();
    match_DsFit_dR_pi_Ds_vec.clear();

    match_dxy_phi_Ds_vec.clear();
    match_dz_phi_Ds_vec.clear();

    match_DsFit_Mconstraint_Ds_invm_vec.clear();

    match_phiFit_vec.clearAll();
    match_DsFit_vec.clearAll();
    match_PVnoDs_vec.clearAll();
    match_PVwithDs_vec.clearAll();

    match_Ds_primvtx_FD_vec.clearAll();
    match_Kp_primvtx_ip_vec.clearAll();
    match_Km_primvtx_ip_vec.clearAll();
    match_pi_primvtx_ip_vec.clearAll();
    match_phi_primvtx_ip_vec.clearAll();
    match_Ds_primvtx_ip_vec.clearAll();

    match_Ds_PVnoDs_FD_vec.clearAll();
    match_Kp_PVnoDs_ip_vec.clearAll();
    match_Km_PVnoDs_ip_vec.clearAll();
    match_pi_PVnoDs_ip_vec.clearAll();
    match_phi_PVnoDs_ip_vec.clearAll();
    match_Ds_PVnoDs_ip_vec.clearAll();

    match_Ds_PVwithDs_FD_vec.clearAll();
    match_Kp_PVwithDs_ip_vec.clearAll();
    match_Km_PVwithDs_ip_vec.clearAll();
    match_pi_PVwithDs_ip_vec.clearAll();
    match_phi_PVwithDs_ip_vec.clearAll();
    match_Ds_PVwithDs_ip_vec.clearAll();

    match_Ds_IsoR03_vec.clearAll();
    match_Ds_IsoR04_vec.clearAll();
 
    match_mu_vec.clearAll();
    match_mu_IsoR03_vec.clearAll();
    match_mu_IsoR04_vec.clearAll();
    match_mu_primvtx_ip_vec.clearAll();
    
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
 
    Kp_match_vec.clear();
    Km_match_vec.clear();
    pi_match_vec.clear();
    match_entry_vec.clear();
    non_match_entry_vec.clear();

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
 
    best_match_entry_vec.clear();
 
    mu_vec.clearAll();
    mu_IsoR03_vec.clearAll();
    mu_IsoR04_vec.clearAll();
    mu_primvtx_ip_vec.clearAll();
    best_mu_match_vec.clear();
}

void PVTree::CreateBranches()
{
    tree->Branch("num_Gen_Kp",  &num_Gen_Kp,  "num_Gen_Kp/I");
    tree->Branch("num_Gen_Km",  &num_Gen_Km,  "num_Gen_Km/I");
    tree->Branch("num_Gen_pi",  &num_Gen_pi,  "num_Gen_pi/I");
    tree->Branch("num_Gen_phi", &num_Gen_phi, "num_Gen_phi/I");
    tree->Branch("num_Gen_Ds",  &num_Gen_Ds,  "num_Gen_Ds/I");
    tree->Branch("num_Gen_mu",  &num_Gen_mu,  "num_Gen_mu/I");
    tree->Branch("num_Gen_nu",  &num_Gen_nu,  "num_Gen_nu/I");
    tree->Branch("num_Gen_W",   &num_Gen_W,   "num_Gen_W/I");
    tree->Branch("num_Gen_H",   &num_Gen_H,   "num_Gen_H/I");

    tree->Branch("Gen_H",   &Gen_H,   "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_Ds",  &Gen_Ds,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_W",   &Gen_W,   "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_phi", &Gen_phi, "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_Kp",  &Gen_Kp,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_Km",  &Gen_Km,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_pi",  &Gen_pi,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_mu",  &Gen_mu,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");
    tree->Branch("Gen_nu",  &Gen_nu,  "eta/D:phi:vx:vy:vz:p:pt:px:py:pz");

    tree->Branch("Gen_dR_Kp_Km",  &Gen_dR_Kp_Km,  "Gen_dR_Kp_Km/D");
    tree->Branch("Gen_dR_Kp_phi", &Gen_dR_Kp_phi, "Gen_dR_Kp_phi/D");
    tree->Branch("Gen_dR_Km_phi", &Gen_dR_Km_phi, "Gen_dR_Km_phi/D");
    tree->Branch("Gen_dR_Kp_pi",  &Gen_dR_Kp_pi,  "Gen_dR_Kp_pi/D");
    tree->Branch("Gen_dR_Km_pi",  &Gen_dR_Km_pi,  "Gen_dR_Km_pi/D");
    tree->Branch("Gen_dR_pi_phi", &Gen_dR_pi_phi, "Gen_dR_pi_phi/D");
    tree->Branch("Gen_dR_Kp_Ds",  &Gen_dR_Kp_Ds,  "Gen_dR_Kp_Ds/D");
    tree->Branch("Gen_dR_Km_Ds",  &Gen_dR_Km_Ds,  "Gen_dR_Km_Ds/D");
    tree->Branch("Gen_dR_phi_Ds", &Gen_dR_phi_Ds, "Gen_dR_phi_Ds/D");
    tree->Branch("Gen_dR_pi_Ds",  &Gen_dR_pi_Ds,  "Gen_dR_pi_Ds/D");
    tree->Branch("Gen_dR_Kp_mu",  &Gen_dR_Kp_mu,  "Gen_dR_Kp_mu/D");
    tree->Branch("Gen_dR_Km_mu",  &Gen_dR_Km_mu,  "Gen_dR_Km_mu/D");
    tree->Branch("Gen_dR_phi_mu", &Gen_dR_phi_mu, "Gen_dR_phi_mu/D");
    tree->Branch("Gen_dR_pi_mu",  &Gen_dR_pi_mu,  "Gen_dR_pi_mu/D");
    tree->Branch("Gen_dR_Ds_mu",  &Gen_dR_Ds_mu,  "Gen_dR_Ds_mu/D");

    tree->Branch("Gen_Ds_dx",   &Gen_Ds_dx,   "Gen_Ds_dx/D");
    tree->Branch("Gen_Ds_dy",   &Gen_Ds_dy,   "Gen_Ds_dy/D");
    tree->Branch("Gen_Ds_dz",   &Gen_Ds_dz,   "Gen_Ds_dz/D");
    tree->Branch("Gen_Ds_FDxy", &Gen_Ds_FDxy, "Gen_Ds_FDxy/D");
    tree->Branch("Gen_Ds_FD",   &Gen_Ds_FD,   "Gen_Ds_FD/D");

    tree->Branch("beamspot", &beamspotInfo, "type/I:x0/D:y0:z0:sigmaZ:dxdz:dydz:BWX:BWY:x0err:y0err:z0err:sigmaZ0err:dxdzerr:dydzerr:BWXerr:BWYerr:emitX:emitY:betaSta");
    tree->Branch("primvtx",  &primvtxInfo,  "isValid/O:isFake:chi2/D:ndof:chi2ndof:vx:vy:vz:vxerr:vyerr:vzerr"); 
    
    tree->Branch("num_match_Kp",       &num_match_Kp,       "num_match_Kp/I");
    tree->Branch("num_match_Km",       &num_match_Km,       "num_match_Km/I");
    tree->Branch("num_match_pi",       &num_match_pi,       "num_match_pi/I");
    tree->Branch("num_match_mu",       &num_match_mu,       "num_match_mu/I");
    tree->Branch("num_tight_match_Kp", &num_tight_match_Kp, "num_tight_match_Kp/I");
    tree->Branch("num_tight_match_Km", &num_tight_match_Km, "num_tight_match_Km/I");
    tree->Branch("num_tight_match_pi", &num_tight_match_pi, "num_tight_match_pi/I");
    tree->Branch("match_dR_Kp",        &match_dR_Kp,        "match_dR_Kp/D");
    tree->Branch("match_dR_Km",        &match_dR_Km,        "match_dR_Km/D");
    tree->Branch("match_dR_pi",        &match_dR_pi,        "match_dR_pi/D");
    tree->Branch("match_dR_mu",        &match_dR_mu,        "match_dR_mu/D");

    tree->Branch("match_Kp_isIsolatedChargedHadron", &match_Kp_vec.isIsolatedChargedHadron);
    tree->Branch("match_Kp_charge",                  &match_Kp_vec.charge);
    tree->Branch("match_Kp_eta",                     &match_Kp_vec.eta);
    tree->Branch("match_Kp_phi",                     &match_Kp_vec.phi);
    tree->Branch("match_Kp_vx",                      &match_Kp_vec.vx);
    tree->Branch("match_Kp_vy",                      &match_Kp_vec.vy);
    tree->Branch("match_Kp_vz",                      &match_Kp_vec.vz);
    tree->Branch("match_Kp_p",                       &match_Kp_vec.p);
    tree->Branch("match_Kp_pt",                      &match_Kp_vec.pt);
    tree->Branch("match_Kp_px",                      &match_Kp_vec.px);
    tree->Branch("match_Kp_py",                      &match_Kp_vec.py);
    tree->Branch("match_Kp_pz",                      &match_Kp_vec.pz);

    tree->Branch("match_Km_isIsolatedChargedHadron", &match_Km_vec.isIsolatedChargedHadron);
    tree->Branch("match_Km_charge",                  &match_Km_vec.charge);
    tree->Branch("match_Km_eta",                     &match_Km_vec.eta);
    tree->Branch("match_Km_phi",                     &match_Km_vec.phi);
    tree->Branch("match_Km_vx",                      &match_Km_vec.vx);
    tree->Branch("match_Km_vy",                      &match_Km_vec.vy);
    tree->Branch("match_Km_vz",                      &match_Km_vec.vz);
    tree->Branch("match_Km_p",                       &match_Km_vec.p);
    tree->Branch("match_Km_pt",                      &match_Km_vec.pt);
    tree->Branch("match_Km_px",                      &match_Km_vec.px);
    tree->Branch("match_Km_py",                      &match_Km_vec.py);
    tree->Branch("match_Km_pz",                      &match_Km_vec.pz);

    tree->Branch("match_pi_isIsolatedChargedHadron", &match_pi_vec.isIsolatedChargedHadron);
    tree->Branch("match_pi_charge",                  &match_pi_vec.charge);
    tree->Branch("match_pi_eta",                     &match_pi_vec.eta);
    tree->Branch("match_pi_phi",                     &match_pi_vec.phi);
    tree->Branch("match_pi_vx",                      &match_pi_vec.vx);
    tree->Branch("match_pi_vy",                      &match_pi_vec.vy);
    tree->Branch("match_pi_vz",                      &match_pi_vec.vz);
    tree->Branch("match_pi_p",                       &match_pi_vec.p);
    tree->Branch("match_pi_pt",                      &match_pi_vec.pt);
    tree->Branch("match_pi_px",                      &match_pi_vec.px);
    tree->Branch("match_pi_py",                      &match_pi_vec.py);
    tree->Branch("match_pi_pz",                      &match_pi_vec.pz);

    tree->Branch("match_phi_eta",  &match_phi_vec.eta);
    tree->Branch("match_phi_phi",  &match_phi_vec.phi);
    tree->Branch("match_phi_p",    &match_phi_vec.p);
    tree->Branch("match_phi_pt",   &match_phi_vec.pt);
    tree->Branch("match_phi_px",   &match_phi_vec.px);
    tree->Branch("match_phi_py",   &match_phi_vec.py);
    tree->Branch("match_phi_pz",   &match_phi_vec.pz);
    tree->Branch("match_phi_invm", &match_phi_vec.invm);

    tree->Branch("match_Ds_eta",  &match_Ds_vec.eta);
    tree->Branch("match_Ds_phi",  &match_Ds_vec.phi);
    tree->Branch("match_Ds_p",    &match_Ds_vec.p);
    tree->Branch("match_Ds_pt",   &match_Ds_vec.pt);
    tree->Branch("match_Ds_px",   &match_Ds_vec.px);
    tree->Branch("match_Ds_py",   &match_Ds_vec.py);
    tree->Branch("match_Ds_pz",   &match_Ds_vec.pz);
    tree->Branch("match_Ds_invm", &match_Ds_vec.invm);

    tree->Branch("match_phiFit_Kp_eta",  &match_phiFit_Kp_vec.eta);
    tree->Branch("match_phiFit_Kp_phi",  &match_phiFit_Kp_vec.phi);
    tree->Branch("match_phiFit_Kp_p",    &match_phiFit_Kp_vec.p);
    tree->Branch("match_phiFit_Kp_pt",   &match_phiFit_Kp_vec.pt);
    tree->Branch("match_phiFit_Kp_px",   &match_phiFit_Kp_vec.px);
    tree->Branch("match_phiFit_Kp_py",   &match_phiFit_Kp_vec.py);
    tree->Branch("match_phiFit_Kp_pz",   &match_phiFit_Kp_vec.pz);
    tree->Branch("match_phiFit_Kp_invm", &match_phiFit_Kp_vec.invm);

    tree->Branch("match_phiFit_Km_eta",  &match_phiFit_Km_vec.eta);
    tree->Branch("match_phiFit_Km_phi",  &match_phiFit_Km_vec.phi);
    tree->Branch("match_phiFit_Km_p",    &match_phiFit_Km_vec.p);
    tree->Branch("match_phiFit_Km_pt",   &match_phiFit_Km_vec.pt);
    tree->Branch("match_phiFit_Km_px",   &match_phiFit_Km_vec.px);
    tree->Branch("match_phiFit_Km_py",   &match_phiFit_Km_vec.py);
    tree->Branch("match_phiFit_Km_pz",   &match_phiFit_Km_vec.pz);
    tree->Branch("match_phiFit_Km_invm", &match_phiFit_Km_vec.invm);

    tree->Branch("match_phiFit_pi_eta",  &match_phiFit_pi_vec.eta);
    tree->Branch("match_phiFit_pi_phi",  &match_phiFit_pi_vec.phi);
    tree->Branch("match_phiFit_pi_p",    &match_phiFit_pi_vec.p);
    tree->Branch("match_phiFit_pi_pt",   &match_phiFit_pi_vec.pt);
    tree->Branch("match_phiFit_pi_px",   &match_phiFit_pi_vec.px);
    tree->Branch("match_phiFit_pi_py",   &match_phiFit_pi_vec.py);
    tree->Branch("match_phiFit_pi_pz",   &match_phiFit_pi_vec.pz);
    tree->Branch("match_phiFit_pi_invm", &match_phiFit_pi_vec.invm);

    tree->Branch("match_phiFit_phi_eta",  &match_phiFit_phi_vec.eta);
    tree->Branch("match_phiFit_phi_phi",  &match_phiFit_phi_vec.phi);
    tree->Branch("match_phiFit_phi_p",    &match_phiFit_phi_vec.p);
    tree->Branch("match_phiFit_phi_pt",   &match_phiFit_phi_vec.pt);
    tree->Branch("match_phiFit_phi_px",   &match_phiFit_phi_vec.px);
    tree->Branch("match_phiFit_phi_py",   &match_phiFit_phi_vec.py);
    tree->Branch("match_phiFit_phi_pz",   &match_phiFit_phi_vec.pz);
    tree->Branch("match_phiFit_phi_invm", &match_phiFit_phi_vec.invm);

    tree->Branch("match_phiFit_Ds_eta",  &match_phiFit_Ds_vec.eta);
    tree->Branch("match_phiFit_Ds_phi",  &match_phiFit_Ds_vec.phi);
    tree->Branch("match_phiFit_Ds_p",    &match_phiFit_Ds_vec.p);
    tree->Branch("match_phiFit_Ds_pt",   &match_phiFit_Ds_vec.pt);
    tree->Branch("match_phiFit_Ds_px",   &match_phiFit_Ds_vec.px);
    tree->Branch("match_phiFit_Ds_py",   &match_phiFit_Ds_vec.py);
    tree->Branch("match_phiFit_Ds_pz",   &match_phiFit_Ds_vec.pz);
    tree->Branch("match_phiFit_Ds_invm", &match_phiFit_Ds_vec.invm);

    tree->Branch("match_DsFit_Kp_eta",  &match_DsFit_Kp_vec.eta);
    tree->Branch("match_DsFit_Kp_phi",  &match_DsFit_Kp_vec.phi);
    tree->Branch("match_DsFit_Kp_p",    &match_DsFit_Kp_vec.p);
    tree->Branch("match_DsFit_Kp_pt",   &match_DsFit_Kp_vec.pt);
    tree->Branch("match_DsFit_Kp_px",   &match_DsFit_Kp_vec.px);
    tree->Branch("match_DsFit_Kp_py",   &match_DsFit_Kp_vec.py);
    tree->Branch("match_DsFit_Kp_pz",   &match_DsFit_Kp_vec.pz);
    tree->Branch("match_DsFit_Kp_invm", &match_DsFit_Kp_vec.invm);

    tree->Branch("match_DsFit_Km_eta",  &match_DsFit_Km_vec.eta);
    tree->Branch("match_DsFit_Km_phi",  &match_DsFit_Km_vec.phi);
    tree->Branch("match_DsFit_Km_p",    &match_DsFit_Km_vec.p);
    tree->Branch("match_DsFit_Km_pt",   &match_DsFit_Km_vec.pt);
    tree->Branch("match_DsFit_Km_px",   &match_DsFit_Km_vec.px);
    tree->Branch("match_DsFit_Km_py",   &match_DsFit_Km_vec.py);
    tree->Branch("match_DsFit_Km_pz",   &match_DsFit_Km_vec.pz);
    tree->Branch("match_DsFit_Km_invm", &match_DsFit_Km_vec.invm);

    tree->Branch("match_DsFit_pi_eta",  &match_DsFit_pi_vec.eta);
    tree->Branch("match_DsFit_pi_phi",  &match_DsFit_pi_vec.phi);
    tree->Branch("match_DsFit_pi_p",    &match_DsFit_pi_vec.p);
    tree->Branch("match_DsFit_pi_pt",   &match_DsFit_pi_vec.pt);
    tree->Branch("match_DsFit_pi_px",   &match_DsFit_pi_vec.px);
    tree->Branch("match_DsFit_pi_py",   &match_DsFit_pi_vec.py);
    tree->Branch("match_DsFit_pi_pz",   &match_DsFit_pi_vec.pz);
    tree->Branch("match_DsFit_pi_invm", &match_DsFit_pi_vec.invm);

    tree->Branch("match_DsFit_phi_eta",  &match_DsFit_phi_vec.eta);
    tree->Branch("match_DsFit_phi_phi",  &match_DsFit_phi_vec.phi);
    tree->Branch("match_DsFit_phi_p",    &match_DsFit_phi_vec.p);
    tree->Branch("match_DsFit_phi_pt",   &match_DsFit_phi_vec.pt);
    tree->Branch("match_DsFit_phi_px",   &match_DsFit_phi_vec.px);
    tree->Branch("match_DsFit_phi_py",   &match_DsFit_phi_vec.py);
    tree->Branch("match_DsFit_phi_pz",   &match_DsFit_phi_vec.pz);
    tree->Branch("match_DsFit_phi_invm", &match_DsFit_phi_vec.invm);

    tree->Branch("match_DsFit_Ds_eta",  &match_DsFit_Ds_vec.eta);
    tree->Branch("match_DsFit_Ds_phi",  &match_DsFit_Ds_vec.phi);
    tree->Branch("match_DsFit_Ds_p",    &match_DsFit_Ds_vec.p);
    tree->Branch("match_DsFit_Ds_pt",   &match_DsFit_Ds_vec.pt);
    tree->Branch("match_DsFit_Ds_px",   &match_DsFit_Ds_vec.px);
    tree->Branch("match_DsFit_Ds_py",   &match_DsFit_Ds_vec.py);
    tree->Branch("match_DsFit_Ds_pz",   &match_DsFit_Ds_vec.pz);
    tree->Branch("match_DsFit_Ds_invm", &match_DsFit_Ds_vec.invm);

    tree->Branch("match_dR_Kp_Km",  &match_dR_Kp_Km_vec);
    tree->Branch("match_dR_Kp_phi", &match_dR_Kp_phi_vec);
    tree->Branch("match_dR_Km_phi", &match_dR_Km_phi_vec);
    tree->Branch("match_dR_Kp_pi",  &match_dR_Kp_pi_vec);
    tree->Branch("match_dR_Km_pi",  &match_dR_Km_pi_vec);
    tree->Branch("match_dR_pi_phi", &match_dR_pi_phi_vec);
    tree->Branch("match_dR_Kp_Ds",  &match_dR_Kp_Ds_vec);
    tree->Branch("match_dR_Km_Ds",  &match_dR_Km_Ds_vec);
    tree->Branch("match_dR_phi_Ds", &match_dR_phi_Ds_vec);
    tree->Branch("match_dR_pi_Ds",  &match_dR_pi_Ds_vec);
    
    tree->Branch("match_phiFit_dR_Kp_Km",  &match_phiFit_dR_Kp_Km_vec);
    tree->Branch("match_phiFit_dR_Kp_phi", &match_phiFit_dR_Kp_phi_vec);
    tree->Branch("match_phiFit_dR_Km_phi", &match_phiFit_dR_Km_phi_vec);
    tree->Branch("match_phiFit_dR_Kp_pi",  &match_phiFit_dR_Kp_pi_vec);
    tree->Branch("match_phiFit_dR_Km_pi",  &match_phiFit_dR_Km_pi_vec);
    tree->Branch("match_phiFit_dR_pi_phi", &match_phiFit_dR_pi_phi_vec);
    tree->Branch("match_phiFit_dR_Kp_Ds",  &match_phiFit_dR_Kp_Ds_vec);
    tree->Branch("match_phiFit_dR_Km_Ds",  &match_phiFit_dR_Km_Ds_vec);
    tree->Branch("match_phiFit_dR_phi_Ds", &match_phiFit_dR_phi_Ds_vec);
    tree->Branch("match_phiFit_dR_pi_Ds",  &match_phiFit_dR_pi_Ds_vec);

    tree->Branch("match_DsFit_dR_Kp_Km",  &match_DsFit_dR_Kp_Km_vec);
    tree->Branch("match_DsFit_dR_Kp_phi", &match_DsFit_dR_Kp_phi_vec);
    tree->Branch("match_DsFit_dR_Km_phi", &match_DsFit_dR_Km_phi_vec);
    tree->Branch("match_DsFit_dR_Kp_pi",  &match_DsFit_dR_Kp_pi_vec);
    tree->Branch("match_DsFit_dR_Km_pi",  &match_DsFit_dR_Km_pi_vec);
    tree->Branch("match_DsFit_dR_pi_phi", &match_DsFit_dR_pi_phi_vec);
    tree->Branch("match_DsFit_dR_Kp_Ds",  &match_DsFit_dR_Kp_Ds_vec);
    tree->Branch("match_DsFit_dR_Km_Ds",  &match_DsFit_dR_Km_Ds_vec);
    tree->Branch("match_DsFit_dR_phi_Ds", &match_DsFit_dR_phi_Ds_vec);
    tree->Branch("match_DsFit_dR_pi_Ds",  &match_DsFit_dR_pi_Ds_vec);

    tree->Branch("match_dxy_phi_Ds", &match_dxy_phi_Ds_vec);
    tree->Branch("match_dz_phi_Ds",  &match_dz_phi_Ds_vec);

    tree->Branch("match_DsFit_Mconstraint_Ds_invm", &match_DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("match_phiFit_isValid",  &match_phiFit_vec.isValid);
    tree->Branch("match_phiFit_isFake",   &match_phiFit_vec.isFake);
    tree->Branch("match_phiFit_chi2",     &match_phiFit_vec.chi2);
    tree->Branch("match_phiFit_ndof",     &match_phiFit_vec.ndof);
    tree->Branch("match_phiFit_chi2ndof", &match_phiFit_vec.chi2ndof);
    tree->Branch("match_phiFit_vx",       &match_phiFit_vec.vx);
    tree->Branch("match_phiFit_vy",       &match_phiFit_vec.vy);
    tree->Branch("match_phiFit_vz",       &match_phiFit_vec.vz);
    tree->Branch("match_phiFit_vxerr",    &match_phiFit_vec.vxerr);
    tree->Branch("match_phiFit_vyerr",    &match_phiFit_vec.vyerr);
    tree->Branch("match_phiFit_vzerr",    &match_phiFit_vec.vzerr);

    tree->Branch("match_DsFit_isValid",  &match_DsFit_vec.isValid);
    tree->Branch("match_DsFit_isFake",   &match_DsFit_vec.isFake);
    tree->Branch("match_DsFit_chi2",     &match_DsFit_vec.chi2);
    tree->Branch("match_DsFit_ndof",     &match_DsFit_vec.ndof);
    tree->Branch("match_DsFit_chi2ndof", &match_DsFit_vec.chi2ndof);
    tree->Branch("match_DsFit_vx",       &match_DsFit_vec.vx);
    tree->Branch("match_DsFit_vy",       &match_DsFit_vec.vy);
    tree->Branch("match_DsFit_vz",       &match_DsFit_vec.vz);
    tree->Branch("match_DsFit_vxerr",    &match_DsFit_vec.vxerr);
    tree->Branch("match_DsFit_vyerr",    &match_DsFit_vec.vyerr);
    tree->Branch("match_DsFit_vzerr",    &match_DsFit_vec.vzerr);

    tree->Branch("match_PVnoDs_isValid",  &match_PVnoDs_vec.isValid);
    tree->Branch("match_PVnoDs_isFake",   &match_PVnoDs_vec.isFake);
    tree->Branch("match_PVnoDs_chi2",     &match_PVnoDs_vec.chi2);
    tree->Branch("match_PVnoDs_ndof",     &match_PVnoDs_vec.ndof);
    tree->Branch("match_PVnoDs_chi2ndof", &match_PVnoDs_vec.chi2ndof);
    tree->Branch("match_PVnoDs_vx",       &match_PVnoDs_vec.vx);
    tree->Branch("match_PVnoDs_vy",       &match_PVnoDs_vec.vy);
    tree->Branch("match_PVnoDs_vz",       &match_PVnoDs_vec.vz);
    tree->Branch("match_PVnoDs_vxerr",    &match_PVnoDs_vec.vxerr);
    tree->Branch("match_PVnoDs_vyerr",    &match_PVnoDs_vec.vyerr);
    tree->Branch("match_PVnoDs_vzerr",    &match_PVnoDs_vec.vzerr);

    tree->Branch("match_PVwithDs_isValid",  &match_PVwithDs_vec.isValid);
    tree->Branch("match_PVwithDs_isFake",   &match_PVwithDs_vec.isFake);
    tree->Branch("match_PVwithDs_chi2",     &match_PVwithDs_vec.chi2);
    tree->Branch("match_PVwithDs_ndof",     &match_PVwithDs_vec.ndof);
    tree->Branch("match_PVwithDs_chi2ndof", &match_PVwithDs_vec.chi2ndof);
    tree->Branch("match_PVwithDs_vx",       &match_PVwithDs_vec.vx);
    tree->Branch("match_PVwithDs_vy",       &match_PVwithDs_vec.vy);
    tree->Branch("match_PVwithDs_vz",       &match_PVwithDs_vec.vz);
    tree->Branch("match_PVwithDs_vxerr",    &match_PVwithDs_vec.vxerr);
    tree->Branch("match_PVwithDs_vyerr",    &match_PVwithDs_vec.vyerr);
    tree->Branch("match_PVwithDs_vzerr",    &match_PVwithDs_vec.vzerr);

    tree->Branch("match_Ds_primvtx_FDxy",       &match_Ds_primvtx_FD_vec.FDxy);
    tree->Branch("match_Ds_primvtx_FDz",        &match_Ds_primvtx_FD_vec.FDz);
    tree->Branch("match_Ds_primvtx_FD",         &match_Ds_primvtx_FD_vec.FD);
    tree->Branch("match_Ds_primvtx_FDxyerr",    &match_Ds_primvtx_FD_vec.FDxyerr);
    tree->Branch("match_Ds_primvtx_FDzerr",     &match_Ds_primvtx_FD_vec.FDzerr);
    tree->Branch("match_Ds_primvtx_FDerr",      &match_Ds_primvtx_FD_vec.FDerr);
    tree->Branch("match_Ds_primvtx_FDxychi2",   &match_Ds_primvtx_FD_vec.FDxychi2);
    tree->Branch("match_Ds_primvtx_FDzchi2",    &match_Ds_primvtx_FD_vec.FDzchi2);
    tree->Branch("match_Ds_primvtx_FDchi2",     &match_Ds_primvtx_FD_vec.FDchi2);
    tree->Branch("match_Ds_primvtx_dira",       &match_Ds_primvtx_FD_vec.dira);
    tree->Branch("match_Ds_primvtx_dira_angle", &match_Ds_primvtx_FD_vec.dira_angle);
    tree->Branch("match_Kp_primvtx_ip",         &match_Kp_primvtx_ip_vec.ip);
    tree->Branch("match_Kp_primvtx_iperr",      &match_Kp_primvtx_ip_vec.iperr);
    tree->Branch("match_Kp_primvtx_ipchi2",     &match_Kp_primvtx_ip_vec.ipchi2);
    tree->Branch("match_Km_primvtx_ip",         &match_Km_primvtx_ip_vec.ip);
    tree->Branch("match_Km_primvtx_iperr",      &match_Km_primvtx_ip_vec.iperr);
    tree->Branch("match_Km_primvtx_ipchi2",     &match_Km_primvtx_ip_vec.ipchi2);
    tree->Branch("match_pi_primvtx_ip",         &match_pi_primvtx_ip_vec.ip);
    tree->Branch("match_pi_primvtx_iperr",      &match_pi_primvtx_ip_vec.iperr);
    tree->Branch("match_pi_primvtx_ipchi2",     &match_pi_primvtx_ip_vec.ipchi2);
    tree->Branch("match_phi_primvtx_ip",        &match_phi_primvtx_ip_vec.ip);
    tree->Branch("match_phi_primvtx_iperr",     &match_phi_primvtx_ip_vec.iperr);
    tree->Branch("match_phi_primvtx_ipchi2",    &match_phi_primvtx_ip_vec.ipchi2);
    tree->Branch("match_Ds_primvtx_ip",         &match_Ds_primvtx_ip_vec.ip);
    tree->Branch("match_Ds_primvtx_iperr",      &match_Ds_primvtx_ip_vec.iperr);
    tree->Branch("match_Ds_primvtx_ipchi2",     &match_Ds_primvtx_ip_vec.ipchi2);

    tree->Branch("match_Ds_PVnoDs_FDxy",       &match_Ds_PVnoDs_FD_vec.FDxy);
    tree->Branch("match_Ds_PVnoDs_FDz",        &match_Ds_PVnoDs_FD_vec.FDz);
    tree->Branch("match_Ds_PVnoDs_FD",         &match_Ds_PVnoDs_FD_vec.FD);
    tree->Branch("match_Ds_PVnoDs_FDxyerr",    &match_Ds_PVnoDs_FD_vec.FDxyerr);
    tree->Branch("match_Ds_PVnoDs_FDzerr",     &match_Ds_PVnoDs_FD_vec.FDzerr);
    tree->Branch("match_Ds_PVnoDs_FDerr",      &match_Ds_PVnoDs_FD_vec.FDerr);
    tree->Branch("match_Ds_PVnoDs_FDxychi2",   &match_Ds_PVnoDs_FD_vec.FDxychi2);
    tree->Branch("match_Ds_PVnoDs_FDzchi2",    &match_Ds_PVnoDs_FD_vec.FDzchi2);
    tree->Branch("match_Ds_PVnoDs_FDchi2",     &match_Ds_PVnoDs_FD_vec.FDchi2);
    tree->Branch("match_Ds_PVnoDs_dira",       &match_Ds_PVnoDs_FD_vec.dira);
    tree->Branch("match_Ds_PVnoDs_dira_angle", &match_Ds_PVnoDs_FD_vec.dira_angle);
    tree->Branch("match_Kp_PVnoDs_ip",         &match_Kp_PVnoDs_ip_vec.ip);
    tree->Branch("match_Kp_PVnoDs_iperr",      &match_Kp_PVnoDs_ip_vec.iperr);
    tree->Branch("match_Kp_PVnoDs_ipchi2",     &match_Kp_PVnoDs_ip_vec.ipchi2);
    tree->Branch("match_Km_PVnoDs_ip",         &match_Km_PVnoDs_ip_vec.ip);
    tree->Branch("match_Km_PVnoDs_iperr",      &match_Km_PVnoDs_ip_vec.iperr);
    tree->Branch("match_Km_PVnoDs_ipchi2",     &match_Km_PVnoDs_ip_vec.ipchi2);
    tree->Branch("match_pi_PVnoDs_ip",         &match_pi_PVnoDs_ip_vec.ip);
    tree->Branch("match_pi_PVnoDs_iperr",      &match_pi_PVnoDs_ip_vec.iperr);
    tree->Branch("match_pi_PVnoDs_ipchi2",     &match_pi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("match_phi_PVnoDs_ip",        &match_phi_PVnoDs_ip_vec.ip);
    tree->Branch("match_phi_PVnoDs_iperr",     &match_phi_PVnoDs_ip_vec.iperr);
    tree->Branch("match_phi_PVnoDs_ipchi2",    &match_phi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("match_Ds_PVnoDs_ip",         &match_Ds_PVnoDs_ip_vec.ip);
    tree->Branch("match_Ds_PVnoDs_iperr",      &match_Ds_PVnoDs_ip_vec.iperr);
    tree->Branch("match_Ds_PVnoDs_ipchi2",     &match_Ds_PVnoDs_ip_vec.ipchi2);

    tree->Branch("match_Ds_PVwithDs_FDxy",       &match_Ds_PVwithDs_FD_vec.FDxy);
    tree->Branch("match_Ds_PVwithDs_FDz",        &match_Ds_PVwithDs_FD_vec.FDz);
    tree->Branch("match_Ds_PVwithDs_FD",         &match_Ds_PVwithDs_FD_vec.FD);
    tree->Branch("match_Ds_PVwithDs_FDxyerr",    &match_Ds_PVwithDs_FD_vec.FDxyerr);
    tree->Branch("match_Ds_PVwithDs_FDzerr",     &match_Ds_PVwithDs_FD_vec.FDzerr);
    tree->Branch("match_Ds_PVwithDs_FDerr",      &match_Ds_PVwithDs_FD_vec.FDerr);
    tree->Branch("match_Ds_PVwithDs_FDxychi2",   &match_Ds_PVwithDs_FD_vec.FDxychi2);
    tree->Branch("match_Ds_PVwithDs_FDzchi2",    &match_Ds_PVwithDs_FD_vec.FDzchi2);
    tree->Branch("match_Ds_PVwithDs_FDchi2",     &match_Ds_PVwithDs_FD_vec.FDchi2);
    tree->Branch("match_Ds_PVwithDs_dira",       &match_Ds_PVwithDs_FD_vec.dira);
    tree->Branch("match_Ds_PVwithDs_dira_angle", &match_Ds_PVwithDs_FD_vec.dira_angle);
    tree->Branch("match_Kp_PVwithDs_ip",         &match_Kp_PVwithDs_ip_vec.ip);
    tree->Branch("match_Kp_PVwithDs_iperr",      &match_Kp_PVwithDs_ip_vec.iperr);
    tree->Branch("match_Kp_PVwithDs_ipchi2",     &match_Kp_PVwithDs_ip_vec.ipchi2);
    tree->Branch("match_Km_PVwithDs_ip",         &match_Km_PVwithDs_ip_vec.ip);
    tree->Branch("match_Km_PVwithDs_iperr",      &match_Km_PVwithDs_ip_vec.iperr);
    tree->Branch("match_Km_PVwithDs_ipchi2",     &match_Km_PVwithDs_ip_vec.ipchi2);
    tree->Branch("match_pi_PVwithDs_ip",         &match_pi_PVwithDs_ip_vec.ip);
    tree->Branch("match_pi_PVwithDs_iperr",      &match_pi_PVwithDs_ip_vec.iperr);
    tree->Branch("match_pi_PVwithDs_ipchi2",     &match_pi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("match_phi_PVwithDs_ip",        &match_phi_PVwithDs_ip_vec.ip);
    tree->Branch("match_phi_PVwithDs_iperr",     &match_phi_PVwithDs_ip_vec.iperr);
    tree->Branch("match_phi_PVwithDs_ipchi2",    &match_phi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("match_Ds_PVwithDs_ip",         &match_Ds_PVwithDs_ip_vec.ip);
    tree->Branch("match_Ds_PVwithDs_iperr",      &match_Ds_PVwithDs_ip_vec.iperr);
    tree->Branch("match_Ds_PVwithDs_ipchi2",     &match_Ds_PVwithDs_ip_vec.ipchi2);

    tree->Branch("match_Ds_IsoR03_sumChargedHadronPt", &match_Ds_IsoR03_vec.sumChargedHadronPt);
    tree->Branch("match_Ds_IsoR03_sumNeutralHadronEt", &match_Ds_IsoR03_vec.sumNeutralHadronEt);
    tree->Branch("match_Ds_IsoR03_sumPhotonEt",        &match_Ds_IsoR03_vec.sumPhotonEt);
    tree->Branch("match_Ds_IsoR03_sumPUPt",            &match_Ds_IsoR03_vec.sumPUPt);
    tree->Branch("match_Ds_IsoR03_PFIso",              &match_Ds_IsoR03_vec.PFIso);

    tree->Branch("match_Ds_IsoR04_sumChargedHadronPt", &match_Ds_IsoR04_vec.sumChargedHadronPt);
    tree->Branch("match_Ds_IsoR04_sumNeutralHadronEt", &match_Ds_IsoR04_vec.sumNeutralHadronEt);
    tree->Branch("match_Ds_IsoR04_sumPhotonEt",        &match_Ds_IsoR04_vec.sumPhotonEt);
    tree->Branch("match_Ds_IsoR04_sumPUPt",            &match_Ds_IsoR04_vec.sumPUPt);
    tree->Branch("match_Ds_IsoR04_PFIso",              &match_Ds_IsoR04_vec.PFIso);

    tree->Branch("match_mu_charge",    &match_mu_vec.charge);
    tree->Branch("match_mu_eta",       &match_mu_vec.eta);
    tree->Branch("match_mu_phi",       &match_mu_vec.phi);
    tree->Branch("match_mu_vx",        &match_mu_vec.vx);
    tree->Branch("match_mu_vy",        &match_mu_vec.vy);
    tree->Branch("match_mu_vz",        &match_mu_vec.vz);
    tree->Branch("match_mu_p",         &match_mu_vec.p);
    tree->Branch("match_mu_pt",        &match_mu_vec.pt);
    tree->Branch("match_mu_px",        &match_mu_vec.px);
    tree->Branch("match_mu_py",        &match_mu_vec.py);
    tree->Branch("match_mu_pz",        &match_mu_vec.pz);
    tree->Branch("match_mu_dxy",       &match_mu_vec.dxy);
    tree->Branch("match_mu_dxyerr",    &match_mu_vec.dxyerr);
    tree->Branch("match_mu_dz",        &match_mu_vec.dz);
    tree->Branch("match_mu_dzerr",     &match_mu_vec.dzerr);
    tree->Branch("match_mu_isHighPt",  &match_mu_vec.isHighPt);
    tree->Branch("match_mu_isLoose",   &match_mu_vec.isLoose);
    tree->Branch("match_mu_isMedium",  &match_mu_vec.isMedium);
    tree->Branch("match_mu_isSoft",    &match_mu_vec.isSoft);
    tree->Branch("match_mu_isTight",   &match_mu_vec.isTight);
    tree->Branch("match_mu_isPF",      &match_mu_vec.isPF);
    tree->Branch("match_mu_isTracker", &match_mu_vec.isTracker);
    tree->Branch("match_mu_isGlobal",  &match_mu_vec.isGlobal);

    tree->Branch("match_mu_IsoR03_sumChargedHadronPt", &match_mu_IsoR03_vec.sumChargedHadronPt);
    tree->Branch("match_mu_IsoR03_sumNeutralHadronEt", &match_mu_IsoR03_vec.sumNeutralHadronEt);
    tree->Branch("match_mu_IsoR03_sumPhotonEt",        &match_mu_IsoR03_vec.sumPhotonEt);
    tree->Branch("match_mu_IsoR03_sumPUPt",            &match_mu_IsoR03_vec.sumPUPt);
    tree->Branch("match_mu_IsoR03_PFIso",              &match_mu_IsoR03_vec.PFIso);

    tree->Branch("match_mu_IsoR04_sumChargedHadronPt", &match_mu_IsoR04_vec.sumChargedHadronPt);
    tree->Branch("match_mu_IsoR04_sumNeutralHadronEt", &match_mu_IsoR04_vec.sumNeutralHadronEt);
    tree->Branch("match_mu_IsoR04_sumPhotonEt",        &match_mu_IsoR04_vec.sumPhotonEt);
    tree->Branch("match_mu_IsoR04_sumPUPt",            &match_mu_IsoR04_vec.sumPUPt);
    tree->Branch("match_mu_IsoR04_PFIso",              &match_mu_IsoR04_vec.PFIso);

    tree->Branch("match_mu_primvtx_ip",         &match_mu_primvtx_ip_vec.ip);
    tree->Branch("match_mu_primvtx_iperr",      &match_mu_primvtx_ip_vec.iperr);
    tree->Branch("match_mu_primvtx_ipchi2",     &match_mu_primvtx_ip_vec.ipchi2);

    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds",  &num_reco_Ds);

    tree->Branch("Kp_isIsolatedChargedHadron", &Kp_vec.isIsolatedChargedHadron);
    tree->Branch("Kp_charge",                  &Kp_vec.charge);
    tree->Branch("Kp_eta",                     &Kp_vec.eta);
    tree->Branch("Kp_phi",                     &Kp_vec.phi);
    tree->Branch("Kp_vx",                      &Kp_vec.vx);
    tree->Branch("Kp_vy",                      &Kp_vec.vy);
    tree->Branch("Kp_vz",                      &Kp_vec.vz);
    tree->Branch("Kp_p",                       &Kp_vec.p);
    tree->Branch("Kp_pt",                      &Kp_vec.pt);
    tree->Branch("Kp_px",                      &Kp_vec.px);
    tree->Branch("Kp_py",                      &Kp_vec.py);
    tree->Branch("Kp_pz",                      &Kp_vec.pz);

    tree->Branch("Km_isIsolatedChargedHadron", &Km_vec.isIsolatedChargedHadron);
    tree->Branch("Km_charge",                  &Km_vec.charge);
    tree->Branch("Km_eta",                     &Km_vec.eta);
    tree->Branch("Km_phi",                     &Km_vec.phi);
    tree->Branch("Km_vx",                      &Km_vec.vx);
    tree->Branch("Km_vy",                      &Km_vec.vy);
    tree->Branch("Km_vz",                      &Km_vec.vz);
    tree->Branch("Km_p",                       &Km_vec.p);
    tree->Branch("Km_pt",                      &Km_vec.pt);
    tree->Branch("Km_px",                      &Km_vec.px);
    tree->Branch("Km_py",                      &Km_vec.py);
    tree->Branch("Km_pz",                      &Km_vec.pz);

    tree->Branch("pi_isIsolatedChargedHadron", &pi_vec.isIsolatedChargedHadron);
    tree->Branch("pi_charge",                  &pi_vec.charge);
    tree->Branch("pi_eta",                     &pi_vec.eta);
    tree->Branch("pi_phi",                     &pi_vec.phi);
    tree->Branch("pi_vx",                      &pi_vec.vx);
    tree->Branch("pi_vy",                      &pi_vec.vy);
    tree->Branch("pi_vz",                      &pi_vec.vz);
    tree->Branch("pi_p",                       &pi_vec.p);
    tree->Branch("pi_pt",                      &pi_vec.pt);
    tree->Branch("pi_px",                      &pi_vec.px);
    tree->Branch("pi_py",                      &pi_vec.py);
    tree->Branch("pi_pz",                      &pi_vec.pz);

    tree->Branch("phi_eta",  &phi_vec.eta);
    tree->Branch("phi_phi",  &phi_vec.phi);
    tree->Branch("phi_p",    &phi_vec.p);
    tree->Branch("phi_pt",   &phi_vec.pt);
    tree->Branch("phi_px",   &phi_vec.px);
    tree->Branch("phi_py",   &phi_vec.py);
    tree->Branch("phi_pz",   &phi_vec.pz);
    tree->Branch("phi_invm", &phi_vec.invm);

    tree->Branch("Ds_eta",  &Ds_vec.eta);
    tree->Branch("Ds_phi",  &Ds_vec.phi);
    tree->Branch("Ds_p",    &Ds_vec.p);
    tree->Branch("Ds_pt",   &Ds_vec.pt);
    tree->Branch("Ds_px",   &Ds_vec.px);
    tree->Branch("Ds_py",   &Ds_vec.py);
    tree->Branch("Ds_pz",   &Ds_vec.pz);
    tree->Branch("Ds_invm", &Ds_vec.invm);

    tree->Branch("phiFit_Kp_eta",  &phiFit_Kp_vec.eta);
    tree->Branch("phiFit_Kp_phi",  &phiFit_Kp_vec.phi);
    tree->Branch("phiFit_Kp_p",    &phiFit_Kp_vec.p);
    tree->Branch("phiFit_Kp_pt",   &phiFit_Kp_vec.pt);
    tree->Branch("phiFit_Kp_px",   &phiFit_Kp_vec.px);
    tree->Branch("phiFit_Kp_py",   &phiFit_Kp_vec.py);
    tree->Branch("phiFit_Kp_pz",   &phiFit_Kp_vec.pz);
    tree->Branch("phiFit_Kp_invm", &phiFit_Kp_vec.invm);

    tree->Branch("phiFit_Km_eta",  &phiFit_Km_vec.eta);
    tree->Branch("phiFit_Km_phi",  &phiFit_Km_vec.phi);
    tree->Branch("phiFit_Km_p",    &phiFit_Km_vec.p);
    tree->Branch("phiFit_Km_pt",   &phiFit_Km_vec.pt);
    tree->Branch("phiFit_Km_px",   &phiFit_Km_vec.px);
    tree->Branch("phiFit_Km_py",   &phiFit_Km_vec.py);
    tree->Branch("phiFit_Km_pz",   &phiFit_Km_vec.pz);
    tree->Branch("phiFit_Km_invm", &phiFit_Km_vec.invm);

    tree->Branch("phiFit_pi_eta",  &phiFit_pi_vec.eta);
    tree->Branch("phiFit_pi_phi",  &phiFit_pi_vec.phi);
    tree->Branch("phiFit_pi_p",    &phiFit_pi_vec.p);
    tree->Branch("phiFit_pi_pt",   &phiFit_pi_vec.pt);
    tree->Branch("phiFit_pi_px",   &phiFit_pi_vec.px);
    tree->Branch("phiFit_pi_py",   &phiFit_pi_vec.py);
    tree->Branch("phiFit_pi_pz",   &phiFit_pi_vec.pz);
    tree->Branch("phiFit_pi_invm", &phiFit_pi_vec.invm);

    tree->Branch("phiFit_phi_eta",  &phiFit_phi_vec.eta);
    tree->Branch("phiFit_phi_phi",  &phiFit_phi_vec.phi);
    tree->Branch("phiFit_phi_p",    &phiFit_phi_vec.p);
    tree->Branch("phiFit_phi_pt",   &phiFit_phi_vec.pt);
    tree->Branch("phiFit_phi_px",   &phiFit_phi_vec.px);
    tree->Branch("phiFit_phi_py",   &phiFit_phi_vec.py);
    tree->Branch("phiFit_phi_pz",   &phiFit_phi_vec.pz);
    tree->Branch("phiFit_phi_invm", &phiFit_phi_vec.invm);

    tree->Branch("phiFit_Ds_eta",  &phiFit_Ds_vec.eta);
    tree->Branch("phiFit_Ds_phi",  &phiFit_Ds_vec.phi);
    tree->Branch("phiFit_Ds_p",    &phiFit_Ds_vec.p);
    tree->Branch("phiFit_Ds_pt",   &phiFit_Ds_vec.pt);
    tree->Branch("phiFit_Ds_px",   &phiFit_Ds_vec.px);
    tree->Branch("phiFit_Ds_py",   &phiFit_Ds_vec.py);
    tree->Branch("phiFit_Ds_pz",   &phiFit_Ds_vec.pz);
    tree->Branch("phiFit_Ds_invm", &phiFit_Ds_vec.invm);

    tree->Branch("DsFit_Kp_eta",  &DsFit_Kp_vec.eta);
    tree->Branch("DsFit_Kp_phi",  &DsFit_Kp_vec.phi);
    tree->Branch("DsFit_Kp_p",    &DsFit_Kp_vec.p);
    tree->Branch("DsFit_Kp_pt",   &DsFit_Kp_vec.pt);
    tree->Branch("DsFit_Kp_px",   &DsFit_Kp_vec.px);
    tree->Branch("DsFit_Kp_py",   &DsFit_Kp_vec.py);
    tree->Branch("DsFit_Kp_pz",   &DsFit_Kp_vec.pz);
    tree->Branch("DsFit_Kp_invm", &DsFit_Kp_vec.invm);

    tree->Branch("DsFit_Km_eta",  &DsFit_Km_vec.eta);
    tree->Branch("DsFit_Km_phi",  &DsFit_Km_vec.phi);
    tree->Branch("DsFit_Km_p",    &DsFit_Km_vec.p);
    tree->Branch("DsFit_Km_pt",   &DsFit_Km_vec.pt);
    tree->Branch("DsFit_Km_px",   &DsFit_Km_vec.px);
    tree->Branch("DsFit_Km_py",   &DsFit_Km_vec.py);
    tree->Branch("DsFit_Km_pz",   &DsFit_Km_vec.pz);
    tree->Branch("DsFit_Km_invm", &DsFit_Km_vec.invm);

    tree->Branch("DsFit_pi_eta",  &DsFit_pi_vec.eta);
    tree->Branch("DsFit_pi_phi",  &DsFit_pi_vec.phi);
    tree->Branch("DsFit_pi_p",    &DsFit_pi_vec.p);
    tree->Branch("DsFit_pi_pt",   &DsFit_pi_vec.pt);
    tree->Branch("DsFit_pi_px",   &DsFit_pi_vec.px);
    tree->Branch("DsFit_pi_py",   &DsFit_pi_vec.py);
    tree->Branch("DsFit_pi_pz",   &DsFit_pi_vec.pz);
    tree->Branch("DsFit_pi_invm", &DsFit_pi_vec.invm);

    tree->Branch("DsFit_phi_eta",  &DsFit_phi_vec.eta);
    tree->Branch("DsFit_phi_phi",  &DsFit_phi_vec.phi);
    tree->Branch("DsFit_phi_p",    &DsFit_phi_vec.p);
    tree->Branch("DsFit_phi_pt",   &DsFit_phi_vec.pt);
    tree->Branch("DsFit_phi_px",   &DsFit_phi_vec.px);
    tree->Branch("DsFit_phi_py",   &DsFit_phi_vec.py);
    tree->Branch("DsFit_phi_pz",   &DsFit_phi_vec.pz);
    tree->Branch("DsFit_phi_invm", &DsFit_phi_vec.invm);

    tree->Branch("DsFit_Ds_eta",  &DsFit_Ds_vec.eta);
    tree->Branch("DsFit_Ds_phi",  &DsFit_Ds_vec.phi);
    tree->Branch("DsFit_Ds_p",    &DsFit_Ds_vec.p);
    tree->Branch("DsFit_Ds_pt",   &DsFit_Ds_vec.pt);
    tree->Branch("DsFit_Ds_px",   &DsFit_Ds_vec.px);
    tree->Branch("DsFit_Ds_py",   &DsFit_Ds_vec.py);
    tree->Branch("DsFit_Ds_pz",   &DsFit_Ds_vec.pz);
    tree->Branch("DsFit_Ds_invm", &DsFit_Ds_vec.invm);

    tree->Branch("dR_Kp_Km",  &dR_Kp_Km_vec);
    tree->Branch("dR_Kp_phi", &dR_Kp_phi_vec);
    tree->Branch("dR_Km_phi", &dR_Km_phi_vec);
    tree->Branch("dR_Kp_pi",  &dR_Kp_pi_vec);
    tree->Branch("dR_Km_pi",  &dR_Km_pi_vec);
    tree->Branch("dR_pi_phi", &dR_pi_phi_vec);
    tree->Branch("dR_Kp_Ds",  &dR_Kp_Ds_vec);
    tree->Branch("dR_Km_Ds",  &dR_Km_Ds_vec);
    tree->Branch("dR_phi_Ds", &dR_phi_Ds_vec);
    tree->Branch("dR_pi_Ds",  &dR_pi_Ds_vec);
    
    tree->Branch("phiFit_dR_Kp_Km",  &phiFit_dR_Kp_Km_vec);
    tree->Branch("phiFit_dR_Kp_phi", &phiFit_dR_Kp_phi_vec);
    tree->Branch("phiFit_dR_Km_phi", &phiFit_dR_Km_phi_vec);
    tree->Branch("phiFit_dR_Kp_pi",  &phiFit_dR_Kp_pi_vec);
    tree->Branch("phiFit_dR_Km_pi",  &phiFit_dR_Km_pi_vec);
    tree->Branch("phiFit_dR_pi_phi", &phiFit_dR_pi_phi_vec);
    tree->Branch("phiFit_dR_Kp_Ds",  &phiFit_dR_Kp_Ds_vec);
    tree->Branch("phiFit_dR_Km_Ds",  &phiFit_dR_Km_Ds_vec);
    tree->Branch("phiFit_dR_phi_Ds", &phiFit_dR_phi_Ds_vec);
    tree->Branch("phiFit_dR_pi_Ds",  &phiFit_dR_pi_Ds_vec);


    tree->Branch("DsFit_dR_Kp_Km",  &DsFit_dR_Kp_Km_vec);
    tree->Branch("DsFit_dR_Kp_phi", &DsFit_dR_Kp_phi_vec);
    tree->Branch("DsFit_dR_Km_phi", &DsFit_dR_Km_phi_vec);
    tree->Branch("DsFit_dR_Kp_pi",  &DsFit_dR_Kp_pi_vec);
    tree->Branch("DsFit_dR_Km_pi",  &DsFit_dR_Km_pi_vec);
    tree->Branch("DsFit_dR_pi_phi", &DsFit_dR_pi_phi_vec);
    tree->Branch("DsFit_dR_Kp_Ds",  &DsFit_dR_Kp_Ds_vec);
    tree->Branch("DsFit_dR_Km_Ds",  &DsFit_dR_Km_Ds_vec);
    tree->Branch("DsFit_dR_phi_Ds", &DsFit_dR_phi_Ds_vec);
    tree->Branch("DsFit_dR_pi_Ds",  &DsFit_dR_pi_Ds_vec);

    tree->Branch("dxy_phi_Ds", &dxy_phi_Ds_vec);
    tree->Branch("dz_phi_Ds",  &dz_phi_Ds_vec);

    tree->Branch("DsFit_Mconstraint_Ds_invm", &DsFit_Mconstraint_Ds_invm_vec);

    tree->Branch("phiFit_isValid",  &phiFit_vec.isValid);
    tree->Branch("phiFit_isFake",   &phiFit_vec.isFake);
    tree->Branch("phiFit_chi2",     &phiFit_vec.chi2);
    tree->Branch("phiFit_ndof",     &phiFit_vec.ndof);
    tree->Branch("phiFit_chi2ndof", &phiFit_vec.chi2ndof);
    tree->Branch("phiFit_vx",       &phiFit_vec.vx);
    tree->Branch("phiFit_vy",       &phiFit_vec.vy);
    tree->Branch("phiFit_vz",       &phiFit_vec.vz);
    tree->Branch("phiFit_vxerr",    &phiFit_vec.vxerr);
    tree->Branch("phiFit_vyerr",    &phiFit_vec.vyerr);
    tree->Branch("phiFit_vzerr",    &phiFit_vec.vzerr);

    tree->Branch("DsFit_isValid",  &DsFit_vec.isValid);
    tree->Branch("DsFit_isFake",   &DsFit_vec.isFake);
    tree->Branch("DsFit_chi2",     &DsFit_vec.chi2);
    tree->Branch("DsFit_ndof",     &DsFit_vec.ndof);
    tree->Branch("DsFit_chi2ndof", &DsFit_vec.chi2ndof);
    tree->Branch("DsFit_vx",       &DsFit_vec.vx);
    tree->Branch("DsFit_vy",       &DsFit_vec.vy);
    tree->Branch("DsFit_vz",       &DsFit_vec.vz);
    tree->Branch("DsFit_vxerr",    &DsFit_vec.vxerr);
    tree->Branch("DsFit_vyerr",    &DsFit_vec.vyerr);
    tree->Branch("DsFit_vzerr",    &DsFit_vec.vzerr);

    tree->Branch("PVnoDs_isValid",  &PVnoDs_vec.isValid);
    tree->Branch("PVnoDs_isFake",   &PVnoDs_vec.isFake);
    tree->Branch("PVnoDs_chi2",     &PVnoDs_vec.chi2);
    tree->Branch("PVnoDs_ndof",     &PVnoDs_vec.ndof);
    tree->Branch("PVnoDs_chi2ndof", &PVnoDs_vec.chi2ndof);
    tree->Branch("PVnoDs_vx",       &PVnoDs_vec.vx);
    tree->Branch("PVnoDs_vy",       &PVnoDs_vec.vy);
    tree->Branch("PVnoDs_vz",       &PVnoDs_vec.vz);
    tree->Branch("PVnoDs_vxerr",    &PVnoDs_vec.vxerr);
    tree->Branch("PVnoDs_vyerr",    &PVnoDs_vec.vyerr);
    tree->Branch("PVnoDs_vzerr",    &PVnoDs_vec.vzerr);

    tree->Branch("PVwithDs_isValid",  &PVwithDs_vec.isValid);
    tree->Branch("PVwithDs_isFake",   &PVwithDs_vec.isFake);
    tree->Branch("PVwithDs_chi2",     &PVwithDs_vec.chi2);
    tree->Branch("PVwithDs_ndof",     &PVwithDs_vec.ndof);
    tree->Branch("PVwithDs_chi2ndof", &PVwithDs_vec.chi2ndof);
    tree->Branch("PVwithDs_vx",       &PVwithDs_vec.vx);
    tree->Branch("PVwithDs_vy",       &PVwithDs_vec.vy);
    tree->Branch("PVwithDs_vz",       &PVwithDs_vec.vz);
    tree->Branch("PVwithDs_vxerr",    &PVwithDs_vec.vxerr);
    tree->Branch("PVwithDs_vyerr",    &PVwithDs_vec.vyerr);
    tree->Branch("PVwithDs_vzerr",    &PVwithDs_vec.vzerr);

    tree->Branch("Ds_primvtx_FDxy",       &Ds_primvtx_FD_vec.FDxy);
    tree->Branch("Ds_primvtx_FDz",        &Ds_primvtx_FD_vec.FDz);
    tree->Branch("Ds_primvtx_FD",         &Ds_primvtx_FD_vec.FD);
    tree->Branch("Ds_primvtx_FDxyerr",    &Ds_primvtx_FD_vec.FDxyerr);
    tree->Branch("Ds_primvtx_FDzerr",     &Ds_primvtx_FD_vec.FDzerr);
    tree->Branch("Ds_primvtx_FDerr",      &Ds_primvtx_FD_vec.FDerr);
    tree->Branch("Ds_primvtx_FDxychi2",   &Ds_primvtx_FD_vec.FDxychi2);
    tree->Branch("Ds_primvtx_FDzchi2",    &Ds_primvtx_FD_vec.FDzchi2);
    tree->Branch("Ds_primvtx_FDchi2",     &Ds_primvtx_FD_vec.FDchi2);
    tree->Branch("Ds_primvtx_dira",       &Ds_primvtx_FD_vec.dira);
    tree->Branch("Ds_primvtx_dira_angle", &Ds_primvtx_FD_vec.dira_angle);
    tree->Branch("Kp_primvtx_ip",         &Kp_primvtx_ip_vec.ip);
    tree->Branch("Kp_primvtx_iperr",      &Kp_primvtx_ip_vec.iperr);
    tree->Branch("Kp_primvtx_ipchi2",     &Kp_primvtx_ip_vec.ipchi2);
    tree->Branch("Km_primvtx_ip",         &Km_primvtx_ip_vec.ip);
    tree->Branch("Km_primvtx_iperr",      &Km_primvtx_ip_vec.iperr);
    tree->Branch("Km_primvtx_ipchi2",     &Km_primvtx_ip_vec.ipchi2);
    tree->Branch("pi_primvtx_ip",         &pi_primvtx_ip_vec.ip);
    tree->Branch("pi_primvtx_iperr",      &pi_primvtx_ip_vec.iperr);
    tree->Branch("pi_primvtx_ipchi2",     &pi_primvtx_ip_vec.ipchi2);
    tree->Branch("phi_primvtx_ip",        &phi_primvtx_ip_vec.ip);
    tree->Branch("phi_primvtx_iperr",     &phi_primvtx_ip_vec.iperr);
    tree->Branch("phi_primvtx_ipchi2",    &phi_primvtx_ip_vec.ipchi2);
    tree->Branch("Ds_primvtx_ip",         &Ds_primvtx_ip_vec.ip);
    tree->Branch("Ds_primvtx_iperr",      &Ds_primvtx_ip_vec.iperr);
    tree->Branch("Ds_primvtx_ipchi2",     &Ds_primvtx_ip_vec.ipchi2);

    tree->Branch("Ds_PVnoDs_FDxy",       &Ds_PVnoDs_FD_vec.FDxy);
    tree->Branch("Ds_PVnoDs_FDz",        &Ds_PVnoDs_FD_vec.FDz);
    tree->Branch("Ds_PVnoDs_FD",         &Ds_PVnoDs_FD_vec.FD);
    tree->Branch("Ds_PVnoDs_FDxyerr",    &Ds_PVnoDs_FD_vec.FDxyerr);
    tree->Branch("Ds_PVnoDs_FDzerr",     &Ds_PVnoDs_FD_vec.FDzerr);
    tree->Branch("Ds_PVnoDs_FDerr",      &Ds_PVnoDs_FD_vec.FDerr);
    tree->Branch("Ds_PVnoDs_FDxychi2",   &Ds_PVnoDs_FD_vec.FDxychi2);
    tree->Branch("Ds_PVnoDs_FDzchi2",    &Ds_PVnoDs_FD_vec.FDzchi2);
    tree->Branch("Ds_PVnoDs_FDchi2",     &Ds_PVnoDs_FD_vec.FDchi2);
    tree->Branch("Ds_PVnoDs_dira",       &Ds_PVnoDs_FD_vec.dira);
    tree->Branch("Ds_PVnoDs_dira_angle", &Ds_PVnoDs_FD_vec.dira_angle);
    tree->Branch("Kp_PVnoDs_ip",         &Kp_PVnoDs_ip_vec.ip);
    tree->Branch("Kp_PVnoDs_iperr",      &Kp_PVnoDs_ip_vec.iperr);
    tree->Branch("Kp_PVnoDs_ipchi2",     &Kp_PVnoDs_ip_vec.ipchi2);
    tree->Branch("Km_PVnoDs_ip",         &Km_PVnoDs_ip_vec.ip);
    tree->Branch("Km_PVnoDs_iperr",      &Km_PVnoDs_ip_vec.iperr);
    tree->Branch("Km_PVnoDs_ipchi2",     &Km_PVnoDs_ip_vec.ipchi2);
    tree->Branch("pi_PVnoDs_ip",         &pi_PVnoDs_ip_vec.ip);
    tree->Branch("pi_PVnoDs_iperr",      &pi_PVnoDs_ip_vec.iperr);
    tree->Branch("pi_PVnoDs_ipchi2",     &pi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("phi_PVnoDs_ip",        &phi_PVnoDs_ip_vec.ip);
    tree->Branch("phi_PVnoDs_iperr",     &phi_PVnoDs_ip_vec.iperr);
    tree->Branch("phi_PVnoDs_ipchi2",    &phi_PVnoDs_ip_vec.ipchi2);
    tree->Branch("Ds_PVnoDs_ip",         &Ds_PVnoDs_ip_vec.ip);
    tree->Branch("Ds_PVnoDs_iperr",      &Ds_PVnoDs_ip_vec.iperr);
    tree->Branch("Ds_PVnoDs_ipchi2",     &Ds_PVnoDs_ip_vec.ipchi2);

    tree->Branch("Ds_PVwithDs_FDxy",       &Ds_PVwithDs_FD_vec.FDxy);
    tree->Branch("Ds_PVwithDs_FDz",        &Ds_PVwithDs_FD_vec.FDz);
    tree->Branch("Ds_PVwithDs_FD",         &Ds_PVwithDs_FD_vec.FD);
    tree->Branch("Ds_PVwithDs_FDxyerr",    &Ds_PVwithDs_FD_vec.FDxyerr);
    tree->Branch("Ds_PVwithDs_FDzerr",     &Ds_PVwithDs_FD_vec.FDzerr);
    tree->Branch("Ds_PVwithDs_FDerr",      &Ds_PVwithDs_FD_vec.FDerr);
    tree->Branch("Ds_PVwithDs_FDxychi2",   &Ds_PVwithDs_FD_vec.FDxychi2);
    tree->Branch("Ds_PVwithDs_FDzchi2",    &Ds_PVwithDs_FD_vec.FDzchi2);
    tree->Branch("Ds_PVwithDs_FDchi2",     &Ds_PVwithDs_FD_vec.FDchi2);
    tree->Branch("Ds_PVwithDs_dira",       &Ds_PVwithDs_FD_vec.dira);
    tree->Branch("Ds_PVwithDs_dira_angle", &Ds_PVwithDs_FD_vec.dira_angle);
    tree->Branch("Kp_PVwithDs_ip",         &Kp_PVwithDs_ip_vec.ip);
    tree->Branch("Kp_PVwithDs_iperr",      &Kp_PVwithDs_ip_vec.iperr);
    tree->Branch("Kp_PVwithDs_ipchi2",     &Kp_PVwithDs_ip_vec.ipchi2);
    tree->Branch("Km_PVwithDs_ip",         &Km_PVwithDs_ip_vec.ip);
    tree->Branch("Km_PVwithDs_iperr",      &Km_PVwithDs_ip_vec.iperr);
    tree->Branch("Km_PVwithDs_ipchi2",     &Km_PVwithDs_ip_vec.ipchi2);
    tree->Branch("pi_PVwithDs_ip",         &pi_PVwithDs_ip_vec.ip);
    tree->Branch("pi_PVwithDs_iperr",      &pi_PVwithDs_ip_vec.iperr);
    tree->Branch("pi_PVwithDs_ipchi2",     &pi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("phi_PVwithDs_ip",        &phi_PVwithDs_ip_vec.ip);
    tree->Branch("phi_PVwithDs_iperr",     &phi_PVwithDs_ip_vec.iperr);
    tree->Branch("phi_PVwithDs_ipchi2",    &phi_PVwithDs_ip_vec.ipchi2);
    tree->Branch("Ds_PVwithDs_ip",         &Ds_PVwithDs_ip_vec.ip);
    tree->Branch("Ds_PVwithDs_iperr",      &Ds_PVwithDs_ip_vec.iperr);
    tree->Branch("Ds_PVwithDs_ipchi2",     &Ds_PVwithDs_ip_vec.ipchi2);

    tree->Branch("Ds_IsoR03_sumChargedHadronPt", &Ds_IsoR03_vec.sumChargedHadronPt);
    tree->Branch("Ds_IsoR03_sumNeutralHadronEt", &Ds_IsoR03_vec.sumNeutralHadronEt);
    tree->Branch("Ds_IsoR03_sumPhotonEt",        &Ds_IsoR03_vec.sumPhotonEt);
    tree->Branch("Ds_IsoR03_sumPUPt",            &Ds_IsoR03_vec.sumPUPt);
    tree->Branch("Ds_IsoR03_PFIso",              &Ds_IsoR03_vec.PFIso);

    tree->Branch("Ds_IsoR04_sumChargedHadronPt", &Ds_IsoR04_vec.sumChargedHadronPt);
    tree->Branch("Ds_IsoR04_sumNeutralHadronEt", &Ds_IsoR04_vec.sumNeutralHadronEt);
    tree->Branch("Ds_IsoR04_sumPhotonEt",        &Ds_IsoR04_vec.sumPhotonEt);
    tree->Branch("Ds_IsoR04_sumPUPt",            &Ds_IsoR04_vec.sumPUPt);
    tree->Branch("Ds_IsoR04_PFIso",              &Ds_IsoR04_vec.PFIso);

    tree->Branch("Kp_match",        &Kp_match_vec);
    tree->Branch("Km_match",        &Km_match_vec);
    tree->Branch("pi_match",        &pi_match_vec);
    tree->Branch("match_entry",     &match_entry_vec);
    tree->Branch("non_match_entry", &non_match_entry_vec);

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

    tree->Branch("best_match_entry", &best_match_entry_vec);

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

    tree->Branch("best_mu_match", &best_mu_match_vec);
}

void PVTree::Gen_Reset()
{
    Gen_Kp.clearAll();
    Gen_Km.clearAll();
    Gen_pi.clearAll();
    Gen_phi.clearAll();
    Gen_Ds.clearAll();
    Gen_mu.clearAll();
    Gen_nu.clearAll();
    Gen_W.clearAll();
    Gen_H.clearAll();

    Gen_dR_Kp_Km  = null;
    Gen_dR_Kp_phi = null;
    Gen_dR_Km_phi = null;
    Gen_dR_Kp_pi  = null;
    Gen_dR_Km_pi  = null;
    Gen_dR_pi_phi = null;
    Gen_dR_Kp_Ds  = null;
    Gen_dR_Km_Ds  = null;
    Gen_dR_phi_Ds = null;
    Gen_dR_pi_Ds  = null;
    Gen_dR_Kp_mu  = null;
    Gen_dR_Km_mu  = null;
    Gen_dR_phi_mu = null;
    Gen_dR_pi_mu  = null;
    Gen_dR_Ds_mu  = null;

    Gen_Ds_dx   = null;
    Gen_Ds_dy   = null;
    Gen_Ds_dz   = null;
    Gen_Ds_FDxy = null;
    Gen_Ds_FD   = null;
}

void PVTree::Match_Reset()
{
    match_Kp.clearAll();
    match_Km.clearAll();
    match_pi.clearAll();
    match_phi.clearAll();
    match_Ds.clearAll();
    match_phiFit_Kp.clearAll();
    match_phiFit_Km.clearAll();
    match_phiFit_pi.clearAll();
    match_phiFit_phi.clearAll();
    match_phiFit_Ds.clearAll();
    match_DsFit_Kp.clearAll();
    match_DsFit_Km.clearAll();
    match_DsFit_pi.clearAll();
    match_DsFit_phi.clearAll();
    match_DsFit_Ds.clearAll();

    match_dR_Kp_Km  = null;
    match_dR_Kp_phi = null;
    match_dR_Km_phi = null;
    match_dR_Kp_pi  = null;
    match_dR_Km_pi  = null;
    match_dR_pi_phi = null;
    match_dR_Kp_Ds  = null;
    match_dR_Km_Ds  = null;
    match_dR_phi_Ds = null;
    match_dR_pi_Ds  = null;

    match_phiFit_dR_Kp_Km  = null;
    match_phiFit_dR_Kp_phi = null;
    match_phiFit_dR_Km_phi = null;
    match_phiFit_dR_Kp_pi  = null;
    match_phiFit_dR_Km_pi  = null;
    match_phiFit_dR_pi_phi = null;
    match_phiFit_dR_Kp_Ds  = null;
    match_phiFit_dR_Km_Ds  = null;
    match_phiFit_dR_phi_Ds = null;
    match_phiFit_dR_pi_Ds  = null;

    match_DsFit_dR_Kp_Km  = null;
    match_DsFit_dR_Kp_phi = null;
    match_DsFit_dR_Km_phi = null;
    match_DsFit_dR_Kp_pi  = null;
    match_DsFit_dR_Km_pi  = null;
    match_DsFit_dR_pi_phi = null;
    match_DsFit_dR_Kp_Ds  = null;
    match_DsFit_dR_Km_Ds  = null;
    match_DsFit_dR_phi_Ds = null;
    match_DsFit_dR_pi_Ds  = null;

    match_dxy_phi_Ds = null;
    match_dz_phi_Ds  = null;

    match_DsFit_Mconstraint_Ds_invm = null;

    match_phiFit.clearAll();
    match_DsFit.clearAll();
    match_PVnoDs.clearAll();
    match_PVwithDs.clearAll();

    match_Ds_primvtx_FD.clearAll();
    match_Kp_primvtx_ip.clearAll();
    match_Km_primvtx_ip.clearAll();
    match_pi_primvtx_ip.clearAll();
    match_phi_primvtx_ip.clearAll();
    match_Ds_primvtx_ip.clearAll();

    match_Ds_PVnoDs_FD.clearAll();
    match_Kp_PVnoDs_ip.clearAll();
    match_Km_PVnoDs_ip.clearAll();
    match_pi_PVnoDs_ip.clearAll();
    match_phi_PVnoDs_ip.clearAll();
    match_Ds_PVnoDs_ip.clearAll();

    match_Ds_PVwithDs_FD.clearAll();
    match_Kp_PVwithDs_ip.clearAll();
    match_Km_PVwithDs_ip.clearAll();
    match_pi_PVwithDs_ip.clearAll();
    match_phi_PVwithDs_ip.clearAll();
    match_Ds_PVwithDs_ip.clearAll();

    match_Ds_IsoR03.clearAll();
    match_Ds_IsoR04.clearAll();

    match_mu.clearAll();
    match_mu_IsoR03.clearAll();
    match_mu_IsoR04.clearAll();
    match_mu_primvtx_ip.clearAll();
}

void PVTree::FillMatchDs()
{
    match_Kp_vec.fillAll( match_Kp );
    match_Km_vec.fillAll( match_Km );
    match_pi_vec.fillAll( match_pi );
    match_phi_vec.fillAll( match_phi );
    match_Ds_vec.fillAll( match_Ds );
    match_phiFit_Kp_vec.fillAll( match_phiFit_Kp );
    match_phiFit_Km_vec.fillAll( match_phiFit_Km );
    match_phiFit_pi_vec.fillAll( match_phiFit_pi );
    match_phiFit_phi_vec.fillAll( match_phiFit_phi );
    match_phiFit_Ds_vec.fillAll( match_phiFit_Ds );
    match_DsFit_Kp_vec.fillAll( match_DsFit_Kp );
    match_DsFit_Km_vec.fillAll( match_DsFit_Km );
    match_DsFit_pi_vec.fillAll( match_DsFit_pi );
    match_DsFit_phi_vec.fillAll( match_DsFit_phi );
    match_DsFit_Ds_vec.fillAll( match_DsFit_Ds );
    
    match_dR_Kp_Km_vec.push_back( match_dR_Kp_Km );
    match_dR_Kp_phi_vec.push_back( match_dR_Kp_phi );
    match_dR_Km_phi_vec.push_back( match_dR_Km_phi );
    match_dR_Kp_pi_vec.push_back( match_dR_Kp_pi );
    match_dR_Km_pi_vec.push_back( match_dR_Km_pi );
    match_dR_pi_phi_vec.push_back( match_dR_pi_phi );
    match_dR_Kp_Ds_vec.push_back( match_dR_Kp_Ds );
    match_dR_Km_Ds_vec.push_back( match_dR_Km_Ds );
    match_dR_phi_Ds_vec.push_back( match_dR_phi_Ds );
    match_dR_pi_Ds_vec.push_back( match_dR_pi_Ds );

    match_phiFit_dR_Kp_Km_vec.push_back( match_phiFit_dR_Kp_Km );
    match_phiFit_dR_Kp_phi_vec.push_back( match_phiFit_dR_Kp_phi );
    match_phiFit_dR_Km_phi_vec.push_back( match_phiFit_dR_Km_phi );
    match_phiFit_dR_Kp_pi_vec.push_back( match_phiFit_dR_Kp_pi );
    match_phiFit_dR_Km_pi_vec.push_back( match_phiFit_dR_Km_pi );
    match_phiFit_dR_pi_phi_vec.push_back( match_phiFit_dR_pi_phi );
    match_phiFit_dR_Kp_Ds_vec.push_back( match_phiFit_dR_Kp_Ds );
    match_phiFit_dR_Km_Ds_vec.push_back( match_phiFit_dR_Km_Ds );
    match_phiFit_dR_phi_Ds_vec.push_back( match_phiFit_dR_phi_Ds );
    match_phiFit_dR_pi_Ds_vec.push_back( match_phiFit_dR_pi_Ds );

    match_DsFit_dR_Kp_Km_vec.push_back( match_DsFit_dR_Kp_Km );
    match_DsFit_dR_Kp_phi_vec.push_back( match_DsFit_dR_Kp_phi );
    match_DsFit_dR_Km_phi_vec.push_back( match_DsFit_dR_Km_phi );
    match_DsFit_dR_Kp_pi_vec.push_back( match_DsFit_dR_Kp_pi );
    match_DsFit_dR_Km_pi_vec.push_back( match_DsFit_dR_Km_pi );
    match_DsFit_dR_pi_phi_vec.push_back( match_DsFit_dR_pi_phi );
    match_DsFit_dR_Kp_Ds_vec.push_back( match_DsFit_dR_Kp_Ds );
    match_DsFit_dR_Km_Ds_vec.push_back( match_DsFit_dR_Km_Ds );
    match_DsFit_dR_phi_Ds_vec.push_back( match_DsFit_dR_phi_Ds );
    match_DsFit_dR_pi_Ds_vec.push_back( match_DsFit_dR_pi_Ds );

    match_dxy_phi_Ds_vec.push_back( match_dxy_phi_Ds );
    match_dz_phi_Ds_vec.push_back( match_dz_phi_Ds );

    match_DsFit_Mconstraint_Ds_invm_vec.push_back( match_DsFit_Mconstraint_Ds_invm );

    match_phiFit_vec.fillAll( match_phiFit );
    match_DsFit_vec.fillAll( match_DsFit );
    match_PVnoDs_vec.fillAll( match_PVnoDs );
    match_PVwithDs_vec.fillAll( match_PVwithDs );

    match_Ds_primvtx_FD_vec.fillAll( match_Ds_primvtx_FD );
    match_Kp_primvtx_ip_vec.fillAll( match_Kp_primvtx_ip );
    match_Km_primvtx_ip_vec.fillAll( match_Km_primvtx_ip );
    match_pi_primvtx_ip_vec.fillAll( match_pi_primvtx_ip );
    match_phi_primvtx_ip_vec.fillAll( match_phi_primvtx_ip );
    match_Ds_primvtx_ip_vec.fillAll( match_Ds_primvtx_ip );

    match_Ds_PVnoDs_FD_vec.fillAll( match_Ds_PVnoDs_FD );
    match_Kp_PVnoDs_ip_vec.fillAll( match_Kp_PVnoDs_ip );
    match_Km_PVnoDs_ip_vec.fillAll( match_Km_PVnoDs_ip );
    match_pi_PVnoDs_ip_vec.fillAll( match_pi_PVnoDs_ip );
    match_phi_PVnoDs_ip_vec.fillAll( match_phi_PVnoDs_ip );
    match_Ds_PVnoDs_ip_vec.fillAll( match_Ds_PVnoDs_ip );

    match_Ds_PVwithDs_FD_vec.fillAll( match_Ds_PVwithDs_FD );
    match_Kp_PVwithDs_ip_vec.fillAll( match_Kp_PVwithDs_ip );
    match_Km_PVwithDs_ip_vec.fillAll( match_Km_PVwithDs_ip );
    match_pi_PVwithDs_ip_vec.fillAll( match_pi_PVwithDs_ip );
    match_phi_PVwithDs_ip_vec.fillAll( match_phi_PVwithDs_ip );
    match_Ds_PVwithDs_ip_vec.fillAll( match_Ds_PVwithDs_ip );

    match_Ds_IsoR03_vec.fillAll( match_Ds_IsoR03 );
    match_Ds_IsoR04_vec.fillAll( match_Ds_IsoR04 );
}

void PVTree::FillMatchmu(){
    match_mu_vec.fillAll( match_mu );
    match_mu_IsoR03_vec.fillAll( match_mu_IsoR03 );
    match_mu_IsoR04_vec.fillAll( match_mu_IsoR04 );
    match_mu_primvtx_ip_vec.fillAll( match_mu_primvtx_ip );
}

void PVTree::FillDs()
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

    Kp_match_vec.push_back(Kp_match);
    Km_match_vec.push_back(Km_match);
    pi_match_vec.push_back(pi_match);
    match_entry_vec.push_back(match_entry);
    non_match_entry_vec.push_back(non_match_entry);
}

void PVTree::FillBestDs(int idxmax)
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

    best_match_entry_vec.push_back( match_entry_vec[idxmax] );
}

void PVTree::FillBestmu()
{
    mu_vec.fillAll( mu );
    mu_IsoR03_vec.fillAll( mu_IsoR03 );
    mu_IsoR04_vec.fillAll( mu_IsoR04 );
    mu_primvtx_ip_vec.fillAll( mu_primvtx_ip );
    best_mu_match_vec.push_back(best_mu_match);
}

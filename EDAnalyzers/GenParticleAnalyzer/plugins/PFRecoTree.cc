#include "EDAnalyzers/GenParticleAnalyzer/interface/PFRecoTree.h"

PFRecoTree::PFRecoTree(TTree *tree_)
{
    tree = tree_;
}

void PFRecoTree::Init()
{
    num_reco_phi = 0;
    num_reco_Ds = 0;

    Kp_ETA.clear();
    Kp_PHI.clear();
    Kp_ORIVX_X.clear();
    Kp_ORIVX_Y.clear();
    Kp_ORIVX_Z.clear();
    Kp_P.clear();
    Kp_PT.clear();
    Kp_PX.clear();
    Kp_PY.clear();
    Kp_PZ.clear();
    Kp_PP.clear();
    Kp_PL.clear();

    Km_ETA.clear();
    Km_PHI.clear();
    Km_ORIVX_X.clear();
    Km_ORIVX_Y.clear();
    Km_ORIVX_Z.clear();
    Km_P.clear();
    Km_PT.clear();
    Km_PX.clear();
    Km_PY.clear();
    Km_PZ.clear();
    Km_PP.clear();
    Km_PL.clear();

    phi_ETA.clear();
    phi_PHI.clear();
    phi_P.clear();
    phi_PT.clear();
    phi_PX.clear();
    phi_PY.clear();
    phi_PZ.clear();
    phi_PP.clear();
    phi_PL.clear();
    phi_CHI2.clear();
    phi_M.clear();
    phi_ENDVX_X.clear();
    phi_ENDVX_Y.clear();
    phi_ENDVX_Z.clear();

    pi_ETA.clear();
    pi_PHI.clear();
    pi_ORIVX_X.clear();
    pi_ORIVX_Y.clear();
    pi_ORIVX_Z.clear();
    pi_P.clear();
    pi_PT.clear();
    pi_PX.clear();
    pi_PY.clear();
    pi_PZ.clear();
    pi_PP.clear();
    pi_PL.clear();

    Ds_ETA.clear();
    Ds_PHI.clear();
    Ds_P.clear();
    Ds_PT.clear();
    Ds_PX.clear();
    Ds_PY.clear();
    Ds_PZ.clear();
    Ds_CHI2.clear();
    Ds_M.clear();
    Ds_ENDVX_X.clear();
    Ds_ENDVX_Y.clear();
    Ds_ENDVX_Z.clear();

    dR_Kp_Km.clear();
    dR_Kp_phi.clear();
    dR_Km_phi.clear();
    dR_Kp_pi.clear();
    dR_Km_pi.clear();
    dR_phi_pi.clear();
    dR_Kp_Ds.clear();
    dR_Km_Ds.clear();
    dR_phi_Ds.clear();
    dR_pi_Ds.clear();

    d_Kp_Km.clear();
    d_Kp_phi.clear();
    d_Km_phi.clear();
    d_Kp_pi.clear();
    d_Km_pi.clear();
    d_phi_pi.clear();
    d_Kp_Ds.clear();
    d_Km_Ds.clear();
    d_phi_Ds.clear();
    d_pi_Ds.clear();

    phiFit_Kp_ETA.clear();
    phiFit_Kp_PHI.clear();
    phiFit_Kp_P.clear();
    phiFit_Kp_PT.clear();
    phiFit_Kp_PX.clear();
    phiFit_Kp_PY.clear();
    phiFit_Kp_PZ.clear();
    phiFit_Kp_PP.clear();
    phiFit_Kp_PL.clear();

    phiFit_Km_ETA.clear();
    phiFit_Km_PHI.clear();
    phiFit_Km_P.clear();
    phiFit_Km_PT.clear();
    phiFit_Km_PX.clear();
    phiFit_Km_PY.clear();
    phiFit_Km_PZ.clear();
    phiFit_Km_PP.clear();
    phiFit_Km_PL.clear();

    DsFit_Kp_ETA.clear();
    DsFit_Kp_PHI.clear();
    DsFit_Kp_P.clear();
    DsFit_Kp_PT.clear();
    DsFit_Kp_PX.clear();
    DsFit_Kp_PY.clear();
    DsFit_Kp_PZ.clear();
    DsFit_Kp_PP.clear();
    DsFit_Kp_PL.clear();

    DsFit_Km_ETA.clear();
    DsFit_Km_PHI.clear();
    DsFit_Km_P.clear();
    DsFit_Km_PT.clear();
    DsFit_Km_PX.clear();
    DsFit_Km_PY.clear();
    DsFit_Km_PZ.clear();
    DsFit_Km_PP.clear();
    DsFit_Km_PL.clear();

    DsFit_phi_ETA.clear();
    DsFit_phi_PHI.clear();
    DsFit_phi_P.clear();
    DsFit_phi_PT.clear();
    DsFit_phi_PX.clear();
    DsFit_phi_PY.clear();
    DsFit_phi_PZ.clear();
    DsFit_phi_PP.clear();
    DsFit_phi_PL.clear();
    DsFit_phi_M.clear();

    DsFit_pi_ETA.clear();
    DsFit_pi_PHI.clear();
    DsFit_pi_P.clear();
    DsFit_pi_PT.clear();
    DsFit_pi_PX.clear();
    DsFit_pi_PY.clear();
    DsFit_pi_PZ.clear();
    DsFit_pi_PP.clear();
    DsFit_pi_PL.clear();

    DsFit_Mconstraint_Ds_M.clear();
}

void PFRecoTree::CreateBranches()
{
    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

    tree->Branch("Kp_ETA", &Kp_ETA);
    tree->Branch("Kp_PHI", &Kp_PHI);
    tree->Branch("Kp_ORIVX_X", &Kp_ORIVX_X);
    tree->Branch("Kp_ORIVX_Y", &Kp_ORIVX_Y);
    tree->Branch("Kp_ORIVX_Z", &Kp_ORIVX_Z);
    tree->Branch("Kp_P", &Kp_P);
    tree->Branch("Kp_PT", &Kp_PT);
    tree->Branch("Kp_PX", &Kp_PX);
    tree->Branch("Kp_PY", &Kp_PY);
    tree->Branch("Kp_PZ", &Kp_PZ);
    tree->Branch("Kp_PP", &Kp_PP);
    tree->Branch("Kp_PL", &Kp_PL);

    tree->Branch("Km_ETA", &Km_ETA);
    tree->Branch("Km_PHI", &Km_PHI);
    tree->Branch("Km_ORIVX_X", &Km_ORIVX_X);
    tree->Branch("Km_ORIVX_Y", &Km_ORIVX_Y);
    tree->Branch("Km_ORIVX_Z", &Km_ORIVX_Z);
    tree->Branch("Km_P", &Km_P);
    tree->Branch("Km_PT", &Km_PT);
    tree->Branch("Km_PX", &Km_PX);
    tree->Branch("Km_PY", &Km_PY);
    tree->Branch("Km_PZ", &Km_PZ);
    tree->Branch("Km_PP", &Km_PP);
    tree->Branch("Km_PL", &Km_PL);

    tree->Branch("phi_ETA", &phi_ETA);
    tree->Branch("phi_PHI", &phi_PHI);
    tree->Branch("phi_P", &phi_P);
    tree->Branch("phi_PT", &phi_PT);
    tree->Branch("phi_PX", &phi_PX);
    tree->Branch("phi_PY", &phi_PY);
    tree->Branch("phi_PZ", &phi_PZ);
    tree->Branch("phi_PP", &phi_PP);
    tree->Branch("phi_PL", &phi_PL);
    tree->Branch("phi_CHI2", &phi_CHI2);
    tree->Branch("phi_M", &phi_M);
    tree->Branch("phi_ENDVX_X", &phi_ENDVX_X);
    tree->Branch("phi_ENDVX_Y", &phi_ENDVX_Y);
    tree->Branch("phi_ENDVX_Z", &phi_ENDVX_Z);

    tree->Branch("pi_ETA", &pi_ETA);
    tree->Branch("pi_PHI", &pi_PHI);
    tree->Branch("pi_ORIVX_X", &pi_ORIVX_X);
    tree->Branch("pi_ORIVX_Y", &pi_ORIVX_Y);
    tree->Branch("pi_ORIVX_Z", &pi_ORIVX_Z);
    tree->Branch("pi_P", &pi_P);
    tree->Branch("pi_PT", &pi_PT);
    tree->Branch("pi_PX", &pi_PX);
    tree->Branch("pi_PY", &pi_PY);
    tree->Branch("pi_PZ", &pi_PZ);
    tree->Branch("pi_PP", &pi_PP);
    tree->Branch("pi_PL", &pi_PL);

    tree->Branch("Ds_ETA", &Ds_ETA);
    tree->Branch("Ds_PHI", &Ds_PHI);
    tree->Branch("Ds_P", &Ds_P);
    tree->Branch("Ds_PT", &Ds_PT);
    tree->Branch("Ds_PX", &Ds_PX);
    tree->Branch("Ds_PY", &Ds_PY);
    tree->Branch("Ds_PZ", &Ds_PZ);
    tree->Branch("Ds_CHI2", &Ds_CHI2);
    tree->Branch("Ds_M", &Ds_M);
    tree->Branch("Ds_ENDVX_X", &Ds_ENDVX_X);
    tree->Branch("Ds_ENDVX_Y", &Ds_ENDVX_Y);
    tree->Branch("Ds_ENDVX_Z", &Ds_ENDVX_Z);

    tree->Branch("dR_Kp_Km", &dR_Kp_Km);
    tree->Branch("dR_Kp_phi", &dR_Kp_phi);
    tree->Branch("dR_Km_phi", &dR_Km_phi);
    tree->Branch("dR_Kp_pi", &dR_Kp_pi);
    tree->Branch("dR_Km_pi", &dR_Km_pi);
    tree->Branch("dR_phi_pi", &dR_phi_pi);
    tree->Branch("dR_Kp_Ds", &dR_Kp_Ds);
    tree->Branch("dR_Km_Ds", &dR_Km_Ds);
    tree->Branch("dR_phi_Ds", &dR_phi_Ds);
    tree->Branch("dR_pi_Ds", &dR_pi_Ds);

    tree->Branch("d_Kp_Km", &d_Kp_Km);
    tree->Branch("d_Kp_phi", &d_Kp_phi);
    tree->Branch("d_Km_phi", &d_Km_phi);
    tree->Branch("d_Kp_pi", &d_Kp_pi);
    tree->Branch("d_Km_pi", &d_Km_pi);
    tree->Branch("d_phi_pi", &d_phi_pi);
    tree->Branch("d_Kp_Ds", &d_Kp_Ds);
    tree->Branch("d_Km_Ds", &d_Km_Ds);
    tree->Branch("d_phi_Ds", &d_phi_Ds);
    tree->Branch("d_pi_Ds", &d_pi_Ds);

    tree->Branch("phiFit_Kp_ETA", &phiFit_Kp_ETA);
    tree->Branch("phiFit_Kp_PHI", &phiFit_Kp_PHI);
    tree->Branch("phiFit_Kp_P", &phiFit_Kp_P);
    tree->Branch("phiFit_Kp_PT", &phiFit_Kp_PT);
    tree->Branch("phiFit_Kp_PX", &phiFit_Kp_PX);
    tree->Branch("phiFit_Kp_PY", &phiFit_Kp_PY);
    tree->Branch("phiFit_Kp_PZ", &phiFit_Kp_PZ);
    tree->Branch("phiFit_Kp_PP", &phiFit_Kp_PP);
    tree->Branch("phiFit_Kp_PL", &phiFit_Kp_PL);

    tree->Branch("phiFit_Km_ETA", &phiFit_Km_ETA);
    tree->Branch("phiFit_Km_PHI", &phiFit_Km_PHI);
    tree->Branch("phiFit_Km_P", &phiFit_Km_P);
    tree->Branch("phiFit_Km_PT", &phiFit_Km_PT);
    tree->Branch("phiFit_Km_PX", &phiFit_Km_PX);
    tree->Branch("phiFit_Km_PY", &phiFit_Km_PY);
    tree->Branch("phiFit_Km_PZ", &phiFit_Km_PZ);
    tree->Branch("phiFit_Km_PP", &phiFit_Km_PP);
    tree->Branch("phiFit_Km_PL", &phiFit_Km_PL);

    tree->Branch("DsFit_Kp_ETA", &DsFit_Kp_ETA);
    tree->Branch("DsFit_Kp_PHI", &DsFit_Kp_PHI);
    tree->Branch("DsFit_Kp_P", &DsFit_Kp_P);
    tree->Branch("DsFit_Kp_PT", &DsFit_Kp_PT);
    tree->Branch("DsFit_Kp_PX", &DsFit_Kp_PX);
    tree->Branch("DsFit_Kp_PY", &DsFit_Kp_PY);
    tree->Branch("DsFit_Kp_PZ", &DsFit_Kp_PZ);
    tree->Branch("DsFit_Kp_PP", &DsFit_Kp_PP);
    tree->Branch("DsFit_Kp_PL", &DsFit_Kp_PL);

    tree->Branch("DsFit_Km_ETA", &DsFit_Km_ETA);
    tree->Branch("DsFit_Km_PHI", &DsFit_Km_PHI);
    tree->Branch("DsFit_Km_P", &DsFit_Km_P);
    tree->Branch("DsFit_Km_PT", &DsFit_Km_PT);
    tree->Branch("DsFit_Km_PX", &DsFit_Km_PX);
    tree->Branch("DsFit_Km_PY", &DsFit_Km_PY);
    tree->Branch("DsFit_Km_PZ", &DsFit_Km_PZ);
    tree->Branch("DsFit_Km_PP", &DsFit_Km_PP);
    tree->Branch("DsFit_Km_PL", &DsFit_Km_PL);

    tree->Branch("DsFit_phi_ETA", &DsFit_phi_ETA);
    tree->Branch("DsFit_phi_PHI", &DsFit_phi_PHI);
    tree->Branch("DsFit_phi_P", &DsFit_phi_P);
    tree->Branch("DsFit_phi_PT", &DsFit_phi_PT);
    tree->Branch("DsFit_phi_PX", &DsFit_phi_PX);
    tree->Branch("DsFit_phi_PY", &DsFit_phi_PY);
    tree->Branch("DsFit_phi_PZ", &DsFit_phi_PZ);
    tree->Branch("DsFit_phi_PP", &DsFit_phi_PP);
    tree->Branch("DsFit_phi_PL", &DsFit_phi_PL);
    tree->Branch("DsFit_phi_M", &DsFit_phi_M);

    tree->Branch("DsFit_pi_ETA", &DsFit_pi_ETA);
    tree->Branch("DsFit_pi_PHI", &DsFit_pi_PHI);
    tree->Branch("DsFit_pi_P", &DsFit_pi_P);
    tree->Branch("DsFit_pi_PT", &DsFit_pi_PT);
    tree->Branch("DsFit_pi_PX", &DsFit_pi_PX);
    tree->Branch("DsFit_pi_PY", &DsFit_pi_PY);
    tree->Branch("DsFit_pi_PZ", &DsFit_pi_PZ);
    tree->Branch("DsFit_pi_PP", &DsFit_pi_PP);
    tree->Branch("DsFit_pi_PL", &DsFit_pi_PL);

    tree->Branch("DsFit_Mconstraint_Ds_M", &DsFit_Mconstraint_Ds_M);
};

#include "EDAnalyzers/GenParticleAnalyzer/interface/GenTree.h"

GenTree::GenTree(TTree *tree_)
{
    tree = tree_;
}

void GenTree::Init()
{
    num_Gen_Kp = 0;
    num_Gen_Km = 0;
    num_Gen_pi = 0;
    num_Gen_phi = 0;
    num_Gen_Ds = 0;

    Gen_Kp_ETA.clear();
    Gen_Kp_PHI.clear();
    Gen_Kp_ORIVX_X.clear();
    Gen_Kp_ORIVX_Y.clear();
    Gen_Kp_ORIVX_Z.clear();
    Gen_Kp_P.clear();
    Gen_Kp_PT.clear();
    Gen_Kp_PX.clear();
    Gen_Kp_PY.clear();
    Gen_Kp_PZ.clear();
    Gen_Kp_PP.clear();
    Gen_Kp_PL.clear();

    Gen_Km_ETA.clear();
    Gen_Km_PHI.clear();
    Gen_Km_ORIVX_X.clear();
    Gen_Km_ORIVX_Y.clear();
    Gen_Km_ORIVX_Z.clear();
    Gen_Km_P.clear();
    Gen_Km_PT.clear();
    Gen_Km_PX.clear();
    Gen_Km_PY.clear();
    Gen_Km_PZ.clear();
    Gen_Km_PP.clear();
    Gen_Km_PL.clear();

    Gen_pi_ETA.clear();
    Gen_pi_PHI.clear();
    Gen_pi_ORIVX_X.clear();
    Gen_pi_ORIVX_Y.clear();
    Gen_pi_ORIVX_Z.clear();
    Gen_pi_P.clear();
    Gen_pi_PT.clear();
    Gen_pi_PX.clear();
    Gen_pi_PY.clear();
    Gen_pi_PZ.clear();
    Gen_pi_PP.clear();
    Gen_pi_PL.clear();

    Gen_phi_ETA.clear();
    Gen_phi_PHI.clear();
    Gen_phi_ORIVX_X.clear();
    Gen_phi_ORIVX_Y.clear();
    Gen_phi_ORIVX_Z.clear();
    Gen_phi_P.clear();
    Gen_phi_PT.clear();
    Gen_phi_PX.clear();
    Gen_phi_PY.clear();
    Gen_phi_PZ.clear();
    Gen_phi_PP.clear();
    Gen_phi_PL.clear();

    Gen_Ds_ETA.clear();
    Gen_Ds_PHI.clear();
    Gen_Ds_ORIVX_X.clear();
    Gen_Ds_ORIVX_Y.clear();
    Gen_Ds_ORIVX_Z.clear();
    Gen_Ds_P.clear();
    Gen_Ds_PT.clear();
    Gen_Ds_PX.clear();
    Gen_Ds_PY.clear();
    Gen_Ds_PZ.clear();

    Gen_dR_Kp_Km.clear();
    Gen_dR_Kp_phi.clear();
    Gen_dR_Km_phi.clear();
    Gen_dR_Kp_pi.clear();
    Gen_dR_Km_pi.clear();
    Gen_dR_phi_pi.clear();
    Gen_dR_Kp_Ds.clear();
    Gen_dR_Km_Ds.clear();
    Gen_dR_phi_Ds.clear();
    Gen_dR_pi_Ds.clear();

    Gen_d_Kp_Km.clear();
    Gen_d_Kp_pi.clear();
    Gen_d_Km_pi.clear();

    // Matched PF Candidates
    num_match_Kp = 0;
    num_match_Km = 0;
    num_match_pi = 0;

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

void GenTree::CreateBranches()
{
    // Gen Particles
    tree->Branch("num_Gen_Kp", &num_Gen_Kp);
    tree->Branch("num_Gen_Km", &num_Gen_Km);
    tree->Branch("num_Gen_pi", &num_Gen_pi);
    tree->Branch("num_Gen_phi", &num_Gen_phi);
    tree->Branch("num_Gen_Ds", &num_Gen_Ds);

    tree->Branch("Gen_Kp_ETA", &Gen_Kp_ETA);
    tree->Branch("Gen_Kp_PHI", &Gen_Kp_PHI);
    tree->Branch("Gen_Kp_ORIVX_X", &Gen_Kp_ORIVX_X);
    tree->Branch("Gen_Kp_ORIVX_Y", &Gen_Kp_ORIVX_Y);
    tree->Branch("Gen_Kp_ORIVX_Z", &Gen_Kp_ORIVX_Z);
    tree->Branch("Gen_Kp_P", &Gen_Kp_P);
    tree->Branch("Gen_Kp_PT", &Gen_Kp_PT);
    tree->Branch("Gen_Kp_PX", &Gen_Kp_PX);
    tree->Branch("Gen_Kp_PY", &Gen_Kp_PY);
    tree->Branch("Gen_Kp_PZ", &Gen_Kp_PZ);
    tree->Branch("Gen_Kp_PP", &Gen_Kp_PP);
    tree->Branch("Gen_Kp_PL", &Gen_Kp_PL);

    tree->Branch("Gen_Km_ETA", &Gen_Km_ETA);
    tree->Branch("Gen_Km_PHI", &Gen_Km_PHI);
    tree->Branch("Gen_Km_ORIVX_X", &Gen_Km_ORIVX_X);
    tree->Branch("Gen_Km_ORIVX_Y", &Gen_Km_ORIVX_Y);
    tree->Branch("Gen_Km_ORIVX_Z", &Gen_Km_ORIVX_Z);
    tree->Branch("Gen_Km_P", &Gen_Km_P);
    tree->Branch("Gen_Km_PT", &Gen_Km_PT);
    tree->Branch("Gen_Km_PX", &Gen_Km_PX);
    tree->Branch("Gen_Km_PY", &Gen_Km_PY);
    tree->Branch("Gen_Km_PZ", &Gen_Km_PZ);
    tree->Branch("Gen_Km_PP", &Gen_Km_PP);
    tree->Branch("Gen_Km_PL", &Gen_Km_PL);

    tree->Branch("Gen_pi_ETA", &Gen_pi_ETA);
    tree->Branch("Gen_pi_PHI", &Gen_pi_PHI);
    tree->Branch("Gen_pi_ORIVX_X", &Gen_pi_ORIVX_X);
    tree->Branch("Gen_pi_ORIVX_Y", &Gen_pi_ORIVX_Y);
    tree->Branch("Gen_pi_ORIVX_Z", &Gen_pi_ORIVX_Z);
    tree->Branch("Gen_pi_P", &Gen_pi_P);
    tree->Branch("Gen_pi_PT", &Gen_pi_PT);
    tree->Branch("Gen_pi_PX", &Gen_pi_PX);
    tree->Branch("Gen_pi_PY", &Gen_pi_PY);
    tree->Branch("Gen_pi_PZ", &Gen_pi_PZ);
    tree->Branch("Gen_pi_PP", &Gen_pi_PP);
    tree->Branch("Gen_pi_PL", &Gen_pi_PL);

    tree->Branch("Gen_phi_ETA", &Gen_phi_ETA);
    tree->Branch("Gen_phi_PHI", &Gen_phi_PHI);
    tree->Branch("Gen_phi_ORIVX_X", &Gen_phi_ORIVX_X);
    tree->Branch("Gen_phi_ORIVX_Y", &Gen_phi_ORIVX_Y);
    tree->Branch("Gen_phi_ORIVX_Z", &Gen_phi_ORIVX_Z);
    tree->Branch("Gen_phi_P", &Gen_phi_P);
    tree->Branch("Gen_phi_PT", &Gen_phi_PT);
    tree->Branch("Gen_phi_PX", &Gen_phi_PX);
    tree->Branch("Gen_phi_PY", &Gen_phi_PY);
    tree->Branch("Gen_phi_PZ", &Gen_phi_PZ);
    tree->Branch("Gen_phi_PP", &Gen_phi_PP);
    tree->Branch("Gen_phi_PL", &Gen_phi_PL);

    tree->Branch("Gen_Ds_ETA", &Gen_Ds_ETA);
    tree->Branch("Gen_Ds_PHI", &Gen_Ds_PHI);
    tree->Branch("Gen_Ds_ORIVX_X", &Gen_Ds_ORIVX_X);
    tree->Branch("Gen_Ds_ORIVX_Y", &Gen_Ds_ORIVX_Y);
    tree->Branch("Gen_Ds_ORIVX_Z", &Gen_Ds_ORIVX_Z);
    tree->Branch("Gen_Ds_P", &Gen_Ds_P);
    tree->Branch("Gen_Ds_PT", &Gen_Ds_PT);
    tree->Branch("Gen_Ds_PX", &Gen_Ds_PX);
    tree->Branch("Gen_Ds_PY", &Gen_Ds_PY);
    tree->Branch("Gen_Ds_PZ", &Gen_Ds_PZ);

    tree->Branch("Gen_dR_Kp_Km", &Gen_dR_Kp_Km);
    tree->Branch("Gen_dR_Kp_phi", &Gen_dR_Kp_phi);
    tree->Branch("Gen_dR_Km_phi", &Gen_dR_Km_phi);
    tree->Branch("Gen_dR_Kp_pi", &Gen_dR_Kp_pi);
    tree->Branch("Gen_dR_Km_pi", &Gen_dR_Km_pi);
    tree->Branch("Gen_dR_phi_pi", &Gen_dR_phi_pi);
    tree->Branch("Gen_dR_Kp_Ds", &Gen_dR_Kp_Ds);
    tree->Branch("Gen_dR_Km_Ds", &Gen_dR_Km_Ds);
    tree->Branch("Gen_dR_phi_Ds", &Gen_dR_phi_Ds);
    tree->Branch("Gen_dR_pi_Ds", &Gen_dR_pi_Ds);

    tree->Branch("Gen_d_Kp_Km", &Gen_d_Kp_Km);
    tree->Branch("Gen_d_Kp_pi", &Gen_d_Kp_pi);
    tree->Branch("Gen_d_Km_pi", &Gen_d_Km_pi);

    // Matched PF Candidates
    tree->Branch("num_match_Kp", &num_match_Kp);
    tree->Branch("num_match_Km", &num_match_Km);
    tree->Branch("num_match_pi", &num_match_pi);

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

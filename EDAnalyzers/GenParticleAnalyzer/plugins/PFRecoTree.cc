#include "EDAnalyzers/GenParticleAnalyzer/interface/PFRecoTree.h"

PFRecoTree::PFRecoTree(TTree *tree_)
{
    tree = tree_;
}

void PFRecoTree::Init()
{
    num_reco_phi = 0;
    num_reco_Ds = 0;

    Kp_ETA_vec.clear();
    Kp_PHI_vec.clear();
    Kp_ORIVX_X_vec.clear();
    Kp_ORIVX_Y_vec.clear();
    Kp_ORIVX_Z_vec.clear();
    Kp_P_vec.clear();
    Kp_PT_vec.clear();
    Kp_PX_vec.clear();
    Kp_PY_vec.clear();
    Kp_PZ_vec.clear();
    Kp_PP_vec.clear();
    Kp_PL_vec.clear();

    Km_ETA_vec.clear();
    Km_PHI_vec.clear();
    Km_ORIVX_X_vec.clear();
    Km_ORIVX_Y_vec.clear();
    Km_ORIVX_Z_vec.clear();
    Km_P_vec.clear();
    Km_PT_vec.clear();
    Km_PX_vec.clear();
    Km_PY_vec.clear();
    Km_PZ_vec.clear();
    Km_PP_vec.clear();
    Km_PL_vec.clear();

    phi_ETA_vec.clear();
    phi_PHI_vec.clear();
    phi_P_vec.clear();
    phi_PT_vec.clear();
    phi_PX_vec.clear();
    phi_PY_vec.clear();
    phi_PZ_vec.clear();
    phi_PP_vec.clear();
    phi_PL_vec.clear();
    phi_CHI2_vec.clear();
    phi_M_vec.clear();
    phi_ENDVX_X_vec.clear();
    phi_ENDVX_Y_vec.clear();
    phi_ENDVX_Z_vec.clear();

    pi_ETA_vec.clear();
    pi_PHI_vec.clear();
    pi_ORIVX_X_vec.clear();
    pi_ORIVX_Y_vec.clear();
    pi_ORIVX_Z_vec.clear();
    pi_P_vec.clear();
    pi_PT_vec.clear();
    pi_PX_vec.clear();
    pi_PY_vec.clear();
    pi_PZ_vec.clear();
    pi_PP_vec.clear();
    pi_PL_vec.clear();

    Ds_ETA_vec.clear();
    Ds_PHI_vec.clear();
    Ds_P_vec.clear();
    Ds_PT_vec.clear();
    Ds_PX_vec.clear();
    Ds_PY_vec.clear();
    Ds_PZ_vec.clear();
    Ds_CHI2_vec.clear();
    Ds_M_vec.clear();
    Ds_ENDVX_X_vec.clear();
    Ds_ENDVX_Y_vec.clear();
    Ds_ENDVX_Z_vec.clear();

    dR_Kp_Km_vec.clear();
    dR_Kp_phi_vec.clear();
    dR_Km_phi_vec.clear();
    dR_Kp_pi_vec.clear();
    dR_Km_pi_vec.clear();
    dR_phi_pi_vec.clear();
    dR_Kp_Ds_vec.clear();
    dR_Km_Ds_vec.clear();
    dR_phi_Ds_vec.clear();
    dR_pi_Ds_vec.clear();

    d_Kp_Km_vec.clear();
    d_Kp_phi_vec.clear();
    d_Km_phi_vec.clear();
    d_Kp_pi_vec.clear();
    d_Km_pi_vec.clear();
    d_phi_pi_vec.clear();
    d_Kp_Ds_vec.clear();
    d_Km_Ds_vec.clear();
    d_phi_Ds_vec.clear();
    d_pi_Ds_vec.clear();

    phiFit_Kp_ETA_vec.clear();
    phiFit_Kp_PHI_vec.clear();
    phiFit_Kp_P_vec.clear();
    phiFit_Kp_PT_vec.clear();
    phiFit_Kp_PX_vec.clear();
    phiFit_Kp_PY_vec.clear();
    phiFit_Kp_PZ_vec.clear();
    phiFit_Kp_PP_vec.clear();
    phiFit_Kp_PL_vec.clear();

    phiFit_Km_ETA_vec.clear();
    phiFit_Km_PHI_vec.clear();
    phiFit_Km_P_vec.clear();
    phiFit_Km_PT_vec.clear();
    phiFit_Km_PX_vec.clear();
    phiFit_Km_PY_vec.clear();
    phiFit_Km_PZ_vec.clear();
    phiFit_Km_PP_vec.clear();
    phiFit_Km_PL_vec.clear();

    DsFit_Kp_ETA_vec.clear();
    DsFit_Kp_PHI_vec.clear();
    DsFit_Kp_P_vec.clear();
    DsFit_Kp_PT_vec.clear();
    DsFit_Kp_PX_vec.clear();
    DsFit_Kp_PY_vec.clear();
    DsFit_Kp_PZ_vec.clear();
    DsFit_Kp_PP_vec.clear();
    DsFit_Kp_PL_vec.clear();

    DsFit_Km_ETA_vec.clear();
    DsFit_Km_PHI_vec.clear();
    DsFit_Km_P_vec.clear();
    DsFit_Km_PT_vec.clear();
    DsFit_Km_PX_vec.clear();
    DsFit_Km_PY_vec.clear();
    DsFit_Km_PZ_vec.clear();
    DsFit_Km_PP_vec.clear();
    DsFit_Km_PL_vec.clear();

    DsFit_phi_ETA_vec.clear();
    DsFit_phi_PHI_vec.clear();
    DsFit_phi_P_vec.clear();
    DsFit_phi_PT_vec.clear();
    DsFit_phi_PX_vec.clear();
    DsFit_phi_PY_vec.clear();
    DsFit_phi_PZ_vec.clear();
    DsFit_phi_PP_vec.clear();
    DsFit_phi_PL_vec.clear();
    DsFit_phi_M_vec.clear();

    DsFit_pi_ETA_vec.clear();
    DsFit_pi_PHI_vec.clear();
    DsFit_pi_P_vec.clear();
    DsFit_pi_PT_vec.clear();
    DsFit_pi_PX_vec.clear();
    DsFit_pi_PY_vec.clear();
    DsFit_pi_PZ_vec.clear();
    DsFit_pi_PP_vec.clear();
    DsFit_pi_PL_vec.clear();

    DsFit_Mconstraint_Ds_M_vec.clear();
}

void PFRecoTree::Kp_Reset()
{
    Kp_ETA = null;
    Kp_PHI = null;
    Kp_ORIVX_X = null;
    Kp_ORIVX_Y = null;
    Kp_ORIVX_Z = null;
    Kp_P = null;
    Kp_PT = null;
    Kp_PX = null;
    Kp_PY = null;
    Kp_PZ = null;
}

void PFRecoTree::Km_Reset()
{
    Km_ETA = null;
    Km_PHI = null;
    Km_ORIVX_X = null;
    Km_ORIVX_Y = null;
    Km_ORIVX_Z = null;
    Km_P = null;
    Km_PT = null;
    Km_PX = null;
    Km_PY = null;
    Km_PZ = null;
    dR_Kp_Km = null;
    d_Kp_Km = null;
}

void PFRecoTree::phi_Reset()
{
    phiFit_Kp_ETA = null;
    phiFit_Kp_PHI = null;
    phiFit_Kp_P = null;
    phiFit_Kp_PT = null;
    phiFit_Kp_PX = null;
    phiFit_Kp_PY = null;
    phiFit_Kp_PZ = null;

    phiFit_Km_ETA = null;
    phiFit_Km_PHI = null;
    phiFit_Km_P = null;
    phiFit_Km_PT = null;
    phiFit_Km_PX = null;
    phiFit_Km_PY = null;
    phiFit_Km_PZ = null;
    
    phi_ETA = null;
    phi_PHI = null;
    phi_P = null;
    phi_PT = null;
    phi_PX = null;
    phi_PY = null;
    phi_PZ = null;
    phi_CHI2 = null;
    phi_M = null;
    phi_ENDVX_X = null;
    phi_ENDVX_Y = null;
    phi_ENDVX_Z = null;
    
    Kp_PP = null;
    Kp_PL = null;
    Km_PP = null;
    Km_PL = null;
    phiFit_Kp_PP = null;
    phiFit_Kp_PL = null;
    phiFit_Km_PP = null;
    phiFit_Km_PL = null;

    dR_Kp_phi = null;
    dR_Km_phi = null;
    d_Kp_phi = null;
    d_Km_phi = null;
}

void PFRecoTree::pi_Reset()
{
    pi_ETA = null;
    pi_PHI = null;
    pi_ORIVX_X = null;
    pi_ORIVX_Y = null;
    pi_ORIVX_Z = null;
    pi_P = null;
    pi_PT = null;
    pi_PX = null;
    pi_PY = null;
    pi_PZ = null;
    
    dR_Kp_pi = null;
    dR_Km_pi = null;
    dR_phi_pi = null;
    d_Kp_pi = null;
    d_Km_pi = null;
    d_phi_pi = null;
}

void PFRecoTree::Ds_Reset()
{
    
    DsFit_Kp_ETA = null;
    DsFit_Kp_PHI = null;
    DsFit_Kp_P = null;
    DsFit_Kp_PT = null;
    DsFit_Kp_PX = null;
    DsFit_Kp_PY = null;
    DsFit_Kp_PZ = null;

    DsFit_Km_ETA = null;
    DsFit_Km_PHI = null;
    DsFit_Km_P = null;
    DsFit_Km_PT = null;
    DsFit_Km_PX = null;
    DsFit_Km_PY = null;
    DsFit_Km_PZ = null;

    DsFit_phi_ETA = null;
    DsFit_phi_PHI = null;
    DsFit_phi_P = null;
    DsFit_phi_PT = null;
    DsFit_phi_PX = null;
    DsFit_phi_PY = null;
    DsFit_phi_PZ = null;
    DsFit_phi_M = null;

    DsFit_pi_ETA = null;
    DsFit_pi_PHI = null;
    DsFit_pi_P = null;
    DsFit_pi_PT = null;
    DsFit_pi_PX = null;
    DsFit_pi_PY = null;
    DsFit_pi_PZ = null;

    Ds_ETA = null;
    Ds_PHI = null;
    Ds_P = null;
    Ds_PT = null;
    Ds_PX = null;
    Ds_PY = null;
    Ds_PZ = null;
    Ds_CHI2 = null;
    Ds_M = null;
    Ds_ENDVX_X = null;
    Ds_ENDVX_Y = null;
    Ds_ENDVX_Z = null;

    DsFit_Mconstraint_Ds_M = null;

    phi_PP = null;
    phi_PL = null;
    pi_PP = null;
    pi_PL = null;
    
    DsFit_phi_PP = null;
    DsFit_phi_PL = null;
    DsFit_pi_PP = null;
    DsFit_pi_PL = null;
    
    DsFit_Kp_PP = null;
    DsFit_Kp_PL = null;
    DsFit_Km_PP = null;
    DsFit_Km_PL = null;

    dR_Kp_Ds = null;
    dR_Km_Ds = null;
    dR_phi_Ds = null;
    dR_pi_Ds = null;
    d_Kp_Ds = null;
    d_Km_Ds = null;
    d_phi_Ds = null;
    d_pi_Ds = null;
}

void PFRecoTree::Fill_Vector()
{
    Kp_ETA_vec.push_back(Kp_ETA);
    Kp_PHI_vec.push_back(Kp_PHI);
    Kp_ORIVX_X_vec.push_back(Kp_ORIVX_X);
    Kp_ORIVX_Y_vec.push_back(Kp_ORIVX_Y);
    Kp_ORIVX_Z_vec.push_back(Kp_ORIVX_Z);
    Kp_P_vec.push_back(Kp_P);
    Kp_PT_vec.push_back(Kp_PT);
    Kp_PX_vec.push_back(Kp_PX);
    Kp_PY_vec.push_back(Kp_PY);
    Kp_PZ_vec.push_back(Kp_PZ);
    Kp_PP_vec.push_back(Kp_PP);
    Kp_PL_vec.push_back(Kp_PL);

    Km_ETA_vec.push_back(Km_ETA);
    Km_PHI_vec.push_back(Km_PHI);
    Km_ORIVX_X_vec.push_back(Km_ORIVX_X);
    Km_ORIVX_Y_vec.push_back(Km_ORIVX_Y);
    Km_ORIVX_Z_vec.push_back(Km_ORIVX_Z);
    Km_P_vec.push_back(Km_P);
    Km_PT_vec.push_back(Km_PT);
    Km_PX_vec.push_back(Km_PX);
    Km_PY_vec.push_back(Km_PY);
    Km_PZ_vec.push_back(Km_PZ);
    Km_PP_vec.push_back(Km_PP);
    Km_PL_vec.push_back(Km_PL);

    phi_ETA_vec.push_back(phi_ETA);
    phi_PHI_vec.push_back(phi_PHI);
    phi_P_vec.push_back(phi_P);
    phi_PT_vec.push_back(phi_PT);
    phi_PX_vec.push_back(phi_PX);
    phi_PY_vec.push_back(phi_PY);
    phi_PZ_vec.push_back(phi_PZ);
    phi_PP_vec.push_back(phi_PP);
    phi_PL_vec.push_back(phi_PL);
    phi_CHI2_vec.push_back(phi_CHI2);
    phi_M_vec.push_back(phi_M);
    phi_ENDVX_X_vec.push_back(phi_ENDVX_X);
    phi_ENDVX_Y_vec.push_back(phi_ENDVX_Y);
    phi_ENDVX_Z_vec.push_back(phi_ENDVX_Z);

    pi_ETA_vec.push_back(pi_ETA);
    pi_PHI_vec.push_back(pi_PHI);
    pi_ORIVX_X_vec.push_back(pi_ORIVX_X);
    pi_ORIVX_Y_vec.push_back(pi_ORIVX_Y);
    pi_ORIVX_Z_vec.push_back(pi_ORIVX_Z);
    pi_P_vec.push_back(pi_P);
    pi_PT_vec.push_back(pi_PT);
    pi_PX_vec.push_back(pi_PX);
    pi_PY_vec.push_back(pi_PY);
    pi_PZ_vec.push_back(pi_PZ);
    pi_PP_vec.push_back(pi_PP);
    pi_PL_vec.push_back(pi_PL);

    Ds_ETA_vec.push_back(Ds_ETA);
    Ds_PHI_vec.push_back(Ds_PHI);
    Ds_P_vec.push_back(Ds_P);
    Ds_PT_vec.push_back(Ds_PT);
    Ds_PX_vec.push_back(Ds_PX);
    Ds_PY_vec.push_back(Ds_PY);
    Ds_PZ_vec.push_back(Ds_PZ);
    Ds_CHI2_vec.push_back(Ds_CHI2);
    Ds_M_vec.push_back(Ds_M);
    Ds_ENDVX_X_vec.push_back(Ds_ENDVX_X);
    Ds_ENDVX_Y_vec.push_back(Ds_ENDVX_Y);
    Ds_ENDVX_Z_vec.push_back(Ds_ENDVX_Z);

    dR_Kp_Km_vec.push_back(dR_Kp_Km);
    dR_Kp_phi_vec.push_back(dR_Kp_phi);
    dR_Km_phi_vec.push_back(dR_Km_phi);
    dR_Kp_pi_vec.push_back(dR_Kp_pi);
    dR_Km_pi_vec.push_back(dR_Km_pi);
    dR_phi_pi_vec.push_back(dR_phi_pi);
    dR_Kp_Ds_vec.push_back(dR_Kp_Ds);
    dR_Km_Ds_vec.push_back(dR_Km_Ds);
    dR_phi_Ds_vec.push_back(dR_phi_Ds);
    dR_pi_Ds_vec.push_back(dR_pi_Ds);

    d_Kp_Km_vec.push_back(d_Kp_Km);
    d_Kp_phi_vec.push_back(d_Kp_phi);
    d_Km_phi_vec.push_back(d_Km_phi);
    d_Kp_pi_vec.push_back(d_Kp_pi);
    d_Km_pi_vec.push_back(d_Km_pi);
    d_phi_pi_vec.push_back(d_phi_pi);
    d_Kp_Ds_vec.push_back(d_Kp_Ds);
    d_Km_Ds_vec.push_back(d_Km_Ds);
    d_phi_Ds_vec.push_back(d_phi_Ds);
    d_pi_Ds_vec.push_back(d_pi_Ds);

    phiFit_Kp_ETA_vec.push_back(phiFit_Kp_ETA);
    phiFit_Kp_PHI_vec.push_back(phiFit_Kp_PHI);
    phiFit_Kp_P_vec.push_back(phiFit_Kp_P);
    phiFit_Kp_PT_vec.push_back(phiFit_Kp_PT);
    phiFit_Kp_PX_vec.push_back(phiFit_Kp_PX);
    phiFit_Kp_PY_vec.push_back(phiFit_Kp_PY);
    phiFit_Kp_PZ_vec.push_back(phiFit_Kp_PZ);
    phiFit_Kp_PP_vec.push_back(phiFit_Kp_PP);
    phiFit_Kp_PL_vec.push_back(phiFit_Kp_PL);

    phiFit_Km_ETA_vec.push_back(phiFit_Km_ETA);
    phiFit_Km_PHI_vec.push_back(phiFit_Km_PHI);
    phiFit_Km_P_vec.push_back(phiFit_Km_P);
    phiFit_Km_PT_vec.push_back(phiFit_Km_PT);
    phiFit_Km_PX_vec.push_back(phiFit_Km_PX);
    phiFit_Km_PY_vec.push_back(phiFit_Km_PY);
    phiFit_Km_PZ_vec.push_back(phiFit_Km_PZ);
    phiFit_Km_PP_vec.push_back(phiFit_Km_PP);
    phiFit_Km_PL_vec.push_back(phiFit_Km_PL);

    DsFit_Kp_ETA_vec.push_back(DsFit_Kp_ETA);
    DsFit_Kp_PHI_vec.push_back(DsFit_Kp_PHI);
    DsFit_Kp_P_vec.push_back(DsFit_Kp_P);
    DsFit_Kp_PT_vec.push_back(DsFit_Kp_PT);
    DsFit_Kp_PX_vec.push_back(DsFit_Kp_PX);
    DsFit_Kp_PY_vec.push_back(DsFit_Kp_PY);
    DsFit_Kp_PZ_vec.push_back(DsFit_Kp_PZ);
    DsFit_Kp_PP_vec.push_back(DsFit_Kp_PP);
    DsFit_Kp_PL_vec.push_back(DsFit_Kp_PL);

    DsFit_Km_ETA_vec.push_back(DsFit_Km_ETA);
    DsFit_Km_PHI_vec.push_back(DsFit_Km_PHI);
    DsFit_Km_P_vec.push_back(DsFit_Km_P);
    DsFit_Km_PT_vec.push_back(DsFit_Km_PT);
    DsFit_Km_PX_vec.push_back(DsFit_Km_PX);
    DsFit_Km_PY_vec.push_back(DsFit_Km_PY);
    DsFit_Km_PZ_vec.push_back(DsFit_Km_PZ);
    DsFit_Km_PP_vec.push_back(DsFit_Km_PP);
    DsFit_Km_PL_vec.push_back(DsFit_Km_PL);

    DsFit_phi_ETA_vec.push_back(DsFit_phi_ETA);
    DsFit_phi_PHI_vec.push_back(DsFit_phi_PHI);
    DsFit_phi_P_vec.push_back(DsFit_phi_P);
    DsFit_phi_PT_vec.push_back(DsFit_phi_PT);
    DsFit_phi_PX_vec.push_back(DsFit_phi_PX);
    DsFit_phi_PY_vec.push_back(DsFit_phi_PY);
    DsFit_phi_PZ_vec.push_back(DsFit_phi_PZ);
    DsFit_phi_PP_vec.push_back(DsFit_phi_PP);
    DsFit_phi_PL_vec.push_back(DsFit_phi_PL);
    DsFit_phi_M_vec.push_back(DsFit_phi_M);

    DsFit_pi_ETA_vec.push_back(DsFit_pi_ETA);
    DsFit_pi_PHI_vec.push_back(DsFit_pi_PHI);
    DsFit_pi_P_vec.push_back(DsFit_pi_P);
    DsFit_pi_PT_vec.push_back(DsFit_pi_PT);
    DsFit_pi_PX_vec.push_back(DsFit_pi_PX);
    DsFit_pi_PY_vec.push_back(DsFit_pi_PY);
    DsFit_pi_PZ_vec.push_back(DsFit_pi_PZ);
    DsFit_pi_PP_vec.push_back(DsFit_pi_PP);
    DsFit_pi_PL_vec.push_back(DsFit_pi_PL);

    DsFit_Mconstraint_Ds_M_vec.push_back(DsFit_Mconstraint_Ds_M);
}

void PFRecoTree::CreateBranches()
{
    tree->Branch("num_reco_phi", &num_reco_phi);
    tree->Branch("num_reco_Ds", &num_reco_Ds);

    tree->Branch("Kp_ETA", &Kp_ETA_vec);
    tree->Branch("Kp_PHI", &Kp_PHI_vec);
    tree->Branch("Kp_ORIVX_X", &Kp_ORIVX_X_vec);
    tree->Branch("Kp_ORIVX_Y", &Kp_ORIVX_Y_vec);
    tree->Branch("Kp_ORIVX_Z", &Kp_ORIVX_Z_vec);
    tree->Branch("Kp_P", &Kp_P_vec);
    tree->Branch("Kp_PT", &Kp_PT_vec);
    tree->Branch("Kp_PX", &Kp_PX_vec);
    tree->Branch("Kp_PY", &Kp_PY_vec);
    tree->Branch("Kp_PZ", &Kp_PZ_vec);
    tree->Branch("Kp_PP", &Kp_PP_vec);
    tree->Branch("Kp_PL", &Kp_PL_vec);

    tree->Branch("Km_ETA", &Km_ETA_vec);
    tree->Branch("Km_PHI", &Km_PHI_vec);
    tree->Branch("Km_ORIVX_X", &Km_ORIVX_X_vec);
    tree->Branch("Km_ORIVX_Y", &Km_ORIVX_Y_vec);
    tree->Branch("Km_ORIVX_Z", &Km_ORIVX_Z_vec);
    tree->Branch("Km_P", &Km_P_vec);
    tree->Branch("Km_PT", &Km_PT_vec);
    tree->Branch("Km_PX", &Km_PX_vec);
    tree->Branch("Km_PY", &Km_PY_vec);
    tree->Branch("Km_PZ", &Km_PZ_vec);
    tree->Branch("Km_PP", &Km_PP_vec);
    tree->Branch("Km_PL", &Km_PL_vec);

    tree->Branch("phi_ETA", &phi_ETA_vec);
    tree->Branch("phi_PHI", &phi_PHI_vec);
    tree->Branch("phi_P", &phi_P_vec);
    tree->Branch("phi_PT", &phi_PT_vec);
    tree->Branch("phi_PX", &phi_PX_vec);
    tree->Branch("phi_PY", &phi_PY_vec);
    tree->Branch("phi_PZ", &phi_PZ_vec);
    tree->Branch("phi_PP", &phi_PP_vec);
    tree->Branch("phi_PL", &phi_PL_vec);
    tree->Branch("phi_CHI2", &phi_CHI2_vec);
    tree->Branch("phi_M", &phi_M_vec);
    tree->Branch("phi_ENDVX_X", &phi_ENDVX_X_vec);
    tree->Branch("phi_ENDVX_Y", &phi_ENDVX_Y_vec);
    tree->Branch("phi_ENDVX_Z", &phi_ENDVX_Z_vec);

    tree->Branch("pi_ETA", &pi_ETA_vec);
    tree->Branch("pi_PHI", &pi_PHI_vec);
    tree->Branch("pi_ORIVX_X", &pi_ORIVX_X_vec);
    tree->Branch("pi_ORIVX_Y", &pi_ORIVX_Y_vec);
    tree->Branch("pi_ORIVX_Z", &pi_ORIVX_Z_vec);
    tree->Branch("pi_P", &pi_P_vec);
    tree->Branch("pi_PT", &pi_PT_vec);
    tree->Branch("pi_PX", &pi_PX_vec);
    tree->Branch("pi_PY", &pi_PY_vec);
    tree->Branch("pi_PZ", &pi_PZ_vec);
    tree->Branch("pi_PP", &pi_PP_vec);
    tree->Branch("pi_PL", &pi_PL_vec);

    tree->Branch("Ds_ETA", &Ds_ETA_vec);
    tree->Branch("Ds_PHI", &Ds_PHI_vec);
    tree->Branch("Ds_P", &Ds_P_vec);
    tree->Branch("Ds_PT", &Ds_PT_vec);
    tree->Branch("Ds_PX", &Ds_PX_vec);
    tree->Branch("Ds_PY", &Ds_PY_vec);
    tree->Branch("Ds_PZ", &Ds_PZ_vec);
    tree->Branch("Ds_CHI2", &Ds_CHI2_vec);
    tree->Branch("Ds_M", &Ds_M_vec);
    tree->Branch("Ds_ENDVX_X", &Ds_ENDVX_X_vec);
    tree->Branch("Ds_ENDVX_Y", &Ds_ENDVX_Y_vec);
    tree->Branch("Ds_ENDVX_Z", &Ds_ENDVX_Z_vec);

    tree->Branch("dR_Kp_Km", &dR_Kp_Km_vec);
    tree->Branch("dR_Kp_phi", &dR_Kp_phi_vec);
    tree->Branch("dR_Km_phi", &dR_Km_phi_vec);
    tree->Branch("dR_Kp_pi", &dR_Kp_pi_vec);
    tree->Branch("dR_Km_pi", &dR_Km_pi_vec);
    tree->Branch("dR_phi_pi", &dR_phi_pi_vec);
    tree->Branch("dR_Kp_Ds", &dR_Kp_Ds_vec);
    tree->Branch("dR_Km_Ds", &dR_Km_Ds_vec);
    tree->Branch("dR_phi_Ds", &dR_phi_Ds_vec);
    tree->Branch("dR_pi_Ds", &dR_pi_Ds_vec);

    tree->Branch("d_Kp_Km", &d_Kp_Km_vec);
    tree->Branch("d_Kp_phi", &d_Kp_phi_vec);
    tree->Branch("d_Km_phi", &d_Km_phi_vec);
    tree->Branch("d_Kp_pi", &d_Kp_pi_vec);
    tree->Branch("d_Km_pi", &d_Km_pi_vec);
    tree->Branch("d_phi_pi", &d_phi_pi_vec);
    tree->Branch("d_Kp_Ds", &d_Kp_Ds_vec);
    tree->Branch("d_Km_Ds", &d_Km_Ds_vec);
    tree->Branch("d_phi_Ds", &d_phi_Ds_vec);
    tree->Branch("d_pi_Ds", &d_pi_Ds_vec);

    tree->Branch("phiFit_Kp_ETA", &phiFit_Kp_ETA_vec);
    tree->Branch("phiFit_Kp_PHI", &phiFit_Kp_PHI_vec);
    tree->Branch("phiFit_Kp_P", &phiFit_Kp_P_vec);
    tree->Branch("phiFit_Kp_PT", &phiFit_Kp_PT_vec);
    tree->Branch("phiFit_Kp_PX", &phiFit_Kp_PX_vec);
    tree->Branch("phiFit_Kp_PY", &phiFit_Kp_PY_vec);
    tree->Branch("phiFit_Kp_PZ", &phiFit_Kp_PZ_vec);
    tree->Branch("phiFit_Kp_PP", &phiFit_Kp_PP_vec);
    tree->Branch("phiFit_Kp_PL", &phiFit_Kp_PL_vec);

    tree->Branch("phiFit_Km_ETA", &phiFit_Km_ETA_vec);
    tree->Branch("phiFit_Km_PHI", &phiFit_Km_PHI_vec);
    tree->Branch("phiFit_Km_P", &phiFit_Km_P_vec);
    tree->Branch("phiFit_Km_PT", &phiFit_Km_PT_vec);
    tree->Branch("phiFit_Km_PX", &phiFit_Km_PX_vec);
    tree->Branch("phiFit_Km_PY", &phiFit_Km_PY_vec);
    tree->Branch("phiFit_Km_PZ", &phiFit_Km_PZ_vec);
    tree->Branch("phiFit_Km_PP", &phiFit_Km_PP_vec);
    tree->Branch("phiFit_Km_PL", &phiFit_Km_PL_vec);

    tree->Branch("DsFit_Kp_ETA", &DsFit_Kp_ETA_vec);
    tree->Branch("DsFit_Kp_PHI", &DsFit_Kp_PHI_vec);
    tree->Branch("DsFit_Kp_P", &DsFit_Kp_P_vec);
    tree->Branch("DsFit_Kp_PT", &DsFit_Kp_PT_vec);
    tree->Branch("DsFit_Kp_PX", &DsFit_Kp_PX_vec);
    tree->Branch("DsFit_Kp_PY", &DsFit_Kp_PY_vec);
    tree->Branch("DsFit_Kp_PZ", &DsFit_Kp_PZ_vec);
    tree->Branch("DsFit_Kp_PP", &DsFit_Kp_PP_vec);
    tree->Branch("DsFit_Kp_PL", &DsFit_Kp_PL_vec);

    tree->Branch("DsFit_Km_ETA", &DsFit_Km_ETA_vec);
    tree->Branch("DsFit_Km_PHI", &DsFit_Km_PHI_vec);
    tree->Branch("DsFit_Km_P", &DsFit_Km_P_vec);
    tree->Branch("DsFit_Km_PT", &DsFit_Km_PT_vec);
    tree->Branch("DsFit_Km_PX", &DsFit_Km_PX_vec);
    tree->Branch("DsFit_Km_PY", &DsFit_Km_PY_vec);
    tree->Branch("DsFit_Km_PZ", &DsFit_Km_PZ_vec);
    tree->Branch("DsFit_Km_PP", &DsFit_Km_PP_vec);
    tree->Branch("DsFit_Km_PL", &DsFit_Km_PL_vec);

    tree->Branch("DsFit_phi_ETA", &DsFit_phi_ETA_vec);
    tree->Branch("DsFit_phi_PHI", &DsFit_phi_PHI_vec);
    tree->Branch("DsFit_phi_P", &DsFit_phi_P_vec);
    tree->Branch("DsFit_phi_PT", &DsFit_phi_PT_vec);
    tree->Branch("DsFit_phi_PX", &DsFit_phi_PX_vec);
    tree->Branch("DsFit_phi_PY", &DsFit_phi_PY_vec);
    tree->Branch("DsFit_phi_PZ", &DsFit_phi_PZ_vec);
    tree->Branch("DsFit_phi_PP", &DsFit_phi_PP_vec);
    tree->Branch("DsFit_phi_PL", &DsFit_phi_PL_vec);
    tree->Branch("DsFit_phi_M", &DsFit_phi_M_vec);

    tree->Branch("DsFit_pi_ETA", &DsFit_pi_ETA_vec);
    tree->Branch("DsFit_pi_PHI", &DsFit_pi_PHI_vec);
    tree->Branch("DsFit_pi_P", &DsFit_pi_P_vec);
    tree->Branch("DsFit_pi_PT", &DsFit_pi_PT_vec);
    tree->Branch("DsFit_pi_PX", &DsFit_pi_PX_vec);
    tree->Branch("DsFit_pi_PY", &DsFit_pi_PY_vec);
    tree->Branch("DsFit_pi_PZ", &DsFit_pi_PZ_vec);
    tree->Branch("DsFit_pi_PP", &DsFit_pi_PP_vec);
    tree->Branch("DsFit_pi_PL", &DsFit_pi_PL_vec);

    tree->Branch("DsFit_Mconstraint_Ds_M", &DsFit_Mconstraint_Ds_M_vec);
};

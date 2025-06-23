#include <TChain.h>

const double m_Ds = 1.96835;
const double m_phi = 1.019460;
const double m_K = 0.493677;
const double m_pi = 0.13957039;
const double pst_phi = m_phi/2;
const double pst_Ds = sqrt((pow(m_Ds,4)+pow(m_phi,4)+pow(m_pi,4)-2*pow(m_Ds*m_phi,2)-2*pow(m_Ds*m_pi,2)-2*pow(m_phi*m_pi,2))/(4*pow(m_Ds,2)));
const double Est_phi = std::sqrt(pow(m_phi,2) + pow(pst_Ds,2));
const double Est_pi = std::sqrt(pow(m_pi,2) + pow(pst_Ds,2));


int new_sel_eff(){

    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SS1.root");

    std::vector<double> *match_Kp_PT = nullptr;
    std::vector<double> *match_Km_PT = nullptr;
    std::vector<double> *match_pi_PT = nullptr;
    
    std::vector<double> *match_Kp_P = nullptr;
    std::vector<double> *match_Km_P = nullptr;
    std::vector<double> *match_pi_P = nullptr;

    std::vector<double> *match_phiFit_Kp_PP = nullptr;
    std::vector<double> *match_phiFit_Km_PP = nullptr;
    std::vector<double> *match_DsFit_phi_PP = nullptr;
    std::vector<double> *match_DsFit_pi_PP = nullptr;

    std::vector<double> *match_phiFit_Kp_PL = nullptr;
    std::vector<double> *match_phiFit_Km_PL = nullptr;
    std::vector<double> *match_DsFit_phi_PL = nullptr;
    std::vector<double> *match_DsFit_pi_PL = nullptr;
    
    std::vector<double> *match_phiFit_phi_P = nullptr;
    std::vector<double> *match_DsFit_Ds_P = nullptr;

    std::vector<double> *match_dR_Kp_Km = nullptr;
    std::vector<double> *match_dR_Kp_pi = nullptr;
    std::vector<double> *match_dR_Km_pi = nullptr;

    std::vector<double> *match_phiFit_dR_Kp_phi = nullptr;
    std::vector<double> *match_phiFit_dR_Km_phi = nullptr;
    std::vector<double> *match_phiFit_dR_pi_phi = nullptr;

    std::vector<double> *match_DsFit_dR_Kp_Ds = nullptr;
    std::vector<double> *match_DsFit_dR_Km_Ds = nullptr;
    std::vector<double> *match_DsFit_dR_pi_Ds = nullptr;
    std::vector<double> *match_DsFit_dR_phi_Ds = nullptr;
    
    std::vector<double> *match_phiFit_CHI2NDOF = nullptr;
    std::vector<double> *match_phiFit_phi_M = nullptr;
    std::vector<double> *match_DsFit_CHI2NDOF = nullptr;
    std::vector<double> *match_DsFit_Ds_M = nullptr;


    mychain->SetBranchAddress("match_Kp_PT", &match_Kp_PT);
    mychain->SetBranchAddress("match_Km_PT", &match_Km_PT);
    mychain->SetBranchAddress("match_pi_PT", &match_pi_PT);
    
    mychain->SetBranchAddress("match_Kp_P", &match_Kp_P);
    mychain->SetBranchAddress("match_Km_P", &match_Km_P);
    mychain->SetBranchAddress("match_pi_P", &match_pi_P);

    mychain->SetBranchAddress("match_phiFit_Kp_PP", &match_phiFit_Kp_PP);
    mychain->SetBranchAddress("match_phiFit_Km_PP", &match_phiFit_Km_PP);
    mychain->SetBranchAddress("match_DsFit_phi_PP", &match_DsFit_phi_PP);
    mychain->SetBranchAddress("match_DsFit_pi_PP", &match_DsFit_pi_PP);

    mychain->SetBranchAddress("match_phiFit_Kp_PL", &match_phiFit_Kp_PL);
    mychain->SetBranchAddress("match_phiFit_Km_PL", &match_phiFit_Km_PL);
    mychain->SetBranchAddress("match_DsFit_phi_PL", &match_DsFit_phi_PL);
    mychain->SetBranchAddress("match_DsFit_pi_PL", &match_DsFit_pi_PL);
    
    mychain->SetBranchAddress("match_phiFit_phi_P", &match_phiFit_phi_P);
    mychain->SetBranchAddress("match_DsFit_Ds_P", &match_DsFit_Ds_P);
    
    mychain->SetBranchAddress("match_dR_Kp_Km", &match_dR_Kp_Km);
    mychain->SetBranchAddress("match_dR_Kp_pi", &match_dR_Kp_pi);
    mychain->SetBranchAddress("match_dR_Km_pi", &match_dR_Km_pi);
    mychain->SetBranchAddress("match_phiFit_dR_Kp_phi", &match_phiFit_dR_Kp_phi);
    mychain->SetBranchAddress("match_phiFit_dR_Km_phi", &match_phiFit_dR_Km_phi);
    mychain->SetBranchAddress("match_phiFit_dR_pi_phi", &match_phiFit_dR_pi_phi);
    mychain->SetBranchAddress("match_DsFit_dR_Kp_Ds", &match_DsFit_dR_Kp_Ds);
    mychain->SetBranchAddress("match_DsFit_dR_Km_Ds", &match_DsFit_dR_Km_Ds);
    mychain->SetBranchAddress("match_DsFit_dR_pi_Ds", &match_DsFit_dR_pi_Ds);
    mychain->SetBranchAddress("match_DsFit_dR_phi_Ds", &match_DsFit_dR_phi_Ds);
    
    mychain->SetBranchAddress("match_phiFit_CHI2NDOF", &match_phiFit_CHI2NDOF);
    mychain->SetBranchAddress("match_phiFit_phi_M", &match_phiFit_phi_M);
    mychain->SetBranchAddress("match_DsFit_CHI2NDOF", &match_DsFit_CHI2NDOF);
    mychain->SetBranchAddress("match_DsFit_Ds_M", &match_DsFit_Ds_M);

    int nevents = mychain->GetEntries();

    std::cout << nevents << std::endl;

    int n_match = 0;
    int n_Kp_PT = 0;
    int n_Km_PT = 0;
    int n_pi_PT = 0;
    int n_Kp_P = 0;
    int n_Km_P = 0;
    int n_pi_P = 0;
    int n_Kp_Km_PT = 0;
    int n_dR_Kp_Km = 0;
    int n_dR_Kp_pi = 0;
    int n_dR_Km_pi = 0;
    int n_phiFit_dR_Kp_phi = 0;
    int n_phiFit_dR_Km_phi = 0;
    int n_phiFit_dR_pi_phi = 0;
    int n_DsFit_dR_Kp_Ds = 0;
    int n_DsFit_dR_Km_Ds = 0;
    int n_DsFit_dR_pi_Ds = 0;
    int n_DsFit_dR_phi_Ds = 0;
    int n_phiFit_CHI2NDOF = 0;
    int n_phiFit_phi_M = 0;
    int n_DsFit_CHI2NDOF = 0;
    int n_DsFit_Ds_M = 0;
    int n_APplot = 0;

    for(int i=0; i<nevents; i++){

        mychain->GetEntry(i);

        if( match_Kp_PT->size() < 1 ) continue; n_match++;

        if( match_Kp_PT->at(0) < 0.5 ) continue; n_Kp_PT++;
        if( match_Km_PT->at(0) < 0.5 ) continue; n_Km_PT++;
        if( match_pi_PT->at(0) < 0.5 ) continue; n_pi_PT++;

        if( match_Kp_P->at(0) < 1 ) continue; n_Kp_P++;
        if( match_Km_P->at(0) < 1 ) continue; n_Km_P++;
        if( match_pi_P->at(0) < 1 ) continue; n_pi_P++;

        if( abs(match_Kp_PT->at(0) - match_Km_PT->at(0)) > 20 ) continue; n_Kp_Km_PT++;

        if( match_dR_Kp_Km->at(0) > 0.1) continue; n_dR_Kp_Km++;
        if( match_dR_Kp_pi->at(0) > 0.4) continue; n_dR_Kp_pi++;
        if( match_dR_Km_pi->at(0) > 0.4) continue; n_dR_Km_pi++;

        if( match_phiFit_dR_Kp_phi->at(0) > 0.05 ) continue; n_phiFit_dR_Kp_phi++;
        if( match_phiFit_dR_Km_phi->at(0) > 0.05 ) continue; n_phiFit_dR_Km_phi++;
        if( match_phiFit_dR_pi_phi->at(0) > 0.4 ) continue; n_phiFit_dR_pi_phi++;

        if( match_DsFit_dR_Kp_Ds->at(0) > 0.15 ) continue; n_DsFit_dR_Kp_Ds++;
        if( match_DsFit_dR_Km_Ds->at(0) > 0.15 ) continue; n_DsFit_dR_Km_Ds++;
        if( match_DsFit_dR_pi_Ds->at(0) > 0.4 ) continue; n_DsFit_dR_pi_Ds++;
        if( match_DsFit_dR_phi_Ds->at(0) > 0.15 ) continue; n_DsFit_dR_phi_Ds++;
        
        if( match_phiFit_CHI2NDOF->at(0) < 0 ) continue;
        if( match_phiFit_CHI2NDOF->at(0) > 10 ) continue; n_phiFit_CHI2NDOF++;
        if( match_phiFit_phi_M->at(0) < 0.99 ) continue;
        if( match_phiFit_phi_M->at(0) > 1.05 ) continue; n_phiFit_phi_M++;

        if( match_DsFit_CHI2NDOF->at(0) < 0 ) continue;
        if( match_DsFit_CHI2NDOF->at(0) > 10 ) continue; n_DsFit_CHI2NDOF++;
        if( match_DsFit_Ds_M->at(0) < 1.85 ) continue;
        if( match_DsFit_Ds_M->at(0) > 2.1 ) continue; n_DsFit_Ds_M++;

        double alpha_phi = (match_phiFit_Kp_PL->at(0) - match_phiFit_Km_PL->at(0)) / (match_phiFit_Kp_PL->at(0) + match_phiFit_Km_PL->at(0));
        double alpha_Ds = (match_DsFit_phi_PL->at(0) - match_DsFit_pi_PL->at(0)) / (match_DsFit_phi_PL->at(0) + match_DsFit_pi_PL->at(0));

        double beta_phi = std::sqrt(pow(match_phiFit_phi_P->at(0),2) / (pow(m_phi,2) + pow(match_phiFit_phi_P->at(0),2)));
        double beta_Ds = std::sqrt(pow(match_DsFit_Ds_P->at(0),2) / (pow(m_Ds,2) + pow(match_DsFit_Ds_P->at(0),2)));

        double var_phi = pow(match_phiFit_Kp_PP->at(0),2) + pow(alpha_phi*beta_phi*m_phi,2)/4;
        double var_Ds = pow(match_DsFit_phi_PP->at(0),2) + pow(beta_Ds*(Est_phi-Est_pi)-alpha_Ds*beta_Ds*(Est_phi+Est_pi),2)/4;

        if(var_phi < 0.006) continue;
        if(var_phi > 0.028) continue;
        if(var_Ds < 0.4) continue;
        if(var_Ds > 0.6) continue; n_APplot++;

    }

    std::cout << "n_match " << n_match << std::endl;
    std::cout << "n_Kp_PT " << n_Kp_PT << std::endl;
    std::cout << "n_Km_PT " << n_Km_PT << std::endl;
    std::cout << "n_pi_PT " << n_pi_PT << std::endl;
    std::cout << "n_Kp_P " << n_Kp_P << std::endl;
    std::cout << "n_Km_P " << n_Km_P << std::endl;
    std::cout << "n_pi_P " << n_pi_P << std::endl;
    std::cout << "n_Kp_Km_PT " << n_Kp_Km_PT << std::endl;
    std::cout << "n_dR_Kp_Km " << n_dR_Kp_Km << std::endl;
    std::cout << "n_dR_Kp_pi " << n_dR_Kp_pi << std::endl;
    std::cout << "n_dR_Km_pi " << n_dR_Km_pi << std::endl;
    std::cout << "n_phiFit_dR_Kp_phi " << n_phiFit_dR_Kp_phi << std::endl;
    std::cout << "n_phiFit_dR_Km_phi " << n_phiFit_dR_Km_phi << std::endl;
    std::cout << "n_phiFit_dR_pi_phi " << n_phiFit_dR_pi_phi << std::endl;
    std::cout << "n_DsFit_dR_Kp_Ds " << n_DsFit_dR_Kp_Ds << std::endl;
    std::cout << "n_DsFit_dR_Km_Ds " << n_DsFit_dR_Km_Ds << std::endl;
    std::cout << "n_DsFit_dR_pi_Ds " << n_DsFit_dR_pi_Ds << std::endl;
    std::cout << "n_DsFit_dR_phi_Ds " << n_DsFit_dR_phi_Ds << std::endl;
    std::cout << "n_phiFit_CHI2NDOF " << n_phiFit_CHI2NDOF << std::endl;
    std::cout << "n_phiFit_phi_M " << n_phiFit_phi_M << std::endl;
    std::cout << "n_DsFit_CHI2NDOF " << n_DsFit_CHI2NDOF << std::endl;
    std::cout << "n_DsFit_Ds_M " << n_DsFit_Ds_M << std::endl;
    std::cout << "n_APplot " << n_APplot << std::endl;

    return 0;
}

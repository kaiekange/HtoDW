#include <TChain.h>

const double m_Ds = 1.96835;
const double m_phi = 1.019460;
const double m_K = 0.493677;
const double m_pi = 0.13957039;
const double pst_phi = m_phi/2;
const double pst_Ds = sqrt((pow(m_Ds,4)+pow(m_phi,4)+pow(m_pi,4)-2*pow(m_Ds*m_phi,2)-2*pow(m_Ds*m_pi,2)-2*pow(m_phi*m_pi,2))/(4*pow(m_Ds,2)));
const double Est_phi = std::sqrt(pow(m_phi,2) + pow(pst_Ds,2));
const double Est_pi = std::sqrt(pow(m_pi,2) + pow(pst_Ds,2));


int APcut()
{
    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    std::vector<double> *best_phiFit_Kp_PP = nullptr;
    std::vector<double> *best_phiFit_Kp_PL = nullptr;
    std::vector<double> *best_phiFit_Km_PP = nullptr;
    std::vector<double> *best_phiFit_Km_PL = nullptr;
    std::vector<double> *best_phiFit_phi_P = nullptr;

    std::vector<double> *best_DsFit_phi_PP = nullptr;
    std::vector<double> *best_DsFit_phi_PL = nullptr;
    std::vector<double> *best_DsFit_pi_PP = nullptr;
    std::vector<double> *best_DsFit_pi_PL = nullptr;
    std::vector<double> *best_DsFit_Ds_P = nullptr;

    std::vector<bool> *best_match_entry = nullptr;

    mychain->SetBranchAddress("best_phiFit_Kp_PP", &best_phiFit_Kp_PP);
    mychain->SetBranchAddress("best_phiFit_Kp_PL", &best_phiFit_Kp_PL);
    mychain->SetBranchAddress("best_phiFit_Km_PP", &best_phiFit_Km_PP);
    mychain->SetBranchAddress("best_phiFit_Km_PL", &best_phiFit_Km_PL);
    mychain->SetBranchAddress("best_phiFit_phi_P", &best_phiFit_phi_P);

    mychain->SetBranchAddress("best_DsFit_phi_PP", &best_DsFit_phi_PP);
    mychain->SetBranchAddress("best_DsFit_phi_PL", &best_DsFit_phi_PL);
    mychain->SetBranchAddress("best_DsFit_pi_PP", &best_DsFit_pi_PP);
    mychain->SetBranchAddress("best_DsFit_pi_PL", &best_DsFit_pi_PL);
    mychain->SetBranchAddress("best_DsFit_Ds_P", &best_DsFit_Ds_P);

    mychain->SetBranchAddress("best_match_entry", &best_match_entry);

    int nentries = mychain->GetEntries();

    int n_match=0, n_nonmatch=0;

    for(int i=0; i<nentries; i++){

        mychain->GetEntry(i);

        if(best_match_entry->size() < 1) continue;

        double alpha_phi = (best_phiFit_Kp_PL->at(0) - best_phiFit_Km_PL->at(0)) / (best_phiFit_Kp_PL->at(0) + best_phiFit_Km_PL->at(0));
        double alpha_Ds = (best_DsFit_phi_PL->at(0) - best_DsFit_pi_PL->at(0)) / (best_DsFit_phi_PL->at(0) + best_DsFit_pi_PL->at(0));

        double beta_phi = std::sqrt(pow(best_phiFit_phi_P->at(0),2) / (pow(m_phi,2) + pow(best_phiFit_phi_P->at(0),2)));
        double beta_Ds = std::sqrt(pow(best_DsFit_Ds_P->at(0),2) / (pow(m_Ds,2) + pow(best_DsFit_Ds_P->at(0),2)));
        
        double var_phi = pow(best_phiFit_Kp_PP->at(0),2) + pow(alpha_phi*beta_phi*m_phi,2)/4;
        double var_Ds = pow(best_DsFit_phi_PP->at(0),2) + pow(beta_Ds*(Est_phi-Est_pi)-alpha_Ds*beta_Ds*(Est_phi+Est_pi),2)/4;

        if(var_phi < 0.006) continue;
        if(var_phi > 0.03) continue;
        if(var_Ds < 0.4) continue;
        if(var_Ds > 0.6) continue;

        if(best_match_entry->at(0) == 1) n_match++;
        else if(best_match_entry->at(0) == 0) n_nonmatch++;
    }

    std::cout << n_match << " " << n_nonmatch << std::endl;

    return 0;
}

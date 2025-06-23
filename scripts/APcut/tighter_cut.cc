#include <TChain.h>

const double m_Ds = 1.96835;
const double m_phi = 1.019460;
const double m_K = 0.493677;
const double m_pi = 0.13957039;
const double pst_phi = m_phi/2;
const double pst_Ds = sqrt((pow(m_Ds,4)+pow(m_phi,4)+pow(m_pi,4)-2*pow(m_Ds*m_phi,2)-2*pow(m_Ds*m_pi,2)-2*pow(m_phi*m_pi,2))/(4*pow(m_Ds,2)));
const double Est_phi = std::sqrt(pow(m_phi,2) + pow(pst_Ds,2));
const double Est_pi = std::sqrt(pow(m_pi,2) + pow(pst_Ds,2));

int checkAP()
{
    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    std::vector<double> *phiFit_Kp_PP = nullptr;
    std::vector<double> *phiFit_Kp_PL = nullptr;
    std::vector<double> *phiFit_Km_PP = nullptr;
    std::vector<double> *phiFit_Km_PL = nullptr;
    std::vector<double> *phiFit_phi_P = nullptr;

    std::vector<double> *DsFit_phi_PP = nullptr;
    std::vector<double> *DsFit_phi_PL = nullptr;
    std::vector<double> *DsFit_pi_PP = nullptr;
    std::vector<double> *DsFit_pi_PL = nullptr;
    std::vector<double> *DsFit_Ds_P = nullptr;
    std::vector<double> *DsFit_Ds_PT = nullptr;

    std::vector<bool> *match_entry = nullptr;
    bool best_match_entry;

    mychain->SetBranchAddress("phiFit_Kp_PP", &phiFit_Kp_PP);
    mychain->SetBranchAddress("phiFit_Kp_PL", &phiFit_Kp_PL);
    mychain->SetBranchAddress("phiFit_Km_PP", &phiFit_Km_PP);
    mychain->SetBranchAddress("phiFit_Km_PL", &phiFit_Km_PL);
    mychain->SetBranchAddress("phiFit_phi_P", &phiFit_phi_P);

    mychain->SetBranchAddress("DsFit_phi_PP", &DsFit_phi_PP);
    mychain->SetBranchAddress("DsFit_phi_PL", &DsFit_phi_PL);
    mychain->SetBranchAddress("DsFit_pi_PP", &DsFit_pi_PP);
    mychain->SetBranchAddress("DsFit_pi_PL", &DsFit_pi_PL);
    mychain->SetBranchAddress("DsFit_Ds_P", &DsFit_Ds_P);
    mychain->SetBranchAddress("DsFit_Ds_PT", &DsFit_Ds_PT);

    mychain->SetBranchAddress("match_entry", &match_entry);

    TFile * outfile = new TFile("./checkAP.root", "RECREATE");
    TTree * outtree = new TTree("Events", "Events");
    outtree->Branch("best_match_entry", & best_match_entry);

    int nentries = mychain->GetEntries();

    int n_match=0, n_nonmatch=0;

    for(int i=0; i<nentries; i++){

        best_match_entry = false;
        mychain->GetEntry(i);

        int vecsize = match_entry->size();

        double maxPT = 0;
        int maxidx = -1;

        for(int j=0; j<vecsize; j++){
            double alpha_phi = (phiFit_Kp_PL->at(j) - phiFit_Km_PL->at(j)) / (phiFit_Kp_PL->at(j) + phiFit_Km_PL->at(j));
            double alpha_Ds = (DsFit_phi_PL->at(j) - DsFit_pi_PL->at(j)) / (DsFit_phi_PL->at(j) + DsFit_pi_PL->at(j));

            double beta_phi = std::sqrt(pow(phiFit_phi_P->at(j),2) / (pow(m_phi,2) + pow(phiFit_phi_P->at(j),2)));
            double beta_Ds = std::sqrt(pow(DsFit_Ds_P->at(j),2) / (pow(m_Ds,2) + pow(DsFit_Ds_P->at(j),2)));

            double var_phi = pow(phiFit_Kp_PP->at(j),2) + pow(alpha_phi*beta_phi*m_phi,2)/4;
            double var_Ds = pow(DsFit_phi_PP->at(j),2) + pow(beta_Ds*(Est_phi-Est_pi)-alpha_Ds*beta_Ds*(Est_phi+Est_pi),2)/4;

            if(var_phi < 0.006) continue;
            if(var_phi > 0.028) continue;
            if(var_Ds < 0.4) continue;
            if(var_Ds > 0.6) continue;

            if(DsFit_Ds_PT->at(j) > maxPT){
                maxPT = DsFit_Ds_PT->at(j);
                maxidx = j;
            }
        }

        if(maxidx > -1){
            best_match_entry = match_entry->at(maxidx);
            outtree->Fill();
        }
    }

    outtree->Write();
    delete outtree;
	outfile->Close();
	delete outfile;
    /* std::cout << n_match << " " << n_nonmatch << std::endl; */

    return 0;
}

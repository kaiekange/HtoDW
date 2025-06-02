#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

int checkpt(int file_num) {
    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add(Form("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/output_%d.root",file_num));

    TFile * outfile = new TFile(Form("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/checkpt_%d.root",file_num), "RECREATE");
    TTree * outtree = new TTree("Events", "Events");

    vector<bool> *match_entry_vec = nullptr;
    vector<double> *phiFit_phi_M_vec = nullptr;
    vector<double> *DsFit_Ds_M_vec = nullptr;
    vector<double> *DsFit_Ds_PT_vec = nullptr;
    vector<double> *DsFit_Ds_P_vec = nullptr;
    vector<double> *Ds_PT_vec = nullptr;
    vector<double> *Kp_PT_vec = nullptr;
    vector<double> *Km_PT_vec = nullptr;
    vector<double> *pi_PT_vec = nullptr;
  
    bool DsFit_Ds_PT_match;
    bool DsFit_Ds_P_match;
    bool Ds_PT_match;
    bool ScalarSum_PT_match;
    bool ScalarSum_PT2_match;

    mychain->SetBranchAddress("match_entry", &match_entry_vec);
    mychain->SetBranchAddress("phiFit_phi_M", &phiFit_phi_M_vec);
    mychain->SetBranchAddress("DsFit_Ds_M", &DsFit_Ds_M_vec);
    mychain->SetBranchAddress("DsFit_Ds_PT", &DsFit_Ds_PT_vec);
    mychain->SetBranchAddress("DsFit_Ds_P", &DsFit_Ds_P_vec);
    mychain->SetBranchAddress("Ds_PT", &Ds_PT_vec);
    mychain->SetBranchAddress("Kp_PT", &Kp_PT_vec);
    mychain->SetBranchAddress("Km_PT", &Km_PT_vec);
    mychain->SetBranchAddress("pi_PT", &pi_PT_vec);
  
    outtree->Branch("DsFit_Ds_PT_match", &DsFit_Ds_PT_match);
    outtree->Branch("DsFit_Ds_P_match", &DsFit_Ds_P_match);
    outtree->Branch("Ds_PT_match", &Ds_PT_match);
    outtree->Branch("ScalarSum_PT_match", &ScalarSum_PT_match);
    outtree->Branch("ScalarSum_PT2_match", &ScalarSum_PT2_match);

    int nentries = mychain->GetEntries();

    for(int i=0; i<nentries; i++){
        mychain->GetEntry(i);

        double max_DsFit_Ds_PT = -1;
        double max_DsFit_Ds_P = -1;
        double max_Ds_PT = -1;
        double max_ScalarSum_PT = -1;
        double max_ScalarSum_PT2 = -1;

        DsFit_Ds_PT_match = false;
        DsFit_Ds_P_match = false;
        Ds_PT_match = false;
        ScalarSum_PT_match = false;
        ScalarSum_PT2_match = false;

        int vecsize = match_entry_vec->size();

        for(int j=0; j<vecsize; j++){

            if( (phiFit_phi_M_vec->at(j) > 0.99) && (phiFit_phi_M_vec->at(j) < 1.05) && (DsFit_Ds_M_vec->at(j) > 1.85) && (DsFit_Ds_M_vec->at(j) < 2.1) ){

                if(DsFit_Ds_PT_vec->at(j) > max_DsFit_Ds_PT){
                    max_DsFit_Ds_PT = DsFit_Ds_PT_vec->at(j);
                    DsFit_Ds_PT_match = match_entry_vec->at(j);
                }

                if(DsFit_Ds_P_vec->at(j) > max_DsFit_Ds_P){
                    max_DsFit_Ds_P = DsFit_Ds_P_vec->at(j);
                    DsFit_Ds_P_match = match_entry_vec->at(j);
                }

                if(Ds_PT_vec->at(j) > max_Ds_PT){
                    max_Ds_PT = Ds_PT_vec->at(j);
                    Ds_PT_match = match_entry_vec->at(j);
                }

                if( (Kp_PT_vec->at(j) + Km_PT_vec->at(j) + pi_PT_vec->at(j)) > max_ScalarSum_PT ){
                    max_ScalarSum_PT = (Kp_PT_vec->at(j) + Km_PT_vec->at(j) + pi_PT_vec->at(j));
                    ScalarSum_PT_match = match_entry_vec->at(j);
                }

                if( (pow(Kp_PT_vec->at(j),2) + pow(Km_PT_vec->at(j),2) + pow(pi_PT_vec->at(j),2)) > max_ScalarSum_PT2 ){
                    max_ScalarSum_PT2 = (pow(Kp_PT_vec->at(j),2) + pow(Km_PT_vec->at(j),2) + pow(pi_PT_vec->at(j),2));
                    ScalarSum_PT2_match = match_entry_vec->at(j);
                }
            }
        }

        outtree->Fill();
        match_entry_vec->clear();
        phiFit_phi_M_vec->clear();
        DsFit_Ds_M_vec->clear();
        DsFit_Ds_PT_vec->clear();
        DsFit_Ds_P_vec->clear();
        Ds_PT_vec->clear();
        Kp_PT_vec->clear();
        Km_PT_vec->clear();
        pi_PT_vec->clear();
    }
   
    outtree->Write();
    delete outtree;
	outfile->Close();
	delete outfile;

    delete mychain;

    return 0;
}

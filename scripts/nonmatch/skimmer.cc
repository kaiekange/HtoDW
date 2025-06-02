#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

int skimmer(int file_num) {
    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add(Form("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/output_%d.root",file_num));

    TFile * outfile = new TFile(Form("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/skimmed_%d.root",file_num), "RECREATE");
    TTree * outtree = new TTree("Events", "Events");

    vector<bool> *match_entry_vec = nullptr;
    vector<bool> *non_match_entry_vec = nullptr;
    vector<bool> *Kp_match_vec = nullptr;
    vector<bool> *Km_match_vec = nullptr;
    vector<bool> *pi_match_vec = nullptr;
    vector<double> *phi_M_vec = nullptr;
    vector<double> *phiFit_phi_M_vec = nullptr;
    vector<double> *DsFit_phi_M_vec = nullptr;
    vector<double> *Ds_M_vec = nullptr;
    vector<double> *phiFit_Ds_M_vec = nullptr;
    vector<double> *DsFit_Ds_M_vec = nullptr;
    vector<double> *DsFit_Mconstraint_Ds_M_vec = nullptr;
  
    bool match_entry = false;
    bool non_match_entry = false;
    bool Kp_match = false;
    bool Km_match = false;
    bool pi_match = false;
    double phi_M = 0;
    double phiFit_phi_M = 0;
    double DsFit_phi_M = 0;
    double Ds_M = 0;
    double phiFit_Ds_M = 0;
    double DsFit_Ds_M = 0;
    double DsFit_Mconstraint_Ds_M = 0;

    mychain->SetBranchAddress("match_entry", &match_entry_vec);
    mychain->SetBranchAddress("non_match_entry", &non_match_entry_vec);
    mychain->SetBranchAddress("Kp_match", &Kp_match_vec);
    mychain->SetBranchAddress("Km_match", &Km_match_vec);
    mychain->SetBranchAddress("pi_match", &pi_match_vec);
    mychain->SetBranchAddress("phi_M", &phi_M_vec);
    mychain->SetBranchAddress("phiFit_phi_M", &phiFit_phi_M_vec);
    mychain->SetBranchAddress("DsFit_phi_M", &DsFit_phi_M_vec);
    mychain->SetBranchAddress("Ds_M", &Ds_M_vec);
    mychain->SetBranchAddress("phiFit_Ds_M", &phiFit_Ds_M_vec);
    mychain->SetBranchAddress("DsFit_Ds_M", &DsFit_Ds_M_vec);
    mychain->SetBranchAddress("DsFit_Mconstraint_Ds_M", &DsFit_Mconstraint_Ds_M_vec);
  
    outtree->Branch("match_entry", &match_entry);
    outtree->Branch("non_match_entry", &non_match_entry);
    outtree->Branch("Kp_match", &Kp_match);
    outtree->Branch("Km_match", &Km_match);
    outtree->Branch("pi_match", &pi_match);
    outtree->Branch("phi_M", &phi_M);
    outtree->Branch("phiFit_phi_M", &phiFit_phi_M);
    outtree->Branch("DsFit_phi_M", &DsFit_phi_M);
    outtree->Branch("Ds_M", &Ds_M);
    outtree->Branch("phiFit_Ds_M", &phiFit_Ds_M);
    outtree->Branch("DsFit_Ds_M", &DsFit_Ds_M);
    outtree->Branch("DsFit_Mconstraint_Ds_M", &DsFit_Mconstraint_Ds_M);

    int nentries = mychain->GetEntries();

    for(int i=0; i<nentries; i++){
        mychain->GetEntry(i);
        int vecsize = match_entry_vec->size();

        for(int j=0; j<vecsize; j++){
            match_entry = match_entry_vec->at(j);
            non_match_entry = non_match_entry_vec->at(j);
            Kp_match = Kp_match_vec->at(j);
            Km_match = Km_match_vec->at(j);
            pi_match = pi_match_vec->at(j);
            phi_M = phi_M_vec->at(j);
            phiFit_phi_M = phiFit_phi_M_vec->at(j);
            DsFit_phi_M = DsFit_phi_M_vec->at(j);
            Ds_M = Ds_M_vec->at(j);
            phiFit_Ds_M = phiFit_Ds_M_vec->at(j);
            DsFit_Ds_M = DsFit_Ds_M_vec->at(j);
            DsFit_Mconstraint_Ds_M = DsFit_Mconstraint_Ds_M_vec->at(j);
            outtree->Fill();
        }
        match_entry_vec->clear();
        non_match_entry_vec->clear();
        Kp_match_vec->clear();
        Km_match_vec->clear();
        pi_match_vec->clear();
        phi_M_vec->clear();
        phiFit_phi_M_vec->clear();
        DsFit_phi_M_vec->clear();
        Ds_M_vec->clear();
        phiFit_Ds_M_vec->clear();
        DsFit_Ds_M_vec->clear();
        DsFit_Mconstraint_Ds_M_vec->clear();
    }
   
    outtree->Write();
    delete outtree;
	outfile->Close();
	delete outfile;

    delete mychain;

    return 0;
}

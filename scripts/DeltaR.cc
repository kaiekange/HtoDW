#include <TChain.h>
#include <TTree.h>

double getDeltaR(double eta1, double phi1, double eta2, double phi2){      
   double DeltaPhi = abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

int DeltaR(){

    TChain * mychain = new TChain("RecoAnalyzer/Events");

    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/RecoAnalyzer/RecoStudy.root");

    double Gen_Kp_ETA;
    double Gen_Kp_PHI;
    double Gen_Km_ETA;
    double Gen_Km_PHI;
    double Gen_pi_ETA;
    double Gen_pi_PHI;

    std::vector<double> * match_Kp_ETA = nullptr;
    std::vector<double> * match_Kp_PHI = nullptr;
    std::vector<double> * match_Km_ETA = nullptr;
    std::vector<double> * match_Km_PHI = nullptr;
    std::vector<double> * match_pi_ETA = nullptr;
    std::vector<double> * match_pi_PHI = nullptr;

    std::vector<double> * match_phiFit_Kp_ETA = nullptr;
    std::vector<double> * match_phiFit_Kp_PHI = nullptr;
    std::vector<double> * match_phiFit_Km_ETA = nullptr;
    std::vector<double> * match_phiFit_Km_PHI = nullptr;
    std::vector<double> * match_phiFit_pi_ETA = nullptr;
    std::vector<double> * match_phiFit_pi_PHI = nullptr;

    std::vector<double> * match_DsFit_Kp_ETA = nullptr;
    std::vector<double> * match_DsFit_Kp_PHI = nullptr;
    std::vector<double> * match_DsFit_Km_ETA = nullptr;
    std::vector<double> * match_DsFit_Km_PHI = nullptr;
    std::vector<double> * match_DsFit_pi_ETA = nullptr;
    std::vector<double> * match_DsFit_pi_PHI = nullptr;
    
    mychain->SetBranchAddress("Gen_Kp_ETA", &Gen_Kp_ETA);
    mychain->SetBranchAddress("Gen_Kp_PHI", &Gen_Kp_PHI);
    mychain->SetBranchAddress("Gen_Km_ETA", &Gen_Km_ETA);
    mychain->SetBranchAddress("Gen_Km_PHI", &Gen_Km_PHI);
    mychain->SetBranchAddress("Gen_pi_ETA", &Gen_pi_ETA);
    mychain->SetBranchAddress("Gen_pi_PHI", &Gen_pi_PHI);

    mychain->SetBranchAddress("match_Kp_ETA", &match_Kp_ETA);
    mychain->SetBranchAddress("match_Kp_PHI", &match_Kp_PHI);
    mychain->SetBranchAddress("match_Km_ETA", &match_Km_ETA);
    mychain->SetBranchAddress("match_Km_PHI", &match_Km_PHI);
    mychain->SetBranchAddress("match_pi_ETA", &match_pi_ETA);
    mychain->SetBranchAddress("match_pi_PHI", &match_pi_PHI);

    mychain->SetBranchAddress("match_phiFit_Kp_ETA", &match_phiFit_Kp_ETA);
    mychain->SetBranchAddress("match_phiFit_Kp_PHI", &match_phiFit_Kp_PHI);
    mychain->SetBranchAddress("match_phiFit_Km_ETA", &match_phiFit_Km_ETA);
    mychain->SetBranchAddress("match_phiFit_Km_PHI", &match_phiFit_Km_PHI);
    mychain->SetBranchAddress("match_phiFit_pi_ETA", &match_phiFit_pi_ETA);
    mychain->SetBranchAddress("match_phiFit_pi_PHI", &match_phiFit_pi_PHI);

    mychain->SetBranchAddress("match_DsFit_Kp_ETA", &match_DsFit_Kp_ETA);
    mychain->SetBranchAddress("match_DsFit_Kp_PHI", &match_DsFit_Kp_PHI);
    mychain->SetBranchAddress("match_DsFit_Km_ETA", &match_DsFit_Km_ETA);
    mychain->SetBranchAddress("match_DsFit_Km_PHI", &match_DsFit_Km_PHI);
    mychain->SetBranchAddress("match_DsFit_pi_ETA", &match_DsFit_pi_ETA);
    mychain->SetBranchAddress("match_DsFit_pi_PHI", &match_DsFit_pi_PHI);

    TFile * outfile = new TFile("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/RecoAnalyzer/DeltaR.root", "RECREATE");
	TTree * outtree = new TTree("Events", "Events");

    double Kp_deltaR;
    double Km_deltaR;
    double pi_deltaR;
    double Kp_phiFit_deltaR;
    double Km_phiFit_deltaR;
    double pi_phiFit_deltaR;
    double Kp_DsFit_deltaR;
    double Km_DsFit_deltaR;
    double pi_DsFit_deltaR;
    
    outtree->Branch("Kp_deltaR", &Kp_deltaR);
    outtree->Branch("Km_deltaR", &Km_deltaR);
    outtree->Branch("pi_deltaR", &pi_deltaR);
    outtree->Branch("Kp_phiFit_deltaR", &Kp_phiFit_deltaR);
    outtree->Branch("Km_phiFit_deltaR", &Km_phiFit_deltaR);
    outtree->Branch("pi_phiFit_deltaR", &pi_phiFit_deltaR);
    outtree->Branch("Kp_DsFit_deltaR", &Kp_DsFit_deltaR);
    outtree->Branch("Km_DsFit_deltaR", &Km_DsFit_deltaR);
    outtree->Branch("pi_DsFit_deltaR", &pi_DsFit_deltaR);

    int nentries = mychain->GetEntries();

    for(int i=0; i<nentries; i++){

        match_Kp_ETA->clear();
        match_Kp_PHI->clear();
        match_Km_ETA->clear();
        match_Km_PHI->clear();
        match_pi_ETA->clear();
        match_pi_PHI->clear();

        match_phiFit_Kp_ETA->clear();
        match_phiFit_Kp_PHI->clear();
        match_phiFit_Km_ETA->clear();
        match_phiFit_Km_PHI->clear();
        match_phiFit_pi_ETA->clear();
        match_phiFit_pi_PHI->clear();

        match_DsFit_Kp_ETA->clear();
        match_DsFit_Kp_PHI->clear();
        match_DsFit_Km_ETA->clear();
        match_DsFit_Km_PHI->clear();
        match_DsFit_pi_ETA->clear();
        match_DsFit_pi_PHI->clear();
        
        mychain->GetEntry(i);
        
        if(match_DsFit_Kp_ETA->size() < 1) continue;
        
        Kp_deltaR = getDeltaR(match_Kp_ETA->at(0), match_Kp_PHI->at(0), Gen_Kp_ETA, Gen_Kp_PHI);
        Km_deltaR = getDeltaR(match_Km_ETA->at(0), match_Km_PHI->at(0), Gen_Km_ETA, Gen_Km_PHI);
        pi_deltaR = getDeltaR(match_pi_ETA->at(0), match_pi_PHI->at(0), Gen_pi_ETA, Gen_pi_PHI);
        Kp_phiFit_deltaR = getDeltaR(match_phiFit_Kp_ETA->at(0), match_phiFit_Kp_PHI->at(0), Gen_Kp_ETA, Gen_Kp_PHI);
        Km_phiFit_deltaR = getDeltaR(match_phiFit_Km_ETA->at(0), match_phiFit_Km_PHI->at(0), Gen_Km_ETA, Gen_Km_PHI);
        pi_phiFit_deltaR = getDeltaR(match_phiFit_pi_ETA->at(0), match_phiFit_pi_PHI->at(0), Gen_pi_ETA, Gen_pi_PHI);
        Kp_DsFit_deltaR = getDeltaR(match_DsFit_Kp_ETA->at(0), match_DsFit_Kp_PHI->at(0), Gen_Kp_ETA, Gen_Kp_PHI);
        Km_DsFit_deltaR = getDeltaR(match_DsFit_Km_ETA->at(0), match_DsFit_Km_PHI->at(0), Gen_Km_ETA, Gen_Km_PHI);
        pi_DsFit_deltaR = getDeltaR(match_DsFit_pi_ETA->at(0), match_DsFit_pi_PHI->at(0), Gen_pi_ETA, Gen_pi_PHI);
        
        outtree->Fill();
    }

    outtree->Write();
	delete outtree;
	outfile->Close();
	delete outfile;


    return 0;
}

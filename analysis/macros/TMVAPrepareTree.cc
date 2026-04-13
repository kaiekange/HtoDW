#include <cstdlib>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"

Int_t TMVAPrepareTree()
{
    //========================================//
    //  Load files into corresponding chains  //
    //========================================//
    TChain *SigTree = new TChain("Events");
    SigTree->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_ReDefine.root");

    TChain *BkgTree = new TChain("Events");
    BkgTree->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_ReDefine.root");


    TCut SigCut = "event_match && num_ss_mu_tight == 0 && num_os_mu_tight == 1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125 && TMath::Finite(Ds_FDsig) && TMath::Finite(Km_ipsig) && TMath::Finite(Kp_ipsig) && TMath::Finite(pi_ipsig)";
    TCut BkgCut = "num_ss_mu == 0 && num_os_mu == 1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125 && TMath::Finite(Ds_FDsig) && TMath::Finite(Km_ipsig) && TMath::Finite(Kp_ipsig) && TMath::Finite(pi_ipsig)";

    TRandom3 rng(12345);
    TTree *SigTreeFiltered = SigTree->CopyTree(SigCut);
    TTree *SigTrain = SigTreeFiltered->CloneTree(0);
    TTree *SigVal = SigTreeFiltered->CloneTree(0);
    TTree *SigTest = SigTreeFiltered->CloneTree(0);
    Long64_t nSignal = SigTreeFiltered->GetEntries();
    for (Long64_t i = 0; i < nSignal; i++)
    {
        SigTreeFiltered->GetEntry(i);
        Double_t r = rng.Uniform();
        if (r < 0.5)
            SigTrain->Fill();
        else if (r < 0.75)
            SigVal->Fill();
        else
            SigTest->Fill();
    }

    TFile *SigTrainFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_Train.root", "RECREATE");
    SigTrain->Write();l
    SigTrainFile->Close();
    delete SigTrainFile;

    TFile *SigValFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_Val.root", "RECREATE");
    SigVal->Write();
    SigValFile->Close();
    delete SigValFile;

    TFile *SigTestFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_Test.root", "RECREATE");
    SigTest->Write();
    SigTestFile->Close();
    delete SigTestFile;

    TTree *BkgTreeFiltered = BkgTree->CopyTree(BkgCut);
    TTree *BkgTrain = BkgTreeFiltered->CloneTree(0);
    TTree *BkgVal = BkgTreeFiltered->CloneTree(0);
    TTree *BkgTest = BkgTreeFiltered->CloneTree(0);
    Long64_t nBackground = BkgTreeFiltered->GetEntries();
    for (Long64_t i = 0; i < nBackground; i++)
    {
        BkgTreeFiltered->GetEntry(i);
        Double_t r = rng.Uniform();
        if (r < 0.5)
            BkgTrain->Fill();
        else if (r < 0.75)
            BkgVal->Fill();
        else
            BkgTest->Fill();
    }

    TFile *BkgTrainFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_Train.root", "RECREATE");
    BkgTrain->Write();
    BkgTrainFile->Close();
    delete BkgTrainFile;

    TFile *BkgValFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_Val.root", "RECREATE");
    BkgVal->Write();
    BkgValFile->Close();
    delete BkgValFile;

    TFile *BkgTestFile = new TFile("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_Test.root", "RECREATE");
    BkgTest->Write();
    BkgTestFile->Close();
    delete BkgTestFile;

    return 0;
}
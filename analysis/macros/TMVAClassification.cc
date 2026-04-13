#include <cstdlib>
#include <iostream>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

Int_t TMVAClassification()
{
    TMVA::Tools::Instance();

    std::cout << "\n==> Start TMVAClassification" << std::endl;

    //========================================//
    //  Load files into corresponding chains  //
    //========================================//
    TFile *SigTrainFile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_Train.root");
    TFile *SigValFile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_Val.root");
    TFile *BkgTrainFile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_Train.root");
    TFile *BkgValFile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_Val.root");

    TTree *SigTrain = (TTree *)SigTrainFile->Get("Events");
    TTree *SigVal = (TTree *)SigValFile->Get("Events");
    TTree *BkgTrain = (TTree *)BkgTrainFile->Get("Events");
    TTree *BkgVal = (TTree *)BkgValFile->Get("Events");

    //====================================================//
    //  Add TMVA classification discriminating variables  //
    //====================================================//
    TString TMVApath = "TMVA.root";
    TFile *outfile = TFile::Open(TMVApath, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outfile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

    TString datasetpath = "dataset";
    TMVA::DataLoader *dataloader = new TMVA::DataLoader(datasetpath);

    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    // dataloader->AddVariable( "myvar1 := var1+var2", 'F' );
    // dataloader->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
    // dataloader->AddVariable( "var3",                "Variable 3", "units", 'F' );
    // dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );

    dataloader->AddVariable("Ds_IsoR04_RelIso", "#it{D_{s}^{+}} IsoR04 RelIso", "units", 'F');
    dataloader->AddVariable("log(dira)", "log(cos #alpha(PV,SV))", "units", 'F');
    dataloader->AddVariable("log(Ds_FDsig)", "log(#it{D_{s}^{+}} FD significance)", "units", 'F');
    // dataloader->AddVariable("Kp_ipsig", "#it{K^{+}} IP significance", "units", 'F');
    // dataloader->AddVariable("Km_ipsig", "#it{K^{-}} IP significance", "units", 'F');
    // dataloader->AddVariable("pi_ipsig", "#it{#pi^{+}} IP significance", "units", 'F');
    dataloader->AddVariable("log(min(min(Kp_ipsig, Km_ipsig), pi_ipsig))", "log(min IP significance of #it{K^{+},K^{-},#pi^{+}})", "units", 'F');

    dataloader->AddVariable("K1_pt", "#it{p_{T}(K_{1})}", "units", 'F');
    // dataloader->AddVariable("K2_pt", "#it{p_{T}(K_{2})}", "units", 'F');
    dataloader->AddVariable("K1_pt-K2_pt", "#it{p_{T}(K_{1})-p_{T}(K_{2})}", "units", 'F');
    dataloader->AddVariable("pi_pt", "#it{p_{T}(#pi^{+})}", "units", 'F');
    dataloader->AddVariable("mu_pt", "#it{p_{T}(#mu^{-})}", "units", 'F');
    dataloader->AddVariable("H_pt", "#it{p_{T}(D_{s}^{+}#mu^{-})}", "units", 'F');
    dataloader->AddVariable("phi_pt/Ds_pt", "#it{p_{T}(#phi)/p_{T}(D_{s}^{+})}", "units", 'F');
    dataloader->AddVariable("Ds_pt/H_pt", "#it{p_{T}(D_{s}^{+})/p_{T}(D_{s}^{+}#mu^{-})}", "units", 'F');

    dataloader->AddVariable("dR_Kp_Km", "#it{#DeltaR(K_{1},K_{2})}", "units", 'F');
    // dataloader->AddVariable("dR_K1_pi", "#it{#DeltaR(K_{1},#pi^{+})}", "units", 'F');
    // dataloader->AddVariable("dR_K1_phi", "#it{#DeltaR(K_{1},#phi)}", "units", 'F');
    // dataloader->AddVariable("dR_K1_Ds", "#it{#DeltaR(K_{1},D_{s}^{+})}", "units", 'F');
    // dataloader->AddVariable("dR_K2_pi", "#it{#DeltaR(K_{2},#pi^{+})}", "units", 'F');
    // dataloader->AddVariable("dR_K2_phi", "#it{#DeltaR(K_{2},#phi)}", "units", 'F');
    // dataloader->AddVariable("dR_K2_Ds", "#it{#DeltaR(K_{2},D_{s}^{+})}", "units", 'F');
    dataloader->AddVariable("dR_pi_phi", "#it{#DeltaR(#pi^{+},#phi)}", "units", 'F');
    // dataloader->AddVariable("dR_pi_Ds", "#it{#DeltaR(#pi^{+},D_{s}^{+})}", "units", 'F');
    // dataloader->AddVariable("dR_phi_Ds", "#it{#DeltaR(#phi,D_{s}^{+})}", "units", 'F');

    // dataloader->AddVariable("dR_K1_mu", "#it{#DeltaR(K_{1},#mu^{-})}", "units", 'F');
    // dataloader->AddVariable("dR_K2_mu", "#it{#DeltaR(K_{2},#mu^{-})}", "units", 'F');
    // dataloader->AddVariable("dR_pi_mu", "#it{#DeltaR(#pi^{+},#mu^{-})}", "units", 'F');
    // dataloader->AddVariable("dR_phi_mu", "#it{#DeltaR(#phi,#mu^{-})}", "units", 'F');
    dataloader->AddVariable("dR_Ds_mu", "#it{#DeltaR(D_{s}^{+},#mu^{-})}", "units", 'F');

    dataloader->AddVariable("met_pt", "MET #it{p_{T}}", "units", 'F');
    dataloader->AddVariable("met_sumEt", "MET #sum#it{E_{T}}", "units", 'F');
    dataloader->AddVariable("num_jets_pt30", "# jets (#it{p_{T}}>30 GeV)", "units", 'I');
    dataloader->AddVariable("W_masst", "#it{M_{T}(W^{-})}", "units", 'F');

    // Assign samples weight
    Double_t sigweight = 1.0;
    Double_t bkgweight = 1.0;

    // You can add an arbitrary number of signal or background trees
    // dataloader->AddSignalTree(SigFiltered, sigweight);
    // dataloader->AddBackgroundTree(BkgFiltered, bkgweight);

    dataloader->AddSignalTree(SigTrain, 1.0, TMVA::Types::kTraining);
    dataloader->AddSignalTree(SigVal, 1.0, TMVA::Types::kTesting);

    dataloader->AddBackgroundTree(BkgTrain, 1.0, TMVA::Types::kTraining);
    dataloader->AddBackgroundTree(BkgVal, 1.0, TMVA::Types::kTesting);

    // Tell the dataloader how to use the training and testing events
    //
    // If no numbers of events are given, half of the events in the tree are used
    // for training, and the other half for testing:
    //
    // dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "SplitMode=random:!V" );
    //
    // To also specify the number of testing events, use:
    //
    //    dataloader->PrepareTrainingAndTestTree( mycut,
    //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );

    // dataloader->PrepareTrainingAndTestTree("", "", "SplitMode=Random:NormMode=NumEvents:!V");
    dataloader->PrepareTrainingAndTestTree("", "", "NormMode=NumEvents:!V");
    // dataloader->PrepareTrainingAndTestTree( SigCut, BkgCut, "nTrain_Signal=67223:nTrain_Background=150000:SplitMode=Random:NormMode=NumEvents:!V" );

    // Use BDT
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
    // factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

    // "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

    // For an example of the category classifier usage, see: TMVAClassificationCategory
    //
    // --------------------------------------------------------------------------------------------------
    //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
    // STILL EXPERIMENTAL and only implemented for BDT's !
    //
    //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
    //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
    //
    // --------------------------------------------------------------------------------------------------

    // factory->OptimizeAllMethods("ROCIntegral","Scan");
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outfile->Close();

    std::cout << "==> Wrote root file: " << outfile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch())
        TMVA::TMVAGui(TMVApath);

    return 0;
}
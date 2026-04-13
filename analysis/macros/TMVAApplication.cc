#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStopwatch.h>

#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>

using namespace TMVA;

const TString ChannelSample[14] = {
	"Signal",
	"WJets",
	"W0J",
	"W1J",
	"W2J",
	"TT",
	"TW",
	"WW",
	"WZ",
	"WZto2L2Q",
	"ZZto2L2Q",
	"WWto2L2Nu_OS",
	"WWto2L2Nu_SS_EW",
	"WWto2L2Nu_SS_QCD"};

Int_t TMVAApplication()
{
	// Load library
	TMVA::Tools::Instance();
	std::cout << "\n==> Start TMVAClassificationApplication" << std::endl;

	TMVA::Reader * reader = new TMVA::Reader( "!Color:!Silent" );

	Float_t dl_Ds_IsoR04_RelIso;
	Float_t dl_log_dira;
	Float_t dl_log_Ds_FDsig;
	Float_t dl_min_ipsig;

	Float_t dl_K1_pt;
	Float_t dl_K1_pt_minus_K2_pt;
	Float_t dl_pi_pt;
	Float_t dl_mu_pt;
	Float_t dl_H_pt;
	Float_t dl_phi_pt_div_Ds_pt;
	Float_t dl_Ds_pt_div_H_pt;

	Float_t dl_dR_Kp_Km;
	Float_t dl_dR_pi_phi;
	Float_t dl_dR_Ds_mu;

	Float_t dl_met_pt;
	Float_t dl_met_sumEt;
	Float_t dl_num_jets_pt30;
	Float_t dl_W_masst;

	reader->AddVariable( "Ds_IsoR04_RelIso", &dl_Ds_IsoR04_RelIso );
	reader->AddVariable("log(dira)", &dl_log_dira);
    reader->AddVariable("log(Ds_FDsig)", &dl_log_Ds_FDsig);
    reader->AddVariable("log(min(min(Kp_ipsig, Km_ipsig), pi_ipsig))", &dl_min_ipsig);

    reader->AddVariable("K1_pt", &dl_K1_pt);
    reader->AddVariable("K1_pt-K2_pt", &dl_K1_pt_minus_K2_pt);
    reader->AddVariable("pi_pt", &dl_pi_pt);
    reader->AddVariable("mu_pt", &dl_mu_pt);
    reader->AddVariable("H_pt", &dl_H_pt);
    reader->AddVariable("phi_pt/Ds_pt", &dl_phi_pt_div_Ds_pt);
    reader->AddVariable("Ds_pt/H_pt", &dl_Ds_pt_div_H_pt);
    
    reader->AddVariable("dR_Kp_Km", &dl_dR_Kp_Km);
    reader->AddVariable("dR_pi_phi", &dl_dR_pi_phi);
    reader->AddVariable("dR_Ds_mu", &dl_dR_Ds_mu);

    reader->AddVariable("met_pt", &dl_met_pt);
    reader->AddVariable("met_sumEt", &dl_met_sumEt);
    reader->AddVariable("num_jets_pt30", &dl_num_jets_pt30);
    reader->AddVariable("W_masst", &dl_W_masst);

	//========== Load in classifier ==========//
	TString methodName = "BDT method";
	TString weightfile = "dataset/weights/TMVAClassification_BDT.weights.xml";
	reader->BookMVA( methodName, weightfile );

	for (Int_t cid = 0; cid < 2; cid++)
    {
        TString datafilepath;
		// if(cid < 2)
		// 	datafilepath = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[cid] + "_Test.root";
		// else
			// datafilepath = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[cid] + "_ReDefine.root";
		datafilepath = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[cid] + "_ReDefine.root";
		
		TString TMVAApppath = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[cid] + "_TMVAapp.root";
		//========== Load in sample ==========//
		TChain * mychain = new TChain("Events");
		std::cout << "--- TMVAClassificationApp    : Using input file: " << datafilepath << std::endl;
		mychain->Add(datafilepath);

		//========== Load discriminating variables ==========//
		Float_t Ds_IsoR04_RelIso;
		Float_t dira;
		Float_t Ds_FDsig;
		Float_t Kp_ipsig;
		Float_t Km_ipsig;
		Float_t pi_ipsig;

		Float_t K1_pt;
		Float_t K2_pt;
		Float_t pi_pt;
		Float_t mu_pt;
		Float_t H_pt;
		Float_t phi_pt;
		Float_t Ds_pt;

		Float_t dR_Kp_Km;
		Float_t dR_pi_phi;
		Float_t dR_Ds_mu;

		Float_t met_pt;
		Float_t met_sumEt;
		Int_t num_jets_pt30;
		Float_t W_masst;

		mychain->SetBranchAddress("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso);
		mychain->SetBranchAddress("dira", &dira);
		mychain->SetBranchAddress("Ds_FDsig", &Ds_FDsig);
		mychain->SetBranchAddress("Kp_ipsig", &Kp_ipsig);
		mychain->SetBranchAddress("Km_ipsig", &Km_ipsig);
		mychain->SetBranchAddress("pi_ipsig", &pi_ipsig);

		mychain->SetBranchAddress("K1_pt", &K1_pt);
		mychain->SetBranchAddress("K2_pt", &K2_pt);
		mychain->SetBranchAddress("pi_pt", &pi_pt);
		mychain->SetBranchAddress("mu_pt", &mu_pt);
		mychain->SetBranchAddress("H_pt", &H_pt);
		mychain->SetBranchAddress("phi_pt", &phi_pt);
		mychain->SetBranchAddress("Ds_pt", &Ds_pt);

		mychain->SetBranchAddress("dR_Kp_Km", &dR_Kp_Km);
		mychain->SetBranchAddress("dR_pi_phi", &dR_pi_phi);
		mychain->SetBranchAddress("dR_Ds_mu", &dR_Ds_mu);

		mychain->SetBranchAddress("met_pt", &met_pt);
		mychain->SetBranchAddress("met_sumEt", &met_sumEt);
		mychain->SetBranchAddress("num_jets_pt30", &num_jets_pt30);
		mychain->SetBranchAddress("W_masst", &W_masst);

		Int_t nentries = mychain->GetEntries();
		std::cout << "--- Processing: " << nentries << " events" << std::endl;
		TStopwatch sw;
		sw.Start();

		TFile *outfile  = new TFile( TMVAApppath, "RECREATE" );
		TTree *outtree= mychain->CloneTree(0);
		Float_t BDT;
		outtree->Branch("BDT", &BDT, "BDT/F");

		for (Long64_t i=0; i<nentries; i++) {
			BDT =-100;
			if (i%10000 == 0) std::cout << "--- ... Processing event: " << i << std::endl;
			mychain->GetEntry(i);

			dl_Ds_IsoR04_RelIso = Ds_IsoR04_RelIso;
			dl_log_dira = log(dira);
			dl_log_Ds_FDsig = (Ds_FDsig > 0) ? log(Ds_FDsig) : -10;
			float kp = std::isfinite(Kp_ipsig) ? Kp_ipsig : 1e6;
			float km = std::isfinite(Km_ipsig) ? Km_ipsig : 1e6;
			float pi = std::isfinite(pi_ipsig) ? pi_ipsig : 1e6;

			dl_min_ipsig = log(std::min(std::min(kp, km), pi));

			// dl_min_ipsig = std::min(std::min(Kp_ipsig, Km_ipsig), pi_ipsig);

			dl_K1_pt = K1_pt;
			dl_K1_pt_minus_K2_pt = K1_pt - K2_pt;
			dl_pi_pt = pi_pt;
			dl_mu_pt = mu_pt;
			dl_H_pt = H_pt;
			dl_phi_pt_div_Ds_pt = phi_pt / Ds_pt;
			dl_Ds_pt_div_H_pt = Ds_pt / H_pt;

			dl_dR_Kp_Km = dR_Kp_Km;
			dl_dR_pi_phi = dR_pi_phi;
			dl_dR_Ds_mu = dR_Ds_mu;

			dl_met_pt = met_pt;
			dl_met_sumEt = met_sumEt;
			dl_num_jets_pt30 = num_jets_pt30;
			dl_W_masst = W_masst;

			BDT =reader->EvaluateMVA("BDT method");
			outtree->Fill();
		}

		outtree->Write();
		sw.Stop();
		std::cout << "--- End of event loop: "; sw.Print();
		outfile->Close();
		std::cout << "--- Created root file: " << TMVAApppath << "containing the MVA output histograms" << std::endl;
	}

	delete reader;
	std::cout << "==> TMVAClassificationApplication is done!\n" << std::endl;
	return 0;
}
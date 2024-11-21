#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCut.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLatex.h>

const double myfontsize = 0.05;

int draw(){

	//================//
	// Use LHCb Style //
	//================//
	TH1::SetDefaultSumw2();
    #include "lhcbStyle.C"
    gROOT->ProcessLine(".x lhcbStyle.C");
    // Initialize chains
	TString treename = "Events";
	TChain * mychain = new TChain(treename); 
    for(Int_t i=0; i<100; i++ ){
        TString inpath = Form("../20241113/samples/process_%d/output/S5_PAT.root", i);
	    mychain->Add(inpath);
    }
    std::cout << mychain->GetEntries() << std::endl;
	// Initialize canvas
	TCanvas* c1 = new TCanvas("c1","",800,600);
	gStyle->SetOptStat(0);
	gPad->SetLeftMargin(0.2);

	Double_t varmax = 150;
	Double_t varmin = 0;
	Int_t nbins = 20;
	Double_t binsize = (varmax - varmin) / Double_t(nbins);
	Int_t precision = std::max(0, -Int_t(std::log10(binsize)) + 2);
	TString variable = "patPackedGenParticles_packedGenParticles__PAT.obj.packedPt_";
	TString partID = "25";
	TString myunit = "GeV/#it{c}";
	TString xtitle = "#it{P_{t}(H)} [" + myunit + "]";
	TString ytitleprefix = Form("# candidates/ (%.*f", precision, binsize);
    TString figdirpath = "../20241113/figures/";

	TH1F * h1 = new TH1F("h1", "h1", nbins, varmin, varmax);

    TCut  mycut;
    mycut += "abs(patPackedGenParticles_packedGenParticles__PAT.obj.pdgId_) == " + partID;
    //mycut += "recoGenParticles_prunedGenParticles__PAT.obj.m_state.status_ == 22";

    mychain->Draw(variable + ">> h1", mycut);

	h1->SetLineColor(1);
	h1->SetMarkerColor(1);
	h1->SetMarkerSize(0.5);
	h1->SetMinimum(0.0);
	h1->GetXaxis()->SetTitle(xtitle);
	h1->GetYaxis()->SetTitle(ytitleprefix + " " + myunit + ")");
	h1->GetXaxis()->SetTitleOffset(1.3);
	h1->GetYaxis()->SetTitleOffset(1.5);
	h1->GetXaxis()->CenterTitle(true);
	h1->GetXaxis()->SetTitleSize(myfontsize);
    h1->GetXaxis()->SetLabelSize(myfontsize); 
	h1->GetYaxis()->CenterTitle(true);
	h1->GetYaxis()->SetTitleSize(myfontsize);
    h1->GetYaxis()->SetLabelSize(myfontsize); 

	// Draw histograms
	h1->Draw("ep");
	// Save figure and clean meomry
	c1->Update();
	TString outplotpath = figdirpath + "plot_" + partID + "_Pt" + ".png";
	c1->SaveAs(outplotpath);
	
	return 0;
}

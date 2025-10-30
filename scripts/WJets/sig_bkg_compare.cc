#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void compare(TChain *mychain1, TChain *mychain2, TString myvar1, TString myvar2, TString leg1, TString leg2, TString vartitle, float varmin, float varmax, TCut cut1, TCut cut2, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", 100, varmin, varmax);
    mychain1->Project("h1", myvar1, cut1);
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlack);
    h1->SetMarkerColor(kBlack);
    h1->SetMarkerSize(0.7);

    TH1F *h2 = new TH1F("h2", "", 100, varmin, varmax);
    mychain2->Project("h2", myvar2, cut2);
    h2->Scale(1./h2->Integral());
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
    h2->SetMarkerSize(0.7);

    float height = std::max(h1->GetMaximum(), h2->GetMaximum());

    h1->SetMaximum(height*1.3);
    h1->SetMinimum(0);
    h1->GetXaxis()->SetTitle(vartitle);
    h1->GetYaxis()->SetTitle("Normalized # candidates");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    h1->Draw("ep");
    h2->Draw("ep same");

    TLegend * mylegend = new TLegend(0.7, 0.75, 0.9, 0.9);
    mylegend->AddEntry(h1, leg1, "l");
    mylegend->AddEntry(h2, leg2, "l");
    mylegend->SetFillColor(0);
    mylegend->SetLineWidth(0);
    mylegend->SetTextSize(0.04);
    mylegend->Draw();

    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/sig_bkg_compare/"+figpath+".png");
    delete h1;
    delete h2;
}

int sig_bkg_compare() {
    setTDRStyle();

    TChain *mychain1 = new TChain("PVStudy/Events");
    TChain *mychain2 = new TChain("RecoAnalyzer/Events");
    mychain1->Add("../../tuples/PVStudy_noAP.root");
    mychain2->Add("../../tuples/WJetsToLNu.root");
  
    compare(mychain1, mychain2, "best_Ds_Iso_R0p3", "Ds_Iso_R0p3", "signal", "W+jets", "IsoR03", 0, 200, "1", "1", "Ds_Iso_R0p3");
    compare(mychain1, mychain2, "best_Ds_Iso_R0p3", "Ds_Iso_R0p3", "signal", "W+jets", "IsoR03", 0, 80, "1", "1", "Ds_Iso_R0p3_small");
    
    compare(mychain1, mychain2, "best_Ds_Iso_R0p4", "Ds_Iso_R0p4", "signal", "W+jets", "IsoR04", 0, 200, "1", "1", "Ds_Iso_R0p4");
    compare(mychain1, mychain2, "best_Ds_Iso_R0p4", "Ds_Iso_R0p4", "signal", "W+jets", "IsoR04", 0, 80, "1", "1", "Ds_Iso_R0p4_small");

    compare(mychain1, mychain2, "best_Ds_IsoRel_R0p3", "Ds_IsoRel_R0p3", "signal", "W+jets", "IsoR03 / #it{p_{T}}(#it{D_{s}^{+}})", 0, 20, "1", "1", "Ds_IsoRel_R0p3");
    compare(mychain1, mychain2, "best_Ds_IsoRel_R0p3", "Ds_IsoRel_R0p3", "signal", "W+jets", "IsoR03 / #it{p_{T}}(#it{D_{s}^{+}})", 0, 5, "1", "1", "Ds_IsoRel_R0p3_small");
    
    compare(mychain1, mychain2, "best_Ds_IsoRel_R0p4", "Ds_IsoRel_R0p4", "signal", "W+jets", "IsoR04 / #it{p_{T}}(#it{D_{s}^{+}})", 0, 20, "1", "1", "Ds_IsoRel_R0p4");
    compare(mychain1, mychain2, "best_Ds_IsoRel_R0p4", "Ds_IsoRel_R0p4", "signal", "W+jets", "IsoR04 / #it{p_{T}}(#it{D_{s}^{+}})", 0, 5, "1", "1", "Ds_IsoRel_R0p4_small");
    
    delete mychain1;
    delete mychain2;

    return 0;
}

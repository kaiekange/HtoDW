#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TCut mycut, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", vartitle, 100, varmin, varmax);

    mychain->Project("hist", myvar, mycut);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw("ep");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs(figpath);
    delete hist;
}

void draw_2d(TChain *mychain, TString myvar, TString vartitle, float xvarmin, float xvarmax, float yvarmin, float yvarmax, TCut mycut, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH2F *hist = new TH2F("hist", vartitle, 100, xvarmin, xvarmax, 100, yvarmin, yvarmax);

    mychain->Project("hist", myvar, mycut);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    /* hist->SetMaximum(1.3*hist->GetMaximum()); */
    hist->SetMinimum(0);
    hist->Draw("COLZ");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs(figpath);
    delete hist;
}

void compare(TChain *mychain, TString myvar1, TString myvar2, TString leg1, TString leg2, TString vartitle, float varmin, float varmax, TCut cut1, TCut cut2, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", 100, varmin, varmax);
    mychain->Project("h1", myvar1, cut1);
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlue);
    h1->SetMarkerColor(kBlue);
    h1->SetMarkerSize(0.7);


    TH1F *h2 = new TH1F("h2", "", 100, varmin, varmax);
    mychain->Project("h2", myvar2, cut2);
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
    canvas->SaveAs(figpath);
    delete h1;
    delete h2;
}


int PV_plots() {
    setTDRStyle();

    TChain *mychain = new TChain("PVStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PVStudy/PVStudy.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/PV_figures")) gSystem->MakeDirectory("./figures/PV_figures");

    compare(mychain, "Gen_H_ORIVX_X", "PV_X", "Gen", "Fit", "#it{x}_{PV} [cm]", -0.1, 0.1, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Gen_Fit_PV_X.png");
    compare(mychain, "Gen_H_ORIVX_Y", "PV_Y", "Gen", "Fit", "#it{y}_{PV} [cm]", -0.1, 0.2, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Gen_Fit_PV_Y.png");
    compare(mychain, "Gen_H_ORIVX_Z", "PV_Z", "Gen", "Fit", "#it{z}_{PV} [cm]", -15, 15, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Gen_Fit_PV_Z.png");

    compare(mychain, "Gen_Kp_ORIVX_X", "PV_X", "Gen #it{D_{s}^{+}} decay", "Fit PV", "#it{x} [cm]", -0.1, 0.1, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Kp_PV_X.png");
    compare(mychain, "Gen_Kp_ORIVX_Y", "PV_Y", "Gen #it{D_{s}^{+}} decay", "Fit PV", "#it{y} [cm]", -0.1, 0.2, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Kp_PV_Y.png");
    compare(mychain, "Gen_Kp_ORIVX_Z", "PV_Z", "Gen #it{D_{s}^{+}} decay", "Fit PV", "#it{z} [cm]", -15, 15, "1", "PV_IsValid && !PV_IsFake", "./figures/PV_figures/Kp_PV_Z.png");

    /* draw_1d(mychain, "(Gen_dR_Km_pi - match_dR_Km_pi)", ";Gen #Delta#it{R(K^{-},#pi^{+})}-reco #Delta#it{R(K^{-},#pi^{+})};# candidates", -0.8, 0.8, "1", "./figures/gen_match_compare/diff_dR_Km_pi.png"); */
    /* draw_2d(mychain, "Gen_Kp_PP:((Gen_Kp_PL - Gen_Km_PL)/(Gen_Kp_PL + Gen_Km_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "1", "./figures/APplot/Gen_APplot_phi.png"); */


    delete mychain;

    return 0;
}

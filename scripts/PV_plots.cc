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
    if(gSystem->AccessPathName("./figures/PVStudy")) gSystem->MakeDirectory("./figures/PVStudy");
    
    // compare gen and match 
    /* compare(mychain, "Gen_H_ORIVX_X", "match_PV_PT1_X", "Gen", "match", "#it{x}_{PV} [cm]", -0.04, -0.01, "1", "1", "./figures/PVStudy/compare_PVX_Gen_PT1.png"); */
    /* compare(mychain, "Gen_H_ORIVX_X", "match_PV_PT2_X", "Gen", "match", "#it{x}_{PV} [cm]", -0.04, -0.01, "1", "1", "./figures/PVStudy/compare_PVX_Gen_PT2.png"); */
    /* compare(mychain, "Gen_H_ORIVX_X", "match_PV_PT5_X", "Gen", "match", "#it{x}_{PV} [cm]", -0.04, -0.01, "1", "1", "./figures/PVStudy/compare_PVX_Gen_PT5.png"); */
    /* compare(mychain, "Gen_H_ORIVX_X", "match_PV_PT10_X", "Gen", "match", "#it{x}_{PV} [cm]", -0.04, -0.01, "1", "1", "./figures/PVStudy/compare_PVX_Gen_PT10.png"); */

    /* compare(mychain, "Gen_H_ORIVX_Y", "match_PV_PT1_Y", "Gen", "match", "#it{y}_{PV} [cm]", 0.058, 0.08, "1", "1", "./figures/PVStudy/compare_PVY_Gen_PT1.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Y", "match_PV_PT2_Y", "Gen", "match", "#it{y}_{PV} [cm]", 0.058, 0.08, "1", "1", "./figures/PVStudy/compare_PVY_Gen_PT2.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Y", "match_PV_PT5_Y", "Gen", "match", "#it{y}_{PV} [cm]", 0.058, 0.08, "1", "1", "./figures/PVStudy/compare_PVY_Gen_PT5.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Y", "match_PV_PT10_Y", "Gen", "match", "#it{y}_{PV} [cm]", 0.058, 0.08, "1", "1", "./figures/PVStudy/compare_PVY_Gen_PT10.png"); */

    /* compare(mychain, "Gen_H_ORIVX_Z", "match_PV_PT1_Z", "Gen", "match", "#it{z}_{PV} [cm]", -20, 20, "1", "1", "./figures/PVStudy/compare_PVZ_Gen_PT1.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Z", "match_PV_PT2_Z", "Gen", "match", "#it{z}_{PV} [cm]", -20, 20, "1", "1", "./figures/PVStudy/compare_PVZ_Gen_PT2.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Z", "match_PV_PT5_Z", "Gen", "match", "#it{z}_{PV} [cm]", -20, 20, "1", "1", "./figures/PVStudy/compare_PVZ_Gen_PT5.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Z", "match_PV_PT10_Z", "Gen", "match", "#it{z}_{PV} [cm]", -20, 20, "1", "1", "./figures/PVStudy/compare_PVZ_Gen_PT10.png"); */

    /* draw_1d(mychain, "match_PV_PT1_X - Gen_H_ORIVX_X", ";fit - Gen #it{x}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVX_PT1.png"); */
    /* draw_1d(mychain, "match_PV_PT2_X - Gen_H_ORIVX_X", ";fit - Gen #it{x}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVX_PT2.png"); */
    /* draw_1d(mychain, "match_PV_PT5_X - Gen_H_ORIVX_X", ";fit - Gen #it{x}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVX_PT5.png"); */
    /* draw_1d(mychain, "match_PV_PT10_X - Gen_H_ORIVX_X", ";fit - Gen #it{x}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVX_PT10.png"); */

    /* draw_1d(mychain, "match_PV_PT1_Y - Gen_H_ORIVX_Y", ";fit - Gen #it{y}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVY_PT1.png"); */
    /* draw_1d(mychain, "match_PV_PT2_Y - Gen_H_ORIVX_Y", ";fit - Gen #it{y}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVY_PT2.png"); */
    /* draw_1d(mychain, "match_PV_PT5_Y - Gen_H_ORIVX_Y", ";fit - Gen #it{y}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVY_PT5.png"); */
    /* draw_1d(mychain, "match_PV_PT10_Y - Gen_H_ORIVX_Y", ";fit - Gen #it{y}_{PV} [cm];# candidates", -1, 1, "1", "./figures/PVStudy/diff_PVY_PT10.png"); */

    /* draw_1d(mychain, "match_PV_PT1_Z - Gen_H_ORIVX_Z", ";fit - Gen #it{z}_{PV} [cm];# candidates", -15, 15, "1", "./figures/PVStudy/diff_PVZ_PT1.png"); */
    /* draw_1d(mychain, "match_PV_PT2_Z - Gen_H_ORIVX_Z", ";fit - Gen #it{z}_{PV} [cm];# candidates", -15, 15, "1", "./figures/PVStudy/diff_PVZ_PT2.png"); */
    /* draw_1d(mychain, "match_PV_PT5_Z - Gen_H_ORIVX_Z", ";fit - Gen #it{z}_{PV} [cm];# candidates", -15, 15, "1", "./figures/PVStudy/diff_PVZ_PT5.png"); */
    /* draw_1d(mychain, "match_PV_PT10_Z - Gen_H_ORIVX_Z", ";fit - Gen #it{z}_{PV} [cm];# candidates", -15, 15, "1", "./figures/PVStudy/diff_PVZ_PT10.png"); */
    
    /* draw_2d(mychain, "match_PV_PT1_XERR:(match_PV_PT1_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{x}_{PV} error [cm]", -1, 1, 0, 0.0015, "1", "./figures/PVStudy/err_PVX_PT1.png"); */
    /* draw_2d(mychain, "match_PV_PT2_XERR:(match_PV_PT2_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{x}_{PV} error [cm]", -1, 1, 0, 0.0025, "1", "./figures/PVStudy/err_PVX_PT2.png"); */
    /* draw_2d(mychain, "match_PV_PT5_XERR:(match_PV_PT5_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{x}_{PV} error [cm]", -1, 1, 0, 0.01, "1", "./figures/PVStudy/err_PVX_PT5.png"); */
    /* draw_2d(mychain, "match_PV_PT10_XERR:(match_PV_PT10_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{x}_{PV} error [cm]", -1, 1, 0, 0.01, "1", "./figures/PVStudy/err_PVX_PT10.png"); */
    
    /* draw_2d(mychain, "match_PV_PT1_YERR:(match_PV_PT1_Y - Gen_H_ORIVX_Y)", ";fit - Gen #it{y}_{PV} [cm];fit #it{y}_{PV} error [cm]", -1, 1, 0, 0.0015, "1", "./figures/PVStudy/err_PVY_PT1.png"); */
    /* draw_2d(mychain, "match_PV_PT2_YERR:(match_PV_PT2_Y - Gen_H_ORIVX_Y)", ";fit - Gen #it{y}_{PV} [cm];fit #it{y}_{PV} error [cm]", -1, 1, 0, 0.0025, "1", "./figures/PVStudy/err_PVY_PT2.png"); */
    /* draw_2d(mychain, "match_PV_PT5_YERR:(match_PV_PT5_Y - Gen_H_ORIVX_Y)", ";fit - Gen #it{y}_{PV} [cm];fit #it{y}_{PV} error [cm]", -1, 1, 0, 0.01, "1", "./figures/PVStudy/err_PVY_PT5.png"); */
    /* draw_2d(mychain, "match_PV_PT10_YERR:(match_PV_PT10_Y - Gen_H_ORIVX_Y)", ";fit - Gen #it{y}_{PV} [cm];fit #it{y}_{PV} error [cm]", -1, 1, 0, 0.01, "1", "./figures/PVStudy/err_PVY_PT10.png"); */
    
    /* draw_2d(mychain, "match_PV_PT1_ZERR:(match_PV_PT1_Z - Gen_H_ORIVX_Z)", ";fit - Gen #it{z}_{PV} [cm];fit #it{z}_{PV} error [cm]", -1, 1, 0, 0.003, "1", "./figures/PVStudy/err_PVZ_PT1.png"); */
    /* draw_2d(mychain, "match_PV_PT2_ZERR:(match_PV_PT2_Z - Gen_H_ORIVX_Z)", ";fit - Gen #it{z}_{PV} [cm];fit #it{z}_{PV} error [cm]", -1, 1, 0, 0.005, "1", "./figures/PVStudy/err_PVZ_PT2.png"); */
    /* draw_2d(mychain, "match_PV_PT5_ZERR:(match_PV_PT5_Z - Gen_H_ORIVX_Z)", ";fit - Gen #it{z}_{PV} [cm];fit #it{z}_{PV} error [cm]", -1, 1, 0, 0.02, "1", "./figures/PVStudy/err_PVZ_PT5.png"); */
    /* draw_2d(mychain, "match_PV_PT10_ZERR:(match_PV_PT10_Z - Gen_H_ORIVX_Z)", ";fit - Gen #it{z}_{PV} [cm];fit #it{z}_{PV} error [cm]", -1, 1, 0, 0.02, "1", "./figures/PVStudy/err_PVZ_PT10.png"); */

    /* draw_2d(mychain, "match_PV_PT1_CHI2NDOF:(match_PV_PT1_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{#chi}^{2}/NDOF", -1, 1, 0, 400000, "1", "./figures/PVStudy/CHI2NDOF_PV_PT1.png"); */
    /* draw_2d(mychain, "match_PV_PT2_CHI2NDOF:(match_PV_PT2_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{#chi}^{2}/NDOF", -1, 1, 0, 700000, "1", "./figures/PVStudy/CHI2NDOF_PV_PT2.png"); */
    /* draw_2d(mychain, "match_PV_PT5_CHI2NDOF:(match_PV_PT5_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{#chi}^{2}/NDOF", -0.4, 0.4, 0, 10, "1", "./figures/PVStudy/CHI2NDOF_PV_PT5.png"); */
    /* draw_2d(mychain, "match_PV_PT10_CHI2NDOF:(match_PV_PT10_X - Gen_H_ORIVX_X)", ";fit - Gen #it{x}_{PV} [cm];fit #it{#chi}^{2}/NDOF", -0.2, 0.2, 0, 10, "1", "./figures/PVStudy/CHI2NDOF_PV_PT10.png"); */

/*     draw_1d(mychain, "match_PV_PT1_CHI2NDOF", ";fit #it{#chi}^{2}/NDOF;# candidates", 0, 400000, "1", "./figures/PVStudy/CHI2NDOF_PT1.png"); */
/*     draw_1d(mychain, "match_PV_PT2_CHI2NDOF", ";fit #it{#chi}^{2}/NDOF;# candidates", 0, 700000, "1", "./figures/PVStudy/CHI2NDOF_PT2.png"); */
/*     draw_1d(mychain, "match_PV_PT5_CHI2NDOF", ";fit #it{#chi}^{2}/NDOF;# candidates", 0, 10, "1", "./figures/PVStudy/CHI2NDOF_PT5.png"); */
/*     draw_1d(mychain, "match_PV_PT10_CHI2NDOF", ";fit #it{#chi}^{2}/NDOF;# candidates", 0, 10, "1", "./figures/PVStudy/CHI2NDOF_PT10.png"); */

    draw_1d(mychain, "match_PV_PT1_XERR", ";fit #it{x}_{PV} error [cm];# candidates", 0, 0.0015, "1", "./figures/PVStudy/XERR_PT1.png");
    draw_1d(mychain, "match_PV_PT2_XERR", ";fit #it{x}_{PV} error [cm];# candidates", 0, 0.0025, "1", "./figures/PVStudy/XERR_PT2.png");
    draw_1d(mychain, "match_PV_PT5_XERR", ";fit #it{x}_{PV} error [cm];# candidates", 0, 0.01, "1", "./figures/PVStudy/XERR_PT5.png");
    draw_1d(mychain, "match_PV_PT10_XERR", ";fit #it{x}_{PV} error [cm];# candidates", 0, 0.01, "1", "./figures/PVStudy/XERR_PT10.png");

    draw_1d(mychain, "match_PV_PT1_YERR", ";fit #it{y}_{PV} error [cm];# candidates", 0, 0.0015, "1", "./figures/PVStudy/YERR_PT1.png");
    draw_1d(mychain, "match_PV_PT2_YERR", ";fit #it{y}_{PV} error [cm];# candidates", 0, 0.0025, "1", "./figures/PVStudy/YERR_PT2.png");
    draw_1d(mychain, "match_PV_PT5_YERR", ";fit #it{y}_{PV} error [cm];# candidates", 0, 0.01, "1", "./figures/PVStudy/YERR_PT5.png");
    draw_1d(mychain, "match_PV_PT10_YERR", ";fit #it{y}_{PV} error [cm];# candidates", 0, 0.01, "1", "./figures/PVStudy/YERR_PT10.png");

    draw_1d(mychain, "match_PV_PT1_ZERR", ";fit #it{z}_{PV} error [cm];# candidates", 0, 0.003, "1", "./figures/PVStudy/ZERR_PT1.png");
    draw_1d(mychain, "match_PV_PT2_ZERR", ";fit #it{z}_{PV} error [cm];# candidates", 0, 0.005, "1", "./figures/PVStudy/ZERR_PT2.png");
    draw_1d(mychain, "match_PV_PT5_ZERR", ";fit #it{z}_{PV} error [cm];# candidates", 0, 0.02, "1", "./figures/PVStudy/ZERR_PT5.png");
    draw_1d(mychain, "match_PV_PT10_ZERR", ";fit #it{z}_{PV} error [cm];# candidates", 0, 0.02, "1", "./figures/PVStudy/ZERR_PT10.png");
    delete mychain;

    return 0;
}

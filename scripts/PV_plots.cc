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
    h1->GetYaxis()->SetTitle("Normalized # Events");

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

    if(gSystem->AccessPathName("./PVfigures")) gSystem->MakeDirectory("./PVfigures");
    if(gSystem->AccessPathName("./PVfigures/noBS")) gSystem->MakeDirectory("./PVfigures/noBS");
    if(gSystem->AccessPathName("./PVfigures/withBS")) gSystem->MakeDirectory("./PVfigures/withBS");

    TCut sel_withBS = "PV_withBS_IsValid && !PV_withBS_IsFake";
    TCut sel_noBS = "PV_noBS_IsValid && !PV_noBS_IsFake";

    /* compare(mychain, "Gen_H_ORIVX_X", "PV_noBS_X", "Gen", "Fit noBS", "#it{x}_{PV} [cm]", -0.035, -0.015, "1", sel_noBS, "./PVfigures/noBS/compare_PV_X.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Y", "PV_noBS_Y", "Gen", "Fit noBS", "#it{y}_{PV} [cm]", 0.06, 0.08, "1", sel_noBS, "./PVfigures/noBS/compare_PV_Y.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Z", "PV_noBS_Z", "Gen", "Fit noBS", "#it{z}_{PV} [cm]", -15, 15, "1", sel_noBS, "./PVfigures/noBS/compare_PV_Z.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_X - PV_noBS_X)", ";noBS #Delta#it{x}_{PV} [cm];# events", -0.01, 0.01, sel_noBS, "./PVfigures/noBS/diff_PV_X.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Y - PV_noBS_Y)", ";noBS #Delta#it{y}_{PV} [cm];# events", -0.01, 0.01, sel_noBS, "./PVfigures/noBS/diff_PV_Y.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Z - PV_noBS_Z)", ";noBS #Delta#it{z}_{PV} [cm];# events", -0.01, 0.01, sel_noBS, "./PVfigures/noBS/diff_PV_Z.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_X - PV_noBS_X)/Gen_H_ORIVX_X", ";noBS #it{#varepsilonx}_{PV};# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_PV_X.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Y - PV_noBS_Y)/Gen_H_ORIVX_Y", ";noBS #it{#varepsilony}_{PV};# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_PV_Y.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Z - PV_noBS_Z)/Gen_H_ORIVX_Z", ";noBS #it{#varepsilonz}_{PV};# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_PV_Z.png"); */
    /* compare(mychain, "sqrt(Gen_H_ORIVX_X*Gen_H_ORIVX_X + Gen_H_ORIVX_Y*Gen_H_ORIVX_Y)", "sqrt(PV_noBS_X*PV_noBS_X + PV_noBS_Y*PV_noBS_Y)", "Gen", "Fit noBS", "#it{d_{xy}}(PV) [cm]", 0.05, 0.1, "1", sel_noBS, "./PVfigures/noBS/compare_dxy.png"); */
    /* draw_1d(mychain, "(sqrt(Gen_H_ORIVX_X*Gen_H_ORIVX_X + Gen_H_ORIVX_Y*Gen_H_ORIVX_Y) - sqrt(PV_noBS_X*PV_noBS_X + PV_noBS_Y*PV_noBS_Y))", ";noBS #Delta#it{d_{xy}}(PV) [cm];# events", -0.01, 0.01, sel_noBS, "./PVfigures/noBS/diff_dxy.png"); */

    /* compare(mychain, "Gen_H_ORIVX_X", "PV_withBS_X", "Gen", "Fit withBS", "#it{x}_{PV} [cm]", -0.035, -0.015, "1", sel_withBS, "./PVfigures/withBS/compare_PV_X.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Y", "PV_withBS_Y", "Gen", "Fit withBS", "#it{y}_{PV} [cm]", 0.06, 0.08, "1", sel_withBS, "./PVfigures/withBS/compare_PV_Y.png"); */
    /* compare(mychain, "Gen_H_ORIVX_Z", "PV_withBS_Z", "Gen", "Fit withBS", "#it{z}_{PV} [cm]", -15, 15, "1", sel_withBS, "./PVfigures/withBS/compare_PV_Z.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_X - PV_withBS_X)", ";withBS #Delta#it{x}_{PV} [cm];# events", -0.01, 0.01, sel_withBS, "./PVfigures/withBS/diff_PV_X.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Y - PV_withBS_Y)", ";withBS #Delta#it{y}_{PV} [cm];# events", -0.01, 0.01, sel_withBS, "./PVfigures/withBS/diff_PV_Y.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Z - PV_withBS_Z)", ";withBS #Delta#it{z}_{PV} [cm];# events", -0.01, 0.01, sel_withBS, "./PVfigures/withBS/diff_PV_Z.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_X - PV_withBS_X)/Gen_H_ORIVX_X", ";withBS #it{#varepsilonx}_{PV};# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_PV_X.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Y - PV_withBS_Y)/Gen_H_ORIVX_Y", ";withBS #it{#varepsilony}_{PV};# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_PV_Y.png"); */
    /* draw_1d(mychain, "(Gen_H_ORIVX_Z - PV_withBS_Z)/Gen_H_ORIVX_Z", ";withBS #it{#varepsilonz}_{PV};# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_PV_Z.png"); */
    /* compare(mychain, "sqrt(Gen_H_ORIVX_X*Gen_H_ORIVX_X + Gen_H_ORIVX_Y*Gen_H_ORIVX_Y)", "sqrt(PV_withBS_X*PV_withBS_X + PV_withBS_Y*PV_withBS_Y)", "Gen", "Fit withBS", "#it{d_{xy}}(PV) [cm]", 0.05, 0.1, "1", sel_withBS, "./PVfigures/withBS/compare_dxy.png"); */
    /* draw_1d(mychain, "(sqrt(Gen_H_ORIVX_X*Gen_H_ORIVX_X + Gen_H_ORIVX_Y*Gen_H_ORIVX_Y) - sqrt(PV_withBS_X*PV_withBS_X + PV_withBS_Y*PV_withBS_Y))", ";withBS #Delta#it{d_{xy}}(PV) [cm];# events", -0.01, 0.01, sel_withBS, "./PVfigures/withBS/diff_dxy.png"); */

    /* compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2))", "sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2))", "Gen", "Fit noBS", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_noBS, "./PVfigures/noBS/compare_FDxy_Ds_PV.png"); */
    /* compare(mychain, "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z)", "abs(match_DsFit_ENDVX_Z-PV_noBS_Z)", "Gen", "Fit noBS", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_noBS, "./PVfigures/noBS/compare_FDz_Ds_PV.png"); */
    /* compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2))", "sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_noBS_Z,2))", "Gen", "Fit noBS", "FD(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_noBS, "./PVfigures/noBS/compare_FD_Ds_PV.png"); */
    
    /* draw_1d(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2))", ";noBS #DeltaFD_{#it{xy}}(#it{D_{s}^{+}}) [cm];# events", -0.2, 0.2, sel_noBS, "./PVfigures/noBS/diff_FDxy_Ds_PV.png"); */
    /* draw_1d(mychain, "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_noBS_Z)", ";noBS #DeltaFD_{#it{z}}(#it{D_{s}^{+}}) [cm];# events", -0.2, 0.2, sel_noBS, "./PVfigures/noBS/diff_FDz_Ds_PV.png"); */
    /* draw_1d(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_noBS_Z,2))", ";noBS #DeltaFD(#it{D_{s}^{+}}) [cm];# events", -0.2, 0.2, sel_noBS, "./PVfigures/noBS/diff_FD_Ds_PV.png"); */
    
    /* draw_1d(mychain, "(sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2)))/sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2))", ";noBS #it{#varepsilon}FD_{#it{xy}}(#it{D_{s}^{+}});# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_FDxy_Ds_PV.png"); */
    /* draw_1d(mychain, "(abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_noBS_Z))/abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z)", ";noBS #it{#varepsilon}FD_{#it{z}}(#it{D_{s}^{+}});# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_FDz_Ds_PV.png"); */
    /* draw_1d(mychain, "(sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_noBS_Z,2)))/sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2))", ";noBS #it{#varepsilon}FD(#it{D_{s}^{+}});# events", -1, 1, sel_noBS, "./PVfigures/noBS/epsilon_FD_Ds_PV.png"); */

    /* compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2))", "sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2))", "Gen", "Fit withBS", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_withBS, "./PVfigures/withBS/compare_FDxy_Ds_PV.png"); */
    /* compare(mychain, "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z)", "abs(match_DsFit_ENDVX_Z-PV_withBS_Z)", "Gen", "Fit withBS", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_withBS, "./PVfigures/withBS/compare_FDz_Ds_PV.png"); */
    /* compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2))", "sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_withBS_Z,2))", "Gen", "Fit withBS", "FD(#it{D_{s}^{+}}) [cm]", 0, 2, "1", sel_withBS, "./PVfigures/withBS/compare_FD_Ds_PV.png"); */
    
    /* draw_1d(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2))", ";withBS #DeltaFD_{#it{xy}}(#it{D_{s}^{+}}) [cm];# events", -0.5, 0.5, sel_withBS, "./PVfigures/withBS/diff_FDxy_Ds_PV.png"); */
    /* draw_1d(mychain, "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_withBS_Z)", ";withBS #DeltaFD_{#it{z}}(#it{D_{s}^{+}}) [cm];# events", -0.5, 0.5, sel_withBS, "./PVfigures/withBS/diff_FDz_Ds_PV.png"); */
    /* draw_1d(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_withBS_Z,2))", ";withBS #DeltaFD(#it{D_{s}^{+}}) [cm];# events", -0.5, 0.5, sel_withBS, "./PVfigures/withBS/diff_FD_Ds_PV.png"); */

    /* draw_1d(mychain, "(sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2)))/sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2))", ";withBS #it{#varepsilon}FD_{#it{xy}}(#it{D_{s}^{+}});# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_FDxy_Ds_PV.png"); */
    /* draw_1d(mychain, "(abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_withBS_Z))/abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z)", ";withBS #it{#varepsilon}FD_{#it{z}}(#it{D_{s}^{+}});# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_FDz_Ds_PV.png"); */
    /* draw_1d(mychain, "(sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_withBS_Z,2)))/sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2))", ";withBS #it{#varepsilon}FD(#it{D_{s}^{+}});# events", -1, 1, sel_withBS, "./PVfigures/withBS/epsilon_FD_Ds_PV.png"); */

    compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2))", "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2))", "Fit noBS", "Fit withBS", "#DeltaFD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", -0.5, 0.5, sel_noBS, sel_withBS, "./PVfigures/compare_FDxy_Ds_PV.png");
    compare(mychain, "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_noBS_Z)", "abs(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z) - abs(match_DsFit_ENDVX_Z-PV_withBS_Z)", "Fit noBS", "Fit withBS", "#DeltaFD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", -0.5, 0.5, sel_noBS, sel_withBS, "./PVfigures/compare_FDz_Ds_PV.png");
    compare(mychain, "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_noBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_noBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_noBS_Z,2))", "sqrt(pow(Gen_Kp_ORIVX_X-Gen_H_ORIVX_X,2)+pow(Gen_Kp_ORIVX_Y-Gen_H_ORIVX_Y,2)+pow(Gen_Kp_ORIVX_Z-Gen_H_ORIVX_Z,2)) - sqrt(pow(match_DsFit_ENDVX_X-PV_withBS_X,2)+pow(match_DsFit_ENDVX_Y-PV_withBS_Y,2)+pow(match_DsFit_ENDVX_Z-PV_withBS_Z,2))", "Fit noBS", "Fit withBS", "#DeltaFD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", -0.5, 0.5, sel_noBS, sel_withBS, "./PVfigures/compare_FD_Ds_PV.png");
    
    compare(mychain, "PV_noBS_XERR", "PV_withBS_XERR", "Fit noBS", "Fit withBS", "#it{#deltax}_{PV} [cm]", 0, 0.003, sel_noBS, sel_withBS, "./PVfigures/compare_PV_XERR.png");
    compare(mychain, "PV_noBS_YERR", "PV_withBS_YERR", "Fit noBS", "Fit withBS", "#it{#deltay}_{PV} [cm]", 0, 0.003, sel_noBS, sel_withBS, "./PVfigures/compare_PV_YERR.png");
    compare(mychain, "PV_noBS_ZERR", "PV_withBS_ZERR", "Fit noBS", "Fit withBS", "#it{#deltaz}_{PV} [cm]", 0, 0.006, sel_noBS, sel_withBS, "./PVfigures/compare_PV_ZERR.png");

    compare(mychain, "abs(PV_noBS_X)/PV_noBS_XERR", "abs(PV_withBS_X)/PV_withBS_XERR", "Fit noBS", "Fit withBS", "sig(#it{x}_{PV}) [cm]", 0, 100, sel_noBS, sel_withBS, "./PVfigures/compare_sig_X.png");
    compare(mychain, "abs(PV_noBS_Y)/PV_noBS_YERR", "abs(PV_withBS_Y)/PV_withBS_YERR", "Fit noBS", "Fit withBS", "sig(#it{y}_{PV}) [cm]", 0, 200, sel_noBS, sel_withBS, "./PVfigures/compare_sig_Y.png");
    compare(mychain, "abs(PV_noBS_Z)/PV_noBS_ZERR", "abs(PV_withBS_Z)/PV_withBS_ZERR", "Fit noBS", "Fit withBS", "sig(#it{z}_{PV}) [cm]", 0, 10000, sel_noBS, sel_withBS, "./PVfigures/compare_sig_Z.png");
    
    delete mychain;

    return 0;
}

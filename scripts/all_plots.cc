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


int all_plots() {
    setTDRStyle();

    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    /* if(gSystem->AccessPathName("./figures/fit_results")) gSystem->MakeDirectory("./figures/fit_results"); */
    if(gSystem->AccessPathName("./figures/gen_match_compare")) gSystem->MakeDirectory("./figures/gen_match_compare");
    if(gSystem->AccessPathName("./figures/APplot")) gSystem->MakeDirectory("./figures/APplot");
    if(gSystem->AccessPathName("./figures/fit_vertex")) gSystem->MakeDirectory("./figures/fit_vertex");
    if(gSystem->AccessPathName("./figures/match_nonmatch_compare")) gSystem->MakeDirectory("./figures/match_nonmatch_compare");
    if(gSystem->AccessPathName("./figures/dxy_dz_compare")) gSystem->MakeDirectory("./figures/dxy_dz_compare");

    /* // compare gen and match */ 
    /* compare(mychain, "Gen_Kp_ETA", "match_Kp_ETA", "Gen", "match", "#it{#eta}(#it{K^{+}})", -6, 6, "1", "1", "./figures/gen_match_compare/Kp_ETA.png"); */
    /* compare(mychain, "Gen_Km_ETA", "match_Km_ETA", "Gen", "match", "#it{#eta}(#it{K^{-}})", -6, 6, "1", "1", "./figures/gen_match_compare/Km_ETA.png"); */
    /* compare(mychain, "Gen_pi_ETA", "match_pi_ETA", "Gen", "match", "#it{#eta}(#it{#pi^{+}})", -6, 6, "1", "1", "./figures/gen_match_compare/pi_ETA.png"); */
    /* compare(mychain, "Gen_phi_ETA", "match_phiFit_phi_ETA", "Gen", "match", "#it{#eta}(#it{#phi})", -6, 6, "1", "1", "./figures/gen_match_compare/phi_ETA.png"); */
    /* compare(mychain, "Gen_Ds_ETA", "match_DsFit_Ds_ETA", "Gen", "match", "#it{#eta}(#it{D_{s}^{+}})", -6, 6, "1", "1", "./figures/gen_match_compare/Ds_ETA.png"); */
    
    /* draw_1d(mychain, "(Gen_Kp_ETA - match_Kp_ETA)/Gen_Kp_ETA", ";#it{#varepsilon#eta}(#it{K^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Kp_ETA.png"); */
    /* draw_1d(mychain, "(Gen_Km_ETA - match_Km_ETA)/Gen_Km_ETA", ";#it{#varepsilon#eta}(#it{K^{-}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Km_ETA.png"); */
    /* draw_1d(mychain, "(Gen_pi_ETA - match_pi_ETA)/Gen_pi_ETA", ";#it{#varepsilon#eta}(#it{#pi^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_pi_ETA.png"); */
    /* draw_1d(mychain, "(Gen_phi_ETA - match_phiFit_phi_ETA)/Gen_phi_ETA", ";#it{#varepsilon#eta}(#it{#phi});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_phi_ETA.png"); */
    /* draw_1d(mychain, "(Gen_Ds_ETA - match_DsFit_Ds_ETA)/Gen_Ds_ETA", ";#it{#varepsilon#eta}(#it{D_{s}^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Ds_ETA.png"); */


    /* compare(mychain, "Gen_Kp_PHI", "match_Kp_PHI", "Gen", "match", "#it{#phi}(#it{K^{+}})", -3.14, 3.14, "1", "1", "./figures/gen_match_compare/Kp_PHI.png"); */
    /* compare(mychain, "Gen_Km_PHI", "match_Km_PHI", "Gen", "match", "#it{#phi}(#it{K^{-}})", -3.14, 3.14, "1", "1", "./figures/gen_match_compare/Km_PHI.png"); */
    /* compare(mychain, "Gen_pi_PHI", "match_pi_PHI", "Gen", "match", "#it{#phi}(#it{#pi^{+}})", -3.14, 3.14, "1", "1", "./figures/gen_match_compare/pi_PHI.png"); */
    /* compare(mychain, "Gen_phi_PHI", "match_phiFit_phi_PHI", "Gen", "match", "#it{#phi}(#it{#phi})", -3.14, 3.14, "1", "1", "./figures/gen_match_compare/phi_PHI.png"); */
    /* compare(mychain, "Gen_Ds_PHI", "match_DsFit_Ds_PHI", "Gen", "match", "#it{#phi}(#it{D_{s}^{+}})", -3.14, 3.14, "1", "1", "./figures/gen_match_compare/Ds_PHI.png"); */

    /* draw_1d(mychain, "(Gen_Kp_PHI - match_Kp_PHI)/Gen_Kp_PHI", ";#it{#varepsilon#phi}(#it{K^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Kp_PHI.png"); */
    /* draw_1d(mychain, "(Gen_Km_PHI - match_Km_PHI)/Gen_Km_PHI", ";#it{#varepsilon#phi}(#it{K^{-}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Km_PHI.png"); */
    /* draw_1d(mychain, "(Gen_pi_PHI - match_pi_PHI)/Gen_pi_PHI", ";#it{#varepsilon#phi}(#it{#pi^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_pi_PHI.png"); */
    /* draw_1d(mychain, "(Gen_phi_PHI - match_phiFit_phi_PHI)/Gen_phi_PHI", ";#it{#varepsilon#phi}(#it{#phi});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_phi_PHI.png"); */
    /* draw_1d(mychain, "(Gen_Ds_PHI - match_DsFit_Ds_PHI)/Gen_Ds_PHI", ";#it{#varepsilon#phi}(#it{D_{s}^{+}});# candidates", -0.5, 0.5, "1", "./figures/gen_match_compare/diff_Ds_PHI.png"); */


    /* compare(mychain, "Gen_Kp_P", "match_Kp_P", "Gen", "match", "#it{p}(#it{K^{+}}) [GeV]", 0, 300, "1", "1", "./figures/gen_match_compare/Kp_P.png"); */
    /* compare(mychain, "Gen_Km_P", "match_Km_P", "Gen", "match", "#it{p}(#it{K^{-}}) [GeV]", 0, 300, "1", "1", "./figures/gen_match_compare/Km_P.png"); */
    /* compare(mychain, "Gen_pi_P", "match_pi_P", "Gen", "match", "#it{p}(#it{#pi^{+}}) [GeV]", 0, 300, "1", "1", "./figures/gen_match_compare/pi_P.png"); */
    /* compare(mychain, "Gen_phi_P", "match_phiFit_phi_P", "Gen", "match", "#it{p}(#it{#phi}) [GeV]", 0, 600, "1", "1", "./figures/gen_match_compare/phi_P.png"); */
    /* compare(mychain, "Gen_Ds_P", "match_DsFit_Ds_P", "Gen", "match", "#it{p}(#it{D_{s}^{+}}) [GeV]", 0, 800, "1", "1", "./figures/gen_match_compare/Ds_P.png"); */

    /* draw_1d(mychain, "(Gen_Kp_P - match_Kp_P)/Gen_Kp_P", ";#it{#varepsilonp}(#it{K^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Kp_P.png"); */
    /* draw_1d(mychain, "(Gen_Km_P - match_Km_P)/Gen_Km_P", ";#it{#varepsilonp}(#it{K^{-}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Km_P.png"); */
    /* draw_1d(mychain, "(Gen_pi_P - match_pi_P)/Gen_pi_P", ";#it{#varepsilonp}(#it{#pi^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_pi_P.png"); */
    /* draw_1d(mychain, "(Gen_phi_P - match_phiFit_phi_P)/Gen_phi_P", ";#it{#varepsilonp}(#it{#phi});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_phi_P.png"); */
    /* draw_1d(mychain, "(Gen_Ds_P - match_DsFit_Ds_P)/Gen_Ds_P", ";#it{#varepsilonp}(#it{D_{s}^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Ds_P.png"); */


    /* compare(mychain, "Gen_Kp_PT", "match_Kp_PT", "Gen", "match", "#it{p_{T}}(#it{K^{+}}) [GeV]", 0, 50, "1", "1", "./figures/gen_match_compare/Kp_PT.png"); */
    /* compare(mychain, "Gen_Km_PT", "match_Km_PT", "Gen", "match", "#it{p_{T}}(#it{K^{-}}) [GeV]", 0, 50, "1", "1", "./figures/gen_match_compare/Km_PT.png"); */
    /* compare(mychain, "Gen_pi_PT", "match_pi_PT", "Gen", "match", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", 0, 80, "1", "1", "./figures/gen_match_compare/pi_PT.png"); */
    /* compare(mychain, "Gen_phi_PT", "match_phiFit_phi_PT", "Gen", "match", "#it{p_{T}}(#it{#phi}) [GeV]", 0, 100, "1", "1", "./figures/gen_match_compare/phi_PT.png"); */
    /* compare(mychain, "Gen_Ds_PT", "match_DsFit_Ds_PT", "Gen", "match", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", 0, 150, "1", "1", "./figures/gen_match_compare/Ds_PT.png"); */

    /* draw_1d(mychain, "(Gen_Kp_PT - match_Kp_PT)/Gen_Kp_PT", ";#it{#varepsilonp_{T}}(#it{K^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Kp_PT.png"); */
    /* draw_1d(mychain, "(Gen_Km_PT - match_Km_PT)/Gen_Km_PT", ";#it{#varepsilonp_{T}}(#it{K^{-}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Km_PT.png"); */
    /* draw_1d(mychain, "(Gen_pi_PT - match_pi_PT)/Gen_pi_PT", ";#it{#varepsilonp_{T}}(#it{#pi^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_pi_PT.png"); */
    /* draw_1d(mychain, "(Gen_phi_PT - match_phiFit_phi_PT)/Gen_phi_PT", ";#it{#varepsilonp_{T}}(#it{#phi});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_phi_PT.png"); */
    /* draw_1d(mychain, "(Gen_Ds_PT - match_DsFit_Ds_PT)/Gen_Ds_PT", ";#it{#varepsilonp_{T}}(#it{D_{s}^{+}});# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_Ds_PT.png"); */


    /* compare(mychain, "Gen_Kp_ORIVX_X", "match_Kp_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{K^{+}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/Kp_ORIVX_X.png"); */
    /* compare(mychain, "Gen_Km_ORIVX_X", "match_Km_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{K^{-}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/Km_ORIVX_X.png"); */
    /* compare(mychain, "Gen_pi_ORIVX_X", "match_pi_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{#pi^{+}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/pi_ORIVX_X.png"); */

    /* draw_1d(mychain, "(Gen_Kp_ORIVX_X - match_Kp_ORIVX_X)", ";#Delta#it{x}_{orig vtx}(#it{K^{+}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_Kp_ORIVX_X.png"); */
    /* draw_1d(mychain, "(Gen_Km_ORIVX_X - match_Km_ORIVX_X)", ";#Delta#it{x}_{orig vtx}(#it{K^{-}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_Km_ORIVX_X.png"); */
    /* draw_1d(mychain, "(Gen_pi_ORIVX_X - match_pi_ORIVX_X)", ";#Delta#it{x}_{orig vtx}(#it{#pi^{+}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_pi_ORIVX_X.png"); */


    /* compare(mychain, "Gen_Kp_ORIVX_Y", "match_Kp_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{K^{+}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/Kp_ORIVX_Y.png"); */
    /* compare(mychain, "Gen_Km_ORIVX_Y", "match_Km_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{K^{-}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/Km_ORIVX_Y.png"); */
    /* compare(mychain, "Gen_pi_ORIVX_Y", "match_pi_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{#pi^{+}}) [cm]",  -0.5, 0.5, "1", "1", "./figures/gen_match_compare/pi_ORIVX_Y.png"); */
    
    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Y - match_Kp_ORIVX_Y)", ";#Delta#it{y}_{orig vtx}(#it{K^{+}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_Kp_ORIVX_Y.png"); */
    /* draw_1d(mychain, "(Gen_Km_ORIVX_Y - match_Km_ORIVX_Y)", ";#Delta#it{y}_{orig vtx}(#it{K^{-}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_Km_ORIVX_Y.png"); */
    /* draw_1d(mychain, "(Gen_pi_ORIVX_Y - match_pi_ORIVX_Y)", ";#Delta#it{y}_{orig vtx}(#it{#pi^{+}});# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_pi_ORIVX_Y.png"); */


    /* compare(mychain, "Gen_Kp_ORIVX_Z", "match_Kp_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{K^{+}})",  -20, 20, "1", "1", "./figures/gen_match_compare/Kp_ORIVX_Z.png"); */
    /* compare(mychain, "Gen_Km_ORIVX_Z", "match_Km_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{K^{-}})",  -20, 20, "1", "1", "./figures/gen_match_compare/Km_ORIVX_Z.png"); */
    /* compare(mychain, "Gen_pi_ORIVX_Z", "match_pi_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{#pi^{+}})",  -20, 20, "1", "1", "./figures/gen_match_compare/pi_ORIVX_Z.png"); */

    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Z - match_Kp_ORIVX_Z)", ";#Delta#it{z}_{orig vtx}(#it{K^{+}});# candidates", -10, 10, "1", "./figures/gen_match_compare/diff_Kp_ORIVX_Z.png"); */
    /* draw_1d(mychain, "(Gen_Km_ORIVX_Z - match_Km_ORIVX_Z)", ";#Delta#it{z}_{orig vtx}(#it{K^{-}});# candidates", -10, 10, "1", "./figures/gen_match_compare/diff_Km_ORIVX_Z.png"); */
    /* draw_1d(mychain, "(Gen_pi_ORIVX_Z - match_pi_ORIVX_Z)", ";#Delta#it{z}_{orig vtx}(#it{#pi^{+}});# candidates", -10, 10, "1", "./figures/gen_match_compare/diff_pi_ORIVX_Z.png"); */


    /* compare(mychain, "Gen_dR_Kp_Km", "match_dR_Kp_Km", "Gen", "match", "#Delta#it{R(K^{+},K^{-})}", 0, 0.2, "1", "1", "./figures/gen_match_compare/dR_Kp_Km.png"); */
    /* compare(mychain, "Gen_dR_Kp_pi", "match_dR_Kp_pi", "Gen", "match", "#Delta#it{R(K^{+},#pi^{+})}", 0, 0.8, "1", "1", "./figures/gen_match_compare/dR_Kp_pi.png"); */
    /* compare(mychain, "Gen_dR_Km_pi", "match_dR_Km_pi", "Gen", "match", "#Delta#it{R(K^{-},#pi^{+})}", 0, 0.8, "1", "1", "./figures/gen_match_compare/dR_Km_pi.png"); */
    
    /* draw_1d(mychain, "(Gen_dR_Kp_Km - match_dR_Kp_Km)", ";Gen #Delta#it{R(K^{+},K^{-})}-reco #Delta#it{R(K^{+},K^{-})};# candidates", -0.2, 0.2, "1", "./figures/gen_match_compare/diff_dR_Kp_Km.png"); */
    /* draw_1d(mychain, "(Gen_dR_Kp_pi - match_dR_Kp_pi)", ";Gen #Delta#it{R(K^{+},#pi^{+})}-reco #Delta#it{R(K^{+},#pi^{+})};# candidates", -0.8, 0.8, "1", "./figures/gen_match_compare/diff_dR_Kp_pi.png"); */
    /* draw_1d(mychain, "(Gen_dR_Km_pi - match_dR_Km_pi)", ";Gen #Delta#it{R(K^{-},#pi^{+})}-reco #Delta#it{R(K^{-},#pi^{+})};# candidates", -0.8, 0.8, "1", "./figures/gen_match_compare/diff_dR_Km_pi.png"); */
   


    
    /* // AP plot */
    /* draw_2d(mychain, "Gen_Kp_PP:((Gen_Kp_PL - Gen_Km_PL)/(Gen_Kp_PL + Gen_Km_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "1", "./figures/APplot/Gen_APplot_phi.png"); */
    /* draw_2d(mychain, "match_phiFit_Kp_PP:((match_phiFit_Kp_PL - match_phiFit_Km_PL)/(match_phiFit_Kp_PL + match_phiFit_Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "1", "./figures/APplot/match_APplot_phi.png"); */
    /* draw_2d(mychain, "phiFit_Kp_PP:((phiFit_Kp_PL - phiFit_Km_PL)/(phiFit_Kp_PL + phiFit_Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "non_match_entry==1", "./figures/APplot/nonmatch_APplot_phi.png"); */
    
    /* draw_2d(mychain, "Gen_phi_PP:((Gen_phi_PL - Gen_pi_PL)/(Gen_phi_PL + Gen_pi_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "1", "./figures/APplot/Gen_APplot_Ds.png"); */
    /* draw_2d(mychain, "match_DsFit_phi_PP:((match_DsFit_phi_PL - match_DsFit_pi_PL)/(match_DsFit_phi_PL + match_DsFit_pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "1", "./figures/APplot/match_APplot_Ds.png"); */
    /* draw_2d(mychain, "DsFit_phi_PP:((DsFit_phi_PL - DsFit_pi_PL)/(DsFit_phi_PL + DsFit_pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "non_match_entry==1", "./figures/APplot/nonmatch_APplot_Ds.png"); */





    // fit vertex
    compare(mychain, "Gen_Kp_ORIVX_X", "match_phiFit_ENDVX_X", "Gen", "match", "vtx_{#it{x}}(#it{#phi}) [cm]", -1, 1, "1", "1", "./figures/fit_vertex/phi_VX.png");
    compare(mychain, "Gen_Kp_ORIVX_Y", "match_phiFit_ENDVX_Y", "Gen", "match", "vtx_{#it{y}}(#it{#phi}) [cm]", -1, 1, "1", "1", "./figures/fit_vertex/phi_VY.png");
    compare(mychain, "Gen_Kp_ORIVX_Z", "match_phiFit_ENDVX_Z", "Gen", "match", "vtx_{#it{z}}(#it{#phi}) [cm]", -15, 15, "1", "1", "./figures/fit_vertex/phi_VZ.png");

    compare(mychain, "Gen_Kp_ORIVX_X", "phiFit_ENDVX_X", "Gen", "non-match", "vtx_{#it{x}}(#it{#phi}) [cm]", -1, 1, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_VX.png");
    compare(mychain, "Gen_Kp_ORIVX_Y", "phiFit_ENDVX_Y", "Gen", "non-match", "vtx_{#it{y}}(#it{#phi}) [cm]", -1, 1, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_VY.png");
    compare(mychain, "Gen_Kp_ORIVX_Z", "phiFit_ENDVX_Z", "Gen", "non-match", "vtx_{#it{z}}(#it{#phi}) [cm]", -15, 15, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_VZ.png");
    
    draw_1d(mychain, "(Gen_Kp_ORIVX_X - match_phiFit_ENDVX_X)", ";fit vtx_{#it{x}}(#it{#phi})-Gen vtx_{#it{x}}(#it{#phi}) [cm];# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_phi_VX.png");
    draw_1d(mychain, "(Gen_Kp_ORIVX_Y - match_phiFit_ENDVX_Y)", ";fit vtx_{#it{y}}(#it{#phi})-Gen vtx_{#it{y}}(#it{#phi}) [cm];# candidates", -1, 1, "1", "./figures/gen_match_compare/diff_phi_VY.png");
    draw_1d(mychain, "(Gen_Kp_ORIVX_Z - match_phiFit_ENDVX_Z)", ";fit vtx_{#it{z}}(#it{#phi})-Gen vtx_{#it{z}}(#it{#phi}) [cm];# candidates", -2, 2, "1", "./figures/gen_match_compare/diff_phi_VZ.png");

    compare(mychain, "phiFit_ENDVX_XERR", "phiFit_ENDVX_XERR", "match", "non-match", "#Deltavtx_{#it{x}}(#it{#phi}) [cm]", 0, 1, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_VXERR.png");
    compare(mychain, "phiFit_ENDVX_YERR", "phiFit_ENDVX_YERR", "match", "non-match", "#Deltavtx_{#it{y}}(#it{#phi}) [cm]", 0, 1, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_VYERR.png");
    compare(mychain, "phiFit_ENDVX_ZERR", "phiFit_ENDVX_ZERR", "match", "non-match", "#Deltavtx_{#it{z}}(#it{#phi}) [cm]", 0, 2, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_VZERR.png");

    draw_2d(mychain, "match_phiFit_ENDVX_XERR:(match_phiFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";fit vtx_{#it{x}}(#it{#phi})-Gen vtx_{#it{x}}(#it{#phi}) [cm];fit #Deltavtx_{#it{x}}(#it{#phi})|", -1, 1, 0, 0.6, "1", "./figures/fit_vertex/phi_err_vtx_x.png");
    draw_2d(mychain, "match_phiFit_ENDVX_YERR:(match_phiFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";fit vtx_{#it{y}}(#it{#phi})-Gen vtx_{#it{y}}(#it{#phi}) [cm];fit #Deltavtx_{#it{y}}(#it{#phi})|", -1, 1, 0, 0.6, "1", "./figures/fit_vertex/phi_err_vtx_y.png");
    draw_2d(mychain, "match_phiFit_ENDVX_ZERR:(match_phiFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";fit vtx_{#it{z}}(#it{#phi})-Gen vtx_{#it{z}}(#it{#phi}) [cm];fit #Deltavtx_{#it{z}}(#it{#phi})|", -2, 2, 0, 2, "1", "./figures/fit_vertex/phi_err_vtx_z.png");
    
    /* draw_2d(mychain, "phiFit_ENDVX_XERR:(phiFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";nonmatch fit vtx_{#it{x}}(#it{#phi})-Gen vtx_{#it{x}}(#it{#phi}) [cm];nonmatch fit #Deltavtx_{#it{x}}(#it{#phi})|", -1, 1, 0, 0.6, "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_err_vtx_x.png"); */
    /* draw_2d(mychain, "phiFit_ENDVX_YERR:(phiFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";nonmatch fit vtx_{#it{y}}(#it{#phi})-Gen vtx_{#it{y}}(#it{#phi}) [cm];nonmatch fit #Deltavtx_{#it{y}}(#it{#phi})|", -1, 1, 0, 0.6, "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_err_vtx_y.png"); */
    /* draw_2d(mychain, "phiFit_ENDVX_ZERR:(phiFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";nonmatch fit vtx_{#it{z}}(#it{#phi})-Gen vtx_{#it{z}}(#it{#phi}) [cm];nonmatch fit #Deltavtx_{#it{z}}(#it{#phi})|", -2, 2, 0, 2, "non_match_entry==1", "./figures/fit_vertex/nonmatch_phi_err_vtx_z.png"); */

    /* draw_1d(mychain, "(Gen_Kp_ORIVX_X - phiFit_ENDVX_X)", ";nonmatch fit vtx_{#it{x}}(#it{#phi})-Gen vtx_{#it{x}}(#it{#phi}) [cm];# candidates", -1, 1, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_phi_VX.png"); */
    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Y - phiFit_ENDVX_Y)", ";nonmatch fit vtx_{#it{y}}(#it{#phi})-Gen vtx_{#it{y}}(#it{#phi}) [cm];# candidates", -1, 1, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_phi_VY.png"); */
    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Z - phiFit_ENDVX_Z)", ";nonmatch fit vtx_{#it{z}}(#it{#phi})-Gen vtx_{#it{z}}(#it{#phi}) [cm];# candidates", -2, 2, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_phi_VZ.png"); */


    compare(mychain, "Gen_Kp_ORIVX_X", "match_DsFit_ENDVX_X", "Gen", "match", "vtx_{#it{x}}(#it{D_{s}^{+}}) [cm]", -1, 1, "1", "1", "./figures/fit_vertex/Ds_VX.png");
    compare(mychain, "Gen_Kp_ORIVX_Y", "match_DsFit_ENDVX_Y", "Gen", "match", "vtx_{#it{y}}(#it{D_{s}^{+}}) [cm]", -1, 1, "1", "1", "./figures/fit_vertex/Ds_VY.png");
    compare(mychain, "Gen_Kp_ORIVX_Z", "match_DsFit_ENDVX_Z", "Gen", "match", "vtx_{#it{z}}(#it{D_{s}^{+}}) [cm]", -15, 15, "1", "1", "./figures/fit_vertex/Ds_VZ.png");

    compare(mychain, "Gen_Kp_ORIVX_X", "DsFit_ENDVX_X", "Gen", "non-match", "vtx_{#it{x}}(#it{D_{s}^{+}}) [cm]", -1, 1, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_VX.png");
    compare(mychain, "Gen_Kp_ORIVX_Y", "DsFit_ENDVX_Y", "Gen", "non-match", "vtx_{#it{y}}(#it{D_{s}^{+}}) [cm]", -1, 1, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_VY.png");
    compare(mychain, "Gen_Kp_ORIVX_Z", "DsFit_ENDVX_Z", "Gen", "non-match", "vtx_{#it{z}}(#it{D_{s}^{+}}) [cm]", -15, 15, "1", "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_VZ.png");
    
    draw_1d(mychain, "(Gen_Kp_ORIVX_X - match_DsFit_ENDVX_X)", ";fit vtx_{#it{x}}(#it{D_{s}^{+}})-Gen vtx_{#it{x}}(#it{D_{s}^{+}}) [cm];# candidates", -0.5, 0.5, "1", "./figures/fit_vertex/diff_Ds_VX.png");
    draw_1d(mychain, "(Gen_Kp_ORIVX_Y - match_DsFit_ENDVX_Y)", ";fit vtx_{#it{y}}(#it{D_{s}^{+}})-Gen vtx_{#it{y}}(#it{D_{s}^{+}}) [cm];# candidates", -0.5, 0.5, "1", "./figures/fit_vertex/diff_Ds_VY.png");
    draw_1d(mychain, "(Gen_Kp_ORIVX_Z - match_DsFit_ENDVX_Z)", ";fit vtx_{#it{z}}(#it{D_{s}^{+}})-Gen vtx_{#it{z}}(#it{D_{s}^{+}}) [cm];# candidates", -1, 1, "1", "./figures/fit_vertex/diff_Ds_VZ.png");

    compare(mychain, "DsFit_ENDVX_XERR", "DsFit_ENDVX_XERR", "match", "non-match", "#Deltavtx_{#it{x}}(#it{D_{s}^{+}}) [cm]", 0, 0.25, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/Ds_VXERR.png");
    compare(mychain, "DsFit_ENDVX_YERR", "DsFit_ENDVX_YERR", "match", "non-match", "#Deltavtx_{#it{y}}(#it{D_{s}^{+}}) [cm]", 0, 0.25, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/Ds_VYERR.png");
    compare(mychain, "DsFit_ENDVX_ZERR", "DsFit_ENDVX_ZERR", "match", "non-match", "#Deltavtx_{#it{z}}(#it{D_{s}^{+}}) [cm]", 0, 0.5, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/Ds_VZERR.png");

    draw_2d(mychain, "match_DsFit_ENDVX_XERR:(match_DsFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";fit vtx_{#it{x}}(#it{D_{s}^{+}})-Gen vtx_{#it{x}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{x}}(#it{D_{s}^{+}})|", -0.3, 0.3, 0, 0.1, "1", "./figures/fit_vertex/Ds_err_vtx_x.png");
    draw_2d(mychain, "match_DsFit_ENDVX_YERR:(match_DsFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";fit vtx_{#it{y}}(#it{D_{s}^{+}})-Gen vtx_{#it{y}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{y}}(#it{D_{s}^{+}})|", -0.3, 0.3, 0, 0.1, "1", "./figures/fit_vertex/Ds_err_vtx_y.png");
    draw_2d(mychain, "match_DsFit_ENDVX_ZERR:(match_DsFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";fit vtx_{#it{z}}(#it{D_{s}^{+}})-Gen vtx_{#it{z}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{z}}(#it{D_{s}^{+}})|", -0.5, 0.5, 0, 0.3, "1", "./figures/fit_vertex/Ds_err_vtx_z.png");

    /* draw_2d(mychain, "DsFit_ENDVX_XERR:(DsFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";nonmatch fit vtx_{#it{x}}(#it{D_{s}^{+}})-Gen vtx_{#it{x}}(#it{D_{s}^{+}}) [cm];nonmatch fit #Deltavtx_{#it{x}}(#it{D_{s}^{+}})|", -1, 1, 0, 0.6, "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_err_vtx_x.png"); */
    /* draw_2d(mychain, "DsFit_ENDVX_YERR:(DsFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";nonmatch fit vtx_{#it{y}}(#it{D_{s}^{+}})-Gen vtx_{#it{y}}(#it{D_{s}^{+}}) [cm];nonmatch fit #Deltavtx_{#it{y}}(#it{D_{s}^{+}})|", -1, 1, 0, 0.6, "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_err_vtx_y.png"); */
    /* draw_2d(mychain, "DsFit_ENDVX_ZERR:(DsFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";nonmatch fit vtx_{#it{z}}(#it{D_{s}^{+}})-Gen vtx_{#it{z}}(#it{D_{s}^{+}}) [cm];nonmatch fit #Deltavtx_{#it{z}}(#it{D_{s}^{+}})|", -2, 2, 0, 2, "non_match_entry==1", "./figures/fit_vertex/nonmatch_Ds_err_vtx_z.png"); */

    /* draw_1d(mychain, "(Gen_Kp_ORIVX_X - DsFit_ENDVX_X)", ";nonmatch fit vtx_{#it{x}}(#it{D_{s}^{+}})-Gen vtx_{#it{x}}(#it{D_{s}^{+}}) [cm];# candidates", -1, 1, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_Ds_VX.png"); */
    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Y - DsFit_ENDVX_Y)", ";nonmatch fit vtx_{#it{y}}(#it{D_{s}^{+}})-Gen vtx_{#it{y}}(#it{D_{s}^{+}}) [cm];# candidates", -1, 1, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_Ds_VY.png"); */
    /* draw_1d(mychain, "(Gen_Kp_ORIVX_Z - DsFit_ENDVX_Z)", ";nonmatch fit vtx_{#it{z}}(#it{D_{s}^{+}})-Gen vtx_{#it{z}}(#it{D_{s}^{+}}) [cm];# candidates", -2, 2, "non_match_entry==1", "./figures/gen_match_compare/nonmatch_diff_Ds_VZ.png"); */


    compare(mychain, "phiFit_ENDVX_X-DsFit_ENDVX_X", "phiFit_ENDVX_X-DsFit_ENDVX_X", "match", "non-match", "vtx_{#it{x}}(#phi)-vtx_{#it{x}}(#it{D_{s}^{+}}) [cm]", -1, 1, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_Ds_VX.png");
    compare(mychain, "phiFit_ENDVX_Y-DsFit_ENDVX_Y", "phiFit_ENDVX_Y-DsFit_ENDVX_Y", "match", "non-match", "vtx_{#it{y}}(#phi)-vtx_{#it{y}}(#it{D_{s}^{+}}) [cm]", -1, 1, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_Ds_VY.png");
    compare(mychain, "phiFit_ENDVX_Z-DsFit_ENDVX_Z", "phiFit_ENDVX_Z-DsFit_ENDVX_Z", "match", "non-match", "vtx_{#it{z}}(#phi)-vtx_{#it{z}}(#it{D_{s}^{+}}) [cm]", -2, 2, "match_entry==1", "non_match_entry==1", "./figures/fit_vertex/phi_Ds_VZ.png");



/*     // match nonmatch compare */
/*     compare(mychain, "dxy_Kp_Km", "dxy_Kp_Km", "match", "non-match", "#it{d_{xy}(K^{+},K^{-})} [cm]", 0, 0.1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Kp_Km.png"); */
/*     compare(mychain, "dxy_Kp_pi", "dxy_Kp_pi", "match", "non-match", "#it{d_{xy}(K^{+},#pi^{+})} [cm]", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Kp_pi.png"); */
/*     compare(mychain, "dxy_Km_pi", "dxy_Km_pi", "match", "non-match", "#it{d_{xy}(K^{-},#pi^{+})} [cm]", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Km_pi.png"); */
/*     compare(mychain, "dxy_Kp_phi", "dxy_Kp_phi", "match", "non-match", "#it{d_{xy}(K^{+},#phi)} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Kp_phi.png"); */
/*     compare(mychain, "dxy_Km_phi", "dxy_Km_phi", "match", "non-match", "#it{d_{xy}(K^{-},#phi)} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Km_phi.png"); */
/*     compare(mychain, "dxy_pi_phi", "dxy_pi_phi", "match", "non-match", "#it{d_{xy}(#pi^{+},#phi)} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_pi_phi.png"); */
/*     compare(mychain, "dxy_Kp_Ds", "dxy_Kp_Ds", "match", "non-match", "#it{d_{xy}(K^{+},D_{s}^{+})} [cm]", 0, 1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Kp_Ds.png"); */
/*     compare(mychain, "dxy_Km_Ds", "dxy_Km_Ds", "match", "non-match", "#it{d_{xy}(K^{-},D_{s}^{+})} [cm]", 0, 1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_Km_Ds.png"); */
/*     compare(mychain, "dxy_pi_Ds", "dxy_pi_Ds", "match", "non-match", "#it{d_{xy}(#pi^{+},D_{s}^{+})} [cm]", 0, 1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_pi_Ds.png"); */
/*     compare(mychain, "dxy_phi_Ds", "dxy_phi_Ds", "match", "non-match", "#it{d_{xy}(#phi,D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dxy_phi_Ds.png"); */

/*     compare(mychain, "dz_Kp_Km", "dz_Kp_Km", "match", "non-match", "#it{d_{z}(K^{+},K^{-})} [cm]", 0, 0.2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Kp_Km.png"); */
/*     compare(mychain, "dz_Kp_pi", "dz_Kp_pi", "match", "non-match", "#it{d_{z}(K^{+},#pi^{+})} [cm]", 0, 0.2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Kp_pi.png"); */
/*     compare(mychain, "dz_Km_pi", "dz_Km_pi", "match", "non-match", "#it{d_{z}(K^{-},#pi^{+})} [cm]", 0, 0.2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Km_pi.png"); */
/*     compare(mychain, "dz_Kp_phi", "dz_Kp_phi", "match", "non-match", "#it{d_{z}(K^{+},#phi)} [cm]", 0, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Kp_phi.png"); */
/*     compare(mychain, "dz_Km_phi", "dz_Km_phi", "match", "non-match", "#it{d_{z}(K^{-},#phi)} [cm]", 0, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Km_phi.png"); */
/*     compare(mychain, "dz_pi_phi", "dz_pi_phi", "match", "non-match", "#it{d_{z}(#pi^{+},#phi)} [cm]", 0, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_pi_phi.png"); */
/*     compare(mychain, "dz_Kp_Ds", "dz_Kp_Ds", "match", "non-match", "#it{d_{z}(K^{+},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Kp_Ds.png"); */
/*     compare(mychain, "dz_Km_Ds", "dz_Km_Ds", "match", "non-match", "#it{d_{z}(K^{-},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_Km_Ds.png"); */
/*     compare(mychain, "dz_pi_Ds", "dz_pi_Ds", "match", "non-match", "#it{d_{z}(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_pi_Ds.png"); */
/*     compare(mychain, "dz_phi_Ds", "dz_phi_Ds", "match", "non-match", "#it{d_{z}(#phi,D_{s}^{+})} [cm]", 0, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dz_phi_Ds.png"); */

/*     compare(mychain, "dR_Kp_Km", "dR_Kp_Km", "match", "non-match", "#it{#DeltaR(K^{+},K^{-})}", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dR_Kp_Km.png"); */
/*     compare(mychain, "dR_Kp_pi", "dR_Kp_pi", "match", "non-match", "#it{#DeltaR(K^{+},#pi^{+})}", 0, 0.6, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dR_Kp_pi.png"); */
/*     compare(mychain, "dR_Km_pi", "dR_Km_pi", "match", "non-match", "#it{#DeltaR(K^{-},#pi^{+})}", 0, 0.6, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/dR_Km_pi.png"); */
/*     compare(mychain, "phiFit_dR_Kp_phi", "phiFit_dR_Kp_phi", "match", "non-match", "#it{#DeltaR(K^{+},#phi)}", 0, 0.1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_dR_Kp_phi.png"); */
/*     compare(mychain, "phiFit_dR_Km_phi", "phiFit_dR_Km_phi", "match", "non-match", "#it{#DeltaR(K^{-},#phi)}", 0, 0.1, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_dR_Km_phi.png"); */
/*     compare(mychain, "phiFit_dR_pi_phi", "phiFit_dR_pi_phi", "match", "non-match", "#it{#DeltaR(#pi^{+},#phi)}", 0, 0.6, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_dR_pi_phi.png"); */
/*     compare(mychain, "DsFit_dR_Kp_Ds", "DsFit_dR_Kp_Ds", "match", "non-match", "#it{#DeltaR(K^{+},D_{s}^{+})}", 0, 0.4, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_dR_Kp_Ds.png"); */
/*     compare(mychain, "DsFit_dR_Km_Ds", "DsFit_dR_Km_Ds", "match", "non-match", "#it{#DeltaR(K^{-},D_{s}^{+})}", 0, 0.4, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_dR_Km_Ds.png"); */
/*     compare(mychain, "DsFit_dR_pi_Ds", "DsFit_dR_pi_Ds", "match", "non-match", "#it{#DeltaR(#pi^{+},D_{s}^{+})}", 0, 0.6, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_dR_pi_Ds.png"); */
/*     compare(mychain, "DsFit_dR_phi_Ds", "DsFit_dR_phi_Ds", "match", "non-match", "#it{#DeltaR(#phi,D_{s}^{+})}", 0, 0.4, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_dR_phi_Ds.png"); */

/*     compare(mychain, "Kp_PT", "Kp_PT", "match", "non-match", "#it{p_{T}(K^{+})} [GeV]", 0, 30, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_PT.png"); */
/*     compare(mychain, "Km_PT", "Km_PT", "match", "non-match", "#it{p_{T}(K^{-})} [GeV]", 0, 30, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_PT.png"); */
/*     compare(mychain, "pi_PT", "pi_PT", "match", "non-match", "#it{p_{T}(#pi^{+})} [GeV]", 0, 30, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_PT.png"); */
/*     compare(mychain, "phiFit_phi_PT", "phiFit_phi_PT", "match", "non-match", "#it{p_{T}(#phi)} [GeV]", 0, 60, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_phi_PT.png"); */
/*     compare(mychain, "DsFit_Ds_PT", "DsFit_Ds_PT", "match", "non-match", "#it{p_{T}(D_{s}^{+})} [GeV]", 0, 80, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_Ds_PT.png"); */

/*     compare(mychain, "Kp_P", "Kp_P", "match", "non-match", "#it{p(K^{+})} [GeV]", 0, 50, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_P.png"); */
/*     compare(mychain, "Km_P", "Km_P", "match", "non-match", "#it{p(K^{-})} [GeV]", 0, 50, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_P.png"); */
/*     compare(mychain, "pi_P", "pi_P", "match", "non-match", "#it{p(#pi^{+})} [GeV]", 0, 50, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_P.png"); */
/*     compare(mychain, "phiFit_phi_P", "phiFit_phi_P", "match", "non-match", "#it{p(#phi)} [GeV]", 0, 100, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_phi_P.png"); */
/*     compare(mychain, "DsFit_Ds_P", "DsFit_Ds_P", "match", "non-match", "#it{p(D_{s}^{+})} [GeV]", 0, 150, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_Ds_P.png"); */

/*     compare(mychain, "Kp_PHI", "Kp_PHI", "match", "non-match", "#it{#phi(K^{+})}", -3.15, 3.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_PHI.png"); */
/*     compare(mychain, "Km_PHI", "Km_PHI", "match", "non-match", "#it{#phi(K^{-})}", -3.15, 3.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_PHI.png"); */
/*     compare(mychain, "pi_PHI", "pi_PHI", "match", "non-match", "#it{#phi(#pi^{+})}", -3.15, 3.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_PHI.png"); */
/*     compare(mychain, "phiFit_phi_PHI", "phiFit_phi_PHI", "match", "non-match", "#it{#phi(#phi)}", -3.15, 3.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_phi_PHI.png"); */
/*     compare(mychain, "DsFit_Ds_PHI", "DsFit_Ds_PHI", "match", "non-match", "#it{#phi(D_{s}^{+})}", -3.15, 3.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_Ds_PHI.png"); */

/*     compare(mychain, "Kp_ETA", "Kp_ETA", "match", "non-match", "#it{#eta(K^{+})}", -5, 5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_ETA.png"); */
/*     compare(mychain, "Km_ETA", "Km_ETA", "match", "non-match", "#it{#eta(K^{-})}", -5, 5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_ETA.png"); */
/*     compare(mychain, "pi_ETA", "pi_ETA", "match", "non-match", "#it{#eta(#pi^{+})}", -5, 5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_ETA.png"); */
/*     compare(mychain, "phiFit_phi_ETA", "phiFit_phi_ETA", "match", "non-match", "#it{#eta(#phi)}", -5, 5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_phi_ETA.png"); */
/*     compare(mychain, "DsFit_Ds_ETA", "DsFit_Ds_ETA", "match", "non-match", "#it{#eta(D_{s}^{+})}", -5, 5, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_Ds_ETA.png"); */

/*     compare(mychain, "Kp_ORIVX_X", "Kp_ORIVX_X", "match", "non-match", "#it{x}_{orig}#it{(K^{+})} [cm]", -0.1, 0.05, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_ORIVX_X.png"); */
/*     compare(mychain, "Kp_ORIVX_Y", "Kp_ORIVX_Y", "match", "non-match", "#it{y}_{orig}#it{(K^{+})} [cm]", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_ORIVX_Y.png"); */
/*     compare(mychain, "Kp_ORIVX_Z", "Kp_ORIVX_Z", "match", "non-match", "#it{z}_{orig}#it{(K^{+})} [cm]", -15, 15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Kp_ORIVX_Z.png"); */

/*     compare(mychain, "Km_ORIVX_X", "Km_ORIVX_X", "match", "non-match", "#it{x}_{orig}#it{(K^{-})} [cm]", -0.1, 0.05, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_ORIVX_X.png"); */
/*     compare(mychain, "Km_ORIVX_Y", "Km_ORIVX_Y", "match", "non-match", "#it{y}_{orig}#it{(K^{-})} [cm]", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_ORIVX_Y.png"); */
/*     compare(mychain, "Km_ORIVX_Z", "Km_ORIVX_Z", "match", "non-match", "#it{z}_{orig}#it{(K^{-})} [cm]", -15, 15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/Km_ORIVX_Z.png"); */

/*     compare(mychain, "pi_ORIVX_X", "pi_ORIVX_X", "match", "non-match", "#it{x}_{orig}#it{(#pi^{+})} [cm]", -0.1, 0.05, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_ORIVX_X.png"); */
/*     compare(mychain, "pi_ORIVX_Y", "pi_ORIVX_Y", "match", "non-match", "#it{y}_{orig}#it{(#pi^{+})} [cm]", 0, 0.15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_ORIVX_Y.png"); */
/*     compare(mychain, "pi_ORIVX_Z", "pi_ORIVX_Z", "match", "non-match", "#it{z}_{orig}#it{(#pi^{+})} [cm]", -15, 15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/pi_ORIVX_Z.png"); */

/*     compare(mychain, "phiFit_ENDVX_X", "phiFit_ENDVX_X", "match", "non-match", "#it{x}_{decay}(#it{#phi}) [cm]", -2, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_ENDVX_X.png"); */
/*     compare(mychain, "phiFit_ENDVX_Y", "phiFit_ENDVX_Y", "match", "non-match", "#it{y}_{decay}(#it{#phi}) [cm]", -2, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_ENDVX_Y.png"); */
/*     compare(mychain, "phiFit_ENDVX_Z", "phiFit_ENDVX_Z", "match", "non-match", "#it{z}_{decay}(#it{#phi}) [cm]", -15, 15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/phiFit_ENDVX_Z.png"); */

/*     compare(mychain, "DsFit_ENDVX_X", "DsFit_ENDVX_X", "match", "non-match", "#it{x}_{decay}(#it{D_{s}^{+}}) [cm]", -2, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_ENDVX_X.png"); */
/*     compare(mychain, "DsFit_ENDVX_Y", "DsFit_ENDVX_Y", "match", "non-match", "#it{y}_{decay}(#it{D_{s}^{+}}) [cm]", -2, 2, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_ENDVX_Y.png"); */
/*     compare(mychain, "DsFit_ENDVX_Z", "DsFit_ENDVX_Z", "match", "non-match", "#it{z}_{decay}(#it{D_{s}^{+}}) [cm]", -15, 15, "match_entry==1", "non_match_entry==1", "./figures/match_nonmatch_compare/DsFit_ENDVX_Z.png"); */







/*     // dxy dz compare */
/*     compare(mychain, "dxy_Kp_Km", "dz_Kp_Km", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},K^{-})} [cm]", 0, 0.2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Kp_Km.png"); */
/*     compare(mychain, "dxy_Kp_pi", "dz_Kp_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#pi^{+})} [cm]", 0, 0.2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Kp_pi.png"); */
/*     compare(mychain, "dxy_Km_pi", "dz_Km_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#pi^{+})} [cm]", 0, 0.2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Km_pi.png"); */
/*     compare(mychain, "dxy_Kp_phi", "dz_Kp_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#phi)} [cm]", 0, 2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Kp_phi.png"); */
/*     compare(mychain, "dxy_Km_phi", "dz_Km_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#phi)} [cm]", 0, 2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Km_phi.png"); */
/*     compare(mychain, "dxy_pi_phi", "dz_pi_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},#phi)} [cm]", 0, 2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_pi_phi.png"); */
/*     compare(mychain, "dxy_Kp_Ds", "dz_Kp_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Kp_Ds.png"); */
/*     compare(mychain, "dxy_Km_Ds", "dz_Km_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_Km_Ds.png"); */
/*     compare(mychain, "dxy_pi_Ds", "dz_pi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_pi_Ds.png"); */
/*     compare(mychain, "dxy_phi_Ds", "dz_phi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#phi,D_{s}^{+})} [cm]", 0, 2, "match_entry==1", "match_entry==1", "./figures/dxy_dz_compare/match_d_phi_Ds.png"); */

/*     compare(mychain, "dxy_Kp_Km", "dz_Kp_Km", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},K^{-})} [cm]", 0, 0.2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Kp_Km.png"); */
/*     compare(mychain, "dxy_Kp_pi", "dz_Kp_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#pi^{+})} [cm]", 0, 0.2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Kp_pi.png"); */
/*     compare(mychain, "dxy_Km_pi", "dz_Km_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#pi^{+})} [cm]", 0, 0.2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Km_pi.png"); */
/*     compare(mychain, "dxy_Kp_phi", "dz_Kp_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#phi)} [cm]", 0, 2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Kp_phi.png"); */
/*     compare(mychain, "dxy_Km_phi", "dz_Km_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#phi)} [cm]", 0, 2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Km_phi.png"); */
/*     compare(mychain, "dxy_pi_phi", "dz_pi_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},#phi)} [cm]", 0, 2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_pi_phi.png"); */
/*     compare(mychain, "dxy_Kp_Ds", "dz_Kp_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},D_{s}^{+})} [cm]", 0, 1.5, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Kp_Ds.png"); */
/*     compare(mychain, "dxy_Km_Ds", "dz_Km_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},D_{s}^{+})} [cm]", 0, 1.5, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_Km_Ds.png"); */
/*     compare(mychain, "dxy_pi_Ds", "dz_pi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_pi_Ds.png"); */
/*     compare(mychain, "dxy_phi_Ds", "dz_phi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#phi,D_{s}^{+})} [cm]", 0, 2, "non_match_entry==1", "non_match_entry==1", "./figures/dxy_dz_compare/nonmatch_d_phi_Ds.png"); */




    // 4 masses
    /* draw_1d(mychain, "phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "all_phi_M", mass_cut); */
    /* draw_1d(mychain, "phiFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "all_phiFit_phi_M", mass_cut); */
    /* draw_1d(mychain, "DsFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "all_DsFit_phi_M", mass_cut); */


    /* draw_1d(mychain, "phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "match_phi_M", match_sel); */
    /* draw_1d(mychain, "phiFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "match_phiFit_phi_M", match_sel); */
    /* draw_1d(mychain, "DsFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "match_DsFit_phi_M", match_sel); */

    /* draw_1d(mychain, "phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "nonmatch_phi_M", nonmatch_sel); */
    /* draw_1d(mychain, "phiFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "nonmatch_phiFit_phi_M", nonmatch_sel); */
    /* draw_1d(mychain, "DsFit_phi_M", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "nonmatch_DsFit_phi_M", nonmatch_sel); */

    /* draw_1d(mychain, "Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "all_Ds_M", mass_cut); */
    /* draw_1d(mychain, "phiFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "all_phiFit_Ds_M", mass_cut); */
    /* draw_1d(mychain, "DsFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "all_DsFit_Ds_M", mass_cut); */
    /* draw_1d(mychain, "DsFit_Mconstraint_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "all_DsFit_Mconstraint_Ds_M", mass_cut); */

    /* draw_1d(mychain, "Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "match_Ds_M", match_sel); */
    /* draw_1d(mychain, "phiFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "match_phiFit_Ds_M", match_sel); */
    /* draw_1d(mychain, "DsFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "match_DsFit_Ds_M", match_sel); */
    /* draw_1d(mychain, "DsFit_Mconstraint_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "match_DsFit_Mconstraint_Ds_M", match_sel); */

    /* draw_1d(mychain, "Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "nonmatch_Ds_M", nonmatch_sel); */
    /* draw_1d(mychain, "phiFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "nonmatch_phiFit_Ds_M", nonmatch_sel); */
    /* draw_1d(mychain, "DsFit_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "nonmatch_DsFit_Ds_M", nonmatch_sel); */
    /* draw_1d(mychain, "DsFit_Mconstraint_Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "nonmatch_DsFit_Mconstraint_Ds_M", nonmatch_sel); */


    delete mychain;

    return 0;
}

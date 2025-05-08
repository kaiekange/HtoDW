#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void draw_fit_results_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", vartitle, 100, varmin, varmax);

    mychain->Project("hist", myvar);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw();
    CMS_lumi(canvas);
	canvas->Update();
	canvas->RedrawAxis();
    canvas->SaveAs("./figures/fit_results/"+figpath+".png");
    delete hist;
}

void draw_fit_results_2d(TChain *mychain, TString myvar, TString vartitle, float xvarmin, float xvarmax, float yvarmin, float yvarmax, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH2F *hist = new TH2F("hist", vartitle, 100, xvarmin, xvarmax, 100, yvarmin, yvarmax);

    mychain->Project("hist", myvar);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    /* hist->SetMaximum(1.3*hist->GetMaximum()); */
    hist->SetMinimum(0);
    hist->Draw("COLZ");
    CMS_lumi(canvas);
	canvas->Update();
	canvas->RedrawAxis();
    canvas->SaveAs("./figures/fit_results/"+figpath+".png");
    delete hist;
}

int fit_results() {
    setTDRStyle();

    /* TChain *mychain1 = new TChain("PackedGenParticle/Events"); */
    /* mychain1->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked1/GenPacked1.root"); */
    TChain *mychain2 = new TChain("PackedGenParticle/Events");
    mychain2->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked2/GenPacked2.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/fit_results")) gSystem->MakeDirectory("./figures/fit_results");
   
    /* draw_fit_results_1d(mychain1, "dR_Kp_vec", ";#Delta#it{R(K^{+})};# events", 0, 10, "dR_Kp"); */
    /* draw_fit_results_1d(mychain1, "dR_Km_vec", ";#Delta#it{R(K^{-})};# events", 0, 10, "dR_Km"); */
    /* draw_fit_results_1d(mychain1, "dR_pi_vec", ";#Delta#it{R(#pi^{+})};# events", 0, 10, "dR_pi"); */
    /* draw_fit_results_1d(mychain1, "dR_Kp_vec", ";#Delta#it{R(K^{+})};# events", 0, 0.1, "dR_Kp_small"); */
    /* draw_fit_results_1d(mychain1, "dR_Km_vec", ";#Delta#it{R(K^{-})};# events", 0, 0.1, "dR_Km_small"); */
    /* draw_fit_results_1d(mychain1, "dR_pi_vec", ";#Delta#it{R(#pi^{+})};# events", 0, 0.1, "dR_pi_small"); */
   
    /* draw_fit_results_1d(mychain2, "Gen_pt_Kp", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 60, "Gen_pt_Kp"); */
    /* draw_fit_results_1d(mychain2, "Gen_pt_Km", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 60, "Gen_pt_Km"); */
    /* draw_fit_results_1d(mychain2, "Gen_pt_pi", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 80, "Gen_pt_pi"); */
    /* draw_fit_results_1d(mychain2, "Gen_pt_Kp", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 5, "Gen_pt_Kp_small"); */
    /* draw_fit_results_1d(mychain2, "Gen_pt_Km", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 5, "Gen_pt_Km_small"); */
    /* draw_fit_results_1d(mychain2, "Gen_pt_pi", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "Gen_pt_pi_small"); */

    /* draw_fit_results_1d(mychain2, "pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 80, "pt_Kp"); */
    /* draw_fit_results_1d(mychain2, "pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 80, "pt_Km"); */
    /* draw_fit_results_1d(mychain2, "pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 100, "pt_pi"); */
    /* draw_fit_results_1d(mychain2, "pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 5, "pt_Kp_small"); */
    /* draw_fit_results_1d(mychain2, "pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 5, "pt_Km_small"); */
    /* draw_fit_results_1d(mychain2, "pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "pt_pi_small"); */

    /* draw_fit_results_1d(mychain2, "Gen_p_Kp", ";Gen #it{p(K^{+})} [GeV];# events", 0, 60, "Gen_p_Kp"); */
    /* draw_fit_results_1d(mychain2, "Gen_p_Km", ";Gen #it{p(K^{-})} [GeV];# events", 0, 60, "Gen_p_Km"); */
    /* draw_fit_results_1d(mychain2, "Gen_p_pi", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 80, "Gen_p_pi"); */
    /* draw_fit_results_1d(mychain2, "Gen_p_Kp", ";Gen #it{p(K^{+})} [GeV];# events", 0, 5, "Gen_p_Kp_small"); */
    /* draw_fit_results_1d(mychain2, "Gen_p_Km", ";Gen #it{p(K^{-})} [GeV];# events", 0, 5, "Gen_p_Km_small"); */
    /* draw_fit_results_1d(mychain2, "Gen_p_pi", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 5, "Gen_p_pi_small"); */

    /* draw_fit_results_1d(mychain2, "p_Kp", ";#it{p(K^{+})} [GeV];# events", 0, 60, "p_Kp"); */
    /* draw_fit_results_1d(mychain2, "p_Km", ";#it{p(K^{-})} [GeV];# events", 0, 60, "p_Km"); */
    /* draw_fit_results_1d(mychain2, "p_pi", ";#it{p(#pi^{+})} [GeV];# events", 0, 80, "p_pi"); */
    /* draw_fit_results_1d(mychain2, "p_Kp", ";#it{p(K^{+})} [GeV];# events", 0, 5, "p_Kp_small"); */
    /* draw_fit_results_1d(mychain2, "p_Km", ";#it{p(K^{-})} [GeV];# events", 0, 5, "p_Km_small"); */
    /* draw_fit_results_1d(mychain2, "p_pi", ";#it{p(#pi^{+})} [GeV];# events", 0, 5, "p_pi_small"); */

    draw_fit_results_1d(mychain2, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 20, "all_pt");
    draw_fit_results_1d(mychain2, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 5, "all_pt_small");
    draw_fit_results_1d(mychain2, "all_p", ";track #it{p} [GeV];# entries", 0, 40, "all_p");
    draw_fit_results_1d(mychain2, "all_p", ";track #it{p} [GeV];# entries", 0, 5, "all_p_small");
    
    /* draw_fit_results_1d(mychain2, "dR_Kp_Km", ";#Delta#it{R(K^{+},K^{-})};# events", 0, 0.2, "dR_Kp_Km"); */
    /* draw_fit_results_1d(mychain2, "dR_Kp_pi", ";#Delta#it{R(K^{+},#pi^{+})};# events", 0, 0.8, "dR_Kp_pi"); */
    /* draw_fit_results_1d(mychain2, "dR_Km_pi", ";#Delta#it{R(K^{-},#pi^{+})};# events", 0, 0.8, "dR_Km_pi"); */

    /* draw_fit_results_1d(mychain2, "m_phi", ";#it{M(K^{+}K^{-})} [GeV];# events", 0.99, 1.05, "m_phi"); */
    /* draw_fit_results_1d(mychain2, "m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# events", 1.85, 2.1, "m_Ds"); */
    /* draw_fit_results_2d(mychain2, "m_Ds:m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{M(K^{+}K^{-}#pi^{+})} [GeV]", 0.99, 1.05, 1.85, 2.1, "m_phi_Ds"); */

    /* draw_fit_results_1d(mychain2, "chi2_phi", ";#it{#chi^{2}(#phi)};# events", 0, 10, "chi2_phi"); */
    /* draw_fit_results_1d(mychain2, "chi2_Ds", ";#it{#chi^{2}(D_{s}^{+})};# events", 0, 25, "chi2_Ds"); */
    /* draw_fit_results_2d(mychain2, "chi2_phi:m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{#chi^{2}(#phi)}", 0.99, 1.05, 0, 10, "m_chi2_phi"); */
    /* draw_fit_results_2d(mychain2, "chi2_Ds:m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];#it{#chi^{2}(D_{s}^{+})}", 1.85, 2.1, 0, 25, "m_chi2_Ds"); */

    /* draw_fit_results_1d(mychain2, "d_Kp_Km", ";#Delta#it{d(K^{+},K^{-})} [cm];# events", 0, 0.5, "d_Kp_Km"); */
    /* draw_fit_results_1d(mychain2, "d_Kp_pi", ";#Delta#it{d(K^{+},#pi^{+})} [cm];# events", 0, 0.5, "d_Kp_pi"); */
    /* draw_fit_results_1d(mychain2, "d_Km_pi", ";#Delta#it{d(K^{-},#pi^{+})} [cm];# events", 0, 0.5, "d_Km_pi"); */

    /* draw_fit_results_1d(mychain2, "d_phi_pi", ";#Delta#it{#phi,#pi^{+})} [cm];# events", 0, 5, "d_phi_pi"); */
    /* draw_fit_results_1d(mychain2, "d_phi_Ds", ";#Delta#it{#phi,D_{s}^{+})} [cm];# events", 0, 5, "d_phi_Ds"); */

    /* delete mychain1; */
    delete mychain2;

    return 0;
}

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

    TChain *mychain = new TChain("PackedGenParticle/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked/GenPacked.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/fit_results")) gSystem->MakeDirectory("./figures/fit_results");

    /* draw_fit_results_1d(mychain, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 20, "all_pt"); */
    /* draw_fit_results_1d(mychain, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 5, "all_pt_small"); */
    /* draw_fit_results_1d(mychain, "all_p", ";track #it{p} [GeV];# entries", 0, 40, "all_p"); */
    /* draw_fit_results_1d(mychain, "all_p", ";track #it{p} [GeV];# entries", 0, 5, "all_p_small"); */

    /* draw_fit_results_1d(mychain, "dR_Kp", ";#Delta#it{R(K^{+})};# events", 0, 10, "dR_Kp"); */
    /* draw_fit_results_1d(mychain, "dR_Km", ";#Delta#it{R(K^{-})};# events", 0, 10, "dR_Km"); */
    /* draw_fit_results_1d(mychain, "dR_pi", ";#Delta#it{R(#pi^{+})};# events", 0, 10, "dR_pi"); */
    /* draw_fit_results_1d(mychain, "dR_Kp", ";#Delta#it{R(K^{+})};# events", 0, 0.1, "dR_Kp_small"); */
    /* draw_fit_results_1d(mychain, "dR_Km", ";#Delta#it{R(K^{-})};# events", 0, 0.1, "dR_Km_small"); */
    /* draw_fit_results_1d(mychain, "dR_pi", ";#Delta#it{R(#pi^{+})};# events", 0, 0.1, "dR_pi_small"); */

    /* draw_fit_results_1d(mychain, "Gen_p_Kp", ";Gen #it{p(K^{+})} [GeV];# events", 0, 60, "Gen_p_Kp"); */
    /* draw_fit_results_1d(mychain, "Gen_p_Km", ";Gen #it{p(K^{-})} [GeV];# events", 0, 60, "Gen_p_Km"); */
    /* draw_fit_results_1d(mychain, "Gen_p_pi", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 80, "Gen_p_pi"); */
    /* draw_fit_results_1d(mychain, "Gen_p_Kp", ";Gen #it{p(K^{+})} [GeV];# events", 0, 5, "Gen_p_Kp_small"); */
    /* draw_fit_results_1d(mychain, "Gen_p_Km", ";Gen #it{p(K^{-})} [GeV];# events", 0, 5, "Gen_p_Km_small"); */
    /* draw_fit_results_1d(mychain, "Gen_p_pi", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 5, "Gen_p_pi_small"); */

    /* draw_fit_results_1d(mychain, "Gen_pt_Kp", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 60, "Gen_pt_Kp"); */
    /* draw_fit_results_1d(mychain, "Gen_pt_Km", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 60, "Gen_pt_Km"); */
    /* draw_fit_results_1d(mychain, "Gen_pt_pi", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 80, "Gen_pt_pi"); */
    /* draw_fit_results_1d(mychain, "Gen_pt_Kp", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 5, "Gen_pt_Kp_small"); */
    /* draw_fit_results_1d(mychain, "Gen_pt_Km", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 5, "Gen_pt_Km_small"); */
    /* draw_fit_results_1d(mychain, "Gen_pt_pi", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "Gen_pt_pi_small"); */

    /* draw_fit_results_1d(mychain, "Gen_eta_Kp", ";Gen #it{#eta(K^{+})};# events", -6, 6, "Gen_eta_Kp"); */
    /* draw_fit_results_1d(mychain, "Gen_eta_Km", ";Gen #it{#eta(K^{-})};# events", -6, 6, "Gen_eta_Km"); */
    /* draw_fit_results_1d(mychain, "Gen_eta_pi", ";Gen #it{#eta(#pi^{+})};# events", -6, 6, "Gen_eta_pi"); */
    /* draw_fit_results_1d(mychain, "Gen_phi_Kp", ";Gen #it{#phi(K^{+})};# events", -3.15, 3.15, "Gen_phi_Kp"); */
    /* draw_fit_results_1d(mychain, "Gen_phi_Km", ";Gen #it{#phi(K^{-})};# events", -3.15, 3.15, "Gen_phi_Km"); */
    /* draw_fit_results_1d(mychain, "Gen_phi_pi", ";Gen #it{#phi(#pi^{+})};# events", -3.15, 3.15, "Gen_phi_pi"); */

    /* draw_fit_results_1d(mychain, "p_Kp", ";#it{p(K^{+})} [GeV];# events", 0, 60, "p_Kp"); */
    /* draw_fit_results_1d(mychain, "p_Km", ";#it{p(K^{-})} [GeV];# events", 0, 60, "p_Km"); */
    /* draw_fit_results_1d(mychain, "p_pi", ";#it{p(#pi^{+})} [GeV];# events", 0, 80, "p_pi"); */
    /* draw_fit_results_1d(mychain, "p_Kp", ";#it{p(K^{+})} [GeV];# events", 0, 5, "p_Kp_small"); */
    /* draw_fit_results_1d(mychain, "p_Km", ";#it{p(K^{-})} [GeV];# events", 0, 5, "p_Km_small"); */
    /* draw_fit_results_1d(mychain, "p_pi", ";#it{p(#pi^{+})} [GeV];# events", 0, 5, "p_pi_small"); */

    /* draw_fit_results_1d(mychain, "pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 60, "pt_Kp"); */
    /* draw_fit_results_1d(mychain, "pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 60, "pt_Km"); */
    /* draw_fit_results_1d(mychain, "pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 80, "pt_pi"); */
    /* draw_fit_results_1d(mychain, "pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 5, "pt_Kp_small"); */
    /* draw_fit_results_1d(mychain, "pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 5, "pt_Km_small"); */
    /* draw_fit_results_1d(mychain, "pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "pt_pi_small"); */

    draw_fit_results_1d(mychain, "pt_Kp-pt_Km", ";#it{p_{T}(K^{+})-p_{T}(K^{-})} [GeV];# events", -30, 30, "pt_Kp_minus_pt_Km");
    draw_fit_results_2d(mychain, "pt_Km:pt_Kp", ";#it{p_{T}(K^{+})} [GeV];#it{p_{T}(K^{-})} [GeV]", 0, 60, 0, 60, "pt_Kp_Km");

    /* draw_fit_results_1d(mychain, "eta_Kp", ";#it{#eta(K^{+})};# events", -6, 6, "eta_Kp"); */
    /* draw_fit_results_1d(mychain, "eta_Km", ";#it{#eta(K^{-})};# events", -6, 6, "eta_Km"); */
    /* draw_fit_results_1d(mychain, "eta_pi", ";#it{#eta(#pi^{+})};# events", -6, 6, "eta_pi"); */
    /* draw_fit_results_1d(mychain, "phi_Kp", ";#it{#phi(K^{+})};# events", -3.15, 3.15, "phi_Kp"); */
    /* draw_fit_results_1d(mychain, "phi_Km", ";#it{#phi(K^{-})};# events", -3.15, 3.15, "phi_Km"); */
    /* draw_fit_results_1d(mychain, "phi_pi", ";#it{#phi(#pi^{+})};# events", -3.15, 3.15, "phi_pi"); */

    /* draw_fit_results_1d(mychain, "Gen_dR_Kp_Km", ";Gen #Delta#it{R(K^{+},K^{-})};# events", 0, 0.2, "Gen_dR_Kp_Km"); */
    /* draw_fit_results_1d(mychain, "Gen_dR_Kp_pi", ";Gen #Delta#it{R(K^{+},#pi^{+})};# events", 0, 0.8, "Gen_dR_Kp_pi"); */
    /* draw_fit_results_1d(mychain, "Gen_dR_Km_pi", ";Gen #Delta#it{R(K^{-},#pi^{+})};# events", 0, 0.8, "Gen_dR_Km_pi"); */

    /* draw_fit_results_1d(mychain, "dR_Kp_Km", ";#Delta#it{R(K^{+},K^{-})};# events", 0, 0.2, "dR_Kp_Km"); */
    /* draw_fit_results_1d(mychain, "dR_Kp_pi", ";#Delta#it{R(K^{+},#pi^{+})};# events", 0, 0.8, "dR_Kp_pi"); */
    /* draw_fit_results_1d(mychain, "dR_Km_pi", ";#Delta#it{R(K^{-},#pi^{+})};# events", 0, 0.8, "dR_Km_pi"); */

    /* draw_fit_results_1d(mychain, "d_Kp_Km", ";#Delta#it{d(K^{+},K^{-})} [cm];# events", 0, 0.5, "d_Kp_Km"); */
    /* draw_fit_results_1d(mychain, "d_Kp_pi", ";#Delta#it{d(K^{+},#pi^{+})} [cm];# events", 0, 0.5, "d_Kp_pi"); */
    /* draw_fit_results_1d(mychain, "d_Km_pi", ";#Delta#it{d(K^{-},#pi^{+})} [cm];# events", 0, 0.5, "d_Km_pi"); */

    /* draw_fit_results_1d(mychain, "chi2_phi", ";#it{#chi^{2}(#phi)};# events", 0, 10, "chi2_phi"); */
    /* draw_fit_results_1d(mychain, "refitted_pt_phi", ";#it{p_{T}(K^{+}K^{-})} [GeV];# events", 0, 150, "refitted_pt_phi"); */
    /* draw_fit_results_1d(mychain, "refitted_p_phi", ";#it{p(K^{+}K^{-})} [GeV];# events", 0, 300, "refitted_p_phi"); */
    /* draw_fit_results_1d(mychain, "refitted_m_phi", ";#it{M(K^{+}K^{-})} [GeV];# events", 0.99, 1.05, "refitted_m_phi"); */
    draw_fit_results_1d(mychain, "dR_phi_pi", ";#DeltaR(#it{#phi,#pi^{+})};# events", 0, 1, "dR_phi_pi");
    /* draw_fit_results_1d(mychain, "d_phi_pi", ";#it{d(#phi,#pi^{+})} [cm];# events", 0, 8, "d_phi_pi"); */

    draw_fit_results_2d(mychain, "pt_pi:refitted_pt_phi", ";#it{p_{T}(#phi)} [GeV];#it{p_{T}(#pi^{+})} [GeV]", 0, 80, 0, 60, "pt_phi_pi");

    /* draw_fit_results_1d(mychain, "chi2_Ds", ";#it{#chi^{2}(D_{s}^{+})};# events", 0, 25, "chi2_Ds"); */
    /* draw_fit_results_1d(mychain, "refitted2_pt_Ds", ";#it{p_{T}(K^{+}K^{-}#pi^{+})} [GeV];# events", 0, 150, "refitted2_pt_Ds"); */
    /* draw_fit_results_1d(mychain, "refitted2_p_Ds", ";#it{p(K^{+}K^{-}#pi^{+})} [GeV];# events", 0, 500, "refitted2_p_Ds"); */
    /* draw_fit_results_1d(mychain, "refitted2_m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# events", 1.85, 2.1, "refitted2_m_Ds"); */
    draw_fit_results_1d(mychain, "dR_phi_Ds", ";#DeltaR(#it{#phi,D_{s}^{+})};# events", 0, 0.5, "dR_phi_Ds");
    /* draw_fit_results_1d(mychain, "d_phi_Ds", ";#it{d(#phi,D_{s}^{+})} [cm];# events", 0, 8, "d_phi_Ds"); */

    /* draw_fit_results_2d(mychain, "refitted2_m_Ds:refitted_m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{M(K^{+}K^{-}#pi^{+})} [GeV]", 0.99, 1.05, 1.85, 2.1, "m_phi_Ds"); */
    /* draw_fit_results_2d(mychain, "chi2_phi:refitted_m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{#chi^{2}(#phi)}", 0.99, 1.05, 0, 10, "m_chi2_phi"); */
    /* draw_fit_results_2d(mychain, "chi2_Ds:refitted2_m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];#it{#chi^{2}(D_{s}^{+})}", 1.85, 2.1, 0, 25, "m_chi2_Ds"); */

    delete mychain;

    return 0;
}

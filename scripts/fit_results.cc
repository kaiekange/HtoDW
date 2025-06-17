#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath){

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

void draw_2d(TChain *mychain, TString myvar, TString vartitle, float xvarmin, float xvarmax, float yvarmin, float yvarmax, TString figpath){

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

    TChain *mychain = new TChain("SelectionStudy/Events");
    /* mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenTrackMatch/GenTrackMatch.root"); */
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/fit_results")) gSystem->MakeDirectory("./figures/fit_results");

    /* draw_1d(mychain, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 20, "all_pt"); */
    /* draw_1d(mychain, "all_pt", ";track #it{p_{T}} [GeV];# entries", 0, 5, "all_pt_small"); */
    /* draw_1d(mychain, "all_p", ";track #it{p} [GeV];# entries", 0, 40, "all_p"); */
    /* draw_1d(mychain, "all_p", ";track #it{p} [GeV];# entries", 0, 5, "all_p_small"); */

    /* draw_1d(mychain, "dR_Kp", ";#Delta#it{R(K^{+})};# events", 0, 10, "dR_Kp"); */
    /* draw_1d(mychain, "dR_Km", ";#Delta#it{R(K^{-})};# events", 0, 10, "dR_Km"); */
    /* draw_1d(mychain, "dR_pi", ";#Delta#it{R(#pi^{+})};# events", 0, 10, "dR_pi"); */
    /* draw_1d(mychain, "dR_Kp", ";#Delta#it{R(K^{+})};# events", 0, 0.1, "dR_Kp_small"); */
    /* draw_1d(mychain, "dR_Km", ";#Delta#it{R(K^{-})};# events", 0, 0.1, "dR_Km_small"); */
    /* draw_1d(mychain, "dR_pi", ";#Delta#it{R(#pi^{+})};# events", 0, 0.1, "dR_pi_small"); */

    /* draw_1d(mychain, "Gen_Kp_P", ";Gen #it{p(K^{+})} [GeV];# events", 0, 60, "Gen_Kp_P"); */
    /* draw_1d(mychain, "Gen_Km_P", ";Gen #it{p(K^{-})} [GeV];# events", 0, 60, "Gen_Km_P"); */
    /* draw_1d(mychain, "Gen_pi_P", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 80, "Gen_pi_P"); */
    /* draw_1d(mychain, "Gen_Kp_P", ";Gen #it{p(K^{+})} [GeV];# events", 0, 5, "Gen_Kp_small_P"); */
    /* draw_1d(mychain, "Gen_Km_P", ";Gen #it{p(K^{-})} [GeV];# events", 0, 5, "Gen_Km_small_P"); */
    /* draw_1d(mychain, "Gen_pi_P", ";Gen #it{p(#pi^{+})} [GeV];# events", 0, 5, "Gen_pi_small_P"); */

    /* draw_1d(mychain, "Gen_Kp_PT", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 60, "Gen_Kp_PT"); */
    /* draw_1d(mychain, "Gen_Km_PT", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 60, "Gen_Km_PT"); */
    /* draw_1d(mychain, "Gen_pi_PT", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 80, "Gen_pi_PT"); */
    /* draw_1d(mychain, "Gen_Kp_PT", ";Gen #it{p_{T}(K^{+})} [GeV];# events", 0, 5, "Gen_small_PT_Kp"); */
    /* draw_1d(mychain, "Gen_Km_PT", ";Gen #it{p_{T}(K^{-})} [GeV];# events", 0, 5, "Gen_small_PT_Km"); */
    /* draw_1d(mychain, "Gen_pi_PT", ";Gen #it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "Gen_small_PT_pi"); */

    /* draw_1d(mychain, "Gen_Kp_ETA", ";Gen #it{#eta(K^{+})};# events", -6, 6, "Gen_Kp_ETA"); */
    /* draw_1d(mychain, "Gen_Km_ETA", ";Gen #it{#eta(K^{-})};# events", -6, 6, "Gen_Km_ETA"); */
    /* draw_1d(mychain, "Gen_pi_ETA", ";Gen #it{#eta(#pi^{+})};# events", -6, 6, "Gen_pi_ETA"); */
    /* draw_1d(mychain, "Gen_Kp_PHI", ";Gen #it{#phi(K^{+})};# events", -3.15, 3.15, "Gen_Kp_PHI"); */
    /* draw_1d(mychain, "Gen_Km_PHI", ";Gen #it{#phi(K^{-})};# events", -3.15, 3.15, "Gen_Km_PHI"); */
    /* draw_1d(mychain, "Gen_pi_PHI", ";Gen #it{#phi(#pi^{+})};# events", -3.15, 3.15, "Gen_pi_PHI"); */

    /* draw_1d(mychain, "Kp_P", ";#it{p(K^{+})} [GeV];# events", 0, 60, "Kp_P"); */
    /* draw_1d(mychain, "Km_P", ";#it{p(K^{-})} [GeV];# events", 0, 60, "Km_P"); */
    /* draw_1d(mychain, "pi_P", ";#it{p(#pi^{+})} [GeV];# events", 0, 80, "pi_P"); */
    /* draw_1d(mychain, "Kp_P", ";#it{p(K^{+})} [GeV];# events", 0, 5, "Kp_small_P"); */
    /* draw_1d(mychain, "Km_P", ";#it{p(K^{-})} [GeV];# events", 0, 5, "Km_small_P"); */
    /* draw_1d(mychain, "pi_P", ";#it{p(#pi^{+})} [GeV];# events", 0, 5, "pi_small_P"); */
    
    /* draw_1d(mychain, "Kp_PT", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 60, "Kp_PT"); */
    /* draw_1d(mychain, "Km_PT", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 60, "Km_PT"); */
    /* draw_1d(mychain, "pi_PT", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 80, "pi_PT"); */
    /* draw_1d(mychain, "Kp_PT", ";#it{p_{T}(K^{+})} [GeV];# events", 0, 5, "small_PT_Kp"); */
    /* draw_1d(mychain, "Km_PT", ";#it{p_{T}(K^{-})} [GeV];# events", 0, 5, "small_PT_Km"); */
    /* draw_1d(mychain, "pi_PT", ";#it{p_{T}(#pi^{+})} [GeV];# events", 0, 5, "small_PT_pi"); */

    /* draw_1d(mychain, "Kp_ETA", ";#it{#eta(K^{+})};# events", -6, 6, "Kp_ETA"); */
    /* draw_1d(mychain, "Km_ETA", ";#it{#eta(K^{-})};# events", -6, 6, "Km_ETA"); */
    /* draw_1d(mychain, "pi_ETA", ";#it{#eta(#pi^{+})};# events", -6, 6, "pi_ETA"); */
    /* draw_1d(mychain, "Kp_PHI", ";#it{#phi(K^{+})};# events", -3.15, 3.15, "Kp_PHI"); */
    /* draw_1d(mychain, "Km_PHI", ";#it{#phi(K^{-})};# events", -3.15, 3.15, "Km_PHI"); */
    /* draw_1d(mychain, "pi_PHI", ";#it{#phi(#pi^{+})};# events", -3.15, 3.15, "pi_PHI"); */

    /* draw_1d(mychain, "Gen_Kp_PT-Gen_Km_PT", ";Gen #it{p_{T}(K^{+})-p_{T}(K^{-})} [GeV];# events", -30, 30, "Gen_Kp_PT_minus_Km_PT"); */
    /* draw_2d(mychain, "Gen_Km_PT:Gen_Kp_PT", ";Gen #it{p_{T}(K^{+})} [GeV];Gen #it{p_{T}(K^{-})} [GeV]", 0, 60, 0, 60, "Gen_Kp_Km_PT"); */
    /* draw_1d(mychain, "Gen_Kp_PT/Gen_Km_PT", ";Gen #it{p_{T}(K^{+})/p_{T}(K^{-})} [GeV];# events", 0, 5, "Gen_Kp_PT_divide_Km_PT"); */
    /* draw_1d(mychain, "Kp_PT-Km_PT", ";#it{p_{T}(K^{+})-p_{T}(K^{-})} [GeV];# events", -30, 30, "Kp_PT_minus_Km_PT"); */
    /* draw_2d(mychain, "Km_PT:Kp_PT", ";#it{p_{T}(K^{+})} [GeV];#it{p_{T}(K^{-})} [GeV]", 0, 60, 0, 60, "Kp_Km_PT"); */
    /* draw_1d(mychain, "Kp_PT/Km_PT", ";#it{p_{T}(K^{+})/p_{T}(K^{-})} [GeV];# events", 0, 5, "Kp_PT_divide_Km_PT"); */
    
    /* draw_1d(mychain, "Gen_phi_PT-Gen_pi_PT", ";Gen #it{p_{T}(#phi)-p_{T}(#pi^{+})} [GeV];# events", -50, 80, "Gen_phi_PT_minus_pi_PT"); */
    /* draw_2d(mychain, "Gen_pi_PT:Gen_phi_PT", ";Gen #it{p_{T}(#phi)} [GeV];Gen #it{p_{T}(#pi^{+})} [GeV]", 0, 80, 0, 60, "Gen_phi_pi_PT"); */
    /* draw_1d(mychain, "Gen_phi_PT/Gen_pi_PT", ";Gen #it{p_{T}(#phi)/p_{T}(#pi^{+})} [GeV];# events", 0, 20, "Gen_phi_PT_divide_pi_PT"); */
    /* draw_1d(mychain, "phi_PT-pi_PT", ";#it{p_{T}(#phi)-p_{T}(#pi^{+})} [GeV];# events", -50, 80, "phi_PT_minus_pi_PT"); */
    /* draw_2d(mychain, "pi_PT:phi_PT", ";#it{p_{T}(#phi)} [GeV];#it{p_{T}(#pi^{+})} [GeV]", 0, 80, 0, 60, "phi_pi_PT"); */
    /* draw_1d(mychain, "phi_PT/pi_PT", ";#it{p_{T}(#phi)/p_{T}(#pi^{+})} [GeV];# events", 0, 20, "phi_PT_divide_pi_PT"); */
    
    /* draw_2d(mychain, "Gen_Kp_PP:((Gen_Kp_PL - Gen_Km_PL)/(Gen_Kp_PL + Gen_Km_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "Gen_APplot_phi"); */
    /* draw_2d(mychain, "Kp_PP:((Kp_PL - Km_PL)/(Kp_PL + Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "APplot_phi"); */
    /* draw_2d(mychain, "phiFit_Kp_PP:((phiFit_Kp_PL - phiFit_Km_PL)/(phiFit_Kp_PL + phiFit_Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "phiFit_APplot_phi"); */
    
    /* draw_2d(mychain, "Gen_phi_PP:((Gen_phi_PL - Gen_pi_PL)/(Gen_phi_PL + Gen_pi_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "Gen_APplot_Ds"); */
    /* draw_2d(mychain, "phi_PP:((phi_PL - pi_PL)/(phi_PL + pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "APplot_Ds"); */
    /* draw_2d(mychain, "DsFit_phi_PP:((DsFit_phi_PL - DsFit_pi_PL)/(DsFit_phi_PL + DsFit_pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "DsFit_APplot_Ds"); */

    /* draw_2d(mychain, "abs(phiFit_ENDVX_XERR/phiFit_ENDVX_X):(phiFit_ENDVX_X - Gen_Kp_ORIVX_X)/Gen_Kp_ORIVX_X", ";#varepsilonvtx_{#it{x}}(#it{#phi});|#Deltavtx_{#it{x}}(#it{#phi})/vtx_{#it{x}}(#it{#phi})|", -3, 3, 0, 1, "phi_sig_vtx_x"); */
    /* draw_2d(mychain, "abs(phiFit_ENDVX_YERR/phiFit_ENDVX_Y):(phiFit_ENDVX_Y - Gen_Kp_ORIVX_Y)/Gen_Kp_ORIVX_Y", ";#varepsilonvtx_{#it{y}}(#it{#phi});|#Deltavtx_{#it{y}}(#it{#phi})/vtx_{#it{y}}(#it{#phi})|", -3, 3, 0, 1, "phi_sig_vtx_y"); */
    /* draw_2d(mychain, "abs(phiFit_ENDVX_ZERR/phiFit_ENDVX_Z):(phiFit_ENDVX_Z - Gen_Kp_ORIVX_Z)/Gen_Kp_ORIVX_Z", ";#varepsilonvtx_{#it{z}}(#it{#phi});|#Deltavtx_{#it{z}}(#it{#phi})/vtx_{#it{z}}(#it{#phi})|", -3, 3, 0, 1, "phi_sig_vtx_z"); */
    
    /* draw_2d(mychain, "abs(DsFit_ENDVX_XERR/DsFit_ENDVX_X):(DsFit_ENDVX_X - Gen_Kp_ORIVX_X)/Gen_Kp_ORIVX_X", ";#varepsilonvtx_{#it{x}}(#it{D_{s}^{+}});|#Deltavtx_{#it{x}}(#it{D_{s}^{+}})/vtx_{#it{x}}(#it{D_{s}^{+}})|", -3, 3, 0, 1, "Ds_sig_vtx_x"); */
    /* draw_2d(mychain, "abs(DsFit_ENDVX_YERR/DsFit_ENDVX_Y):(DsFit_ENDVX_Y - Gen_Kp_ORIVX_Y)/Gen_Kp_ORIVX_Y", ";#varepsilonvtx_{#it{y}}(#it{D_{s}^{+}});|#Deltavtx_{#it{y}}(#it{D_{s}^{+}})/vtx_{#it{y}}(#it{D_{s}^{+}})|", -3, 3, 0, 1, "Ds_sig_vtx_y"); */
    /* draw_2d(mychain, "abs(DsFit_ENDVX_ZERR/DsFit_ENDVX_Z):(DsFit_ENDVX_Z - Gen_Kp_ORIVX_Z)/Gen_Kp_ORIVX_Z", ";#varepsilonvtx_{#it{z}}(#it{D_{s}^{+}});|#Deltavtx_{#it{z}}(#it{D_{s}^{+}})/vtx_{#it{z}}(#it{D_{s}^{+}})|", -3, 3, 0, 1, "Ds_sig_vtx_z"); */
    
    draw_2d(mychain, "phiFit_ENDVX_XERR:(phiFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";fit vtx_{#it{x}}(#it{#phi})-Gen vtx_{#it{x}}(#it{#phi}) [cm];fit #Deltavtx_{#it{x}}(#it{#phi})|", -1, 1, 0, 0.6, "phi_err_vtx_x");
    draw_2d(mychain, "phiFit_ENDVX_YERR:(phiFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";fit vtx_{#it{y}}(#it{#phi})-Gen vtx_{#it{y}}(#it{#phi}) [cm];fit #Deltavtx_{#it{y}}(#it{#phi})|", -1, 1, 0, 0.6, "phi_err_vtx_y");
    draw_2d(mychain, "phiFit_ENDVX_ZERR:(phiFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";fit vtx_{#it{z}}(#it{#phi})-Gen vtx_{#it{z}}(#it{#phi}) [cm];fit #Deltavtx_{#it{z}}(#it{#phi})|", -2, 2, 0, 2, "phi_err_vtx_z");
    
    draw_2d(mychain, "DsFit_ENDVX_XERR:(DsFit_ENDVX_X - Gen_Kp_ORIVX_X)", ";fit vtx_{#it{x}}(#it{D_{s}^{+}})-Gen vtx_{#it{x}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{x}}(#it{D_{s}^{+}})|", -0.3, 0.3, 0, 0.1, "Ds_err_vtx_x");
    draw_2d(mychain, "DsFit_ENDVX_YERR:(DsFit_ENDVX_Y - Gen_Kp_ORIVX_Y)", ";fit vtx_{#it{y}}(#it{D_{s}^{+}})-Gen vtx_{#it{y}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{y}}(#it{D_{s}^{+}})|", -0.3, 0.3, 0, 0.1, "Ds_err_vtx_y");
    draw_2d(mychain, "DsFit_ENDVX_ZERR:(DsFit_ENDVX_Z - Gen_Kp_ORIVX_Z)", ";fit vtx_{#it{z}}(#it{D_{s}^{+}})-Gen vtx_{#it{z}}(#it{D_{s}^{+}}) [cm];fit #Deltavtx_{#it{z}}(#it{D_{s}^{+}})|", -0.5, 0.5, 0, 0.3, "Ds_err_vtx_z");


    /* draw_1d(mychain, "Gen_dR_Kp_Km", ";Gen #Delta#it{R(K^{+},K^{-})};# events", 0, 0.2, "Gen_dR_Kp_Km"); */
    /* draw_1d(mychain, "Gen_dR_Kp_pi", ";Gen #Delta#it{R(K^{+},#pi^{+})};# events", 0, 0.8, "Gen_dR_Kp_pi"); */
    /* draw_1d(mychain, "Gen_dR_Km_pi", ";Gen #Delta#it{R(K^{-},#pi^{+})};# events", 0, 0.8, "Gen_dR_Km_pi"); */

    /* draw_1d(mychain, "dR_Kp_Km", ";#Delta#it{R(K^{+},K^{-})};# events", 0, 0.2, "dR_Kp_Km"); */
    /* draw_1d(mychain, "dR_Kp_pi", ";#Delta#it{R(K^{+},#pi^{+})};# events", 0, 0.8, "dR_Kp_pi"); */
    /* draw_1d(mychain, "dR_Km_pi", ";#Delta#it{R(K^{-},#pi^{+})};# events", 0, 0.8, "dR_Km_pi"); */

    /* draw_1d(mychain, "d_Kp_Km", ";#Delta#it{d(K^{+},K^{-})} [cm];# events", 0, 0.5, "d_Kp_Km"); */
    /* draw_1d(mychain, "d_Kp_pi", ";#Delta#it{d(K^{+},#pi^{+})} [cm];# events", 0, 0.5, "d_Kp_pi"); */
    /* draw_1d(mychain, "d_Km_pi", ";#Delta#it{d(K^{-},#pi^{+})} [cm];# events", 0, 0.5, "d_Km_pi"); */

    /* draw_1d(mychain, "phi_CHI2", ";#it{#chi^{2}(#phi)};# events", 0, 10, "phi_CHI2"); */
    /* draw_1d(mychain, "phi_PT", ";#it{p_{T}(K^{+}K^{-})} [GeV];# events", 0, 150, "phi_PT"); */
    /* draw_1d(mychain, "phi_P", ";#it{p(K^{+}K^{-})} [GeV];# events", 0, 300, "phi_P"); */
    /* draw_1d(mychain, "phi_M", ";#it{M(K^{+}K^{-})} [GeV];# events", 0.99, 1.05, "phi_M"); */
    /* draw_1d(mychain, "dR_phi_pi", ";#DeltaR(#it{#phi,#pi^{+})};# events", 0, 1, "dR_phi_pi"); */
    /* draw_1d(mychain, "d_phi_pi", ";#it{d(#phi,#pi^{+})} [cm];# events", 0, 8, "d_phi_pi"); */


    /* draw_1d(mychain, "Ds_CHI2", ";#it{#chi^{2}(D_{s}^{+})};# events", 0, 25, "Ds_CHI2"); */
    /* draw_1d(mychain, "Ds_PT", ";#it{p_{T}(K^{+}K^{-}#pi^{+})} [GeV];# events", 0, 150, "Ds_pt"); */
    /* draw_1d(mychain, "Ds_P", ";#it{p(K^{+}K^{-}#pi^{+})} [GeV];# events", 0, 500, "Ds_p"); */
    /* draw_1d(mychain, "Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# events", 1.85, 2.1, "Ds_M"); */
    /* draw_1d(mychain, "dR_phi_Ds", ";#DeltaR(#it{#phi,D_{s}^{+})};# events", 0, 0.5, "dR_phi_Ds"); */
    /* draw_1d(mychain, "d_phi_Ds", ";#it{d(#phi,D_{s}^{+})} [cm];# events", 0, 8, "d_phi_Ds"); */

    /* draw_2d(mychain, "Ds_M:phi_M", ";#it{M(K^{+}K^{-})} [GeV];#it{M(K^{+}K^{-}#pi^{+})} [GeV]", 0.99, 1.05, 1.85, 2.1, "phi_Ds_M"); */
    /* draw_2d(mychain, "phi_CHI2:phi_M", ";#it{M(K^{+}K^{-})} [GeV];#it{#chi^{2}(#phi)}", 0.99, 1.05, 0, 10, "phi_M_CHI2"); */
    /* draw_2d(mychain, "Ds_CHI2:Ds_M", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];#it{#chi^{2}(D_{s}^{+})}", 1.85, 2.1, 0, 25, "Ds_M_CHI2"); */

    delete mychain;

    return 0;
}

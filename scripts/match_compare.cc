#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath, TCut cutstr="1"){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", vartitle, 100, varmin, varmax);

    mychain->Project("hist", myvar, cutstr);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw("ep");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
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
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
    delete h1;
    delete h2;
}

int match_compare() {
    setTDRStyle();

    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/match_compare")) gSystem->MakeDirectory("./figures/match_compare");

    TCut mass_cut = "DsFit_Ds_M > 1.85 && DsFit_Ds_M < 2.1 && phiFit_phi_M > 0.99 && phiFit_phi_M < 1.05";
    TCut match_sel = "match_entry == 1";
    TCut nonmatch_sel = "non_match_entry == 1";
    match_sel += mass_cut;
    nonmatch_sel += mass_cut;
  
    /* compare(mychain, "dxy_Kp_Km", "dxy_Kp_Km", "match", "non-match", "#it{d_{xy}(K^{+},K^{-})} [cm]", 0, 0.1, match_sel, nonmatch_sel, "dxy_Kp_Km"); */ 
    /* compare(mychain, "dxy_Kp_pi", "dxy_Kp_pi", "match", "non-match", "#it{d_{xy}(K^{+},#pi^{+})} [cm]", 0, 0.15, match_sel, nonmatch_sel, "dxy_Kp_pi"); */ 
    /* compare(mychain, "dxy_Km_pi", "dxy_Km_pi", "match", "non-match", "#it{d_{xy}(K^{-},#pi^{+})} [cm]", 0, 0.15, match_sel, nonmatch_sel, "dxy_Km_pi"); */ 
    /* compare(mychain, "dxy_Kp_phi", "dxy_Kp_phi", "match", "non-match", "#it{d_{xy}(K^{+},#phi)} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dxy_Kp_phi"); */ 
    /* compare(mychain, "dxy_Km_phi", "dxy_Km_phi", "match", "non-match", "#it{d_{xy}(K^{-},#phi)} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dxy_Km_phi"); */ 
    /* compare(mychain, "dxy_pi_phi", "dxy_pi_phi", "match", "non-match", "#it{d_{xy}(#pi^{+},#phi)} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dxy_pi_phi"); */ 
    /* compare(mychain, "dxy_Kp_Ds", "dxy_Kp_Ds", "match", "non-match", "#it{d_{xy}(K^{+},D_{s}^{+})} [cm]", 0, 1, match_sel, nonmatch_sel, "dxy_Kp_Ds"); */ 
    /* compare(mychain, "dxy_Km_Ds", "dxy_Km_Ds", "match", "non-match", "#it{d_{xy}(K^{-},D_{s}^{+})} [cm]", 0, 1, match_sel, nonmatch_sel, "dxy_Km_Ds"); */ 
    /* compare(mychain, "dxy_pi_Ds", "dxy_pi_Ds", "match", "non-match", "#it{d_{xy}(#pi^{+},D_{s}^{+})} [cm]", 0, 1, match_sel, nonmatch_sel, "dxy_pi_Ds"); */ 
    /* compare(mychain, "dxy_phi_Ds", "dxy_phi_Ds", "match", "non-match", "#it{d_{xy}(#phi,D_{s}^{+})} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dxy_phi_Ds"); */ 

    /* compare(mychain, "dz_Kp_Km", "dz_Kp_Km", "match", "non-match", "#it{d_{z}(K^{+},K^{-})} [cm]", 0, 0.2, match_sel, nonmatch_sel, "dz_Kp_Km"); */ 
    /* compare(mychain, "dz_Kp_pi", "dz_Kp_pi", "match", "non-match", "#it{d_{z}(K^{+},#pi^{+})} [cm]", 0, 0.2, match_sel, nonmatch_sel, "dz_Kp_pi"); */ 
    /* compare(mychain, "dz_Km_pi", "dz_Km_pi", "match", "non-match", "#it{d_{z}(K^{-},#pi^{+})} [cm]", 0, 0.2, match_sel, nonmatch_sel, "dz_Km_pi"); */ 
    /* compare(mychain, "dz_Kp_phi", "dz_Kp_phi", "match", "non-match", "#it{d_{z}(K^{+},#phi)} [cm]", 0, 2, match_sel, nonmatch_sel, "dz_Kp_phi"); */ 
    /* compare(mychain, "dz_Km_phi", "dz_Km_phi", "match", "non-match", "#it{d_{z}(K^{-},#phi)} [cm]", 0, 2, match_sel, nonmatch_sel, "dz_Km_phi"); */ 
    /* compare(mychain, "dz_pi_phi", "dz_pi_phi", "match", "non-match", "#it{d_{z}(#pi^{+},#phi)} [cm]", 0, 2, match_sel, nonmatch_sel, "dz_pi_phi"); */ 
    /* compare(mychain, "dz_Kp_Ds", "dz_Kp_Ds", "match", "non-match", "#it{d_{z}(K^{+},D_{s}^{+})} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dz_Kp_Ds"); */ 
    /* compare(mychain, "dz_Km_Ds", "dz_Km_Ds", "match", "non-match", "#it{d_{z}(K^{-},D_{s}^{+})} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dz_Km_Ds"); */ 
    /* compare(mychain, "dz_pi_Ds", "dz_pi_Ds", "match", "non-match", "#it{d_{z}(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, match_sel, nonmatch_sel, "dz_pi_Ds"); */ 
    /* compare(mychain, "dz_phi_Ds", "dz_phi_Ds", "match", "non-match", "#it{d_{z}(#phi,D_{s}^{+})} [cm]", 0, 2, match_sel, nonmatch_sel, "dz_phi_Ds"); */ 

    /* compare(mychain, "dxy_Kp_Km", "dz_Kp_Km", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},K^{-})} [cm]", 0, 0.2, match_sel, match_sel, "match_d_Kp_Km"); */ 
    /* compare(mychain, "dxy_Kp_pi", "dz_Kp_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#pi^{+})} [cm]", 0, 0.2, match_sel, match_sel, "match_d_Kp_pi"); */ 
    /* compare(mychain, "dxy_Km_pi", "dz_Km_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#pi^{+})} [cm]", 0, 0.2, match_sel, match_sel, "match_d_Km_pi"); */ 
    /* compare(mychain, "dxy_Kp_phi", "dz_Kp_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#phi)} [cm]", 0, 2, match_sel, match_sel, "match_d_Kp_phi"); */ 
    /* compare(mychain, "dxy_Km_phi", "dz_Km_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#phi)} [cm]", 0, 2, match_sel, match_sel, "match_d_Km_phi"); */ 
    /* compare(mychain, "dxy_pi_phi", "dz_pi_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},#phi)} [cm]", 0, 2, match_sel, match_sel, "match_d_pi_phi"); */ 
    /* compare(mychain, "dxy_Kp_Ds", "dz_Kp_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},D_{s}^{+})} [cm]", 0, 1.5, match_sel, match_sel, "match_d_Kp_Ds"); */ 
    /* compare(mychain, "dxy_Km_Ds", "dz_Km_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},D_{s}^{+})} [cm]", 0, 1.5, match_sel, match_sel, "match_d_Km_Ds"); */ 
    /* compare(mychain, "dxy_pi_Ds", "dz_pi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, match_sel, match_sel, "match_d_pi_Ds"); */ 
    /* compare(mychain, "dxy_phi_Ds", "dz_phi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#phi,D_{s}^{+})} [cm]", 0, 2, match_sel, match_sel, "match_d_phi_Ds"); */ 

    /* compare(mychain, "dxy_Kp_Km", "dz_Kp_Km", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},K^{-})} [cm]", 0, 0.2, nonmatch_sel, nonmatch_sel, "nonmatch_d_Kp_Km"); */ 
    /* compare(mychain, "dxy_Kp_pi", "dz_Kp_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#pi^{+})} [cm]", 0, 0.2, nonmatch_sel, nonmatch_sel, "nonmatch_d_Kp_pi"); */ 
    /* compare(mychain, "dxy_Km_pi", "dz_Km_pi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#pi^{+})} [cm]", 0, 0.2, nonmatch_sel, nonmatch_sel, "nonmatch_d_Km_pi"); */ 
    /* compare(mychain, "dxy_Kp_phi", "dz_Kp_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},#phi)} [cm]", 0, 2, nonmatch_sel, nonmatch_sel, "nonmatch_d_Kp_phi"); */ 
    /* compare(mychain, "dxy_Km_phi", "dz_Km_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},#phi)} [cm]", 0, 2, nonmatch_sel, nonmatch_sel, "nonmatch_d_Km_phi"); */ 
    /* compare(mychain, "dxy_pi_phi", "dz_pi_phi", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},#phi)} [cm]", 0, 2, nonmatch_sel, nonmatch_sel, "nonmatch_d_pi_phi"); */ 
    /* compare(mychain, "dxy_Kp_Ds", "dz_Kp_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{+},D_{s}^{+})} [cm]", 0, 1.5, nonmatch_sel, nonmatch_sel, "nonmatch_d_Kp_Ds"); */ 
    /* compare(mychain, "dxy_Km_Ds", "dz_Km_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(K^{-},D_{s}^{+})} [cm]", 0, 1.5, nonmatch_sel, nonmatch_sel, "nonmatch_d_Km_Ds"); */ 
    /* compare(mychain, "dxy_pi_Ds", "dz_pi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#pi^{+},D_{s}^{+})} [cm]", 0, 1.5, nonmatch_sel, nonmatch_sel, "nonmatch_d_pi_Ds"); */ 
    /* compare(mychain, "dxy_phi_Ds", "dz_phi_Ds", "#it{d_{xy}}", "#it{d_{z}}", "#it{d(#phi,D_{s}^{+})} [cm]", 0, 2, nonmatch_sel, nonmatch_sel, "nonmatch_d_phi_Ds"); */ 


    /* compare(mychain, "dR_Kp_Km", "dR_Kp_Km", "match", "non-match", "#it{#DeltaR(K^{+},K^{-})}", 0, 0.15, match_sel, nonmatch_sel, "dR_Kp_Km"); */ 
    /* compare(mychain, "dR_Kp_pi", "dR_Kp_pi", "match", "non-match", "#it{#DeltaR(K^{+},#pi^{+})}", 0, 0.6, match_sel, nonmatch_sel, "dR_Kp_pi"); */ 
    /* compare(mychain, "dR_Km_pi", "dR_Km_pi", "match", "non-match", "#it{#DeltaR(K^{-},#pi^{+})}", 0, 0.6, match_sel, nonmatch_sel, "dR_Km_pi"); */ 
    /* compare(mychain, "phiFit_dR_Kp_phi", "phiFit_dR_Kp_phi", "match", "non-match", "#it{#DeltaR(K^{+},#phi)}", 0, 0.1, match_sel, nonmatch_sel, "phiFit_dR_Kp_phi"); */ 
    /* compare(mychain, "phiFit_dR_Km_phi", "phiFit_dR_Km_phi", "match", "non-match", "#it{#DeltaR(K^{-},#phi)}", 0, 0.1, match_sel, nonmatch_sel, "phiFit_dR_Km_phi"); */ 
    /* compare(mychain, "phiFit_dR_pi_phi", "phiFit_dR_pi_phi", "match", "non-match", "#it{#DeltaR(#pi^{+},#phi)}", 0, 0.6, match_sel, nonmatch_sel, "phiFit_dR_pi_phi"); */ 
    /* compare(mychain, "DsFit_dR_Kp_Ds", "DsFit_dR_Kp_Ds", "match", "non-match", "#it{#DeltaR(K^{+},D_{s}^{+})}", 0, 0.4, match_sel, nonmatch_sel, "DsFit_dR_Kp_Ds"); */ 
    /* compare(mychain, "DsFit_dR_Km_Ds", "DsFit_dR_Km_Ds", "match", "non-match", "#it{#DeltaR(K^{-},D_{s}^{+})}", 0, 0.4, match_sel, nonmatch_sel, "DsFit_dR_Km_Ds"); */ 
    /* compare(mychain, "DsFit_dR_pi_Ds", "DsFit_dR_pi_Ds", "match", "non-match", "#it{#DeltaR(#pi^{+},D_{s}^{+})}", 0, 0.6, match_sel, nonmatch_sel, "DsFit_dR_pi_Ds"); */ 
    /* compare(mychain, "DsFit_dR_phi_Ds", "DsFit_dR_phi_Ds", "match", "non-match", "#it{#DeltaR(#phi,D_{s}^{+})}", 0, 0.4, match_sel, nonmatch_sel, "DsFit_dR_phi_Ds"); */ 

    /* compare(mychain, "Kp_PT", "Kp_PT", "match", "non-match", "#it{p_{T}(K^{+})} [GeV]", 0, 30, match_sel, nonmatch_sel, "Kp_PT"); */ 
    /* compare(mychain, "Km_PT", "Km_PT", "match", "non-match", "#it{p_{T}(K^{-})} [GeV]", 0, 30, match_sel, nonmatch_sel, "Km_PT"); */ 
    /* compare(mychain, "pi_PT", "pi_PT", "match", "non-match", "#it{p_{T}(#pi^{+})} [GeV]", 0, 30, match_sel, nonmatch_sel, "pi_PT"); */ 
    /* compare(mychain, "phiFit_phi_PT", "phiFit_phi_PT", "match", "non-match", "#it{p_{T}(#phi)} [GeV]", 0, 60, match_sel, nonmatch_sel, "phiFit_phi_PT"); */ 
    /* compare(mychain, "DsFit_Ds_PT", "DsFit_Ds_PT", "match", "non-match", "#it{p_{T}(D_{s}^{+})} [GeV]", 0, 80, match_sel, nonmatch_sel, "DsFit_Ds_PT"); */ 

    /* compare(mychain, "Kp_P", "Kp_P", "match", "non-match", "#it{p(K^{+})} [GeV]", 0, 50, match_sel, nonmatch_sel, "Kp_P"); */ 
    /* compare(mychain, "Km_P", "Km_P", "match", "non-match", "#it{p(K^{-})} [GeV]", 0, 50, match_sel, nonmatch_sel, "Km_P"); */ 
    /* compare(mychain, "pi_P", "pi_P", "match", "non-match", "#it{p(#pi^{+})} [GeV]", 0, 50, match_sel, nonmatch_sel, "pi_P"); */ 
    /* compare(mychain, "phiFit_phi_P", "phiFit_phi_P", "match", "non-match", "#it{p(#phi)} [GeV]", 0, 100, match_sel, nonmatch_sel, "phiFit_phi_P"); */ 
    /* compare(mychain, "DsFit_Ds_P", "DsFit_Ds_P", "match", "non-match", "#it{p(D_{s}^{+})} [GeV]", 0, 150, match_sel, nonmatch_sel, "DsFit_Ds_P"); */ 

    /* compare(mychain, "Kp_PHI", "Kp_PHI", "match", "non-match", "#it{#phi(K^{+})}", -3.15, 3.15, match_sel, nonmatch_sel, "Kp_PHI"); */ 
    /* compare(mychain, "Km_PHI", "Km_PHI", "match", "non-match", "#it{#phi(K^{-})}", -3.15, 3.15, match_sel, nonmatch_sel, "Km_PHI"); */ 
    /* compare(mychain, "pi_PHI", "pi_PHI", "match", "non-match", "#it{#phi(#pi^{+})}", -3.15, 3.15, match_sel, nonmatch_sel, "pi_PHI"); */ 
    /* compare(mychain, "phiFit_phi_PHI", "phiFit_phi_PHI", "match", "non-match", "#it{#phi(#phi)}", -3.15, 3.15, match_sel, nonmatch_sel, "phiFit_phi_PHI"); */ 
    /* compare(mychain, "DsFit_Ds_PHI", "DsFit_Ds_PHI", "match", "non-match", "#it{#phi(D_{s}^{+})}", -3.15, 3.15, match_sel, nonmatch_sel, "DsFit_Ds_PHI"); */ 

    /* compare(mychain, "Kp_ETA", "Kp_ETA", "match", "non-match", "#it{#eta(K^{+})}", -5, 5, match_sel, nonmatch_sel, "Kp_ETA"); */ 
    /* compare(mychain, "Km_ETA", "Km_ETA", "match", "non-match", "#it{#eta(K^{-})}", -5, 5, match_sel, nonmatch_sel, "Km_ETA"); */ 
    /* compare(mychain, "pi_ETA", "pi_ETA", "match", "non-match", "#it{#eta(#pi^{+})}", -5, 5, match_sel, nonmatch_sel, "pi_ETA"); */ 
    /* compare(mychain, "phiFit_phi_ETA", "phiFit_phi_ETA", "match", "non-match", "#it{#eta(#phi)}", -5, 5, match_sel, nonmatch_sel, "phiFit_phi_ETA"); */ 
    /* compare(mychain, "DsFit_Ds_ETA", "DsFit_Ds_ETA", "match", "non-match", "#it{#eta(D_{s}^{+})}", -5, 5, match_sel, nonmatch_sel, "DsFit_Ds_ETA"); */ 
   
    /* draw_1d(mychain, "Kp_ORIVX_X", ";#it{x}_{orig}#it{(K^{+})} [cm];# candidates", -0.1, 0.05, "match_Kp_ORIVX_X", match_sel); */
    /* draw_1d(mychain, "Kp_ORIVX_X", ";#it{x}_{orig}#it{(K^{+})} [cm];# candidates", -0.1, 0.05, "nonmatch_Kp_ORIVX_X", nonmatch_sel); */
    /* draw_1d(mychain, "Kp_ORIVX_Y", ";#it{y}_{orig}#it{(K^{+})} [cm];# candidates", 0, 0.15, "match_Kp_ORIVX_Y", match_sel); */
    /* draw_1d(mychain, "Kp_ORIVX_Y", ";#it{y}_{orig}#it{(K^{+})} [cm];# candidates", 0, 0.15, "nonmatch_Kp_ORIVX_Y", nonmatch_sel); */
    /* draw_1d(mychain, "Kp_ORIVX_Z", ";#it{z}_{orig}#it{(K^{+})} [cm];# candidates", -15, 15, "match_Kp_ORIVX_Z", match_sel); */
    /* draw_1d(mychain, "Kp_ORIVX_Z", ";#it{z}_{orig}#it{(K^{+})} [cm];# candidates", -15, 15, "nonmatch_Kp_ORIVX_Z", nonmatch_sel); */

    /* draw_1d(mychain, "phiFit_ENDVX_X", ";#it{x}_{decay}#it{(#phi)} [cm];# candidates", -2, 2, "match_phiFit_ENDVX_X", match_sel); */
    /* draw_1d(mychain, "phiFit_ENDVX_X", ";#it{x}_{decay}#it{(#phi)} [cm];# candidates", -2, 2, "nonmatch_phiFit_ENDVX_X", nonmatch_sel); */
    /* draw_1d(mychain, "phiFit_ENDVX_Y", ";#it{y}_{decay}#it{(#phi)} [cm];# candidates", -2, 2, "match_phiFit_ENDVX_Y", match_sel); */
    /* draw_1d(mychain, "phiFit_ENDVX_Y", ";#it{y}_{decay}#it{(#phi)} [cm];# candidates", -2, 2, "nonmatch_phiFit_ENDVX_Y", nonmatch_sel); */
    /* draw_1d(mychain, "phiFit_ENDVX_Z", ";#it{z}_{decay}#it{(#phi)} [cm];# candidates", -15, 15, "match_phiFit_ENDVX_Z", match_sel); */
    /* draw_1d(mychain, "phiFit_ENDVX_Z", ";#it{z}_{decay}#it{(#phi)} [cm];# candidates", -15, 15, "nonmatch_phiFit_ENDVX_Z", nonmatch_sel); */

    compare(mychain, "Kp_ORIVX_X", "Kp_ORIVX_X", "match", "non-match", "#it{x}_{orig}#it{(K^{+})} [cm]", -0.1, 0.05, match_sel, nonmatch_sel, "Kp_ORIVX_X"); 
    compare(mychain, "Kp_ORIVX_Y", "Kp_ORIVX_Y", "match", "non-match", "#it{y}_{orig}#it{(K^{+})} [cm]", 0, 0.15, match_sel, nonmatch_sel, "Kp_ORIVX_Y"); 
    compare(mychain, "Kp_ORIVX_Z", "Kp_ORIVX_Z", "match", "non-match", "#it{z}_{orig}#it{(K^{+})} [cm]", -15, 15, match_sel, nonmatch_sel, "Kp_ORIVX_Z"); 

    compare(mychain, "phiFit_ENDVX_X", "phiFit_ENDVX_X", "match", "non-match", "#it{x}_{decay}(#it{#phi}) [cm]", -2, 2, match_sel, nonmatch_sel, "phiFit_ENDVX_X"); 
    compare(mychain, "phiFit_ENDVX_Y", "phiFit_ENDVX_Y", "match", "non-match", "#it{y}_{decay}(#it{#phi}) [cm]", -2, 2, match_sel, nonmatch_sel, "phiFit_ENDVX_Y"); 
    compare(mychain, "phiFit_ENDVX_Z", "phiFit_ENDVX_Z", "match", "non-match", "#it{z}_{decay}(#it{#phi}) [cm]", -15, 15, match_sel, nonmatch_sel, "phiFit_ENDVX_Z"); 

    compare(mychain, "DsFit_ENDVX_X", "DsFit_ENDVX_X", "match", "non-match", "#it{x}_{decay}(#it{D_{s}^{+}}) [cm]", -2, 2, match_sel, nonmatch_sel, "DsFit_ENDVX_X"); 
    compare(mychain, "DsFit_ENDVX_Y", "DsFit_ENDVX_Y", "match", "non-match", "#it{y}_{decay}(#it{D_{s}^{+}}) [cm]", -2, 2, match_sel, nonmatch_sel, "DsFit_ENDVX_Y"); 
    compare(mychain, "DsFit_ENDVX_Z", "DsFit_ENDVX_Z", "match", "non-match", "#it{z}_{decay}(#it{D_{s}^{+}}) [cm]", -15, 15, match_sel, nonmatch_sel, "DsFit_ENDVX_Z"); 

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

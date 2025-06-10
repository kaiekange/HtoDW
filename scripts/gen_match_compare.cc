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
    canvas->SaveAs("./figures/gen_match_compare/"+figpath+".png");
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
    canvas->SaveAs("./figures/gen_match_compare/"+figpath+".png");
    delete h1;
    delete h2;
}

int gen_match_compare() {
    setTDRStyle();

    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/gen_match_compare")) gSystem->MakeDirectory("./figures/gen_match_compare");

    /* compare(mychain, "Gen_Kp_ETA", "match_Kp_ETA", "Gen", "match", "#it{#eta}(#it{K^{+}})", -6, 6, "1", "1", "Kp_ETA"); */
    /* compare(mychain, "Gen_Km_ETA", "match_Km_ETA", "Gen", "match", "#it{#eta}(#it{K^{-}})", -6, 6, "1", "1", "Km_ETA"); */
    /* compare(mychain, "Gen_pi_ETA", "match_pi_ETA", "Gen", "match", "#it{#eta}(#it{#pi^{+}})", -6, 6, "1", "1", "pi_ETA"); */
    /* compare(mychain, "Gen_phi_ETA", "match_phiFit_phi_ETA", "Gen", "match", "#it{#eta}(#it{#phi})", -6, 6, "1", "1", "phi_ETA"); */
    /* compare(mychain, "Gen_Ds_ETA", "match_DsFit_Ds_ETA", "Gen", "match", "#it{#eta}(#it{D_{s}^{+}})", -6, 6, "1", "1", "Ds_ETA"); */

    /* compare(mychain, "Gen_Kp_PHI", "match_Kp_PHI", "Gen", "match", "#it{#phi}(#it{K^{+}})", -3.14, 3.14, "1", "1", "Kp_PHI"); */
    /* compare(mychain, "Gen_Km_PHI", "match_Km_PHI", "Gen", "match", "#it{#phi}(#it{K^{-}})", -3.14, 3.14, "1", "1", "Km_PHI"); */
    /* compare(mychain, "Gen_pi_PHI", "match_pi_PHI", "Gen", "match", "#it{#phi}(#it{#pi^{+}})", -3.14, 3.14, "1", "1", "pi_PHI"); */
    /* compare(mychain, "Gen_phi_PHI", "match_phiFit_phi_PHI", "Gen", "match", "#it{#phi}(#it{#phi})", -3.14, 3.14, "1", "1", "phi_PHI"); */
    /* compare(mychain, "Gen_Ds_PHI", "match_DsFit_Ds_PHI", "Gen", "match", "#it{#phi}(#it{D_{s}^{+}})", -3.14, 3.14, "1", "1", "Ds_PHI"); */

    /* compare(mychain, "Gen_Kp_ORIVX_X", "match_Kp_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{K^{+}}) [cm]",  -0.5, 0.5, "1", "1", "Kp_ORIVX_X"); */
    /* compare(mychain, "Gen_Km_ORIVX_X", "match_Km_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{K^{-}}) [cm]",  -0.5, 0.5, "1", "1", "Km_ORIVX_X"); */
    /* compare(mychain, "Gen_pi_ORIVX_X", "match_pi_ORIVX_X", "Gen", "match", "#it{x}_{orig}(#it{#pi^{+}}) [cm]",  -0.5, 0.5, "1", "1", "pi_ORIVX_X"); */

    /* compare(mychain, "Gen_Kp_ORIVX_Y", "match_Kp_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{K^{+}}) [cm]",  -0.5, 0.5, "1", "1", "Kp_ORIVX_Y"); */
    /* compare(mychain, "Gen_Km_ORIVX_Y", "match_Km_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{K^{-}}) [cm]",  -0.5, 0.5, "1", "1", "Km_ORIVX_Y"); */
    /* compare(mychain, "Gen_pi_ORIVX_Y", "match_pi_ORIVX_Y", "Gen", "match", "#it{y}_{orig}(#it{#pi^{+}}) [cm]",  -0.5, 0.5, "1", "1", "pi_ORIVX_Y"); */

    /* compare(mychain, "Gen_Kp_ORIVX_Z", "match_Kp_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{K^{+}})",  -20, 20, "1", "1", "Kp_ORIVX_Z"); */
    /* compare(mychain, "Gen_Km_ORIVX_Z", "match_Km_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{K^{-}})",  -20, 20, "1", "1", "Km_ORIVX_Z"); */
    /* compare(mychain, "Gen_pi_ORIVX_Z", "match_pi_ORIVX_Z", "Gen", "match", "#it{z}_{orig}(#it{#pi^{+}})",  -20, 20, "1", "1", "pi_ORIVX_Z"); */

    /* compare(mychain, "Gen_Kp_P", "match_Kp_P", "Gen", "match", "#it{p}(#it{K^{+}}) [GeV]", 0, 300, "1", "1", "Kp_P"); */
    /* compare(mychain, "Gen_Km_P", "match_Km_P", "Gen", "match", "#it{p}(#it{K^{-}}) [GeV]", 0, 300, "1", "1", "Km_P"); */
    /* compare(mychain, "Gen_pi_P", "match_pi_P", "Gen", "match", "#it{p}(#it{#pi^{+}}) [GeV]", 0, 300, "1", "1", "pi_P"); */
    /* compare(mychain, "Gen_phi_P", "match_phiFit_phi_P", "Gen", "match", "#it{p}(#it{#phi}) [GeV]", 0, 600, "1", "1", "phi_P"); */
    /* compare(mychain, "Gen_Ds_P", "match_DsFit_Ds_P", "Gen", "match", "#it{p}(#it{D_{s}^{+}}) [GeV]", 0, 800, "1", "1", "Ds_P"); */

    /* compare(mychain, "Gen_Kp_PT", "match_Kp_PT", "Gen", "match", "#it{p_{T}}(#it{K^{+}}) [GeV]", 0, 50, "1", "1", "Kp_PT"); */
    /* compare(mychain, "Gen_Km_PT", "match_Km_PT", "Gen", "match", "#it{p_{T}}(#it{K^{-}}) [GeV]", 0, 50, "1", "1", "Km_PT"); */
    /* compare(mychain, "Gen_pi_PT", "match_pi_PT", "Gen", "match", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", 0, 80, "1", "1", "pi_PT"); */
    /* compare(mychain, "Gen_phi_PT", "match_phiFit_phi_PT", "Gen", "match", "#it{p_{T}}(#it{#phi}) [GeV]", 0, 100, "1", "1", "phi_PT"); */
    /* compare(mychain, "Gen_Ds_PT", "match_DsFit_Ds_PT", "Gen", "match", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", 0, 150, "1", "1", "Ds_PT"); */

    /* draw_1d(mychain, "(Gen_Kp_ETA - match_Kp_ETA)/Gen_Kp_ETA", ";#it{#epsilon#eta}(#it{K^{+}});# candidates", -0.5, 0.5, "diff_Kp_ETA", "1"); */
    /* draw_1d(mychain, "(Gen_Km_ETA - match_Km_ETA)/Gen_Km_ETA", ";#it{#epsilon#eta}(#it{K^{-}});# candidates", -0.5, 0.5, "diff_Km_ETA", "1"); */
    /* draw_1d(mychain, "(Gen_pi_ETA - match_pi_ETA)/Gen_pi_ETA", ";#it{#epsilon#eta}(#it{#pi^{+}});# candidates", -0.5, 0.5, "diff_pi_ETA", "1"); */
    /* draw_1d(mychain, "(Gen_phi_ETA - match_phiFit_phi_ETA)/Gen_phi_ETA", ";#it{#epsilon#eta}(#it{#phi});# candidates", -0.5, 0.5, "diff_phi_ETA", "1"); */
    /* draw_1d(mychain, "(Gen_Ds_ETA - match_DsFit_Ds_ETA)/Gen_Ds_ETA", ";#it{#epsilon#eta}(#it{D_{s}^{+}});# candidates", -0.5, 0.5, "diff_Ds_ETA", "1"); */
  
    /* draw_1d(mychain, "(Gen_Kp_PHI - match_Kp_PHI)/Gen_Kp_PHI", ";#it{#epsilon#phi}(#it{K^{+}});# candidates", -0.5, 0.5, "diff_Kp_PHI", "1"); */
    /* draw_1d(mychain, "(Gen_Km_PHI - match_Km_PHI)/Gen_Km_PHI", ";#it{#epsilon#phi}(#it{K^{-}});# candidates", -0.5, 0.5, "diff_Km_PHI", "1"); */
    /* draw_1d(mychain, "(Gen_pi_PHI - match_pi_PHI)/Gen_pi_PHI", ";#it{#epsilon#phi}(#it{#pi^{+}});# candidates", -0.5, 0.5, "diff_pi_PHI", "1"); */
    /* draw_1d(mychain, "(Gen_phi_PHI - match_phiFit_phi_PHI)/Gen_phi_PHI", ";#it{#epsilon#phi}(#it{#phi});# candidates", -0.5, 0.5, "diff_phi_PHI", "1"); */
    /* draw_1d(mychain, "(Gen_Ds_PHI - match_DsFit_Ds_PHI)/Gen_Ds_PHI", ";#it{#epsilon#phi}(#it{D_{s}^{+}});# candidates", -0.5, 0.5, "diff_Ds_PHI", "1"); */
    
    /* draw_1d(mychain, "(Gen_Kp_P - match_Kp_P)/Gen_Kp_P", ";#it{#epsilonp}(#it{K^{+}});# candidates", -1, 1, "diff_Kp_P", "1"); */
    /* draw_1d(mychain, "(Gen_Km_P - match_Km_P)/Gen_Km_P", ";#it{#epsilonp}(#it{K^{-}});# candidates", -1, 1, "diff_Km_P", "1"); */
    /* draw_1d(mychain, "(Gen_pi_P - match_pi_P)/Gen_pi_P", ";#it{#epsilonp}(#it{#pi^{+}});# candidates", -1, 1, "diff_pi_P", "1"); */
    /* draw_1d(mychain, "(Gen_phi_P - match_phiFit_phi_P)/Gen_phi_P", ";#it{#epsilonp}(#it{#phi});# candidates", -1, 1, "diff_phi_P", "1"); */
    /* draw_1d(mychain, "(Gen_Ds_P - match_DsFit_Ds_P)/Gen_Ds_P", ";#it{#epsilonp}(#it{D_{s}^{+}});# candidates", -1, 1, "diff_Ds_P", "1"); */
    
    /* draw_1d(mychain, "(Gen_Kp_PT - match_Kp_PT)/Gen_Kp_PT", ";#it{#epsilonp_{T}}(#it{K^{+}});# candidates", -1, 1, "diff_Kp_PT", "1"); */
    /* draw_1d(mychain, "(Gen_Km_PT - match_Km_PT)/Gen_Km_PT", ";#it{#epsilonp_{T}}(#it{K^{-}});# candidates", -1, 1, "diff_Km_PT", "1"); */
    /* draw_1d(mychain, "(Gen_pi_PT - match_pi_PT)/Gen_pi_PT", ";#it{#epsilonp_{T}}(#it{#pi^{+}});# candidates", -1, 1, "diff_pi_PT", "1"); */
    /* draw_1d(mychain, "(Gen_phi_PT - match_phiFit_phi_PT)/Gen_phi_PT", ";#it{#epsilonp_{T}}(#it{#phi});# candidates", -1, 1, "diff_phi_PT", "1"); */
    /* draw_1d(mychain, "(Gen_Ds_PT - match_DsFit_Ds_PT)/Gen_Ds_PT", ";#it{#epsilonp_{T}}(#it{D_{s}^{+}});# candidates", -1, 1, "diff_Ds_PT", "1"); */
    
    draw_1d(mychain, "(Gen_Kp_ORIVX_X - match_Kp_ORIVX_X)/Gen_Kp_ORIVX_X", ";#it{#epsilonx}_{orig}(#it{K^{+}});# candidates", -0.5, 2, "diff_Kp_ORIVX_X", "1");
    draw_1d(mychain, "(Gen_Km_ORIVX_X - match_Km_ORIVX_X)/Gen_Km_ORIVX_X", ";#it{#epsilonx}_{orig}(#it{K^{-}});# candidates", -0.5, 2, "diff_Km_ORIVX_X", "1");
    draw_1d(mychain, "(Gen_pi_ORIVX_X - match_pi_ORIVX_X)/Gen_pi_ORIVX_X", ";#it{#epsilonx}_{orig}(#it{#pi^{+}});# candidates", -0.5, 2, "diff_pi_ORIVX_X", "1");
    
    draw_1d(mychain, "(Gen_Kp_ORIVX_Y - match_Kp_ORIVX_Y)/Gen_Kp_ORIVX_Y", ";#it{#epsilony}_{orig}(#it{K^{+}});# candidates", -0.5, 2, "diff_Kp_ORIVX_Y", "1");
    draw_1d(mychain, "(Gen_Km_ORIVX_Y - match_Km_ORIVX_Y)/Gen_Km_ORIVX_Y", ";#it{#epsilony}_{orig}(#it{K^{-}});# candidates", -0.5, 2, "diff_Km_ORIVX_Y", "1");
    draw_1d(mychain, "(Gen_pi_ORIVX_Y - match_pi_ORIVX_Y)/Gen_pi_ORIVX_Y", ";#it{#epsilony}_{orig}(#it{#pi^{+}});# candidates", -0.5, 2, "diff_pi_ORIVX_Y", "1");
    
    draw_1d(mychain, "(Gen_Kp_ORIVX_Z - match_Kp_ORIVX_Z)/Gen_Kp_ORIVX_Z", ";#it{#epsilonz}_{orig}(#it{K^{+}});# candidates", -1, 1, "diff_Kp_ORIVX_Z", "1");
    draw_1d(mychain, "(Gen_Km_ORIVX_Z - match_Km_ORIVX_Z)/Gen_Km_ORIVX_Z", ";#it{#epsilonz}_{orig}(#it{K^{-}});# candidates", -1, 1, "diff_Km_ORIVX_Z", "1");
    draw_1d(mychain, "(Gen_pi_ORIVX_Z - match_pi_ORIVX_Z)/Gen_pi_ORIVX_Z", ";#it{#epsilonz}_{orig}(#it{#pi^{+}});# candidates", -1, 1, "diff_pi_ORIVX_Z", "1");
    
    delete mychain;

    return 0;
}

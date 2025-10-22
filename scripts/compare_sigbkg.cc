#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void compare(TChain *sigchain, TChain *bkgchain, TString sigvar, TString bkgvar, TString sigleg, TString bkgleg, TString vartitle, float varmin, float varmax, TCut sigcut, TCut bkgcut, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", 100, varmin, varmax);
    sigchain->Project("h1", sigvar, sigcut);
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlue);
    h1->SetMarkerColor(kBlue);
    h1->SetMarkerSize(0.7);


    TH1F *h2 = new TH1F("h2", "", 100, varmin, varmax);
    bkgchain->Project("h2", bkgvar, bkgcut);
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

    TLegend * mylegend = new TLegend(0.65, 0.75, 0.85, 0.9);
    mylegend->AddEntry(h1, sigleg, "l");
    mylegend->AddEntry(h2, bkgleg, "l");
    mylegend->SetBorderSize(0);
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


int compare_sigbkg() {
    setTDRStyle();

    TChain *sigchain = new TChain("RecoAnalyzer/Events");
    TChain *bkgchain = new TChain("RecoBestAnalyzer/Events");
    sigchain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/RecoAnalyzer/RecoStudy.root");
    bkgchain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Reconstruction/20250417/2017UL/WJetsBKG/tuples/WJetsBKG.root");

    if(gSystem->AccessPathName("./SBcompare")) gSystem->MakeDirectory("./SBcompare");

    TCut sel_withBS = "PV_withBS_IsValid && !PV_withBS_IsFake";
    TCut sel_noBS = "PV_noBS_IsValid && !PV_noBS_IsFake";

    compare(sigchain, bkgchain, "match_DsFit_Ds_M", "best_DsFit_Ds_M", "signal", "W+jets", "#it{M}(#it{K^{+}K^{-}#pi^{+}}) [GeV]", 1.85, 2.1, "1", "1", "./SBcompare/Ds_M.png");
    compare(sigchain, bkgchain, "match_phiFit_phi_M", "best_phiFit_phi_M", "signal", "W+jets", "#it{M}(#it{K^{+}K^{-}}) [GeV]", 0.99, 1.05, "1", "1", "./SBcompare/phi_M.png");
    
    compare(sigchain, bkgchain, "match_Kp_ETA", "best_Kp_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", -3, 3, "1", "1", "./SBcompare/Kp_ETA.png");
    compare(sigchain, bkgchain, "match_Km_ETA", "best_Km_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", -3, 3, "1", "1", "./SBcompare/Km_ETA.png");
    compare(sigchain, bkgchain, "match_pi_ETA", "best_pi_ETA", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", -3, 3, "1", "1", "./SBcompare/pi_ETA.png");
    compare(sigchain, bkgchain, "match_phi_ETA", "best_phi_ETA", "signal", "W+jets", "#it{#eta}(#it{#phi})", -3, 3, "1", "1", "./SBcompare/phi_ETA.png");
    compare(sigchain, bkgchain, "match_Ds_ETA", "best_Ds_ETA", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", -3, 3, "1", "1", "./SBcompare/Ds_ETA.png");

    compare(sigchain, bkgchain, "match_phiFit_Kp_ETA", "best_phiFit_Kp_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", -3, 3, "1", "1", "./SBcompare/Kp_ETA_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Km_ETA", "best_phiFit_Km_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", -3, 3, "1", "1", "./SBcompare/Km_ETA_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_pi_ETA", "best_phiFit_pi_ETA", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", -3, 3, "1", "1", "./SBcompare/pi_ETA_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_phi_ETA", "best_phiFit_phi_ETA", "signal", "W+jets", "#it{#eta}(#it{#phi})", -3, 3, "1", "1", "./SBcompare/phi_ETA_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Ds_ETA", "best_phiFit_Ds_ETA", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", -3, 3, "1", "1", "./SBcompare/Ds_ETA_phiFit.png");

    compare(sigchain, bkgchain, "match_DsFit_Kp_ETA", "best_DsFit_Kp_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", -3, 3, "1", "1", "./SBcompare/Kp_ETA_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Km_ETA", "best_DsFit_Km_ETA", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", -3, 3, "1", "1", "./SBcompare/Km_ETA_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_pi_ETA", "best_DsFit_pi_ETA", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", -3, 3, "1", "1", "./SBcompare/pi_ETA_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_phi_ETA", "best_DsFit_phi_ETA", "signal", "W+jets", "#it{#eta}(#it{#phi})", -3, 3, "1", "1", "./SBcompare/phi_ETA_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Ds_ETA", "best_DsFit_Ds_ETA", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", -3, 3, "1", "1", "./SBcompare/Ds_ETA_DsFit.png");

    compare(sigchain, bkgchain, "match_Kp_PHI", "best_Kp_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Kp_PHI.png");
    compare(sigchain, bkgchain, "match_Km_PHI", "best_Km_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", -3.15, 3.15, "1", "1", "./SBcompare/Km_PHI.png");
    compare(sigchain, bkgchain, "match_pi_PHI", "best_pi_PHI", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/pi_PHI.png");
    compare(sigchain, bkgchain, "match_phi_PHI", "best_phi_PHI", "signal", "W+jets", "#it{#phi}(#it{#phi})", -3.15, 3.15, "1", "1", "./SBcompare/phi_PHI.png");
    compare(sigchain, bkgchain, "match_Ds_PHI", "best_Ds_PHI", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Ds_PHI.png");

    compare(sigchain, bkgchain, "match_phiFit_Kp_PHI", "best_phiFit_Kp_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Kp_PHI_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Km_PHI", "best_phiFit_Km_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", -3.15, 3.15, "1", "1", "./SBcompare/Km_PHI_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_pi_PHI", "best_phiFit_pi_PHI", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/pi_PHI_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_phi_PHI", "best_phiFit_phi_PHI", "signal", "W+jets", "#it{#phi}(#it{#phi})", -3.15, 3.15, "1", "1", "./SBcompare/phi_PHI_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Ds_PHI", "best_phiFit_Ds_PHI", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Ds_PHI_phiFit.png");

    compare(sigchain, bkgchain, "match_DsFit_Kp_PHI", "best_DsFit_Kp_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Kp_PHI_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Km_PHI", "best_DsFit_Km_PHI", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", -3.15, 3.15, "1", "1", "./SBcompare/Km_PHI_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_pi_PHI", "best_DsFit_pi_PHI", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/pi_PHI_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_phi_PHI", "best_DsFit_phi_PHI", "signal", "W+jets", "#it{#phi}(#it{#phi})", -3.15, 3.15, "1", "1", "./SBcompare/phi_PHI_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Ds_PHI", "best_DsFit_Ds_PHI", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", -3.15, 3.15, "1", "1", "./SBcompare/Ds_PHI_DsFit.png");

    compare(sigchain, bkgchain, "match_Kp_PT", "best_Kp_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Kp_PT.png");
    compare(sigchain, bkgchain, "match_Km_PT", "best_Km_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Km_PT.png");
    compare(sigchain, bkgchain, "match_pi_PT", "best_pi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/pi_PT.png");
    compare(sigchain, bkgchain, "match_phi_PT", "best_phi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", 0, 80, "1", "1", "./SBcompare/phi_PT.png");
    compare(sigchain, bkgchain, "match_Ds_PT", "best_Ds_PT", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Ds_PT.png");

    compare(sigchain, bkgchain, "match_phiFit_Kp_PT", "best_phiFit_Kp_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Kp_PT_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Km_PT", "best_phiFit_Km_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Km_PT_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_pi_PT", "best_phiFit_pi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/pi_PT_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_phi_PT", "best_phiFit_phi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", 0, 80, "1", "1", "./SBcompare/phi_PT_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Ds_PT", "best_phiFit_Ds_PT", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Ds_PT_phiFit.png");

    compare(sigchain, bkgchain, "match_DsFit_Kp_PT", "best_DsFit_Kp_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Kp_PT_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Km_PT", "best_DsFit_Km_PT", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", 0, 50, "1", "1", "./SBcompare/Km_PT_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_pi_PT", "best_DsFit_pi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", 0, 50, "1", "1", "./SBcompare/pi_PT_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_phi_PT", "best_DsFit_phi_PT", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", 0, 80, "1", "1", "./SBcompare/phi_PT_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Ds_PT", "best_DsFit_Ds_PT", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Ds_PT_DsFit.png");

    compare(sigchain, bkgchain, "match_Kp_P", "best_Kp_P", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Kp_P.png");
    compare(sigchain, bkgchain, "match_Km_P", "best_Km_P", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Km_P.png");
    compare(sigchain, bkgchain, "match_pi_P", "best_pi_P", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/pi_P.png");
    compare(sigchain, bkgchain, "match_phi_P", "best_phi_P", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", 0, 200, "1", "1", "./SBcompare/phi_P.png");
    compare(sigchain, bkgchain, "match_Ds_P", "best_Ds_P", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", 0, 250, "1", "1", "./SBcompare/Ds_P.png");

    compare(sigchain, bkgchain, "match_phiFit_Kp_P", "best_phiFit_Kp_P", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Kp_P_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Km_P", "best_phiFit_Km_P", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Km_P_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_pi_P", "best_phiFit_pi_P", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/pi_P_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_phi_P", "best_phiFit_phi_P", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", 0, 200, "1", "1", "./SBcompare/phi_P_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_Ds_P", "best_phiFit_Ds_P", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", 0, 250, "1", "1", "./SBcompare/Ds_P_phiFit.png");

    compare(sigchain, bkgchain, "match_DsFit_Kp_P", "best_DsFit_Kp_P", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Kp_P_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Km_P", "best_DsFit_Km_P", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", 0, 100, "1", "1", "./SBcompare/Km_P_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_pi_P", "best_DsFit_pi_P", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", 0, 100, "1", "1", "./SBcompare/pi_P_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_phi_P", "best_DsFit_phi_P", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", 0, 200, "1", "1", "./SBcompare/phi_P_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_Ds_P", "best_DsFit_Ds_P", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", 0, 250, "1", "1", "./SBcompare/Ds_P_DsFit.png");

    compare(sigchain, bkgchain, "match_dR_Kp_Km", "best_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", 0, 0.1, "1", "1", "./SBcompare/dR_Kp_Km.png");
    compare(sigchain, bkgchain, "match_dR_Kp_pi", "best_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Kp_pi.png");
    compare(sigchain, bkgchain, "match_dR_Kp_phi", "best_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Kp_phi.png");
    compare(sigchain, bkgchain, "match_dR_Kp_Ds", "best_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Kp_Ds.png");
    compare(sigchain, bkgchain, "match_dR_Km_pi", "best_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Km_pi.png");
    compare(sigchain, bkgchain, "match_dR_Km_phi", "best_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Km_phi.png");
    compare(sigchain, bkgchain, "match_dR_Km_Ds", "best_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Km_Ds.png");
    compare(sigchain, bkgchain, "match_dR_pi_phi", "best_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_phi.png");
    compare(sigchain, bkgchain, "match_dR_pi_Ds", "best_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_Ds.png");
    compare(sigchain, bkgchain, "match_dR_phi_Ds", "best_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_phi_Ds.png");

    compare(sigchain, bkgchain, "match_phiFit_dR_Kp_Km", "best_phiFit_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", 0, 0.1, "1", "1", "./SBcompare/dR_Kp_Km_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Kp_pi", "best_phiFit_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Kp_pi_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Kp_phi", "best_phiFit_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Kp_phi_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Kp_Ds", "best_phiFit_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Kp_Ds_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Km_pi", "best_phiFit_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Km_pi_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Km_phi", "best_phiFit_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Km_phi_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_Km_Ds", "best_phiFit_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Km_Ds_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_pi_phi", "best_phiFit_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_phi_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_pi_Ds", "best_phiFit_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_Ds_phiFit.png");
    compare(sigchain, bkgchain, "match_phiFit_dR_phi_Ds", "best_phiFit_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_phi_Ds_phiFit.png");

    compare(sigchain, bkgchain, "match_DsFit_dR_Kp_Km", "best_DsFit_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", 0, 0.1, "1", "1", "./SBcompare/dR_Kp_Km_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Kp_pi", "best_DsFit_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Kp_pi_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Kp_phi", "best_DsFit_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Kp_phi_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Kp_Ds", "best_DsFit_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Kp_Ds_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Km_pi", "best_DsFit_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_Km_pi_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Km_phi", "best_DsFit_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", 0, 0.05, "1", "1", "./SBcompare/dR_Km_phi_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_Km_Ds", "best_DsFit_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_Km_Ds_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_pi_phi", "best_DsFit_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_phi_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_pi_Ds", "best_DsFit_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", 0, 0.4, "1", "1", "./SBcompare/dR_pi_Ds_DsFit.png");
    compare(sigchain, bkgchain, "match_DsFit_dR_phi_Ds", "best_DsFit_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", 0, 0.15, "1", "1", "./SBcompare/dR_phi_Ds_DsFit.png");

    compare(sigchain, bkgchain, "match_dxy_Kp_Km", "best_dxy_Kp_Km", "signal", "W+jets", "#it{d_{xy}}(#it{K^{+}, K^{-}}) [cm]", 0, 0.05, "1", "1", "./SBcompare/dxy_Kp_Km.png");
    compare(sigchain, bkgchain, "match_dxy_Kp_pi", "best_dxy_Kp_pi", "signal", "W+jets", "#it{d_{xy}}(#it{K^{+}, #pi^{+}}) [cm]", 0, 0.1, "1", "1", "./SBcompare/dxy_Kp_pi.png");
    compare(sigchain, bkgchain, "match_dxy_Kp_phi", "best_dxy_Kp_phi", "signal", "W+jets", "#it{d_{xy}}(#it{K^{+}, #phi}) [cm]", 0, 1.5, "1", "1", "./SBcompare/dxy_Kp_phi.png");
    compare(sigchain, bkgchain, "match_dxy_Kp_Ds", "best_dxy_Kp_Ds", "signal", "W+jets", "#it{d_{xy}}(#it{K^{+}, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dxy_Kp_Ds.png");
    compare(sigchain, bkgchain, "match_dxy_Km_pi", "best_dxy_Km_pi", "signal", "W+jets", "#it{d_{xy}}(#it{K^{-}, #pi^{+}}) [cm]", 0, 0.1, "1", "1", "./SBcompare/dxy_Km_pi.png");
    compare(sigchain, bkgchain, "match_dxy_Km_phi", "best_dxy_Km_phi", "signal", "W+jets", "#it{d_{xy}}(#it{K^{-}, #phi}) [cm]", 0, 1.5, "1", "1", "./SBcompare/dxy_Km_phi.png");
    compare(sigchain, bkgchain, "match_dxy_Km_Ds", "best_dxy_Km_Ds", "signal", "W+jets", "#it{d_{xy}}(#it{K^{-}, D_{s}^{+}}) [cm]", 0, 0.1, "1", "1", "./SBcompare/dxy_Km_Ds.png");
    compare(sigchain, bkgchain, "match_dxy_pi_phi", "best_dxy_pi_phi", "signal", "W+jets", "#it{d_{xy}}(#it{#pi^{+}, #phi}) [cm]", 0, 1.5, "1", "1", "./SBcompare/dxy_pi_phi.png");
    compare(sigchain, bkgchain, "match_dxy_pi_Ds", "best_dxy_pi_Ds", "signal", "W+jets", "#it{d_{xy}}(#it{#pi^{+}, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dxy_pi_Ds.png");
    compare(sigchain, bkgchain, "match_dxy_phi_Ds", "best_dxy_phi_Ds", "signal", "W+jets", "#it{d_{xy}}(#it{#phi, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dxy_phi_Ds.png");

    compare(sigchain, bkgchain, "match_dz_Kp_Km", "best_dz_Kp_Km", "signal", "W+jets", "#it{d_{z}}(#it{K^{+}, K^{-}}) [cm]", 0, 0.1, "1", "1", "./SBcompare/dz_Kp_Km.png");
    compare(sigchain, bkgchain, "match_dz_Kp_pi", "best_dz_Kp_pi", "signal", "W+jets", "#it{d_{z}}(#it{K^{+}, #pi^{+}}) [cm]", 0, 0.2, "1", "1", "./SBcompare/dz_Kp_pi.png");
    compare(sigchain, bkgchain, "match_dz_Kp_phi", "best_dz_Kp_phi", "signal", "W+jets", "#it{d_{z}}(#it{K^{+}, #phi}) [cm]", 0, 2, "1", "1", "./SBcompare/dz_Kp_phi.png");
    compare(sigchain, bkgchain, "match_dz_Kp_Ds", "best_dz_Kp_Ds", "signal", "W+jets", "#it{d_{z}}(#it{K^{+}, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dz_Kp_Ds.png");
    compare(sigchain, bkgchain, "match_dz_Km_pi", "best_dz_Km_pi", "signal", "W+jets", "#it{d_{z}}(#it{K^{-}, #pi^{+}}) [cm]", 0, 0.2, "1", "1", "./SBcompare/dz_Km_pi.png");
    compare(sigchain, bkgchain, "match_dz_Km_phi", "best_dz_Km_phi", "signal", "W+jets", "#it{d_{z}}(#it{K^{-}, #phi}) [cm]", 0, 2, "1", "1", "./SBcompare/dz_Km_phi.png");
    compare(sigchain, bkgchain, "match_dz_Km_Ds", "best_dz_Km_Ds", "signal", "W+jets", "#it{d_{z}}(#it{K^{-}, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dz_Km_Ds.png");
    compare(sigchain, bkgchain, "match_dz_pi_phi", "best_dz_pi_phi", "signal", "W+jets", "#it{d_{z}}(#it{#pi^{+}, #phi}) [cm]", 0, 2, "1", "1", "./SBcompare/dz_pi_phi.png");
    compare(sigchain, bkgchain, "match_dz_pi_Ds", "best_dz_pi_Ds", "signal", "W+jets", "#it{d_{z}}(#it{#pi^{+}, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dz_pi_Ds.png");
    compare(sigchain, bkgchain, "match_dz_phi_Ds", "best_dz_phi_Ds", "signal", "W+jets", "#it{d_{z}}(#it{#phi, D_{s}^{+}}) [cm]", 0, 1, "1", "1", "./SBcompare/dz_phi_Ds.png");
    
    compare(sigchain, bkgchain, "match_Ds_FDxy", "best_Ds_FDxy", "signal", "W+jets", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", "1", "./SBcompare/Ds_FDxy.png");
    compare(sigchain, bkgchain, "match_Ds_FDz", "best_Ds_FDz", "signal", "W+jets", "FD_{#it{z}}(#it{D_{s}^{+}}) [cm]", 0, 2, "1", "1", "./SBcompare/Ds_FDz.png");
    compare(sigchain, bkgchain, "match_Ds_FD", "best_Ds_FD", "signal", "W+jets", "FD(#it{D_{s}^{+}}) [cm]", 0, 3, "1", "1", "./SBcompare/Ds_FD.png");

    compare(sigchain, bkgchain, "match_Ds_FDxy_Chi2", "best_Ds_FDxy_Chi2", "signal", "W+jets", "#chi^{2} FD_{#it{xy}}(#it{D_{s}^{+}})", 0, 60, "1", "1", "./SBcompare/Ds_FDxy_Chi2.png");
    compare(sigchain, bkgchain, "match_Ds_FDz_Chi2", "best_Ds_FDz_Chi2", "signal", "W+jets", "#chi^{2} FD_{#it{z}}(#it{D_{s}^{+}})", 0, 60, "1", "1", "./SBcompare/Ds_FDz_Chi2.png");
    compare(sigchain, bkgchain, "match_Ds_FD_Chi2", "best_Ds_FD_Chi2", "signal", "W+jets", "#chi^{2} FD(#it{D_{s}^{+}})", 0, 60, "1", "1", "./SBcompare/Ds_FD_Chi2.png");

    compare(sigchain, bkgchain, "match_Ds_DIRA_angle", "best_Ds_DIRA_angle", "signal", "W+jets", "#it{D_{s}^{+}} direction angle #it{#theta}", -1, 1, "1", "1", "./SBcompare/Ds_DIRA_angle.png");
    compare(sigchain, bkgchain, "match_Ds_DIRA", "best_Ds_DIRA", "signal", "W+jets", "cos(#it{#theta})", -1, 1, "1", "1", "./SBcompare/Ds_DIRA.png");

    compare(sigchain, bkgchain, "match_Kp_IP", "best_Kp_IP", "signal", "W+jets", "IP(#it{K^{+}}) [cm]", 0, 0.05, "1", "1", "./SBcompare/Kp_IP.png");
    compare(sigchain, bkgchain, "match_Km_IP", "best_Km_IP", "signal", "W+jets", "IP(#it{K^{-}}) [cm]", 0, 0.05, "1", "1", "./SBcompare/Km_IP.png");
    compare(sigchain, bkgchain, "match_pi_IP", "best_pi_IP", "signal", "W+jets", "IP(#it{#pi^{+}}) [cm]", 0, 0.1, "1", "1", "./SBcompare/pi_IP.png");

    compare(sigchain, bkgchain, "match_Kp_IP_Chi2", "best_Kp_IP_Chi2", "signal", "W+jets", "#chi^{2} IP(#it{K^{+}})", 0, 20, "1", "1", "./SBcompare/Kp_IP_Chi2.png");
    compare(sigchain, bkgchain, "match_Km_IP_Chi2", "best_Km_IP_Chi2", "signal", "W+jets", "#chi^{2} IP(#it{K^{-}})", 0, 20, "1", "1", "./SBcompare/Km_IP_Chi2.png");
    compare(sigchain, bkgchain, "match_pi_IP_Chi2", "best_pi_IP_Chi2", "signal", "W+jets", "#chi^{2} IP(#it{#pi^{+}})", 0, 40, "1", "1", "./SBcompare/pi_IP_Chi2.png");

    delete sigchain;
    delete bkgchain;

    return 0;
}

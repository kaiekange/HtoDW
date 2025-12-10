#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void compare(TChain *mychain1, TChain *mychain2, TString myvar1, TString myvar2, TString leg1, TString leg2, TString vartitle, bool logorno, int nbins, float varmin, float varmax, TCut cut1, TCut cut2, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", nbins, varmin, varmax);
    mychain1->Project("h1", myvar1, cut1);
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlack);
    h1->SetMarkerColor(kBlack);
    h1->SetMarkerSize(0.7);

    TH1F *h2 = new TH1F("h2", "", nbins, varmin, varmax);
    mychain2->Project("h2", myvar2, cut2);
    h2->Scale(1./h2->Integral());
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
    h2->SetMarkerSize(0.7);

    float height = std::max(h1->GetMaximum(), h2->GetMaximum());

    if(logorno) h1->SetMaximum(height*100);
    else{
        h1->SetMaximum(height*1.3);
        h1->SetMinimum(0.0);
    }
    h1->GetXaxis()->SetTitle(vartitle);
    h1->GetYaxis()->SetTitle("Normalized # candidates");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    if(logorno) canvas->SetLogy(1);
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
    canvas->SaveAs("./figures/sig_bkg_compare_muptcut/"+figpath+".png");
    /* canvas->SaveAs("./figures/sig_bkg_compare/"+figpath+".png"); */
    delete h1;
    delete h2;
}

int sig_bkg_compare() {
    setTDRStyle();

    /* ROOT::EnableImplicitMT(); */

    TChain *mychain1 = new TChain("PVStudy/Events");
    TChain *mychain2 = new TChain("RecoAnalyzer/Events");
    mychain1->Add("../../tuples/PVStudy.root");
    mychain2->Add("../../tuples/WJets_double_skim.root");

    TCut sigcut = "best_match_entry && best_mu_match && best_mu_pt > 15";
    TCut bkgcut = "best_mu_pt > 15";

    /* TCut sigcut = "best_match_entry && best_mu_match"; */
    /* TCut bkgcut = "1"; */

    compare(mychain1, mychain2, "best_Ds_IsoR03_sumChargedHadronPt", "best_Ds_IsoR03_sumChargedHadronPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumChargedHadronPt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumChargedHadronPt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumNeutralHadronEt", "best_Ds_IsoR03_sumNeutralHadronEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumNeutralHadronEt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumNeutralHadronEt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumPhotonEt", "best_Ds_IsoR03_sumPhotonEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumPhotonEt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumPhotonEt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumPUPt", "best_Ds_IsoR03_sumPUPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumPUPt", true, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR03_sumPUPt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR03_PFIso", "best_Ds_IsoR03_PFIso", "signal", "W+jets", "#it{D_{s}^{+}} RelIsoR03", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_PFIsoRel_log");
    compare(mychain1, mychain2, "best_Ds_IsoR03_PFIso * best_DsFit_Ds_pt", "best_Ds_IsoR03_PFIso * best_DsFit_Ds_pt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03", true, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR03_PFIso_log");
 
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumChargedHadronPt", "best_Ds_IsoR04_sumChargedHadronPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumChargedHadronPt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumChargedHadronPt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumNeutralHadronEt", "best_Ds_IsoR04_sumNeutralHadronEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumNeutralHadronEt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumNeutralHadronEt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumPhotonEt", "best_Ds_IsoR04_sumPhotonEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumPhotonEt", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumPhotonEt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumPUPt", "best_Ds_IsoR04_sumPUPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumPUPt", true, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR04_sumPUPt_log");
    compare(mychain1, mychain2, "best_Ds_IsoR04_PFIso", "best_Ds_IsoR04_PFIso", "signal", "W+jets", "#it{D_{s}^{+}} RelIsoR04", true, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_PFIsoRel_log");
    compare(mychain1, mychain2, "best_Ds_IsoR04_PFIso * best_DsFit_Ds_pt", "best_Ds_IsoR04_PFIso * best_DsFit_Ds_pt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04", true, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR04_PFIso_log");

    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumChargedHadronPt", "best_mu_IsoR03_sumChargedHadronPt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumChargedHadronPt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumChargedHadronPt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumNeutralHadronEt", "best_mu_IsoR03_sumNeutralHadronEt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumNeutralHadronEt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumNeutralHadronEt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumPhotonEt", "best_mu_IsoR03_sumPhotonEt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumPhotonEt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumPhotonEt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumPUPt", "best_mu_IsoR03_sumPUPt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumPUPt", true, 100, 0, 10, sigcut, bkgcut, "mu_IsoR03_sumPUPt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_PFIso", "best_mu_IsoR03_PFIso", "signal", "W+jets", "#it{#mu^{-}} RelIsoR03", true, 100, 0, 5, sigcut, bkgcut, "mu_IsoR03_PFIsoRel_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_PFIso * best_mu_pt", "best_mu_IsoR03_PFIso * best_mu_pt", "signal", "W+jets", "#it{#mu^{-}} IsoR03", true, 100, 0, 5, sigcut, bkgcut, "mu_IsoR03_PFIso_log"); */

    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumChargedHadronPt", "best_mu_IsoR04_sumChargedHadronPt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumChargedHadronPt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumChargedHadronPt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumNeutralHadronEt", "best_mu_IsoR04_sumNeutralHadronEt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumNeutralHadronEt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumNeutralHadronEt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumPhotonEt", "best_mu_IsoR04_sumPhotonEt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumPhotonEt", true, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumPhotonEt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumPUPt", "best_mu_IsoR04_sumPUPt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumPUPt", true, 100, 0, 10, sigcut, bkgcut, "mu_IsoR04_sumPUPt_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_PFIso", "best_mu_IsoR04_PFIso", "signal", "W+jets", "#it{#mu^{-}} IsoR04", true, 100, 0, 5, sigcut, bkgcut, "mu_IsoR04_PFIso_log"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_PFIso * best_mu_pt", "best_mu_IsoR03_PFIso * best_mu_pt", "signal", "W+jets", "#it{#mu^{-}} IsoR03", true, 100, 0, 5, sigcut, bkgcut, "mu_IsoR03_PFIso_log"); */

    /* compare(mychain1, mychain2, "best_Ds_primvtx_dira_angle", "best_Ds_primvtx_dira_angle", "signal", "W+jets", "#it{D_{s}^{+}} direction angle #it{#theta}", true, 100, 0, 3.1416, sigcut, bkgcut, "Ds_primvtx_dira_angle_log"); */
    /* compare(mychain1, mychain2, "best_Ds_primvtx_dira", "best_Ds_primvtx_dira", "signal", "W+jets", "cos(#it{#theta})", true, 100, -1, 1, sigcut, bkgcut, "Ds_primvtx_dira_log"); */
    /* compare(mychain1, mychain2, "best_Ds_PVnoDs_dira_angle", "best_Ds_PVnoDs_dira_angle", "signal", "W+jets", "#it{D_{s}^{+}} direction angle #it{#theta}", true, 100, 0, 3.1416, sigcut, bkgcut, "Ds_PVnoDs_dira_angle_log"); */
    /* compare(mychain1, mychain2, "best_Ds_PVnoDs_dira", "best_Ds_PVnoDs_dira", "signal", "W+jets", "cos(#it{#theta})", true, 100, -1, 1, sigcut, bkgcut, "Ds_PVnoDs_dira_log"); */

 
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumChargedHadronPt", "best_Ds_IsoR03_sumChargedHadronPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumChargedHadronPt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumChargedHadronPt");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumNeutralHadronEt", "best_Ds_IsoR03_sumNeutralHadronEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumNeutralHadronEt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumNeutralHadronEt");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumPhotonEt", "best_Ds_IsoR03_sumPhotonEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumPhotonEt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_sumPhotonEt");
    compare(mychain1, mychain2, "best_Ds_IsoR03_sumPUPt", "best_Ds_IsoR03_sumPUPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03 sumPUPt", false, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR03_sumPUPt");
    compare(mychain1, mychain2, "best_Ds_IsoR03_PFIso", "best_Ds_IsoR03_PFIso", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR03_PFIsoRel");
    compare(mychain1, mychain2, "best_Ds_IsoR03_PFIso * best_DsFit_Ds_pt", "best_Ds_IsoR03_PFIso * best_DsFit_Ds_pt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR03", false, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR03_PFIso");
 
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumChargedHadronPt", "best_Ds_IsoR04_sumChargedHadronPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumChargedHadronPt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumChargedHadronPt");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumNeutralHadronEt", "best_Ds_IsoR04_sumNeutralHadronEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumNeutralHadronEt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumNeutralHadronEt");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumPhotonEt", "best_Ds_IsoR04_sumPhotonEt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumPhotonEt", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_sumPhotonEt");
    compare(mychain1, mychain2, "best_Ds_IsoR04_sumPUPt", "best_Ds_IsoR04_sumPUPt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04 sumPUPt", false, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR04_sumPUPt");
    compare(mychain1, mychain2, "best_Ds_IsoR04_PFIso", "best_Ds_IsoR04_PFIso", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04", false, 100, 0, 5, sigcut, bkgcut, "Ds_IsoR04_PFIsoRel");
    compare(mychain1, mychain2, "best_Ds_IsoR04_PFIso * best_DsFit_Ds_pt", "best_Ds_IsoR04_PFIso * best_DsFit_Ds_pt", "signal", "W+jets", "#it{D_{s}^{+}} IsoR04", false, 100, 0, 10, sigcut, bkgcut, "Ds_IsoR04_PFIso");

    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumChargedHadronPt", "best_mu_IsoR03_sumChargedHadronPt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumChargedHadronPt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumChargedHadronPt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumNeutralHadronEt", "best_mu_IsoR03_sumNeutralHadronEt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumNeutralHadronEt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumNeutralHadronEt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumPhotonEt", "best_mu_IsoR03_sumPhotonEt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumPhotonEt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR03_sumPhotonEt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_sumPUPt", "best_mu_IsoR03_sumPUPt", "signal", "W+jets", "#it{#mu^{-}} IsoR03 sumPUPt", false, 100, 0, 10, sigcut, bkgcut, "mu_IsoR03_sumPUPt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_PFIso", "best_mu_IsoR03_PFIso", "signal", "W+jets", "#it{#mu^{-}} RelIsoR03", false, 100, 0, 5, sigcut, bkgcut, "mu_IsoR03_PFIsoRel"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR03_PFIso * best_mu_pt", "best_mu_IsoR03_PFIso * best_mu_pt", "signal", "W+jets", "#it{#mu^{-}} IsoR03", false, 100, 0, 5, sigcut, bkgcut, "mu_IsoR03_PFIso"); */

    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumChargedHadronPt", "best_mu_IsoR04_sumChargedHadronPt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumChargedHadronPt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumChargedHadronPt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumNeutralHadronEt", "best_mu_IsoR04_sumNeutralHadronEt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumNeutralHadronEt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumNeutralHadronEt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumPhotonEt", "best_mu_IsoR04_sumPhotonEt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumPhotonEt", false, 100, 0, 2, sigcut, bkgcut, "mu_IsoR04_sumPhotonEt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_sumPUPt", "best_mu_IsoR04_sumPUPt", "signal", "W+jets", "#it{#mu^{-}} IsoR04 sumPUPt", false, 100, 0, 10, sigcut, bkgcut, "mu_IsoR04_sumPUPt"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_PFIso", "best_mu_IsoR04_PFIso", "signal", "W+jets", "#it{#mu^{-}} RelIsoR04", false, 100, 0, 5, sigcut, bkgcut, "mu_IsoR04_PFIsoRel"); */
    /* compare(mychain1, mychain2, "best_mu_IsoR04_PFIso * best_mu_pt", "best_mu_IsoR04_PFIso * best_mu_pt", "signal", "W+jets", "#it{#mu^{-}} IsoR04", false, 100, 0, 5, sigcut, bkgcut, "mu_IsoR04_PFIso"); */

    /* compare(mychain1, mychain2, "best_Ds_primvtx_dira_angle", "best_Ds_primvtx_dira_angle", "signal", "W+jets", "#it{D_{s}^{+}} direction angle #it{#theta}", false, 100, 0, 3.1416, sigcut, bkgcut, "Ds_primvtx_dira_angle"); */
    /* compare(mychain1, mychain2, "best_Ds_primvtx_dira", "best_Ds_primvtx_dira", "signal", "W+jets", "cos(#it{#theta})", false, 100, -1, 1, sigcut, bkgcut, "Ds_primvtx_dira"); */
    /* compare(mychain1, mychain2, "best_Ds_PVnoDs_dira_angle", "best_Ds_PVnoDs_dira_angle", "signal", "W+jets", "#it{D_{s}^{+}} direction angle #it{#theta}", false, 100, 0, 3.1416, sigcut, bkgcut, "Ds_PVnoDs_dira_angle"); */
    /* compare(mychain1, mychain2, "best_Ds_PVnoDs_dira", "best_Ds_PVnoDs_dira", "signal", "W+jets", "cos(#it{#theta})", false, 100, -1, 1, sigcut, bkgcut, "Ds_PVnoDs_dira"); */

    /* compare(mychain1, mychain2, "best_DsFit_Ds_invm", "best_DsFit_Ds_invm", "signal", "W+jets", "#it{M}(#it{K^{+}K^{-}#pi^{+}}) [GeV]", false, 100, 1.85, 2.1, sigcut, bkgcut, "Ds_invm"); */
    /* compare(mychain1, mychain2, "best_phiFit_phi_invm", "best_phiFit_phi_invm", "signal", "W+jets", "#it{M}(#it{K^{+}K^{-}}) [GeV]", false, 100, 0.99, 1.05, sigcut, bkgcut, "phi_invm"); */

    /* compare(mychain1, mychain2, "best_Kp_eta", "best_Kp_eta", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Kp_eta"); */
    /* compare(mychain1, mychain2, "best_Km_eta", "best_Km_eta", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", false, 100, -3, 3, sigcut, bkgcut, "Km_eta"); */
    /* compare(mychain1, mychain2, "best_pi_eta", "best_pi_eta", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", false, 100, -3, 3, sigcut, bkgcut, "pi_eta"); */
    /* compare(mychain1, mychain2, "best_phi_eta", "best_phi_eta", "signal", "W+jets", "#it{#eta}(#it{#phi})", false, 100, -3, 3, sigcut, bkgcut, "phi_eta"); */
    /* compare(mychain1, mychain2, "best_Ds_eta", "best_Ds_eta", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Ds_eta"); */

    /* compare(mychain1, mychain2, "best_phiFit_Kp_eta", "best_phiFit_Kp_eta", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Kp_eta_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Km_eta", "best_phiFit_Km_eta", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", false, 100, -3, 3, sigcut, bkgcut, "Km_eta_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_pi_eta", "best_phiFit_pi_eta", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", false, 100, -3, 3, sigcut, bkgcut, "pi_eta_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_phi_eta", "best_phiFit_phi_eta", "signal", "W+jets", "#it{#eta}(#it{#phi})", false, 100, -3, 3, sigcut, bkgcut, "phi_eta_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Ds_eta", "best_phiFit_Ds_eta", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Ds_eta_phiFit"); */

    /* compare(mychain1, mychain2, "best_DsFit_Kp_eta", "best_DsFit_Kp_eta", "signal", "W+jets", "#it{#eta}(#it{K^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Kp_eta_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Km_eta", "best_DsFit_Km_eta", "signal", "W+jets", "#it{#eta}(#it{K^{-}})", false, 100, -3, 3, sigcut, bkgcut, "Km_eta_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_pi_eta", "best_DsFit_pi_eta", "signal", "W+jets", "#it{#eta}(#it{#pi^{+}})", false, 100, -3, 3, sigcut, bkgcut, "pi_eta_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_phi_eta", "best_DsFit_phi_eta", "signal", "W+jets", "#it{#eta}(#it{#phi})", false, 100, -3, 3, sigcut, bkgcut, "phi_eta_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Ds_eta", "best_DsFit_Ds_eta", "signal", "W+jets", "#it{#eta}(#it{D_{s}^{+}})", false, 100, -3, 3, sigcut, bkgcut, "Ds_eta_DsFit"); */

    /* compare(mychain1, mychain2, "best_Kp_phi", "best_Kp_phi", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Kp_phi"); */
    /* compare(mychain1, mychain2, "best_Km_phi", "best_Km_phi", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Km_phi"); */
    /* compare(mychain1, mychain2, "best_pi_phi", "best_pi_phi", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "pi_phi"); */
    /* compare(mychain1, mychain2, "best_phi_phi", "best_phi_phi", "signal", "W+jets", "#it{#phi}(#it{#phi})", false, 100, -3.15, 3.15, sigcut, bkgcut, "phi_phi"); */
    /* compare(mychain1, mychain2, "best_Ds_phi", "best_Ds_phi", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Ds_phi"); */

    /* compare(mychain1, mychain2, "best_phiFit_Kp_phi", "best_phiFit_Kp_phi", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Kp_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Km_phi", "best_phiFit_Km_phi", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Km_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_pi_phi", "best_phiFit_pi_phi", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "pi_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_phi_phi", "best_phiFit_phi_phi", "signal", "W+jets", "#it{#phi}(#it{#phi})", false, 100, -3.15, 3.15, sigcut, bkgcut, "phi_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Ds_phi", "best_phiFit_Ds_phi", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Ds_phi_phiFit"); */

    /* compare(mychain1, mychain2, "best_DsFit_Kp_phi", "best_DsFit_Kp_phi", "signal", "W+jets", "#it{#phi}(#it{K^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Kp_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Km_phi", "best_DsFit_Km_phi", "signal", "W+jets", "#it{#phi}(#it{K^{-}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Km_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_pi_phi", "best_DsFit_pi_phi", "signal", "W+jets", "#it{#phi}(#it{#pi^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "pi_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_phi_phi", "best_DsFit_phi_phi", "signal", "W+jets", "#it{#phi}(#it{#phi})", false, 100, -3.15, 3.15, sigcut, bkgcut, "phi_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Ds_phi", "best_DsFit_Ds_phi", "signal", "W+jets", "#it{#phi}(#it{D_{s}^{+}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "Ds_phi_DsFit"); */

    /* compare(mychain1, mychain2, "best_Kp_pt", "best_Kp_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Kp_pt"); */
    /* compare(mychain1, mychain2, "best_Km_pt", "best_Km_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Km_pt"); */
    /* compare(mychain1, mychain2, "best_pi_pt", "best_pi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "pi_pt"); */
    /* compare(mychain1, mychain2, "best_phi_pt", "best_phi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", false, 100, 0, 80, sigcut, bkgcut, "phi_pt"); */
    /* compare(mychain1, mychain2, "best_Ds_pt", "best_Ds_pt", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Ds_pt"); */

    /* compare(mychain1, mychain2, "best_phiFit_Kp_pt", "best_phiFit_Kp_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Kp_pt_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Km_pt", "best_phiFit_Km_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Km_pt_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_pi_pt", "best_phiFit_pi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "pi_pt_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_phi_pt", "best_phiFit_phi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", false, 100, 0, 80, sigcut, bkgcut, "phi_pt_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Ds_pt", "best_phiFit_Ds_pt", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Ds_pt_phiFit"); */

    /* compare(mychain1, mychain2, "best_DsFit_Kp_pt", "best_DsFit_Kp_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Kp_pt_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Km_pt", "best_DsFit_Km_pt", "signal", "W+jets", "#it{p_{T}}(#it{K^{-}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "Km_pt_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_pi_pt", "best_DsFit_pi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#pi^{+}}) [GeV]", false, 100, 0, 50, sigcut, bkgcut, "pi_pt_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_phi_pt", "best_DsFit_phi_pt", "signal", "W+jets", "#it{p_{T}}(#it{#phi}) [GeV]", false, 100, 0, 80, sigcut, bkgcut, "phi_pt_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Ds_pt", "best_DsFit_Ds_pt", "signal", "W+jets", "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Ds_pt_DsFit"); */

    /* compare(mychain1, mychain2, "best_Kp_p", "best_Kp_p", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Kp_p"); */
    /* compare(mychain1, mychain2, "best_Km_p", "best_Km_p", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Km_p"); */
    /* compare(mychain1, mychain2, "best_pi_p", "best_pi_p", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "pi_p"); */
    /* compare(mychain1, mychain2, "best_phi_p", "best_phi_p", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", false, 100, 0, 200, sigcut, bkgcut, "phi_p"); */
    /* compare(mychain1, mychain2, "best_Ds_p", "best_Ds_p", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 250, sigcut, bkgcut, "Ds_p"); */

    /* compare(mychain1, mychain2, "best_phiFit_Kp_p", "best_phiFit_Kp_p", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Kp_p_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Km_p", "best_phiFit_Km_p", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Km_p_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_pi_p", "best_phiFit_pi_p", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "pi_p_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_phi_p", "best_phiFit_phi_p", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", false, 100, 0, 200, sigcut, bkgcut, "phi_p_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_Ds_p", "best_phiFit_Ds_p", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 250, sigcut, bkgcut, "Ds_p_phiFit"); */

    /* compare(mychain1, mychain2, "best_DsFit_Kp_p", "best_DsFit_Kp_p", "signal", "W+jets", "#it{p}(#it{K^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Kp_p_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Km_p", "best_DsFit_Km_p", "signal", "W+jets", "#it{p}(#it{K^{-}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "Km_p_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_pi_p", "best_DsFit_pi_p", "signal", "W+jets", "#it{p}(#it{#pi^{+}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "pi_p_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_phi_p", "best_DsFit_phi_p", "signal", "W+jets", "#it{p}(#it{#phi}) [GeV]", false, 100, 0, 200, sigcut, bkgcut, "phi_p_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_Ds_p", "best_DsFit_Ds_p", "signal", "W+jets", "#it{p}(#it{D_{s}^{+}}) [GeV]", false, 100, 0, 250, sigcut, bkgcut, "Ds_p_DsFit"); */

    /* compare(mychain1, mychain2, "best_dR_Kp_Km", "best_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", false, 100, 0, 0.1, sigcut, bkgcut, "dR_Kp_Km"); */
    /* compare(mychain1, mychain2, "best_dR_Kp_pi", "best_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Kp_pi"); */
    /* compare(mychain1, mychain2, "best_dR_Kp_phi", "best_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Kp_phi"); */
    /* compare(mychain1, mychain2, "best_dR_Kp_Ds", "best_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Kp_Ds"); */
    /* compare(mychain1, mychain2, "best_dR_Km_pi", "best_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Km_pi"); */
    /* compare(mychain1, mychain2, "best_dR_Km_phi", "best_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Km_phi"); */
    /* compare(mychain1, mychain2, "best_dR_Km_Ds", "best_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Km_Ds"); */
    /* compare(mychain1, mychain2, "best_dR_pi_phi", "best_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_phi"); */
    /* compare(mychain1, mychain2, "best_dR_pi_Ds", "best_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_Ds"); */
    /* compare(mychain1, mychain2, "best_dR_phi_Ds", "best_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_phi_Ds"); */

    /* compare(mychain1, mychain2, "best_phiFit_dR_Kp_Km", "best_phiFit_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", false, 100, 0, 0.1, sigcut, bkgcut, "dR_Kp_Km_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Kp_pi", "best_phiFit_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Kp_pi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Kp_phi", "best_phiFit_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Kp_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Kp_Ds", "best_phiFit_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Kp_Ds_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Km_pi", "best_phiFit_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Km_pi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Km_phi", "best_phiFit_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Km_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_Km_Ds", "best_phiFit_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Km_Ds_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_pi_phi", "best_phiFit_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_phi_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_pi_Ds", "best_phiFit_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_Ds_phiFit"); */
    /* compare(mychain1, mychain2, "best_phiFit_dR_phi_Ds", "best_phiFit_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_phi_Ds_phiFit"); */

    /* compare(mychain1, mychain2, "best_DsFit_dR_Kp_Km", "best_DsFit_dR_Kp_Km", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, K^{-}})", false, 100, 0, 0.1, sigcut, bkgcut, "dR_Kp_Km_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Kp_pi", "best_DsFit_dR_Kp_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Kp_pi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Kp_phi", "best_DsFit_dR_Kp_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Kp_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Kp_Ds", "best_DsFit_dR_Kp_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{+}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Kp_Ds_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Km_pi", "best_DsFit_dR_Km_pi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #pi^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_Km_pi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Km_phi", "best_DsFit_dR_Km_phi", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, #phi})", false, 100, 0, 0.05, sigcut, bkgcut, "dR_Km_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_Km_Ds", "best_DsFit_dR_Km_Ds", "signal", "W+jets", "#Delta#it{R}(#it{K^{-}, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_Km_Ds_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_pi_phi", "best_DsFit_dR_pi_phi", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, #phi})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_phi_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_pi_Ds", "best_DsFit_dR_pi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#pi^{+}, D_{s}^{+}})", false, 100, 0, 0.4, sigcut, bkgcut, "dR_pi_Ds_DsFit"); */
    /* compare(mychain1, mychain2, "best_DsFit_dR_phi_Ds", "best_DsFit_dR_phi_Ds", "signal", "W+jets", "#Delta#it{R}(#it{#phi, D_{s}^{+}})", false, 100, 0, 0.15, sigcut, bkgcut, "dR_phi_Ds_DsFit"); */

    /* compare(mychain1, mychain2, "best_dxy_phi_Ds", "best_dxy_phi_Ds", "signal", "W+jets", "#it{d_{xy}}(#it{#phi, D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "dxy_phi_Ds"); */
    /* compare(mychain1, mychain2, "best_dz_phi_Ds", "best_dz_phi_Ds", "signal", "W+jets", "#it{d_{z}}(#it{#phi, D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "dz_phi_Ds"); */

    /* compare(mychain1, mychain2, "best_phiFit_chi2ndof", "best_phiFit_chi2ndof", "signal", "W+jets", "#it{#phi} fit #it{#chi^{2}}/NDOF", false, 100, 0, 10, sigcut, bkgcut, "phiFit_chi2ndof"); */
    /* compare(mychain1, mychain2, "best_DsFit_chi2ndof", "best_DsFit_chi2ndof", "signal", "W+jets", "#it{D_{s}^{+}} fit #it{#chi^{2}}/NDOF", false, 100, 0, 10, sigcut, bkgcut, "DsFit_chi2ndof"); */

    /* compare(mychain1, mychain1, "best_Ds_primvtx_FDxy", "best_Ds_primvtx_FDxy", "signal", "W+jets", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_primvtx_FDxy"); */
    /* compare(mychain1, mychain1, "best_Ds_primvtx_FDz", "best_Ds_primvtx_FDz", "signal", "W+jets", "FD_{#it{z}}(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_primvtx_FDz"); */
    /* compare(mychain1, mychain1, "best_Ds_primvtx_FD", "best_Ds_primvtx_FD", "signal", "W+jets", "FD(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_primvtx_FD"); */

    /* compare(mychain1, mychain1, "best_Ds_primvtx_FDxychi2", "best_Ds_primvtx_FDxychi2", "signal", "W+jets", "FD_{#it{xy}}(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_primvtx_FDxychi2"); */
    /* compare(mychain1, mychain1, "best_Ds_primvtx_FDzchi2", "best_Ds_primvtx_FDzchi2", "signal", "W+jets", "FD_{#it{z}}(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_primvtx_FDzchi2"); */
    /* compare(mychain1, mychain1, "best_Ds_primvtx_FDchi2", "best_Ds_primvtx_FDchi2", "signal", "W+jets", "FD(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_primvtx_FDchi2"); */

    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FDxy", "best_Ds_PVwithDs_FDxy", "signal", "W+jets", "FD_{#it{xy}}(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_PVwithDs_FDxy"); */
    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FDz", "best_Ds_PVwithDs_FDz", "signal", "W+jets", "FD_{#it{z}}(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_PVwithDs_FDz"); */
    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FD", "best_Ds_PVwithDs_FD", "signal", "W+jets", "FD(#it{D_{s}^{+}}) [cm]", false, 100, 0, 1, sigcut, bkgcut, "Ds_PVwithDs_FD"); */

    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FDxychi2", "best_Ds_PVwithDs_FDxychi2", "signal", "W+jets", "FD_{#it{xy}}(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_PVwithDs_FDxychi2"); */
    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FDzchi2", "best_Ds_PVwithDs_FDzchi2", "signal", "W+jets", "FD_{#it{z}}(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_PVwithDs_FDzchi2"); */
    /* compare(mychain1, mychain1, "best_Ds_PVwithDs_FDchi2", "best_Ds_PVwithDs_FDchi2", "signal", "W+jets", "FD(#it{D_{s}^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Ds_PVwithDs_FDchi2"); */

    /* compare(mychain1, mychain2, "best_Kp_primvtx_ip", "best_Kp_primvtx_ip", "signal", "W+jets", "IP(#it{K^{+}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Kp_primvtx_ip"); */
    /* compare(mychain1, mychain2, "best_Km_primvtx_ip", "best_Km_primvtx_ip", "signal", "W+jets", "IP(#it{K^{-}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Km_primvtx_ip"); */
    /* compare(mychain1, mychain2, "best_pi_primvtx_ip", "best_pi_primvtx_ip", "signal", "W+jets", "IP(#it{#pi^{+}}) [cm]", false, 100, 0, 0.1, sigcut, bkgcut, "pi_primvtx_ip"); */
    /* compare(mychain1, mychain2, "best_phi_primvtx_ip", "best_phi_primvtx_ip", "signal", "W+jets", "IP(#it{#phi}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "phi_primvtx_ip"); */
    /* compare(mychain1, mychain2, "best_Ds_primvtx_ip", "best_Ds_primvtx_ip", "signal", "W+jets", "IP(#it{D_{s}^{+}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Ds_primvtx_ip"); */

    /* compare(mychain1, mychain2, "best_Kp_primvtx_ipchi2", "best_Kp_primvtx_ipchi2", "signal", "W+jets", "IP(#it{K^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Kp_primvtx_ipchi2"); */
    /* compare(mychain1, mychain2, "best_Km_primvtx_ipchi2", "best_Km_primvtx_ipchi2", "signal", "W+jets", "IP(#it{K^{-}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Km_primvtx_ipchi2"); */
    /* compare(mychain1, mychain2, "best_pi_primvtx_ipchi2", "best_pi_primvtx_ipchi2", "signal", "W+jets", "IP(#it{#pi^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "pi_primvtx_ipchi2"); */
    /* compare(mychain1, mychain2, "best_phi_primvtx_ipchi2", "best_phi_primvtx_ipchi2", "signal", "W+jets", "IP(#it{#phi}) significance", false, 100, 0, 0.05, sigcut, bkgcut, "phi_primvtx_ipchi2"); */
    /* compare(mychain1, mychain2, "best_Ds_primvtx_ipchi2", "best_Ds_primvtx_ipchi2", "signal", "W+jets", "IP(#it{D_{s}^{+}}) significance", false, 100, 0, 0.02, sigcut, bkgcut, "Ds_primvtx_ipchi2"); */

    /* compare(mychain1, mychain2, "best_Kp_PVwithDs_ip", "best_Kp_PVwithDs_ip", "signal", "W+jets", "IP(#it{K^{+}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Kp_PVwithDs_ip"); */
    /* compare(mychain1, mychain2, "best_Km_PVwithDs_ip", "best_Km_PVwithDs_ip", "signal", "W+jets", "IP(#it{K^{-}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Km_PVwithDs_ip"); */
    /* compare(mychain1, mychain2, "best_pi_PVwithDs_ip", "best_pi_PVwithDs_ip", "signal", "W+jets", "IP(#it{#pi^{+}}) [cm]", false, 100, 0, 0.1, sigcut, bkgcut, "pi_PVwithDs_ip"); */
    /* compare(mychain1, mychain2, "best_phi_PVwithDs_ip", "best_phi_PVwithDs_ip", "signal", "W+jets", "IP(#it{#phi}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "phi_PVwithDs_ip"); */
    /* compare(mychain1, mychain2, "best_Ds_PVwithDs_ip", "best_Ds_PVwithDs_ip", "signal", "W+jets", "IP(#it{D_{s}^{+}}) [cm]", false, 100, 0, 0.05, sigcut, bkgcut, "Ds_PVwithDs_ip"); */

    /* compare(mychain1, mychain2, "best_Kp_PVwithDs_ipchi2", "best_Kp_PVwithDs_ipchi2", "signal", "W+jets", "IP(#it{K^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Kp_PVwithDs_ipchi2"); */
    /* compare(mychain1, mychain2, "best_Km_PVwithDs_ipchi2", "best_Km_PVwithDs_ipchi2", "signal", "W+jets", "IP(#it{K^{-}}) significance", false, 100, 0, 20, sigcut, bkgcut, "Km_PVwithDs_ipchi2"); */
    /* compare(mychain1, mychain2, "best_pi_PVwithDs_ipchi2", "best_pi_PVwithDs_ipchi2", "signal", "W+jets", "IP(#it{#pi^{+}}) significance", false, 100, 0, 20, sigcut, bkgcut, "pi_PVwithDs_ipchi2"); */
    /* compare(mychain1, mychain2, "best_phi_PVwithDs_ipchi2", "best_phi_PVwithDs_ipchi2", "signal", "W+jets", "IP(#it{#phi}) significance", false, 100, 0, 0.05, sigcut, bkgcut, "phi_PVwithDs_ipchi2"); */
    /* compare(mychain1, mychain2, "best_Ds_PVwithDs_ipchi2", "best_Ds_PVwithDs_ipchi2", "signal", "W+jets", "IP(#it{D_{s}^{+}}) significance", false, 100, 0, 0.02, sigcut, bkgcut, "Ds_PVwithDs_ipchi2"); */

    /* compare(mychain1, mychain2, "best_Ds_primvtx_FDxy * best_DsFit_Ds_invm * 20 / best_DsFit_Ds_p", "best_Ds_primvtx_FDxy * best_DsFit_Ds_invm * 20 / best_DsFit_Ds_p", "signal", "W+jets", "#tau(#it{D_{s}^{+}}) [ps]", false, 100, 0, 2, sigcut, bkgcut, "Ds_primvtx_tau"); */
    /* compare(mychain1, mychain2, "best_Ds_PVwithDs_FDxy * best_DsFit_Ds_invm * 20 / best_DsFit_Ds_p", "best_Ds_PVwithDs_FDxy * best_DsFit_Ds_invm * 20 / best_DsFit_Ds_p", "signal", "W+jets", "#tau(#it{D_{s}^{+}}) [ps]", false, 100, 0, 2, sigcut, bkgcut, "Ds_PVwithDs_tau"); */

    /* compare(mychain1, mychain2, "best_mu_eta", "best_mu_eta", "signal", "W+jets", "#it{#eta}(#it{#mu^{-}})", false, 100, -3, 3, sigcut, bkgcut, "mu_eta"); */
    /* compare(mychain1, mychain2, "best_mu_phi", "best_mu_phi", "signal", "W+jets", "#it{#phi}(#it{#mu^{-}})", false, 100, -3.15, 3.15, sigcut, bkgcut, "mu_phi"); */
    /* compare(mychain1, mychain2, "best_mu_pt", "best_mu_pt", "signal", "W+jets", "#it{p_{T}}(#it{#mu^{-}}) [GeV]", false, 100, 0, 100, sigcut, bkgcut, "mu_pt"); */
    /* compare(mychain1, mychain2, "best_mu_p", "best_mu_p", "signal", "W+jets", "#it{p}(#it{#mu^{-}}) [GeV]", false, 100, 0, 200, sigcut, bkgcut, "mu_p"); */
    /* compare(mychain1, mychain2, "best_mu_primvtx_ip", "best_mu_primvtx_ip", "signal", "W+jets", "IP(#it{#mu^{-}}) [cm]", false, 100, 0, 0.02, sigcut, bkgcut, "mu_primvtx_ip"); */
    /* compare(mychain1, mychain2, "best_mu_primvtx_ipchi2", "best_mu_primvtx_ipchi2", "signal", "W+jets", "IP(#it{#mu^{-}}) significance", false, 100, 0, 5, sigcut, bkgcut, "mu_primvtx_ipchi2"); */

    /* compare(mychain1, mychain2, "best_mu_dxy", "best_mu_dxy", "signal", "W+jets", "#it{d_{xy}}(#it{#mu^{-}}) [cm]", false, 100, 0, 0.005, sigcut, bkgcut, "mu_dxy"); */
    /* compare(mychain1, mychain2, "best_mu_dz", "best_mu_dz", "signal", "W+jets", "#it{d_{z}}(#it{#mu^{-}}) [cm]", false, 100, 0, 0.02, sigcut, bkgcut, "mu_dz"); */

    /* compare(mychain1, mychain2, "best_mu_isHighPt", "best_mu_isHighPt", "signal", "W+jets", "isHighPtMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isHighPt"); */
    /* compare(mychain1, mychain2, "best_mu_isLoose", "best_mu_isLoose", "signal", "W+jets", "isLooseMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isLoose"); */
    /* compare(mychain1, mychain2, "best_mu_isMedium", "best_mu_isMedium", "signal", "W+jets", "isMediumMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isMedium"); */
    /* compare(mychain1, mychain2, "best_mu_isSoft", "best_mu_isSoft", "signal", "W+jets", "isSoftMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isSoft"); */
    /* compare(mychain1, mychain2, "best_mu_isPF", "best_mu_isPF", "signal", "W+jets", "isPFMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isPF"); */
    /* compare(mychain1, mychain2, "best_mu_isTracker", "best_mu_isTracker", "signal", "W+jets", "isTrackerMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isTracker"); */
    /* compare(mychain1, mychain2, "best_mu_isGlobal", "best_mu_isGlobal", "signal", "W+jets", "isGlobalMuon", false, 2, 0, 2, sigcut, bkgcut, "mu_isGlobal"); */

    delete mychain1;
    delete mychain2;

    return 0;
}

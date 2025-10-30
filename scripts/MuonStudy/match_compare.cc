#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void compare(TChain *mychain, TString myvar, TString vartitle, int nbins, float varmin, float varmax, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", nbins, varmin, varmax);
    mychain->Project("h1", myvar, "mu_match");
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlack);
    h1->SetMarkerColor(kBlack);
    h1->SetMarkerSize(0.7);

    TH1F *h2 = new TH1F("h2", "", nbins, varmin, varmax);
    mychain->Project("h2", myvar, "!mu_match");
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
    mylegend->AddEntry(h1, "match", "l");
    mylegend->AddEntry(h2, "non-match", "l");
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

    TChain *mychain = new TChain("MuonStudy/Events");
    /* mychain->Add("../../tuples/MuonStudy.root"); */
    mychain->Add("../../tuples/Muon2.root");

    /* compare(mychain, "mu_eta", "#it{#eta}(#it{#mu^{-}})", 100, -3, 3, "mu_eta"); */
    /* compare(mychain, "mu_phi", "#it{#phi}(#it{#mu^{-}})", 100, -3.1416, 3.1416, "mu_phi"); */
    compare(mychain, "mu_p", "#it{p}(#it{#mu^{-}}) [GeV]", 100, 0, 200, "mu_p");
    compare(mychain, "mu_pt", "#it{p_{T}}(#it{#mu^{-}}) [GeV]", 100, 0, 100, "mu_pt");
    /* compare(mychain, "mu_isHighPt", "isHighPtMuon", 2, 0, 2, "mu_isHighPt"); */
    /* compare(mychain, "mu_isMedium", "isMediumMuon", 2, 0, 2, "mu_isMedium"); */
    /* compare(mychain, "mu_isSoft", "isSoftMuon", 2, 0, 2, "mu_isSoft"); */
    /* compare(mychain, "mu_isTight", "isTightMuon", 2, 0, 2, "mu_isTight"); */
    /* compare(mychain, "mu_isPF", "isPFMuon", 2, 0, 2, "mu_isPF"); */
    /* compare(mychain, "mu_isTracker", "isTrackerMuon", 2, 0, 2, "mu_isTracker"); */
    /* compare(mychain, "mu_isGlobal", "isGlobalMuon", 2, 0, 2, "mu_isGlobal"); */
    /* compare(mychain, "mu_IsoR03_sumChargedHadronPt", "IsoR03 sumChargedHadronPt", 100, 0, 3, "mu_IsoR03_sumChargedHadronPt"); */
    /* compare(mychain, "mu_IsoR03_sumChargedParticlePt", "IsoR03 sumChargedParticlePt", 100, 0, 3, "mu_IsoR03_sumChargedParticlePt"); */
    /* compare(mychain, "mu_IsoR03_sumNeutralHadronEt", "IsoR03 sumNeutralHadronEt", 100, 0, 3, "mu_IsoR03_sumNeutralHadronEt"); */
    /* compare(mychain, "mu_IsoR03_sumPhotonEt", "IsoR03 sumPhotonEt", 100, 0, 3, "mu_IsoR03_sumPhotonEt"); */
    /* compare(mychain, "mu_IsoR03_sumPUPt", "IsoR03 sumPUPt", 100, 0, 10, "mu_IsoR03_sumPUPt"); */
    compare(mychain, "mu_PFIsoR03", "PFIsoR03", 100, 0, 0.3, "mu_IsoR03");
    /* compare(mychain, "mu_IsoR04_sumChargedHadronPt", "IsoR04 sumChargedHadronPt", 100, 0, 3, "mu_IsoR04_sumChargedHadronPt"); */
    /* compare(mychain, "mu_IsoR04_sumChargedParticlePt", "IsoR04 sumChargedParticlePt", 100, 0, 3, "mu_IsoR04_sumChargedParticlePt"); */
    /* compare(mychain, "mu_IsoR04_sumNeutralHadronEt", "IsoR04 sumNeutralHadronEt", 100, 0, 3, "mu_IsoR04_sumNeutralHadronEt"); */
    /* compare(mychain, "mu_IsoR04_sumPhotonEt", "IsoR04 sumPhotonEt", 100, 0, 3, "mu_IsoR04_sumPhotonEt"); */
    /* compare(mychain, "mu_IsoR04_sumPUPt", "IsoR04 sumPUPt", 100, 0, 10, "mu_IsoR04_sumPUPt"); */
    compare(mychain, "mu_PFIsoR04", "PFIsoR04", 100, 0, 0.3, "mu_IsoR04");
    /* compare(mychain, "mu_primvtx_dxy", "#it{d_{xy}}(#it{#mu^{-}}) [cm]", 100, -0.05, 0.05, "mu_primvtx_dxy"); */
    /* compare(mychain, "mu_primvtx_dz", "#it{d_{z}}(#it{#mu^{-}}) [cm]", 100, -0.2, 0.2, "mu_primvtx_dz"); */
    /* compare(mychain, "mu_primvtx_ip", "IP(#it{#mu^{-}}) [cm]", 100, 0, 0.05, "mu_primvtx_ip"); */
    /* compare(mychain, "mu_primvtx_ipchi2", "#it{#chi^{2}} IP(#it{#mu^{-}})", 100, 0, 5, "mu_primvtx_ipchi2"); */

    delete mychain;

    return 0;
}

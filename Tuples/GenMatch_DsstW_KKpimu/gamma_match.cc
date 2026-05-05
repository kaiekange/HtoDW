#include "TFile.h"
#include "TTree.h"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/tdrStyle.cc"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/CMS_lumi.cc"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/draw_funcs.cc"

Int_t gamma_match()
{
    setTDRStyle();
    lumi_sqrtS = "13.6 TeV, 2022";

    ROOT::EnableImplicitMT();
    TFile *myfile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/GenMatch_DsstW_KKpimu/GenMatch_DsstW_KKpimu.root");
    TTree *mytree = (TTree *)myfile->Get("GenMatchDsstWKKpimu/Events");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    TH1F *h1 = new TH1F("h1", ";#gamma #DeltaR(Reco,Gen);# events", 100, 0, 0.3);
    mytree->Project("h1", "gamma_deltaR");
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(2);
    h1->Draw("hist");
    h1->SetMaximum(1.3*h1->GetMaximum());
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("figures/dR.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvas_setup(c2);
    TH2F *h2 = new TH2F("h2", ";#gamma #Delta#it{p_{T}}/Gen #it{p_{T}};#gamma #DeltaR(Reco,Gen)", 100, 0, 2, 100, 0, 0.1);
    mytree->Project("h2", "gamma_deltaR:abs(Gen_gamma_pt - gamma_pt)/Gen_gamma_pt");
    h2->Draw("colz");
    CMS_lumi(c2);
    c2->Update();
    c2->RedrawAxis();
    c2->SaveAs("figures/dR_pt.png");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    canvas_setup(c3);
    TH1F *h3 = new TH1F("h3", ";#it{m(D_{s}^{*})} [GeV];# events", 100, 1.85, 2.35);
    mytree->Project("h3", "Dsst_mass");
    h3->SetLineColor(kBlack);
    h3->SetLineWidth(2);
    h3->Draw("hist");
    h3->SetMaximum(1.3*h3->GetMaximum());
    CMS_lumi(c3);
    c3->Update();
    c3->RedrawAxis();
    c3->SaveAs("figures/Dsst_mass.png");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    canvas_setup(c4);
    TH1F *h4 = new TH1F("h4", ";#it{m(D_{s}^{*}) - m(D_{s})} [GeV];# events", 100, 0, 0.3);
    mytree->Project("h4", "Dsst_mass - Ds_mass");
    h4->SetLineColor(kBlack);
    h4->SetLineWidth(2);
    h4->Draw("hist");
    h4->SetMaximum(1.3*h4->GetMaximum());
    CMS_lumi(c4);
    c4->Update();
    c4->RedrawAxis();
    c4->SaveAs("figures/mass_diff.png");

    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    canvas_setup(c5);
    TH1F *h5 = new TH1F("h5", ";#it{m(D_{s}^{*})} [GeV];# events", 100, 1.85, 2.35);
    mytree->Project("h5", "Dsst_mass", "gamma_deltaR < 0.02 && abs(Gen_gamma_pt - gamma_pt)/Gen_gamma_pt < 0.5");
    h5->SetLineColor(kBlack);
    h5->SetLineWidth(2);
    h5->Draw("hist");
    h5->SetMaximum(1.3*h5->GetMaximum());
    CMS_lumi(c5);
    c5->Update();
    c5->RedrawAxis();
    c5->SaveAs("figures/Dsst_mass_tight.png");

    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
    canvas_setup(c6);
    TH1F *h6 = new TH1F("h6", ";#it{m(D_{s}^{*}) - m(D_{s})} [GeV];# events", 100, 0, 0.3);
    mytree->Project("h6", "Dsst_mass - Ds_mass", "gamma_deltaR < 0.02 && abs(Gen_gamma_pt - gamma_pt)/Gen_gamma_pt < 0.5");
    h6->SetLineColor(kBlack);
    h6->SetLineWidth(2);
    h6->Draw("hist");
    h6->SetMaximum(1.3*h6->GetMaximum());
    CMS_lumi(c6);
    c6->Update();
    c6->RedrawAxis();
    c6->SaveAs("figures/mass_diff_tight.png");


    return 0;
}
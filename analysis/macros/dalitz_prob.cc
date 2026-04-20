#include <vector>
#include <iostream>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>

#include "../functions/tdrStyle.cc"
#include "../functions/CMS_lumi.cc"
#include "../functions/draw_funcs.cc"

Int_t dalitz_prob()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    TChain *mychain = new TChain("GenMatchDsWKKpimu/Events");

    mychain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/GenMatch_DsW_KKpimu/GenMatch_DsW_KKpimu.root");

    TH2F * h1 = new TH2F("h1", ";SV fit #it{#chi}^{2} probablity;Gen #it{M}(#it{K^{+}K^{-}}) [GeV]", 100, 0, 1, 100, 0.95, 2.2);
    TH2F * h2 = new TH2F("h2", ";SV fit #it{#chi}^{2} probablity;Gen #it{M}(#it{K^{+}#pi^{+}}) [GeV]", 100, 0, 1, 100, 0.6, 2.2);
    TH2F * h3 = new TH2F("h3", ";SV fit #it{#chi}^{2} probablity;Gen #it{M}(#it{K^{-}#pi^{+}}) [GeV]", 100, 0, 1, 100, 0.6, 2.2);

    TH2F * h4 = new TH2F("h4", ";SV fit #it{#chi}^{2}/ndf;Gen #it{M}(#it{K^{+}K^{-}}) [GeV]", 100, 0, 20, 100, 0.95, 2.2);
    TH2F * h5 = new TH2F("h5", ";SV fit #it{#chi}^{2}/ndf;Gen #it{M}(#it{K^{+}#pi^{+}}) [GeV]", 100, 0, 20, 100, 0.6, 2.2);
    TH2F * h6 = new TH2F("h6", ";SV fit #it{#chi}^{2}/ndf;Gen #it{M}(#it{K^{-}#pi^{+}}) [GeV]", 100, 0, 20, 100, 0.6, 2.2);

    TH1F * h7 = new TH1F("h7", ";Gen #it{M^{2}}(#it{K^{+}K^{-}}) [GeV]; # events", 100, 0.9, 3.5);
    TH1F * h8 = new TH1F("h8", ";Gen #it{M^{2}}(#it{K^{+}#pi^{+}}) [GeV]; # events", 100, 0.3, 2.4);
    TH1F * h9 = new TH1F("h9", ";Gen #it{M^{2}}(#it{K^{-}#pi^{+}}) [GeV]; # events", 100, 0.3, 2.4);

    mychain->Project("h1", "KK_mass:Ds_prob");
    mychain->Project("h2", "Kppi_mass:Ds_prob");
    mychain->Project("h3", "Kmpi_mass:Ds_prob");

    mychain->Project("h4", "KK_mass:Ds_chi2ndf");
    mychain->Project("h5", "Kppi_mass:Ds_chi2ndf");
    mychain->Project("h6", "Kmpi_mass:Ds_chi2ndf");

    mychain->Project("h7", "KK_mass*KK_mass");
    mychain->Project("h8", "Kppi_mass*Kppi_mass");
    mychain->Project("h9", "Kmpi_mass*Kmpi_mass");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    h1->Draw("colz");
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/KK_vs_prob.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvas_setup(c2);
    h2->Draw("colz");
    CMS_lumi(c2);
    c2->Update();
    c2->RedrawAxis();
    c2->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kppi_vs_prob.png");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    canvas_setup(c3);
    h3->Draw("colz");
    CMS_lumi(c3);
    c3->Update();
    c3->RedrawAxis();
    c3->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kmpi_vs_prob.png");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    canvas_setup(c4);
    h4->Draw("colz");
    CMS_lumi(c4);
    c4->Update();
    c4->RedrawAxis();
    c4->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/KK_vs_chi2ndf.png");


    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    canvas_setup(c5);
    h5->Draw("colz");
    CMS_lumi(c5);
    c5->Update();
    c5->RedrawAxis();
    c5->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kppi_vs_chi2ndf.png");


    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
    canvas_setup(c6);
    h6->Draw("colz");
    CMS_lumi(c6);
    c6->Update();
    c6->RedrawAxis();
    c6->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kmpi_vs_chi2ndf.png");

    TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
    canvas_setup(c7);
    h7->Draw("hist");
    CMS_lumi(c7);
    c7->Update();
    c7->RedrawAxis();
    c7->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/KK_mass2.png");

    TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
    canvas_setup(c8);
    h8->Draw("hist");
    CMS_lumi(c8);
    c8->Update();
    c8->RedrawAxis();
    c8->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kppi_mass2.png");

    TCanvas *c9 = new TCanvas("c9", "c9", 800, 600);
    canvas_setup(c9);
    h9->Draw("hist");
    CMS_lumi(c9);
    c9->Update();
    c9->RedrawAxis();
    c9->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsW_KKpimu/figures/Kmpi_mass2.png");
    
    return 0;
}

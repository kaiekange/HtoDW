#include "TFile.h"
#include "TTree.h"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/tdrStyle.cc"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/CMS_lumi.cc"
#include "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/functions/draw_funcs.cc"

void draw_compare(TTree *mytree, TString varname, Float_t varmin, Float_t varmax, TString htitle, TString figname)
{
    setTDRStyle();
    lumi_sqrtS = "13.6 TeV, 2022";

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    canvas->SetRightMargin(0.2);
    TLegend *HistLegend = new TLegend(0.8, 0.7, 0.93, 0.9, NULL, "brNDC");
    HistLegend->SetTextFont(42);
    HistLegend->SetTextColor(kBlack);
    HistLegend->SetBorderSize(0);
    HistLegend->SetFillColor(0);
    HistLegend->SetLineWidth(0);
    HistLegend->SetTextSize(0.04);
    
    TH1F *h1 = new TH1F("h1", ";"+htitle+";# normalized events", 100, varmin, varmax);
    TH1F *h2 = new TH1F("h2", "", 100, varmin, varmax);

    mytree->Project("h1", varname, "event_match");
    mytree->Project("h2", varname, "!event_match");

    h1->Scale(1. / h1->Integral());
    h2->Scale(1. / h2->Integral());

    HistLegend->AddEntry(h1, "Match", "l");
    HistLegend->AddEntry(h2, "Non-match", "l");

    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);

    h1->Draw("hist");
    h2->Draw("hist same");
    h1->SetMaximum(1.5*std::max(h1->GetMaximum(), h2->GetMaximum()));

    HistLegend->Draw();

    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("figures/" + figname);

}

Int_t compare()
{
    ROOT::EnableImplicitMT();
    TFile *myfile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/RecoAll_DsW_KKpimu/RecoAll_DsW_KKpimu.root");
    TTree *mytree = (TTree *)myfile->Get("RecoAllDsWKKpimu/Events");

    draw_compare(mytree, "Ds_FDxysig", 0, 10, "#it{D_{s}} FD_{#it{xy}} significance", "Ds_FDxysig.png");
    draw_compare(mytree, "Ds_IsoR04_RelIso", 0, 10, "#it{D_{s}} IsoR04 RelIso", "Ds_IsoR04_RelIso.png");
    draw_compare(mytree, "dR_Kp_Km", 0, 1, "#it{#DeltaR(K^{+},K^{-})}", "dR_Kp_Km.png");
    draw_compare(mytree, "dR_Kp_pi", 0, 1, "#it{#DeltaR(K^{+},#pi^{+})}", "dR_Kp_pi.png");
    draw_compare(mytree, "dR_Km_pi", 0, 1, "#it{#DeltaR(K^{-},#pi^{+})}", "dR_Km_pi.png");
    draw_compare(mytree, "dR_Kp_Ds", 0, 1, "#it{#DeltaR(K^{+},D_{s}^{+})}", "dR_Kp_Ds.png");
    draw_compare(mytree, "dR_Km_Ds", 0, 1, "#it{#DeltaR(K^{-},D_{s}^{+})}", "dR_Km_Ds.png");
    draw_compare(mytree, "dR_pi_Ds", 0, 1, "#it{#DeltaR(#pi^{+},D_{s}^{+})}", "dR_pi_Ds.png");

    return 0;
}
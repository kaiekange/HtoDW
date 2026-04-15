#include <vector>
#include <iostream>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>

#include "../functions/tdrStyle.cc"
#include "../functions/CMS_lumi.cc"
#include "../functions/draw_funcs.cc"

int compare_Dsst()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    TChain *chain1 = new TChain("GenRecoDsW/Events");
    TChain *chain2 = new TChain("GenRecoDsstW/Events");

    chain1->Add("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenRecoDsW/output_.root");
    chain2->Add("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenRecoDsstW/output_.root");

    TCanvas *canvas = new TCanvas("c1", "canvas", 800, 600);
    canvas_setup(canvas);
    canvas->SetRightMargin(0.13);
    TLegend *HistLegend = new TLegend(0.85, 0.5, 0.95, 0.9, NULL, "brNDC");
    HistLegend->SetTextFont(42);
    HistLegend->SetTextColor(kBlack);
    HistLegend->SetBorderSize(0);
    HistLegend->SetFillColor(0);
    HistLegend->SetLineWidth(0);
    HistLegend->SetTextSize(0.04);

    TH1F * h1 = new TH1F("h1", ";Gen #it{M}_{vis} [GeV];# events", 100, 0, 125);
    TH1F * h2 = new TH1F("h2", "", 100, 0, 125);

    chain1->Project("h1", "Gen_vis_mass");
    chain2->Project("h2", "Gen_vis_mass");

    h1->Scale(1. / h1->Integral());
    h2->Scale(1. / h2->Integral());

    HistLegend->AddEntry(h1, "#it{M(D_{s}l)}", "l");
    HistLegend->AddEntry(h2, "#it{M(D_{s}*l)}", "l");

    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);

    h1->Draw("hist");
    h2->Draw("hist same");

    HistLegend->Draw();

    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/vis_mass.png");

    return 0;
}

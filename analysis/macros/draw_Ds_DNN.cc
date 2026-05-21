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

int draw_Ds_DNN()
{
    lumi_sqrtS = "2022 (13.6 TeV)";
    setTDRStyle();

    TChain *chain1 = new TChain("DsWKKpimu/Events");

    chain1->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsW_KKpimu/SIG_DsW_KKpimu_scored.root");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    TH1F * h1 = new TH1F("h1", ";#it{M(D_{s}^{+})} [GeV];# events", 100, 1.8, 2.2);
    chain1->Project("h1", "Ds_mass", "!event_match && score > 0.8");
    h1->SetMaximum(1.3*h1->GetMaximum());
    h1->Draw("hist");
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsW_KKpimu/figures/Ds_mass.png");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    canvas_setup(c4);
    TH1F * h4 = new TH1F("h4", ";score;# events", 100, 0., 1.);
    chain1->Project("h4", "score", "!event_match");
    h4->SetMaximum(1.3*h4->GetMaximum());
    h4->Draw("hist");
    CMS_lumi(c4);
    c4->Update();
    c4->RedrawAxis();
    c4->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsW_KKpimu/figures/match_score.png");

    return 0;
}

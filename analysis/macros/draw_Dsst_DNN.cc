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

int draw_Dsst_DNN()
{
    lumi_sqrtS = "2022 (13.6 TeV)";
    setTDRStyle();

    TChain *chain1 = new TChain("DsstWKKpimu/Events");

    chain1->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsstW_KKpimu/SIG_DsstW_KKpimu_scored.root");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    TH1F * h1 = new TH1F("h1", ";#it{M(D_{s}^{+})} [GeV];# events", 100, 1.8, 2.2);
    chain1->Project("h1", "Ds_mass", "!event_match && score > 0.8");
    h1->SetMaximum(1.3*h1->GetMaximum());
    h1->Draw("hist");
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu/figures/Ds_mass.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvas_setup(c2);
    TH1F * h2 = new TH1F("h2", ";#it{M(D_{s}^{*+})} [GeV];# events", 100, 1.9, 2.3);
    chain1->Project("h2", "Dsst_mass", "!event_match && score > 0.8");
    h2->SetMaximum(1.3*h2->GetMaximum());
    h2->Draw("hist");
    CMS_lumi(c2);
    c2->Update();
    c2->RedrawAxis();
    c2->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu/figures/Dsst_mass.png");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    canvas_setup(c3);
    TH1F * h3 = new TH1F("h3", ";#it{M(D_{s}^{*+}) - M(D_{s}^{+})} [GeV];# events", 100, 0., 0.3);
    chain1->Project("h3", "Dsst_mass - Ds_mass", "!event_match && score > 0.8");
    h3->SetMaximum(1.3*h3->GetMaximum());
    h3->Draw("hist");
    CMS_lumi(c3);
    c3->Update();
    c3->RedrawAxis();
    c3->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu/figures/Ds_mass_diff.png");

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    canvas_setup(c4);
    TH1F * h4 = new TH1F("h4", ";score;# events", 100, 0., 1.);
    chain1->Project("h4", "score", "!event_match");
    h4->SetMaximum(1.3*h4->GetMaximum());
    h4->Draw("hist");
    CMS_lumi(c4);
    c4->Update();
    c4->RedrawAxis();
    c4->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu/figures/match_score.png");

    return 0;
}

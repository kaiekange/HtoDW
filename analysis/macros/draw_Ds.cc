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

int draw_Dsst()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    TChain *chain1 = new TChain("GenMatchDsstWKKpimu/Events");

    chain1->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/GenMatch_DsstW_KKpimu/GenMatch_DsstW_KKpimu.root");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    TH1F * h1 = new TH1F("h1", ";#it{M(D_{s}^{*+})} [GeV];# events", 100, 1.85, 3.0);
    chain1->Project("h1", "Dsst_mass");
    h1->SetMaximum(1.3*h1->GetMaximum());
    h1->Draw("hist");
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsstW_KKpimu/figures/Dsst_mass.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvas_setup(c2);
    TH1F * h2 = new TH1F("h2", ";#it{M(D_{s}^{*+}) - M(D_{s}^{+})} [GeV];# events", 100, 0., 1.);
    chain1->Project("h2", "Dsst_mass - Ds_mass");
    h2->SetMaximum(1.3*h2->GetMaximum());
    h2->Draw("hist");
    CMS_lumi(c2);
    c2->Update();
    c2->RedrawAxis();
    c2->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/GenMatch_DsstW_KKpimu/figures/Ds_mass_diff.png");

    return 0;
}

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

int draw_Ds()
{
    lumi_sqrtS = "2022 (13.6 TeV)";
    setTDRStyle();

    TChain *chain1 = new TChain("DsstWKKpimuonDsW/Events");

    chain1->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsstW_KKpimu_onDsW/SIG_DsstW_KKpimu_onDsW.root");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    TH1F * h1 = new TH1F("h1", ";#it{M(D_{s}^{*+})} [GeV];# events", 100, 1.9, 2.3);
    chain1->Project("h1", "Dsst_mass", "Ds_match && mu_match");
    h1->SetMaximum(1.3*h1->GetMaximum());
    h1->Draw("hist");
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu_onDsW/figures/Dsst_mass.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvas_setup(c2);
    TH1F * h2 = new TH1F("h2", ";#it{M(D_{s}^{*+}) - M(D_{s}^{+})} [GeV];# events", 100, 0., 0.3);
    chain1->Project("h2", "Dsst_mass - Ds_mass", "Ds_match && mu_match");
    h2->SetMaximum(1.3*h2->GetMaximum());
    h2->Draw("hist");
    CMS_lumi(c2);
    c2->Update();
    c2->RedrawAxis();
    c2->SaveAs("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Signal_DsstW_KKpimu_onDsW/figures/Ds_mass_diff.png");

    return 0;
}

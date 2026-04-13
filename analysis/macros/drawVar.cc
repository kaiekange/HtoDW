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

Int_t drawVar()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();
    TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/Signal/";

    TChain *mychain = new TChain("MatchRecoKineFit/Events");
    mychain->Add("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/MatchRecoKineFit/MatchRecoKineFit.root");

    // TChain *mychain = new TChain("Events");
    // mychain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_ReDefine.root");
    // TCut mycut = "event_match && num_ss_mu_tight == 0 && num_os_mu_tight == 1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";
    TCut mycut = "";

    // TCanvas *canvas_pi = new TCanvas("canvas_pi", "", 800, 600);
    // canvas_setup(canvas_pi);
    // TH1F *hist_pi = new TH1F("hist_pi", "", 100, 0, 5);
    // mychain->Project("hist_pi", "pi_pt", mycut);
    // hist_pi->GetXaxis()->SetTitle("#it{p_{T}(#pi^{+})} [GeV]");
    // hist_pi->GetYaxis()->SetTitle("# events");
    // hist_pi->SetMaximum(hist_pi->GetMaximum()*1.3);
    // hist_pi->Draw("ep");
    // CMS_lumi(canvas_pi);
    // canvas_pi->Update();
    // canvas_pi->RedrawAxis();
    // canvas_pi->SaveAs(figdir + "pi_pt.png");

    // TCanvas *canvas_phi = new TCanvas("canvas_phi", "", 800, 600);
    // canvas_setup(canvas_phi);
    // TH1F *hist_phi = new TH1F("hist_phi", "", 100, 1.019460-0.03, 1.019460+0.03);
    // mychain->Project("hist_phi", "phi_mass", mycut);
    // hist_phi->GetXaxis()->SetTitle("#it{M(K^{+}K^{-})} [GeV]");
    // hist_phi->GetYaxis()->SetTitle("# events");
    // hist_phi->SetMaximum(hist_phi->GetMaximum()*1.3);
    // hist_phi->Draw("ep");
    // CMS_lumi(canvas_phi);
    // canvas_phi->Update();
    // canvas_phi->RedrawAxis();
    // canvas_phi->SaveAs(figdir + "phi_mass.png");

    // TCanvas *canvas_Ds = new TCanvas("canvas_Ds", "", 800, 600);
    // canvas_setup(canvas_Ds);
    // TH1F *hist_Ds = new TH1F("hist_Ds", "", 100, 1.96835-0.15, 1.96835+0.15);
    // mychain->Project("hist_Ds", "Ds_mass", mycut);
    // hist_Ds->GetXaxis()->SetTitle("#it{M(K^{+}K^{-}#pi^{+})} [GeV]");
    // hist_Ds->GetYaxis()->SetTitle("# events");
    // hist_Ds->SetMaximum(hist_Ds->GetMaximum()*1.3);
    // hist_Ds->Draw("ep");
    // CMS_lumi(canvas_Ds);
    // canvas_Ds->Update();
    // canvas_Ds->RedrawAxis();
    // canvas_Ds->SaveAs(figdir + "Ds_mass.png");

    TCanvas *canvas_Ds_pt = new TCanvas("canvas_Ds_pt", "", 800, 600);
    canvas_setup(canvas_Ds_pt);
    TH1F *hist_Ds_pt = new TH1F("hist_Ds_pt", "", 100, 0, 80);
    mychain->Project("hist_Ds_pt", "Ds_pt", mycut);
    hist_Ds_pt->GetXaxis()->SetTitle("#it{p_{T}(D_{s}^{+})} [GeV]");
    hist_Ds_pt->GetYaxis()->SetTitle("# events");
    hist_Ds_pt->SetMaximum(hist_Ds_pt->GetMaximum()*1.3);
    hist_Ds_pt->Draw("ep");
    CMS_lumi(canvas_Ds_pt);
    canvas_Ds_pt->Update();
    canvas_Ds_pt->RedrawAxis();
    canvas_Ds_pt->SaveAs(figdir + "Ds_pt.png");


    return 0;
}

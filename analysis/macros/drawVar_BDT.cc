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

Int_t drawVar_BDT()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();
    TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/BDT/";

    TChain *mychain = new TChain("Events");
    mychain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/TT_TMVAapp.root");
    TCut mycut = "num_ss_mu==0 && num_os_mu==1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125 && BDT > 0.24";

    TCanvas *canvas_phi = new TCanvas("canvas_phi", "", 800, 600);
    canvas_setup(canvas_phi);
    TH1F *hist_phi = new TH1F("hist_phi", "", 30, 1.019460-0.03, 1.019460+0.03);
    mychain->Project("hist_phi", "phi_mass", mycut);
    hist_phi->GetXaxis()->SetTitle("#it{M(K^{+}K^{-})} [GeV]");
    hist_phi->GetYaxis()->SetTitle("# events");
    hist_phi->Draw("hist");
    CMS_lumi(canvas_phi);
    canvas_phi->Update();
    canvas_phi->RedrawAxis();
    canvas_phi->SaveAs(figdir + "TT_phi_mass.png");

    TCanvas *canvas_Ds = new TCanvas("canvas_Ds", "", 800, 600);
    canvas_setup(canvas_Ds);
    TH1F *hist_Ds = new TH1F("hist_Ds", "", 30, 1.96835-0.15, 1.96835+0.15);
    mychain->Project("hist_Ds", "Ds_mass", mycut);
    hist_Ds->GetXaxis()->SetTitle("#it{M(K^{+}K^{-}#pi^{+})} [GeV]");
    hist_Ds->GetYaxis()->SetTitle("# events");
    hist_Ds->Draw("hist");
    CMS_lumi(canvas_Ds);
    canvas_Ds->Update();
    canvas_Ds->RedrawAxis();
    canvas_Ds->SaveAs(figdir + "TT_Ds_mass.png");


    return 0;
}

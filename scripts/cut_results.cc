#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath, TString cutstr="1"){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", vartitle, 100, varmin, varmax);

    mychain->Project("hist", myvar, cutstr);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw();
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/cut_results/"+figpath+".png");
    delete hist;
}

void draw_2d(TChain *mychain, TString myvar, TString vartitle, float xvarmin, float xvarmax, float yvarmin, float yvarmax, TString figpath, TString cutstr="1"){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH2F *hist = new TH2F("hist", vartitle, 100, xvarmin, xvarmax, 100, yvarmin, yvarmax);

    mychain->Project("hist", myvar, cutstr);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    /* hist->SetMaximum(1.3*hist->GetMaximum()); */
    hist->SetMinimum(0);
    hist->Draw("COLZ");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/cut_results/"+figpath+".png");
    delete hist;
}

int cut_results() {
    setTDRStyle();

    TChain *mychain = new TChain("PreliminaryCut/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PreliminaryCut/PreliminaryCut2.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/cut_results")) gSystem->MakeDirectory("./figures/cut_results");

    TString cut_loose = "m_Ds > 1.85 && m_Ds < 2.1 && chi2_Ds > 0 && m_phi > 0.99 && m_phi < 1.05 && chi2_phi > 0";
    TString cut_tight = "m_Ds > 1.85 && m_Ds < 2.1 && chi2_Ds > 0 && m_phi > 0.99 && m_phi < 1.05 && chi2_phi > 0 && dR_Kp_Km < 0.05 && dR_Kp_pi < 0.2 && dR_Km_pi < 0.2 && dR_phi_pi < 0.2 && dR_phi_Ds < 0.1 && pt_Ds > 15 && pt_phi > 7";
   
    draw_1d(mychain, "p_Kp", ";#it{p(K^{+})} [GeV];# candidates", 1, 60, "p_Kp", cut_loose);
    draw_1d(mychain, "p_Km", ";#it{p(K^{-})} [GeV];# candidates", 1, 60, "p_Km", cut_loose);
    draw_1d(mychain, "p_pi", ";#it{p(#pi^{+})} [GeV];# candidates", 1, 80, "p_pi", cut_loose);

    draw_1d(mychain, "pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# candidates", 0.5, 60, "pt_Kp", cut_loose);
    draw_1d(mychain, "pt_Km", ";#it{p_{T}(K^{-})} [GeV];# candidates", 0.5, 60, "pt_Km", cut_loose);
    draw_1d(mychain, "pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# candidates", 0.5, 80, "pt_pi", cut_loose);

    draw_1d(mychain, "pt_Kp-pt_Km", ";#it{p_{T}(K^{+})-p_{T}(K^{-})} [GeV];# candidates", -20, 20, "pt_Kp_minus_pt_Km", cut_loose);
    draw_2d(mychain, "pt_Km:pt_Kp", ";#it{p_{T}(K^{+})} [GeV];#it{p_{T}(K^{-})} [GeV]", 0, 60, 0, 60, "pt_Kp_Km", cut_loose);

    draw_1d(mychain, "dR_Kp_Km", ";#Delta#it{R(K^{+},K^{-})};# candidates", 0, 0.15, "dR_Kp_Km", cut_loose);
    draw_1d(mychain, "dR_Kp_pi", ";#Delta#it{R(K^{+},#pi^{+})};# candidates", 0, 0.6, "dR_Kp_pi", cut_loose);
    draw_1d(mychain, "dR_Km_pi", ";#Delta#it{R(K^{-},#pi^{+})};# candidates", 0, 0.6, "dR_Km_pi", cut_loose);
    draw_1d(mychain, "d_Kp_Km", ";#Delta#it{d(K^{+},K^{-})} [cm];# candidates", 0, 0.2, "d_Kp_Km", cut_loose);
    draw_1d(mychain, "d_Kp_pi", ";#Delta#it{d(K^{+},#pi^{+})} [cm];# candidates", 0, 0.4, "d_Kp_pi", cut_loose);
    draw_1d(mychain, "d_Km_pi", ";#Delta#it{d(K^{-},#pi^{+})} [cm];# candidates", 0, 0.4, "d_Km_pi", cut_loose);

    draw_1d(mychain, "chi2_phi", ";#it{#chi^{2}(#phi)};# candidates", 0, 10, "chi2_phi", cut_loose);
    draw_1d(mychain, "pt_phi", ";#it{p_{T}(K^{+}K^{-})} [GeV];# candidates", 0, 70, "pt_phi", cut_loose);
    draw_1d(mychain, "p_phi", ";#it{p(K^{+}K^{-})} [GeV];# candidates", 0, 200, "p_phi", cut_loose);
    draw_1d(mychain, "m_phi", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "m_phi", cut_loose);
    draw_1d(mychain, "dR_phi_pi", ";#DeltaR(#it{#phi,#pi^{+})};# candidates", 0, 1, "dR_phi_pi", cut_loose);
    draw_1d(mychain, "d_phi_pi", ";#it{d(#phi,#pi^{+})} [cm];# candidates", 0, 6, "d_phi_pi", cut_loose);

    draw_1d(mychain, "chi2_Ds", ";#it{#chi^{2}(D_{s}^{+})};# candidates", 0, 25, "chi2_Ds", cut_loose);
    draw_1d(mychain, "pt_Ds", ";#it{p_{T}(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 0, 100, "pt_Ds", cut_loose);
    draw_1d(mychain, "p_Ds", ";#it{p(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 0, 300, "p_Ds", cut_loose);
    draw_1d(mychain, "m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "m_Ds", cut_loose);
    draw_1d(mychain, "dR_phi_Ds", ";#DeltaR(#it{#phi,D_{s}^{+})};# candidates", 0, 0.4, "dR_phi_Ds", cut_loose);
    draw_1d(mychain, "d_phi_Ds", ";#it{d(#phi,D_{s}^{+})} [cm];# candidates", 0, 4, "d_phi_Ds", cut_loose);

    draw_2d(mychain, "m_Ds:m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{M(K^{+}K^{-}#pi^{+})} [GeV]", 0.99, 1.05, 1.85, 2.1, "m_phi_Ds", cut_loose);
    draw_2d(mychain, "chi2_phi:m_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{#chi^{2}(#phi)}", 0.99, 1.05, 0, 10, "m_chi2_phi", cut_loose);
    draw_2d(mychain, "chi2_Ds:m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];#it{#chi^{2}(D_{s}^{+})}", 1.85, 2.1, 0, 25, "m_chi2_Ds", cut_loose);

    draw_1d(mychain, "m_phi", ";#it{M(K^{+}K^{-})} [GeV];# candidates", 0.99, 1.05, "m_phi_tight", cut_tight);
    draw_1d(mychain, "m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# candidates", 1.85, 2.1, "m_Ds_tight", cut_tight);
    
    
    delete mychain;

    return 0;
}

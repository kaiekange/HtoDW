#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"


int fit_results() {
    setTDRStyle();

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TChain *mychain1 = new TChain("PackedGenParticle/Events");
    mychain1->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked1/output_*.root");
    /* TChain *mychain2 = new TChain("PackedGenParticle/Events"); */
    /* mychain2->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked2/output_*.root"); */
    TChain *mychain3 = new TChain("PackedGenParticle/Events");
    mychain3->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/GenPacked3/output_*.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/fit_results")) gSystem->MakeDirectory("./figures/fit_results");
    
    TH1F *h_dR_Kp = new TH1F("h_dR_Kp", ";#Delta#it{R(K^{+})};# events", 100, 0, 10);
    TH1F *h_dR_Km = new TH1F("h_dR_Km", ";#Delta#it{R(K^{-})};# events", 100, 0, 10);
    TH1F *h_dR_pi = new TH1F("h_dR_pi", ";#Delta#it{R(#pi^{+})};# events", 100, 0, 10);
    TH1F *h1_dR_Kp = new TH1F("h1_dR_Kp", ";#Delta#it{R(K^{+})};# events", 100, 0, 0.1);
    TH1F *h1_dR_Km = new TH1F("h1_dR_Km", ";#Delta#it{R(K^{-})};# events", 100, 0, 0.1);
    TH1F *h1_dR_pi = new TH1F("h1_dR_pi", ";#Delta#it{R(#pi^{+})};# events", 100, 0, 0.1);

    TH1F *h_pt_Kp = new TH1F("h_pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 100, 0, 60);
    TH1F *h_pt_Km = new TH1F("h_pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 100, 0, 60);
    TH1F *h_pt_pi = new TH1F("h_pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 100, 0, 80);
    TH1F *h_dR_Kp_Km = new TH1F("h_dR_Kp_Km", ";#Delta#it{R(K^{+},K^{-})};# events", 100, 0, 0.2);
    TH1F *h_dR_Kp_pi = new TH1F("h_dR_Kp_pi", ";#Delta#it{R(K^{+},#pi^{+})};# events", 100, 0, 0.8);
    TH1F *h_dR_Km_pi = new TH1F("h_dR_Km_pi", ";#Delta#it{R(K^{-},#pi^{+})};# events", 100, 0, 0.8);
    TH1F *h1_pt_Kp = new TH1F("h1_pt_Kp", ";#it{p_{T}(K^{+})} [GeV];# events", 100, 0, 5);
    TH1F *h1_pt_Km = new TH1F("h1_pt_Km", ";#it{p_{T}(K^{-})} [GeV];# events", 100, 0, 5);
    TH1F *h1_pt_pi = new TH1F("h1_pt_pi", ";#it{p_{T}(#pi^{+})} [GeV];# events", 100, 0, 5);

    mychain1->Project("h_dR_Kp", "dR_Kp_vec");
    mychain1->Project("h_dR_Km", "dR_Km_vec");
    mychain1->Project("h_dR_pi", "dR_pi_vec");
    mychain1->Project("h1_dR_Kp", "dR_Kp_vec");
    mychain1->Project("h1_dR_Km", "dR_Km_vec");
    mychain1->Project("h1_dR_pi", "dR_pi_vec");

    mychain3->Project("h_pt_Kp", "pt_Kp");
    mychain3->Project("h_pt_Km", "pt_Km");
    mychain3->Project("h_pt_pi", "pt_pi");
    mychain3->Project("h_dR_Kp_Km", "dR_Kp_Km");
    mychain3->Project("h_dR_Kp_pi", "dR_Kp_pi");
    mychain3->Project("h_dR_Km_pi", "dR_Km_pi");
    mychain3->Project("h1_pt_Kp", "pt_Kp");
    mychain3->Project("h1_pt_Km", "pt_Km");
    mychain3->Project("h1_pt_pi", "pt_pi");
    
    TCanvas *c_dR_Kp = new TCanvas("c_dR_Kp", "c_dR_Kp", 800, 600);
    canvas_setup(c_dR_Kp);
    h_dR_Kp->SetMaximum(1.3*h_dR_Kp->GetMaximum());
    h_dR_Kp->SetMinimum(0);
    h_dR_Kp->Draw();
    CMS_lumi(c_dR_Kp);
	c_dR_Kp->Update();
	c_dR_Kp->RedrawAxis();
    c_dR_Kp->SaveAs("./figures/fit_results/dR_Kp.png");
    delete h_dR_Kp;

    TCanvas *c_dR_Km = new TCanvas("c_dR_Km", "c_dR_Km", 800, 600);
    canvas_setup(c_dR_Km);
    h_dR_Km->SetMaximum(1.3*h_dR_Km->GetMaximum());
    h_dR_Km->SetMinimum(0);
    h_dR_Km->Draw();
    CMS_lumi(c_dR_Km);
	c_dR_Km->Update();
	c_dR_Km->RedrawAxis();
    c_dR_Km->SaveAs("./figures/fit_results/dR_Km.png");
    delete h_dR_Km;

    TCanvas *c_dR_pi = new TCanvas("c_dR_pi", "c_dR_pi", 800, 600);
    canvas_setup(c_dR_pi);
    h_dR_pi->SetMaximum(1.3*h_dR_pi->GetMaximum());
    h_dR_pi->SetMinimum(0);
    h_dR_pi->Draw();
    CMS_lumi(c_dR_pi);
	c_dR_pi->Update();
	c_dR_pi->RedrawAxis();
    c_dR_pi->SaveAs("./figures/fit_results/dR_pi.png");
    delete h_dR_pi;

    TCanvas *c1_dR_Kp = new TCanvas("c1_dR_Kp", "c1_dR_Kp", 800, 600);
    canvas_setup(c1_dR_Kp);
    h1_dR_Kp->SetMaximum(1.3*h1_dR_Kp->GetMaximum());
    h1_dR_Kp->SetMinimum(0);
    h1_dR_Kp->Draw();
    CMS_lumi(c1_dR_Kp);
	c1_dR_Kp->Update();
	c1_dR_Kp->RedrawAxis();
    c1_dR_Kp->SaveAs("./figures/fit_results/dR_Kp_small.png");
    delete h1_dR_Kp;

    TCanvas *c1_dR_Km = new TCanvas("c1_dR_Km", "c1_dR_Km", 800, 600);
    canvas_setup(c1_dR_Km);
    h1_dR_Km->SetMaximum(1.3*h1_dR_Km->GetMaximum());
    h1_dR_Km->SetMinimum(0);
    h1_dR_Km->Draw();
    CMS_lumi(c1_dR_Km);
	c1_dR_Km->Update();
	c1_dR_Km->RedrawAxis();
    c1_dR_Km->SaveAs("./figures/fit_results/dR_Km_small.png");
    delete h1_dR_Km;

    TCanvas *c1_dR_pi = new TCanvas("c1_dR_pi", "c1_dR_pi", 800, 600);
    canvas_setup(c1_dR_pi);
    h1_dR_pi->SetMaximum(1.3*h1_dR_pi->GetMaximum());
    h1_dR_pi->SetMinimum(0);
    h1_dR_pi->Draw();
    CMS_lumi(c1_dR_pi);
	c1_dR_pi->Update();
	c1_dR_pi->RedrawAxis();
    c1_dR_pi->SaveAs("./figures/fit_results/dR_pi_small.png");
    delete h1_dR_pi;

//////
    TCanvas *c_pt_Kp = new TCanvas("c_pt_Kp", "c_pt_Kp", 800, 600);
    canvas_setup(c_pt_Kp);
    h_pt_Kp->SetMaximum(1.3*h_pt_Kp->GetMaximum());
    h_pt_Kp->SetMinimum(0);
    h_pt_Kp->Draw();
    CMS_lumi(c_pt_Kp);
	c_pt_Kp->Update();
	c_pt_Kp->RedrawAxis();
    c_pt_Kp->SaveAs("./figures/fit_results/pt_Kp.png");
    delete h_pt_Kp;

    TCanvas *c_pt_Km = new TCanvas("c_pt_Km", "c_pt_Km", 800, 600);
    canvas_setup(c_pt_Km);
    h_pt_Km->SetMaximum(1.3*h_pt_Km->GetMaximum());
    h_pt_Km->SetMinimum(0);
    h_pt_Km->Draw();
    CMS_lumi(c_pt_Km);
	c_pt_Km->Update();
	c_pt_Km->RedrawAxis();
    c_pt_Km->SaveAs("./figures/fit_results/pt_Km.png");
    delete h_pt_Km;

    TCanvas *c_pt_pi = new TCanvas("c_pt_pi", "c_pt_pi", 800, 600);
    canvas_setup(c_pt_pi);
    h_pt_pi->SetMaximum(1.3*h_pt_pi->GetMaximum());
    h_pt_pi->SetMinimum(0);
    h_pt_pi->Draw();
    CMS_lumi(c_pt_pi);
	c_pt_pi->Update();
	c_pt_pi->RedrawAxis();
    c_pt_pi->SaveAs("./figures/fit_results/pt_pi.png");
    delete h_pt_pi;

    TCanvas *c_dR_Kp_Km = new TCanvas("c_dR_Kp_Km", "c_dR_Kp_Km", 800, 600);
    canvas_setup(c_dR_Kp_Km);
    h_dR_Kp_Km->SetMaximum(1.3*h_dR_Kp_Km->GetMaximum());
    h_dR_Kp_Km->SetMinimum(0);
    h_dR_Kp_Km->Draw();
    CMS_lumi(c_dR_Kp_Km);
	c_dR_Kp_Km->Update();
	c_dR_Kp_Km->RedrawAxis();
    c_dR_Kp_Km->SaveAs("./figures/fit_results/dR_Kp_Km.png");
    delete h_dR_Kp_Km;

    TCanvas *c_dR_Kp_pi = new TCanvas("c_dR_Kp_pi", "c_dR_Kp_pi", 800, 600);
    canvas_setup(c_dR_Kp_pi);
    h_dR_Kp_pi->SetMaximum(1.3*h_dR_Kp_pi->GetMaximum());
    h_dR_Kp_pi->SetMinimum(0);
    h_dR_Kp_pi->Draw();
    CMS_lumi(c_dR_Kp_pi);
	c_dR_Kp_pi->Update();
	c_dR_Kp_pi->RedrawAxis();
    c_dR_Kp_pi->SaveAs("./figures/fit_results/dR_Kp_pi.png");
    delete h_dR_Kp_pi;

    TCanvas *c_dR_Km_pi = new TCanvas("c_dR_Km_pi", "c_dR_Km_pi", 800, 600);
    canvas_setup(c_dR_Km_pi);
    h_dR_Km_pi->SetMaximum(1.3*h_dR_Km_pi->GetMaximum());
    h_dR_Km_pi->SetMinimum(0);
    h_dR_Km_pi->Draw();
    CMS_lumi(c_dR_Km_pi);
	c_dR_Km_pi->Update();
	c_dR_Km_pi->RedrawAxis();
    c_dR_Km_pi->SaveAs("./figures/fit_results/dR_Km_pi.png");
    delete h_dR_Km_pi;
    
    TCanvas *c1_pt_Kp = new TCanvas("c1_pt_Kp", "c1_pt_Kp", 800, 600);
    canvas_setup(c1_pt_Kp);
    h1_pt_Kp->SetMaximum(1.3*h1_pt_Kp->GetMaximum());
    h1_pt_Kp->SetMinimum(0);
    h1_pt_Kp->Draw();
    CMS_lumi(c1_pt_Kp);
	c1_pt_Kp->Update();
	c1_pt_Kp->RedrawAxis();
    c1_pt_Kp->SaveAs("./figures/fit_results/pt_Kp_small.png");
    delete h1_pt_Kp;

    TCanvas *c1_pt_Km = new TCanvas("c1_pt_Km", "c1_pt_Km", 800, 600);
    canvas_setup(c1_pt_Km);
    h1_pt_Km->SetMaximum(1.3*h1_pt_Km->GetMaximum());
    h1_pt_Km->SetMinimum(0);
    h1_pt_Km->Draw();
    CMS_lumi(c1_pt_Km);
	c1_pt_Km->Update();
	c1_pt_Km->RedrawAxis();
    c1_pt_Km->SaveAs("./figures/fit_results/pt_Km_small.png");
    delete h1_pt_Km;

    TCanvas *c1_pt_pi = new TCanvas("c1_pt_pi", "c1_pt_pi", 800, 600);
    canvas_setup(c1_pt_pi);
    h1_pt_pi->SetMaximum(1.3*h1_pt_pi->GetMaximum());
    h1_pt_pi->SetMinimum(0);
    h1_pt_pi->Draw();
    CMS_lumi(c1_pt_pi);
	c1_pt_pi->Update();
	c1_pt_pi->RedrawAxis();
    c1_pt_pi->SaveAs("./figures/fit_results/pt_pi_small.png");
    delete h1_pt_pi;

    TH1F *h_m_phi = new TH1F("h_m_phi", ";#it{M(K^{+}K^{-})} [GeV];# events", 100, 0.99, 1.05);
    TH1F *h_m_Ds = new TH1F("h_m_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];# events", 100, 1.85, 2.1);
    TH2F *h_m_phi_Ds = new TH2F("h_m_phi_Ds", ";#it{M(K^{+}K^{-})} [GeV];#it{M(K^{+}K^{-}#pi^{+})} [GeV]", 100, 0.99, 1.05, 100, 1.85, 2.1);

    mychain3->Project("h_m_phi", "m_phi");
    mychain3->Project("h_m_Ds", "m_Ds");
    mychain3->Project("h_m_phi_Ds", "m_Ds:m_phi");

    TCanvas *c_m_phi = new TCanvas("c_m_phi", "c_m_phi", 800, 600);
    canvas_setup(c_m_phi);
    h_m_phi->SetMaximum(1.3*h_m_phi->GetMaximum());
    h_m_phi->SetMinimum(0);
    h_m_phi->Draw();
    CMS_lumi(c_m_phi);
	c_m_phi->Update();
	c_m_phi->RedrawAxis();
    c_m_phi->SaveAs("./figures/fit_results/m_phi.png");
    delete h_m_phi;

    TCanvas *c_m_Ds = new TCanvas("c_m_Ds", "c_m_Ds", 800, 600);
    canvas_setup(c_m_Ds);
    h_m_Ds->SetMaximum(1.3*h_m_Ds->GetMaximum());
    h_m_Ds->SetMinimum(0);
    h_m_Ds->Draw();
    CMS_lumi(c_m_Ds);
	c_m_Ds->Update();
	c_m_Ds->RedrawAxis();
    c_m_Ds->SaveAs("./figures/fit_results/m_Ds.png");
    delete h_m_Ds;
    
    TCanvas *c_m_phi_Ds = new TCanvas("c_m_phi_Ds", "c_m_phi_Ds", 800, 600);
    canvas_setup(c_m_phi_Ds);
    h_m_phi_Ds->SetMaximum(1.3*h_m_phi_Ds->GetMaximum());
    h_m_phi_Ds->SetMinimum(0);
    h_m_phi_Ds->Draw("COLZ");
    CMS_lumi(c_m_phi_Ds);
	c_m_phi_Ds->Update();
	c_m_phi_Ds->RedrawAxis();
    c_m_phi_Ds->SaveAs("./figures/fit_results/m_phi_Ds.png");
    delete h_m_phi_Ds;

    TH1F *h_chi2_phi = new TH1F("h_chi2_phi", ";#it{#chi^{2}(#phi)};# events", 100, 0, 10);
    TH1F *h_chi2_Ds = new TH1F("h_chi2_Ds", ";#it{#chi^{2}(D_{s}^{+})};# events", 100, 0, 25);
    TH2F *h_m_chi2_phi = new TH2F("h_m_chi2_phi", ";#it{M(K^{+}K^{-})} [GeV];#it{#chi^{2}(#phi)}", 100, 0.99, 1.05, 100, 0, 10);
    TH2F *h_m_chi2_Ds = new TH2F("h_m_chi2_Ds", ";#it{M(K^{+}K^{-}#pi^{+})} [GeV];#it{#chi^{2}(D_{s}^{+})}", 100, 1.85, 2.1, 100, 0, 25);
    mychain3->Project("h_chi2_phi", "chi2_phi");
    mychain3->Project("h_chi2_Ds", "chi2_Ds");
    mychain3->Project("h_m_chi2_phi", "chi2_phi:m_phi");
    mychain3->Project("h_m_chi2_Ds", "chi2_Ds:m_Ds");
    
    TCanvas *c_chi2_phi = new TCanvas("c_chi2_phi", "c_chi2_phi", 800, 600);
    canvas_setup(c_chi2_phi);
    h_chi2_phi->SetMaximum(1.3*h_chi2_phi->GetMaximum());
    h_chi2_phi->SetMinimum(0);
    h_chi2_phi->Draw();
    CMS_lumi(c_chi2_phi);
	c_chi2_phi->Update();
	c_chi2_phi->RedrawAxis();
    c_chi2_phi->SaveAs("./figures/fit_results/chi2_phi.png");
    delete h_chi2_phi;

    TCanvas *c_chi2_Ds = new TCanvas("c_chi2_Ds", "c_chi2_Ds", 800, 600);
    canvas_setup(c_chi2_Ds);
    h_chi2_Ds->SetMaximum(1.3*h_chi2_Ds->GetMaximum());
    h_chi2_Ds->SetMinimum(0);
    h_chi2_Ds->Draw();
    CMS_lumi(c_chi2_Ds);
	c_chi2_Ds->Update();
	c_chi2_Ds->RedrawAxis();
    c_chi2_Ds->SaveAs("./figures/fit_results/chi2_Ds.png");
    delete h_chi2_Ds;

    TCanvas *c_m_chi2_phi = new TCanvas("c_m_chi2_phi", "c_m_chi2_phi", 800, 600);
    canvas_setup(c_m_chi2_phi);
    h_m_chi2_phi->SetMaximum(1.3*h_m_chi2_phi->GetMaximum());
    h_m_chi2_phi->SetMinimum(0);
    h_m_chi2_phi->Draw("COLZ");
    CMS_lumi(c_m_chi2_phi);
	c_m_chi2_phi->Update();
	c_m_chi2_phi->RedrawAxis();
    c_m_chi2_phi->SaveAs("./figures/fit_results/m_chi2_phi.png");
    delete h_m_chi2_phi;

    TCanvas *c_m_chi2_Ds = new TCanvas("c_m_chi2_Ds", "c_m_chi2_Ds", 800, 600);
    canvas_setup(c_m_chi2_Ds);
    h_m_chi2_Ds->SetMaximum(1.3*h_m_chi2_Ds->GetMaximum());
    h_m_chi2_Ds->SetMinimum(0);
    h_m_chi2_Ds->Draw("COLZ");
    CMS_lumi(c_m_chi2_Ds);
	c_m_chi2_Ds->Update();
	c_m_chi2_Ds->RedrawAxis();
    c_m_chi2_Ds->SaveAs("./figures/fit_results/m_chi2_Ds.png");
    delete h_m_chi2_Ds;


    TH1F *h_d_Kp_Km = new TH1F("h_d_Kp_Km", ";#Delta#it{d(K^{+},K^{-})} [cm];# events", 100, 0, 0.5);
    TH1F *h_d_Kp_pi = new TH1F("h_d_Kp_pi", ";#Delta#it{d(K^{+},#pi^{+})} [cm];# events", 100, 0, 0.5);
    TH1F *h_d_Km_pi = new TH1F("h_d_Km_pi", ";#Delta#it{d(K^{-},#pi^{+})} [cm];# events", 100, 0, 0.5);
    TH1F *h_d_phi_pi = new TH1F("h_d_phi_pi", ";#Delta#it{d(#phi,#pi^{+})} [cm];# events", 100, 0, 5);
    TH1F *h_d_phi_Ds = new TH1F("h_d_phi_Ds", ";#Delta#it{d(#phi,D_{s}^{+})} [cm];# events", 100, 0, 5);
    mychain3->Project("h_d_Kp_Km", "d_Kp_Km");
    mychain3->Project("h_d_Kp_pi", "d_Kp_pi");
    mychain3->Project("h_d_Km_pi", "d_Km_pi");
    mychain3->Project("h_d_phi_pi", "d_phi_pi");
    mychain3->Project("h_d_phi_Ds", "d_phi_Ds");

    TCanvas *c_d_Kp_Km = new TCanvas("c_d_Kp_Km", "c_d_Kp_Km", 800, 600);
    canvas_setup(c_d_Kp_Km);
    h_d_Kp_Km->SetMaximum(1.3*h_d_Kp_Km->GetMaximum());
    h_d_Kp_Km->SetMinimum(0);
    h_d_Kp_Km->Draw();
    CMS_lumi(c_d_Kp_Km);
	c_d_Kp_Km->Update();
	c_d_Kp_Km->RedrawAxis();
    c_d_Kp_Km->SaveAs("./figures/fit_results/d_Kp_Km.png");
    delete h_d_Kp_Km;

    TCanvas *c_d_Kp_pi = new TCanvas("c_d_Kp_pi", "c_d_Kp_pi", 800, 600);
    canvas_setup(c_d_Kp_pi);
    h_d_Kp_pi->SetMaximum(1.3*h_d_Kp_pi->GetMaximum());
    h_d_Kp_pi->SetMinimum(0);
    h_d_Kp_pi->Draw();
    CMS_lumi(c_d_Kp_pi);
	c_d_Kp_pi->Update();
	c_d_Kp_pi->RedrawAxis();
    c_d_Kp_pi->SaveAs("./figures/fit_results/d_Kp_pi.png");
    delete h_d_Kp_pi;

    TCanvas *c_d_Km_pi = new TCanvas("c_d_Km_pi", "c_d_Km_pi", 800, 600);
    canvas_setup(c_d_Km_pi);
    h_d_Km_pi->SetMaximum(1.3*h_d_Km_pi->GetMaximum());
    h_d_Km_pi->SetMinimum(0);
    h_d_Km_pi->Draw();
    CMS_lumi(c_d_Km_pi);
	c_d_Km_pi->Update();
	c_d_Km_pi->RedrawAxis();
    c_d_Km_pi->SaveAs("./figures/fit_results/d_Km_pi.png");
    delete h_d_Km_pi;

    TCanvas *c_d_phi_pi = new TCanvas("c_d_phi_pi", "c_d_phi_pi", 800, 600);
    canvas_setup(c_d_phi_pi);
    h_d_phi_pi->SetMaximum(1.3*h_d_phi_pi->GetMaximum());
    h_d_phi_pi->SetMinimum(0);
    h_d_phi_pi->Draw();
    CMS_lumi(c_d_phi_pi);
	c_d_phi_pi->Update();
	c_d_phi_pi->RedrawAxis();
    c_d_phi_pi->SaveAs("./figures/fit_results/d_phi_pi.png");
    delete h_d_phi_pi;

    TCanvas *c_d_phi_Ds = new TCanvas("c_d_phi_Ds", "c_d_phi_Ds", 800, 600);
    canvas_setup(c_d_phi_Ds);
    h_d_phi_Ds->SetMaximum(1.3*h_d_phi_Ds->GetMaximum());
    h_d_phi_Ds->SetMinimum(0);
    h_d_phi_Ds->Draw();
    CMS_lumi(c_d_phi_Ds);
	c_d_phi_Ds->Update();
	c_d_phi_Ds->RedrawAxis();
    c_d_phi_Ds->SaveAs("./figures/fit_results/d_phi_Ds.png");
    delete h_d_phi_Ds;

    delete mychain1;
    delete mychain3;

    return 0;
}

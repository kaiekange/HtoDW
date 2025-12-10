#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath, TCut cutstr="1"){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", "", 100, varmin, varmax);

    mychain->Project("hist", myvar, cutstr);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->GetXaxis()->SetTitle(vartitle);
    hist->GetYaxis()->SetTitle("Normalized # candidates");
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw("ep");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
    delete hist;
}

void compare(TChain *mychain, TString myvar1, TString myvar2, TString leg1, TString leg2, TString vartitle, bool logorno, float varmin, float varmax, TCut cut1, TCut cut2, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h1 = new TH1F("h1", "", 100, varmin, varmax);
    mychain->Project("h1", myvar1, cut1);
    h1->Scale(1./h1->Integral());
    h1->SetLineColor(kBlack);
    h1->SetMarkerColor(kBlack);
    h1->SetMarkerSize(0.7);

    TH1F *h2 = new TH1F("h2", "", 100, varmin, varmax);
    mychain->Project("h2", myvar2, cut2);
    h2->Scale(1./h2->Integral());
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
    h2->SetMarkerSize(0.7);

    float height = std::max(h1->GetMaximum(), h2->GetMaximum());

    if (logorno) h1->SetMaximum(height*100);
    else {
        h1->SetMaximum(height*1.3);
        h1->SetMinimum(0);
    }
    h1->GetXaxis()->SetTitle(vartitle);
    h1->GetYaxis()->SetTitle("Normalized # candidates");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    if (logorno) canvas->SetLogy(1);
    canvas_setup(canvas);
    h1->Draw("ep");
    h2->Draw("ep same");

    TLegend * mylegend = new TLegend(0.7, 0.75, 0.9, 0.9);
    mylegend->AddEntry(h1, leg1, "l");
    mylegend->AddEntry(h2, leg2, "l");
    mylegend->SetFillColor(0);
    mylegend->SetLineWidth(0);
    mylegend->SetTextSize(0.04);
    mylegend->Draw();

    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
    delete h1;
    delete h2;
}

int match_compare() {
    setTDRStyle();

    TChain *mychain = new TChain("PVStudy/Events");
    mychain->Add("../../tuples/PVStudy.root");

    /* compare(mychain, "match_PVwithDs_vx - Gen_H.vx", "primvtx.vx - Gen_H.vx", "refitted PV", "MiniAOD PV", "Reco #it{x}_{PV} - Gen #it{x}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_x"); */ 
    /* compare(mychain, "match_PVwithDs_vy - Gen_H.vy", "primvtx.vy - Gen_H.vy", "refitted PV", "MiniAOD PV", "Reco #it{y}_{PV} - Gen #it{y}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_y"); */ 
    /* compare(mychain, "match_PVwithDs_vz - Gen_H.vz", "primvtx.vz - Gen_H.vz", "refitted PV", "MiniAOD PV", "Reco #it{z}_{PV} - Gen #it{z}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_z"); */ 

    /* draw_1d(mychain, "Gen_Kp.vx", "Gen #it{x}_{SV} [cm]", false, -1, 1, "Gen_SVx"); */
    /* draw_1d(mychain, "Gen_Kp.vy", "Gen #it{y}_{SV} [cm]", false, -1, 1, "Gen_SVy"); */
    /* draw_1d(mychain, "Gen_Kp.vz", "Gen #it{z}_{SV} [cm]", false, -15, 15, "Gen_SVz"); */

    /* draw_1d(mychain, "match_DsFit_vx - Gen_Kp.vx", "Residual #it{x}_{SV} [cm]", false, -0.2, 0.2, "res_SVx"); */
    /* draw_1d(mychain, "match_DsFit_vy - Gen_Kp.vy", "Residual #it{y}_{SV} [cm]", false, -0.2, 0.2, "res_SVy"); */
    /* draw_1d(mychain, "match_DsFit_vz - Gen_Kp.vz", "Residual #it{z}_{SV} [cm]", false, -0.2, 0.2, "res_SVz"); */
    
    /* draw_1d(mychain, "match_Kp_pt - Gen_Kp.pt", "Residual #it{p_{T}(K^{+})} [GeV]", false, -5, 5, "res_Kp_pt"); */
    /* draw_1d(mychain, "match_Km_pt - Gen_Km.pt", "Residual #it{p_{T}(K^{-})} [GeV]", false, -5, 5, "res_Km_pt"); */
    /* draw_1d(mychain, "match_pi_pt - Gen_pi.pt", "Residual #it{p_{T}(#pi^{+})} [GeV]", false, -5, 5, "res_pi_pt"); */
    /* draw_1d(mychain, "match_phiFit_phi_pt - Gen_phi.pt", "Residual #it{p_{T}(#phi)} [GeV]", false, -5, 5, "res_phi_pt"); */
    /* draw_1d(mychain, "match_DsFit_Ds_pt - Gen_Ds.pt", "Residual #it{p_{T}(D_{s}^{+})} [GeV]", false, -5, 5, "res_Ds_pt"); */
    
    /* draw_1d(mychain, "match_Kp_p - Gen_Kp.p", "Residual #it{p(K^{+})} [GeV]", false, -5, 5, "res_Kp_p"); */
    /* draw_1d(mychain, "match_Km_p - Gen_Km.p", "Residual #it{p(K^{-})} [GeV]", false, -5, 5, "res_Km_p"); */
    /* draw_1d(mychain, "match_pi_p - Gen_pi.p", "Residual #it{p(#pi^{+})} [GeV]", false, -5, 5, "res_pi_p"); */
    /* draw_1d(mychain, "match_phiFit_phi_p - Gen_phi.p", "Residual #it{p(#phi)} [GeV]", false, -5, 5, "res_phi_p"); */
    /* draw_1d(mychain, "match_DsFit_Ds_p - Gen_Ds.p", "Residual #it{p(D_{s}^{+})} [GeV]", false, -5, 5, "res_Ds_p"); */
    
    /* draw_1d(mychain, "match_Kp_eta - Gen_Kp.eta", "Residual #it{#eta(K^{+})}", false, -0.005, 0.005, "res_Kp_eta"); */
    /* draw_1d(mychain, "match_Km_eta - Gen_Km.eta", "Residual #it{#eta(K^{-})}", false, -0.005, 0.005, "res_Km_eta"); */
    /* draw_1d(mychain, "match_pi_eta - Gen_pi.eta", "Residual #it{#eta(#pi^{+})}", false, -0.005, 0.005, "res_pi_eta"); */
    /* draw_1d(mychain, "match_phiFit_phi_eta - Gen_phi.eta", "Residual #it{#eta(#phi)}", false, -0.005, 0.005, "res_phi_eta"); */
    /* draw_1d(mychain, "match_DsFit_Ds_eta - Gen_Ds.eta", "Residual #it{#eta(D_{s}^{+})}", false, -0.005, 0.005, "res_Ds_eta"); */
    
    /* draw_1d(mychain, "match_Kp_phi - Gen_Kp.phi", "Residual #it{#phi(K^{+})}", false, -0.005, 0.005, "res_Kp_phi"); */
    /* draw_1d(mychain, "match_Km_phi - Gen_Km.phi", "Residual #it{#phi(K^{-})}", false, -0.005, 0.005, "res_Km_phi"); */
    /* draw_1d(mychain, "match_pi_phi - Gen_pi.phi", "Residual #it{#phi(#pi^{+})}", false, -0.005, 0.005, "res_pi_phi"); */
    /* draw_1d(mychain, "match_phiFit_phi_phi - Gen_phi.phi", "Residual #it{#phi(#phi)}", false, -0.005, 0.005, "res_phi_phi"); */
    /* draw_1d(mychain, "match_DsFit_Ds_phi - Gen_Ds.phi", "Residual #it{#phi(D_{s}^{+})}", false, -0.005, 0.005, "res_Ds_phi"); */

    /* compare(mychain, "match_PVwithDs_vx - Gen_H.vx", "primvtx.vx - Gen_H.vx", "refitted PV", "MiniAOD PV", "Reco #it{x}_{PV} - Gen #it{x}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_x"); */ 
    /* compare(mychain, "match_PVwithDs_vy - Gen_H.vy", "primvtx.vy - Gen_H.vy", "refitted PV", "MiniAOD PV", "Reco #it{y}_{PV} - Gen #it{y}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_y"); */ 
    /* compare(mychain, "match_PVwithDs_vz - Gen_H.vz", "primvtx.vz - Gen_H.vz", "refitted PV", "MiniAOD PV", "Reco #it{z}_{PV} - Gen #it{z}_{PV} [cm]", false, -0.01, 0.01, "1", "1", "PV_z"); */ 

    /* compare(mychain, "match_Ds_PVwithDs_FDxy - sqrt(pow(Gen_Kp.vx-Gen_H.vx,2)+pow(Gen_Kp.vy-Gen_H.vy,2))", "match_Ds_primvtx_FDxy - sqrt(pow(Gen_Kp.vx-Gen_H.vx,2)+pow(Gen_Kp.vy-Gen_H.vy,2))", "refitted PV", "MiniAOD PV", "Reco FD_{#it{xy}}(#it{D_{s}^{+}}) - Gen FD_{#it{xy}}(#it{D_{s}^{+}})", false, -0.3, 0.5, "1", "1", "FDxy"); */ 
    /* compare(mychain, "match_Ds_PVwithDs_FDz - abs(Gen_Kp.vz-Gen_H.vz)", "match_Ds_primvtx_FDz - abs(Gen_Kp.vz-Gen_H.vz)", "refitted PV", "MiniAOD PV", "Reco FD_{#it{z}}(#it{D_{s}^{+}}) - Gen FD_{#it{z}}(#it{D_{s}^{+}})", false, -0.3, 0.5, "1", "1", "FDz"); */ 
    /* compare(mychain, "match_Ds_PVwithDs_FD - sqrt(pow(Gen_Kp.vx-Gen_H.vx,2)+pow(Gen_Kp.vy-Gen_H.vy,2)+pow(Gen_Kp.vz-Gen_H.vz,2))", "match_Ds_primvtx_FD - sqrt(pow(Gen_Kp.vx-Gen_H.vx,2)+pow(Gen_Kp.vy-Gen_H.vy,2)+pow(Gen_Kp.vz-Gen_H.vz,2))", "refitted PV", "MiniAOD PV", "Reco FD(#it{D_{s}^{+}}) - Gen FD(#it{D_{s}^{+}})", false, -0.3, 0.5, "1", "1", "FD"); */ 

    /* compare(mychain, "match_Ds_PVwithDs_FDxychi2", "match_Ds_primvtx_FDxychi2", "refitted PV", "MiniAOD PV", "FD_{#it{xy}}(#it{D_{s}^{+}}) significance", false, 0, 60, "1", "1", "FDxychi2"); */ 
    /* compare(mychain, "match_Ds_PVwithDs_FDzchi2", "match_Ds_primvtx_FDzchi2", "refitted PV", "MiniAOD PV", "FD_{#it{z}}(#it{D_{s}^{+}}) significance", false, 0, 60, "1", "1", "FDzchi2"); */ 
    /* compare(mychain, "match_Ds_PVwithDs_FDchi2", "match_Ds_primvtx_FDchi2", "refitted PV", "MiniAOD PV", "FD(#it{D_{s}^{+}}) significance", false, 0, 60, "1", "1", "FDchi2"); */ 

    compare(mychain, "match_Kp_PVwithDs_ip", "match_Kp_primvtx_ip", "refitted PV", "MiniAOD PV", "IP(#it{K^{+}}) [cm]", false, 0, 0.05, "1", "1", "Kp_ip"); 
    compare(mychain, "match_Km_PVwithDs_ip", "match_Km_primvtx_ip", "refitted PV", "MiniAOD PV", "IP(#it{K^{-}}) [cm]", false, 0, 0.05, "1", "1", "Km_ip"); 
    compare(mychain, "match_pi_PVwithDs_ip", "match_pi_primvtx_ip", "refitted PV", "MiniAOD PV", "IP(#it{#pi^{+}}) [cm]", false, 0, 0.1, "1", "1", "pi_ip"); 
    compare(mychain, "match_phi_PVwithDs_ip", "match_phi_primvtx_ip", "refitted PV", "MiniAOD PV", "IP(#it{#phi}) [cm]", false, 0, 0.05, "1", "1", "phi_ip"); 
    compare(mychain, "match_Ds_PVwithDs_ip", "match_Ds_primvtx_ip", "refitted PV", "MiniAOD PV", "IP(#it{D_{s}^{+}}) [cm]", false, 0, 0.05, "1", "1", "Ds_ip"); 

    compare(mychain, "match_Kp_PVwithDs_ipchi2", "match_Kp_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "IP(#it{K^{+}}) significance", false, 0, 20, "1", "1", "Kp_ipchi2"); 
    compare(mychain, "match_Km_PVwithDs_ipchi2", "match_Km_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "IP(#it{K^{-}}) significance", false, 0, 20, "1", "1", "Km_ipchi2"); 
    compare(mychain, "match_pi_PVwithDs_ipchi2", "match_pi_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "IP(#it{#pi^{+}}) significance", false, 0, 20, "1", "1", "pi_ipchi2"); 
    compare(mychain, "match_phi_PVwithDs_ipchi2", "match_phi_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "IP(#it{#phi}) significance", false, 0, 0.05, "1", "1", "phi_ipchi2"); 
    compare(mychain, "match_Ds_PVwithDs_ipchi2", "match_Ds_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "IP(#it{D_{s}^{+}}) significance", false, 0, 0.02, "1", "1", "Ds_ipchi2"); 

    compare(mychain, "match_Ds_PVwithDs_dira", "match_Ds_primvtx_dira", "refitted PV", "MiniAOD PV", "cos(#it{#theta})", false, -1, 1, "1", "1", "dira"); 
    compare(mychain, "match_Ds_PVwithDs_dira_angle", "match_Ds_primvtx_dira_angle", "refitted PV", "MiniAOD PV", "#it{D_{s}^{+}} direction angle #it{#theta}", false, 0, 3.1416, "1", "1", "dira_angle"); 

    compare(mychain, "match_Ds_PVwithDs_dira", "match_Ds_primvtx_dira", "refitted PV", "MiniAOD PV", "cos(#it{#theta})", true, -1, 1, "1", "1", "dira_log"); 
    compare(mychain, "match_Ds_PVwithDs_dira_angle", "match_Ds_primvtx_dira_angle", "refitted PV", "MiniAOD PV", "#it{D_{s}^{+}} direction angle #it{#theta}", true, 0, 3.1416, "1", "1", "dira_angle_log"); 


    delete mychain;

    return 0;
}

#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void draw_1d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, TString figpath, TCut cutstr="1"){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *hist = new TH1F("hist", vartitle, 100, varmin, varmax);

    mychain->Project("hist", myvar, cutstr);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    hist->SetMaximum(1.3*hist->GetMaximum());
    hist->SetMinimum(0);
    hist->Draw("ep");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
    delete hist;
}

void compare(TChain *mychain, TString myvar1, TString myvar2, TString leg1, TString leg2, TString vartitle, float varmin, float varmax, TCut cut1, TCut cut2, TString figpath){

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

    h1->SetMaximum(height*1.3);
    h1->SetMinimum(0);
    h1->GetXaxis()->SetTitle(vartitle);
    h1->GetYaxis()->SetTitle("Normalized # candidates");

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
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
    mychain->Add("../../tuples/PVStudy_noAP.root");

    TCut match_sel = "match_entry == 1";
    TCut nonmatch_sel = "non_match_entry == 1";
  
    /* compare(mychain, "match_PV_noDs_noBS_x - Gen_H_vx", "PV_vx - Gen_H_vx", "refitted PV", "MiniAOD PV", "Reco #it{x}_{PV} - Gen #it{x}_{PV} [cm]", -0.01, 0.01, "1", "1", "PV_x"); */ 
    /* compare(mychain, "match_PV_noDs_noBS_y - Gen_H_vy", "PV_vy - Gen_H_vy", "refitted PV", "MiniAOD PV", "Reco #it{y}_{PV} - Gen #it{y}_{PV} [cm]", -0.01, 0.01, "1", "1", "PV_y"); */ 
    /* compare(mychain, "match_PV_noDs_noBS_z - Gen_H_vz", "PV_vz - Gen_H_vz", "refitted PV", "MiniAOD PV", "Reco #it{z}_{PV} - Gen #it{z}_{PV} [cm]", -0.01, 0.01, "1", "1", "PV_z"); */ 
   
    /* compare(mychain, "match_Ds_PVnoDs_FDxy - sqrt(pow(Gen_Kp_vx-Gen_H_vx,2)+pow(Gen_Kp_vy-Gen_H_vy,2))", "match_Ds_primvtx_FDxy - sqrt(pow(Gen_Kp_vx-Gen_H_vx,2)+pow(Gen_Kp_vy-Gen_H_vy,2))", "refitted PV", "MiniAOD PV", "Reco FD_{#it{xy}}(#it{D_{s}^{+}}) - Gen FD_{#it{xy}}(#it{D_{s}^{+}})", -0.3, 0.5, "1", "1", "FDxy"); */ 
    /* compare(mychain, "match_Ds_PVnoDs_FDz - abs(Gen_Kp_vz-Gen_H_vz)", "match_Ds_primvtx_FDz - abs(Gen_Kp_vz-Gen_H_vz)", "refitted PV", "MiniAOD PV", "Reco FD_{#it{z}}(#it{D_{s}^{+}}) - Gen FD_{#it{z}}(#it{D_{s}^{+}})", -1, 2, "1", "1", "FDz"); */ 
    /* compare(mychain, "match_Ds_PVnoDs_FD - sqrt(pow(Gen_Kp_vx-Gen_H_vx,2)+pow(Gen_Kp_vy-Gen_H_vy,2)+pow(Gen_Kp_vz-Gen_H_vz,2))", "match_Ds_primvtx_FD - sqrt(pow(Gen_Kp_vx-Gen_H_vx,2)+pow(Gen_Kp_vy-Gen_H_vy,2)+pow(Gen_Kp_vz-Gen_H_vz,2))", "refitted PV", "MiniAOD PV", "Reco FD(#it{D_{s}^{+}}) - Gen FD(#it{D_{s}^{+}})", -1, 2, "1", "1", "FD"); */ 
    
    /* compare(mychain, "match_Ds_PVnoDs_FDxychi2", "match_Ds_primvtx_FDxychi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} FD_{#it{xy}}(#it{D_{s}^{+}})", 0, 60, "1", "1", "FDxychi2"); */ 
    /* compare(mychain, "match_Ds_PVnoDs_FDzchi2", "match_Ds_primvtx_FDzchi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} FD_{#it{z}}(#it{D_{s}^{+}})", 0, 60, "1", "1", "FDzchi2"); */ 
    /* compare(mychain, "match_Ds_PVnoDs_FDchi2", "match_Ds_primvtx_FDchi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} FD(#it{D_{s}^{+}})", 0, 60, "1", "1", "FDchi2"); */ 
    
    /* compare(mychain, "match_Kp_PVnoDs_ipchi2", "match_Kp_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} IP(#it{K^{+}})", 0, 30, "1", "1", "Kp_ipchi2"); */ 
    /* compare(mychain, "match_Km_PVnoDs_ipchi2", "match_Km_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} IP(#it{K^{-}})", 0, 30, "1", "1", "Km_ipchi2"); */ 
    /* compare(mychain, "match_pi_PVnoDs_ipchi2", "match_pi_primvtx_ipchi2", "refitted PV", "MiniAOD PV", "#it{#chi^{2}} IP(#it{#pi^{+}})", 0, 30, "1", "1", "pi_ipchi2"); */ 
    
    /* compare(mychain, "match_Ds_PVnoDs_dira", "match_Ds_primvtx_dira", "refitted PV", "MiniAOD PV", "cos(#it{#theta})", -1, 1, "1", "1", "dira"); */ 
    /* compare(mychain, "match_Ds_PVnoDs_dira_angle", "match_Ds_primvtx_dira_angle", "refitted PV", "MiniAOD PV", "#it{D_{s}^{+}} direction angle #it{#theta}", 0, 3.1416, "1", "1", "dira_angle"); */ 

    compare(mychain, "match_Ds_PVnoDs_dira", "match_Ds_primvtx_dira", "refitted PV", "MiniAOD PV", "cos(#it{#theta})", -1, 1, "1", "1", "dira"); 

    delete mychain;

    return 0;
}

#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "../CMSplots/tdrStyle.c"
#include "../CMSplots/CMS_lumi.c"
#include "../CMSplots/draw_funcs.c"

void draw_2d(TChain *mychain, TString myvar, TString vartitle, float varmin, float varmax, float yvarmin, float yvarmax, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH2F *hist = new TH2F("hist", vartitle, 100, varmin, varmax, 100, yvarmin, yvarmax);

    mychain->Project("hist", myvar);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    /* hist->SetMaximum(1.3*hist->GetMaximum()); */
    /* hist->SetMinimum(0); */
    hist->Draw("COLZ");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/match_compare/"+figpath+".png");
    delete hist;
}

int APvar() {
    setTDRStyle();

    TChain *mychain = new TChain("PVStudy/Events");
    mychain->Add("../../tuples/PVStudy_noAP.root");

    draw_2d(mychain, "APvar_phi:phiFit_phi_invm", ";#it{M}(#it{phi}) [GeV]; constant for #it{#phi}", 0.99, 1.05, 0, 0.04, "AP_phi"); 
    draw_2d(mychain, "APvar_Ds:DsFit_Ds_invm", ";#it{M}(#it{D_{s}^{+}}) [GeV]; constant for #it{D_{s}^{+}}", 1.85, 2.1, 0.4, 0.6, "AP_Ds"); 
    
    delete mychain;

    return 0;
}

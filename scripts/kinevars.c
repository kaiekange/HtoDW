#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void drawkinevar(TChain *mychain, TString var, TString vartitle, TString pdgid, float varmin, float varmax){
    TString fig_dir = "./figures";
    if(gSystem->AccessPathName(fig_dir)) gSystem->MakeDirectory(fig_dir);

    TH1F *h1 = new TH1F("h1", ";"+vartitle+";# events", 20, varmin, varmax);

    TString figpath1 = fig_dir + "/"+var+"_"+pdgid+".png";

    mychain->Project("h1", "match_part.m_state.p4Polar_.fCoordinates.f"+var, "match_part.m_state.pdgId_=="+pdgid);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    h1->SetMaximum(1.3*h1->GetMaximum());
    h1->SetMinimum(0);
    h1->Draw();
	c1->Update();
	c1->RedrawAxis();
    c1->SaveAs(figpath1);
    delete h1;
}

int kinevars() {
    setTDRStyle();

    TChain *mychain = new TChain("prunedGenParticles/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Higgs_charm/Simulation/20250202/2017/tuples/GenTM/output.root");

    TString partlist[] = {"H", "W^{-}", "D_{s}^{+}", "#phi", "#pi^{+}", "K^{+}", "K^{-}", "#mu^{-}", "#bar{#nu}_{#mu}"};
    TString partpdg[] = {"25", "-24", "431", "333", "211", "321", "-321", "13", "-14"};

    for(int i=0; i<9; i++){
        drawkinevar(mychain, "Pt", "#it{p_{T}("+partlist[i]+")} [GeV/#it{c}]", partpdg[i], 0, 150);
        drawkinevar(mychain, "Eta", "#it{#eta("+partlist[i]+")}", partpdg[i], -5, 5);
        drawkinevar(mychain, "Phi", "#it{#phi("+partlist[i]+")}", partpdg[i], -3.15, 3.15);
    }

    delete mychain;
    
    return 0;
}

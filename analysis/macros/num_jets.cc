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

enum ChannelId
{
    Signal = 0,
    WJets,
    // W0J,
    // W1J,
    // W2J,
    TT,
    TW,
    WW,
    WZ,
    NCOMB
};

const TString ChannelSample[NCOMB] = {
    "Signal",
    "WJets",
    // "W0J",
    // "W1J",
    // "W2J",
    "TT",
    "TW",
    "WW",
    "WZ"};

const TString ChannelLegend[NCOMB] = {
    "Signal",
    "W+jets",
    // "W+0j",
    // "W+1j",
    // "W+2j",
    "tt",
    "tW",
    "WW",
    "WZ"};

Int_t num_jets()
{

    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    const TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/";
    if (gSystem->AccessPathName(figdir))
        gSystem->MakeDirectory(figdir);

    TChain *ChannelChain[NCOMB] = {nullptr};

    for (Int_t cid = 0; cid < NCOMB; ++cid)
    {
        ChannelChain[cid] = new TChain("Events");
        ChannelChain[cid]->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[cid] + "_ReDefine.root");
    }

    TH1F *ChannelHist[NCOMB] = {nullptr};

    for (Int_t cid = 0; cid < NCOMB; ++cid)
    {
        TCut mycut;
        if(cid == 0)
            mycut = "event_match && num_ss_mu_tight==0 && num_os_mu_tight==1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";
        else
            mycut = "num_ss_mu==0 && num_os_mu==1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";

        TCanvas *canvas = new TCanvas(Form("canvas_%d", cid), "canvas", 800, 600);
        canvas_setup(canvas);
        ChannelHist[cid] = new TH1F(Form("hist_%s", ChannelSample[cid].Data()), "", 80, 0, 8);
        ChannelChain[cid]->Project(Form("hist_%s", ChannelSample[cid].Data()), "num_jets_pt30", mycut);
        ChannelHist[cid]->SetLineWidth(2);

        ChannelHist[cid]->GetXaxis()->SetTitle("# jets (#it{p_{T}} > 30 GeV)");
        ChannelHist[cid]->GetYaxis()->SetTitle("# events");
        ChannelHist[cid]->Draw("hist");

        CMS_lumi(canvas);
        canvas->Update();
        canvas->RedrawAxis();
        canvas->SaveAs(figdir + "num_jet_pt30_" + ChannelSample[cid] + ".png");
    }

    return 0;
}

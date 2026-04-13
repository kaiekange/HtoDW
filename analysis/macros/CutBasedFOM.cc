#include <vector>
#include <iostream>
#include <TChain.h>
#include <TString.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "../functions/tdrStyle.cc"
#include "../functions/CMS_lumi.cc"
#include "../functions/draw_funcs.cc"

enum CutId
{
    Ds_pt = 0,
    dR_Kp_Km,
    dR_Kp_pi,
    dR_Kp_phi,
    dR_Kp_Ds,
    dR_pi_phi,
    dR_pi_Ds,
    dR_phi_Ds,
    NCUT
};

const TString cutbasestr[NCUT] = {
    "Ds_pt > ",
    "dR_Kp_Km < ",
    "dR_Kp_pi < ",
    "dR_Kp_phi < ",
    "dR_Kp_Ds < ",
    "dR_pi_phi < ",
    "dR_pi_Ds < ",
    "dR_phi_Ds < "};

const TString grXtitle[NCUT] = {
    "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]",
    "#it{#DeltaR}(#it{K^{+}, #it{K^{-}}})",
    "#it{#DeltaR}(#it{K^{+}, #it{#pi^{+}}})",
    "#it{#DeltaR}(#it{K^{+}, #it{#phi}})",
    "#it{#DeltaR}(#it{K^{+}, #it{D_{s}^{+}}})",
    "#it{#DeltaR}(#it{#pi^{+}, #it{#phi}})",
    "#it{#DeltaR}(#it{#pi^{+}, #it{D_{s}^{+}}})",
    "#it{#DeltaR}(#it{#phi, #it{D_{s}^{+}}})"};

const TString vartitle[NCUT] = {
    "Ds_pt",
    "dR_Kp_Km",
    "dR_Kp_pi",
    "dR_Kp_phi",
    "dR_Kp_Ds",
    "dR_pi_phi",
    "dR_pi_Ds",
    "dR_phi_Ds"};

const Float_t varmin[NCUT] = {
    7,
    0,
    0,
    0,
    0,
    0,
    0,
    0};

const Float_t varmax[NCUT] = {
    30,
    0.2,
    0.8,
    0.1,
    0.3,
    0.8,
    0.8,
    0.2};

const Float_t stepsize[NCUT] = {
    0.1,
    0.001,
    0.01,
    0.001,
    0.001,
    0.005,
    0.005,
    0.001};

Int_t CutBasedFOM(Int_t idx)
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    const TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/CutBasedFOM/";
    if (gSystem->AccessPathName(figdir))
        gSystem->MakeDirectory(figdir);

    TChain *SigChain = new TChain("Events");
    TChain *BkgChain = new TChain("Events");

    SigChain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/tuples/Signal_ReDefine.root");
    BkgChain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/tuples/WJets_ReDefine.root");

    TCut sigbasecut = "event_match && num_ss_mu_tight == 0 && num_os_mu_tight == 1 && Ds_pt > 15 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";
    TCut bkgbasecut = "num_ss_mu == 0 && num_os_mu == 1 && Ds_pt > 15 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";

    Int_t nSignal = SigChain->GetEntries(sigbasecut);
    Int_t nBackground = BkgChain->GetEntries(bkgbasecut);

    Int_t nsteps = Int_t((varmax[idx] - varmin[idx]) / stepsize[idx]);

    Float_t cutvalue[nsteps];
    Float_t punzifom[nsteps];

    Float_t maxpunzi = -1;
    Float_t bestcutvalue;
    Float_t bestsigeff;
    Float_t bestbkgeff;

    for (Int_t i = 0; i < nsteps; i++)
    {
        cutvalue[i] = varmin[idx] + (i + 1) * stepsize[idx];
        TString varcutstr = cutbasestr[idx] + Form("%f", cutvalue[i]);
        Int_t nSig = SigChain->GetEntries(sigbasecut + TCut(varcutstr));
        Int_t nBkg = BkgChain->GetEntries(bkgbasecut + TCut(varcutstr));

        Float_t sigeff = Float_t(nSig) / Float_t(nSignal);

        punzifom[i] = sigeff / (1.5 + std::sqrt(nBkg));

        if (punzifom[i] > maxpunzi)
        {
            maxpunzi = punzifom[i];
            bestcutvalue = cutvalue[i];
            bestsigeff = sigeff;
            bestbkgeff = Float_t(nBkg) / Float_t(nBackground);
        }
    }

    std::cout << "max punzi fom: " << maxpunzi << ", at " + cutbasestr[idx] << bestcutvalue << ", signal effiicency: " << bestsigeff << ", background effiicency: " << bestbkgeff << std::endl;

    TGraph *gr = new TGraph(nsteps, cutvalue, punzifom);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    canvas->SetLeftMargin(0.2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle(grXtitle[idx]);
    gr->GetYaxis()->SetTitle("Punzi figure of merit");
    gr->GetYaxis()->SetNdivisions(810);
    gr->GetYaxis()->SetTitleOffset(1.6);
    CMS_lumi(canvas);
    canvas->Update();
    canvas->SaveAs(figdir + vartitle[idx] + ".png");

    return 0;
}
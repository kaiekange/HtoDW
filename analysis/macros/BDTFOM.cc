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

Int_t BDTFOM()
{
    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    const TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/BDT/";
    if (gSystem->AccessPathName(figdir))
        gSystem->MakeDirectory(figdir);

    TChain *SigChain = new TChain("Events");
    TChain *BkgChain = new TChain("Events");

    SigChain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_TMVAapp.root");
    BkgChain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/WJets_TMVAapp.root");

    Int_t nSignal = SigChain->GetEntries();
    Int_t nBackground = BkgChain->GetEntries();

    Float_t BDTmin = 0.;
    Float_t BDTmax = 0.29;
    Float_t stepsize = 0.01;

    const Int_t nsteps = Int_t((BDTmax - BDTmin) / stepsize) + 1;

    Float_t BDTvalue[nsteps];
    Float_t punzifom[nsteps];

    Float_t maxpunzi = -1;
    Float_t bestBDTvalue;
    Float_t bestsigeff;
    Float_t bestbkgeff;

    for (Int_t i = 0; i < nsteps; i++)
    {
        BDTvalue[i] = BDTmin + i * stepsize;
        TString BDTcutstr = Form("BDT > %f", BDTvalue[i]);
        Int_t nSig = SigChain->GetEntries(TCut(BDTcutstr));
        Int_t nBkg = BkgChain->GetEntries(TCut(BDTcutstr));
        Float_t nBkg_equiv = Float_t(nBkg) * 4 / 76566809. * 67710.0 * 1000 * 300;

        Float_t sigeff = Float_t(nSig) / Float_t(nSignal);

        punzifom[i] = sigeff / (1.5 + std::sqrt(nBkg_equiv));

        std::cout << "BDT > " << BDTvalue[i] << ", nsig " << nSig << ", nbkg " << nBkg << ", nbkg lumi-scaled " << nBkg_equiv << ", F.o.M: " << punzifom[i] << std::endl;

        if (punzifom[i] > maxpunzi)
        {
            maxpunzi = punzifom[i];
            bestBDTvalue = BDTvalue[i];
            bestsigeff = sigeff;
            bestbkgeff = Float_t(nBkg) / Float_t(nBackground);
        }
    }

    std::cout << "max punzi fom: " << maxpunzi << ", at BDT > " << bestBDTvalue << ", signal effiicency: " << bestsigeff << ", background effiicency: " << bestbkgeff << std::endl;

    TGraph *gr = new TGraph(nsteps, BDTvalue, punzifom);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    canvas->SetLeftMargin(0.2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("BDT cut");
    gr->GetYaxis()->SetTitle("Punzi figure of merit");
    gr->GetYaxis()->SetNdivisions(810);
    gr->GetYaxis()->SetTitleOffset(1.6);
    CMS_lumi(canvas);
    canvas->Update();
    canvas->SaveAs(figdir + "FOM.png");

    return 0;
}
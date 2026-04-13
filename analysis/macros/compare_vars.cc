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
    W0J,
    W1J,
    W2J,
    TT,
    TW,
    WW,
    WZ,
    NCOMB
};

const TString ChannelSample[NCOMB] = {
    "JMEcorr",
    "WJets",
    "W0J",
    "W1J",
    "W2J",
    "TT",
    "TW",
    "WW",
    "WZ"};

const TString ChannelLegend[NCOMB] = {
    "Signal",
    "W+jets",
    "W+0j",
    "W+1j",
    "W+2j",
    "tt",
    "tW",
    "WW",
    "WZ"};

const Color_t LineColors[NCOMB] = {
    kBlack,       // Signal
    kRed + 1,     // WJets
    kAzure + 1,   // W0J
    kBlue + 1,    // W1J
    kCyan + 1,    // W2J
    kMagenta + 1, // TT
    kOrange + 7,  // TW
    kGreen + 2,   // WW
    kViolet + 1   // WZ
};

enum HistId
{
    Ds_IsoR04_RelIso = 0,
    Ds_IsoR04_PFIso,
    Ds_IsoR03_RelIso,
    Ds_IsoR03_PFIso,

    mu_IsoR04_RelIso,
    mu_IsoR04_PFIso,
    mu_IsoR03_RelIso,
    mu_IsoR03_PFIso,

    phi_mass,
    Ds_mass,

    dira,
    Ds_chi2ndf,
    Ds_prob,

    met_pt,
    met_sumEt,

    num_jets_pt30,

    Kp_pt,
    Km_pt,
    pi_pt,
    phi_pt,
    Ds_pt,
    mu_pt,

    Kp_ip,
    Km_ip,
    pi_ip,
    Ds_ip,
    mu_ip,

    Kp_ipsig,
    Km_ipsig,
    pi_ipsig,
    Ds_ipsig,
    mu_ipsig,

    Ds_FDxy,
    Ds_FD,

    Ds_FDxysig,
    Ds_FDsig,

    dR_Kp_Km,
    dR_Kp_pi,
    dR_Kp_phi,
    dR_Kp_Ds,
    dR_Kp_mu,
    dR_Km_pi,
    dR_Km_phi,
    dR_Km_Ds,
    dR_Km_mu,
    dR_pi_phi,
    dR_pi_Ds,
    dR_pi_mu,
    dR_phi_Ds,
    dR_phi_mu,
    dR_Ds_mu,
    NHIST
};

const TString HistVar[NHIST] = {
    "Ds_IsoR04_RelIso",
    "Ds_IsoR04_PFIso",
    "Ds_IsoR03_RelIso",
    "Ds_IsoR03_PFIso",

    "mu_IsoR04_RelIso",
    "mu_IsoR04_PFIso",
    "mu_IsoR03_RelIso",
    "mu_IsoR03_PFIso",

    "phi_mass",
    "Ds_mass",

    "dira",
    "Ds_chi2ndf",
    "Ds_prob",

    "met_pt",
    "met_sumEt",

    "num_jets_pt30",

    "Kp_pt",
    "Km_pt",
    "pi_pt",
    "phi_pt",
    "Ds_pt",
    "mu_pt",

    "Kp_ip",
    "Km_ip",
    "pi_ip",
    "Ds_ip",
    "mu_ip",

    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "Ds_ipsig",
    "mu_ipsig",

    "Ds_FDxy",
    "Ds_FD",

    "Ds_FDxysig",
    "Ds_FDsig",

    "dR_Kp_Km",
    "dR_Kp_pi",
    "dR_Kp_phi",
    "dR_Kp_Ds",
    "dR_Kp_mu",
    "dR_Km_pi",
    "dR_Km_phi",
    "dR_Km_Ds",
    "dR_Km_mu",
    "dR_pi_phi",
    "dR_pi_Ds",
    "dR_pi_mu",
    "dR_phi_Ds",
    "dR_phi_mu",
    "dR_Ds_mu"};

const Int_t nbins[NHIST] = {
    100,  // Ds_IsoR04_RelIso
    100,  // Ds_IsoR04_PFIso
    100,  // Ds_IsoR03_RelIso
    100,  // Ds_IsoR03_PFIso

    100,  // mu_IsoR04_RelIso
    100,  // mu_IsoR04_PFIso
    100,  // mu_IsoR03_RelIso
    100,  // mu_IsoR03_PFIso

    100,  // phi_mass
    100,  // Ds_mass

    100,  // dira
    100,  // Ds_chi2ndf
    100,  // Ds_prob

    100,  // met_pt
    100,  // met_sumEt

    8,    // num_jets_pt30

    100,  // Kp_pt
    100,  // Km_pt
    100,  // pi_pt
    100,  // phi_pt
    100,  // Ds_pt
    100,  // mu_pt

    100,  // Kp_ip
    100,  // Km_ip
    100,  // pi_ip
    100,  // Ds_ip
    100,  // mu_ip

    100,  // Kp_ipsig
    100,  // Km_ipsig
    100,  // pi_ipsig
    100,  // Ds_ipsig
    100,  // mu_ipsig

    100,  // Ds_FDxy
    100,  // Ds_FD

    100,  // Ds_FDxysig
    100,  // Ds_FDsig

    100,  // dR_Kp_Km
    100,  // dR_Kp_pi
    100,  // dR_Kp_phi
    100,  // dR_Kp_Ds
    100,  // dR_Kp_mu
    100,  // dR_Km_pi
    100,  // dR_Km_phi
    100,  // dR_Km_Ds
    100,  // dR_Km_mu
    100,  // dR_pi_phi
    100,  // dR_pi_Ds
    100,  // dR_pi_mu
    100,  // dR_phi_Ds
    100,  // dR_phi_mu
    100}; // dR_Ds_mu

const Float_t varmin[NHIST] = {
    0,      // Ds_IsoR04_RelIso
    0,      // Ds_IsoR04_PFIso
    0,      // Ds_IsoR03_RelIso
    0,      // Ds_IsoR03_PFIso

    0,      // mu_IsoR04_RelIso
    0,      // mu_IsoR04_PFIso
    0,      // mu_IsoR03_RelIso
    0,      // mu_IsoR03_PFIso

    0.99,   // phi_mass
    1.82,   // Ds_mass

    0.9,    // dira
    0,      // Ds_chi2ndf
    0.01,   // Ds_prob

    0,      // met_pt
    0,      // met_sumEt

    0,      // num_jets_pt30

    1,      // Kp_pt
    1,      // Km_pt
    0.5,      // pi_pt
    3,      // phi_pt
    7,      // Ds_pt
    24,      // mu_pt

    0,      // Kp_ip
    0,      // Km_ip
    0,      // pi_ip
    0,      // Ds_ip
    0,      // mu_ip

    0,      // Kp_ipsig
    0,      // Km_ipsig
    0,      // pi_ipsig
    0,      // Ds_ipsig
    0,      // mu_ipsig

    0,      // Ds_FDxy
    0,      // Ds_FD

    0,      // Ds_FDxysig
    0,      // Ds_FDsig

    0,      // dR_Kp_Km
    0,      // dR_Kp_pi
    0,      // dR_Kp_phi
    0,      // dR_Kp_Ds
    0,      // dR_Kp_mu
    0,      // dR_Km_pi
    0,      // dR_Km_phi
    0,      // dR_Km_Ds
    0,      // dR_Km_mu
    0,      // dR_pi_phi
    0,      // dR_pi_Ds
    0,      // dR_pi_mu
    0,      // dR_phi_Ds
    0,      // dR_phi_mu
    0};     // dR_Ds_mu

const Float_t varmax[NHIST] = {
    0.5,    // Ds_IsoR04_RelIso
    30,     // Ds_IsoR04_PFIso
    0.5,    // Ds_IsoR03_RelIso
    30,     // Ds_IsoR03_PFIso

    0.2,    // mu_IsoR04_RelIso
    8,      // mu_IsoR04_PFIso
    0.2,    // mu_IsoR03_RelIso
    8,      // mu_IsoR03_PFIso

    1.048,  // phi_mass
    2.1,    // Ds_mass

    1,      // dira
    3.7,    // Ds_chi2ndf
    1,      // Ds_prob

    200,    // met_pt
    1000,   // met_sumEt

    8,      // num_jets_pt30

    40,     // Kp_pt
    40,     // Km_pt
    40,     // pi_pt
    60,    // phi_pt
    80,    // Ds_pt
    200,    // mu_pt

    0.05,    // Kp_ip
    0.05,    // Km_ip
    0.1,    // pi_ip
    0.02,   // Ds_ip
    0.01,   // mu_ip

    20,     // Kp_ipsig
    20,     // Km_ipsig
    30,     // pi_ipsig
    10,     // Ds_ipsig
    8,      // mu_ipsig

    1,      // Ds_FDxy
    2,      // Ds_FD

    40,     // Ds_FDxysig
    40,     // Ds_FDsig

    0.2,    // dR_Kp_Km
    0.8,    // dR_Kp_pi
    0.1,    // dR_Kp_phi
    0.3,    // dR_Kp_Ds
    5,      // dR_Kp_mu
    0.8,    // dR_Km_pi
    0.1,    // dR_Km_phi
    0.3,    // dR_Km_Ds
    5,      // dR_Km_mu
    1,      // dR_pi_phi
    1,      // dR_pi_Ds
    5,      // dR_pi_mu
    0.3,    // dR_phi_Ds
    5,      // dR_phi_mu
    5};     // dR_Ds_mu

const TString HistXTitle[NHIST] = {
    "#it{D_{s}^{+}} IsoR04 RelIso",  // Ds_IsoR04_RelIso
    "#it{D_{s}^{+}} IsoR04 PFIso",   // Ds_IsoR04_PFIso
    "#it{D_{s}^{+}} IsoR03 RelIso",  // Ds_IsoR03_RelIso
    "#it{D_{s}^{+}} IsoR03 PFIso",   // Ds_IsoR03_PFIso

    "#it{#mu^{-}} IsoR04 RelIso",    // mu_IsoR04_RelIso
    "#it{#mu^{-}} IsoR04 PFIso",     // mu_IsoR04_PFIso
    "#it{#mu^{-}} IsoR03 RelIso",    // mu_IsoR03_RelIso
    "#it{#mu^{-}} IsoR03 PFIso",     // mu_IsoR03_PFIso

    "#it{M}(#it{K^{+}K^{-}}) [GeV]",      // phi_mass
    "#it{M}(#it{K^{+}K^{-}#pi^{+}}) [GeV]", // Ds_mass

    "cos #alpha(PV, SV)",                // dira
    "#it{D_{s}^{+}} fit #chi^{2}/ndf",   // Ds_chi2ndf
    "#it{D_{s}^{+}} fit probability",    // Ds_prob

    "MET #it{p_T} [GeV]",           // met_pt
    "MET #it{#sum E_{T}} [GeV]",    // met_sumEt

    "# jets (#it{p_{T}} > 30 GeV)", // num_jets_pt30

    "#it{p_{T}}(#it{K^{+}}) [GeV]",  // Kp_pt
    "#it{p_{T}}(#it{K^{-}}) [GeV]",  // Km_pt
    "#it{p_{T}}(#it{#pi^{+}}) [GeV]", // pi_pt
    "#it{p_{T}}(#it{#phi}) [GeV]",   // phi_pt
    "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", // Ds_pt
    "#it{p_{T}}(#it{#mu^{-}}) [GeV]", // mu_pt

    "IP(#it{K^{+}})",    // Kp_ip
    "IP(#it{K^{-}})",    // Km_ip
    "IP(#it{#pi^{+}})",  // pi_ip
    "IP(#it{D_{s}^{+}})", // Ds_ip
    "IP(#it{#mu^{-}})",  // mu_ip

    "IP significance(#it{K^{+}})",  // Kp_ipsig
    "IP significance(#it{K^{-}})",  // Km_ipsig
    "IP significance(#it{#pi^{+}})", // pi_ipsig
    "IP significance(#it{D_{s}^{+}})", // Ds_ipsig
    "IP significance(#it{#mu^{-}})", // mu_ipsig

    "#it{D_{s}^{+}} FD in xy [cm]",  // Ds_FDxy
    "#it{D_{s}^{+}} FD [cm]",        // Ds_FD

    "#it{D_{s}^{+}} FD significance in xy", // Ds_FDxysig
    "#it{D_{s}^{+}} FD significance",       // Ds_FDsig

    "#it{#DeltaR}(#it{K^{+}}, #it{K^{-}})", // dR_Kp_Km
    "#it{#DeltaR}(#it{K^{+}}, #it{#pi^{+}})", // dR_Kp_pi
    "#it{#DeltaR}(#it{K^{+}}, #it{#phi})",   // dR_Kp_phi
    "#it{#DeltaR}(#it{K^{+}}, #it{D_{s}^{+}})", // dR_Kp_Ds
    "#it{#DeltaR}(#it{K^{+}}, #it{#mu^{-}})", // dR_Kp_mu

    "#it{#DeltaR}(#it{K^{-}}, #it{#pi^{+}})", // dR_Km_pi
    "#it{#DeltaR}(#it{K^{-}}, #it{#phi})",   // dR_Km_phi
    "#it{#DeltaR}(#it{K^{-}}, #it{D_{s}^{+}})", // dR_Km_Ds
    "#it{#DeltaR}(#it{K^{-}}, #it{#mu^{-}})", // dR_Km_mu

    "#it{#DeltaR}(#it{#pi^{+}}, #it{#phi})",  // dR_pi_phi
    "#it{#DeltaR}(#it{#pi^{+}}, #it{D_{s}^{+}})", // dR_pi_Ds
    "#it{#DeltaR}(#it{#pi^{+}}, #it{#mu^{-}})", // dR_pi_mu

    "#it{#DeltaR}(#it{#phi}, #it{D_{s}^{+}})", // dR_phi_Ds
    "#it{#DeltaR}(#it{#phi}, #it{#mu^{-}})",  // dR_phi_mu
    
    "#it{#DeltaR}(#it{D_{s}^{+}}, #it{#mu^{-}})" // dR_Ds_mu
};

int compare_vars()
{

    lumi_sqrtS = "13.6 TeV, 2022";
    setTDRStyle();

    const TString figdir = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/figures/";
    if (gSystem->AccessPathName(figdir))
        gSystem->MakeDirectory(figdir);

    TChain *ChannelChain[NCOMB] = {nullptr};

    for (Int_t cid = 0; cid < NCOMB; ++cid)
    {
        ChannelChain[cid] = new TChain("FullReco/Events");
        ChannelChain[cid]->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/tuples/" + ChannelSample[cid] + ".root");
    }

    for (Int_t hid = 0; hid < NHIST; ++hid)
    {
        TH1F *ChannelHist[NCOMB] = {nullptr};
        TCanvas *canvas = new TCanvas(Form("canvas_%d", hid), "canvas", 800, 600);
        canvas_setup(canvas);
        canvas->SetRightMargin(0.13);
        TLegend *HistLegend = new TLegend(0.88, 0.5, 1, 0.9, NULL, "brNDC");
        HistLegend->SetTextFont(42);
        HistLegend->SetTextColor(kBlack);
        HistLegend->SetBorderSize(0);
        HistLegend->SetFillColor(0);
        HistLegend->SetLineWidth(0);
        HistLegend->SetTextSize(0.04);

        Double_t height = 0;
        Double_t ground = 100;

        for (Int_t cid = 0; cid < NCOMB; ++cid)
        {
            ChannelHist[cid] = new TH1F(Form("hist_%s_%d", ChannelSample[cid].Data(), hid), "", nbins[hid], varmin[hid], varmax[hid]);
            ChannelChain[cid]->Project(Form("hist_%s_%d", ChannelSample[cid].Data(), hid), HistVar[hid], "num_ss_mu_tight==0 && num_os_mu_tight==1");
            ChannelHist[cid]->Scale(1. / ChannelHist[cid]->Integral());
            HistLegend->AddEntry(ChannelHist[cid], ChannelLegend[cid], "l");
            ChannelHist[cid]->SetLineColor(LineColors[cid]);
            ChannelHist[cid]->SetLineWidth(2);

            height = std::max(ChannelHist[cid]->GetMaximum(), height);
            ground = std::min(ChannelHist[cid]->GetMinimum(), ground);

            if (cid <= NCOMB / 2)
                ChannelHist[cid]->SetLineStyle(1);
            else
                ChannelHist[cid]->SetLineStyle(2);

            if (cid == 0)
            {
                ChannelHist[cid]->GetXaxis()->SetTitle(HistXTitle[hid]);
                ChannelHist[cid]->GetYaxis()->SetTitle("# normalized events");
                ChannelHist[cid]->Draw("hist");
            }
            else
                ChannelHist[cid]->Draw("hist same");
        }

        ChannelHist[0]->SetMaximum(height * 1.5);
        ChannelHist[0]->SetMinimum(0);
        HistLegend->Draw();

        CMS_lumi(canvas);
        canvas->Update();
        canvas->RedrawAxis();
        canvas->SaveAs(figdir + HistVar[hid] + ".png");

        ChannelHist[0]->SetMaximum(height * 50);
        ChannelHist[0]->SetMinimum(ground / 50);
        canvas->SetLogy(1);
        canvas->Update();
        canvas->RedrawAxis();
        canvas->SaveAs(figdir + HistVar[hid] + "_log.png");
    }

    return 0;
}

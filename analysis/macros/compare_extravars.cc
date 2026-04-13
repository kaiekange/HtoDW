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

const Color_t LineColors[] = {
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
    K1_pt = 0, // K1_pt
    K2_pt,     // K2_pt
    pi_pt,     // pi_pt
    phi_pt,    // phi_pt
    Ds_pt,     // Ds_pt
    mu_pt,     // mu_pt

    K1_pt_minus_K2_pt, // K1_pt - K2_pt
    phi_pt_div_Ds_pt,  // phi_pt/Ds_pt
    Ds_pt_div_H_pt,    // Ds_pt/H_pt

    Ds_IsoR03_RelIso,       // Ds_IsoR03_RelIso
    Ds_IsoR03_RelIso_large, // Ds_IsoR03_RelIso
    Ds_IsoR04_RelIso,       // Ds_IsoR04_RelIso
    Ds_IsoR04_RelIso_large, // Ds_IsoR04_RelIso
    mu_IsoR03_RelIso,       // mu_IsoR03_RelIso
    mu_IsoR04_RelIso,       // mu_IsoR04_RelIso
    mu_IsoR03_PFIso,        // mu_IsoR03_PFIso
    mu_IsoR04_PFIso,        // mu_IsoR04_PFIso
    dira,                   // dira
    log_dira,               // log(dira)
    Ds_FD,                  // Ds_FD
    Ds_FDsig,               // Ds_FDsig
    log_Ds_FDsig,           // log(Ds_FDsig)

    Kp_ip,         // Kp_ip
    Km_ip,         // Km_ip
    pi_ip,         // pi_ip
    Ds_ip,         // Ds_ip
    mu_ip,         // mu_ip
    Kp_ipsig,      // Kp_ipsig
    Km_ipsig,      // Km_ipsig
    pi_ipsig,      // pi_ipsig
    Ds_ipsig,      // Ds_ipsig
    mu_ipsig,      // mu_ipsig
    min_ipsig,     // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)
    log_min_ipsig, // log min(min(Kp_ipsig, Km_ipsig), pi_ipsig)

    dR_Kp_Km,  // dR_Kp_Km
    dR_K1_pi,  // dR_K1_pi
    dR_K1_phi, // dR_K1_phi
    dR_K1_Ds,  // dR_K1_Ds
    dR_K1_mu,  // dR_K1_mu

    dR_K2_pi,  // dR_K2_pi
    dR_K2_phi, // dR_K2_phi
    dR_K2_Ds,  // dR_K2_Ds
    dR_K2_mu,  // dR_K2_mu

    dR_pi_phi, // dR_pi_phi
    dR_pi_Ds,  // dR_pi_Ds
    dR_pi_mu,  // dR_pi_mu

    dR_phi_Ds, // dR_phi_Ds
    dR_phi_mu, // dR_phi_mu

    dR_Ds_mu, // dR_Ds_mu

    met_pt,        // met_pt
    met_sumEt,     // met_sumEt
    num_jets_pt30, // num_jets_pt30

    W_masst, // W_masst
    H_mass,  // H_mass
    H_pt,    // H_pt

    phi_mass, // phi_mass
    Ds_mass,  // Ds_mass
    NHIST
};

const TString HistVar[NHIST] = {
    "K1_pt",
    "K2_pt",
    "pi_pt",
    "phi_pt",
    "Ds_pt",
    "mu_pt",

    "K1_pt - K2_pt",
    "phi_pt/Ds_pt",
    "Ds_pt/H_pt",

    "Ds_IsoR03_RelIso",
    "Ds_IsoR03_RelIso",
    "Ds_IsoR04_RelIso",
    "Ds_IsoR04_RelIso",
    "mu_IsoR03_RelIso",
    "mu_IsoR04_RelIso",
    "mu_IsoR03_PFIso",
    "mu_IsoR04_PFIso",
    "dira",
    "log(dira)",
    "Ds_FD",
    "Ds_FDsig",
    "log(Ds_FDsig)",

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
    "min(min(Kp_ipsig, Km_ipsig), pi_ipsig)",
    "log(min(min(Kp_ipsig, Km_ipsig), pi_ipsig))",

    "dR_Kp_Km",
    "dR_K1_pi",
    "dR_K1_phi",
    "dR_K1_Ds",
    "dR_K1_mu",

    "dR_K2_pi",
    "dR_K2_phi",
    "dR_K2_Ds",
    "dR_K2_mu",

    "dR_pi_phi",
    "dR_pi_Ds",
    "dR_pi_mu",

    "dR_phi_Ds",
    "dR_phi_mu",

    "dR_Ds_mu",

    "met_pt",
    "met_sumEt",
    "num_jets_pt30",

    "W_masst",
    "H_mass",
    "H_pt",

    "phi_mass", // phi_mass
    "Ds_mass",  // Ds_mass
};

const TString VarFigName[NHIST] = {
    "K1_pt",  // K1_pt
    "K2_pt",  // K2_pt
    "pi_pt",  // pi_pt
    "phi_pt", // phi_pt
    "Ds_pt",  // Ds_pt
    "mu_pt",  // mu_pt

    "K1_pt_minus_K2_pt", // K1_pt - K2_pt
    "phi_pt_div_Ds_pt",  // phi_pt/Ds_pt
    "Ds_pt_div_H_pt",    // Ds_pt/H_pt

    "Ds_IsoR03_RelIso",       // Ds_IsoR03_RelIso
    "Ds_IsoR03_RelIso_large", // Ds_IsoR03_RelIso
    "Ds_IsoR04_RelIso",       // Ds_IsoR04_RelIso
    "Ds_IsoR04_RelIso_large", // Ds_IsoR04_RelIso
    "mu_IsoR03_RelIso",       // mu_IsoR03_RelIso
    "mu_IsoR04_RelIso",       // mu_IsoR04_RelIso
    "mu_IsoR03_PFIso",        // mu_IsoR03_PFIso
    "mu_IsoR04_PFIso",        // mu_IsoR04_PFIso
    "dira",                   // dira
    "log_dira",               // log(dira)
    "Ds_FD",                  // Ds_FD
    "Ds_FDsig",               // Ds_FDsig
    "log_Ds_FDsig",           // log(Ds_FDsig)

    "Kp_ip",         // Kp_ip
    "Km_ip",         // Km_ip
    "pi_ip",         // pi_ip
    "Ds_ip",         // Ds_ip
    "mu_ip",         // mu_ip
    "Kp_ipsig",      // Kp_ipsig
    "Km_ipsig",      // Km_ipsig
    "pi_ipsig",      // pi_ipsig
    "Ds_ipsig",      // Ds_ipsig
    "mu_ipsig",      // mu_ipsig
    "min_ipsig",     // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)
    "log_min_ipsig", // log min(min(Kp_ipsig, Km_ipsig), pi_ipsig)

    "dR_Kp_Km",  // dR_Kp_Km
    "dR_K1_pi",  // dR_K1_pi
    "dR_K1_phi", // dR_K1_phi
    "dR_K1_Ds",  // dR_K1_Ds
    "dR_K1_mu",  // dR_K1_mu

    "dR_K2_pi",  // dR_K2_pi
    "dR_K2_phi", // dR_K2_phi
    "dR_K2_Ds",  // dR_K2_Ds
    "dR_K2_mu",  // dR_K2_mu

    "dR_pi_phi", // dR_pi_phi
    "dR_pi_Ds",  // dR_pi_Ds
    "dR_pi_mu",  // dR_pi_mu

    "dR_phi_Ds", // dR_phi_Ds
    "dR_phi_mu", // dR_phi_mu

    "dR_Ds_mu", // dR_Ds_mu

    "met_pt",        // met_pt
    "met_sumEt",     // met_sumEt
    "num_jets_pt30", // num_jets_pt30

    "W_masst", // W_masst
    "H_mass",  // H_mass
    "H_pt",    // H_pt

    "phi_mass", // phi_mass
    "Ds_mass",  // Ds_mass
};

const Float_t varmin[NHIST] = {
    0,  // K1_pt
    0,  // K2_pt
    0,  // pi_pt
    0,  // phi_pt
    0,  // Ds_pt
    0, // mu_pt

    0, // K1_pt - K2_pt
    0, // phi_pt/Ds_pt
    0, // Ds_pt/H_pt

    0,      // Ds_IsoR03_RelIso
    0,      // Ds_IsoR03_RelIso
    0,      // Ds_IsoR04_RelIso
    0,      // Ds_IsoR04_RelIso
    0,      // mu_IsoR03_RelIso
    0,      // mu_IsoR04_RelIso
    0,      // mu_IsoR03_PFIso
    0,      // mu_IsoR04_PFIso
    0.9,    // dira
    -0.106, // log(dira)
    0,      // Ds_FD
    0,      // Ds_FDsig
    -5,     // log(Ds_FDsig)

    0,  // Kp_ip
    0,  // Km_ip
    0,  // pi_ip
    0,  // Ds_ip
    0,  // mu_ip
    0,  // Kp_ipsig
    0,  // Km_ipsig
    0,  // pi_ipsig
    0,  // Ds_ipsig
    0,  // mu_ipsig
    0,  // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)
    -5, // log min(min(Kp_ipsig, Km_ipsig), pi_ipsig)

    0, // dR_Kp_Km
    0, // dR_K1_pi
    0, // dR_K1_phi
    0, // dR_K1_Ds
    0, // dR_K1_mu

    0, // dR_K2_pi
    0, // dR_K2_phi
    0, // dR_K2_Ds
    0, // dR_K2_mu

    0, // dR_pi_phi
    0, // dR_pi_Ds
    0, // dR_pi_mu

    0, // dR_phi_Ds
    0, // dR_phi_mu

    0, // dR_Ds_mu

    0, // met_pt
    0, // met_sumEt
    0, // num_jets_pt30

    0, // W_masst
    0, // H_mass
    0, // H_pt

    1.019460-0.03,
    1.96835-0.15,
};

const Float_t varmax[NHIST] = {
    30,  // K1_pt
    30,  // K2_pt
    40,  // pi_pt
    50,  // phi_pt
    70,  // Ds_pt
    150, // mu_pt

    15, // K1_pt - K2_pt
    1,  // phi_pt/Ds_pt
    5,  // Ds_pt/H_pt

    0.5,    // Ds_IsoR03_RelIso
    10,     // Ds_IsoR03_RelIso
    0.5,    // Ds_IsoR04_RelIso
    10,     // Ds_IsoR04_RelIso
    10,     // mu_IsoR03_RelIso
    10,     // mu_IsoR04_RelIso
    10,     // mu_IsoR03_PFIso
    10,     // mu_IsoR04_PFIso
    0.9,    // dira
    -0.106, // log(dira)
    3,      // Ds_FD
    50,      // Ds_FDsig
    8,     // log(Ds_FDsig)

    0.05, // Kp_ip
    0.05, // Km_ip
    0.1,  // pi_ip
    0.02, // Ds_ip
    0.01, // mu_ip
    20,   // Kp_ipsig
    20,   // Km_ipsig
    30,   // pi_ipsig
    10,   // Ds_ipsig
    8,    // mu_ipsig
    20,   // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)
    8,    // log min(min(Kp_ipsig, Km_ipsig), pi_ipsig)

    0.2, // dR_Kp_Km
    0.8, // dR_K1_pi
    0.1, // dR_K1_phi
    0.3, // dR_K1_Ds
    5,   // dR_K1_mu

    0.8, // dR_K2_pi
    0.1, // dR_K2_phi
    0.3, // dR_K2_Ds
    5,   // dR_K2_mu

    1, // dR_pi_phi
    1, // dR_pi_Ds
    5, // dR_pi_mu

    0.3, // dR_phi_Ds
    5,   // dR_phi_mu

    5, // dR_Ds_mu

    250,  // met_pt
    1000, // met_sumEt
    8,    // num_jets_pt30

    250, // W_masst
    125, // H_mass
    250, // H_pt

    1.019460+0.03,
    1.96835+0.15,
};

const TString HistXTitle[NHIST] = {
    "#it{p_{T}}(#it{K_{1}}) [GeV]",     // K1_pt
    "#it{p_{T}}(#it{K_{2}}) [GeV]",     // K2_pt
    "#it{p_{T}}(#it{#pi^{+}}) [GeV]",   // pi_pt
    "#it{p_{T}}(#it{#phi}) [GeV]",      // phi_pt
    "#it{p_{T}}(#it{D_{s}^{+}}) [GeV]", // Ds_pt
    "#it{p_{T}}(#it{#mu^{-}}) [GeV]",   // mu_pt

    "#it{p_{T}(K_{1})-p_{T}(K_{2})}",                // K1_pt - K2_pt
    "#it{p_{T}(#phi)/p_{T}(D_{s}^{+})}",             // phi_pt/Ds_pt
    "#it{p_{T}(D_{s}^{+})/p_{T}(D_{s}^{+}#mu^{-})}", // Ds_pt/H_pt

    "#it{D_{s}^{+}} IsoR03 RelIso",        // Ds_IsoR04_RelIso
    "#it{D_{s}^{+}} IsoR03 RelIso",        // Ds_IsoR04_RelIso
    "#it{D_{s}^{+}} IsoR04 RelIso",        // Ds_IsoR04_RelIso
    "#it{D_{s}^{+}} IsoR04 RelIso",        // Ds_IsoR04_RelIso
    "#it{#mu^{-}} IsoR03 RelIso",          // mu_IsoR03_RelIso
    "#it{#mu^{-}} IsoR04 RelIso",          // mu_IsoR04_RelIso
    "#it{#mu^{-}} IsoR03 PFIso",           // mu_IsoR03_PFIso
    "#it{#mu^{-}} IsoR04 PFIso",           // mu_IsoR04_PFIso
    "cos #alpha(PV, SV)",                  // dira
    "log(cos #alpha(PV, SV))",             // log(dira)
    "#it{D_{s}^{+}} FD [cm]",              // Ds_FD
    "#it{D_{s}^{+}} FD significance",      // Ds_FDsig
    "log(#it{D_{s}^{+}} FD significance)", // log(Ds_FDsig)

    "IP(#it{K^{+}})",                                       // Kp_ip
    "IP(#it{K^{-}})",                                       // Km_ip
    "IP(#it{#pi^{+}})",                                     // pi_ip
    "IP(#it{D_{s}^{+}})",                                   // Ds_ip
    "IP(#it{#mu^{-}})",                                     // mu_ip
    "IP significance(#it{K^{+}})",                          // Kp_ipsig
    "IP significance(#it{K^{-}})",                          // Km_ipsig
    "IP significance(#it{#pi^{+}})",                        // pi_ipsig
    "IP significance(#it{D_{s}^{+}})",                      // Ds_ipsig
    "IP significance(#it{#mu^{-}})",                        // mu_ipsig
    "min IP significance of #it{K^{+},K^{-},#pi^{+}}",      // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)
    "log(min IP significance of #it{K^{+},K^{-},#pi^{+}})", // min(min(Kp_ipsig, Km_ipsig), pi_ipsig)

    "#it{#DeltaR}(#it{K_{1}}, #it{K_{2}})",     // dR_Kp_Km
    "#it{#DeltaR}(#it{K_{1}}, #it{#pi^{+}})",   // dR_K1_pi
    "#it{#DeltaR}(#it{K_{1}}, #it{#phi})",      // dR_K1_phi
    "#it{#DeltaR}(#it{K_{1}}, #it{D_{s}^{+}})", // dR_K1_Ds
    "#it{#DeltaR}(#it{K_{1}}, #it{#mu^{-}})",   // dR_K1_mu

    "#it{#DeltaR}(#it{K_{2}}, #it{#pi^{+}})",   // dR_K2_pi
    "#it{#DeltaR}(#it{K_{2}}, #it{#phi})",      // dR_K2_phi
    "#it{#DeltaR}(#it{K_{2}}, #it{D_{s}^{+}})", // dR_K2_Ds
    "#it{#DeltaR}(#it{K_{2}}, #it{#mu^{-}})",   // dR_K2_mu

    "#it{#DeltaR}(#it{#pi^{+}}, #it{#phi})",      // dR_pi_phi
    "#it{#DeltaR}(#it{#pi^{+}}, #it{D_{s}^{+}})", // dR_pi_Ds
    "#it{#DeltaR}(#it{#pi^{+}}, #it{#mu^{-}})",   // dR_pi_mu

    "#it{#DeltaR}(#it{#phi}, #it{D_{s}^{+}})", // dR_phi_Ds
    "#it{#DeltaR}(#it{#phi}, #it{#mu^{-}})",   // dR_phi_mu

    "#it{#DeltaR}(#it{D_{s}^{+}}, #it{#mu^{-}})", // dR_Ds_mu

    "MET #it{p_{T}} [GeV]",         // met_pt
    "#it{#sum E_{T}} [GeV]",    // met_sumEt
    "# jets (#it{p_{T}} > 30 GeV)", // num_jets_pt30

    "#it{M_{T}(W^{-})} [GeV]",            // W_masst
    "#it{M(D_{s}^{+}#mu^{-})} [GeV]",     // H_mass
    "#it{p_{T}(D_{s}^{+}#mu^{-})} [GeV]", // H_pt

    "#it{M(K^{+}K^{-})} [GeV]", // phi_mass
    "#it{M(K^{+}K^{-}#pi^{+})} [GeV]", // Ds_mass
};

Int_t compare_extravars()
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

    for (Int_t hid = 5; hid < 6; ++hid)
    // for (Int_t hid = NHIST - 7; hid < NHIST - 6; ++hid)
    {
        Int_t nbins = 100;

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
            TCut mycut;
            if (cid == 0)
                mycut = "event_match && num_ss_mu_tight==0 && num_os_mu_tight==1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";
            else
                mycut = "num_ss_mu==0 && num_os_mu==1 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125";

            ChannelHist[cid] = new TH1F(Form("hist_%s_%d", ChannelSample[cid].Data(), hid), "", nbins, varmin[hid], varmax[hid]);
            ChannelChain[cid]->Project(Form("hist_%s_%d", ChannelSample[cid].Data(), hid), HistVar[hid], mycut);
            ChannelHist[cid]->Scale(1. / ChannelHist[cid]->Integral());
            HistLegend->AddEntry(ChannelHist[cid], ChannelLegend[cid], "l");
            ChannelHist[cid]->SetLineColor(LineColors[cid]);
            ChannelHist[cid]->SetLineWidth(2);

            height = std::max(ChannelHist[cid]->GetMaximum(), height);
            ground = std::min(ChannelHist[cid]->GetMinimum(), ground);

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
        canvas->SaveAs(figdir + VarFigName[hid] + ".png");

        // ChannelHist[0]->SetMaximum(height * 50);
        // ChannelHist[0]->SetMinimum(0.0001);
        // canvas->SetLogy(1);
        // canvas->Update();
        // canvas->RedrawAxis();
        // canvas->SaveAs(figdir + VarFigName[hid] + "_log.png");
    }

    return 0;
}

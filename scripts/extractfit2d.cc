#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

int extractfit2d() {
    setTDRStyle();
    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TChain *mychain = new TChain("Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/all_skimmed.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/fit_compare")) gSystem->MakeDirectory("./figures/fit_compare");

    TCut mass_cut = "DsFit_Ds_M > 1.85 && DsFit_Ds_M < 2.1 && phiFit_phi_M > 0.99 && phiFit_phi_M < 1.05";
    TCut match_sel = "match_entry == 1";
    TCut nonmatch_sel = "non_match_entry == 1";
    match_sel += mass_cut;
    nonmatch_sel += mass_cut;

    TTree *match_tree = mychain->CopyTree(match_sel);
    TTree *nonmatch_tree = mychain->CopyTree(nonmatch_sel);

    double nentries = nonmatch_tree->GetEntries();

    RooRealVar DsFit_Ds_M("DsFit_Ds_M", "DsFit_Ds_M", 1.85, 2.1);
    RooRealVar phiFit_phi_M("phiFit_phi_M", "phiFit_phi_M", 0.99, 1.05);
    RooDataSet *fit_dataset = new RooDataSet("fit_dataset", "fit_dataset", nonmatch_tree, RooArgSet(DsFit_Ds_M, phiFit_phi_M));

    RooDataSet *Ds_dataset = new RooDataSet("Ds_dataset", "Ds_dataset", match_tree, RooArgSet(DsFit_Ds_M));
    RooDataSet *phi_dataset = new RooDataSet("phi_dataset", "phi_dataset", match_tree, RooArgSet(phiFit_phi_M));

    RooKeysPdf Ds_peak("Ds_peak", "Ds_peak", DsFit_Ds_M, *Ds_dataset);
    RooKeysPdf phi_peak("phi_peak", "phi_peak", phiFit_phi_M, *phi_dataset);

    RooRealVar Ds_c0("Ds_c0", "Ds_c0", 0, -1, 1) ;
    RooRealVar Ds_c1("Ds_c1", "Ds_c1", 0, -1, 1) ;
    RooRealVar Ds_c2("Ds_c2", "Ds_c2", 0, -1, 1) ;
    RooChebychev Ds_chebychev("Ds_chebychev", "Ds_chebychev", DsFit_Ds_M, RooArgList(Ds_c0,Ds_c1,Ds_c2));

    RooRealVar phi_c0("phi_c0", "phi_c0", 0, -1, 1) ;
    RooRealVar phi_c1("phi_c1", "phi_c1", 0, -1, 1) ;
    RooRealVar phi_c2("phi_c2", "phi_c2", 0, -1, 1) ;
    RooChebychev phi_chebychev("phi_chebychev", "phi_chebychev", phiFit_phi_M, RooArgList(phi_c0,phi_c1,phi_c2));

    RooRealVar n_Ds_peak_phi_peak("n_Ds_peak_phi_peak", "n_Ds_peak_phi_peak", nentries*0.5, 0, nentries*1.1);
    RooRealVar n_Ds_chebychev_phi_peak("n_Ds_chebychev_phi_peak", "n_Ds_chebychev_phi_peak", nentries*0.5, 0, nentries*1.1);
    RooRealVar n_Ds_peak_phi_chebychev("n_Ds_peak_phi_chebychev", "n_Ds_peak_phi_chebychev", nentries*0.5, 0, nentries*1.1);
    RooRealVar n_Ds_chebychev_phi_chebychev("n_Ds_chebychev_phi_chebychev", "n_Ds_chebychev_phi_chebychev", nentries*0.5, 0, nentries*1.1);

    RooProdPdf Ds_peak_phi_peak("Ds_peak_phi_peak", "Ds_peak_phi_peak", RooArgList(Ds_peak, phi_peak));
    RooProdPdf Ds_chebychev_phi_peak("Ds_chebychev_phi_peak", "Ds_chebychev_phi_peak", RooArgList(Ds_chebychev, phi_peak));
    RooProdPdf Ds_peak_phi_chebychev("Ds_peak_phi_chebychev", "Ds_peak_phi_chebychev", RooArgList(Ds_peak, phi_chebychev));
    RooProdPdf Ds_chebychev_phi_chebychev("Ds_chebychev_phi_chebychev", "Ds_chebychev_phi_chebychev", RooArgList(Ds_chebychev, phi_chebychev));

    RooAddPdf model("model", "model", RooArgList(Ds_peak_phi_peak, Ds_chebychev_phi_peak, Ds_peak_phi_chebychev, Ds_chebychev_phi_chebychev), RooArgList(n_Ds_peak_phi_peak, n_Ds_chebychev_phi_peak, n_Ds_peak_phi_chebychev, n_Ds_chebychev_phi_chebychev));

    RooFitResult *fitResult = model.fitTo(*fit_dataset, RooFit::Save(true), RooFit::Extended(true));
    fitResult->Print();

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 700);
    canvas_setup(c1);
    c1->Divide(1,2);
    c1->cd(1);
    TVirtualPad* c11 = c1->GetPad(1);
    canvas_setup_sub(c11);
    c11->SetPad(0,0.3,1,1);
    c11->SetBottomMargin(0.0);
    c11->SetRightMargin(0.05);
    c11->SetLogy(0);
    c11->SetFillColor(0);
    c11->SetFrameFillColor(0);

    RooPlot *Ds_frame = new RooPlot("", "", DsFit_Ds_M, 1.85, 2.1, 100);
    fit_dataset->plotOn(Ds_frame, RooFit::Name("nonmatch Ds"), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(1.1), RooFit::Binning(100), RooFit::DrawOption("ep"));
    model.plotOn(Ds_frame, RooFit::Name("model"), RooFit::Components("model"), RooFit::LineStyle(9), RooFit::LineColor(kBlue), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(Ds_frame, RooFit::Name("Ds_peak_phi_peak"), RooFit::Components("Ds_peak_phi_peak"), RooFit::LineStyle(9), RooFit::LineColor(kRed), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(Ds_frame, RooFit::Name("Ds_chebychev_phi_peak"), RooFit::Components("Ds_chebychev_phi_peak"), RooFit::LineStyle(9), RooFit::LineColor(kViolet), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(Ds_frame, RooFit::Name("Ds_peak_phi_chebychev"), RooFit::Components("Ds_peak_phi_chebychev"), RooFit::LineStyle(9), RooFit::LineColor(kGreen+1), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(Ds_frame, RooFit::Name("Ds_chebychev_phi_chebychev"), RooFit::Components("Ds_chebychev_phi_chebychev"), RooFit::LineStyle(9), RooFit::LineColor(kOrange-3), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));

    Ds_frame->Draw("");
    Ds_frame->SetYTitle("# non-matched candidates");
    Ds_frame->SetMinimum(0.01);
    Ds_frame->SetMaximum(Ds_frame->GetMaximum()*1.3);

    TLegend *Ds_legend = new TLegend(0.2, 0.2, 0.4, 0.4, NULL, "brNDC");
    Ds_legend->AddEntry(Ds_frame->findObject("Non-matched #it{D_{s}^{+}} candidates"), "nonmatch Ds", "ep");
    Ds_legend->AddEntry(Ds_frame->findObject("model"), "Fit model", "l");
    Ds_legend->AddEntry(Ds_frame->findObject("Ds_peak_phi_peak"), "#it{D_{s}^{+}} peak, #it{#phi} peak", "l");
    Ds_legend->AddEntry(Ds_frame->findObject("Ds_chebychev_phi_peak"), "#it{D_{s}^{+}} chebychev, #it{#phi} peak", "l");
    Ds_legend->AddEntry(Ds_frame->findObject("Ds_peak_phi_chebychev"), "#it{D_{s}^{+}} peak, #it{#phi} chebychev", "l");
    Ds_legend->AddEntry(Ds_frame->findObject("Ds_chebychev_phi_chebychev"), "#it{D_{s}^{+}} chebychev, #it{#phi} chebychev", "l");
    Ds_legend->SetBorderSize(0);
    Ds_legend->SetFillColor(0);
    Ds_legend->SetLineWidth(0);
    Ds_legend->SetTextSize(0.03);
    Ds_legend->Draw();
    
    write_text(0.4, 0.8, Form("Fit Results:"));
    write_text(0.4, 0.75, Form("# #it{D_{s}^{+}} peak, #it{#phi} peak = %.0f #pm %.0f", n_Ds_peak_phi_peak.getVal(), n_Ds_peak_phi_peak.getError()));
    write_text(0.4, 0.7, Form("# #it{D_{s}^{+}} chebychev, #it{#phi} peak = %.0f #pm %.0f", n_Ds_chebychev_phi_peak.getVal(), n_Ds_chebychev_phi_peak.getError()));
    write_text(0.4, 0.65, Form("# #it{D_{s}^{+}} peak, #it{#phi} chebychev = %.0f #pm %.0f", n_Ds_peak_phi_chebychev.getVal(), n_Ds_peak_phi_chebychev.getError()));
    write_text(0.4, 0.6, Form("# #it{D_{s}^{+}} chebychev, #it{#phi} chebychev = %.0f #pm %.0f", n_Ds_chebychev_phi_chebychev.getVal(), n_Ds_chebychev_phi_chebychev.getError()));
    CMS_lumi_sub(c11);

    c1->cd(2);
    TVirtualPad* c12 = c1->GetPad(2);
    canvas_setup_sub(c12);
    c12->SetPad(0,0,1,0.3);
    c12->SetRightMargin(0.05);
    c12->SetTopMargin(0.05);
    c12->SetBottomMargin(0.4);
    RooHist* Ds_hpull = Ds_frame->pullHist("nonmatch Ds","model");
    RooPlot* Ds_pull = new RooPlot("", "", DsFit_Ds_M, 1.85, 2.1, 100);
    Ds_hpull->SetFillColor(1);
    Ds_pull->addPlotable(Ds_hpull, "BX");
    Ds_pull->Draw();
    TAxis* Ds_xpull = Ds_pull->GetXaxis();
    TAxis* Ds_ypull = Ds_pull->GetYaxis();
    Ds_xpull->SetTitleOffset(3.2);
    Ds_xpull->SetTickLength(0.06);
    Ds_ypull->SetNdivisions(504);
    Ds_pull->SetXTitle("#it{M}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]");
    /* Ds_pull->GetXaxis()->CenterTitle(true); */
    Ds_pull->GetXaxis()->SetTitleOffset(1.);
    Ds_pull->SetYTitle("pull");
    Ds_pull->GetYaxis()->CenterTitle(true);
    Ds_pull->GetYaxis()->SetRangeUser(-5, 5);
    Ds_pull->GetYaxis()->SetTitleOffset(0.3);
    Ds_pull->GetXaxis()->SetTitleSize(0.1);
    Ds_pull->GetXaxis()->SetLabelSize(0.12); 
    Ds_pull->GetYaxis()->SetTitleSize(0.12);
    Ds_pull->GetYaxis()->SetLabelSize(0.12);
    c1->Update();
    c1->SaveAs("./figures/fit_compare/2D_fit_Ds.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 700);
    canvas_setup(c2);
    c2->Divide(1,2);
    c2->cd(1);
    TVirtualPad* c21 = c2->GetPad(1);
    canvas_setup_sub(c21);
    c21->SetPad(0,0.3,1,1);
    c21->SetBottomMargin(0.0);
    c21->SetRightMargin(0.05);
    c21->SetLogy(0);
    c21->SetFillColor(0);
    c21->SetFrameFillColor(0);

    RooPlot *phi_frame = new RooPlot("", "", phiFit_phi_M, 0.99, 1.05, 100);
    fit_dataset->plotOn(phi_frame, RooFit::Name("nonmatch phi"), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(1.1), RooFit::Binning(100), RooFit::DrawOption("ep"));
    model.plotOn(phi_frame, RooFit::Name("model"), RooFit::Components("model"), RooFit::LineStyle(9), RooFit::LineColor(kBlue), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(phi_frame, RooFit::Name("Ds_peak_phi_peak"), RooFit::Components("Ds_peak_phi_peak"), RooFit::LineStyle(9), RooFit::LineColor(kRed), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(phi_frame, RooFit::Name("Ds_chebychev_phi_peak"), RooFit::Components("Ds_chebychev_phi_peak"), RooFit::LineStyle(9), RooFit::LineColor(kViolet), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(phi_frame, RooFit::Name("Ds_peak_phi_chebychev"), RooFit::Components("Ds_peak_phi_chebychev"), RooFit::LineStyle(9), RooFit::LineColor(kGreen+1), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(phi_frame, RooFit::Name("Ds_chebychev_phi_chebychev"), RooFit::Components("Ds_chebychev_phi_chebychev"), RooFit::LineStyle(9), RooFit::LineColor(kOrange-3), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));

    phi_frame->Draw("");
    phi_frame->SetYTitle("# non-matched candidates");
    phi_frame->SetMinimum(0.01);
    phi_frame->SetMaximum(phi_frame->GetMaximum()*1.3);

    TLegend *phi_legend = new TLegend(0.2, 0.2, 0.4, 0.4, NULL, "brNDC");
    phi_legend->AddEntry(phi_frame->findObject("Non-matched #it{#phi} candidates"), "nonmatch phi", "ep");
    phi_legend->AddEntry(phi_frame->findObject("model"), "Fit model", "l");
    phi_legend->AddEntry(phi_frame->findObject("Ds_peak_phi_peak"), "#it{D_{s}^{+}} peak, #it{#phi} peak", "l");
    phi_legend->AddEntry(phi_frame->findObject("Ds_chebychev_phi_peak"), "#it{D_{s}^{+}} chebychev, #it{#phi} peak", "l");
    phi_legend->AddEntry(phi_frame->findObject("Ds_peak_phi_chebychev"), "#it{D_{s}^{+}} peak, #it{#phi} chebychev", "l");
    phi_legend->AddEntry(phi_frame->findObject("Ds_chebychev_phi_chebychev"), "#it{D_{s}^{+}} chebychev, #it{#phi} chebychev", "l");
    phi_legend->SetBorderSize(0);
    phi_legend->SetFillColor(0);
    phi_legend->SetLineWidth(0);
    phi_legend->SetTextSize(0.03);
    phi_legend->Draw();
    
    write_text(0.4, 0.8, Form("Fit Results:"));
    write_text(0.4, 0.75, Form("# #it{D_{s}^{+}} peak, #it{#phi} peak = %.0f #pm %.0f", n_Ds_peak_phi_peak.getVal(), n_Ds_peak_phi_peak.getError()));
    write_text(0.4, 0.7, Form("# #it{D_{s}^{+}} chebychev, #it{#phi} peak = %.0f #pm %.0f", n_Ds_chebychev_phi_peak.getVal(), n_Ds_chebychev_phi_peak.getError()));
    write_text(0.4, 0.65, Form("# #it{D_{s}^{+}} peak, #it{#phi} chebychev = %.0f #pm %.0f", n_Ds_peak_phi_chebychev.getVal(), n_Ds_peak_phi_chebychev.getError()));
    write_text(0.4, 0.6, Form("# #it{D_{s}^{+}} chebychev, #it{#phi} chebychev = %.0f #pm %.0f", n_Ds_chebychev_phi_chebychev.getVal(), n_Ds_chebychev_phi_chebychev.getError()));
    CMS_lumi_sub(c21);

    c2->cd(2);
    TVirtualPad* c22 = c2->GetPad(2);
    canvas_setup_sub(c22);
    c22->SetPad(0,0,1,0.3);
    c22->SetRightMargin(0.05);
    c22->SetTopMargin(0.05);
    c22->SetBottomMargin(0.4);
    RooHist* phi_hpull = phi_frame->pullHist("nonmatch phi","model");
    RooPlot* phi_pull = new RooPlot("", "", phiFit_phi_M, 0.99, 1.05, 100);
    phi_hpull->SetFillColor(1);
    phi_pull->addPlotable(phi_hpull, "BX");
    phi_pull->Draw();
    TAxis* phi_xpull = phi_pull->GetXaxis();
    TAxis* phi_ypull = phi_pull->GetYaxis();
    phi_xpull->SetTitleOffset(3.2);
    phi_xpull->SetTickLength(0.06);
    phi_ypull->SetNdivisions(504);
    phi_pull->SetXTitle("#it{M}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]");
    /* phi_pull->GetXaxis()->CenterTitle(true); */
    phi_pull->GetXaxis()->SetTitleOffset(1.);
    phi_pull->SetYTitle("pull");
    phi_pull->GetYaxis()->CenterTitle(true);
    phi_pull->GetYaxis()->SetRangeUser(-5, 5);
    phi_pull->GetYaxis()->SetTitleOffset(0.3);
    phi_pull->GetXaxis()->SetTitleSize(0.1);
    phi_pull->GetXaxis()->SetLabelSize(0.12); 
    phi_pull->GetYaxis()->SetTitleSize(0.12);
    phi_pull->GetYaxis()->SetLabelSize(0.12);
    c2->Update();
    c2->SaveAs("./figures/fit_compare/2D_fit_phi.png");
    delete mychain;

    return 0;
}

#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

void fit_draw_match(TTree *mytree, TString varname, double varmin, double varmax, TString xtitle, TString figpath) {

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h_temp = new TH1F("h_temp", "", 100, varmin, varmax);
    mytree->Project("h_temp", varname);
    double nentries = h_temp->GetEntries();

    RooRealVar fit_var(varname, varname, varmin, varmax);
    RooDataSet *dataset = new RooDataSet("dataset", "dataset", mytree, RooArgSet(fit_var));

    RooRealVar mu("mu", "mu", varmin, varmax);
    RooRealVar kuan("width", "width", h_temp->GetRMS(), 0., h_temp->GetRMS()*3);
    RooRealVar sigma("sigma", "sigma", h_temp->GetRMS(), 0., h_temp->GetRMS()*3);
    RooRealVar nsig("nsig", "nsig", nentries, 0., nentries*1.1);
    RooVoigtian model("signal", "signal", fit_var, mu, kuan, sigma);

    RooFitResult *fitResult = model.fitTo(*dataset, RooFit::Save(true));
    fitResult->Print();

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 700);
    canvas_setup(canvas);
    canvas->Divide(1,2);
    canvas->cd(1);
    TVirtualPad* canvas1 = canvas->GetPad(1);
    canvas_setup_sub(canvas1);
    canvas1->SetPad(0,0.3,1,1);
    canvas1->SetRightMargin(0.05);
    canvas1->SetBottomMargin(0.0);
    canvas1->SetLogy(0);
    canvas1->SetFillColor(0);
    canvas1->SetFrameFillColor(0);

    RooPlot *frame = fit_var.frame();
    dataset->plotOn(frame, RooFit::Name("Candidates"), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(1.1), RooFit::Binning(100), RooFit::DrawOption("ep"));
    model.plotOn(frame, RooFit::Name("Fit"), RooFit::Components("signal"), RooFit::LineStyle(9), RooFit::LineColor(kRed), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));

    frame->Draw("");
    frame->SetYTitle("# candidates");
	/* frame->GetYaxis()->SetTitleOffset(1.); */ 
	/* frame->GetYaxis()->CenterTitle(true); */
    frame->SetMinimum(0.01);
    frame->SetMaximum(h_temp->GetMaximum()*1.3);

    TLegend *legend = new TLegend(0.2, 0.5, 0.4, 0.7, NULL, "brNDC");
    legend->AddEntry(frame->findObject("Candidates"), "Candidates", "ep");
    legend->AddEntry(frame->findObject("Fit"), "Voigtian", "l");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.04);
    legend->Draw();
    
    write_text(0.6, 0.8, Form("Fit Results:"));
    write_text(0.6, 0.75, Form("#mu = %.3f #pm %.3f", mu.getVal(), mu.getError()));
    write_text(0.6, 0.7, Form("width = %.4f #pm %.4f", kuan.getVal(), kuan.getError()));
    write_text(0.6, 0.65, Form("#sigma = %.4f #pm %.4f", sigma.getVal(), sigma.getError()));
    write_text(0.6, 0.6, Form("# voigtian = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
    CMS_lumi_sub(canvas1);

    canvas->cd(2);
    TVirtualPad* canvas2 = canvas->GetPad(2);
    canvas_setup_sub(canvas2);
    canvas2->SetPad(0,0,1,0.3);
    canvas2->SetRightMargin(0.05);
    canvas2->SetTopMargin(0.05);
    canvas2->SetBottomMargin(0.4);
    RooHist* hpull = frame->pullHist("Candidates","Fit");
    RooPlot* pull = new RooPlot("", "", fit_var, varmin, varmax, 100);
    hpull->SetFillColor(1);
    pull->addPlotable(hpull, "BX");
    pull->Draw();
    TAxis* xpull = pull->GetXaxis();
    TAxis* ypull = pull->GetYaxis();
    xpull->SetTitleOffset(3.2);
    xpull->SetTickLength(0.06);
    ypull->SetNdivisions(504);
    pull->SetXTitle(xtitle);
    /* pull->GetXaxis()->CenterTitle(true); */
    pull->GetXaxis()->SetTitleOffset(1.);
    pull->SetYTitle("Pull");
    pull->GetYaxis()->CenterTitle(true);
    pull->GetYaxis()->SetRangeUser(-5, 5);
    pull->GetYaxis()->SetTitleOffset(0.3);
    pull->GetXaxis()->SetTitleSize(0.1);
    pull->GetXaxis()->SetLabelSize(0.12); 
    pull->GetYaxis()->SetTitleSize(0.12);
    pull->GetYaxis()->SetLabelSize(0.12);
    canvas->Update();
    canvas->SaveAs(figpath);
}

void fit_draw_nonmatch(TTree *mytree, TString varname, double varmin, double varmax, TString xtitle, double nsig_max, double nbkg_max, TString figpath) {

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH1F *h_temp = new TH1F("h_temp", "", 100, varmin, varmax);
    mytree->Project("h_temp", varname);
    double nentries = h_temp->GetEntries();

    RooRealVar fit_var(varname, varname, varmin, varmax);
    RooDataSet *dataset = new RooDataSet("dataset", "dataset", mytree, RooArgSet(fit_var));

    double kuan_max = std::min(h_temp->GetRMS(),0.02);
    double sigma_max = std::min(h_temp->GetRMS(),0.02);

    RooRealVar mu("mu", "mu", varmin, varmax);
    RooRealVar kuan("width", "width", kuan_max, 0., kuan_max*3);
    RooRealVar sigma("sigma", "sigma", sigma_max, 0., sigma_max*3);
    RooVoigtian sig("signal", "signal", fit_var, mu, kuan, sigma);

    RooRealVar c0("c0", "c0", 0, -1, 1) ;
    RooRealVar c1("c1", "c1", 0, -1, 1) ;
    RooRealVar c2("c2", "c2", 0, -1, 1) ;
    RooChebychev bkg("background", "background", fit_var, RooArgList(c0,c1,c2));

    RooRealVar nsig("nsig", "nsig", nsig_max/2, 0., nsig_max*1.1);
    RooRealVar nbkg("nbkg", "nbkg", nbkg_max/2, 0., nbkg_max*1.1);
    RooAddPdf model("model", "model", RooArgList(sig,bkg), RooArgList(nsig,nbkg));

    RooFitResult *fitResult = model.fitTo(*dataset, RooFit::Save(true), RooFit::Extended(true));
    fitResult->Print();

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 700);
    canvas_setup(canvas);
    canvas->Divide(1,2);
    canvas->cd(1);
    TVirtualPad* canvas1 = canvas->GetPad(1);
    canvas_setup_sub(canvas1);
    canvas1->SetPad(0,0.3,1,1);
    canvas1->SetBottomMargin(0.0);
    canvas1->SetRightMargin(0.05);
    canvas1->SetLogy(0);
    canvas1->SetFillColor(0);
    canvas1->SetFrameFillColor(0);

    RooPlot *frame = fit_var.frame();
    dataset->plotOn(frame, RooFit::Name("Candidates"), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(1.1), RooFit::Binning(100), RooFit::DrawOption("ep"));
    model.plotOn(frame, RooFit::Name("Voigtian"), RooFit::Components("signal"), RooFit::LineStyle(9), RooFit::LineColor(kOrange-3), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(frame, RooFit::Name("ChebyChev"), RooFit::Components("background"), RooFit::LineStyle(9), RooFit::LineColor(kBlue), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));
    model.plotOn(frame, RooFit::Name("Fit"), RooFit::Components("model"), RooFit::LineStyle(9), RooFit::LineColor(kRed), RooFit::LineWidth(2.0), RooFit::DrawOption("L"));

    frame->Draw("");
    frame->SetYTitle("# candidates");
	/* frame->GetYaxis()->SetTitleOffset(1.); */ 
	/* frame->GetYaxis()->CenterTitle(true); */
    frame->SetMinimum(0.01);
    frame->SetMaximum(h_temp->GetMaximum()*1.3);

    TLegend *legend = new TLegend(0.2, 0.2, 0.4, 0.4, NULL, "brNDC");
    legend->AddEntry(frame->findObject("Candidates"), "Candidates", "ep");
    legend->AddEntry(frame->findObject("Voigtian"), "Voigtian", "l");
    legend->AddEntry(frame->findObject("ChebyChev"), "ChebyChev", "l");
    legend->AddEntry(frame->findObject("Fit"), "Fit model", "l");
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.04);
    legend->Draw();
    
    write_text(0.6, 0.8, Form("Fit Results:"));
    write_text(0.6, 0.75, Form("#mu = %.3f #pm %.3f", mu.getVal(), mu.getError()));
    write_text(0.6, 0.7, Form("width = %.4f #pm %.4f", kuan.getVal(), kuan.getError()));
    write_text(0.6, 0.65, Form("#sigma = %.4f #pm %.4f", sigma.getVal(), sigma.getError()));
    write_text(0.6, 0.6, Form("# voigtian = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
    write_text(0.6, 0.55, Form("# chebychev = %.0f #pm %.0f", nbkg.getVal(), nbkg.getError()));
    CMS_lumi_sub(canvas1);

    canvas->cd(2);
    TVirtualPad* canvas2 = canvas->GetPad(2);
    canvas_setup_sub(canvas2);
    canvas2->SetPad(0,0,1,0.3);
    canvas2->SetRightMargin(0.05);
    canvas2->SetTopMargin(0.05);
    canvas2->SetBottomMargin(0.4);
    RooHist* hpull = frame->pullHist("Candidates","Fit");
    RooPlot* pull = new RooPlot("", "", fit_var, varmin, varmax, 100);
    hpull->SetFillColor(1);
    pull->addPlotable(hpull, "BX");
    pull->Draw();
    TAxis* xpull = pull->GetXaxis();
    TAxis* ypull = pull->GetYaxis();
    xpull->SetTitleOffset(3.2);
    xpull->SetTickLength(0.06);
    ypull->SetNdivisions(504);
    pull->SetXTitle(xtitle);
    /* pull->GetXaxis()->CenterTitle(true); */
    pull->GetXaxis()->SetTitleOffset(1.);
    pull->SetYTitle("Pull");
    pull->GetYaxis()->CenterTitle(true);
    pull->GetYaxis()->SetRangeUser(-5, 5);
    pull->GetYaxis()->SetTitleOffset(0.3);
    pull->GetXaxis()->SetTitleSize(0.1);
    pull->GetXaxis()->SetLabelSize(0.12); 
    pull->GetYaxis()->SetTitleSize(0.12);
    pull->GetYaxis()->SetLabelSize(0.12);
    canvas->Update();
    canvas->SaveAs(figpath);
}

int nonmatchphi() {
    setTDRStyle();

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

    /* fit_draw_nonmatch(nonmatch_tree, "phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", nonmatch_tree->GetEntries(), nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_phi_m1.png"); */
    /* fit_draw_nonmatch(nonmatch_tree, "phiFit_phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", nonmatch_tree->GetEntries(), nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_phi_m2.png"); */
    /* fit_draw_nonmatch(nonmatch_tree, "DsFit_phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", nonmatch_tree->GetEntries(), nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_phi_m3.png"); */
    fit_draw_match(match_tree, "phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_phi_m1.png");
    fit_draw_match(match_tree, "phiFit_phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_phi_m2.png");
    fit_draw_match(match_tree, "DsFit_phi_M", 0.99, 1.05, "#it{m}(#it{K^{+}K^{-}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_phi_m3.png");

    /* fit_draw_nonmatch(nonmatch_tree, "Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", 100, nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_Ds_m1.png"); */
    /* fit_draw_nonmatch(nonmatch_tree, "phiFit_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", 100, nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_Ds_m2.png"); */
    /* fit_draw_nonmatch(nonmatch_tree, "DsFit_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", 100, nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_Ds_m3.png"); */
    /* fit_draw_nonmatch(nonmatch_tree, "DsFit_Mconstraint_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", 100, nonmatch_tree->GetEntries(), "./figures/fit_compare/nonmatch_Ds_m4.png"); */
    fit_draw_match(match_tree, "Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_Ds_m1.png");
    fit_draw_match(match_tree, "phiFit_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_Ds_m2.png");
    fit_draw_match(match_tree, "DsFit_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_Ds_m3.png");
    fit_draw_match(match_tree, "DsFit_Mconstraint_Ds_M", 1.85, 2.1, "#it{m}(#it{K^{+}K^{-}#pi^{+}}) [GeV/#it{c}^{2}]", "./figures/fit_compare/match_Ds_m4.png");




    delete mychain;

    return 0;
}

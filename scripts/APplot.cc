#include <iostream>
#include <vector>
#include <TString.h>
#include <TChain.h>

#include "CMSplots/tdrStyle.c"
#include "CMSplots/CMS_lumi.c"
#include "CMSplots/draw_funcs.c"

const double m_Ds = 1.96835;
const double m_phi = 1.019460;
const double m_K = 0.493677;
const double m_pi = 0.13957039;
const double pst = sqrt((pow(m_Ds,4)+pow(m_phi,4)+pow(m_pi,4)-2*pow(m_Ds*m_phi,2)-2*pow(m_Ds*m_pi,2)-2*pow(m_phi*m_pi,2))/(4*pow(m_Ds,2)));

void draw_2d(TChain *mychain, TString myvar, TString vartitle, float xvarmin, float xvarmax, float yvarmin, float yvarmax, TCut mycut, TString figpath){

    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TH2F *hist = new TH2F("hist", vartitle, 100, xvarmin, xvarmax, 100, yvarmin, yvarmax);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas_setup(canvas);
    mychain->Project("hist", myvar, mycut);

	TBox *cutbox = new TBox(-0.4, 0, 0.4, 0.18);
	cutbox->SetFillColor(2);
	cutbox->SetFillStyle(3944);
	/* cutbox->SetLineWidth(0.05); */
	cutbox->Draw("same");
	gPad->Update();
	TLine *cutline1 = new TLine(-0.4, 0, -0.4, 0.18);
	TLine *cutline2 = new TLine(0.4, 0, 0.4, 0.18);
	TLine *cutline3 = new TLine(-0.4, 0, 0.4, 0);
	TLine *cutline4 = new TLine(-0.4, 0.18, 0.4, 0.18);
	cutline1->SetLineColor(2);
	cutline2->SetLineColor(2);
	cutline3->SetLineColor(2);
	cutline4->SetLineColor(2);
    cutline1->Draw("same");
    cutline2->Draw("same");
    cutline3->Draw("same");
    cutline4->Draw("same");

    /* hist->SetMaximum(1.3*hist->GetMaximum()); */
    hist->SetMinimum(0);
    hist->Draw("COLZ");
    CMS_lumi(canvas);
    canvas->Update();
    canvas->RedrawAxis();
    canvas->SaveAs("./figures/APplot/"+figpath+".png");
    delete hist;
}

int APplot() {

    setTDRStyle();
    lumi_sqrtS = "13 TeV, 2017 Simulation";

    TChain *mychain = new TChain("SelectionStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/SelectionStudy.root");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/SelectionStudy/output_1.root");

    if(gSystem->AccessPathName("./figures")) gSystem->MakeDirectory("./figures");
    if(gSystem->AccessPathName("./figures/APplot")) gSystem->MakeDirectory("./figures/APplot");

    draw_2d(mychain, "Gen_Kp_PP:((Gen_Kp_PL - Gen_Km_PL)/(Gen_Kp_PL + Gen_Km_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "1", "Gen_APplot_phi");
    draw_2d(mychain, "phiFit_Kp_PP:((phiFit_Kp_PL - phiFit_Km_PL)/(phiFit_Kp_PL + phiFit_Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "match_entry==1", "match_APplot_phi");
    draw_2d(mychain, "phiFit_Kp_PP:((phiFit_Kp_PL - phiFit_Km_PL)/(phiFit_Kp_PL + phiFit_Km_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 0.3, "non_match_entry==1", "nonmatch_APplot_phi");

    /* draw_2d(mychain, "Gen_phi_PP:((Gen_phi_PL - Gen_pi_PL)/(Gen_phi_PL + Gen_pi_PL))", ";Gen #it{#alpha};Gen #it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "1", "Gen_APplot_Ds"); */
    /* draw_2d(mychain, "DsFit_phi_PP:((DsFit_phi_PL - DsFit_pi_PL)/(DsFit_phi_PL + DsFit_pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "match_entry==1", "match_APplot_Ds"); */
    /* draw_2d(mychain, "DsFit_phi_PP:((DsFit_phi_PL - DsFit_pi_PL)/(DsFit_phi_PL + DsFit_pi_PL))", ";#it{#alpha};#it{p_{#perp}} [GeV]", -1, 1, 0, 1.5, "non_match_entry==1", "nonmatch_APplot_Ds"); */


    TH1F *h1 = new TH1F("h1", ";#it{p}^{*2}_{#perp}+#it{#alpha}^{2}#it{#beta}^{2}#it{E}^{*2} [GeV^{2}];# entries", 100, 0, 0.04);
    TH1F *h2 = new TH1F("h2", ";#it{p}^{*2}_{#perp}+#it{#alpha}^{2}#it{#beta}^{2}#it{E}^{*2} [GeV^{2}];# entries", 100, 0, 0.04);
    TH1F *h3 = new TH1F("h3", ";#it{p}^{*2}_{#perp}+[-#it{#alpha#beta}(#it{E}_{1}*+#it{E}_{2}*)/2+#it{#beta}(#it{E}_{1}*-#it{E}_{2}*)/2]^{2} [GeV^{2}];# entries", 100, 0, 1);
    TH1F *h4 = new TH1F("h4", ";#it{p}^{*2}_{#perp}+[-#it{#alpha#beta}(#it{E}_{1}*+#it{E}_{2}*)/2+#it{#beta}(#it{E}_{1}*-#it{E}_{2}*)/2]^{2} [GeV^{2}];# entries", 100, 0, 1);


    std::vector<bool> *match_entry_vec = nullptr;
    std::vector<bool> *non_match_entry_vec = nullptr;

    std::vector<double> *phiFit_Kp_PP_vec = nullptr;
    std::vector<double> *phiFit_Kp_PL_vec = nullptr;
    std::vector<double> *phiFit_Km_PP_vec = nullptr;
    std::vector<double> *phiFit_Km_PL_vec = nullptr;
    std::vector<double> *phiFit_phi_P_vec = nullptr;

    std::vector<double> *DsFit_phi_PP_vec = nullptr;
    std::vector<double> *DsFit_phi_PL_vec = nullptr;
    std::vector<double> *DsFit_pi_PP_vec = nullptr;
    std::vector<double> *DsFit_pi_PL_vec = nullptr;
    std::vector<double> *DsFit_Ds_P_vec = nullptr;

    mychain->SetBranchAddress("match_entry", &match_entry_vec);
    mychain->SetBranchAddress("non_match_entry", &non_match_entry_vec);

    mychain->SetBranchAddress("phiFit_Kp_PP", &phiFit_Kp_PP_vec);
    mychain->SetBranchAddress("phiFit_Kp_PL", &phiFit_Kp_PL_vec);
    mychain->SetBranchAddress("phiFit_Km_PP", &phiFit_Km_PP_vec);
    mychain->SetBranchAddress("phiFit_Km_PL", &phiFit_Km_PL_vec);
    mychain->SetBranchAddress("phiFit_phi_P", &phiFit_phi_P_vec);

    mychain->SetBranchAddress("DsFit_phi_PP", &DsFit_phi_PP_vec);
    mychain->SetBranchAddress("DsFit_phi_PL", &DsFit_phi_PL_vec);
    mychain->SetBranchAddress("DsFit_pi_PP", &DsFit_pi_PP_vec);
    mychain->SetBranchAddress("DsFit_pi_PL", &DsFit_pi_PL_vec);
    mychain->SetBranchAddress("DsFit_Ds_P", &DsFit_Ds_P_vec);
    
    int nentries = mychain->GetEntries();

    for(int i=0; i<nentries; i++){
        mychain->GetEntry(i);

        double phialpha;
        double phibeta;
        double Dsalpha;
        double Dsbeta;

        int vecsize = match_entry_vec->size();

        for(int j=0; j<vecsize; j++){

            phialpha = (phiFit_Kp_PL_vec->at(j)-phiFit_Km_PL_vec->at(j))/(phiFit_Kp_PL_vec->at(j)+phiFit_Km_PL_vec->at(j));
            phibeta = sqrt(pow(phiFit_phi_P_vec->at(j),2) / (pow(phiFit_phi_P_vec->at(j),2) + pow(m_phi,2)));

            double E1st = sqrt(pow(m_phi,2)+pow(pst,2));
            double E2st = sqrt(pow(m_pi,2)+pow(pst,2));

            Dsalpha = (DsFit_phi_PL_vec->at(j) - DsFit_pi_PL_vec->at(j))/(DsFit_phi_PL_vec->at(j) + DsFit_pi_PL_vec->at(j));
            Dsbeta = sqrt(pow(DsFit_Ds_P_vec->at(j),2) / (pow(DsFit_Ds_P_vec->at(j),2) + pow(m_Ds,2)));

            if(match_entry_vec->at(j)==1){
                h1->Fill(pow(phiFit_Kp_PP_vec->at(j),2)+pow(phialpha*phibeta,2)*pow(m_phi/2,2));
                h3->Fill(pow(DsFit_phi_PP_vec->at(j),2)+pow(Dsbeta*(E1st-E2st)/2-Dsalpha*Dsbeta*(E1st+E2st)/2,2));
            }
            if(non_match_entry_vec->at(j)==1){
                h2->Fill(pow(phiFit_Kp_PP_vec->at(j),2)+pow(phialpha*phibeta,2)*pow(m_phi/2,2));
                h4->Fill(pow(DsFit_phi_PP_vec->at(j),2)+pow(Dsbeta*(E1st-E2st)/2-Dsalpha*Dsbeta*(E1st+E2st)/2,2));
            }
        }

    }

    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
    h3->Scale(1./h3->Integral());
    h4->Scale(1./h4->Integral());

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvas_setup(c1);
    h1->SetMaximum(1.3*std::max(h1->GetMaximum(),h2->GetMaximum()));
    h1->SetMinimum(0);
    h1->SetMarkerColor(kBlue);
    gPad->SetBottomMargin(0.15);
    h1->Draw("ep");
    h2->SetMarkerColor(kRed);
    h2->Draw("ep same");
    TLegend * mylegend1 = new TLegend(0.7, 0.75, 0.9, 0.9);
    mylegend1->AddEntry(h1, "match", "ep");
    mylegend1->AddEntry(h2, "nonmatch", "ep");
    mylegend1->SetFillColor(0);
    mylegend1->SetLineWidth(0);
    mylegend1->SetTextSize(0.04);
    mylegend1->Draw();
    CMS_lumi(c1);
    c1->Update();
    c1->RedrawAxis();
    c1->SaveAs("./figures/APplot/phi_pst2.png");
    
    /* TCanvas *c2 = new TCanvas("c2", "c2", 800, 600); */
    /* canvas_setup(c2); */
    /* h2->SetMaximum(1.3*h2->GetMaximum()); */
    /* h2->SetMinimum(0); */
    /* h2->Draw("ep"); */
    /* CMS_lumi(c2); */
    /* c2->Update(); */
    /* c2->RedrawAxis(); */
    
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    canvas_setup(c3);
    h3->SetMaximum(1.3*std::max(h3->GetMaximum(),h4->GetMaximum()));
    h3->SetMinimum(0);
    h3->SetMarkerColor(kBlue);
    gPad->SetBottomMargin(0.15);
    h3->Draw("ep");
    h4->SetMarkerColor(kRed);
    h4->Draw("ep same");
    TLegend * mylegend2 = new TLegend(0.7, 0.75, 0.9, 0.9);
    mylegend2->AddEntry(h3, "match", "ep");
    mylegend2->AddEntry(h4, "nonmatch", "ep");
    mylegend2->SetFillColor(0);
    mylegend2->SetLineWidth(0);
    mylegend2->SetTextSize(0.04);
    mylegend2->Draw();
    CMS_lumi(c3);
    c3->Update();
    c3->RedrawAxis();
    c3->SaveAs("./figures/APplot/Ds_pst2.png");
    
    /* TCanvas *c4 = new TCanvas("c4", "c4", 800, 600); */
    /* canvas_setup(c4); */
    /* h4->SetMaximum(1.3*h4->GetMaximum()); */
    /* h4->SetMinimum(0); */
    /* h4->Draw("ep"); */
    /* CMS_lumi(c4); */
    /* c4->Update(); */
    /* c4->RedrawAxis(); */

    delete mychain;

    return 0;
}

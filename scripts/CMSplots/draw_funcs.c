#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>

void canvas_setup(TPad* pad) {
    pad->SetFillColor(0);
    pad->SetBorderMode(0);
    pad->SetFrameFillStyle(0);
    pad->SetFrameBorderMode(0);
    pad->SetLeftMargin( 0.16 );
    pad->SetRightMargin( 0.12 );
    pad->SetTopMargin( 0.08 );
    pad->SetBottomMargin( 0.12 );
    pad->SetTickx(0);
    pad->SetTicky(0);
}

void write_text(float pos_x, float pos_y, TString the_text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(12);
    latex.SetTextSize(0.04);
    latex.DrawLatex(pos_x,pos_y,the_text);
}

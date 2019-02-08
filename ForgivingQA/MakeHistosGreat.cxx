/*
 * MakeHistosGreat.cxx
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#include "MakeHistosGreat.h"
#include "TCanvas.h"
#include "TStyle.h"
std::vector<int> fFillColors = { kGray + 1, kRed - 10, kBlue - 9, kGreen - 8,
    kMagenta - 9, kOrange - 9, kCyan - 8, kYellow - 7 };
std::vector<int> fColors = { kBlack, kRed + 1, kBlue + 2, kGreen + 3, kMagenta
    + 1, kOrange - 1, kCyan + 2, kYellow + 2 };
std::vector<int> fMarkers = { kFullCircle, kFullSquare, kOpenCircle,
    kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar,
    kOpenStar };

MakeHistosGreat::MakeHistosGreat() {
  // TODO Auto-generated constructor stub

}

MakeHistosGreat::~MakeHistosGreat() {
  // TODO Auto-generated destructor stub
}

void MakeHistosGreat::FormatHistogram(TH1* hist, unsigned int marker,
                                      unsigned int color) {
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetLabelOffset(0.01);

  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.25);

  hist->SetMarkerStyle(fMarkers[marker]);
  hist->SetMarkerColor(fColors[color]);
  hist->SetLineColor(fColors[color]);

  hist->SetLineWidth(2);
}

void MakeHistosGreat::DrawAndStore(TH1* hist, const char* outname,
                                   const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  c1->cd();
  hist->Draw(drawOption);
  c1->SaveAs(Form("%s.pdf", outname));
  delete c1;
  return;
}

void MakeHistosGreat::DrawLogYAndStore(TH1* hist, const char* outname,
                                   const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  c1->cd();
  c1->SetLogy();
  hist->Draw(drawOption);
  c1->SaveAs(Form("%s.pdf", outname));
  delete c1;
  return;
}

void MakeHistosGreat::DrawAndStore(TH2* hist, const char* outname,
                                   const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  c1->cd();
  hist->Draw(drawOption);
  c1->SaveAs(Form("%s.pdf", outname));
  delete c1;
  return;
}

void MakeHistosGreat::SetStyle(bool graypalette, bool title) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if (graypalette)
    gStyle->SetPalette(8, 0);
  else
    gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetLabelOffset(0.01, "x");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
}

/*
 * MakeHistosGreat.cxx
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#include "MakeHistosGreat.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
std::vector<int> fFillColors = { kGray + 1, kRed - 10, kBlue - 9, kGreen - 8,
    kMagenta - 9, kOrange - 9, kCyan - 8, kYellow - 7 };
std::vector<int> fColors = { kBlack, kRed + 1, kBlue + 2, kGreen + 3, kMagenta
    + 1, kOrange - 1, kCyan + 2, kYellow + 2, kBlue + 3 };
std::vector<int> fMarkers = { kFullCircle, kFullSquare, kOpenCircle,
    kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar,
    kOpenStar };

MakeHistosGreat::MakeHistosGreat()
    : fTightMargin(false),
      fOutfile(nullptr) {
}

MakeHistosGreat::MakeHistosGreat(const char* outname)
    : fTightMargin(false),
      fOutfile(nullptr) {
  fOutfile = new TFile(Form("%s.root", outname), "RECREATE");
}

MakeHistosGreat::~MakeHistosGreat() {
  delete fOutfile;
  // TODO Auto-generated destructor stub
}

void MakeHistosGreat::FormatHistogram(TH1* hist, unsigned int marker,
                                      unsigned int color, float size) {
  hist->GetXaxis()->SetLabelSize(28);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetLabelSize(28);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetTitleFont(43);

  hist->GetXaxis()->SetTitleSize(28);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetTitleSize(28);
  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetTitleFont(43);

  hist->SetMarkerStyle(fMarkers[marker]);
  hist->SetMarkerColor(fColors[color]);
  hist->SetMarkerSize(size);
  hist->SetLineColor(fColors[color]);
  hist->SetLineWidth(2);
}

void MakeHistosGreat::FormatSmallHistogram(TH1* hist, unsigned int marker,
                                           unsigned int color, float size) {
  hist->GetXaxis()->SetLabelSize(15);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetTitleSize(15);
  hist->GetXaxis()->SetTitleOffset(2.2);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetTitleFont(43);

  hist->GetYaxis()->SetMaxDigits(2);
  hist->GetYaxis()->SetLabelSize(15);
  hist->GetYaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetTitleSize(15);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetTitleFont(43);

  hist->SetMarkerStyle(fMarkers[marker]);
  hist->SetMarkerColor(fColors[color]);
  hist->SetMarkerSize(size);
  hist->SetLineColor(fColors[color]);
  hist->SetLineWidth(2);
}

void MakeHistosGreat::FormatHistogram(TH2 *histo) {
  histo->GetXaxis()->SetLabelSize(28);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetLabelSize(28);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->GetXaxis()->SetTitleSize(28);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleSize(28);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
}

void MakeHistosGreat::DrawAndStore(std::vector<TH1*> hist, const char* outname,
                                   const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  if (fTightMargin) {
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.01);
  }
  c1->cd();
  TString DrawOpt = Form("%s", drawOption);
  bool oneTime = false;
  for (auto it : hist) {
    it->Draw(DrawOpt.Data());
    if (!oneTime) {
      oneTime = true;
      DrawOpt += "same";
    }
  }
  c1->SaveAs(Form("%s.pdf", outname));

  if(fOutfile) {
    fOutfile->cd();
    for (auto it : hist) {
      it->Write();
    }
    c1->Write(outname);
  }

  delete c1;
  return;
}

void MakeHistosGreat::DrawOnPad(std::vector<TH1*> hist, TPad* TPain,
                                const char* drawOption) {
  TPain->cd();
  TString DrawOpt = Form("%s", drawOption);
  bool oneTime = false;
  for (auto it : hist) {
    it->Draw(DrawOpt.Data());
    if (!oneTime) {
      oneTime = true;
      DrawOpt += "same";
    }
  }
  return;
}

void MakeHistosGreat::DrawLogYAndStore(std::vector<TH1*> hist,
                                       const char* outname,
                                       const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  if (fTightMargin) {
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.01);
  }
  c1->cd();
  c1->SetLogy();
  TString DrawOpt = Form("%s", drawOption);
  bool oneTime = false;
  for (auto it : hist) {
    it->Draw(DrawOpt.Data());
    if (!oneTime) {
      oneTime = true;
      DrawOpt += "same";
    }
  }
  c1->SaveAs(Form("%s.pdf", outname));

  if(fOutfile) {
    fOutfile->cd();
    for (auto it : hist) {
      it->Write();
    }
    c1->Write(outname);
  }

  delete c1;
  return;
}

void MakeHistosGreat::DrawAndStore(std::vector<TH2*> hist, const char* outname,
                                   const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  if (fTightMargin) {
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.01);
  } else {
    c1->SetRightMargin(0.14);
  }
  c1->cd();
  TString DrawOpt = Form("%s", drawOption);
  bool oneTime = false;
  for (auto it : hist) {
    it->Draw(DrawOpt.Data());
    it->GetZaxis()->SetTitleFont(43);
    it->GetZaxis()->SetTitleSize(28);
    it->GetZaxis()->SetLabelFont(43);
    it->GetZaxis()->SetLabelSize(28);
    if (!oneTime) {
      oneTime = true;
      DrawOpt += "same";
    }
  }
  c1->SaveAs(Form("%s.pdf", outname));

  if(fOutfile) {
    fOutfile->cd();
    for (auto it : hist) {
      it->Write();
    }
    c1->Write(outname);
  }

  delete c1;
  return;
}

void MakeHistosGreat::DrawLogZAndStore(std::vector<TH2*> hist,
                                       const char* outname,
                                       const char* drawOption) {
  auto c1 = new TCanvas(Form("%s", outname), Form("%s", outname));
  if (fTightMargin) {
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.01);
  } else {
    c1->SetRightMargin(0.14);
  }
  c1->cd();
  c1->SetLogz();
  TString DrawOpt = Form("%s", drawOption);
  bool oneTime = false;
  for (auto it : hist) {
    it->Draw(DrawOpt.Data());
    it->GetZaxis()->SetTitleFont(43);
    it->GetZaxis()->SetTitleSize(28);
    it->GetZaxis()->SetLabelFont(43);
    it->GetZaxis()->SetLabelSize(28);
    if (!oneTime) {
      oneTime = true;
      DrawOpt += "same";
    }
  }
  c1->SaveAs(Form("%s.pdf", outname));

  if(fOutfile) {
    fOutfile->cd();
    for (auto it : hist) {
      it->Write();
    }
    c1->Write(outname);
  }

  delete c1;
  return;
}

void MakeHistosGreat::SetStyle(bool title) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.075);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetLabelSize(28, "xyz");
  gStyle->SetTitleSize(28, "xyz");
  gStyle->SetLabelOffset(0.01, "xy");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.25, "x");
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);
  gStyle->SetPalette(kCividis);
}

void MakeHistosGreat::DrawLatexLabel(float pTMin, float pTMax,
                                     ForgivingFitter* fit, TPad* pad,
                                     const char* part, float xPos, float yPos, float offset) {
  const int signal = fit->GetSignalCounts();
  const float purity = fit->GetPurity();
  const float meanMass = fit->GetMeanMass();
  const float meanWidthActual = fit->GetMeanWidth();
  pad->cd();
  TLatex Label;
  Label.SetNDC(kTRUE);
  Label.SetTextSize(gStyle->GetTextSize() * 1.2);
  int counter = 0;

  float signalC = signal;
  int decimal = 0;
  while ( int(signalC) % 1000 != 0 ) {
    signalC /= 1000;
    decimal += 3;
  }
  signalC *= 1000;
  decimal -= 3;

  Label.DrawLatex(gPad->GetUxmax() - xPos,
                  gPad->GetUymax() - yPos + offset * counter++,
                  Form("Purity = %.1f %%", purity * 100.f));
  Label.DrawLatex(
      gPad->GetUxmax() - xPos,
      gPad->GetUymax() - yPos + offset * counter++,
      Form("#sigma_{%s} = %.1f MeV/#it{c}^{2}", part,
           meanWidthActual * 1000.f));
  Label.DrawLatex(
      gPad->GetUxmax() - xPos, gPad->GetUymax() - yPos + offset * counter++,
      Form("#LT#it{M}_{%s}#GT = %.1f MeV/#it{c}^{2}", part, meanMass * 1000.f));
  Label.DrawLatex(gPad->GetUxmax() - xPos,
                  gPad->GetUymax() - yPos + offset * counter++,
                  Form("%s: %.1f #times 10^{%d}", part, signalC, decimal));
  Label.DrawLatex(gPad->GetUxmax() - xPos,
                  gPad->GetUymax() - yPos + offset * counter++,
                  Form("%.2f < #it{p}_{T} < %.2f GeV/#it{c}", pTMin, pTMax));
}

void MakeHistosGreat::DrawPerformance(ForgivingFitter* fit, TPad* pad,
                                      const char* part, float xPos, float yPos,
                                      float pTmin, float pTmax) {
  float signal = (float) fit->GetSignalCounts();
  float background = (float) fit->GetBackgroundCounts();
  pad->cd();
  const float offset = 0.06;
  float counter = 0;
  TLatex Label;
  Label.SetNDC(kTRUE);
  Label.SetTextSize(gStyle->GetTextSize() * 0.9);
  Label.DrawLatex(xPos, yPos - offset * counter++, "#bf{ALICE Performance}");
  Label.DrawLatex(xPos, yPos - offset * counter++, "pp #sqrt{s} = 13 TeV");
  Label.DrawLatex(xPos, yPos - offset * counter++, "High Mult. (0-0.072% INEL)");
  if (pTmin > 0 && pTmax > 0) {
    Label.DrawLatex(xPos, yPos - offset * counter++,
                    Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTmin, pTmax));
  }
  Label.DrawLatex(xPos, yPos - offset * counter++,
                  Form("%s: %.0f", part, signal));
  Label.DrawLatex(
      xPos, yPos - offset * counter++,
      Form("Purity = %.1f %%", signal / (signal + background) * 100.f));
}

void MakeHistosGreat::DrawPublication(ForgivingFitter* fit, TPad* pad,
                                      const char* part, float xPos, float yPos,
                                      float pTmin, float pTmax) {
  TH1F* hist = nullptr;
  if(pad->GetListOfPrimitives()->At(0)) {
    hist = (TH1F*)pad->GetListOfPrimitives()->At(0);
  }

  TString partString = TString::Format("%s", part);
  float signal = (float) fit->GetSignalCounts();
  float background = (float) fit->GetBackgroundCounts();
  pad->cd();
  const float offset = 0.065;
  float counter = 0;
  TLatex Label;
  Label.SetNDC(kTRUE);
  Label.SetTextSize(gStyle->GetTextSize() * 0.9);
  Label.DrawLatex(xPos, yPos - offset * counter++, "ALICE pp #sqrt{#it{s}} = 13 TeV");
  Label.DrawLatex(xPos, yPos - offset * counter++, "High-mult. (0#kern[-0.95]{ }#minus#kern[-0.05]{ }0.072#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
  if (pTmin > 0 && pTmax > 0) {
    Label.DrawLatex(xPos, yPos - offset * counter++,
                    Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", pTmin, pTmax));
  }
  if(partString == TString("#Sigma^{0} + #bar{#Sigma^{0}}")) {
    Label.DrawLatex(xPos, yPos - offset * counter++, "#Sigma^{0} #rightarrow #Lambda#gamma, #bar{#Sigma^{0}} #rightarrow #bar{#Lambda}#gamma");
  } else if(partString == TString("#bar{#Sigma^{0}}")) {
    Label.DrawLatex(xPos, yPos - offset * counter++, "#bar{#Sigma^{0}} #rightarrow #bar{#Lambda}#gamma");
  } else if(partString == TString("#Sigma^{0}")) {
    Label.DrawLatex(xPos, yPos - offset * counter++, "#Sigma^{0} #rightarrow #Lambda#gamma");
  }
  auto leg = new TLegend(xPos - 0.005, yPos - offset * (counter - 0.75), xPos + 0.2, yPos - offset * (counter + 2));
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg->AddEntry(hist, "Data", "P");
  leg->AddEntry(fit->GetFullFitFunction(), "Total fit", "l");
  leg->AddEntry(fit->GetBackgroundFunction(), "Background", "l");
  leg->Draw("same");
}

void MakeHistosGreat::DrawLine(TPad* pad, float xMin, float xMax, float yMin, float yMax, int color) {
  TLine one;
  one.SetLineColor(color);
  one.SetLineWidth(2);
  one.SetLineStyle(3);
  pad->cd();
  one.DrawLine(xMin,yMin,xMax,yMax);
}

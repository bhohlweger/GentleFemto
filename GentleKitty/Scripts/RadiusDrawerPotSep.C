#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "DreamPlot.h"
#include "TFile.h"
#include "TDatabasePDG.h"

int main(int argc, char* argv[]) {
  if(!argv[1]) {
    std::cout << "pp RadFile missing\n";
    return -1;
  }
  if(!argv[2]) {
    std::cout << "pL NLO RadFile missing\n";
    return -1;
  }
  if(!argv[3]) {
    std::cout << "pL LO RadFile missing\n";
    return -1;
  }
  if(!argv[4]) {
    std::cout << "Source Name\n";
    return -1;
  }
  const char* ppFile = argv[1];
  const char* pLNLOFile = argv[2];
  const char* pLLOFile = argv[3];
  const char* sourceName = argv[4];
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile* ppHMFile =
      TFile::Open(
          ppFile,
          "read");
  TGraphErrors* mTppHMSys = (TGraphErrors*) ppHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTppHMStat = (TGraphErrors*) ppHMFile->Get("mTRadiusStat");

  TFile* pLNLOHMFile =
      TFile::Open(
          pLNLOFile,
          "read");
  TGraphErrors* mTpLNLOHMSys = (TGraphErrors*) pLNLOHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLNLOHMStat = (TGraphErrors*) pLNLOHMFile->Get("mTRadiusStat");

  TFile* pLLOHMFile =
      TFile::Open(
          pLLOFile,
          "read");
  TGraphErrors* mTpLLOHMSys = (TGraphErrors*) pLLOHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLLOHMStat = (TGraphErrors*) pLLOHMFile->Get("mTRadiusStat");

  double yMin = 1234567;
  double yMax = 0;
  double x,y;
  for (int iBin = 0; iBin < mTppHMSys->GetN(); iBin++) {
    mTppHMSys->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    mTppHMSys->SetPointError(iBin, 0.4 * mTppHMSys->GetErrorX(iBin),
                             mTppHMSys->GetErrorY(iBin));
  }

  for (int iBin = 0; iBin < mTpLNLOHMSys->GetN(); iBin++) {
    mTpLNLOHMSys->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    mTpLNLOHMSys->SetPointError(iBin, 0.4 * mTpLNLOHMSys->GetErrorX(iBin),
                             mTpLNLOHMSys->GetErrorY(iBin));
  }

  for (int iBin = 0; iBin < mTpLLOHMSys->GetN(); iBin++) {
    mTpLLOHMSys->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    mTpLLOHMSys->SetPointError(iBin, 0.4 * mTpLLOHMSys->GetErrorX(iBin),
                             mTpLLOHMSys->GetErrorY(iBin));
  }
  TFile* out = TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  TLegend* leg = new TLegend(0.55, 0.46, 0.826, 0.6);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetNColumns(2);
  leg->SetTextSizePixels(40);

  mTpLNLOHMSys->SetPoint(mTpLNLOHMSys->GetN(), 0.95, 0.5);
  mTpLNLOHMSys->SetPointError(mTpLNLOHMSys->GetN(), 0., 0.);

  mTpLNLOHMSys->SetPoint(mTpLNLOHMSys->GetN(), 2.7, 1.3);
  mTpLNLOHMSys->SetPointError(mTpLNLOHMSys->GetN(), 0., 0.);

  mTpLLOHMSys->SetPoint(mTpLLOHMSys->GetN(), 0.95, 0.5);
  mTpLLOHMSys->SetPointError(mTpLLOHMSys->GetN(), 0., 0.);

  mTpLLOHMSys->SetPoint(mTpLLOHMSys->GetN(), 2.7, 1.3);
  mTpLLOHMSys->SetPointError(mTpLLOHMSys->GetN(), 0., 0.);

  mTpLNLOHMSys->SetTitle(
      "; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); #it{r}_{core} (fm)");
  mTpLNLOHMSys->GetXaxis()->SetTitleSize(40);
  mTpLNLOHMSys->GetYaxis()->SetTitleSize(40);
  mTpLNLOHMSys->GetXaxis()->SetTitleOffset(1.35);
  mTpLNLOHMSys->GetYaxis()->SetTitleOffset(1.4);
  mTpLNLOHMSys->GetXaxis()->SetLabelSize(40);
  mTpLNLOHMSys->GetYaxis()->SetLabelSize(40);
  mTpLNLOHMSys->GetXaxis()->SetLabelOffset(.02);
  mTpLNLOHMSys->GetYaxis()->SetLabelOffset(.02);
  mTpLNLOHMSys->GetXaxis()->SetRangeUser(0.95, 2.7);
//  mTpLHMSys->GetYaxis()->SetRangeUser(0.8, 2.0);
  mTpLNLOHMSys->GetYaxis()->SetRangeUser(0.8*yMin, 1.2*yMax);

  mTpLNLOHMSys->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLNLOHMSys->SetFillStyle(3225);
  mTpLNLOHMSys->Draw("2AP");

  mTpLNLOHMStat->SetMarkerColor(kRed + 1);
  mTpLNLOHMStat->SetLineColor(kRed + 1);
  mTpLNLOHMStat->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLNLOHMStat->SetLineWidth(1);
  mTpLNLOHMStat->SetFillStyle(3225);
  mTpLNLOHMStat->SetMarkerStyle(41);
  mTpLNLOHMStat->SetMarkerSize(2.0);
  mTpLNLOHMStat->Draw("pez same");

  mTpLLOHMSys->SetFillColorAlpha(kGreen - 7, 0.7);
  mTpLLOHMSys->SetFillStyle(3225);
  mTpLLOHMSys->Draw("2PSAME");

  mTpLLOHMStat->SetMarkerColor(kGreen + 2);
  mTpLLOHMStat->SetLineColor(kGreen + 2);
  mTpLLOHMStat->SetFillColorAlpha(kGreen - 7, 0.7);
  mTpLLOHMStat->SetLineWidth(1);
  mTpLLOHMStat->SetFillStyle(3225);
  mTpLLOHMStat->SetMarkerStyle(34);
  mTpLLOHMStat->SetMarkerSize(2.0);
  mTpLLOHMStat->Draw("pez same");

  mTppHMSys->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppHMSys->SetFillStyle(3225);
  mTppHMSys->Draw("2PSAME");

  mTppHMStat->SetMarkerColor(kBlue + 2);
  mTppHMStat->SetLineColor(kBlue + 2);
  mTppHMStat->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppHMStat->SetLineWidth(1);
  mTppHMStat->SetFillStyle(3225);
  mTppHMStat->SetMarkerStyle(34);
  mTppHMStat->SetMarkerSize(2.0);
  mTppHMStat->Draw("pez same");

//  leg->AddEntry((TObject*)0, "        ", "");
//  leg->AddEntry((TObject*)0, " MB ", "");
//  leg->AddEntry((TObject*)0, " HM ", "");
//  leg->AddEntry((TObject*)0, "p#minus#kern[0.4]{p}", "");
//  leg->AddEntry(mTppMBStat, "    ", "pef");
//  leg->AddEntry(mTppHMStat, "    ", "pef");
//  leg->AddEntry((TObject*)0, "p#minus#kern[0.2]{#Lambda}", "");
//  leg->AddEntry(mTpLMBStat, "    ", "pef");
//  leg->AddEntry(mTpLHMStat, "    ", "pef");
  leg->AddEntry(mTppHMStat, "p#minus#kern[0.4]{p}", "pef");
  leg->AddEntry(mTpLNLOHMStat, "p#minus#kern[0.2]{#Lambda} (NLO)", "pef");
  leg->AddEntry(mTpLLOHMStat, "p#minus#kern[0.2]{#Lambda} (LO)", "pef");

  TLatex BeamText;
  BeamText.SetTextFont(43);
  BeamText.SetTextSize(40);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.55, 0.83,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.55, 0.76, "High-mult.");
  BeamText.DrawLatex(
      0.55,
      0.69,
      "(0#kern[-0.95]{ }#minus#kern[-0.05]{ }0.072#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
//  BeamText.DrawLatex(0.55, 0.62, "Gaussian Source");
  BeamText.DrawLatex(0.17, 0.2, TString::Format("%s",sourceName).Data());

  leg->Draw("same");
  c4->SaveAs(Form("%s/mTvsRad.pdf", gSystem->pwd()));
  c4->Write();
  out->Write();
  out->Close();
  ppHMFile->Close();
  pLNLOHMFile->Close();
  pLLOHMFile->Close();
}


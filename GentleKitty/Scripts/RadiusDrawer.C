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
    std::cout << "pL RadFile missing\n";
    return -1;
  }
  if(!argv[3]) {
    std::cout << "Source Name\n";
    return -1;
  }
  const char* ppFile = argv[1];
  const char* pLFile = argv[2];
  const char* sourceName = argv[3];
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile* ppHMFile =
      TFile::Open(
          ppFile,
          "read");
  TGraphErrors* mTppHMSys = (TGraphErrors*) ppHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTppHMStat = (TGraphErrors*) ppHMFile->Get("mTRadiusStat");

  TFile* pLHMFile =
      TFile::Open(
          pLFile,
          "read");
  TGraphErrors* mTpLHMSys = (TGraphErrors*) pLHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLHMStat = (TGraphErrors*) pLHMFile->Get("mTRadiusStat");
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

  for (int iBin = 0; iBin < mTpLHMSys->GetN(); iBin++) {
    mTpLHMSys->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    mTpLHMSys->SetPointError(iBin, 0.4 * mTpLHMSys->GetErrorX(iBin),
                             mTpLHMSys->GetErrorY(iBin));
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

  mTpLHMSys->SetPoint(mTpLHMSys->GetN(), 0.95, 0.5);
  mTpLHMSys->SetPointError(mTpLHMSys->GetN(), 0., 0.);

  mTpLHMSys->SetPoint(mTpLHMSys->GetN(), 2.7, 1.3);
  mTpLHMSys->SetPointError(mTpLHMSys->GetN(), 0., 0.);

  mTpLHMSys->SetTitle(
      "; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); #it{r}_{core} (fm)");
  mTpLHMSys->GetXaxis()->SetTitleSize(40);
  mTpLHMSys->GetYaxis()->SetTitleSize(40);
  mTpLHMSys->GetXaxis()->SetTitleOffset(1.35);
  mTpLHMSys->GetYaxis()->SetTitleOffset(1.4);
  mTpLHMSys->GetXaxis()->SetLabelSize(40);
  mTpLHMSys->GetYaxis()->SetLabelSize(40);
  mTpLHMSys->GetXaxis()->SetLabelOffset(.02);
  mTpLHMSys->GetYaxis()->SetLabelOffset(.02);
  mTpLHMSys->GetXaxis()->SetRangeUser(0.95, 2.7);
//  mTpLHMSys->GetYaxis()->SetRangeUser(0.8, 2.0);
  mTpLHMSys->GetYaxis()->SetRangeUser(0.8*yMin, 1.2*yMax);

  mTpLHMSys->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLHMSys->SetFillStyle(3225);
  mTpLHMSys->Draw("2AP");

  mTpLHMStat->SetMarkerColor(kRed + 1);
  mTpLHMStat->SetLineColor(kRed + 1);
  mTpLHMStat->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLHMStat->SetLineWidth(1);
  mTpLHMStat->SetFillStyle(3225);
  mTpLHMStat->SetMarkerStyle(41);
  mTpLHMStat->SetMarkerSize(2.0);
  mTpLHMStat->Draw("pez same");

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
  leg->AddEntry(mTpLHMStat, "p#minus#kern[0.2]{#Lambda}", "pef");

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
      "(0#kern[-0.95]{ }#minus#kern[-0.05]{ }0.17#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
//  BeamText.DrawLatex(0.55, 0.62, "Gaussian Source");
  BeamText.DrawLatex(0.17, 0.2, TString::Format("%s",sourceName).Data());

  leg->Draw("same");
  c4->SaveAs(Form("%s/mTvsRad.pdf", gSystem->pwd()));
  c4->Write();
  out->Write();
  out->Close();
  ppHMFile->Close();
  pLHMFile->Close();
}


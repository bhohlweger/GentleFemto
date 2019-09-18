#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "DreamPlot.h"
#include "TFile.h"

int main(int argc, char* argv[]) {
  const char* sourceName = "Gauss";
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.9);
  TFile* ppMBFile =
      TFile::Open(
          "~/cernbox/ppMinimumBias/180824_FullSys/Fits/pp/180828_Gauss/RadppvsmT.root",
          "read");
  TGraphErrors* mTppMBSys = (TGraphErrors*) ppMBFile->Get("mTRadiusSyst");
  TGraphErrors* mTppMBStat = (TGraphErrors*) ppMBFile->Get("mTRadiusStat");

  TFile* ppHMFile =
      TFile::Open(
          "~/cernbox/HM13TeV/AnalysisData/RandomSystematics_mT/fits/pp/190719_mTGauss/RadppvsmT.root",
          "read");
  TGraphErrors* mTppHMSys = (TGraphErrors*) ppHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTppHMStat = (TGraphErrors*) ppHMFile->Get("mTRadiusStat");

  TFile* pLMBFile =
      TFile::Open(
          "~/cernbox/ppMinimumBias/180824_FullSys/Fits/pL/180828_Gauss/RadpLvsmT.root",
          "read");
  TGraphErrors* mTpLMBSys = (TGraphErrors*) pLMBFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLMBStat = (TGraphErrors*) pLMBFile->Get("mTRadiusStat");

  TFile* pLHMFile =
      TFile::Open(
          "~/cernbox/HM13TeV/AnalysisData/RandomSystematics_mT/fits/pL/190725_Gauss/RadpLvsmT.root",
          "read");
  TGraphErrors* mTpLHMSys = (TGraphErrors*) pLHMFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLHMStat = (TGraphErrors*) pLHMFile->Get("mTRadiusStat");
  for (int iBin = 0; iBin < mTppMBSys->GetN(); iBin++) {
    mTppMBSys->SetPointError(iBin, 0.4 * mTppMBSys->GetErrorX(iBin),
                             mTppMBSys->GetErrorY(iBin));
  }
  for (int iBin = 0; iBin < mTppHMSys->GetN(); iBin++) {
    mTppHMSys->SetPointError(iBin, 0.4 * mTppHMSys->GetErrorX(iBin),
                             mTppHMSys->GetErrorY(iBin));
  }
  for (int iBin = 0; iBin < mTpLMBSys->GetN(); iBin++) {
    mTpLMBSys->SetPointError(iBin, 0.4 * mTpLMBSys->GetErrorX(iBin),
                             mTpLMBSys->GetErrorY(iBin));
  }
  for (int iBin = 0; iBin < mTpLHMSys->GetN(); iBin++) {
    mTpLHMSys->SetPointError(iBin, 0.4 * mTpLHMSys->GetErrorX(iBin),
                             mTpLHMSys->GetErrorY(iBin));
  }
  TFile* out = TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  TLegend* leg = new TLegend(0.526, 0.48, 0.826, 0.68);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetNColumns(3);
  leg->SetTextSizePixels(40);

  mTpLHMSys->SetPoint(mTpLHMSys->GetN(), 0.95, 0.5);
  mTpLHMSys->SetPointError(mTpLHMSys->GetN(), 0., 0.);

  mTpLHMSys->SetPoint(mTpLHMSys->GetN(), 2.7, 1.3);
  mTpLHMSys->SetPointError(mTpLHMSys->GetN(), 0., 0.);

  mTpLHMSys->SetTitle("; < #it{m}_{T} >  (MeV/#it{c}^{2}); #it{r}_{Gauss} (fm)");
  mTpLHMSys->GetXaxis()->SetTitleSize(40);
  mTpLHMSys->GetYaxis()->SetTitleSize(40);
  mTpLHMSys->GetXaxis()->SetTitleOffset(1.4);
  mTpLHMSys->GetYaxis()->SetTitleOffset(1.4);
  mTpLHMSys->GetXaxis()->SetLabelSize(40);
  mTpLHMSys->GetYaxis()->SetLabelSize(40);
  mTpLHMSys->GetXaxis()->SetLabelOffset(.02);
  mTpLHMSys->GetYaxis()->SetLabelOffset(.02);
  mTpLHMSys->GetXaxis()->SetRangeUser(0.95, 2.7);
  mTpLHMSys->GetYaxis()->SetRangeUser(0.6, 2.3);

  mTpLHMSys->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLHMSys->SetFillStyle(3225);
  mTpLHMSys->Draw("2AP");

  mTpLHMStat->SetMarkerColor(kRed + 1);
  mTpLHMStat->SetLineColor(kRed + 1);
  mTpLHMStat->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLHMStat->SetLineWidth(1);
  mTpLHMStat->SetFillStyle(3225);
  mTpLHMStat->SetMarkerStyle(41);
  mTpLHMStat->SetMarkerSize(1.8);
  mTpLHMStat->Draw("pez same");

  mTpLMBSys->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLMBSys->SetFillStyle(3252);
  mTpLMBSys->Draw("2PSAME");

  mTpLMBStat->SetMarkerColor(kRed + 1);
  mTpLMBStat->SetLineColor(kRed + 1);
  mTpLMBStat->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLMBStat->SetLineWidth(1);
  mTpLMBStat->SetFillStyle(3252);
  mTpLMBStat->SetMarkerStyle(20);
  mTpLMBStat->SetMarkerSize(1.8);
  mTpLMBStat->Draw("pez same");

  mTppHMSys->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppHMSys->SetFillStyle(3225);
  mTppHMSys->Draw("2PSAME");

  mTppHMStat->SetMarkerColor(kBlue + 2);
  mTppHMStat->SetLineColor(kBlue + 2);
  mTppHMStat->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppHMStat->SetLineWidth(1);
  mTppHMStat->SetFillStyle(3225);
  mTppHMStat->SetMarkerStyle(34);
  mTppHMStat->SetMarkerSize(1.8);
  mTppHMStat->Draw("pez same");

  mTppMBSys->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppMBSys->SetFillStyle(3252);
  mTppMBSys->Draw("2PSAME");

  mTppMBStat->SetMarkerColor(kBlue + 2);
  mTppMBStat->SetLineColor(kBlue + 2);
  mTppMBStat->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppMBStat->SetLineWidth(1);
  mTppMBStat->SetFillStyle(3252);
  mTppMBStat->SetMarkerStyle(33);
  mTppMBStat->SetMarkerSize(1.8);
  mTppMBStat->Draw("pez same");

//  leg->AddEntry((TObject*)0, "        ", "");
//  leg->AddEntry((TObject*)0, " MB ", "");
//  leg->AddEntry((TObject*)0, " HM ", "");
  leg->AddEntry((TObject*)0, "p#minus#kern[0.4]{p}", "");
  leg->AddEntry(mTppMBStat, "    ", "pef");
  leg->AddEntry(mTppHMStat, "    ", "pef");
  leg->AddEntry((TObject*)0, "p#minus#kern[0.2]{#Lambda}", "");
  leg->AddEntry(mTpLMBStat, "    ", "pef");
  leg->AddEntry(mTpLHMStat, "    ", "pef");

  TLatex BeamText;
  BeamText.SetTextFont(43);
  BeamText.SetTextSize(40);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.55, 0.83,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.55, 0.77, "Gaussian Source");

  TLatex BeamText2;
  BeamText2.SetTextFont(43);
  BeamText2.SetTextSize(40);
  BeamText2.SetNDC(kTRUE);
  BeamText2.DrawLatex(0.63, 0.68, "MB");
  BeamText2.DrawLatex(0.725, 0.68, "HM");

  leg->Draw("same");
  c4->SaveAs(Form("%s/mTvsRad%s.pdf", gSystem->pwd(), sourceName));
  c4->Write();
  out->Write();
  out->Close();
  ppMBFile->Close();
  ppHMFile->Close();
  pLMBFile->Close();
  pLHMFile->Close();
}


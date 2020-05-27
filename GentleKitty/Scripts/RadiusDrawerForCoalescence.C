#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "DreamPlot.h"
#include "TFile.h"
#include "TDatabasePDG.h"

int main(int argc, char *argv[]) {
  if (!argv[1]) {
    std::cout << "pp RadFile from latest task missing\n";
    return -1;
  }
  if (!argv[2]) {
    std::cout << "pp  RadFile  from Bernie missing\n";
    return -1;
  }
  if (!argv[2]) {
    std::cout << "Source Name\n";
    return -1;
  }

  const char *ppFile = argv[1];
  const char *pLNLOFile = argv[2];
  const char *sourceName = argv[3];
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile *ppHMFile = TFile::Open(ppFile, "read");
  TGraphErrors *mTppHMSys = (TGraphErrors*) ppHMFile->Get("mTRadiusSyst");
  TGraphErrors *mTppHMStat = (TGraphErrors*) ppHMFile->Get("mTRadiusStat");

  TFile *pLNLOHMFile = TFile::Open(pLNLOFile, "read");
  TGraphErrors *mTpLNLOHMSys = (TGraphErrors*) pLNLOHMFile->Get("mTRadiusSyst");
  TGraphErrors *mTpLNLOHMStat = (TGraphErrors*) pLNLOHMFile->Get(
      "mTRadiusStat");

  double yMin = 1234567;
  double yMax = 0;
  double x, y;
  for (int iBin = 0; iBin < mTppHMSys->GetN(); iBin++) {
    mTppHMSys->GetPoint(iBin, x, y);
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
    mTpLNLOHMSys->GetPoint(iBin, x, y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    mTpLNLOHMSys->SetPointError(iBin, 0.4 * mTpLNLOHMSys->GetErrorX(iBin),
                                mTpLNLOHMSys->GetErrorY(iBin));
  }

  TFile *out = TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  TLegend *leg = new TLegend(0.18, 0.19, 0.53, 0.48);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
//  leg->SetNColumns(2);
  leg->SetTextSizePixels(40);

  mTppHMSys->SetTitle(
      "; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); #it{r}_{core} (fm)");
  mTppHMSys->GetXaxis()->SetTitleSize(40);
  mTppHMSys->GetYaxis()->SetTitleSize(40);
  mTppHMSys->GetXaxis()->SetTitleOffset(1.35);
  mTppHMSys->GetYaxis()->SetTitleOffset(1.4);
  mTppHMSys->GetXaxis()->SetLabelSize(40);
  mTppHMSys->GetYaxis()->SetLabelSize(40);
  mTppHMSys->GetXaxis()->SetLabelOffset(.02);
  mTppHMSys->GetYaxis()->SetLabelOffset(.02);
  mTppHMSys->SetFillColorAlpha(65 + 2, 0.7);
  mTppHMSys->SetFillStyle(3225);
  mTppHMSys->GetXaxis()->SetRangeUser(0.95, 1.8);
  mTppHMSys->GetXaxis()->SetLimits(0.95, 1.9);
  mTppHMSys->GetYaxis()->SetRangeUser(0.9 * yMin, 1.1 * yMax);
  mTppHMSys->Draw("2AP");

  mTppHMStat->SetMarkerColor(kRed + 1);
  mTppHMStat->SetLineColor(kRed + 1);
  mTppHMStat->SetMarkerColorAlpha(65 + 2, 1);
  mTppHMStat->SetLineColorAlpha(65 + 2, 1);
  mTppHMStat->SetFillColorAlpha(65 + 2, 0.7);
  mTppHMStat->SetLineWidth(1);
  mTppHMStat->SetFillStyle(3225);
  mTppHMStat->SetMarkerStyle(41);
  mTppHMStat->SetMarkerSize(2.0);
  mTppHMStat->Draw("pez same");

  mTpLNLOHMSys->SetFillColorAlpha(kBlue + 2, 0.7);
  mTpLNLOHMSys->SetFillStyle(3225);
  mTpLNLOHMSys->Draw("2PSAME");

  mTpLNLOHMStat->SetMarkerColor(kBlue + 2);
  mTpLNLOHMStat->SetLineColor(kBlue + 2);
  mTpLNLOHMStat->SetFillColorAlpha(kBlue + 2, 0.7);
  mTpLNLOHMStat->SetLineWidth(1);
  mTpLNLOHMStat->SetFillStyle(3225);
  mTpLNLOHMStat->SetMarkerStyle(34);
  mTpLNLOHMStat->SetMarkerSize(2.0);
  mTpLNLOHMStat->Draw("pez same");

  leg->AddEntry(mTppHMStat, "p#minus#kern[0.4]{p} bernie", "pef");
  leg->AddEntry(mTpLNLOHMStat, "p#minus#kern[0.4]{p} fixed bug", "pef");

  TLatex BeamText;
  BeamText.SetTextFont(43);
  BeamText.SetTextSize(40);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.48, 0.86,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp",(int) 13));
  BeamText.DrawLatex(0.48, 0.79, "Min-Bias");
  BeamText.DrawLatex(0.48, 0.72, TString::Format("%s", sourceName).Data());
// BeamText.DrawLatex(0.48, 0.65, "Avg Ref. Mult_{#eta < 0.8} ~23");

  leg->Draw("same");
  c4->SaveAs(Form("%s/ComparisionofRadiusWithBenie.pdf", gSystem->pwd()));
  c4->Write();
  out->Write();
  out->Close();
  ppHMFile->Close();
}

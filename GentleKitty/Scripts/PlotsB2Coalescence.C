#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "DreamPlot.h"
#include "TFile.h"
#include "TDatabasePDG.h"

int main(int argc, char* argv[]) {
  if (!argv[1]) {
    std::cout << "B2 data missing\n";
    return -1;
  }
  if (!argv[2]) {
    std::cout << "B2 theory missing\n";
    return -1;
  }
  if (!argv[3]) {
    std::cout << "Source Name\n";
    return -1;
  }
  const char* ppFile = argv[1];
  const char* pLFile = argv[2];
  const char* sourceName = argv[3];
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile* B2mTFile =
    TFile::Open(
      ppFile,
      "read");
  //TGraphAsymmErrors* B2mTSys = (TGraphAsymmErrors*) B2mTFile->Get("grB2stat_mt");
  //TGraphAsymmErrors* B2mTStat = (TGraphAsymmErrors*) B2mTFile->Get("grB2syst_mt");
  TGraphErrors* B2mTSys = (TGraphErrors*) B2mTFile->Get("grB2stat_mt");
  TGraphErrors* B2mTStat = (TGraphErrors*) B2mTFile->Get("grB2syst_mt");
  std::cout << "Here reached -1 \n";
  TFile* B2mTModelFile =
    TFile::Open(
      pLFile,
      "read");
  TGraphErrors* B2mTGauss_syst = (TGraphErrors*) B2mTModelFile->Get("B2eff_Gauss");
  TGraphErrors* B2mTTwoGauss_syst = (TGraphErrors*) B2mTModelFile->Get("B2eff_TwoGauss");
  TGraphErrors* B2mTHulthen_syst = (TGraphErrors*) B2mTModelFile->Get("B2eff_Hulthen");
  TGraphErrors* B2mTChiEFT_syst = (TGraphErrors*) B2mTModelFile->Get("B2eff_ChiEFT");
  /*double yMin = 1234567;
  //double yMax = 0;
  // double x,y;
  //std::cout << "Here reached -2 \n";
  //for (int iBin = 0; iBin < B2mTSys->GetN(); iBin++) {
    B2mTSys->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    B2mTSys->SetPointError(iBin, 0.4 * B2mTSys->GetErrorX(iBin),
                             B2mTSys->GetErrorY(iBin));
  }

  for (int iBin = 0; iBin < B2mTGauss_syst->GetN(); iBin++) {
    B2mTGauss_syst->GetPoint(iBin, x,y);
    if (y < yMin) {
      yMin = y;
    }
    if (yMax < y) {
      yMax = y;
    }
    B2mTGauss_syst->SetPointError(iBin, 0.4 * B2mTGauss_syst->GetErrorX(iBin),
                             B2mTGauss_syst->GetErrorY(iBin));
  }*/
  TFile* out = TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();

  TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.87);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetNColumns(1);
  leg->SetTextSizePixels(40);

  //B2mTSys->SetPoint(B2mTSys->GetN(), 0.95, 0.5);
  //B2mTSys->SetPointError(B2mTSys->GetN(), 0., 0.);

// B2mTSys->SetPoint(B2mTSys->GetN(), 2.7, 1.3);
// B2mTSys->SetPointError(B2mTSys->GetN(), 0., 0.);
  B2mTStat->SetTitle(
    "; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); #it{B}_{2} (GeV^{2}/c^{3})");
  B2mTStat->GetXaxis()->SetTitleSize(40);
  B2mTStat->GetYaxis()->SetTitleSize(40);
  B2mTStat->GetXaxis()->SetTitleOffset(1.35);
  B2mTStat->GetYaxis()->SetTitleOffset(1.4);
  B2mTStat->GetXaxis()->SetLabelSize(40);
  B2mTStat->GetYaxis()->SetLabelSize(40);
  B2mTStat->GetXaxis()->SetLabelOffset(.02);
  B2mTStat->GetYaxis()->SetLabelOffset(.02);
  B2mTStat->GetXaxis()->SetLimits(0.9, 1.9);
  B2mTStat->GetYaxis()->SetRangeUser(0.001, 0.03);


  B2mTStat->SetMarkerColor(kRed + 1);
  B2mTStat->SetLineColor(kRed + 1);
  B2mTStat->SetFillColorAlpha(kRed - 7, 0.8);
  B2mTStat->SetLineWidth(1);
  B2mTStat->SetFillStyle(3001);
  B2mTStat->SetMarkerStyle(41);
  B2mTStat->SetMarkerSize(2.0);
  B2mTStat->Draw("2AP");


  B2mTSys->SetMarkerColor(kRed + 1);
  B2mTSys->SetLineColor(kRed + 1);
  B2mTSys->SetFillColorAlpha(kRed - 7, 0.8);
  B2mTSys->SetFillStyle(3001);

  B2mTSys->Draw("pez same");

  B2mTGauss_syst->SetFillColorAlpha(kYellow + 2, 0.7);
  B2mTGauss_syst->SetFillStyle(3225);
  B2mTGauss_syst->SetMarkerColor(kYellow + 2);
  B2mTGauss_syst->SetLineColor(kYellow + 2);
  B2mTGauss_syst->Draw("2PSAME");

  B2mTTwoGauss_syst->SetFillColorAlpha(kGreen - 2, 0.7);
  B2mTTwoGauss_syst->SetFillStyle(3225);
  B2mTTwoGauss_syst->SetMarkerColor(kGreen + 2);
  B2mTTwoGauss_syst->SetLineColor(kGreen + 2);
  B2mTTwoGauss_syst->Draw("2PSAME");

  B2mTHulthen_syst->SetFillColorAlpha(65 + 1, 0.7);
  B2mTHulthen_syst->SetFillStyle(3225);
  B2mTHulthen_syst->SetMarkerColor(65 + 2);
  B2mTHulthen_syst->SetLineColor(65 + 2);
  B2mTHulthen_syst->Draw("2PSAME");

  B2mTChiEFT_syst->SetFillColorAlpha(kBlue + 2, 0.7);
  B2mTChiEFT_syst->SetFillStyle(3225);
  B2mTChiEFT_syst->SetMarkerColor(kBlue + 2);
  B2mTChiEFT_syst->SetLineColor(kBlue + 2);
  B2mTChiEFT_syst->Draw("2PSAME");
//  leg->AddEntry((TObject*)0, "        ", "");
//  leg->AddEntry((TObject*)0, " MB ", "");
//  leg->AddEntry((TObject*)0, " HM ", "");
//  leg->AddEntry((TObject*)0, "p#minus#kern[0.4]{p}", "");
//  leg->AddEntry(mTppMBStat, "    ", "pef");
//  leg->AddEntry(B2mTStat, "    ", "pef");
//  leg->AddEntry((TObject*)0, "p#minus#kern[0.2]{#Lambda}", "");
//  leg->AddEntry(mTpLMBStat, "    ", "pef");
//  leg->AddEntry(mTpLHMStat, "    ", "pef");
  leg->AddEntry(B2mTStat, "#it{B}_{2} ALICE", "pef");
  leg->AddEntry(B2mTChiEFT_syst, "Chiral EFT", "pef");
  leg->AddEntry(B2mTHulthen_syst, "Hulthen", "pef");
  leg->AddEntry(B2mTTwoGauss_syst, "Two gaussians", "pef");
  leg->AddEntry(B2mTGauss_syst, "Gaussian", "pef");

  TLatex BeamText;
  BeamText.SetTextFont(43);
  BeamText.SetTextSize(40);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.54, 0.83,
                     Form("ALICE %s #sqrt{#it{s}} = %0.2f TeV", "p-Pb", (float) 5.02));
  //BeamText.DrawLatex(0.55, 0.76, "INEL");
  //BeamText.DrawLatex(
  //    0.55,
  ////   0.69,
  //  "(0#kern[-0.95]{ }#minus#kern[-0.05]{ }0.17#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
//  BeamText.DrawLatex(0.55, 0.62, "Gaussian Source");
  //BeamText.DrawLatex(0.17, 0.2, TString::Format("%s",sourceName).Data());

  leg->Draw("same");
  c4->SaveAs(Form("%s/B2plotpPb.pdf", gSystem->pwd()));
  c4->Write();
  out->Write();
  out->Close();
  B2mTFile->Close();
  B2mTModelFile->Close();
}


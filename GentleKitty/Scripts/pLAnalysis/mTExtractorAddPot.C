#include "VariationmTAnalysis.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

void DrawRadiusVsMT(const char* filepp, const char* filepL,
                    const char* sourceName);

int main(int argc, char* argv[]) {
//  gROOT->ProcessLine("gErrorIgnoreLevel = 2001");
  const char* DataDir = argv[1];
  const char* FitDirResonance = argv[2];
  const char* FitDirGauss = argv[3];
  const char* ppFileResonance = (argv[4]) ? argv[4] : "";
  const char* ppFileGauss = (argv[5]) ? argv[5] : "";
  const int nMTBins = 6;
  TString avgmTFile = TString::Format("%s/AveragemT.root", DataDir);
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  if (!mTFile) {
    std::cout << "No mT File \n";
    return -1;
  }
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT_pLVar0");
  if (!avgmT) {
    std::cout << "no average mT graph " << std::endl;
    mTFile->ls();
    return -1;
  }

  VariationmTAnalysis* analyser = new VariationmTAnalysis(2, 42, 162);
  analyser->SetHistName("hCk_RebinnedpLVar");
  analyser->SetFileName("OutFileVarpL.root");
  analyser->SetmTAverage(avgmT);
  analyser->SetLegData("p-#Lambda #oplus #bar{p}-#bar{#Lambda}", "fpe");
  analyser->SetLegModel("Fit Gaussian + Resonance Source", "f", 1, 3252);
  analyser->SetLegModel("Fit Gaussian Source", "f", 1, 3225);
  analyser->SetTextXMin(0.35);
  analyser->SetPlottingRange(0, 230);
  //this implies a folder structure following mTBin_[1,2,3,....]
  for (int imt = 1; imt <= nMTBins; ++imt) {
    std::cout << "========================== \n";
    std::cout << "mTBin: " << imt << std::endl;
    std::cout << "========================== \n";
    TString mTDataDir = Form("%s/mTBin_%u", DataDir, imt);
    analyser->SetSystematic(mTDataDir.Data());
    TString mTFitDirResonance = Form("%s/mTBin_%u/", FitDirResonance, imt);
    analyser->SetVariation(mTFitDirResonance.Data(), 0);
    TString mTFitDirGauss = Form("%s/mTBin_%u/", FitDirGauss, imt);
    analyser->SetVariation(mTFitDirGauss.Data(), 1);
  }
  std::vector<float> mTBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
  analyser->SetmTBins(mTBins);
  analyser->MakeCFPlotsSingleBand();
  analyser->StoreRadvsmT("RadpLvsmTGaussReso.root", 0);
  analyser->StoreRadvsmT("RadpLvsmTGauss.root", 1);

  if (ppFileResonance != "" && ppFileGauss != "") {
    std::cout << "pp file Name: " << ppFileGauss << " and " << ppFileResonance
              << " passed\n";
    DrawRadiusVsMT(ppFileResonance, "RadpLvsmTGaussReso.root",
                   "Gauss+Resonances");
    DrawRadiusVsMT(ppFileGauss, "RadpLvsmTGauss.root", "Gauss");
  }

  return 1;
}

void DrawRadiusVsMT(const char* filepp, const char* filepL,
                    const char* sourceName) {
  DreamPlot::SetStyle();
  TFile* ppFile = TFile::Open(filepp, "read");
  TFile* pLFile = TFile::Open(filepL, "read");
  TGraphErrors* mTppSys = (TGraphErrors*) ppFile->Get("mTRadiusSyst");
  TGraphErrors* mTppStat = (TGraphErrors*) ppFile->Get("mTRadiusStat");
  TGraphErrors* mTpLSys = (TGraphErrors*) pLFile->Get("mTRadiusSyst");
  TGraphErrors* mTpLStat = (TGraphErrors*) pLFile->Get("mTRadiusStat");
  for (int iBin = 0; iBin < mTppSys->GetN(); iBin++) {
    mTppSys->SetPointError(iBin, 0.4 * mTppSys->GetErrorX(iBin),
                           mTppSys->GetErrorY(iBin));
  }
  for (int iBin = 0; iBin < mTpLSys->GetN(); iBin++) {
    mTpLSys->SetPointError(iBin, 0.4 * mTpLSys->GetErrorX(iBin),
                           mTpLSys->GetErrorY(iBin));
  }
  TFile* out = TFile::Open(Form("%s.root",sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  TLegend* leg = new TLegend(0.55, 0.5, 0.85, 0.7);
  leg->SetFillStyle(4000);
  leg->SetTextFont(43);
  leg->SetTextSizePixels(32);
  mTpLSys->SetPoint(mTpLSys->GetN(), 0.95, 0.5);
  mTpLSys->SetPointError(mTpLSys->GetN(), 0., 0.);

  mTpLSys->SetPoint(mTpLSys->GetN(), 2.7, 1.3);
  mTpLSys->SetPointError(mTpLSys->GetN(), 0., 0.);

  mTpLSys->SetTitle("; < m_{T} >  (MeV/#it{c}^{2}); r_{Core} (fm)");
  mTpLSys->GetXaxis()->SetTitleSize(40);
  mTpLSys->GetYaxis()->SetTitleSize(40);
  mTpLSys->GetXaxis()->SetTitleOffset(1.4);
  mTpLSys->GetYaxis()->SetTitleOffset(1.4);
  mTpLSys->GetXaxis()->SetLabelSize(40);
  mTpLSys->GetYaxis()->SetLabelSize(40);
  mTpLSys->GetXaxis()->SetLabelOffset(.02);
  mTpLSys->GetYaxis()->SetLabelOffset(.02);
  mTpLSys->GetXaxis()->SetRangeUser(0.95, 2.7);

  if (strcmp(sourceName, "Gauss") == 0) {
    mTpLSys->GetYaxis()->SetRangeUser(0.8, 1.9);
  } else {
    mTpLSys->GetYaxis()->SetRangeUser(0.45, 1.4);
  }

  mTpLSys->SetFillColorAlpha(kRed - 7, 0.7);
  mTpLSys->Draw("2AP");

  mTpLStat->SetMarkerColor(kRed + 1);
  mTpLStat->SetLineColor(kRed + 1);
  mTpLStat->SetFillColor(kRed + 1);
  mTpLStat->SetLineWidth(1);
  mTpLStat->SetMarkerStyle(23);
  mTpLStat->SetMarkerSize(1.2);
  mTpLStat->Draw("pez same");

  mTppSys->SetFillColorAlpha(kBlue - 7, 0.7);
  mTppSys->Draw("2PSAME");

  mTppStat->SetMarkerColor(kBlue + 2);
  mTppStat->SetLineColor(kBlue + 2);
  mTppStat->SetFillColor(kBlue + 2);
  mTppStat->SetLineWidth(1);
  mTppStat->SetMarkerStyle(20);
  mTppStat->SetMarkerSize(1.2);
  mTppStat->Draw("pez same");

  leg->AddEntry(mTppStat, "p#minus p", "pl");
  leg->AddEntry(mTpLStat, "p#minus #Lambda", "pl");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize() * .85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.55, 0.83,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.55, 0.77, "High Mult. (0-0.072% INEL)");
  BeamText.DrawLatex(0.55, 0.72, Form("%s", sourceName));

  leg->Draw("same");
  c4->SaveAs(Form("%s/mTvsRad%s.pdf", gSystem->pwd(), sourceName));
  c4->Write();
  out->Write();
  out->Close();
  ppFile->Close();
  pLFile->Close();
}

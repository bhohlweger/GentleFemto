#include "TString.h"
#include "TFile.h"
#include "TColor.h"
#include "TNtuple.h"
#include "DreamPlot.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

void PlotSigmaExclusion(TString inputDir) {
  DreamPlot::SetStyle();

  const int NRGBs = 5;
  Double_t stops[NRGBs];
  for (int i = 0; i < NRGBs; ++i) {
    stops[i] = float(i) / (NRGBs - 1);
  }
  Double_t red[NRGBs] = { 1., 29. / 255., 25. / 255., 27. / 255., 32. / 255. };
  Double_t green[NRGBs] = { 1., 221. / 255., 160. / 255., 113. / 255., 74.
      / 255. };
  Double_t blue[NRGBs] = { 1., 221. / 255., 184. / 255., 154. / 255., 129.
      / 255. };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);

  double contours[3];
  contours[0] = 2.3;
  contours[1] = 6.18;
  contours[2] = 11.83;

  auto file = TFile::Open(TString::Format("%s/Result.root", inputDir.Data()));
  auto tuple = (TNtuple*) file->Get("exclusion");
  Float_t bestChi2, defaultChi2, worstChi2, d0, f0inv, CFneg;
  tuple->SetBranchAddress("CFneg", &CFneg);
  tuple->SetBranchAddress("d0", &d0);
  tuple->SetBranchAddress("f0inv", &f0inv);
  tuple->SetBranchAddress("bestChi2", &bestChi2);
  tuple->SetBranchAddress("defChi2", &defaultChi2);
  tuple->SetBranchAddress("worstChi2", &worstChi2);

  // Get the best value of chi2
  auto histBestChi2 = new TH1F("histBestChi2", "; #chi^{2}_{best}; Entries",
                               10000, 0, 100);
  auto histDefChi2 = new TH1F("histDefChi2", "; #chi^{2}_{default}; Entries",
                              10000, 0, 100);
  auto histWorstChi2 = new TH1F("histWorstChi2", "; #chi^{2}_{worst}; Entries",
                                10000, 0, 100);
  DreamPlot::SetStyleHisto(histBestChi2);
  DreamPlot::SetStyleHisto(histDefChi2);
  DreamPlot::SetStyleHisto(histWorstChi2);

  tuple->Draw("bestChi2>>histBestChi2");
  tuple->Draw("defChi2>>histDefChi2");
  tuple->Draw("worstChi2>>histWorstChi2");

  float bestbestChi2 = histBestChi2->GetXaxis()->GetBinLowEdge(
      histBestChi2->FindFirstBinAbove(0.1, 1));
  float bestDefChi2 = histDefChi2->GetXaxis()->GetBinLowEdge(
      histDefChi2->FindFirstBinAbove(0.1, 1));
  float bestworstChi2 = histWorstChi2->GetXaxis()->GetBinLowEdge(
      histWorstChi2->FindFirstBinAbove(0.1, 1));

  std::cout << "Best chi2  : " << bestbestChi2 << "\n";
  std::cout << "Def chi2   : " << bestDefChi2 << "\n";
  std::cout << "Worst chi2 : " << bestworstChi2 << "\n";

  auto mapRelWorst = new TH2F(
      "mapRelWorst", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
      200, -5, 5, 200, 0, 10);
  auto mapWorstChi2 = new TH2F(
      "mapWorstChi2", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #chi^{2}", 200,
      -5, 5, 200, 0, 10);
  auto mapRelDef = new TH2F(
      "mapRelDef", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
      200, -5, 5, 200, 0, 10);
  auto mapDefChi2 = new TH2F(
      "mapDefChi2", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #chi^{2}", 200,
      -5, 5, 200, 0, 10);
  auto mapRelBest = new TH2F(
      "mapRelBest", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
      200, -5, 5, 200, 0, 10);
  auto mapBestChi2 = new TH2F(
      "mapBestChi2", ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #chi^{2}", 200,
      -5, 5, 200, 0, 10);
  auto mapCFneg = new TH2F("mapCFneg",
                           ";1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); CFneg",
                           200, -5, 5, 200, 0, 10);
  const float delta = 1E-6;

  // Iterate over the data to obtain the results
  for (int i = 0; i < tuple->GetEntriesFast(); ++i) {
    tuple->GetEntry(i);
    if(CFneg > 0) mapCFneg->Fill(f0inv + delta, d0 + delta, 1);
    mapWorstChi2->Fill(f0inv + delta, d0 + delta, worstChi2);
    mapRelWorst->Fill(f0inv + delta, d0 + delta, worstChi2 - bestworstChi2);
    mapDefChi2->Fill(f0inv + delta, d0 + delta, defaultChi2);
    mapRelDef->Fill(f0inv + delta, d0 + delta, defaultChi2 - bestDefChi2);
    mapBestChi2->Fill(f0inv + delta, d0 + delta, bestChi2);
    mapRelBest->Fill(f0inv + delta, d0 + delta, bestChi2 - bestbestChi2);
  }

  // Get the lednicky sucks contour
  int nBinsX = mapCFneg->GetXaxis()->GetNbins();
  float xMin = mapCFneg->GetXaxis()->GetBinLowEdge(1);
  float xMax = mapCFneg->GetXaxis()->GetBinUpEdge(nBinsX);
  int nBinsY = mapCFneg->GetYaxis()->GetNbins();
  float yMin = mapCFneg->GetYaxis()->GetBinLowEdge(1);
  float yMax = mapCFneg->GetYaxis()->GetBinUpEdge(nBinsY);
  auto lednickyReallySucks = new TGraphAsymmErrors();
  int iPoint = 0;
  for (int iXbin = 1; iXbin <= nBinsX; ++iXbin) {
    float yBinMin = 0;
    float yBinMax = yMax;
    float binWidth = 0;
    for (int iYbin = 1; iYbin <= nBinsY; ++iYbin) {
      binWidth = mapCFneg->GetYaxis()->GetBinWidth(iYbin);
      if (mapCFneg->GetBinContent(iXbin, iYbin) != 0) {
        yBinMin = mapCFneg->GetYaxis()->GetBinLowEdge(iYbin);
        break;
      }
    }
    if (yBinMin == 0) {
      continue;
    }
    float err = yBinMax - yBinMin;
    float xValue = mapCFneg->GetXaxis()->GetBinLowEdge(iXbin);
    lednickyReallySucks->SetPoint(iPoint, xValue, yBinMin);
    lednickyReallySucks->SetPointError(iPoint++, 0, 0, 0, 20);
  }
  lednickyReallySucks->SetFillStyle(3004);
  lednickyReallySucks->SetLineWidth(0);
  lednickyReallySucks->SetFillColor(kBlack);
  lednickyReallySucks->SetLineColorAlpha(kBlack, 0.7);

  auto fileOut = new TFile(Form("%s/exclusionOutput.root", inputDir.Data()), "RECREATE");
  fileOut->cd();

  auto cBestChi2 = new TCanvas();
  histBestChi2->Rebin(10);
  histBestChi2->Draw();
  cBestChi2->Print(Form("%s/Chi2Best.pdf", inputDir.Data()));
  cBestChi2->Print(Form("%s/Chi2Best.png", inputDir.Data()));
  histBestChi2->Write("Chi2Best");

  auto cDefChi2 = new TCanvas();
  histDefChi2->Rebin(10);
  histDefChi2->Draw();
  cDefChi2->Print(Form("%s/Chi2Def.pdf", inputDir.Data()));
  cDefChi2->Print(Form("%s/Chi2Def.png", inputDir.Data()));
  histDefChi2->Write("Chi2Def");

  auto cWorstChi2 = new TCanvas();
  histWorstChi2->Rebin(10);
  histWorstChi2->Draw();
  cWorstChi2->Print(Form("%s/Chi2Worst.pdf", inputDir.Data()));
  cWorstChi2->Print(Form("%s/Chi2Worst.png", inputDir.Data()));
  histWorstChi2->Write("Chi2Worst");

  auto cWorst = new TCanvas();
  cWorst->SetRightMargin(0.16);
  mapWorstChi2->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  cWorst->Print(Form("%s/MapChi2Worst.pdf", inputDir.Data()));
  cWorst->Print(Form("%s/MapChi2Worst.png", inputDir.Data()));
  mapWorstChi2->Write("MapChi2Worst");

  auto dWorst = new TCanvas();
  dWorst->SetRightMargin(0.16);
  mapRelWorst->SetMaximum(25);
  mapRelWorst->SetContour(3, contours);
  mapRelWorst->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  dWorst->Print(Form("%s/ExclusionWorst.pdf", inputDir.Data()));
  dWorst->Print(Form("%s/ExclusionWorst.png", inputDir.Data()));
  mapRelWorst->Write("ExclusionWorst");
  dWorst->Write("ExclusionWorstPlot");

  auto cDef = new TCanvas();
  cDef->SetRightMargin(0.16);
  mapDefChi2->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  cDef->Print(Form("%s/MapChi2Def.pdf", inputDir.Data()));
  cDef->Print(Form("%s/MapChi2Def.png", inputDir.Data()));
  mapDefChi2->Write("MapChi2Def");

  auto dDef = new TCanvas();
  dDef->SetRightMargin(0.16);
  mapRelDef->SetMaximum(25);
  mapRelDef->SetContour(3, contours);
  mapRelDef->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  dDef->Print(Form("%s/ExclusionDef.pdf", inputDir.Data()));
  dDef->Print(Form("%s/ExclusionDef.png", inputDir.Data()));
  mapRelDef->Write("ExclusionDef");
  dDef->Write("ExclusionDefPlot");

  auto cBest = new TCanvas();
  cBest->SetRightMargin(0.16);
  mapBestChi2->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  cBest->Print(Form("%s/MapChi2Best.pdf", inputDir.Data()));
  cBest->Print(Form("%s/MapChi2Best.png", inputDir.Data()));
  mapBestChi2->Write("MapChi2Best");

  auto dBest = new TCanvas();
  dBest->SetRightMargin(0.16);
  mapRelBest->SetMaximum(25);
  mapRelBest->SetContour(3, contours);
  mapRelBest->Draw("colz");
  lednickyReallySucks->Draw("3 SAME");
  dBest->Print(Form("%s/ExclusionBest.pdf", inputDir.Data()));
  dBest->Print(Form("%s/ExclusionBest.png", inputDir.Data()));
  mapRelBest->Write("ExclusionBest");
  dBest->Write("ExclusionBestPlot");


  auto lednicky = new TCanvas();
  dBest->SetRightMargin(0.16);
  mapCFneg->Draw("colz");
  lednicky->Print(Form("%s/LednickyFails.pdf", inputDir.Data()));
  lednicky->Print(Form("%s/LednickyFails.png", inputDir.Data()));
  mapCFneg->Write("LednickyFails");

  auto lednickyCont = new TCanvas();
  dBest->SetRightMargin(0.16);
  lednickyReallySucks->Draw("A3");
  lednickyCont->Print(Form("%s/LednickyFailsCountour.pdf", inputDir.Data()));
  lednickyCont->Print(Form("%s/LednickyFailsCountour.png", inputDir.Data()));
  lednickyReallySucks->Write("LednickyFailsCountour");

  fileOut->Close();
  file->Close();
}

int main(int argc, char *argv[]) {
  PlotSigmaExclusion(argv[1]);
  return 0;
}

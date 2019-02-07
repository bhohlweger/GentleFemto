#include "TString.h"
#include "TFile.h"
#include "TColor.h"
#include "TNtuple.h"
#include "DreamPlot.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TDatabasePDG.h"
#include "TLatex.h"
#include "TStyle.h"

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

  const float massSigma = TDatabasePDG::Instance()->GetParticle(3212)->Mass();
  const float massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const float effMass = 0.5f * 1000.f * (massSigma + massProton);  // in MeV

  auto file = TFile::Open(TString::Format("%s/Result.root", inputDir.Data()));
  auto tuple = (TNtuple*) file->Get("exclusion");
  Float_t bestChi2, defaultChi2, worstChi2, d0, REf0inv, IMf0inv, CFneg;
  tuple->SetBranchAddress("CFneg", &CFneg);
  tuple->SetBranchAddress("d0", &d0);
  tuple->SetBranchAddress("REf0inv", &REf0inv);
  tuple->SetBranchAddress("IMf0inv", &IMf0inv);
  tuple->SetBranchAddress("bestChi2", &bestChi2);
  tuple->SetBranchAddress("defChi2", &defaultChi2);
  tuple->SetBranchAddress("worstChi2", &worstChi2);

  // Get the best value of chi2
  auto histBestChi2 = new TH1F("histBestChi2", "; #chi^{2}_{best}; Entries",
                               100000, 0, 100);
  auto histDefChi2 = new TH1F("histDefChi2", "; #chi^{2}_{default}; Entries",
                              100000, 0, 100);
  auto histWorstChi2 = new TH1F("histWorstChi2", "; #chi^{2}_{worst}; Entries",
                                100000, 0, 100);
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

  std::vector<float> IMf0invCont = { { -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5,
      -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, } };
  const int nHists = IMf0invCont.size();
  TH2F* mapRelWorst[nHists];
  TH2F* mapWorstChi2[nHists];
  TH2F* mapRelDef[nHists];
  TH2F* mapDefChi2[nHists];
  TH2F* mapRelBest[nHists];
  TH2F* mapBestChi2[nHists];
  TH2F* mapCFneg[nHists];
  TH2F* mapBinding[nHists];

  const int nBinsX = 100;
  const int nBinsY = 200;
  for (size_t i = 0; i < IMf0invCont.size(); ++i) {

    mapRelWorst[i] =
        new TH2F(
            Form("mapRelWorst_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapWorstChi2[i] =
        new TH2F(
            Form("mapWorstChi2_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapRelDef[i] =
        new TH2F(
            Form("mapRelDef_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapDefChi2[i] =
        new TH2F(
            Form("mapDefChi2_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapRelBest[i] =
        new TH2F(
            Form("mapRelBest_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapBestChi2[i] =
        new TH2F(
            Form("mapBestChi2_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); #chi^{2}",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapCFneg[i] =
        new TH2F(
            Form("mapCFneg_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;1/#Rgothic(#it{f}_{0}) (fm^{-1}); #it{d}_{0} (fm); CFneg",
                IMf0invCont[i]),
            nBinsX, -5, 5, nBinsY, -10, 10);
    mapBinding[i] =
        new TH2F(
            Form("mapBinding_%.1f", IMf0invCont[i]),
            Form(
                "#Jgothic(f_{0}^{-1}) = %.1f fm^{-1} ;B_{p#Sigma^{0}} (MeV); #it{d}_{0} (fm); #Delta#chi^{2}",
                IMf0invCont[i]),
            50, 0, 15, 50, 0, 5);
  }
  const float delta = 1E-6;

// Iterate over the data to obtain the results
  const float IMf00 = IMf0invCont[0];
  int imf0 = 0;
  float kappa = 0.f;
  float binding = 0.f;
  const float hbarc = 197.327;
  for (int i = 0; i < tuple->GetEntriesFast(); ++i) {
    tuple->GetEntry(i);
    imf0 = (IMf0inv - IMf00) * 2.f;

    if (CFneg > 0) {
      mapCFneg[imf0]->Fill(REf0inv + delta, d0 + delta, 1);
    }
    mapWorstChi2[imf0]->Fill(REf0inv + delta, d0 + delta, worstChi2);
    mapRelWorst[imf0]->Fill(REf0inv + delta, d0 + delta,
                            worstChi2 - bestworstChi2);
    mapDefChi2[imf0]->Fill(REf0inv + delta, d0 + delta, defaultChi2);
    mapRelDef[imf0]->Fill(REf0inv + delta, d0 + delta,
                          defaultChi2 - bestDefChi2);
    mapBestChi2[imf0]->Fill(REf0inv + delta, d0 + delta, bestChi2);
    mapRelBest[imf0]->Fill(REf0inv + delta, d0 + delta,
                           bestChi2 - bestbestChi2);

    // compute binding energy
    if (REf0inv > -1.5 && REf0inv < 0 && d0 < 5 && d0 > 0
        && 2. * d0 * REf0inv < 1) {
      kappa = (1.f - std::sqrt(1.f - 2.f * d0 * REf0inv)) / d0 * hbarc;
      binding = kappa * kappa / effMass;

      if (defaultChi2 - bestDefChi2 < contours[0]) {
        mapBinding[imf0]->Fill(binding, d0 + delta, 1);
      }
    }
  }

  /// OUTPUT AND PLOTTING
  auto fileOut = new TFile(Form("%s/exclusionOutput.root", inputDir.Data()),
                           "RECREATE");
  TGraphAsymmErrors *lednickyReallySucks[nHists];

  auto cBestChi2 = new TCanvas();
  cBestChi2->SetLogy();
  histBestChi2->Rebin(100);
  histBestChi2->Draw();
  cBestChi2->Print(Form("%s/Chi2Best.pdf", inputDir.Data()));
  cBestChi2->Print(Form("%s/Chi2Best.png", inputDir.Data()));
  histBestChi2->Write("Chi2Best");

  auto cDefChi2 = new TCanvas();
  cDefChi2->SetLogy();
  histDefChi2->Rebin(100);
  histDefChi2->Draw();
  cDefChi2->Print(Form("%s/Chi2Def.pdf", inputDir.Data()));
  cDefChi2->Print(Form("%s/Chi2Def.png", inputDir.Data()));
  histDefChi2->Write("Chi2Def");

  auto cWorstChi2 = new TCanvas();
  cWorstChi2->SetLogy();
  histWorstChi2->Rebin(100);
  histWorstChi2->Draw();
  cWorstChi2->Print(Form("%s/Chi2Worst.pdf", inputDir.Data()));
  cWorstChi2->Print(Form("%s/Chi2Worst.png", inputDir.Data()));
  histWorstChi2->Write("Chi2Worst");
  auto cWorst = new TCanvas();
  auto dWorst = new TCanvas();
  auto cDef = new TCanvas();
  auto dDef = new TCanvas();
  auto cBest = new TCanvas();
  auto dBest = new TCanvas();
  auto lednicky = new TCanvas();
  auto dBinding = new TCanvas();

  TLatex BeamText;
  BeamText.SetNDC(kTRUE);
  BeamText.SetTextSize(0.85 * gStyle->GetTextSize());

  for (size_t imf0 = 0; imf0 < IMf0invCont.size(); ++imf0) {

    // Get the lednicky sucks contour
    int nBinsX = mapCFneg[imf0]->GetXaxis()->GetNbins();
    float xMin = mapCFneg[imf0]->GetXaxis()->GetBinLowEdge(1);
    float xMax = mapCFneg[imf0]->GetXaxis()->GetBinUpEdge(nBinsX);
    int nBinsY = mapCFneg[imf0]->GetYaxis()->GetNbins();
    float yMin = mapCFneg[imf0]->GetYaxis()->GetBinLowEdge(1);
    float yMax = mapCFneg[imf0]->GetYaxis()->GetBinUpEdge(nBinsY);
    lednickyReallySucks[imf0] = new TGraphAsymmErrors();
    int iPoint = 0;
    for (int iXbin = 1; iXbin <= nBinsX; ++iXbin) {
      float yBinMin = 0;
      float yBinMax = yMax;
      float binWidth = 0;
      for (int iYbin = 1; iYbin <= nBinsY; ++iYbin) {
        binWidth = mapCFneg[imf0]->GetYaxis()->GetBinWidth(iYbin);
        if (mapCFneg[imf0]->GetBinContent(iXbin, iYbin) != 0) {
          yBinMin = mapCFneg[imf0]->GetYaxis()->GetBinLowEdge(iYbin);
          break;
        }
      }
      if (yBinMin == 0) {
        continue;
      }
      float err = yBinMax - yBinMin;
      float xValue = mapCFneg[imf0]->GetXaxis()->GetBinLowEdge(iXbin);
      lednickyReallySucks[imf0]->SetPoint(iPoint, xValue, yBinMin);
      lednickyReallySucks[imf0]->SetPointError(iPoint++, 0, 0, 0, 20);
    }
    lednickyReallySucks[imf0]->SetFillStyle(3004);
    lednickyReallySucks[imf0]->SetLineWidth(0);
    lednickyReallySucks[imf0]->SetFillColor(kBlack);
    lednickyReallySucks[imf0]->SetLineColorAlpha(kBlack, 0.7);

    fileOut->cd();

    cWorst->Clear();
    cWorst->SetRightMargin(0.16);
    mapWorstChi2[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    cWorst->Print(
        Form("%s/MapChi2Worst_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    cWorst->Print(
        Form("%s/MapChi2Worst_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapWorstChi2[imf0]->Write(Form("MapChi2Worst_%.1f", IMf0invCont[imf0]));

    dWorst->Clear();
    dWorst->SetRightMargin(0.16);
    mapRelWorst[imf0]->SetMinimum(0);
    mapRelWorst[imf0]->SetMaximum(25);
    mapRelWorst[imf0]->SetContour(3, contours);
    mapRelWorst[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    dWorst->Print(
        Form("%s/ExclusionWorst_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    dWorst->Print(
        Form("%s/ExclusionWorst_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapRelWorst[imf0]->Write(Form("ExclusionWorst_%.1f", IMf0invCont[imf0]));
    dWorst->Write(Form("ExclusionWorstPlot_%.1f", IMf0invCont[imf0]));

    BeamText.DrawLatex(
        0.35,
        0.925,
        TString::Format("#Jgothic(f_{0}^{-1}) = %.3f fm^{-1}",
                        IMf0invCont[imf0]));
    dWorst->Print(Form("%s/ExclusionWorst.gif+", inputDir.Data()));
    dWorst->Print(Form("%s/ExclusionWorst.gif+", inputDir.Data()));
    dWorst->Print(Form("%s/ExclusionWorst.gif+", inputDir.Data()));

    cDef->Clear();
    cDef->SetRightMargin(0.16);
    mapDefChi2[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    cDef->Print(
        Form("%s/MapChi2Def_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    cDef->Print(
        Form("%s/MapChi2Def_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapDefChi2[imf0]->Write(Form("MapChi2Def_%.1f", IMf0invCont[imf0]));

    dDef->Clear();
    dDef->SetRightMargin(0.16);
    mapRelDef[imf0]->SetMinimum(0);
    mapRelDef[imf0]->SetMaximum(25);
    mapRelDef[imf0]->SetContour(3, contours);
    mapRelDef[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    dDef->Print(Form("%s/ExclusionDef_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    dDef->Print(Form("%s/ExclusionDef_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapRelDef[imf0]->Write(Form("ExclusionDef_%.1f", IMf0invCont[imf0]));
    dDef->Write(Form("ExclusionDefPlot_%.1f", IMf0invCont[imf0]));

    BeamText.DrawLatex(
        0.35,
        0.925,
        TString::Format("#Jgothic(f_{0}^{-1}) = %.3f fm^{-1}",
                        IMf0invCont[imf0]));
    dDef->Print(Form("%s/ExclusionDefault.gif+", inputDir.Data()));
    dDef->Print(Form("%s/ExclusionDefault.gif+", inputDir.Data()));
    dDef->Print(Form("%s/ExclusionDefault.gif+", inputDir.Data()));

    cBest->Clear();
    cBest->SetRightMargin(0.16);
    mapBestChi2[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    cBest->Print(
        Form("%s/MapChi2Best_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    cBest->Print(
        Form("%s/MapChi2Best_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapBestChi2[imf0]->Write(Form("MapChi2Best_%.1f", IMf0invCont[imf0]));

    dBest->Clear();
    dBest->SetRightMargin(0.16);
    mapRelBest[imf0]->SetMinimum(0);
    mapRelBest[imf0]->SetMaximum(25);
    mapRelBest[imf0]->SetContour(3, contours);
    mapRelBest[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 SAME");
    dBest->Print(
        Form("%s/ExclusionBest_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    dBest->Print(
        Form("%s/ExclusionBest_%.1f.png", inputDir.Data(), IMf0invCont[imf0]));
    mapRelBest[imf0]->Write(Form("ExclusionBest_%.1f", IMf0invCont[imf0]));
    dBest->Write(Form("ExclusionBestPlot_%.1f", IMf0invCont[imf0]));

    BeamText.DrawLatex(
        0.35,
        0.925,
        TString::Format("#Jgothic(f_{0}^{-1}) = %.3f fm^{-1}",
                        IMf0invCont[imf0]));
    dDef->Print(Form("%s/ExclusionBest.gif+", inputDir.Data()));
    dDef->Print(Form("%s/ExclusionBest.gif+", inputDir.Data()));
    dDef->Print(Form("%s/ExclusionBest.gif+", inputDir.Data()));

    lednicky->Clear();
    dBest->SetRightMargin(0.16);
    mapCFneg[imf0]->Draw("colz");
    lednickyReallySucks[imf0]->Draw("3 same");
    lednicky->Print(
        Form("%s/LednickyFails_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    mapCFneg[imf0]->Write(Form("LednickyFails_%.1f", IMf0invCont[imf0]));
    lednickyReallySucks[imf0]->Write(
        Form("LednickyFailsCountour_%.1f", IMf0invCont[imf0]));

    dBinding->Clear();
    dBinding->SetRightMargin(0.16);
    mapBinding[imf0]->Draw("colz");
    dBinding->Print(
        Form("%s/Binding_%.1f.pdf", inputDir.Data(), IMf0invCont[imf0]));
    mapBinding[imf0]->Write(Form("Binding_%.1f", IMf0invCont[imf0]));
    dBinding->Write(Form("BindingPlot_%.1f", IMf0invCont[imf0]));
  }
  dWorst->Print(Form("%s/ExclusionWorst.gif++", inputDir.Data()));
  dDef->Print(Form("%s/ExclusionDefault.gif++", inputDir.Data()));
  dBest->Print(Form("%s/ExclusionBest.gif++", inputDir.Data()));

  fileOut->Close();
  file->Close();
}

int main(int argc, char *argv[]) {
  PlotSigmaExclusion(argv[1]);
  return 0;
}

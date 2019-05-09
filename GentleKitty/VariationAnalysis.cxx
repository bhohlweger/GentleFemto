/*
 * DrawPP.cxx
 *
 *  Created on: May 8, 2019
 *      Author: schmollweger
 */

#include "VariationAnalysis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TError.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include "TCanvas.h"
VariationAnalysis::VariationAnalysis(const char* histname, const int nVars,
                                     const int nFitVars)
    : fInFile(nullptr),
      fHistname(histname),
      fnDataVars(nVars),
      fnFitVars(nFitVars),
      fRadMean(0),
      fRadSystUp(0),
      fRadSystDown(0),
      fRadStat(0),
      fCk(),
      fModel(nullptr),
      fRadiusDist(nullptr) {
  DreamPlot::SetStyle();
}

VariationAnalysis::~VariationAnalysis() {

}

void VariationAnalysis::ReadFitFile(TString FileName) {
  fInFile = TFile::Open(FileName, "READ");
  if (!fInFile) {
    Error("ReadFitFile", "No input file");
    return;
  }
  TNtuple *resultTuple = (TNtuple*) fInFile->Get("ntResult");
  if (!resultTuple) {
    Error("ReadFitFile", "No Result tuple rip. \n");
    return;
  } else {
    if (fRadiusDist) {
      delete fRadiusDist;
    }
    resultTuple->Draw("Radius_pp>>RadDist");
    fRadiusDist = (TH1D*) gROOT->FindObject("RadDist");
    float radMin = 0.9 * (fRadiusDist->GetMean() - fRadiusDist->GetRMS());
    float radMax = 1.1 * (fRadiusDist->GetMean() + fRadiusDist->GetRMS());
    delete fRadiusDist;
    fRadiusDist = new TH1D("RadDist", "RadDist", 200, radMin, radMax);
    resultTuple->Draw("Radius_pp>>RadDist");
    resultTuple->Draw("RadiusErr_pp>>RadStat");
    TH1F* statErr = (TH1F*) gROOT->FindObject("RadStat");
    fRadStat = statErr->GetMean();
  }
  TFile* tmpFile = TFile::Open(TString::Format("%s/tmp.root", gSystem->pwd()),
                               "RECREATE");
  if (!tmpFile) {
    Error("ReadFitFile", "No Tmp file");
    return;
  }
  TNtuple* Fits = new TNtuple("fitCurves", "fitCurves", "kstar:modelValue");
  Fits->Write();
  TGraph* refGraph;
  for (int iVars = 0; iVars < fnDataVars + 1; ++iVars) {
    TString histname = TString::Format("%s%iMeV_0", fHistname, iVars);
    TH1F* histo = (TH1F*) fInFile->Get(histname.Data());
    if (!histo) {
      TString OutputError = TString::Format("Histogram (%s) missing, rip",
                                            histname.Data());
      Error("ReadFitFile", OutputError.Data());
    } else {
      fCk.push_back(histo);
    }
    TList* outList = (TList*) fInFile->Get(TString::Format("Out%i", iVars));
    if (!outList) {
      TString OutputError = TString::Format(
          "Outlist %s not available", TString::Format("Out%i", iVars).Data())
          .Data();
      Error("ReadFitFile", OutputError.Data());
    }
    //loop over all variations for on fit
    for (int iFitVar = 1; iFitVar < fnFitVars; iFitVar++) {
      TString folderName = TString::Format("Graph_Var_%i_iter_%i", iVars,
                                           iFitVar);
      TList* GraphList = (TList*) outList->FindObject(folderName.Data());
      if (!GraphList) {
        TString OutputError = TString::Format("GraphList %s not available",
                                              folderName.Data()).Data();
        Error("ReadFitFile", OutputError.Data());
        return;
      } else {
        TString GraphName = TString::Format("FitResult_%i", iFitVar);
        TGraph* Graph = (TGraph*) GraphList->FindObject(GraphName.Data());
        if (!Graph) {
          Error("ReadFitFile", GraphName.Data());
          return;
        } else {
          double x, y;
          if (iVars == 0 && iFitVar == 1) {
            refGraph = Graph;
          }
          for (int iPnt = 0; iPnt < Graph->GetN(); ++iPnt) {
            Graph->GetPoint(iPnt, x, y);
            Fits->Fill(x, y);
          }
        }
      }
    }
  }
  fModel = EvaluateCurves(Fits, refGraph);
  fModel->SetName("Model");
  tmpFile->cd();
  fModel->Write();
  tmpFile->Write();
  tmpFile->Close();
}

TGraphErrors * VariationAnalysis::EvaluateCurves(TNtuple * tuple,
                                                 TGraph * RefGraph) {
  //Ref Graph is just any fit graph to have the correct x values for the tuple.
  //user needs to delete grOut.
  TGraphErrors* grOut = new TGraphErrors();
  double kVal, Ck;
  for (int ikstar = 0; ikstar < RefGraph->GetN(); ++ikstar) {
    RefGraph->GetPoint(ikstar, kVal, Ck);
    tuple->Draw(Form("modelValue >> h%i", ikstar),
                Form("std::abs(kstar - %.3f) < 1e-3", kVal));
    TH1F* hist = (TH1F*) gROOT->FindObject(Form("h%i", ikstar));

    double binLow = hist->GetXaxis()->GetBinLowEdge(
        hist->FindFirstBinAbove(0.1, 1));
    double binUp = hist->GetXaxis()->GetBinUpEdge(
        hist->FindLastBinAbove(0.1, 1));
    double DeltaCoulomb = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    double DefaultVal = (binUp + binLow) / 2.;
    grOut->SetPoint(ikstar, kVal, DefaultVal);
    grOut->SetPointError(ikstar, 0, DeltaCoulomb);
    delete hist;
  }
  return grOut;
}
void VariationAnalysis::EvalRadius() {
  fRadMean = fRadiusDist->GetMean();
  int n = fRadiusDist->GetXaxis()->GetNbins();

  auto histRadCumulative = fRadiusDist->GetCumulative();
  histRadCumulative->Scale(1. / (double) fRadiusDist->GetEntries());
  auto c1 = new TCanvas("c4", "c5");
  histRadCumulative->Draw("");
  c1->SaveAs(Form("%s/cumulative.pdf", gSystem->pwd()));
  auto medianBin = histRadCumulative->FindFirstBinAbove(0.5, 1);
  std::cout << "medianBin: " << medianBin << " Median: "
            << histRadCumulative->GetBinCenter(medianBin) << std::endl;
  int binMin = 0;
  int binMax = 0;
  for (int iBin = 0; iBin < histRadCumulative->GetNbinsX(); iBin++) {
    if (binMin == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) - 0.34)) {
      binMin = iBin;
    }
    if (binMax == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) + 0.34)) {
      binMax = iBin;
      break;
    }
  }
  auto radMin = histRadCumulative->GetXaxis()->GetBinCenter(binMin);
  auto radMax = histRadCumulative->GetXaxis()->GetBinCenter(binMax);
  std::cout << "radMin: " << radMin << " radMax: " << radMax << std::endl;
  fRadSystUp = radMax - fRadMean;
  fRadSystDown = fRadMean - radMin;

  auto *canRad2 = new TCanvas();
  canRad2->cd();
  fRadiusDist->Rebin(2);
  DreamPlot::SetStyleHisto(fRadiusDist, 20, 1);
  fRadiusDist->GetXaxis()->SetTitle("Core Radius (fm)");
  fRadiusDist->GetYaxis()->SetTitle("Number of Entries");
  fRadiusDist->Draw();
  auto histRadLimits = (TH1F*) fRadiusDist->Clone("histRadLimits");
  histRadLimits->Reset();
  for (int i = 0; i < fRadiusDist->GetNbinsX(); ++i) {
    if (fRadiusDist->GetBinCenter(i) < radMin
        || fRadiusDist->GetBinCenter(i) > radMax)
      continue;
    histRadLimits->SetBinContent(i, fRadiusDist->GetBinContent(i));
  }
  histRadLimits->SetFillColor(kGray + 1);
  histRadLimits->Draw("same");

  auto lineDefault = new TLine(fRadMean, 0, fRadMean,
                               fRadiusDist->GetMaximum());
  lineDefault->SetLineColor(kRed + 2);
  lineDefault->SetLineWidth(2);
  lineDefault->Draw("same");

  auto lineLow = new TLine(radMin, 0, radMin, fRadiusDist->GetMaximum());
  lineLow->SetLineColor(kGreen + 2);
  lineLow->SetLineWidth(2);
  lineLow->Draw("same");

  auto lineUp = new TLine(radMax, 0, radMax, fRadiusDist->GetMaximum());
  lineUp->SetLineColor(kGreen + 2);
  lineUp->SetLineWidth(2);
  lineUp->Draw("same");

  canRad2->SaveAs(TString::Format("%s/radius.pdf", gSystem->pwd()));
}

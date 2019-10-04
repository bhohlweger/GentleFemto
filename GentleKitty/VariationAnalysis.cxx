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
#include "TEntryList.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include "TCanvas.h"
VariationAnalysis::VariationAnalysis(const char* histname)
    : fInFile(nullptr),
      fSelector(""),
      fHistname(histname),
      fnDataVarsStart(0),
      fnDataVarsEnd(0),
      fnFitVarsStart(0),
      fnFitVarsEnd(0),
      fRadMean(0),
      fRadSystUp(0),
      fRadSystDown(0),
      fRadStat(0),
      fCk(),
      fModel(nullptr),
      fRadiusDist(nullptr),
      fRadiusErrDist(nullptr),
      fDeviationByBin(nullptr) {
  DreamPlot::SetStyle();
}

VariationAnalysis::~VariationAnalysis() {

}

void VariationAnalysis::ReadFitFile(TString FileName) {
  static int imT = 0;
  TFile* tmpFile = TFile::Open(
      TString::Format("%s/tmp_%u.root", gSystem->pwd(), imT++), "RECREATE");
  if (!tmpFile) {
    Error("ReadFitFile", "No Tmp file");
    return;
  }
  TNtuple* Fits = new TNtuple("fitsCurves", "fitsCurves", "kstar:modelValue");
  Fits->Write();

  fInFile = TFile::Open(FileName, "READ");
  if (!fInFile) {
    Error("ReadFitFile", "No input file");
    return;
  }
  fInFile->cd();
  TTree *resultTree = (TTree*) fInFile->Get("ppTree");
  if (!resultTree) {
    Error(
        "ReadFitFile",
        TString::Format("No Result tuple in %s .. .rip. \n", FileName.Data()));
    return;
  } else {
    if (fRadiusDist) {
      delete fRadiusDist;
    }
    if (fRadiusErrDist) {
      delete fRadiusErrDist;
    }
    resultTree->Draw("Radius>>RadDist");
    fRadiusDist = (TH1D*) gROOT->FindObject("RadDist");
    float radMin = fRadiusDist->GetXaxis()->GetXmin(); //0.9 * (fRadiusDist->GetMean() - fRadiusDist->GetRMS());
    float radMax = fRadiusDist->GetXaxis()->GetXmax(); //1.1 * (fRadiusDist->GetMean() + fRadiusDist->GetRMS());
    int nRadBims = (radMax - radMin)/(float)(0.01);
    delete fRadiusDist;
    fRadiusDist = new TH1D("RadDist", "RadDist", nRadBims, radMin, radMax);

    resultTree->Draw("RadiusErr>>RadErrDist");
    fRadiusErrDist = (TH1D*) gROOT->FindObject("RadErrDist");
    float radErrMin = fRadiusErrDist->GetXaxis()->GetXmin(); //0.9 * (fRadiusErrDist->GetMean() - fRadiusErrDist->GetRMS());
    float radErrMax = fRadiusErrDist->GetXaxis()->GetXmax(); //1.1* (fRadiusErrDist->GetMean() + fRadiusErrDist->GetRMS());
    int nRadErrBims = (radErrMax-radErrMin)/(float)0.01;
    delete fRadiusErrDist;
    fRadiusErrDist = new TH1D("RadErrDist", "RadErrDist", nRadErrBims, radErrMin, radErrMax);

    int before = resultTree->GetEntriesFast();
    resultTree->Draw(">>myList", fSelector, "entrylist");
    std::cout << fSelector.GetName() << '\t' << fSelector.GetTitle() << std::endl;
    TEntryList *entrList=(TEntryList*)gDirectory->Get("myList");
    if (!entrList) {
      Error("ReadFitFile", "Entry list missing \n");
      return;
    }

    const int nEntries = entrList->GetN();
    std::cout << "Entries Before: " << before << " and after "
              << nEntries << " meaning a delta of " << before - nEntries
              << std::endl;
    if (nEntries == 0) {
      return;
    }
    int lastDataVar = -1;

    unsigned int dataVar;
    float chisq;
    float radius;
    float radiusErr;
    TH1F* CFHisto = nullptr;
    TGraph* FitCurve = nullptr;
    TGraph* refGraph = nullptr;
    fInFile->cd();
    resultTree->SetBranchAddress("DataVarID", &dataVar);
    resultTree->SetBranchAddress("chiSqNDF", &chisq);
    resultTree->SetBranchAddress("Radius", &radius);
    resultTree->SetBranchAddress("RadiusErr", &radiusErr);
    resultTree->SetBranchAddress("CorrHist", &CFHisto);
    resultTree->SetBranchAddress("FitResult", &FitCurve);
    for (int iEntr = 0; iEntr < nEntries; ++iEntr) {
      resultTree->GetEntry(entrList->GetEntry(iEntr));
      fRadiusDist->Fill(radius);
      fRadiusErrDist->Fill(radiusErr);
      if (!refGraph) {
        refGraph = new TGraph(FitCurve->GetN(), FitCurve->GetX(),
                              FitCurve->GetY());
      }
      if (dataVar != lastDataVar) {
        fCk.push_back(CFHisto);
        lastDataVar = dataVar;
      }
      double x, y;
      for (int iPnt = 0; iPnt < FitCurve->GetN(); ++iPnt) {
        FitCurve->GetPoint(iPnt, x, y);
        Fits->Fill(x, y);
      }
    }
    fModel = EvaluateCurves(Fits, refGraph);
    fModel->SetName("Model");
    fDeviationByBin = DeviationByBin(fCk.at(0), fModel);
    fDeviationByBin->SetName("DeviationPerBin");
    tmpFile->cd();
    fModel->Write();
    fDeviationByBin->Write();
    tmpFile->Write();
    tmpFile->Close();
  }
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

TGraphErrors* VariationAnalysis::DeviationByBin(TH1F* RefHist,
                                                TGraphErrors* model) {
  TGraphErrors* grOut = new TGraphErrors();
  double kVal, Ck;
  for (int ikstar = 0; ikstar < model->GetN(); ++ikstar) {
    model->GetPoint(ikstar, kVal, Ck);
    int iDataBin = RefHist->FindBin(kVal);
    if (std::abs(kVal - RefHist->GetBinCenter(iDataBin)) > 1e-3) {
      TString OutError =
          TString::Format(
              "Deviation between Graph & Histogram of %.3f. Someone should look into this \n",
              std::abs(kVal - RefHist->GetBinCenter(iDataBin)));
      Error("DeviationByBin", OutError.Data());
      return nullptr;
    }
    double CkErr = model->GetErrorY(ikstar);
    double CkData = RefHist->GetBinContent(iDataBin);
    double CkErrStatData = RefHist->GetBinError(iDataBin);

    double deviation = (Ck - CkData) / CkErrStatData;
    double err = CkErr / CkErrStatData;

    grOut->SetPoint(ikstar, kVal, deviation);
    grOut->SetPointError(ikstar, 0, err);
  }
  return grOut;
}

void VariationAnalysis::EvalRadius(const char* bin) {
  fRadMean = fRadiusDist->GetMean();
  fRadStat = fRadiusErrDist->GetMean();
  int n = fRadiusDist->GetXaxis()->GetNbins();

  auto histRadCumulative = fRadiusDist->GetCumulative();
  histRadCumulative->Scale(1. / (double) fRadiusDist->GetEntries());
  auto c1 = new TCanvas("c4", "c5");
  histRadCumulative->Draw("");
  c1->SaveAs(Form("%s/cumulative%s.pdf", gSystem->pwd(), bin));
  delete c1;
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

  canRad2->SaveAs(TString::Format("%s/radius%s.pdf", gSystem->pwd(), bin));
}

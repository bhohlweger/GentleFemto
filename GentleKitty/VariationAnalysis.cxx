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
      fRadiusDist(nullptr),
      fDeviationByBin(nullptr) {
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
  TGraph* refGraph = nullptr;
  for (int iVars = 0; iVars < fnDataVars + 1; ++iVars) {
    TString dirName = TString::Format("Out%i", iVars);
    TDirectoryFile* dir = (TDirectoryFile*) fInFile->FindObjectAny(
        dirName.Data());
    if (!dir) {
      TString OutputError = TString::Format("No directory named %s",
                                            dirName.Data());
      Error("ReadFitFile", OutputError.Data());
    }
    TString histname = TString::Format("%s%iMeV_0", fHistname, iVars);
    TH1F* histo = (TH1F*) dir->FindObjectAny(histname.Data());
    if (!histo) {
      TString OutputError = TString::Format("Histogram (%s) missing, rip",
                                            histname.Data());
      Error("ReadFitFile", OutputError.Data());
    } else {
      fCk.push_back(histo);
    }
    //loop over all variations for on fit
    for (int iFitVar = 1; iFitVar < fnFitVars; iFitVar++) {
//      resultTuple->Draw(
//          "PolBaseLine>>myBaseLine",
//          Form("std::abs(NumIter-%u)<1e-3&&std::abs(IterID-%u)<1e-3", iVars,
//               iFitVar));
//      TH1F* polBL = (TH1F*) gROOT->FindObject("myBaseLine");
//      if (std::abs(polBL->GetMean()-2)<0.1) {
//        continue;
//        delete polBL;
//      }
//      delete polBL;
      resultTuple->Draw(
          "chisqPerndf>>chisq",
          Form("std::abs(NumIter-%u)<1e-3&&std::abs(IterID-%u)<1e-3", iVars,
               iFitVar));
      TH1F* chiSq = (TH1F*) gROOT->FindObject("chisq");
      if (chiSq->GetMean() > 15.) {
        Warning(
            "ReadFitFile",
            Form("Chisq (%.1f) larger than 15, ignoring fit",
                 chiSq->GetMean()));
        delete chiSq;
        continue;
      }
      delete chiSq;
      TString folderName = TString::Format("Graph_Var_%i_iter_%i", iVars,
                                           iFitVar);
      TList* GraphList = (TList*) dir->FindObjectAny(folderName.Data());
      if (!GraphList) {
        TString OutputError = TString::Format("GraphList %s not available",
                                              folderName.Data()).Data();
        Error("ReadFitFile", OutputError.Data());
        return;
      } else {
        TString GraphName = TString::Format("Graph_Var_%i_Iter_%i", iVars,
                                            iFitVar);
        TGraph* Graph = (TGraph*) GraphList->FindObject(GraphName.Data());
        if (!Graph) {
          GraphList->ls();
          Error("ReadFitFile", folderName.Data());
          return;
        } else {
          double x, y;
          if (!refGraph) {
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
  fDeviationByBin = DeviationByBin(fCk.at(0), fModel);
  fDeviationByBin->SetName("DeviationPerBin");
  tmpFile->cd();
  fModel->Write();
  fDeviationByBin->Write();
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
    double DeltaCoulomb = TMath::Abs((binLow - binUp))/TMath::Sqrt(12);
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
  int n = fRadiusDist->GetXaxis()->GetNbins();

  auto histRadCumulative = fRadiusDist->GetCumulative();
  histRadCumulative->Scale(1. / (double) fRadiusDist->GetEntries());
  auto c1 = new TCanvas("c4", "c5");
  histRadCumulative->Draw("");
  c1->SaveAs(Form("%s/cumulative%s.pdf", gSystem->pwd(),bin));
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

  canRad2->SaveAs(TString::Format("%s/radius%s.pdf", gSystem->pwd(),bin));
}

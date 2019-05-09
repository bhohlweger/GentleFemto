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

VariationAnalysis::VariationAnalysis(const char* histname, const int nVars,
                                     const int nFitVars)
    : fInFile(nullptr),
      fHistname(histname),
      fnDataVars(nVars),
      fnFitVars(nFitVars),
      fnkMin(0),
      fnModelBins(0),
      fdkstar(0),
      fCk(nullptr),
      fFits(nullptr) {

}

VariationAnalysis::~VariationAnalysis() {

}

void VariationAnalysis::ReadFitFile(TString FileName) {
  fInFile = TFile::Open(FileName, "READ");
  if (fFits) {
    delete fFits;
  }
  fFits = new TNtuple("fitCurves", "fitCurves", "kstar:modelValue");
  for (int iVars = 0; iVars < fnDataVars + 1; ++iVars) {
    TString histname = TString::Format("%s%uMeV_0");
    TH1F* histo = (TH1F*) fInFile->Get(histname.Data());
    if (!histo) {
      Error("ReadFitFile",
            TString::Format("Histogram (%s) missing, rip", histname.Data()));
    } else {
      fCk.push_back(histo);
    }
    if (iVars == 0) {
      fnkMin = histo->GetXaxis()->GetXmin();
      fdkstar = histo->GetBinWidth(1);
    }
    TList* outList = (TList*) fInFile->Get(TString::Format("Out%u", iVars));
    if (!outList) {
      Error(
          "ReadFitFile",
          TString::Format("Outlist %s not available",
                          TString::Format("Out%u", iVars)));
    }
    //loop over all variations for on fit
    for (int iFitVar = 1; iFitVar < fnFitVars; iFitVar++) {
      TString folderName = TString::Format("Graph_Var%u_iter_%u");
      TList* GraphList = (TList*) outList->FindObject(folderName.Data());
      if (!GraphList) {
        Error("ReadFitFile",
              TString::Format("GraphList %s not available", folderName.Data()));
      } else {
        TString GraphName = TString("FitResult_%u", iFitVar);
        TGraph* Graph = (TGraph*) GraphList->FindObject(GraphName.Data());
        double x, y;
        if (iVars == 0 && iFitVar == 1) {
          fnModelBins = Graph->GetN();
        }
        for (int iPnt = 0; iPnt < Graph->GetN(); ++iPnt) {
          Graph->GetPoint(iPnt, x, y);
          fFits->Fill(x, y);
        }
      }
    }
  }
}

TGraphErrors* VariationAnalysis::ModelFitBands() {
  return EvaluateCurves(fFits, fnModelBins, fnkMin, fdkstar);
}

TGraphErrors* VariationAnalysis::EvaluateCurves(TNtuple* tuple, const int nBins,
                                                const int kMin,
                                                const int dkStar) {
  //user needs to delete grOut.
  TGraphErrors* grOut = new TGraphErrors();
  for (int ikstar = 0; ikstar < nBins; ++ikstar) {
    double kVal = kMin + ikstar * dkStar;
    tuple->Draw("modelValue >> h", Form("std::abs(kstar - %.3f) < 1e-3", kVal));
    TH1F* hist = (TH1F*) gROOT->FindObject("h");
    double binLow = hist->GetXaxis()->GetBinLowEdge(
        hist->FindFirstBinAbove(0.1, 1));
    double binUp = hist->GetXaxis()->GetBinUpEdge(
        hist->FindLastBinAbove(0.1, 1));
    double DeltaCoulomb = (binLow - binUp) / TMath::Sqrt(12);
    double DefaultVal = (binUp + binLow) / 2.;
    //std::cout << ikstar << " " << kVal << " " << CFval << " " << DefaultVal << "\n";
    grOut->SetPoint(ikstar, kVal, DefaultVal);
    grOut->SetPointError(ikstar, 0, DeltaCoulomb);
    delete hist;
  }
  return grOut;
}

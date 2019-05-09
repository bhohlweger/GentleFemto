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
#include <iostream>
VariationAnalysis::VariationAnalysis(const char* histname, const int nVars,
                                     const int nFitVars)
    : fInFile(nullptr),
      fHistname(histname),
      fnDataVars(nVars),
      fnFitVars(nFitVars),
      fnRadBins(0),
      fRadMin(0),
      fRadMax(0),
      fRadStat(0),
      fCk(),
      fModel(nullptr),
      fRadiusDist(nullptr) {

}

VariationAnalysis::~VariationAnalysis() {

}

void VariationAnalysis::ReadFitFile(TString FileName) {
  fInFile = TFile::Open(FileName, "READ");
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
  TNtuple *resultTuple = (TNtuple*) fInFile->Get("ntResult");
  if (!resultTuple) {
    Error("ReadFitFile", "No Result tuple rip. \n");
    return;
  } else {
    if (fRadiusDist) {
      delete fRadiusDist;
    }
    fRadiusDist = new TH1F("RadDist", "RadDist", fnRadBins, fRadMin, fRadMax);
    resultTuple->Draw("Radius_pp>>RadDist");
    TH1F* statErr = new TH1F("RadStat", "RadStat", 500, 0, 2.);
    resultTuple->Draw("RadiusErr_pp>>RadStat");
    fRadStat = statErr->GetMean();
  }
}

TGraphErrors* VariationAnalysis::EvaluateCurves(TNtuple* tuple,
                                                TGraph* RefGraph) {
  //Ref Graph is just any fit graph to have the correct x values for the tuple.
  //user needs to delete grOut.
  TGraphErrors* grOut = new TGraphErrors();
  double kVal, Ck;
  for (int ikstar = 0; ikstar < RefGraph->GetN(); ++ikstar) {
    RefGraph->GetPoint(ikstar, kVal, Ck);
    tuple->Draw("modelValue >> h", Form("std::abs(kstar - %.3f) < 1e-3", kVal));
    TH1F* hist = (TH1F*) gROOT->FindObject("h");
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

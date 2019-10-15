/*
 * DrawPP.cxx
 *
 *  Created on: May 8, 2019
 *      Author: schmollweger
 */

#include "VariationAnalysispAL.h"
#include "TFile.h"
#include "TGraph.h"
#include "TError.h"
#include "TEntryList.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include "TCanvas.h"
VariationAnalysispAL::VariationAnalysispAL(const char* histname)
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

VariationAnalysispAL::~VariationAnalysispAL() {

}

void VariationAnalysispAL::ReadFitFile(TString FileName) {
  TFile* tmpFile = TFile::Open(
      TString::Format("%s/tmp_pAL_0.root", gSystem->pwd()), "RECREATE");
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
  TTree *resultTree = (TTree*) fInFile->Get("SysTree");
  if (!resultTree) {
    Error(
        "ReadFitFile",
        TString::Format("No Result tuple in %s .. .rip. \n", FileName.Data()));
    return;
  }


    int before = resultTree->GetEntriesFast();
    resultTree->Draw(">>myList", fSelector, "entrylist");
    std::cout << fSelector.GetName() << '\t' << fSelector.GetTitle() << std::endl;

    TEntryList *entrList=(TEntryList*)gDirectory->Get("myList");
    if (!entrList) {
      Error("ReadFitFile", "Entry list missing \n");
      return;
    }else {

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
    TH1F* CFHisto_pAL = nullptr;
    TGraph* FitCurve_pAL = nullptr;
    TGraph* refGraph_pAL = nullptr;


    fInFile->cd();
    resultTree->SetBranchAddress("Data_pAL_VarID", &dataVar);
    resultTree->SetBranchAddress("chiSqNDF", &chisq);
    resultTree->SetBranchAddress("CorrHist_pAL", &CFHisto_pAL);
    resultTree->SetBranchAddress("FitResultTotal_pAL", &FitCurve_pAL);


    for (int iEntr = 0; iEntr < nEntries; ++iEntr) {
      resultTree->GetEntry(entrList->GetEntry(iEntr));

      if (!refGraph_pAL) {
        refGraph_pAL = new TGraph(FitCurve_pAL->GetN(), FitCurve_pAL->GetX(),
                              FitCurve_pAL->GetY());
      }
      if (dataVar != lastDataVar) {
        fCk.push_back(CFHisto_pAL);
        lastDataVar = dataVar;
      }
      double x, y;
      for (int iPnt = 0; iPnt < FitCurve_pAL->GetN(); ++iPnt) {
        FitCurve_pAL->GetPoint(iPnt, x, y);
        Fits->Fill(x, y);
      }
    }


    fModel = EvaluateCurves(Fits, refGraph_pAL);
     fModel->SetName("Model_pAL");
     fDeviationByBin = DeviationByBin(fCk.at(0), fModel);
     fDeviationByBin->SetName("DeviationPerBin_pAL");

    tmpFile->cd();
    refGraph_pAL->Write();
    fModel->Write();
    fDeviationByBin->Write();
    tmpFile->Write();
    tmpFile->Close();
  }
}

TGraphErrors * VariationAnalysispAL::EvaluateCurves(TNtuple * tuple,
                                                 TGraph * RefGraph) {
  //Ref Graph is just any fit graph to have the correct x values for the tuple.
  //user needs to delete grOut.
  TGraphErrors* grOut = new TGraphErrors();
  double kVal, Ck;
  for (int ikstar = 0; ikstar < RefGraph->GetN(); ++ikstar) {
    RefGraph->GetPoint(ikstar, kVal, Ck);
    tuple->Draw("modelValue >> h",Form("TMath::Abs(kstar - %.3f) < 1e-3", kVal));

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

TGraphErrors* VariationAnalysispAL::DeviationByBin(TH1F* RefHist,
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

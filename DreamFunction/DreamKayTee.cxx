/*
 * DreamKayTee.cxx
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#include "DreamKayTee.h"
#include "TCanvas.h"
#include <iostream>
DreamKayTee::DreamKayTee()
    : fIskT(true),
      fKayTeeBins(),
      fNKayTeeBins(0),
      fAveragekT(0),
      fCFPart(nullptr),
      fSum(nullptr),
      fNormleft(0),
      fNormright(0) {
  fSEkT[0] = nullptr;
  fSEkT[1] = nullptr;
  fMEkT[0] = nullptr;
  fMEkT[1] = nullptr;
}

DreamKayTee::~DreamKayTee() {
  // TODO Auto-generated destructor stub
}

void DreamKayTee::ObtainTheCorrelationFunction(const char* outFolder,
                                               const char* prefix,
                                               const char* pair) {
  const char* variable = (fIskT) ? "kT" : "mT";
  const int nBins = (int) fKayTeeBins.size();
  if (fSEkT[0]) {
    fKayTeeBins.push_back(fSEkT[0]->GetYaxis()->GetXmax());
    fCFPart = new DreamPair**[2];
    TString PartName[2] = { "Part", "AntiPart" };
    for (int iPart = 0; iPart < 2; ++iPart) {
      fCFPart[iPart] = new DreamPair*[fNKayTeeBins];
      for (int ikT = 0; ikT < fNKayTeeBins; ++ikT) {
        TString PairName = Form("%s_%s_%i", PartName[iPart].Data(), variable,
                                ikT);
        fCFPart[iPart][ikT] = new DreamPair(PairName.Data(), fNormleft,
                                            fNormright);
        int kTminBin = fSEkT[iPart]->GetYaxis()->FindBin(fKayTeeBins[ikT]);
        int kTmaxBin = fSEkT[iPart]->GetYaxis()->FindBin(fKayTeeBins[ikT + 1]);
//        std::cout << kTminBin << '\t' << kTmaxBin << std::endl;
        DreamDist* kTDist = new DreamDist();
        TString SEkTBinName = Form("SE%s", PairName.Data());
        kTDist->SetSEDist(
            (TH1F*) fSEkT[iPart]->ProjectionX(SEkTBinName.Data(), kTminBin + 1,
                                              kTmaxBin),
            "");
        TString MEkTBinName = Form("ME%s", PairName.Data());
        kTDist->SetMEDist(
            (TH1F*) fMEkT[iPart]->ProjectionX(MEkTBinName.Data(), kTminBin + 1,
                                              kTmaxBin),
            "");
        fCFPart[iPart][ikT]->SetPair(kTDist);
        fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(), 2);
      }
    }
    this->AveragekT();
    fSum = new DreamCF*[fNKayTeeBins];
    TString outname = outFolder;
    outname += "/CFOutputALL_";
    outname += variable;
    outname += "_";
    outname += pair;
    outname += "_";
    outname += prefix;
    outname += ".root";
    TFile* allCFsOut = TFile::Open(outname.Data(), "RECREATE");
    if (fAveragekT) {
      fAveragekT->Write(Form("Average%s", variable));
    }
    for (int ikT = 0; ikT < fNKayTeeBins; ++ikT) {
      fSum[ikT] = new DreamCF();
      fSum[ikT]->SetPairs(fCFPart[0][ikT], fCFPart[1][ikT]);
      fSum[ikT]->GetCorrelations();
      std::vector<TH1F*> CFs = fSum[ikT]->GetCorrelationFunctions();
      TString outfileName = outFolder;
      outfileName += "/CFOutput_";
      outfileName += variable;
      outfileName += "_";
      outfileName += pair;
      outfileName += "_";
      outfileName += prefix;
      outfileName += "_";
      outfileName += ikT;
      outfileName += ".root";
      allCFsOut->cd();
      for (auto &it : CFs) {
        TString OutName = Form("%s_%sBin_%i", it->GetName(), variable, ikT);
//        std::cout << OutName.Data() << std::endl;
        it->Write(OutName.Data());
      }
      fSum[ikT]->WriteOutput(outfileName.Data());
    }
  } else {
    std::cout << "No SE " << variable << " histogram with index 0 \n";
  }
  return;
}

void DreamKayTee::AveragekT() {
  const char* variable = (fIskT) ? "kT" : "mT";
  fAveragekT = new TGraphErrors();
  TH2F* kTkStar = (TH2F*) fSEkT[0]->Clone(Form("%skStarForAverage", variable));
  kTkStar->Add(fSEkT[1]);
  TH1F* kTProjection = (TH1F*) kTkStar->ProjectionY(Form("%sDist", variable), 0,
                                                    -1, "e");
  for (int ikT = 0; ikT < fNKayTeeBins; ++ikT) {
    int binLow = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT));
    int binUp = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT + 1));
    kTProjection->GetXaxis()->SetRange(binLow, binUp);
    fAveragekT->SetPoint(ikT, ikT + 1, kTProjection->GetMean());
    fAveragekT->SetPointError(ikT, 0, kTProjection->GetRMS());
  }
}

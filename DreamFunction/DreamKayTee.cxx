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
    : fKayTeeBins(),
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

void DreamKayTee::ObtainTheCorrelationFunction(const char* outFolder) {
  const int nBins = (int) fKayTeeBins.size();
  if (fSEkT[0]) {
    fKayTeeBins.push_back(fSEkT[0]->GetYaxis()->GetXmax());
    fCFPart = new DreamPair**[2];
    TString PartName[2] = { "Part", "AntiPart" };
    for (int iPart = 0; iPart < 2; ++iPart) {
      fCFPart[iPart] = new DreamPair*[fNKayTeeBins];
      for (int ikT = 0; ikT < fNKayTeeBins; ++ikT) {
        TString PairName = Form("%s_kT_%i", PartName[iPart].Data(), ikT);
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
        fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(),2);
      }
    }
    this->AveragekT();
    fSum = new DreamCF*[fNKayTeeBins];
    TFile* allCFsOut = TFile::Open(Form("%s/CFOutputALLkT_pp.root", outFolder),"RECREATE");
    if (fAveragekT) {
      fAveragekT->Write("AveragekT");
    }
    for (int ikT = 0; ikT<fNKayTeeBins; ++ikT) {
      fSum[ikT] = new DreamCF();
      fSum[ikT]->SetPairs(fCFPart[0][ikT],fCFPart[1][ikT]);
      fSum[ikT]->GetCorrelations();
      std::vector<TH1F*> CFs = fSum[ikT]->GetCorrelationFunctions();
      TString outfileName = Form("%s/CFOutput_pp_ikT%i.root", outFolder, ikT);
      allCFsOut->cd();
      for (auto &it : CFs) {
        TString OutName = Form("%s_kTBin_%i",it->GetName(),ikT);
//        std::cout << OutName.Data() << std::endl;
        it->Write(OutName.Data());
      }
      fSum[ikT]->WriteOutput(outfileName.Data());

    }
  } else {
    std::cout << "No SE kT histogram with index 0 \n";
  }
  return;
}

void DreamKayTee::AveragekT() {
  fAveragekT=  new TGraphErrors();
  TH2F* kTkStar = (TH2F*)fSEkT[0]->Clone("kTkStarForAverage");
  kTkStar->Add(fSEkT[1]);
  TH1F* kTProjection = (TH1F*)kTkStar->ProjectionY("kDist",0,-1,"e");
  for (int ikT = 0; ikT < fNKayTeeBins; ++ikT) {
    int binLow = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT));
    int binUp = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT+1));
    kTProjection->GetXaxis()->SetRange(binLow,binUp);
    fAveragekT->SetPoint(ikT,ikT+1,kTProjection->GetMean());
    fAveragekT->SetPointError(ikT,0,kTProjection->GetRMS());
  }
}

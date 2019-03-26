/*
 * DreamKayTee.cxx
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#include "DreamKayTee.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
DreamKayTee::DreamKayTee()
    : fIskT(true),
      fKayTeeBins(),
      fNKayTeeBins(0),
      fAveragekT(0),
      fCFPart(nullptr),
      fSum(nullptr),
      fNormleft(0),
      fNormright(0),
      fSEMEReweighting(nullptr),
      fSEMEReweightingMeV(nullptr) {
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
  if (!fSEMEReweighting) {
    std::cout << "No Reweighting set. Call SetSEMEReweightingRatio first! \n";
  } else {
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
          int kTmaxBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT + 1]);
//        std::cout << kTminBin << '\t' << kTmaxBin << std::endl;
          DreamDist* kTDist = new DreamDist();
          TString SEkTBinName = Form("SE%s", PairName.Data());
          kTDist->SetSEDist(
              (TH1F*) fSEkT[iPart]->ProjectionX(SEkTBinName.Data(),
                                                kTminBin + 1, kTmaxBin, "e"),
              "");
          TString MEkTBinName = Form("ME%s", PairName.Data());
          kTDist->SetMEDist(
              (TH1F*) fMEkT[iPart]->ProjectionX(MEkTBinName.Data(),
                                                kTminBin + 1, kTmaxBin, "e"),
              "");
          fCFPart[iPart][ikT]->SetPair(kTDist);
          fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(), 2);
        }
      }
      this->AveragekT(pair);
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
      for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
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
          TString HistName = it->GetName();
          for (int iBin = 1; iBin < it->GetNbinsX(); ++iBin) {
            float binCont = it->GetBinContent(iBin);
            float binErr = it->GetBinError(iBin);
            float binCent = it->GetBinCenter(iBin);
            float corrFactor = 0;
            if (HistName.Contains("MeV")) {
              corrFactor = fSEMEReweightingMeV->GetBinContent(
                  fSEMEReweightingMeV->FindBin(binCent));
            } else {
              corrFactor = fSEMEReweighting->GetBinContent(
                  fSEMEReweighting->FindBin(binCent));
            }
            if (corrFactor != 0) {
              it->SetBinContent(iBin, binCont / corrFactor);
              it->SetBinError(iBin, binErr);
            } else {
              //dont worry if iBin seems to change, its due to the rebinning!
              std::cout << "======================================\n";
              std::cout << it->GetName() << std::endl;
              std::cout << "!for ikT = " << ikT << " and iBin = " << iBin
                        << "!\n";
              std::cout << "Correction factor 0, CFs meaningless!!! \n";
              std::cout << "======================================\n";
              std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            }
          }
          TString OutName = Form("%s_%sBin_%i", it->GetName(), variable, ikT);
          it->Write(OutName.Data());
        }
        fSum[ikT]->WriteOutput(outfileName.Data());
      }
    } else {
      std::cout << "No SE " << variable << " histogram with index 0 \n";
    }
  }
  return;
}

void DreamKayTee::AveragekT(const char *pair) {
  const char* variable = (fIskT) ? "kT" : "mT";
  fAveragekT = new TGraphErrors();
  TH2F* kTkStar = (TH2F*) fSEkT[0]->Clone(
      Form("%s%skStarForAverage", variable, pair));
  kTkStar->Add(fSEkT[1]);
  TH1F* kTProjection = (TH1F*) kTkStar->ProjectionY(
      Form("%s%sDist", variable, pair), kTkStar->GetXaxis()->FindBin(0.),
      kTkStar->GetXaxis()->FindBin(0.2), "e");
  auto *c1 = new TCanvas(Form("c%s", pair), Form("c%s", pair));
  c1->Divide(2, 1);
  c1->cd(1);
  kTkStar->Draw("COLZ");
  c1->cd(2);
  kTProjection->Draw();
  c1->SaveAs(Form("kTProjection%s.pdf", pair));
  for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
    int binLow = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT));
    int binUp = kTProjection->GetXaxis()->FindBin(fKayTeeBins.at(ikT + 1));
    kTProjection->GetXaxis()->SetRange(binLow + 1, binUp);
    float binErr = kTProjection->GetRMS();
    if (binLow + 1 == binUp) {
      std::cout << "Low Edge: " << kTProjection->GetBinLowEdge(binLow)
                << " Up Edge: " << kTProjection->GetBinLowEdge(binLow + 1)
                << std::endl;
      binErr = (kTProjection->GetBinLowEdge(binLow + 1)
          - kTProjection->GetBinLowEdge(binLow)) / 2.;
    }
    std::cout << ikT << '\t' << "Bin Low: " << binLow + 1 << "(="
              << fKayTeeBins.at(ikT) << ")" << '\t' << "Bin Up: " << binUp
              << "(=" << fKayTeeBins.at(ikT + 1) << ")" << '\t'
              << kTProjection->GetMean() << '\t' << binErr << std::endl;
    fAveragekT->SetPoint(ikT, ikT + 1, kTProjection->GetMean());
    fAveragekT->SetPointError(ikT, 0, binErr);
  }
}

void DreamKayTee::SetSEMEReweightingRatio(const char* pathToFile,
                                          TString pair) {
  TFile* CalibFile = TFile::Open(pathToFile, "read");
  TString HistName = "SEMEReweightingRatio_";
  if (pair == TString("pp")) {
    TString HistNamePair1 = HistName;
    HistNamePair1 += "Particle0_Particle0";
    TString HistNamePair2 = HistName;
    HistNamePair2 += "Particle1_Particle1";
    TH1F* RatioPair1 = (TH1F*) CalibFile->Get(HistNamePair1.Data());
    if (!RatioPair1) {
      std::cout << "Ratio Pair 1 missing \n";
    }
    TH1F* RatioPair2 = (TH1F*) CalibFile->Get(HistNamePair2.Data());
    if (!RatioPair2) {
      std::cout << "Ratio Pair 2 missing \n";
    }
    TF1* tfconstant = new TF1("myConst", "pol0", 0, 3000);
    tfconstant->SetParameter(0, 1.);
    RatioPair1->Divide(tfconstant, 2);
    RatioPair2->Divide(tfconstant, 2);
    fSEMEReweighting = (TH1F*) RatioPair1->Clone("ppReweightingFactor");
    fSEMEReweighting->Add(RatioPair2);
    TString MeVName = Form("%sMeV", fSEMEReweighting->GetName());
    fSEMEReweightingMeV = DreamCF::ConvertToOtherUnit(fSEMEReweighting, 1000,
                                                      MeVName.Data());
  } else if (pair == TString("pL")) {
    TString HistNamePair1 = HistName;
    HistNamePair1 += "Particle0_Particle2_Rebinned_5";
    TString HistNamePair2 = HistName;
    HistNamePair2 += "Particle1_Particle3_Rebinned_5";
    TH1F* RatioPair1 = (TH1F*) CalibFile->Get(HistNamePair1.Data());
    if (!RatioPair1) {
      std::cout << "Ratio Pair 1 missing \n";
    }
    TH1F* RatioPair2 = (TH1F*) CalibFile->Get(HistNamePair2.Data());
    if (!RatioPair2) {
      std::cout << "Ratio Pair 2 missing \n";
    }
    TF1* tfconstant = new TF1("myConst", "pol0", 0, 3000);
    tfconstant->SetParameter(0, 1.);
    RatioPair1->Divide(tfconstant, 2);
    RatioPair2->Divide(tfconstant, 2);
    fSEMEReweighting = (TH1F*) RatioPair1->Clone("pLReweightingFactor");
    fSEMEReweighting->Add(RatioPair2);
    TString MeVName = Form("%sMeV", fSEMEReweighting->GetName());
    fSEMEReweightingMeV = DreamCF::ConvertToOtherUnit(fSEMEReweighting, 1000,
                                                      MeVName.Data());
  } else if (pair == TString("pXi")) {
    TString HistNamePair1 = HistName;
    HistNamePair1 += "Particle0_Particle4_Rebinned_5";
    TString HistNamePair2 = HistName;
    HistNamePair2 += "Particle1_Particle5_Rebinned_5";
    TH1F* RatioPair1 = (TH1F*) CalibFile->Get(HistNamePair1.Data());
    if (!RatioPair1) {
      std::cout << "Ratio Pair 1 missing \n";
    }
    TH1F* RatioPair2 = (TH1F*) CalibFile->Get(HistNamePair2.Data());
    if (!RatioPair2) {
      std::cout << "Ratio Pair 2 missing \n";
    }
    TF1* tfconstant = new TF1("myConst", "pol0", 0, 3000);
    tfconstant->SetParameter(0, 1.);
    RatioPair1->Divide(tfconstant, 2);
    RatioPair2->Divide(tfconstant, 2);
    fSEMEReweighting = (TH1F*) RatioPair1->Clone("pXiReweightingFactor");
    fSEMEReweighting->Add(RatioPair2);
    TString MeVName = Form("%sMeV", fSEMEReweighting->GetName());
    fSEMEReweightingMeV = DreamCF::ConvertToOtherUnit(fSEMEReweighting, 1000,
                                                      MeVName.Data());
  } else {
    std::cout << "======================================\n";
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    std::cout << pair.Data() << "Pair not implemented \n";
    std::cout << "======================================\n";
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }
  return;
}

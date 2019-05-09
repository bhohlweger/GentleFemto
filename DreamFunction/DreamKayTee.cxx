/*
 * DreamKayTee.cxx
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#include "DreamKayTee.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include <iostream>
DreamKayTee::DreamKayTee()
    : fIskT(true),
      fFixShift(),
      fFixShiftValue(),
      fKayTeeBins(),
      fNKayTeeBins(0),
      fRebin(),
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
        for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
          TString PairName = Form("%s_%s_%i", PartName[iPart].Data(), variable,
                                  ikT);
          fCFPart[iPart][ikT] = new DreamPair(PairName.Data(), fNormleft,
                                              fNormright);
          //multiplications due in order not to hit the border of the bin!
          int kTminBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT] * 1.0001);
          int kTmaxBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT + 1] * 0.9999);

          DreamDist* kTDist = new DreamDist();
          TString SEkTBinName = Form("SE%s", PairName.Data());
          kTDist->SetSEDist(
              (TH1F*) fSEkT[iPart]->ProjectionX(SEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          TString MEkTBinName = Form("ME%s", PairName.Data());
          kTDist->SetMEDist(
              (TH1F*) fMEkT[iPart]->ProjectionX(MEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          //little hack to make things work.
          TH2F* fake = new TH2F("fakeMult", "fakeMult", 10, 0, 1, 10, 0, 1);
          kTDist->SetSEMultDist(fake, "SE");
          kTDist->SetMEMultDist(fake, "ME");
          fCFPart[iPart][ikT]->SetPair(kTDist);
          if (fFixShift.size() > 0 && fFixShift.at(ikT)) {
            fCFPart[iPart][ikT]->ShiftForEmpty(fCFPart[iPart][ikT]->GetPair());
            fCFPart[iPart][ikT]->FixShift(fCFPart[iPart][ikT]->GetPair(),
                                          nullptr, fFixShiftValue.at(ikT),
                                          true);
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(
                  fCFPart[iPart][ikT]->GetFixShifted().at(0), it);
            }
          } else {
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(), it);
            }
          }
        }
      }
      this->AveragekT(pair);
      fSum = new DreamCF*[fNKayTeeBins];
      TString outname = TString::Format("%s/AveragemT.root",
                                        outFolder);
      TFile* allCFsOut = TFile::Open(outname.Data(), "UPDATE");
      allCFsOut->cd();
      if (fAveragekT) {
        TString name = TString::Format("Average%s_%s", variable,pair);
        TGraphErrors* copy = new TGraphErrors(*fAveragekT);
        copy->SetName(name.Data());
        copy->Write(name.Data());
        allCFsOut->Write();
      }
      allCFsOut->Close();
      for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
        fSum[ikT] = new DreamCF();
        fSum[ikT]->SetPairs(fCFPart[0][ikT], fCFPart[1][ikT]);
        fSum[ikT]->GetCorrelations(pair);
        std::vector<TH1F*> CFs = fSum[ikT]->GetCorrelationFunctions();
        TString outfileName = TString::Format("%s/CFOutput_%s_%s_%s_%u.root",
                                              outFolder, variable, pair, prefix,
                                              ikT);
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
              it->SetBinContent(iBin, binCont * corrFactor);
              it->SetBinError(iBin, binErr);
            }
//            else {
//              //dont worry if iBin seems to change, its due to the rebinning!
//              std::cout << "======================================\n";
//              std::cout << it->GetName() << std::endl;
//              std::cout << "!for ikT = " << ikT << " and iBin = " << iBin
//                        << "!\n";
//              std::cout << "Correction factor 0, CFs meaningless!!! \n";
//              std::cout << "======================================\n";
//              std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
//            }
          }
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
  if (fAveragekT) {
    delete fAveragekT;
  }
  fAveragekT = new TGraphErrors();
  TH2F* kTkStar = (TH2F*) fSEkT[0]->Clone(
      Form("%s%skStarForAverage", variable, pair));
  kTkStar->Add(fSEkT[1]);
  TH1F* kTProjection = (TH1F*) kTkStar->ProjectionY(
      Form("%s%sDist", variable, pair), kTkStar->GetXaxis()->FindBin(0.),
      kTkStar->GetXaxis()->FindBin(0.2), "e");
  auto *c1 = new TCanvas(Form("c1%s", pair), Form("c1%s", pair));
  auto *c2 = new TCanvas(Form("c2%s", pair), Form("c2%s", pair));
  c1->cd();
  kTkStar->Draw("COLZ");
  c1->SaveAs(Form("mTvskStar%s.pdf", pair));
  c2->cd();
  float totalPairs = kTProjection->Integral();
  kTProjection->Draw();
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(kTProjection);
  kTProjection->GetXaxis()->SetTitle("<m_{T} > (GeV/c^{2})");
  c2->SaveAs(Form("mTDistribution%s.pdf", pair));
  for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
    int binLow = kTProjection->GetXaxis()->FindBin(
        fKayTeeBins.at(ikT) * 1.0001);
    int binUp = kTProjection->GetXaxis()->FindBin(
        fKayTeeBins.at(ikT + 1) * 0.9999);
    kTProjection->GetXaxis()->SetRange(binLow, binUp);
    float binErr = kTProjection->GetMeanError();
    if (binLow == binUp) {
      std::cout << "Low Edge: " << kTProjection->GetBinLowEdge(binLow)
                << " Up Edge: " << kTProjection->GetBinLowEdge(binLow + 1)
                << std::endl;
      binErr = (kTProjection->GetBinLowEdge(binLow + 1)
          - kTProjection->GetBinLowEdge(binLow)) / 2.;
      std::cout << "Percent of total: "
                << kTProjection->GetBinContent(binLow) * 100. / totalPairs
                << std::endl;
    } else {
      std::cout << "Percent of total: "
                << kTProjection->Integral(binLow, binUp) * 100. / totalPairs
                << std::endl;
    }
    std::cout << ikT << '\t' << "Bin Low: " << binLow << "(="
              << fKayTeeBins.at(ikT) << ")" << '\t' << "Bin Up: " << binUp
              << "(=" << fKayTeeBins.at(ikT + 1) << ")" << '\t'
              << kTProjection->GetMean() << '\t' << binErr << std::endl;
    fAveragekT->SetPoint(ikT, ikT + 1, kTProjection->GetMean());
    fAveragekT->SetPointError(ikT, 0, binErr);
  }
}

void DreamKayTee::SetSEMEReweightingRatio(const char* pathToFile,
                                          const char* HistNum,
                                          bool useRebinned) {
  TFile* CalibFile = TFile::Open(pathToFile, "read");
  if (!CalibFile) {
    std::cout << "No Calibration file found in " << pathToFile << std::endl;
  }
  TString HistName = useRebinned ? "Rebinned" : "FixShifted";
  SetSEMEReweightingRatio(
      (TH1F*) CalibFile->Get(
          TString::Format("hCk_%s_%s", HistName.Data(), HistNum)),
      (TH1F*) CalibFile->Get(TString::Format("hCk_Reweighted_%s", HistNum)),
      (TH1F*) CalibFile->Get(
          TString::Format("hCk_%sMeV_%s", HistName.Data(), HistNum)),
      (TH1F*) CalibFile->Get(TString::Format("hCk_ReweightedMeV_%s", HistNum)));
  return;
}

void DreamKayTee::SetSEMEReweightingRatio(TH1F* Rebinned, TH1F* Reweighted,
                                          TH1F* RebinnedMeV,
                                          TH1F* ReweightedMeV) {
  fSEMEReweighting = (TH1F*) Reweighted->Clone(
      TString::Format("hCk_Reweighting"));
  fSEMEReweighting->Divide(Rebinned);
  fSEMEReweightingMeV = (TH1F*) ReweightedMeV->Clone(
      TString::Format("hCk_ReweightingMeV"));
  fSEMEReweightingMeV->Divide(RebinnedMeV);
  return;
}

/*
 * DreamCF.cxx
 *
 *  Created on: Aug 22, 2018
 *      Author: hohlweger
 */

#include "DreamCF.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include "TObject.h"

DreamCF::DreamCF()
    : fCF(),
      fGrCF(),
      fRatio(),
      fPairOne(nullptr),
      fPairTwo(nullptr) {

}

DreamCF::~DreamCF() {
  for (auto it : fCF) {
    delete it;
  }
  for (auto it : fGrCF) {
    delete it;
  }
  for (auto it : fRatio) {
    delete it;
  }
  delete fPairOne;
  delete fPairTwo;
}

void DreamCF::GetCorrelations(const char* pairName) {
  if (fPairOne && fPairOne->GetPair()) {
    if (fPairTwo && fPairTwo->GetPair()) {
      TH1F* CFSum = AddCF(fPairOne->GetPair()->GetCF(),
                          fPairTwo->GetPair()->GetCF(),
                          Form("hCkTotNormWeight%s", pairName));
      if (CFSum) {
        fCF.push_back(CFSum);
        TString CFSumMeVName = Form("%sMeV", CFSum->GetName());
        TH1F* CFMeVSum = ConvertToOtherUnit(CFSum, 1000, CFSumMeVName.Data());
        if (CFMeVSum) {
          fCF.push_back(CFMeVSum);
        }
      } else {
        Warning("DreamCF", "No Pair 2 Set, only setting Pair 1!");
        //existence already checked
        fCF.push_back(fPairOne->GetPair()->GetCF());
        fGrCF.push_back(fPairOne->GetPair()->GetGrCF());
      }
    }
  } else {
    Warning("DreamCF", "No Pair 1 Set, only setting Pair 2");
    if (fPairTwo && fPairTwo->GetPair()) {
      fCF.push_back(fPairTwo->GetPair()->GetCF());
      fGrCF.push_back(fPairTwo->GetPair()->GetGrCF());
    } else {
      Warning("DreamCF", "No Pair 2 set either");
    }
  }
  if (fPairOne && fPairTwo) {
    if (fPairOne->GetNDists() == fPairTwo->GetNDists()) {
      LoopCorrelations(fPairOne->GetShiftedEmpty(), fPairTwo->GetShiftedEmpty(),
                       Form("hCk_Shifted%s", pairName));
      LoopCorrelations(fPairOne->GetFixShifted(), fPairTwo->GetFixShifted(),
                       Form("hCk_FixShifted%s", pairName));
      LoopCorrelations(fPairOne->GetRebinned(), fPairTwo->GetRebinned(),
                       Form("hCk_Rebinned%s", pairName));
      LoopCorrelations(fPairOne->GetReweighted(), fPairTwo->GetReweighted(),
                       Form("hCk_Reweighted%s", pairName));
      LoopCorrelations(fPairOne->GetUnfolded(), fPairTwo->GetUnfolded(),
                       Form("hCk_Unfolded%s", pairName));
    } else {
      Warning("DreamCF", "Pair 1 with %i histograms, Pair 2 with %i histograms",
              fPairOne->GetNDists(), fPairTwo->GetNDists());
    }
  } else if (fPairOne && !fPairTwo) {
    LoopCorrelations(fPairOne->GetShiftedEmpty(),
                     Form("hCk_Shifted%s", pairName));
    LoopCorrelations(fPairOne->GetFixShifted(),
                     Form("hCk_FixShifted%s", pairName));
    LoopCorrelations(fPairOne->GetRebinned(), Form("hCk_Rebinned%s", pairName));
    LoopCorrelations(fPairOne->GetReweighted(),
                     Form("hCk_Reweighted%s", pairName));
    LoopCorrelations(fPairOne->GetUnfolded(), Form("hCk_Unfolded%s", pairName));
  } else if (fPairTwo && !fPairOne) {
    LoopCorrelations(fPairTwo->GetShiftedEmpty(),
                     Form("hCk_Shifted%s", pairName));
    LoopCorrelations(fPairTwo->GetFixShifted(),
                     Form("hCk_FixShifted%s", pairName));
    LoopCorrelations(fPairTwo->GetRebinned(), Form("hCk_Rebinned%s", pairName));
    LoopCorrelations(fPairTwo->GetReweighted(),
                     Form("hCk_Reweighted%s", pairName));
    LoopCorrelations(fPairTwo->GetUnfolded(), Form("hCk_Unfolded%s", pairName));
  } else {
    Error("DreamCF", "Pair 1 and Pair 2 missing");
  }
  return;
}

void DreamCF::LoopCorrelations(std::vector<DreamDist*> PairOne,
                               std::vector<DreamDist*> PairTwo,
                               const char* name) {
  if (PairOne.size() != PairTwo.size()) {
    Warning("DreamCF", "Different size of pair (%i) and antiparticle pair (%i)",
            PairOne.size(), PairTwo.size());
  } else {
    unsigned int iIter = 0;
    while (iIter < PairOne.size()) {
      DreamDist* PartPair = PairOne.at(iIter);
      DreamDist* AntiPartPair = PairTwo.at(iIter);
      TString CFSumName = Form("%s_%i", name, iIter);
      TH1F* CFSum = AddCF(PartPair->GetCF(), AntiPartPair->GetCF(),
                          CFSumName.Data());
      if (CFSum) {
        fCF.push_back(CFSum);
        auto CFgrSum = AddCF(CFSum, { { fPairOne, fPairTwo } },
                             Form("Gr%s", CFSumName.Data()));
        TString CFSumMeVName = Form("%sMeV_%i", name, iIter);
        TH1F* CFMeVSum = ConvertToOtherUnit(CFSum, 1000, CFSumMeVName.Data());
        if (CFMeVSum) {
          fCF.push_back(CFMeVSum);
        }
        if (CFgrSum) {
          fGrCF.push_back(CFgrSum);
          TString CFSumMeVName = Form("%sMeV", CFgrSum->GetName());
          auto CFMeVSum = ConvertToOtherUnit(CFgrSum, 1000,
                                             CFSumMeVName.Data());
          if (CFMeVSum) {
            fGrCF.push_back(CFMeVSum);
          }
        }
      } else {
        if (PartPair->GetCF()) {
          Warning("DreamCF",
                  "For iteration %i Particle Pair CF (%s) is missing", iIter,
                  PartPair->GetSEDist()->GetName());
        } else if (AntiPartPair->GetCF()) {
          Warning("DreamCF",
                  "For iteration %i AntiParticle Pair CF (%s) is missing",
                  iIter, AntiPartPair->GetSEDist()->GetName());
        }
      }
      iIter++;
    }
  }
}

void DreamCF::LoopCorrelations(std::vector<DreamDist*> Pair, const char* name) {
  unsigned int iIter = 0;
  for (auto it : Pair) {
    TString CFSumName = Form("%s_%i", name, iIter++);
    TH1F* CFClone = (TH1F*) it->GetCF()->Clone(CFSumName.Data());
    if (CFClone) {
      fCF.push_back(CFClone);
      TString CFSumMeVName = Form("%sMeV_%i", name, iIter);
      TH1F* CFMeVSum = ConvertToOtherUnit(CFClone, 1000, CFSumMeVName.Data());
      if (CFMeVSum) {
        fCF.push_back(CFMeVSum);
      }
    } else {
      Warning("DreamCF", "For iteration %i Particle Pair CF (%s) is missing \n",
              iIter, it->GetSEDist()->GetName());
    }
  }
}

void DreamCF::WriteOutput(const char* name) {
  TFile* output = TFile::Open(name, "RECREATE");
  WriteOutput(output, true);
}

void DreamCF::WriteOutput(TFile* output, bool closeFile) {
  output->cd();
  for (auto& it : fCF) {
    it->Write();
    delete it;
  }
  for (auto& it : fGrCF) {
    it->Write();
    delete it;
  }
  if (fRatio.size() > 0) {
    TList *RatioList = new TList();
    RatioList->SetOwner();
    RatioList->SetName("Ratios");
    for (auto& it : fRatio) {
      RatioList->Add(it);
    }
    RatioList->Write("RatioList", 1);
  }
  if (fPairOne) {
    TList *PairDist = new TList();
    PairDist->SetOwner();
    PairDist->SetName("PairDist");
    fPairOne->WriteOutput(PairDist);
    PairDist->Write("PairDist", 1);
  } else {
    Warning("DreamCF", "not writing Pair 1");
  }
  if (fPairTwo) {
    TList *AntiPairDist = new TList();
    AntiPairDist->SetOwner();
    AntiPairDist->SetName("AntiPairDist");
    fPairTwo->WriteOutput(AntiPairDist);
    AntiPairDist->Write("AntiPairDist", 1);
  } else {
    Warning("DreamCF", "not writing Pair 2");
  }
  if (closeFile)
    output->Close();
  return;
}

TH1F* DreamCF::AddCF(TH1F* CF1, TH1F* CF2, const char* name) {
  TH1F* hist_CF_sum = nullptr;
  if (CF1 && CF2) {
    if (CF1->GetXaxis()->GetXmin() == CF2->GetXaxis()->GetXmin()) {
      TH1F* Ratio = (TH1F*) CF1->Clone(TString::Format("%sRatio", name));
      Ratio->Divide(CF2);
      fRatio.push_back(Ratio);
      //Calculate CFs with error weighting
      hist_CF_sum = (TH1F*) CF1->Clone(name);

      int NBins = hist_CF_sum->GetNbinsX();

      for (int i = 0; i < NBins; i++) {
        double CF1_val = CF1->GetBinContent(i + 1);
        double CF1_err = CF1->GetBinError(i + 1);
        double CF2_val = CF2->GetBinContent(i + 1);
        double CF2_err = CF2->GetBinError(i + 1);
        //average for bin i:
        if (CF1_val != 0. && CF2_val != 0.) {
          double CF1_err_weight = 1. / TMath::Power(CF1_err, 2.);
          double CF2_err_weight = 1. / TMath::Power(CF2_err, 2.);

          double CF_sum_average = (CF1_err_weight * CF1_val
              + CF2_err_weight * CF2_val) / (CF1_err_weight + CF2_err_weight);
          double CF_sum_err = 1. / TMath::Sqrt(CF1_err_weight + CF2_err_weight);

          hist_CF_sum->SetBinContent(i + 1, CF_sum_average);
          hist_CF_sum->SetBinError(i + 1, CF_sum_err);
        } else if (CF1_val == 0. && CF2_val != 0.) {
          hist_CF_sum->SetBinContent(i + 1, CF2_val);
          hist_CF_sum->SetBinError(i + 1, CF2_err);
        } else if (CF2_val == 0 && CF1_val != 0.) {
          hist_CF_sum->SetBinContent(i + 1, CF1_val);
          hist_CF_sum->SetBinError(i + 1, CF1_err);
        }
      }
    } else {
      Warning(
          "DreamCF",
          "Skipping %s and %s due to uneven beginning of binning (%.3f) and (%.3f)",
          CF1->GetName(), CF2->GetName(), CF1->GetXaxis()->GetXmin(),
          CF2->GetXaxis()->GetXmin());
    }
  }
  return hist_CF_sum;
}

TGraphAsymmErrors* DreamCF::AddCF(TH1F* histSum, std::vector<DreamPair*> pairs,
                                  const char* name) {
  TGraphAsymmErrors* hist_CF_sum = new TGraphAsymmErrors(histSum);
  hist_CF_sum->Set(0);

  TList *list = new TList;
  for (auto it : pairs) {
    list->Add(it->GetPair()->GetSEDist());
  }
  TH1F *histSE = (TH1F*) pairs.at(0)->GetPair()->GetSEDist()->Clone(
      Form("%s_cloneForMixed", pairs.at(0)->GetPair()->GetSEDist()->GetName()));
  histSE->Reset();
  histSE->Merge(list);

  // Shift the center of the bin in x according to the same event distribution
  // In case of very large bins, this can have a sizeable effect!
  int counter = 0;
  float xVal, xErrRight, xErrLeft, centrVal;
  for (int i = 1; i <= histSum->GetNbinsX(); ++i) {
    if (histSum->GetBinContent(i) == 0)
      continue;
    // In case we have a SE distribution available we use the mean in that k* bin
    centrVal = histSum->GetBinCenter(i);
    histSE->GetXaxis()->SetRangeUser(
        histSum->GetBinLowEdge(i),
        histSum->GetBinLowEdge(i) + histSum->GetBinWidth(i));
    xVal = histSE->GetMean();
    xErrLeft = xVal - centrVal + histSum->GetBinWidth(i) / 2.;
    xErrRight = centrVal - xVal + histSum->GetBinWidth(i) / 2.;
    hist_CF_sum->SetPoint(counter, xVal, histSum->GetBinContent(i));
    hist_CF_sum->SetPointError(counter++, xErrLeft, xErrRight,
                               histSum->GetBinError(i),
                               histSum->GetBinError(i));
  }
  delete list;
  return hist_CF_sum;
}

TGraphAsymmErrors* DreamCF::ConvertToOtherUnit(TGraphAsymmErrors* HistCF,
                                               int Scale, const char* name) {
  int nBins = HistCF->GetN();
  TGraphAsymmErrors* hist_CF_MeV = new TGraphAsymmErrors();
  hist_CF_MeV->SetName(name);
  double x, y;
  for (int i = 0; i < nBins; ++i) {
    HistCF->GetPoint(i, x, y);
    hist_CF_MeV->SetPoint(i, x * 1000.f, y);
    hist_CF_MeV->SetPointError(i, HistCF->GetErrorXlow(i) * 1000.f,
                               HistCF->GetErrorXhigh(i) * 1000.f,
                               HistCF->GetErrorYlow(i),
                               HistCF->GetErrorYhigh(i));
  }
  return hist_CF_MeV;
}

TH1F* DreamCF::ConvertToOtherUnit(TH1F* HistCF, int Scale, const char* name) {
  int nBins = HistCF->GetNbinsX();
  float kMin = HistCF->GetXaxis()->GetXmin();
  float kMax = HistCF->GetXaxis()->GetXmax();
  TH1F* HistScaled = new TH1F(name, name, nBins, kMin * Scale, kMax * Scale);
  for (int iBin = 1; iBin <= nBins; ++iBin) {
    HistScaled->SetBinContent(iBin, HistCF->GetBinContent(iBin));
    HistScaled->SetBinError(iBin, HistCF->GetBinError(iBin));
  }
  return HistScaled;
}

TH1F* DreamCF::FindCorrelationFunction(TString name) {
  TH1F* output = nullptr;
  for (auto it : fCF) {
    TString itName = it->GetName();
    if (itName.Contains(name.Data())) {
      Warning("DreamCF", "For Histo: %s \t we use the %s", name.Data(),
              itName.Data());
      output = it;
    }
  }
  if (!output) {
    Warning("DreamCF", "Output Histogram not found for %s", name.Data());
    Warning("DreamCF", "What we offer is the following:");
    for (auto it : fCF) {
      Warning("DreamCF", "%s", it->GetName());
    }
  }

  return output;
}

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

DreamCF::DreamCF()
    : fCF(),
      fPartPair(nullptr),
      fAntiPartPair(nullptr) {

}

DreamCF::~DreamCF() {
}

void DreamCF::GetCorrelations() {
  if (fPartPair->GetPair()) {
    if (fAntiPartPair->GetPair()) {
      TH1F* CFSum = AddCF(fPartPair->GetPair()->GetCF(),
                          fAntiPartPair->GetPair()->GetCF(),
                          "hCkTotNormWeight");
      if (CFSum) {
        fCF.push_back(CFSum);
      }
    } else {
      std::cout << "No Anti-Particle Pair Set! \n";
    }
  } else {
    std::cout << "No Particle Pair Set! \n";
  }
  if (fPartPair->GetNDists() == fAntiPartPair->GetNDists()) {
    LoopCorrelations(fPartPair->GetShiftedEmpty(),
                     fAntiPartPair->GetShiftedEmpty(), "hCk_Shifted");
    LoopCorrelations(fPartPair->GetFixShifted(), fAntiPartPair->GetFixShifted(),
                     "hCk_FixShifted");
    LoopCorrelations(fPartPair->GetRebinned(), fAntiPartPair->GetRebinned(),
                     "hCk_Rebinned");
    LoopCorrelations(fPartPair->GetReweighted(), fAntiPartPair->GetReweighted(),
                     "hCk_Reweighted");
  } else {
    std::cout << "Part Pair with " << fPartPair->GetNDists()
              << "Distributions, Anti Part Pair with "
              << fAntiPartPair->GetNDists() << std::endl;
  }
  return;
}

void DreamCF::LoopCorrelations(std::vector<DreamDist*> partPair,
                               std::vector<DreamDist*> antipartPair,
                               const char* name) {
  if (partPair.size() != antipartPair.size()) {
    std::cout << "Different size of pair(" << partPair.size()
              << ") and antiparticle pair(" << antipartPair.size() << ") ! \n";
  } else {
    unsigned int iIter = 0;
    while (iIter < partPair.size()) {
      DreamDist* PartPair = partPair.at(iIter);
      DreamDist* AntiPartPair = antipartPair.at(iIter);
      TString CFSumName = Form("%s_%i", name, iIter);
      TH1F* CFSum = AddCF(PartPair->GetCF(), AntiPartPair->GetCF(),
                          CFSumName.Data());
      if (CFSum) {
        fCF.push_back(CFSum);
      } else {
        if (PartPair->GetCF()) {
          std::cout << "For iteration " << iIter << " Particle Pair CF ("
                    << PartPair->GetSEDist()->GetName() << ")is missing \n";
        } else if (AntiPartPair->GetCF()) {
          std::cout << "For iteration " << iIter << " AntiParticle Pair CF ("
                    << AntiPartPair->GetSEDist()->GetName() << ")is missing \n";
        }
      }
      iIter++;
    }
  }
}

void DreamCF::WriteOutput(const char* name) {
  TFile* output = TFile::Open(name, "RECREATE");
  for (auto& it : fCF) {
    it->Write();
  }
  TList *PairDist = new TList();
  PairDist->SetOwner();
  PairDist->SetName("PairDist");
  fPartPair->WriteOutput(PairDist);
  PairDist->Write("PairDist", 1);

  TList *AntiPairDist = new TList();
  AntiPairDist->SetOwner();
  AntiPairDist->SetName("AntiPairDist");
  fAntiPartPair->WriteOutput(AntiPairDist);
  AntiPairDist->Write("AntiPairDist", 1);

  output->Close();
  return;
}

TH1F* DreamCF::AddCF(TH1F* CF1, TH1F* CF2, const char* name) {
  TH1F* hist_CF_sum = nullptr;
  if (CF1 && CF2) {
    if (CF1->GetXaxis()->GetXmin() == CF2->GetXaxis()->GetXmin()) {
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
      std::cout << "Skipping " << CF1->GetName() << " and " << CF2->GetName()
                << " due to uneven beginning of binning \n";
    }
  }
  return hist_CF_sum;
}


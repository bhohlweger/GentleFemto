/*
 * SideBandFit.cxx
 *
 *  Created on: Nov 6, 2018
 *      Author: hohlweger
 */

#include <SideBandFit.h>
#include "ReadDreamFile.h"
#include "DreamPair.h"
#include "DreamCF.h"
#include "TMath.h"
#include <iostream>
#include "TF1.h"
SideBandFit::SideBandFit()
    : fAnalysisFileUp(),
      fAnalysisFileDown(),
      fnormleft(200),
      fnormright(400),
      fRebin(5) {
  // TODO Auto-generated constructor stub

}

SideBandFit::~SideBandFit() {
  // TODO Auto-generated destructor stub
}

void SideBandFit::SetSideBandFile(const char* path, const char* suffixUp,
                                  const char* suffixDown) {

  TString filename = Form("%s/AnalysisResults.root", path);
  fAnalysisFileUp = new ReadDreamFile(6, 6);
  fAnalysisFileUp->SetAnalysisFile(filename.Data(), "MB", suffixUp);

  fAnalysisFileDown = new ReadDreamFile(6, 6);
  fAnalysisFileDown->SetAnalysisFile(filename.Data(), "MB", suffixDown);
}

void SideBandFit::SideBandCFs(bool doQA) {
  for (auto it : fSideBands) {
    delete it;
  }
  fSideBands.clear();

  for (auto it : fSideBandCFs) {
    delete it;
  }
  fSideBandCFs.clear();

  std::vector<const char*> sideNames = { "SideUp", "AntiSideUp", "SideDown",
      "AntiSideDown" };
  for (int iSB = 0; iSB < 4; ++iSB) {
    fSideBands.emplace_back(
        new DreamPair(sideNames[iSB], fnormleft / 1000., fnormright / 1000.));
  }
  fSideBands[0]->SetPair(fAnalysisFileUp->GetPairDistributions(0, 4, "Up"));
  fSideBands[1]->SetPair(fAnalysisFileUp->GetPairDistributions(1, 5, "AntiUp"));
  fSideBands[2]->SetPair(fAnalysisFileDown->GetPairDistributions(0, 4, "Down"));
  fSideBands[3]->SetPair(
      fAnalysisFileDown->GetPairDistributions(1, 5, "AntiDown"));
  for (auto it : fSideBands) {
    it->ShiftForEmpty(it->GetPair());
  }
  int iPos = 0;
  int itCounter = 0;
  float minVal = 3001;
  for (auto it : fSideBands) {
    if (it->GetFirstBin() < minVal) {
      minVal = it->GetFirstBin();
      iPos = itCounter;
    }
    itCounter++;
  }
  for (auto it : fSideBands) {
    it->FixShift(it->GetPair(), fSideBands[iPos]->GetPair(), minVal, true);
    it->Rebin(it->GetPair(), fRebin);
    it->Rebin(it->GetPairFixShifted(0), fRebin);
    it->ReweightMixedEvent(it->GetPairRebinned(0), 0.2, 0.9);
  }

  DreamCF* unitConv = new DreamCF();
  TH1F* CF1Reweighted = fSideBands[0]->GetReweighted().at(0)->GetCF();
  TH1F* CF2Reweighted = fSideBands[1]->GetReweighted().at(0)->GetCF();
  TH1F* CF3Reweighted = fSideBands[2]->GetReweighted().at(0)->GetCF();
  TH1F* CF4Reweighted = fSideBands[3]->GetReweighted().at(0)->GetCF();

  TH1F* CF1Rebinned_NoShift = fSideBands[0]->GetRebinned().at(0)->GetCF();
  TH1F* CF2Rebinned_NoShift = fSideBands[1]->GetRebinned().at(0)->GetCF();
  TH1F* CF3Rebinned_NoShift = fSideBands[2]->GetRebinned().at(0)->GetCF();
  TH1F* CF4Rebinned_NoShift = fSideBands[3]->GetRebinned().at(0)->GetCF();
  TH1F* SideBandSumRebinNoShift = AddCF(CF1Rebinned_NoShift,
                                        CF2Rebinned_NoShift,
                                        CF3Rebinned_NoShift,
                                        CF4Rebinned_NoShift, "RebinnedNoShift");
  fSideBandCFs.push_back(SideBandSumRebinNoShift);
  fSideBandCFs.push_back(
      unitConv->ConvertToOtherUnit(SideBandSumRebinNoShift, 1000,
                                   "RebinnedNoShift_MeV"));
  TH1F* CF1Rebinned_Shift = fSideBands[0]->GetRebinned().at(1)->GetCF();
  TH1F* CF2Rebinned_Shift = fSideBands[1]->GetRebinned().at(1)->GetCF();
  TH1F* CF3Rebinned_Shift = fSideBands[2]->GetRebinned().at(1)->GetCF();
  TH1F* CF4Rebinned_Shift = fSideBands[3]->GetRebinned().at(1)->GetCF();
  TH1F* SideBandSumRebinShift = AddCF(CF1Rebinned_Shift, CF2Rebinned_Shift,
                                      CF3Rebinned_Shift, CF4Rebinned_Shift,
                                      "RebinnedShift");
  fSideBandCFs.push_back(SideBandSumRebinShift);
  fSideBandCFs.push_back(
      unitConv->ConvertToOtherUnit(SideBandSumRebinShift, 1000,
                                   "RebinnedShift_MeV"));

  TH1F* SideBandSumReweighted = AddCF(CF1Reweighted, CF2Reweighted,
                                      CF3Reweighted, CF4Reweighted,
                                      "Reweighted");
  fSideBandCFs.push_back(SideBandSumReweighted);
  fSideBandCFs.push_back(
      unitConv->ConvertToOtherUnit(SideBandSumReweighted, 1000,
                                   "Reweighted_MeV"));
  if (doQA) {
    TH1F* SidebandUpSum = unitConv->AddCF(
        fSideBands[0]->GetRebinned().at(1)->GetCF(),
        fSideBands[1]->GetRebinned().at(1)->GetCF(), "UpSum");
    fSideBandCFs.push_back(
        unitConv->ConvertToOtherUnit(SidebandUpSum, 1000, "SideBandSumUpMeV"));
    TH1F* SidebandDownSum = unitConv->AddCF(
        fSideBands[2]->GetRebinned().at(1)->GetCF(),
        fSideBands[3]->GetRebinned().at(1)->GetCF(), "DownSum");
    fSideBandCFs.push_back(
        unitConv->ConvertToOtherUnit(SidebandDownSum, 1000,
                                     "SideBandSumDownMeV"));
    TH1F* SideBandRatio;
    for (auto it : fSideBandCFs) {
      if (it->GetName() == TString("SideBandSumUpMeV")) {
        SideBandRatio = (TH1F*)it->Clone("RatioSideBandMeV");
      }
      if (it->GetName() == TString("SideBandSumDownMeV")) {
        SideBandRatio->Divide(it);
      }
    }
    fSideBandCFs.push_back(SideBandRatio);
  }
  delete unitConv;
}

TH1F* SideBandFit::AddCF(TH1F* CF1, TH1F* CF2, TH1F* CF3, TH1F* CF4,
                         const char* name) {
  TH1F* hist_CF_sum = nullptr;
  if (CF1 && CF2 && CF3 && CF4) {
    double xMinCF1 = CF1->GetXaxis()->GetXmin();
    double xMinCF2 = CF2->GetXaxis()->GetXmin();
    double xMinCF3 = CF3->GetXaxis()->GetXmin();
    double xMinCF4 = CF4->GetXaxis()->GetXmin();

    if (xMinCF1 == xMinCF2 && xMinCF1 == xMinCF3 && xMinCF1 == xMinCF4) {
      //Calculate CFs with error weighting
      hist_CF_sum = (TH1F*) CF1->Clone(name);

      int NBins = hist_CF_sum->GetNbinsX();

      for (int i = 0; i < NBins; i++) {
        double CF1_val = CF1->GetBinContent(i + 1);
        double CF1_err = CF1->GetBinError(i + 1);
        double CF2_val = CF2->GetBinContent(i + 1);
        double CF2_err = CF2->GetBinError(i + 1);
        double CF3_val = CF3->GetBinContent(i + 1);
        double CF3_err = CF3->GetBinError(i + 1);
        double CF4_val = CF4->GetBinContent(i + 1);
        double CF4_err = CF4->GetBinError(i + 1);
        //average for bin i:
        double CF1_err_weight =
            CF1_val != 0. ? 1. / TMath::Power(CF1_err, 2.) : 0;
        double CF2_err_weight =
            CF2_val != 0. ? 1. / TMath::Power(CF2_err, 2.) : 0;
        double CF3_err_weight =
            CF3_val != 0. ? 1. / TMath::Power(CF3_err, 2.) : 0;
        double CF4_err_weight =
            CF4_val != 0. ? 1. / TMath::Power(CF4_err, 2.) : 0;

        if (CF1_val != 0. || CF2_val != 0. || CF3_val != 0. || CF4_val != 0.) {
          double CF_sum_average = (CF1_err_weight * CF1_val
              + CF2_err_weight * CF2_val + CF3_err_weight * CF3_val
              + CF4_err_weight * CF4_val)
              / (CF1_err_weight + CF2_err_weight + CF3_err_weight
                  + CF4_err_weight);
          double CF_sum_err = 1.
              / TMath::Sqrt(
                  CF1_err_weight + CF2_err_weight + CF3_err_weight
                      + CF4_err_weight);

          hist_CF_sum->SetBinContent(i + 1, CF_sum_average);
          hist_CF_sum->SetBinError(i + 1, CF_sum_err);
        } else {
          hist_CF_sum->SetBinContent(i + 1, 0);
          hist_CF_sum->SetBinError(i + 1, 0);
        }
      }
    } else {
      std::cout << "Skipping " << CF1->GetName() << " and " << CF2->GetName()
                << " due to uneven beginning of binning ("
                << CF1->GetXaxis()->GetXmin() << " and "
                << CF2->GetXaxis()->GetXmin() << ") \n";
    }
  }
  return hist_CF_sum;
}

void SideBandFit::WriteOutput(const char* outputPath) {
  TString name = Form("%s/SideBandCFs%2.0f_%2.0f.root", outputPath, fnormleft,
                      fnormright);
  TFile* output = TFile::Open(name, "RECREATE");
  output->cd();

  std::vector<const char*> sideNames = { "SideUp", "AntiSideUp", "SideDown",
      "AntiSideDown" , "SumSideDown" , "SumSideUp" };
  unsigned int iOut = 0;
  for (auto it : fSideBands) {
    TList *outlist = new TList();
    outlist->SetOwner();
    outlist->SetName(sideNames.at(iOut));
    it->WriteOutput(outlist);
    outlist->Write(sideNames.at(iOut++), 1);
  }
  output->cd();
  for (auto it : fSideBandCFs) {
    it->Write();
  }
  output->Close();
  return;
}

double SideBandFit::Parameterization(const double& Momentum,
                                     const double* SourcePar,
                                     const double* PotPar) {
  return PotPar[0] + PotPar[1] * Momentum
      + TMath::Exp(PotPar[2] + PotPar[3] * Momentum);
}

double SideBandFit::ParameterizationROOT(double* Momentum, double* PotPar) {
  return PotPar[0] + PotPar[1] * Momentum[0]
      + TMath::Exp(PotPar[2] + PotPar[3] * Momentum[0]);
}

void SideBandFit::FitSideBands(TH1F* cfSide, double* potPar) {
  TF1* funct = new TF1("SideBandFit", ParameterizationROOT, 0, 1000, 4);
  cfSide->Fit("SideBandFit", "Q, S, N, R, M");
  funct->GetParameters(potPar);
  cfSide->GetListOfFunctions()->Add(funct);
  return;
}

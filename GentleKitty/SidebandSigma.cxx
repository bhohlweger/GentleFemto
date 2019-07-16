#include "SidebandSigma.h"
#include "DreamPair.h"
#include "DreamCF.h"
#include "TMath.h"
#include <iostream>
#include "TF1.h"

SidebandSigma::SidebandSigma()
    : fAnalysisFile(),
      fSidebandCFGr() {
}

SidebandSigma::~SidebandSigma() {
  delete fAnalysisFile;
  for (auto it : fSidebandCFGr) {
    delete it;
  }
  fSidebandCFGr.resize(0);
}

void SidebandSigma::SetSideBandFile(const char* path, const char* trigger,
                                    const char* suffix) {
  TString filename = Form("%s/AnalysisResults.root", path);
  fAnalysisFile = new ReadDreamFile(8, 8);
  fAnalysisFile->SetAnalysisFile(filename.Data(), trigger, suffix);
}

void SidebandSigma::SideBandCFs() {
  for (auto it : fSideBands) {
    delete it;
  }
  fSideBands.clear();

  for (auto it : fSideBandCFs) {
    delete it;
  }
  fSideBandCFs.clear();

  std::vector<const char*> sideNames = { "Part_SB_up", "AntiPart_SB_up",
      "Part_SB_low", "AntiPart_SB_low" };
  for (int iSB = 0; iSB < 4; ++iSB) {
    fSideBands.emplace_back(
        new DreamPair(sideNames[iSB], fnormleft / 1000., fnormright / 1000.));
  }

  fSideBands[0]->SetPair(fAnalysisFile->GetPairDistributions(0, 4, "Up"));
  fSideBands[1]->SetPair(fAnalysisFile->GetPairDistributions(1, 5, "AntiUp"));

  fSideBands[2]->SetPair(fAnalysisFile->GetPairDistributions(0, 6, "Down"));
  fSideBands[3]->SetPair(fAnalysisFile->GetPairDistributions(1, 7, "AntiDown"));

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
  int iPair = 0;
  for (auto it : fSideBands) {
    it->FixShift(it->GetPair(), fSideBands[iPos]->GetPair(), minVal, true);
    it->Rebin(it->GetPair(), fRebin, true);
    it->Rebin(it->GetPairFixShifted(0), fRebin, true);
    it->ReweightMixedEvent(it->GetPairRebinned(0), 0.2, 0.9, fSideBands[iPair++]->GetPair());
  }

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
      DreamCF::ConvertToOtherUnit(SideBandSumRebinNoShift, 1000,
                                   "RebinnedNoShift_MeV"));
  TGraphAsymmErrors * grSideBandSumRebinNoShift = DreamCF::AddCF(
      SideBandSumRebinNoShift,
      { { fSideBands[0], fSideBands[1], fSideBands[2], fSideBands[3]} },
      Form("Gr%s", SideBandSumRebinNoShift->GetName()));
  fSidebandCFGr.push_back(grSideBandSumRebinNoShift);
  fSidebandCFGr.push_back(DreamCF::ConvertToOtherUnit(grSideBandSumRebinNoShift, 1000,
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
      DreamCF::ConvertToOtherUnit(SideBandSumRebinShift, 1000,
                                   "RebinnedShift_MeV"));
  TGraphAsymmErrors *grSideBandSumRebinShift = DreamCF::AddCF(
      SideBandSumRebinShift,
      { { fSideBands[0], fSideBands[1], fSideBands[2], fSideBands[3]} },
      Form("Gr%s", SideBandSumRebinShift->GetName()));
  fSidebandCFGr.push_back(grSideBandSumRebinShift);
  fSidebandCFGr.push_back(DreamCF::ConvertToOtherUnit(grSideBandSumRebinShift, 1000,
                                                      "RebinnedShift_MeV"));

  TH1F* SideBandSumReweighted = AddCF(CF1Reweighted, CF2Reweighted,
                                      CF3Reweighted, CF4Reweighted,
                                      "Reweighted");
  fSideBandCFs.push_back(SideBandSumReweighted);
  fSideBandCFs.push_back(
      DreamCF::ConvertToOtherUnit(SideBandSumReweighted, 1000,
                                   "Reweighted_MeV"));
  TGraphAsymmErrors *grSideBandSumReweighted = DreamCF::AddCF(
      SideBandSumReweighted,
      { { fSideBands[0], fSideBands[1], fSideBands[2], fSideBands[3]} },
      Form("Gr%s", SideBandSumReweighted->GetName()));
  fSidebandCFGr.push_back(grSideBandSumReweighted);
  fSidebandCFGr.push_back(DreamCF::ConvertToOtherUnit(grSideBandSumReweighted, 1000,
                                                      "Reweighted_MeV"));
}

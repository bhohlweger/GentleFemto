#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"
#include <iostream>

int main(int argc, char* argv[]) {
//  HM
 std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
  // std::vector<float> mTppBins = { 1.08, 1.28, 1.4, 1.58, 1.88, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
  double norm1 = 0.;
  double norm2 = 6.;
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(norm1, norm2);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  // p-antiL
  DreamDist* pALCommon = DreamFile->GetPairDistributionsCommon(0, 3, "");
  DreamDist* pALNonCommon = DreamFile->GetPairDistributionsNonCommon(0, 3, "");

  DreamDist* ApLCommon = DreamFile->GetPairDistributionsCommon(1, 2, "");
  DreamDist* ApLNonCommon = DreamFile->GetPairDistributionsNonCommon(1, 2, "");

  DreamCF* CFpALDefCommon = CATSinput->ObtainCFSyst(4, "pALDefCommon", pALCommon, ApLCommon);//rebin = 4
  DreamCF* CFpALDefNonCommon = CATSinput->ObtainCFSyst(4, "pALDefNonCommon", pALNonCommon, ApLNonCommon);//rebin = 4

  const int pairCountsDefaultCommon = CFpALDefCommon->GetFemtoPairs(0, 0.2);
  const int pairCountsDefaultNonCommon = CFpALDefNonCommon->GetFemtoPairs(0, 0.2);

    TH1F* DefRebinCommon = CFpALDefCommon->FindCorrelationFunction("hCk_RebinnedpALDefCommon_0");
    TH1F* DefRebinMeVCommon = CFpALDefCommon->FindCorrelationFunction(
        "hCk_RebinnedpALDefCommonMeV_0");
    TH1F* DefReweightedCommon = CFpALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpALDefCommon_1");
    TH1F* DefReweightedMeVCommon = CFpALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpALDefCommonMeV_1");

    TH1F* DefRebinNonCommon = CFpALDefNonCommon->FindCorrelationFunction("hCk_RebinnedpALDefNonCommon_0");
    TH1F* DefRebinMeVNonCommon = CFpALDefNonCommon->FindCorrelationFunction(
        "hCk_RebinnedpALDefNonCommonMeV_0");
    TH1F* DefReweightedNonCommon = CFpALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpALDefNonCommon_1");
    TH1F* DefReweightedMeVNonCommon = CFpALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpALDefNonCommonMeV_1");

    DreamFile->SetQuite();
    DreamKayTee* mTpALDistsCommon;
    DreamKayTee* mTpALDistsNonCommon;

    DreamFile->ReadmTHistosAncestors(filename, prefix, "8");
    mTpALDistsCommon = DreamFile->GetmTPairDistributionsCommon(0, 3, 1, 2);
    mTpALDistsNonCommon = DreamFile->GetmTPairDistributionsNonCommon(0, 3, 1, 2);

    mTpALDistsCommon->SetSEMEReweightingRatio(DefRebinCommon, DefReweightedCommon, DefRebinMeVCommon,
                                       DefReweightedMeVCommon);
    mTpALDistsCommon->SetKayTeeBins(mTppBins);
    mTpALDistsCommon->SetNormalization(norm1, norm2);
    mTpALDistsCommon->SetRebin( { 4 });
    mTpALDistsCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpALDistsCommon->ObtainTheCorrelationFunctionAncestors(
        gSystem->pwd(), prefix, TString::Format("pALDefCommon"), TString::Format("Common"));

    mTpALDistsNonCommon->SetSEMEReweightingRatio(DefRebinNonCommon, DefReweightedNonCommon, DefRebinMeVNonCommon,
                                       DefReweightedMeVNonCommon);
    mTpALDistsNonCommon->SetKayTeeBins(mTppBins);
    mTpALDistsNonCommon->SetNormalization(norm1, norm2);
    mTpALDistsNonCommon->SetRebin( { 4 });
    mTpALDistsNonCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpALDistsNonCommon->ObtainTheCorrelationFunctionAncestors(
        gSystem->pwd(), prefix, TString::Format("pALDefNonCommon"), TString::Format("NonCommon"));

  CFpALDefCommon->WriteOutput("CFPALCommon.root");
  CFpALDefNonCommon->WriteOutput("CFPALNonCommon.root");

}

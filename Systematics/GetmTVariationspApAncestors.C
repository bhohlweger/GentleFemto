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
//  std::vector<float> mTppBins = { 1.08, 1.32, 1.65, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.180, 0.280);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");

  DreamDist* pApCommon = DreamFile->GetPairDistributionsCommon(0, 1, "");
  DreamDist* pApNonCommon = DreamFile->GetPairDistributionsNonCommon(0, 1, "");

  DreamCF* CFpApDefCommon = CATSinput->ObtainCFSystBBar(1, "pApDefCommon", pApCommon, nullptr);//rebin = 1
  DreamCF* CFpApDefNonCommon = CATSinput->ObtainCFSystBBar(1, "pApDefNonCommon", pApNonCommon, nullptr);//rebin = 1


  const int pairCountsDefaultCommon = CFpApDefCommon->GetFemtoPairs(0, 0.2);
  const int pairCountsDefaultNonCommon = CFpApDefNonCommon->GetFemtoPairs(0, 0.2);

    TH1F* DefRebinCommon = CFpApDefCommon->FindCorrelationFunction("hCk_RebinnedpApDefCommon_0");
    TH1F* DefRebinMeVCommon = CFpApDefCommon->FindCorrelationFunction(
        "hCk_RebinnedpApDefCommonMeV_0");
    TH1F* DefReweightedCommon = CFpApDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefCommon_0");
    TH1F* DefReweightedMeVCommon = CFpApDefCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefCommonMeV_0");

    TH1F* DefRebinNonCommon = CFpApDefNonCommon->FindCorrelationFunction("hCk_RebinnedpApDefNonCommon_0");
    TH1F* DefRebinMeVNonCommon = CFpApDefNonCommon->FindCorrelationFunction(
        "hCk_RebinnedpApDefNonCommonMeV_0");
    TH1F* DefReweightedNonCommon = CFpApDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefNonCommon_0");
    TH1F* DefReweightedMeVNonCommon = CFpApDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedpApDefNonCommonMeV_0");

    DreamFile->SetQuite();
    DreamKayTee* mTpApDistsCommon;
    DreamKayTee* mTpApDistsNonCommon;


    DreamFile->ReadmTHistosAncestors(filename, prefix, "8");
    mTpApDistsCommon = DreamFile->GetmTPairDistributionsCommon(0, 1);
    mTpApDistsNonCommon = DreamFile->GetmTPairDistributionsNonCommon(0, 1);

    mTpApDistsCommon->SetSEMEReweightingRatio(DefRebinCommon, DefReweightedCommon, DefRebinMeVCommon,
                                       DefReweightedMeVCommon);
    mTpApDistsCommon->SetKayTeeBins(mTppBins);
    mTpApDistsCommon->SetNormalization(0.18, 0.28);
    mTpApDistsCommon->SetRebin( { 1 });
    mTpApDistsCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpApDistsCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("pApDefCommon"), TString::Format("Common"));

    mTpApDistsNonCommon->SetSEMEReweightingRatio(DefRebinNonCommon, DefReweightedNonCommon, DefRebinMeVNonCommon,
                                       DefReweightedMeVNonCommon);
    mTpApDistsNonCommon->SetKayTeeBins(mTppBins);
    mTpApDistsNonCommon->SetNormalization(0.18, 0.28);
    mTpApDistsNonCommon->SetRebin( { 1 });
    mTpApDistsNonCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpApDistsNonCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("pApDefNonCommon"), TString::Format("NonCommon"));

  CFpApDefCommon->WriteOutput("CFpApCommon.root");
  CFpApDefNonCommon->WriteOutput("CFpApNonCommon.root");

}

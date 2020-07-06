#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"
#include <iostream>

int main(int argc, char* argv[]) {
//  HM
//  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
 std::vector<float> mTppBins = { 1.08, 1.32, 1.65, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  double norm1 = 0.;
  double norm2 = 6.;
  CATSinput->SetNormalization(norm1, norm2);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  // p-antiL
  DreamDist* LALCommon = DreamFile->GetPairDistributionsCommon(2, 3, "");
  DreamDist* LALNonCommon = DreamFile->GetPairDistributionsNonCommon(2, 3, "");

  DreamCF* CFLALDefCommon = CATSinput->ObtainCFSystBBar(5, "LALDefCommon", LALCommon, nullptr);//rebin = 5
  DreamCF* CFLALDefNonCommon = CATSinput->ObtainCFSystBBar(5, "LALDefNonCommon", LALNonCommon, nullptr);//rebin = 5


  const int pairCountsDefaultCommon = CFLALDefCommon->GetFemtoPairs(0, 0.2);
  const int pairCountsDefaultNonCommon = CFLALDefNonCommon->GetFemtoPairs(0, 0.2);

    TH1F* DefRebinCommon = CFLALDefCommon->FindCorrelationFunction("hCk_RebinnedLALDefCommon_0");
    TH1F* DefRebinMeVCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_RebinnedLALDefCommonMeV_0");
    TH1F* DefReweightedCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedLALDefCommon_1");
    TH1F* DefReweightedMeVCommon = CFLALDefCommon->FindCorrelationFunction(
        "hCk_ReweightedLALDefCommonMeV_1");

    TH1F* DefRebinNonCommon = CFLALDefNonCommon->FindCorrelationFunction("hCk_RebinnedLALDefNonCommon_0");
    TH1F* DefRebinMeVNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_RebinnedLALDefNonCommonMeV_0");
    TH1F* DefReweightedNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedLALDefNonCommon_1");
    TH1F* DefReweightedMeVNonCommon = CFLALDefNonCommon->FindCorrelationFunction(
        "hCk_ReweightedLALDefNonCommonMeV_1");

    DreamFile->SetQuite();
    DreamKayTee* mTLALDistsCommon;
    DreamKayTee* mTLALDistsNonCommon;


    DreamFile->ReadmTHistosAncestors(filename, prefix, "8");
    mTLALDistsCommon = DreamFile->GetmTPairDistributionsCommon(2, 3);
    mTLALDistsNonCommon = DreamFile->GetmTPairDistributionsNonCommon(2, 3);

    mTLALDistsCommon->SetSEMEReweightingRatio(DefRebinCommon, DefReweightedCommon, DefRebinMeVCommon,
                                       DefReweightedMeVCommon);
    mTLALDistsCommon->SetKayTeeBins(mTppBins);
    mTLALDistsCommon->SetNormalization(norm1, norm2);
    mTLALDistsCommon->SetRebin( { 5 });
    mTLALDistsCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLALDistsCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("LALDefCommon"), TString::Format("Common"));

    mTLALDistsNonCommon->SetSEMEReweightingRatio(DefRebinNonCommon, DefReweightedNonCommon, DefRebinMeVNonCommon,
                                       DefReweightedMeVNonCommon);
    mTLALDistsNonCommon->SetKayTeeBins(mTppBins);
    mTLALDistsNonCommon->SetNormalization(norm1, norm2);
    mTLALDistsNonCommon->SetRebin( { 5 });
    mTLALDistsNonCommon->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLALDistsNonCommon->ObtainTheCorrelationFunctionAncestorsSingle(
        gSystem->pwd(), prefix, TString::Format("LALDefNonCommon"), TString::Format("NonCommon"));

  CFLALDefCommon->WriteOutput("CFLALCommon.root");
  CFLALDefNonCommon->WriteOutput("CFLALNonCommon.root");

}

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
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.180, 0.280);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  // p-antiL
  DreamDist* LAL = DreamFile->GetPairDistributions(2, 3, "");
  DreamCF* CFLALDef = CATSinput->ObtainCFSystBBar(5, "LALDef", LAL, nullptr);//rebin = 5
  // p-L
  DreamDist* LL = DreamFile->GetPairDistributions(2, 2, "");
  DreamDist* ALAL = DreamFile->GetPairDistributions(3, 3, "");
  DreamCF* CFLLDef = CATSinput->ObtainCFSyst(5, "LLDef", LL, ALAL);//rebin = 5
  const int pairCountsDefault = CFLALDef->GetFemtoPairs(0, 0.2);
  const int pairCountsDefault_LL = CFLLDef->GetFemtoPairs(0, 0.2);

    TH1F* DefRebin = CFLALDef->FindCorrelationFunction("hCk_RebinnedLALDef_0");
    TH1F* DefRebinMeV = CFLALDef->FindCorrelationFunction(
        "hCk_RebinnedLALDefMeV_0");
    TH1F* DefReweighted = CFLALDef->FindCorrelationFunction(
        "hCk_ReweightedLALDef_1");
    TH1F* DefReweightedMeV = CFLALDef->FindCorrelationFunction(
        "hCk_ReweightedLALDefMeV_1");

    TH1F* DefRebin_LL = CFLLDef->FindCorrelationFunction("hCk_RebinnedLLDef_0");
    TH1F* DefRebinMeV_LL = CFLLDef->FindCorrelationFunction(
        "hCk_RebinnedLLDefMeV_0");
    TH1F* DefReweighted_LL = CFLLDef->FindCorrelationFunction(
        "hCk_ReweightedLLDef_1");
    TH1F* DefReweightedMeV_LL = CFLLDef->FindCorrelationFunction(
        "hCk_ReweightedLLDefMeV_1");

    DreamFile->SetQuite();
    DreamKayTee* mTLLDists;
    DreamKayTee* mTLALDists;

    DreamFile->ReadmTHistos(filename, prefix, "8");
    mTLALDists = DreamFile->GetmTPairDistributionsBBar(2, 3);
    mTLLDists = DreamFile->GetmTPairDistributions(2, 2, 3, 3);

    mTLALDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTLALDists->SetKayTeeBins(mTppBins);
    mTLALDists->SetNormalization(0.18, 0.28);
    mTLALDists->SetRebin( { 5 });
    mTLALDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLALDists->ObtainTheCorrelationFunctionBBar(
        gSystem->pwd(), prefix, TString::Format("LALDef"));


    mTLLDists->SetSEMEReweightingRatio(DefRebin_LL, DefReweighted_LL, DefRebinMeV_LL,
                                       DefReweightedMeV_LL);
    mTLLDists->SetKayTeeBins(mTppBins);
    mTLLDists->SetNormalization(0.18, 0.28);
    mTLLDists->SetRebin( { 5 });
    mTLLDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTLLDists->ObtainTheCorrelationFunction(
        gSystem->pwd(), prefix, TString::Format("LLDef"));

  CFLALDef->WriteOutput("CFLAL.root");
  CFLLDef->WriteOutput("CFLL.root");

}

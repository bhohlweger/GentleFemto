#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"

int main(int argc, char* argv[]) {
//  HM
//  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
  std::vector<float> mTppBins = { 1.08, 1.28, 1.4, 1.58, 1.88, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.180, 0.280);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, "8");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "8");
  DreamDist* pAL = DreamFile->GetPairDistributions(0, 3, "");
  DreamDist* ApL = DreamFile->GetPairDistributions(1, 2, "");
  DreamCF* CFpALDef = CATSinput->ObtainCFSyst(4, "pALDef", pAL, ApL);//rebin = 4
  const int pairCountsDefault = CFpALDef->GetFemtoPairs(0, 0.2);

    TH1F* DefRebin = CFpALDef->FindCorrelationFunction("hCk_RebinnedpALDef_0");
    TH1F* DefRebinMeV = CFpALDef->FindCorrelationFunction(
        "hCk_RebinnedpALDefMeV_0");
    TH1F* DefReweighted = CFpALDef->FindCorrelationFunction(
        "hCk_ReweightedpALDef_1");
    TH1F* DefReweightedMeV = CFpALDef->FindCorrelationFunction(
        "hCk_ReweightedpALDefMeV_1");

    DreamFile->SetQuite();
    DreamKayTee* mTpLDists;
    DreamFile->ReadmTHistos(filename, prefix, "8");
    mTpLDists = DreamFile->GetmTPairDistributions(0, 3, 1, 2);
    mTpLDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTpLDists->SetKayTeeBins(mTppBins);
    mTpLDists->SetNormalization(0.18, 0.28);
    mTpLDists->SetRebin( { 4 });
    mTpLDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpLDists->ObtainTheCorrelationFunction(
        gSystem->pwd(), prefix, TString::Format("pALDef"));

  CFpALDef->WriteOutput("CFPAL.root");
}

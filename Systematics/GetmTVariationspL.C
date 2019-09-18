#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"

int main(int argc, char* argv[]) {
  std::vector<float> mTppBins = { 1.08, 1.3,1.5, 4.5 };
//  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
//  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.240, 0.340);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "0");
  DreamDist* pL = DreamFile->GetPairDistributions(0, 2, "");
  DreamDist* ApAL = DreamFile->GetPairDistributions(1, 3, "");
  DreamCF* CFpLDef = CATSinput->ObtainCFSyst(5, "pLDef", pL, ApAL);
  const int pairCountsDefault = CFpLDef->GetFemtoPairs(0, 0.2);
  int OutCounter = 0;  // 0 is going to be the DEFAULT
  for (int iVar = 0; iVar < 45; ++iVar) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename, prefix, Form("%u", iVar));
    TString VarName = TString::Format("pLVar%u", iVar);
    DreamCF* CFpLVar = CATSinput->ObtainCFSyst(
        5, VarName.Data(), DreamVarFile->GetPairDistributions(0, 2, ""),
        DreamVarFile->GetPairDistributions(1, 3, ""));
    int femtoPairVar = CFpLVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    TH1F* DefRebin = CFpLDef->FindCorrelationFunction("hCk_RebinnedpLDef_0");
    TH1F* DefRebinMeV = CFpLDef->FindCorrelationFunction(
        "hCk_RebinnedpLDefMeV_0");
    TH1F* DefReweighted = CFpLDef->FindCorrelationFunction(
        "hCk_ReweightedpLDef_1");
    TH1F* DefReweightedMeV = CFpLDef->FindCorrelationFunction(
        "hCk_ReweightedpLDefMeV_1");
    ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
    DreamFile->SetQuite();
    DreamKayTee* mTpLDists;
    DreamFile->ReadmTHistos(filename, prefix, TString::Format("%u", iVar));
    mTpLDists = DreamFile->GetmTPairDistributions(0, 2, 1, 3);
    mTpLDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTpLDists->SetKayTeeBins(mTppBins);
    mTpLDists->SetNormalization(0.24, 0.34);
    mTpLDists->SetRebin( { 5 });
    mTpLDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpLDists->ObtainTheCorrelationFunction(
        gSystem->pwd(), prefix, TString::Format("pLVar%u", OutCounter));
    OutCounter++;
  }
  CFpLDef->WriteOutput("CFPL.root");
}

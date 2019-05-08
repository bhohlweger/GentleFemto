#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"
int main(int argc, char* argv[]) {
  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5 };
  const char* filename = argv[1];
  const char* prefix = argv[2];
  auto CATSinput = new CATSInput();
  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.240, 0.340);

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "0");
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(1, "ppDef", pp, ApAp);

  const int pairCountsDefault = CFppDef->GetFemtoPairs(0, 0.2);
  int OutCounter = 0; // 0 is going to be the DEFAULT
  for (int iVar = 0; iVar < 45; ++iVar) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename, prefix, Form("%u", iVar));
    TString VarName = TString::Format("ppVar%u", iVar);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        1, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    int femtoPairVar = CFppVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    TH1F* DefRebin = CFppDef->FindCorrelationFunction("hCk_FixShiftedppDef_0");
    TH1F* DefRebinMeV = CFppDef->FindCorrelationFunction(
        "hCk_FixShiftedppDefMeV_0");
    TH1F* DefReweighted = CFppDef->FindCorrelationFunction(
        "hCk_ReweightedppDef_0");
    TH1F* DefReweightedMeV = CFppDef->FindCorrelationFunction(
        "hCk_ReweightedppDefMeV_0");
    ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
    DreamFile->SetQuite();
    DreamKayTee* mTppDists;
    DreamFile->ReadmTHistos(filename, prefix, TString::Format("%u", iVar));
    mTppDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

    mTppDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTppDists->SetKayTeeBins(mTppBins);
    mTppDists->SetNormalization(0.24, 0.34);
//    mTppDists->SetRebin( { 2 });
    mTppDists->FixShift( { true, true, true, true, true, true, true }, { 0.004,
                            0.004, 0.004, 0.004, 0.004, 0.004, 0.004 });
    mTppDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix,
                                            TString::Format("pp_%u", OutCounter));
    OutCounter++;
  }
}

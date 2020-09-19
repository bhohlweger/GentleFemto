#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"

#include <iostream>
int main(int argc, char *argv[]) {
// MB
  std::vector<float> mTppBins = { 1.02, 1.14, 1.26, 4.5 };
//  HM
  // std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5 };
//  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 4.5 };
  const char *filename = argv[1];
  const char *prefix = argv[2];
  auto CATSinput = new CATSInput();
  CATSinput->SetFixedkStarMinBin(true, 0.004);
  CATSinput->SetNormalization(0.240, 0.340);

  ReadDreamFile *DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, "0");

  ForgivingReader *ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "0");
  DreamDist *pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist *ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF *CFppDef = CATSinput->ObtainCFSyst(1, "ppDef", pp, ApAp);

  const int pairCountsDefault = CFppDef->GetFemtoPairs(0, 0.2);
  int OutCounter = 0;  // 0 is going to be the DEFAULT
  for (int iVar = 0; iVar < 44; ++iVar) {
    ReadDreamFile *DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename, prefix, Form("%u", iVar));
    TString VarName = TString::Format("ppVar%u", iVar);
    DreamCF *CFppVar = CATSinput->ObtainCFSyst(
        1, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    int femtoPairVar = CFppVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    std::cout << "Hey buddy look at the variation = " << TMath::Abs(relDiff)
        << std::endl;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    std::cout << "THIS VARIATION HAS BEEN SELECTED = " << iVar << std::endl;

    DreamFile->ReadmTMultHistos(filename, prefix, TString::Format("%u", iVar),
                                3);
    DreamKayTee *mTBinner = DreamFile->GetmTMultPairDistributions(0, 0, 1, 1,
                                                                  3);
    mTBinner->SetKayTeeBins(mTppBins);
    mTBinner->SetMultBins( { 8 });
    mTBinner->FixShift( { true, true, true }, { 0.008, 0.008, 0.012 });
    std::vector<DreamCF*> mTBin_0 = mTBinner->GetmTMultBinned(0, OutCounter);
    std::vector<DreamCF*> mTBin_1 = mTBinner->GetmTMultBinned(1, OutCounter);
    std::vector<DreamCF*> mTBin_2 = mTBinner->GetmTMultBinned(2, OutCounter);

    OutCounter++;
  }
}

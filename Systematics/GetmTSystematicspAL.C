#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "ForgivingReader.h"

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
  DreamFile->SetAnalysisFile(filename, prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename, prefix,
                                                          "0");
  // p-antiL
  DreamDist* pAL = DreamFile->GetPairDistributions(0, 3, "");
  DreamDist* ApL = DreamFile->GetPairDistributions(1, 2, "");
  DreamCF* CFpALDef = CATSinput->ObtainCFSyst(4, "pALDef", pAL, ApL);//rebin = 4
  // p-L
  DreamDist* pL = DreamFile->GetPairDistributions(0, 2, "");
  DreamDist* ApAL = DreamFile->GetPairDistributions(1, 3, "");
  DreamCF* CFpLDef = CATSinput->ObtainCFSyst(4, "pLDef", pL, ApAL);//rebin = 4
  const int pairCountsDefault = CFpALDef->GetFemtoPairs(0, 0.2);
  const int pairCountsDefault_pL = CFpLDef->GetFemtoPairs(0, 0.2);

//Running Systematics only for p-antiL:
  int OutCounter = 0;  // 0 is going to be the DEFAULT
  for (int iVar = 0; iVar < 45; ++iVar) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename, prefix, Form("%u", iVar));
    TString VarName = TString::Format("pALVar%u", iVar);
    DreamCF* CFpALVar = CATSinput->ObtainCFSyst(
        4, VarName.Data(), DreamVarFile->GetPairDistributions(0, 3, ""),
        DreamVarFile->GetPairDistributions(1, 2, ""));
    int femtoPairVar = CFpALVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    TH1F* DefRebin = CFpALDef->FindCorrelationFunction("hCk_RebinnedpALDef_0");
    TH1F* DefRebinMeV = CFpALDef->FindCorrelationFunction(
        "hCk_RebinnedpALDefMeV_0");
    TH1F* DefReweighted = CFpALDef->FindCorrelationFunction(
        "hCk_ReweightedpALDef_1");
    TH1F* DefReweightedMeV = CFpALDef->FindCorrelationFunction(
        "hCk_ReweightedpALDefMeV_1");
    ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
    DreamFile->SetQuite();
    DreamKayTee* mTpALDists;
    DreamFile->ReadmTHistos(filename, prefix, TString::Format("%u", iVar));
    mTpALDists = DreamFile->GetmTPairDistributions(0, 3, 1, 2);
    mTpALDists->SetSEMEReweightingRatio(DefRebin, DefReweighted, DefRebinMeV,
                                       DefReweightedMeV);
    mTpALDists->SetKayTeeBins(mTppBins);
    mTpALDists->SetNormalization(0.18, 0.28);
    mTpALDists->SetRebin( { 4 });
    mTpALDists->FixShift( { true, true, true, true, true, true, true }, { 0., 0.,
                            0., 0., 0., 0., 0. });
    mTpALDists->ObtainTheCorrelationFunction(
        gSystem->pwd(), prefix, TString::Format("pALVar%u", OutCounter));
    OutCounter++;
  }

  CFpALDef->WriteOutput("CFPAL.root");

}

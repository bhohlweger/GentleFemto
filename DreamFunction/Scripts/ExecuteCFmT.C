#include "DreamKayTee.h"
#include "ReadDreamFile.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* CalibName = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamKayTee* mTDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

  std::vector<float> mTBins = { 1.01,1.07,1.13,1.19,1.25,1.37,1.55,1.97,4.5 };
//  std::vector<float> mTBins = { 1.01,1.13,1.25,1.55,4.5 };

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  mTDists->SetSEMEReweightingRatio(CalibName,"pp");
  mTDists->SetKayTeeBins(mTBins);
  mTDists->SetNormalization(0.2, 0.4);
  mTDists->ObtainTheCorrelationFunction(foldername.Data(), prefix, "pp");

  return 1;
}

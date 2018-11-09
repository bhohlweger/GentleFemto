#include "DreamKayTee.h"
#include "ReadDreamFile.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamKayTee* kTDists;
  DreamFile->ReadkTHistos(filename, prefix, addon);
  kTDists = DreamFile->GetkTPairDistributions(0, 0, 1, 1);

  std::vector<float> kTBins = { 0.48, 0.69, 1., 1.5 };

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  kTDists->SetKayTeeBins(kTBins);
  kTDists->SetNormalization(0.2, 0.4);
  kTDists->ObtainTheCorrelationFunction(foldername.Data(), prefix, "pp");

  return 1;
}

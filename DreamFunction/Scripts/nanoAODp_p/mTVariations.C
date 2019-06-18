#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* CalibFile = argv[2];
  const char* prefix = argv[3];
  const char* addon = argv[4];
  const char* outname = argv[5];
  const char* outpath = argv[6];

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetQuite();
  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  DreamKayTee* mTppDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTppDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5};
  TString CalibPP = Form("%s", CalibFile);
  mTppDists->SetSEMEReweightingRatio(CalibPP, "1", false);
  mTppDists->SetKayTeeBins(mTppBins);
  mTppDists->SetNormalization(0.24, 0.34);
  mTppDists->SetRebin( { 2 });
  mTppDists->FixShift({true,true,true,true,true,true,true},{0.004,0.004,0.008,0.004,0.008,0.008,0.012});
  mTppDists->ObtainTheCorrelationFunction(outpath, prefix, "pp");
  return 1;
}

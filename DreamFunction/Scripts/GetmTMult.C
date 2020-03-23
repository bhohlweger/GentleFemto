#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "DreamKayTee.h"
#include <iostream>

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  //gStyle->SetOptStat(0);

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReadmTMultHistos(filename, prefix, addon,3);
  DreamKayTee* mTBinner = DreamFile->GetmTMultPairDistributions(0, 0, 1, 1, 3); 
  mTBinner->SetMultBins({9});
  std::vector<DreamCF*> mTBin_0 = mTBinner->GetmTMultBinned(0);
  std::vector<DreamCF*> mTBin_1 = mTBinner->GetmTMultBinned(1);
  std::vector<DreamCF*> mTBin_2 = mTBinner->GetmTMultBinned(2); 
  
  return 0; 
}





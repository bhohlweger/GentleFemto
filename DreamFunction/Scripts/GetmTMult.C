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
  int cout = 0; //counter associated to variations in systematics
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReadmTMultHistos(filename, prefix, addon,3);
  DreamKayTee* mTBinner = DreamFile->GetmTMultPairDistributions(0, 0, 1, 1, 3); 
  mTBinner->SetMultBins({8});
  mTBinner->FixShift({true,true,true},{0.008, 0.008, 0.012});
  std::vector<DreamCF*> mTBin_0 = mTBinner->GetmTMultBinned(0, cout);
  std::vector<DreamCF*> mTBin_1 = mTBinner->GetmTMultBinned(1, cout);
  std::vector<DreamCF*> mTBin_2 = mTBinner->GetmTMultBinned(2, cout );
  
  return 0; 
}





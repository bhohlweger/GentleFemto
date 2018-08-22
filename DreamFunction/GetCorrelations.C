#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
//  gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile=new ReadDreamFile(6,6);
  DreamFile->SetAnalysisFile(filename,prefix);

  DreamCF* pXi=new DreamCF();
  pXi->SetPair(DreamFile->GetPairDistributions(0,4,""));
  pXi->ShiftForEmpty();
//  DreamCF* ApAXi=new DreamCF();
//  ApAXi->SetPair(DreamFile->GetPairDistributions(1,5,""));
//  ApAXi->ShiftForEmpty();
}

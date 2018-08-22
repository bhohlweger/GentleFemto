#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
//  gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile=new ReadDreamFile(6,6);
  DreamFile->SetAnalysisFile(filename,prefix);

  DreamPair* pXi=new DreamPair();
  DreamPair* ApAXi=new DreamPair();
  pXi->SetPair(DreamFile->GetPairDistributions(0,4,""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1,5,""));
  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());
  for (int iReb=3;iReb<6;++iReb)
  {
    pXi->Rebin(pXi->GetPairShifted(0),iReb);
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb-3),0.,3.0);
    ApAXi->Rebin(ApAXi->GetPairShifted(0),iReb);
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb-3),0.,3.0);
  }
}

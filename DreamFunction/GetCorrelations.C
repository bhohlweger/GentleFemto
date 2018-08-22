#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
//  gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part",0.2,0.4);
  DreamPair* ApAXi = new DreamPair("AntiPart",0.2,0.4);
  pXi->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 5, ""));

  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());

  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                  pXi->GetFirstBin());

  for (int iReb = 3; iReb < 6; ++iReb) {
    pXi->Rebin(pXi->GetPairFixShifted(0), iReb);
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb - 3), 0., 3.0);
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), iReb);
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb - 3), 0., 3.0);
  }
  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations("CFOutput_pXi.root");
  //Write to file
}

#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
//  gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part",0.2,0.4);
  DreamPair* ApAXi = new DreamPair("AntiPart",0.2,0.4);
  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pXi->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 5, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());
  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                  pXi->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
  for (int iReb = 3; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pXi->Rebin(pXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb - 3), 0., 3.0);
    std::cout << "==Rebinning==" << std::endl;
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb - 3), 0., 3.0);
  }
  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;
  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  CF_pXi->WriteOutput("CFOutput_pXi.root");
}

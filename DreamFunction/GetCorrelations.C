#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.2, 0.4);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.2, 0.4);

  DreamCF* CF_pL = new DreamCF();
  DreamPair* pL = new DreamPair("Part", 0.2, 0.4);
  DreamPair* ApAL = new DreamPair("AntiPart", 0.2, 0.4);

  DreamCF* CF_LL = new DreamCF();
  DreamPair* LL = new DreamPair("Part", 0.4, 0.6);
  DreamPair* ALAL = new DreamPair("AntiPart", 0.4, 0.6);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part", 0.2, 0.4);
  DreamPair* ApAXi = new DreamPair("AntiPart", 0.2, 0.4);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  pL->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApAL->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

  LL->SetPair(DreamFile->GetPairDistributions(2, 2, ""));
  ALAL->SetPair(DreamFile->GetPairDistributions(3, 3, ""));

  pXi->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 5, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pL->ShiftForEmpty(pL->GetPair());
  ApAL->ShiftForEmpty(ApAL->GetPair());

  LL->ShiftForEmpty(LL->GetPair());
  ALAL->ShiftForEmpty(ALAL->GetPair());

  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());

  pL->FixShift(pL->GetPairShiftedEmpty(0), ApAL->GetPairShiftedEmpty(0),
               ApAL->GetFirstBin());
  ApAL->FixShift(ApAL->GetPairShiftedEmpty(0), pL->GetPairShiftedEmpty(0),
                 pL->GetFirstBin());
  LL->FixShift(LL->GetPairShiftedEmpty(0), ALAL->GetPairShiftedEmpty(0),
               ALAL->GetFirstBin());
  ALAL->FixShift(ALAL->GetPairShiftedEmpty(0), LL->GetPairShiftedEmpty(0),
                 LL->GetFirstBin());

  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                  pXi->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 3; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pL->Rebin(pL->GetPairFixShifted(0), iReb);
    LL->Rebin(LL->GetPairFixShifted(0), iReb);
    pXi->Rebin(pXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pL->ReweightMixedEvent(pL->GetPairRebinned(iReb - 3), 0., 3.0);
    LL->ReweightMixedEvent(LL->GetPairRebinned(iReb - 3), 0., 3.0);
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb - 3), 0., 3.0);
    std::cout << "==Rebinning==" << std::endl;
    ApAL->Rebin(ApAL->GetPairFixShifted(0), iReb);
    ALAL->Rebin(ALAL->GetPairFixShifted(0), iReb);
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    ApAL->ReweightMixedEvent(ApAL->GetPairRebinned(iReb - 3), 0., 3.0);
    ALAL->ReweightMixedEvent(ALAL->GetPairRebinned(iReb - 3), 0., 3.0);
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb - 3), 0., 3.0);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0., 3.0);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0., 3.0);

  pL->Rebin(pL->GetPair(), 4);
  pL->Rebin(pL->GetPair(), 5);
  ApAL->Rebin(ApAL->GetPair(), 4);
  ApAL->Rebin(ApAL->GetPair(), 5);
  LL->Rebin(LL->GetPair(), 4);
  LL->Rebin(LL->GetPair(), 5);
  ALAL->Rebin(ALAL->GetPair(), 4);
  ALAL->Rebin(ALAL->GetPair(), 5);
  pXi->Rebin(pXi->GetPair(), 4);
  pXi->Rebin(pXi->GetPair(), 5);
  ApAXi->Rebin(pXi->GetPair(), 4);
  ApAXi->Rebin(pXi->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;
  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput("CFOutput_pp.root");

  CF_pL->SetPairs(pL, ApAL);
  CF_pL->GetCorrelations();
  CF_pL->WriteOutput("CFOutput_pL.root");

  CF_LL->SetPairs(LL, ALAL);
  CF_LL->GetCorrelations();
  CF_LL->WriteOutput("CFOutput_LL.root");

  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  CF_pXi->WriteOutput("CFOutput_pXi.root");
}

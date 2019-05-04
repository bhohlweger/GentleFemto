#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelations(const float fixShift, const char* filename,
                     const char* prefix, const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pL = new DreamCF();
  DreamPair* pL = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAL = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_LL = new DreamCF();
  DreamPair* LL = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ALAL = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAXi = new DreamPair("AntiPart", 0.24, 0.34);

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
  std::cout << "The old one \n";
  std::cout << "p-p \n";
  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  std::cout << "Ap-Ap \n";
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());
  std::cout << "The new one \n";
  std::cout << "p-p \n";
  pp->FixShift(pp->GetPair(), ApAp->GetPair(), 0.004, true);
  std::cout << "Ap-Ap \n";
  ApAp->FixShift(ApAp->GetPair(), pp->GetPair(), 0.004, true);
  std::cout << "over \n";
  pL->FixShift(pL->GetPairShiftedEmpty(0), ApAL->GetPairShiftedEmpty(0),
               ApAL->GetFirstBin());
  ApAL->FixShift(ApAL->GetPairShiftedEmpty(0), pL->GetPairShiftedEmpty(0),
                 pL->GetFirstBin());
  LL->FixShift(LL->GetPairShiftedEmpty(0), ALAL->GetPairShiftedEmpty(0),
               ALAL->GetFirstBin());
  ALAL->FixShift(ALAL->GetPairShiftedEmpty(0), LL->GetPairShiftedEmpty(0),
                 LL->GetFirstBin());
  if (fixShift < 1e-6) {
    pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                  ApAXi->GetFirstBin());
    ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                    pXi->GetFirstBin());
  } else {
    pXi->FixShift(pXi->GetPair(), nullptr, fixShift, true);
    ApAXi->FixShift(ApAXi->GetPair(), nullptr, fixShift, true);
  }
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pL->Rebin(pL->GetPairFixShifted(0), rebinVec[iReb]);
    LL->Rebin(LL->GetPairFixShifted(0), rebinVec[iReb]);
    pXi->Rebin(pXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    pL->ReweightMixedEvent(pL->GetPairRebinned(iReb), 0.2, 0.9);
    LL->ReweightMixedEvent(LL->GetPairRebinned(iReb), 0.2, 0.9);
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApAL->Rebin(ApAL->GetPairFixShifted(0), rebinVec[iReb]);
    ALAL->Rebin(ALAL->GetPairFixShifted(0), rebinVec[iReb]);
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    ApAL->ReweightMixedEvent(ApAL->GetPairRebinned(iReb), 0.2, 0.9);
    ALAL->ReweightMixedEvent(ALAL->GetPairRebinned(iReb), 0.2, 0.9);
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pp->ReweightMixedEvent(pp->GetPairFixShifted(1), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(1), 0.2, 0.9);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = TString::Format("%s",gSystem->pwd());

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%s/CFOutput_pp.root", foldername.Data()));

  CF_pL->SetPairs(pL, ApAL);
  CF_pL->GetCorrelations();
  CF_pL->WriteOutput(Form("%s/CFOutput_pL.root", foldername.Data()));

  CF_LL->SetPairs(LL, ALAL);
  CF_LL->GetCorrelations();
  CF_LL->WriteOutput(Form("%s/CFOutput_LL.root", foldername.Data()));

  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  CF_pXi->WriteOutput(Form("%s/CFOutput_pXi.root", foldername.Data()));
}

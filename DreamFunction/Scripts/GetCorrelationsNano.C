#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelations(const char* filename,
                     const char* prefix, const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAXi = new DreamPair("AntiPart", 0.24, 0.34);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  pXi->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

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

  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                  pXi->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pXi->Rebin(pXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pp->ReweightMixedEvent(pp->GetPairFixShifted(1), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(1), 0.2, 0.9);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%s/CFOutput_pp.root", gSystem->pwd()));

  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  CF_pXi->WriteOutput(Form("%s/CFOutput_pSomething.root", gSystem->pwd()));
}

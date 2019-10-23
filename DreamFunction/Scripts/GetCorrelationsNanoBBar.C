#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoBBar(const char* filename,
                     const char* prefix, const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  Double_t norm1=0.2;
  Double_t norm2=0.4;

  DreamCF* CF_pAp = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart", norm1, norm2);

  DreamCF* CF_pAL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart",norm1,norm2);
  DreamPair* ApL = new DreamPair("AntiPartPart",norm1,norm2);

  DreamCF* CF_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;

  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->ShiftForEmpty(pAp->GetPair());
  pAL->ShiftForEmpty(pAL->GetPair());
  ApL->ShiftForEmpty(ApL->GetPair());
  LAL->ShiftForEmpty(LAL->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->FixShift(pAp->GetPairShiftedEmpty(0), pAp->GetPairShiftedEmpty(0),
               pAp->GetFirstBin());

  pAL->FixShift(pAL->GetPairShiftedEmpty(0), ApL->GetPairShiftedEmpty(0),
		  	   ApL->GetFirstBin());
  ApL->FixShift(ApL->GetPairShiftedEmpty(0), pAL->GetPairShiftedEmpty(0),
		  	   pAL->GetFirstBin());

  LAL->FixShift(LAL->GetPairShiftedEmpty(0), LAL->GetPairShiftedEmpty(0),
		  	   LAL->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
// To be applied only to p-antiL and L-antiL
  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pAL->Rebin(pAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    pAL->ReweightMixedEvent(pAL->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApL->Rebin(ApL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    ApL->ReweightMixedEvent(ApL->GetPairRebinned(iReb), 0.2, 0.9);

    std::cout << "==Rebinning==" << std::endl;
    LAL->Rebin(LAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    LAL->ReweightMixedEvent(LAL->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pAp->ReweightMixedEvent(pAp->GetPairFixShifted(0), 0.2, 0.9);

//  pAL->Rebin(pAL->GetPair(), 4);
//  pAL->Rebin(pAL->GetPair(), 5);
//  ApL->Rebin(ApL->GetPair(), 4);
//  ApL->Rebin(ApL->GetPair(), 5);
//  LAL->Rebin(LAL->GetPair(), 4);
//  LAL->Rebin(LAL->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  CF_pAp->SetPairs(pAp, nullptr);
  CF_pAp->GetCorrelations();
  CF_pAp->WriteOutput(Form("%s/CFOutput_pAp_%s.root", gSystem->pwd(),addon));

  CF_pAL->SetPairs(pAL, ApL);
  CF_pAL->GetCorrelations();
  CF_pAL->WriteOutput(Form("%s/CFOutput_pAL_%s.root", gSystem->pwd(),addon));

  CF_LAL->SetPairs(LAL, nullptr);
  CF_LAL->GetCorrelations();
  CF_LAL->WriteOutput(Form("%s/CFOutput_LAL_%s.root", gSystem->pwd(),addon));
}

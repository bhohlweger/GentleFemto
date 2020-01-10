#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoBBar(const char* filename,
                     const char* prefix, const char* addon = "", double_t norm1 = 0.18, double_t norm2 = 0.28) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
    std::cout<< "debug 1" <<std::endl;

  DreamFile->SetAnalysisFile(filename, prefix, addon);

  // Double_t norm1=0.18;//default 0.2-0.4, 0.18-0.28 where CFs look flatter
  // Double_t norm2=0.28;

  DreamCF* CF_pAp = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart", norm1, norm2);

  DreamCF* CF_pAL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart",norm1,norm2);
  DreamPair* ApL = new DreamPair("AntiPartPart",norm1,norm2);

  DreamCF* CF_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);

  DreamCF* CF_pAXi = new DreamCF();
  DreamPair* pAXi = new DreamPair("PartAntiPart",norm1,norm2);
  DreamPair* ApXi = new DreamPair("AntiPartPart",norm1,norm2);

  DreamCF* CF_XiAXi = new DreamCF();
  DreamPair* XiAXi = new DreamPair("PartAntiPart", norm1, norm2);

  DreamCF* CF_LAXi = new DreamCF();
  DreamPair* LAXi = new DreamPair("PartAntiPart",norm1,norm2);
  DreamPair* ALXi = new DreamPair("AntiPartPart",norm1,norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;

  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));

  std::cout << "========Pair Set pAXi=========" << std::endl;
  pAXi->SetPair(DreamFile->GetPairDistributions(0, 5, ""));
  std::cout << "========Pair Set ApXi=========" << std::endl;
  ApXi->SetPair(DreamFile->GetPairDistributions(1, 4, ""));
  XiAXi->SetPair(DreamFile->GetPairDistributions(4, 5, ""));
  LAXi->SetPair(DreamFile->GetPairDistributions(2, 5, ""));
  ALXi->SetPair(DreamFile->GetPairDistributions(3, 4, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->ShiftForEmpty(pAp->GetPair());
  pAL->ShiftForEmpty(pAL->GetPair());
  ApL->ShiftForEmpty(ApL->GetPair());
  LAL->ShiftForEmpty(LAL->GetPair());

  pAXi->ShiftForEmpty(pAXi->GetPair());
  ApXi->ShiftForEmpty(ApXi->GetPair());
  XiAXi->ShiftForEmpty(XiAXi->GetPair());
  LAXi->ShiftForEmpty(LAXi->GetPair());
  ALXi->ShiftForEmpty(ALXi->GetPair());
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

  pAXi->FixShift(pAXi->GetPairShiftedEmpty(0), ApXi->GetPairShiftedEmpty(0),
		  	   ApXi->GetFirstBin());
  ApXi->FixShift(ApXi->GetPairShiftedEmpty(0), pAXi->GetPairShiftedEmpty(0),
		  	   pAXi->GetFirstBin());

  XiAXi->FixShift(XiAXi->GetPairShiftedEmpty(0), XiAXi->GetPairShiftedEmpty(0),
		  	   XiAXi->GetFirstBin());

  LAXi->FixShift(LAXi->GetPairShiftedEmpty(0), ALXi->GetPairShiftedEmpty(0),
		  	   ALXi->GetFirstBin());
  ALXi->FixShift(ALXi->GetPairShiftedEmpty(0), LAXi->GetPairShiftedEmpty(0),
		  	   LAXi->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
// To be applied only to p-antiL and L-antiL
  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning pAL==" << std::endl;
    pAL->Rebin(pAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting pAL==" << std::endl;
    pAL->ReweightMixedEvent(pAL->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ApL==" << std::endl;
    ApL->Rebin(ApL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ApL==" << std::endl;
    ApL->ReweightMixedEvent(ApL->GetPairRebinned(iReb), 0.2, 0.9);

    std::cout << "==Rebinning LAL==" << std::endl;
    LAL->Rebin(LAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LAL==" << std::endl;
    LAL->ReweightMixedEvent(LAL->GetPairRebinned(iReb), 0.2, 0.9);


    std::cout << "==Rebinning pAXi==" << std::endl;
    pAXi->Rebin(pAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting pAXi==" << std::endl;
    pAXi->ReweightMixedEvent(pAXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ApXi==" << std::endl;
    ApXi->Rebin(ApXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ApXi==" << std::endl;
    ApXi->ReweightMixedEvent(ApXi->GetPairRebinned(iReb), 0.2, 0.9);


    std::cout << "==Rebinning XiAXi==" << std::endl;
    XiAXi->Rebin(XiAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting XiAXi==" << std::endl;
    XiAXi->ReweightMixedEvent(XiAXi->GetPairRebinned(iReb), 0.2, 0.9);

    std::cout << "==Rebinning LAXi==" << std::endl;
    LAXi->Rebin(LAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LAXi==" << std::endl;
    LAXi->ReweightMixedEvent(LAXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ALXi==" << std::endl;
    ALXi->Rebin(ALXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALXi==" << std::endl;
    ALXi->ReweightMixedEvent(ALXi->GetPairRebinned(iReb), 0.2, 0.9);
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

  std::cout << "==========CF pAp===============" << std::endl;
  CF_pAp->SetPairs(pAp, nullptr);
  CF_pAp->GetCorrelations();
  CF_pAp->WriteOutput(Form("%s/CFOutput_pAp_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF pAL===============" << std::endl;
  CF_pAL->SetPairs(pAL, ApL);
  CF_pAL->GetCorrelations();
  CF_pAL->WriteOutput(Form("%s/CFOutput_pAL_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF LAL===============" << std::endl;
  CF_LAL->SetPairs(LAL, nullptr);
  CF_LAL->GetCorrelations();
  CF_LAL->WriteOutput(Form("%s/CFOutput_LAL_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF pAXi===============" << std::endl;
  CF_pAXi->SetPairs(pAXi, ApXi);
  CF_pAXi->GetCorrelations();
  CF_pAXi->WriteOutput(Form("%s/CFOutput_pAXi_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF XiAXi===============" << std::endl;
  CF_XiAXi->SetPairs(XiAXi, nullptr);
  CF_XiAXi->GetCorrelations();
  CF_XiAXi->WriteOutput(Form("%s/CFOutput_XiAXi_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF LAXi===============" << std::endl;
  CF_LAXi->SetPairs(LAXi, ALXi);
  CF_LAXi->GetCorrelations();
  CF_LAXi->WriteOutput(Form("%s/CFOutput_LAXi_%s.root", gSystem->pwd(),addon));
}

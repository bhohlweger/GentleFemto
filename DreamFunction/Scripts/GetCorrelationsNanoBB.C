#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoBB(const char* filename,
                     const char* prefix, const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  Double_t norm1=0.24;//default 0.2-0.4, 0.18-0.28 where CFs look flatter
  Double_t norm2=0.34;// 0.24-0.34 from BB analysis

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", norm1, norm2);
  DreamPair* ApAp = new DreamPair("AntiPart", norm1, norm2);

  DreamCF* CF_pL = new DreamCF();
  DreamPair* pL = new DreamPair("Part",norm1,norm2);
  DreamPair* ApAL = new DreamPair("AntiPart",norm1,norm2);

  DreamCF* CF_LL = new DreamCF();
  DreamPair* LL = new DreamPair("Part", norm1, norm2);
  DreamPair* ALAL = new DreamPair("AntiPart", norm1, norm2);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part",norm1,norm2);
  DreamPair* ApAXi = new DreamPair("AntiPart",norm1,norm2);

  DreamCF* CF_XiXi = new DreamCF();
  DreamPair* XiXi = new DreamPair("Part", norm1, norm2);
  DreamPair* AXiAXi = new DreamPair("AntiPart", norm1, norm2);

  DreamCF* CF_LXi = new DreamCF();
  DreamPair* LXi = new DreamPair("Part",norm1,norm2);
  DreamPair* ALAXi = new DreamPair("AntiPart",norm1,norm2);

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

  XiXi->SetPair(DreamFile->GetPairDistributions(4, 4, ""));
  AXiAXi->SetPair(DreamFile->GetPairDistributions(5, 5, ""));

  LXi->SetPair(DreamFile->GetPairDistributions(2, 4, ""));
  ALAXi->SetPair(DreamFile->GetPairDistributions(3, 5, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==pp==" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  std::cout << "==pL==" << std::endl;
  pL->ShiftForEmpty(pL->GetPair());
  ApAL->ShiftForEmpty(ApAL->GetPair());

  std::cout << "==LL==" << std::endl;
  LL->ShiftForEmpty(LL->GetPair());
  std::cout << "==ALAL==" << std::endl;
  ALAL->ShiftForEmpty(ALAL->GetPair());

  std::cout << "==pXi==" << std::endl;
  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());

  std::cout << "==XiXi==" << std::endl;
  XiXi->ShiftForEmpty(XiXi->GetPair());
  AXiAXi->ShiftForEmpty(AXiAXi->GetPair());

  std::cout << "==LXi==" << std::endl;
  LXi->ShiftForEmpty(LXi->GetPair());
  ALAXi->ShiftForEmpty(ALAXi->GetPair());

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

  XiXi->FixShift(XiXi->GetPairShiftedEmpty(0), AXiAXi->GetPairShiftedEmpty(0),
		  	   AXiAXi->GetFirstBin());
  AXiAXi->FixShift(AXiAXi->GetPairShiftedEmpty(0), XiXi->GetPairShiftedEmpty(0),
		  	   XiXi->GetFirstBin());

  LXi->FixShift(LXi->GetPairShiftedEmpty(0), ALAXi->GetPairShiftedEmpty(0),
		  	   ALAXi->GetFirstBin());
  ALAXi->FixShift(ALAXi->GetPairShiftedEmpty(0), LXi->GetPairShiftedEmpty(0),
		  	   LXi->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
// To be applied only to p-antiL and L-antiL
  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning pL==" << std::endl;
    pL->Rebin(pL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting pL==" << std::endl;
    pL->ReweightMixedEvent(pL->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ApAL==" << std::endl;
    ApAL->Rebin(ApAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ApAL==" << std::endl;
    ApAL->ReweightMixedEvent(ApAL->GetPairRebinned(iReb), 0.2, 0.9);

    std::cout << "==Rebinning LL==" << std::endl;
    LL->Rebin(LL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LL==" << std::endl;
    LL->ReweightMixedEvent(LL->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ALAL==" << std::endl;
    ALAL->Rebin(ALAL->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALAL==" << std::endl;
    ALAL->ReweightMixedEvent(ALAL->GetPairRebinned(iReb), 0.2, 0.9);


    std::cout << "==Rebinning pXi==" << std::endl;
    pXi->Rebin(pXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting pXi==" << std::endl;
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ApXi==" << std::endl;
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ApAXi==" << std::endl;
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb), 0.2, 0.9);


    std::cout << "==Rebinning XiXi==" << std::endl;
    XiXi->Rebin(XiXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting XiXi==" << std::endl;
    XiXi->ReweightMixedEvent(XiXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning AXiAXi==" << std::endl;
    AXiAXi->Rebin(AXiAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting AXiAXi==" << std::endl;
    AXiAXi->ReweightMixedEvent(AXiAXi->GetPairRebinned(iReb), 0.2, 0.9);

    std::cout << "==Rebinning LXi==" << std::endl;
    LXi->Rebin(LXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting LXi==" << std::endl;
    LXi->ReweightMixedEvent(LXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning ALAXi==" << std::endl;
    ALAXi->Rebin(ALAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting ALAXi==" << std::endl;
    ALAXi->ReweightMixedEvent(ALAXi->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

//  pL->Rebin(pL->GetPair(), 4);
//  pL->Rebin(pL->GetPair(), 5);
//  ApAL->Rebin(ApAL->GetPair(), 4);
//  ApAL->Rebin(ApAL->GetPair(), 5);
//  LAL->Rebin(LAL->GetPair(), 4);
//  LAL->Rebin(LAL->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF pAp===============" << std::endl;
  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%s/CFOutput_pp_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF pL===============" << std::endl;
  CF_pL->SetPairs(pL, ApAL);
  CF_pL->GetCorrelations();
  CF_pL->WriteOutput(Form("%s/CFOutput_pL_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF LL===============" << std::endl;
  CF_LL->SetPairs(LL, ALAL);
  CF_LL->GetCorrelations();
  CF_LL->WriteOutput(Form("%s/CFOutput_LL_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF pXi===============" << std::endl;
  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  CF_pXi->WriteOutput(Form("%s/CFOutput_pXi_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF XiXi===============" << std::endl;
  CF_XiXi->SetPairs(XiXi, AXiAXi);
  CF_XiXi->GetCorrelations();
  CF_XiXi->WriteOutput(Form("%s/CFOutput_XiXi_%s.root", gSystem->pwd(),addon));

  std::cout << "==========CF LXi===============" << std::endl;
  CF_LXi->SetPairs(LXi, ALAXi);
  CF_LXi->GetCorrelations();
  CF_LXi->WriteOutput(Form("%s/CFOutput_LXi_%s.root", gSystem->pwd(),addon));
}

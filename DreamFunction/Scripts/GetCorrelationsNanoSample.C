#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoSample(const char* filename,
                     const char* prefix, const char* addon = "", double_t norm1 = 0.24, double_t norm2 = 0.34,
                     double_t norm3 = 0.18, double_t norm4 = 0.28) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFileSample(filename, prefix, addon);

  //default 0.2-0.4, 0.18-0.28 where CFs look flatter
  // 0.24-0.34 from BB analysis
//Baryons
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

//========================
//AntiBaryons
  DreamCF* CF_pAp = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart", norm3, norm4);

  DreamCF* CF_pAL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart",norm3, norm4);
  DreamPair* ApL = new DreamPair("AntiPartPart",norm3, norm4);

  DreamCF* CF_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm3, norm4);

  DreamCF* CF_pAXi = new DreamCF();
  DreamPair* pAXi = new DreamPair("PartAntiPart",norm3, norm4);
  DreamPair* ApXi = new DreamPair("AntiPartPart",norm3, norm4);

  DreamCF* CF_XiAXi = new DreamCF();
  DreamPair* XiAXi = new DreamPair("PartAntiPart",norm3, norm4);

  DreamCF* CF_LAXi = new DreamCF();
  DreamPair* LAXi = new DreamPair("PartAntiPart",norm3, norm4);
  DreamPair* ALXi = new DreamPair("AntiPartPart",norm3, norm4);

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

  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));

  pAXi->SetPair(DreamFile->GetPairDistributions(0, 5, ""));
  ApXi->SetPair(DreamFile->GetPairDistributions(1, 4, ""));
  XiAXi->SetPair(DreamFile->GetPairDistributions(4, 5, ""));
  LAXi->SetPair(DreamFile->GetPairDistributions(2, 5, ""));
  ALXi->SetPair(DreamFile->GetPairDistributions(3, 4, ""));

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

  std::cout << "==pAp==" << std::endl;
  pAp->ShiftForEmpty(pAp->GetPair());
  std::cout << "==pAL==" << std::endl;
  pAL->ShiftForEmpty(pAL->GetPair());
  ApL->ShiftForEmpty(ApL->GetPair());
  std::cout << "==LAL==" << std::endl;
  LAL->ShiftForEmpty(LAL->GetPair());

  std::cout << "==pAXi==" << std::endl;
  pAXi->ShiftForEmpty(pAXi->GetPair());
  ApXi->ShiftForEmpty(ApXi->GetPair());
  std::cout << "==XiAXi==" << std::endl;
  XiAXi->ShiftForEmpty(XiAXi->GetPair());
  std::cout << "==LAXi==" << std::endl;
  LAXi->ShiftForEmpty(LAXi->GetPair());
  ALXi->ShiftForEmpty(ALXi->GetPair());

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
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pAp->ReweightMixedEvent(pAp->GetPairFixShifted(0), 0.2, 0.9);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF pp===============" << std::endl;
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

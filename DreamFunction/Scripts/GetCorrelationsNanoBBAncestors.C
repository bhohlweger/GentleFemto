#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <iostream>

void GetCorrelationsNanoBBAncestors(const char* filename,
                     const char* prefix, const char* addon = "",
                     double_t norm1BB = 0.24, double_t norm2BB = 0.34,
                     double_t norm1BBar = 0.18, double_t norm2BBar = 0.28) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFileAncestors(filename, prefix, addon);

  std::cout << "norm1BB = " << norm1BB << std::endl;
  std::cout << "norm2BB = " << norm2BB << std::endl;

  // Double_t norm1=0.24;//default 0.2-0.4, 0.18-0.28 where CFs look flatter
  // Double_t norm2=0.34;// 0.24-0.34 from BB analysis

//Baryon - Baryon
  DreamCF* CF_ppCommon = new DreamCF();
  DreamPair* ppCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ApApCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

  DreamCF* CF_ppNonCommon = new DreamCF();
  DreamPair* ppNonCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ApApNonCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

  DreamCF* CF_pLCommon = new DreamCF();
  DreamPair* pLCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ApALCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

  DreamCF* CF_pLNonCommon = new DreamCF();
  DreamPair* pLNonCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ApALNonCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

  DreamCF* CF_LLCommon = new DreamCF();
  DreamPair* LLCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ALALCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

  DreamCF* CF_LLNonCommon = new DreamCF();
  DreamPair* LLNonCommon = new DreamPair("Part", norm1BB, norm2BB);
  DreamPair* ALALNonCommon = new DreamPair("AntiPart", norm1BB, norm2BB);

//Baryon - AntiBaryon
  DreamCF* CF_pApCommon = new DreamCF();
  DreamPair* pApCommon = new DreamPair("PartAntiPart", norm1BBar, norm2BBar);

  DreamCF* CF_pApNonCommon = new DreamCF();
  DreamPair* pApNonCommon = new DreamPair("PartAntiPart", norm1BBar, norm2BBar);

  DreamCF* CF_pALCommon = new DreamCF();
  DreamPair* pALCommon = new DreamPair("PartAntiPart",norm1BBar,norm2BBar);
  DreamPair* ApLCommon = new DreamPair("AntiPartPart",norm1BBar,norm2BBar);

  DreamCF* CF_pALNonCommon = new DreamCF();
  DreamPair* pALNonCommon = new DreamPair("PartAntiPart",norm1BBar,norm2BBar);
  DreamPair* ApLNonCommon = new DreamPair("AntiPartPart",norm1BBar,norm2BBar);

  DreamCF* CF_LALCommon = new DreamCF();
  DreamPair* LALCommon = new DreamPair("PartAntiPart", norm1BBar, norm2BBar);

  DreamCF* CF_LALNonCommon = new DreamCF();
  DreamPair* LALNonCommon = new DreamPair("PartAntiPart", norm1BBar, norm2BBar);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;

  ppCommon->SetPair(DreamFile->GetPairDistributionsCommon(0, 0, ""));
  ApApCommon->SetPair(DreamFile->GetPairDistributionsCommon(1, 1, ""));

  ppNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 0, ""));
  ApApNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(1, 1, ""));

  pLCommon->SetPair(DreamFile->GetPairDistributionsCommon(0, 2, ""));
  ApALCommon->SetPair(DreamFile->GetPairDistributionsCommon(1, 3, ""));

  pLNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 2, ""));
  ApALNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(1, 3, ""));

  LLCommon->SetPair(DreamFile->GetPairDistributionsCommon(2, 2, ""));
  ALALCommon->SetPair(DreamFile->GetPairDistributionsCommon(3, 3, ""));

  LLNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(2, 2, ""));
  ALALNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(3, 3, ""));

  pApCommon->SetPair(DreamFile->GetPairDistributionsCommon(0, 1, ""));
  pApNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 1, ""));

  pALCommon->SetPair(DreamFile->GetPairDistributionsCommon(0, 3, ""));
  ApLCommon->SetPair(DreamFile->GetPairDistributionsCommon(1, 2, ""));

  pALNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(0, 3, ""));
  ApLNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(1, 2, ""));

  LALCommon->SetPair(DreamFile->GetPairDistributionsCommon(2, 3, ""));
  LALNonCommon->SetPair(DreamFile->GetPairDistributionsNonCommon(2, 3, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;

  std::cout << "==pp==" << std::endl;
  ppCommon->ShiftForEmpty(ppCommon->GetPair());
  ApApCommon->ShiftForEmpty(ApApCommon->GetPair());
  ppNonCommon->ShiftForEmpty(ppNonCommon->GetPair());
  ApApNonCommon->ShiftForEmpty(ApApNonCommon->GetPair());

  std::cout << "==pL==" << std::endl;
  pLCommon->ShiftForEmpty(pLCommon->GetPair());
  ApALCommon->ShiftForEmpty(ApALCommon->GetPair());

  pLNonCommon->ShiftForEmpty(pLNonCommon->GetPair());
  ApALNonCommon->ShiftForEmpty(ApALNonCommon->GetPair());

  std::cout << "==LL==" << std::endl;
  LLCommon->ShiftForEmpty(LLCommon->GetPair());
  ALALCommon->ShiftForEmpty(ALALCommon->GetPair());

  LLNonCommon->ShiftForEmpty(LLNonCommon->GetPair());
  ALALNonCommon->ShiftForEmpty(ALALNonCommon->GetPair());

  pApCommon->ShiftForEmpty(pApCommon->GetPair());
  pApNonCommon->ShiftForEmpty(pApNonCommon->GetPair());

  pALCommon->ShiftForEmpty(pALCommon->GetPair());
  ApLCommon->ShiftForEmpty(ApLCommon->GetPair());
  pALNonCommon->ShiftForEmpty(pALNonCommon->GetPair());
  ApLNonCommon->ShiftForEmpty(ApLNonCommon->GetPair());

  LALCommon->ShiftForEmpty(LALCommon->GetPair());
  LALNonCommon->ShiftForEmpty(LALNonCommon->GetPair());


  ppCommon->ShiftForEmptyAncestors(ppCommon->GetPair());
  ApApCommon->ShiftForEmptyAncestors(ApApCommon->GetPair());

  ppCommon->ShiftForEmpty(ppCommon->GetPair());
  ApApCommon->ShiftForEmpty(ApApCommon->GetPair());

  ppNonCommon->ShiftForEmpty(ppNonCommon->GetPair());
  ApApNonCommon->ShiftForEmpty(ApApNonCommon->GetPair());

  std::cout << "==pL==" << std::endl;
  pLCommon->ShiftForEmpty(pLCommon->GetPair());
  ApALCommon->ShiftForEmpty(ApALCommon->GetPair());

  pLNonCommon->ShiftForEmpty(pLNonCommon->GetPair());
  ApALNonCommon->ShiftForEmpty(ApALNonCommon->GetPair());

  std::cout << "==LL==" << std::endl;
  LLCommon->ShiftForEmpty(LLCommon->GetPair());
  ALALCommon->ShiftForEmpty(ALALCommon->GetPair());

  LLNonCommon->ShiftForEmpty(LLNonCommon->GetPair());
  ALALNonCommon->ShiftForEmpty(ALALNonCommon->GetPair());

  pApCommon->ShiftForEmpty(pApCommon->GetPair());
  pApNonCommon->ShiftForEmpty(pApNonCommon->GetPair());

  pALCommon->ShiftForEmpty(pALCommon->GetPair());
  ApLCommon->ShiftForEmpty(ApLCommon->GetPair());
  pALNonCommon->ShiftForEmpty(pALNonCommon->GetPair());
  ApLNonCommon->ShiftForEmpty(ApLNonCommon->GetPair());

  LALCommon->ShiftForEmpty(LALCommon->GetPair());
  LALNonCommon->ShiftForEmpty(LALNonCommon->GetPair());


  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  ppCommon->FixShift(ppCommon->GetPairShiftedEmpty(0), ApApCommon->GetPairShiftedEmpty(0),
               ApApCommon->GetFirstBin());
  ApApCommon->FixShift(ApApCommon->GetPairShiftedEmpty(0), ppCommon->GetPairShiftedEmpty(0),
                 ppCommon->GetFirstBin());
  ppNonCommon->FixShift(ppNonCommon->GetPairShiftedEmpty(0), ApApNonCommon->GetPairShiftedEmpty(0),
               ApApNonCommon->GetFirstBin());
  ApApNonCommon->FixShift(ApApNonCommon->GetPairShiftedEmpty(0), ppNonCommon->GetPairShiftedEmpty(0),
                 ppNonCommon->GetFirstBin());

  pLCommon->FixShift(pLCommon->GetPairShiftedEmpty(0), ApALCommon->GetPairShiftedEmpty(0),
               ApALCommon->GetFirstBin());
  ApALCommon->FixShift(ApALCommon->GetPairShiftedEmpty(0), pLCommon->GetPairShiftedEmpty(0),
                 pLCommon->GetFirstBin());
  pLNonCommon->FixShift(pLNonCommon->GetPairShiftedEmpty(0), ApALNonCommon->GetPairShiftedEmpty(0),
               ApALNonCommon->GetFirstBin());
  ApALNonCommon->FixShift(ApALNonCommon->GetPairShiftedEmpty(0), pLNonCommon->GetPairShiftedEmpty(0),
                 pLNonCommon->GetFirstBin());

  LLCommon->FixShift(LLCommon->GetPairShiftedEmpty(0), ALALCommon->GetPairShiftedEmpty(0),
		  	   ALALCommon->GetFirstBin());
  ALALCommon->FixShift(ALALCommon->GetPairShiftedEmpty(0), LLCommon->GetPairShiftedEmpty(0),
		  	   LLCommon->GetFirstBin());

  LLNonCommon->FixShift(LLNonCommon->GetPairShiftedEmpty(0), ALALNonCommon->GetPairShiftedEmpty(0),
		  	   ALALNonCommon->GetFirstBin());
  ALALNonCommon->FixShift(ALALNonCommon->GetPairShiftedEmpty(0), LLNonCommon->GetPairShiftedEmpty(0),
		  	   LLNonCommon->GetFirstBin());

  pApCommon->FixShift(pApCommon->GetPairShiftedEmpty(0), pApCommon->GetPairShiftedEmpty(0),
               pApCommon->GetFirstBin());
  pApNonCommon->FixShift(pApNonCommon->GetPairShiftedEmpty(0), pApNonCommon->GetPairShiftedEmpty(0),
               pApNonCommon->GetFirstBin());

  pALCommon->FixShift(pALCommon->GetPairShiftedEmpty(0), ApLCommon->GetPairShiftedEmpty(0),
		  	   ApLCommon->GetFirstBin());
  ApLCommon->FixShift(ApLCommon->GetPairShiftedEmpty(0), pALCommon->GetPairShiftedEmpty(0),
		  	   pALCommon->GetFirstBin());

  pALNonCommon->FixShift(pALNonCommon->GetPairShiftedEmpty(0), ApLNonCommon->GetPairShiftedEmpty(0),
		  	   ApLNonCommon->GetFirstBin());
  ApLNonCommon->FixShift(ApLNonCommon->GetPairShiftedEmpty(0), pALNonCommon->GetPairShiftedEmpty(0),
		  	   pALNonCommon->GetFirstBin());

  LALCommon->FixShift(LALCommon->GetPairShiftedEmpty(0), LALCommon->GetPairShiftedEmpty(0),
		  	   LALCommon->GetFirstBin());
  LALNonCommon->FixShift(LALNonCommon->GetPairShiftedEmpty(0), LALNonCommon->GetPairShiftedEmpty(0),
		  	   LALNonCommon->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
// To be applied only to p-antiL and L-antiL

  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    pLCommon->Rebin(pLCommon->GetPairFixShifted(0), rebinVec[iReb]);
  std::cout << "===========Reweigthing pL common==============" << std::endl;
    pLCommon->ReweightMixedEvent(pLCommon->GetPairRebinned(iReb), 0.2, 0.9);

    pLNonCommon->Rebin(pLNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    pLNonCommon->ReweightMixedEvent(pLNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ApALCommon->Rebin(ApALCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ApALCommon->ReweightMixedEvent(ApALCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ApALNonCommon->Rebin(ApALNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ApALNonCommon->ReweightMixedEvent(ApALNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    LLCommon->Rebin(LLCommon->GetPairFixShifted(0), rebinVec[iReb]);
    LLCommon->ReweightMixedEvent(LLCommon->GetPairRebinned(iReb), 0.2, 0.9);

    LLNonCommon->Rebin(LLNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    LLNonCommon->ReweightMixedEvent(LLNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ALALCommon->Rebin(ALALCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ALALCommon->ReweightMixedEvent(ALALCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ALALNonCommon->Rebin(ALALNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ALALNonCommon->ReweightMixedEvent(ALALNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    pALCommon->Rebin(pALCommon->GetPairFixShifted(0), rebinVec[iReb]);
    pALCommon->ReweightMixedEvent(pALCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ApLCommon->Rebin(ApLCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ApLCommon->ReweightMixedEvent(ApLCommon->GetPairRebinned(iReb), 0.2, 0.9);

    pALNonCommon->Rebin(pALNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    pALNonCommon->ReweightMixedEvent(pALNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    ApLNonCommon->Rebin(ApLNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    ApLNonCommon->ReweightMixedEvent(ApLNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

    LALCommon->Rebin(LALCommon->GetPairFixShifted(0), rebinVec[iReb]);
    LALCommon->ReweightMixedEvent(LALCommon->GetPairRebinned(iReb), 0.2, 0.9);

    LALNonCommon->Rebin(LALNonCommon->GetPairFixShifted(0), rebinVec[iReb]);
    LALNonCommon->ReweightMixedEvent(LALNonCommon->GetPairRebinned(iReb), 0.2, 0.9);

  }

  ppCommon->ReweightMixedEvent(ppCommon->GetPairFixShifted(0), 0.2, 0.9);
  ApApCommon->ReweightMixedEvent(ApApCommon->GetPairFixShifted(0), 0.2, 0.9);
  pApCommon->ReweightMixedEvent(pApCommon->GetPairFixShifted(0), 0.2, 0.9);
  ppNonCommon->ReweightMixedEvent(ppNonCommon->GetPairFixShifted(0), 0.2, 0.9);
  ApApNonCommon->ReweightMixedEvent(ApApNonCommon->GetPairFixShifted(0), 0.2, 0.9);
  pApNonCommon->ReweightMixedEvent(pApNonCommon->GetPairFixShifted(0), 0.2, 0.9);

  // pLCommon->Rebin(pLCommon->GetPair(), 4);
  // pLCommon->Rebin(pLCommon->GetPair(), 5);
  // ApALCommon->Rebin(ApALCommon->GetPair(), 4);
  // ApALCommon->Rebin(ApALCommon->GetPair(), 5);

  // pLNonCommon->Rebin(pLNonCommon->GetPair(), 4);
  // pLNonCommon->Rebin(pLNonCommon->GetPair(), 5);
  // ApALNonCommon->Rebin(ApALNonCommon->GetPair(), 4);
  // ApALNonCommon->Rebin(ApALNonCommon->GetPair(), 5);


  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  std::cout << "==========CF pAp===============" << std::endl;
  CF_ppCommon->SetPairs(ppCommon, ApApCommon);
  CF_ppCommon->GetCorrelations();
  CF_ppCommon->WriteOutput(Form("%s/CFOutput_pp_%s_Common.root", gSystem->pwd(),addon));

  CF_ppNonCommon->SetPairs(ppNonCommon, ppNonCommon);
  CF_ppNonCommon->GetCorrelations();
  CF_ppNonCommon->WriteOutput(Form("%s/CFOutput_pp_%s_NonCommon.root", gSystem->pwd(),addon));

  CF_pLCommon->SetPairs(pLCommon, ApALCommon);
  CF_pLCommon->GetCorrelations();
  CF_pLCommon->WriteOutput(Form("%s/CFOutput_pL_%s_Common.root", gSystem->pwd(),addon));

  CF_pLNonCommon->SetPairs(pLNonCommon, ApALNonCommon);
  CF_pLNonCommon->GetCorrelations();
  CF_pLNonCommon->WriteOutput(Form("%s/CFOutput_pL_%s_NonCommon.root", gSystem->pwd(),addon));

  CF_LLCommon->SetPairs(LLCommon, ALALCommon);
  CF_LLCommon->GetCorrelations();
  CF_LLCommon->WriteOutput(Form("%s/CFOutput_LL_%s_Common.root", gSystem->pwd(),addon));

  CF_LLNonCommon->SetPairs(LLNonCommon, ALALNonCommon);
  CF_LLNonCommon->GetCorrelations();
  CF_LLNonCommon->WriteOutput(Form("%s/CFOutput_LL_%s_NonCommon.root", gSystem->pwd(),addon));

  CF_pApCommon->SetPairs(pApCommon, nullptr);
  CF_pApCommon->GetCorrelations();
  CF_pApCommon->WriteOutput(Form("%s/CFOutput_pAp_%s_Common.root", gSystem->pwd(),addon));

  CF_pApNonCommon->SetPairs(pApNonCommon, nullptr);
  CF_pApNonCommon->GetCorrelations();
  CF_pApNonCommon->WriteOutput(Form("%s/CFOutput_pAp_%s_NonCommon.root", gSystem->pwd(),addon));

  CF_pALCommon->SetPairs(pALCommon, ApLCommon);
  CF_pALCommon->GetCorrelations();
  CF_pALCommon->WriteOutput(Form("%s/CFOutput_pAL_%s_Common.root", gSystem->pwd(),addon));

  CF_pALNonCommon->SetPairs(pALNonCommon, ApLNonCommon);
  CF_pALNonCommon->GetCorrelations();
  CF_pALNonCommon->WriteOutput(Form("%s/CFOutput_pAL_%s_NonCommon.root", gSystem->pwd(),addon));

  CF_LALCommon->SetPairs(LALCommon, nullptr);
  CF_LALCommon->GetCorrelations();
  CF_LALCommon->WriteOutput(Form("%s/CFOutput_LAL_%s_Common.root", gSystem->pwd(),addon));

  CF_LALNonCommon->SetPairs(LALNonCommon, nullptr);
  CF_LALNonCommon->GetCorrelations();
  CF_LALNonCommon->WriteOutput(Form("%s/CFOutput_LAL_%s_NonCommon.root", gSystem->pwd(),addon));


}

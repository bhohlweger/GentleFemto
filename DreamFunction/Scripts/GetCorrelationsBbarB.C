#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void GetCorrelationsBbarB(const char* filename, const char* prefix,
                     const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pAp_App = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart", 0.2, 0.4);
//  DreamPair* App = new DreamPair("AntiPartPart", 0.2, 0.4);

  DreamCF* CF_pAL_ApL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart", 0.2, 0.4);
  DreamPair* ApL = new DreamPair("AntiPartPart", 0.2, 0.4);

  DreamCF* CF_ALL_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", 0.2, 0.4);
//  DreamPair* ALL = new DreamPair("AntiPartPart", 0.2, 0.4);

  DreamCF* CF_pAXi_ApXi = new DreamCF();
  DreamPair* pAXi = new DreamPair("PartAntiPart", 0.2, 0.4);
  DreamPair* ApXi = new DreamPair("AntiPartPart", 0.2, 0.4);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
//  App->SetPair(DreamFile->GetPairDistributions(1, 0, ""));

  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));
//  ALL->SetPair(DreamFile->GetPairDistributions(3, 2, ""));

  pAXi->SetPair(DreamFile->GetPairDistributions(0, 5, ""));
  ApXi->SetPair(DreamFile->GetPairDistributions(1, 4, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->ShiftForEmpty(pAp->GetPair());
//  App->ShiftForEmpty(App->GetPair());

  pAL->ShiftForEmpty(pAL->GetPair());
  ApL->ShiftForEmpty(ApL->GetPair());

  LAL->ShiftForEmpty(LAL->GetPair());
//  ALL->ShiftForEmpty(ALL->GetPair());

  pAXi->ShiftForEmpty(pAXi->GetPair());
  ApXi->ShiftForEmpty(ApXi->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
//  pAp->FixShift(pAp->GetPairShiftedEmpty(0), App->GetPairShiftedEmpty(0),
//                App->GetFirstBin());
//  App->FixShift(App->GetPairShiftedEmpty(0), pAp->GetPairShiftedEmpty(0),
//                 pAp->GetFirstBin());

  pAL->FixShift(pAL->GetPairShiftedEmpty(0), ApL->GetPairShiftedEmpty(0),
               ApL->GetFirstBin());
  ApL->FixShift(ApL->GetPairShiftedEmpty(0), pAL->GetPairShiftedEmpty(0),
                 pAL->GetFirstBin());
//  LAL->FixShift(LAL->GetPairShiftedEmpty(0), ALL->GetPairShiftedEmpty(0),
//               ALL->GetFirstBin());
//  ALL->FixShift(ALL->GetPairShiftedEmpty(0), LAL->GetPairShiftedEmpty(0),
//                 LAL->GetFirstBin());

  pAXi->FixShift(pAXi->GetPairShiftedEmpty(0), ApXi->GetPairShiftedEmpty(0),
                ApXi->GetFirstBin());
  ApXi->FixShift(ApXi->GetPairShiftedEmpty(0), pAXi->GetPairShiftedEmpty(0),
                  pAXi->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pAL->Rebin(pAL->GetPairFixShifted(0), iReb);
    pAXi->Rebin(pAXi->GetPairFixShifted(0), iReb);
    //for LAL there is no fix shift to be had!
    LAL->Rebin(LAL->GetPairShiftedEmpty(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pAL->ReweightMixedEvent(pAL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    LAL->ReweightMixedEvent(LAL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pAXi->ReweightMixedEvent(pAXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApL->Rebin(ApL->GetPairFixShifted(0), iReb);
    ApXi->Rebin(ApXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    ApL->ReweightMixedEvent(ApL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    ApXi->ReweightMixedEvent(ApXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
  }
  pAp->ReweightMixedEvent(pAp->GetPairShiftedEmpty(0), 0.2, 0.9);

  pAL->Rebin(pAL->GetPair(), 4);
  pAL->Rebin(pAL->GetPair(), 5);
  ApL->Rebin(ApL->GetPair(), 4);
  ApL->Rebin(ApL->GetPair(), 5);
  LAL->Rebin(LAL->GetPair(), 4);
  LAL->Rebin(LAL->GetPair(), 5);
  pAXi->Rebin(pAXi->GetPair(), 4);
  pAXi->Rebin(pAXi->GetPair(), 5);
  ApXi->Rebin(ApXi->GetPair(), 4);
  ApXi->Rebin(ApXi->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  std::cout << "$PWD " << foldername.Data() << std::endl;
  std::cout << "pp CF \n";
  std::cout << "Set Pair \n";
  CF_pAp_App->SetPairs(pAp, nullptr);
  std::cout << "Get CF \n";
  CF_pAp_App->GetCorrelations();
  std::cout << "Write Output \n";
  CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App.root", foldername.Data()));

  std::cout << "pL CF \n";
  CF_pAL_ApL->SetPairs(pAL, ApL);
  CF_pAL_ApL->GetCorrelations();
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL.root", foldername.Data()));

  std::cout << "LL CF \n";
  CF_ALL_LAL->SetPairs(LAL, nullptr);
  CF_ALL_LAL->GetCorrelations();
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL.root", foldername.Data()));

  std::cout << "pXi CF \n";
  CF_pAXi_ApXi->SetPairs(pAXi, ApXi);
  CF_pAXi_ApXi->GetCorrelations();
  CF_pAXi_ApXi->WriteOutput(Form("%sCFOutput_pAXi_ApXi.root", foldername.Data()));
}

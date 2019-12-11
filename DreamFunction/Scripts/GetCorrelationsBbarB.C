#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void GetCorrelationsBbarB(const char* filename, const char* prefix,
                     const char* addon = "", const bool isMC = false,
                     const bool epos = false) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  TString add1="1";
  TString add2="2";
  TString add3="3";
  TString add4="8";
  TString add5="5";

  Double_t norm1=0.2;
  Double_t norm2=0.4;

  DreamCF* CF_pAp_App = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart", norm1, norm2);
//  DreamPair* App = new DreamPair("AntiPartPart", 0.2, 0.4);

  DreamCF* CF_pAL_ApL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair* ApL = new DreamPair("AntiPartPart", norm1, norm2);

  DreamCF* CF_ALL_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);
//  DreamPair* ALL = new DreamPair("AntiPartPart", 0.2, 0.4);

  DreamCF* CF_pAXi_ApXi = new DreamCF();
  DreamPair* pAXi = new DreamPair("PartAntiPart", norm1, norm2);
  DreamPair* ApXi = new DreamPair("AntiPartPart", norm1, norm2);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
//  App->SetPair(DreamFile->GetPairDistributions(1, 0, ""));

  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));
//  ALL->SetPair(DreamFile->GetPairDistributions(3, 2, ""));

  pAXi->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  ApXi->SetPair(DreamFile->GetPairDistributions(2, 2, ""));

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
  //Fix shift singe pair cfs anyway to ensure compatibility!
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

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pAL->Rebin(pAL->GetPairFixShifted(0), iReb);
    pAXi->Rebin(pAXi->GetPairFixShifted(0), iReb);
    LAL->Rebin(LAL->GetPairFixShifted(0), iReb);
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

if(!isMC)
{
  TString foldername = gSystem->pwd();

  std::cout << "$PWD " << foldername << std::endl;
  std::cout << "pp CF \n";
  std::cout << "Set Pair \n";
  CF_pAp_App->SetPairs(pAp, nullptr);
  std::cout << "Get CF \n";
  CF_pAp_App->GetCorrelations();
  std::cout << "Write Output \n";
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st1.root", gSystem->pwd()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st2.root", gSystem->pwd()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st3.root", gSystem->pwd()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_full.root", gSystem->pwd()));
  }
  else if(strcmp(addon, add5)==0){
    std::cout << "Sphericity Cuts [0.9-1]\n";
    CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st5.root", gSystem->pwd()));
    }
;
  std::cout << "pL CF \n";
  CF_pAL_ApL->SetPairs(pAL, ApL);
  CF_pAL_ApL->GetCorrelations();
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st1.root", gSystem->pwd()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st2.root", gSystem->pwd()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st3.root", gSystem->pwd()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_full.root", gSystem->pwd()));
}else if(strcmp(addon, add5)==0){
  std::cout << "Sphericity Cuts [0.9-1]\n";
  CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st5.root", gSystem->pwd()));
}
;

  std::cout << "LL CF \n";
  CF_ALL_LAL->SetPairs(LAL, nullptr);
  CF_ALL_LAL->GetCorrelations();
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st1.root", gSystem->pwd()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st2.root", gSystem->pwd()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st3.root", gSystem->pwd()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_full.root", gSystem->pwd()));
}else if(strcmp(addon, add5)==0){
  std::cout << "Sphericity Cuts [0.9-1]\n";
  CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st5.root", gSystem->pwd()));
}
;
} else if (isMC)
{
  TString foldername;
  if(!epos)
  {
  foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/mc/Pythia/";
  }else if(epos){
  foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/mc/epos/";
  }
  std::cout << "$PWD " << foldername << std::endl;
  std::cout << "pp CF \n";
  std::cout << "Set Pair \n";
  CF_pAp_App->SetPairs(pAp, nullptr);
  std::cout << "Get CF \n";
  CF_pAp_App->GetCorrelations();
  std::cout << "Write Output \n";
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st1.root", foldername.Data()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st2.root", foldername.Data()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st3.root", foldername.Data()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_full.root", foldername.Data()));
  }
  else if(strcmp(addon, add5)==0){
    std::cout << "Sphericity Cuts [0.9-1]\n";
    CF_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st5.root", foldername.Data()));
    }
;
  std::cout << "pL CF \n";
  CF_pAL_ApL->SetPairs(pAL, ApL);
  CF_pAL_ApL->GetCorrelations();
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st1.root", foldername.Data()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st2.root", foldername.Data()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st3.root", foldername.Data()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_full.root", foldername.Data()));
}else if(strcmp(addon, add5)==0){
  std::cout << "Sphericity Cuts [0.9-1]\n";
  CF_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st5.root", foldername.Data()));
}
;

  std::cout << "LL CF \n";
  CF_ALL_LAL->SetPairs(LAL, nullptr);
  CF_ALL_LAL->GetCorrelations();
  if(strcmp(addon, add1)==0){
  std::cout << "Sphericity [0.,0.3] \n";
  CF_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st1.root", foldername.Data()));
}else if(strcmp(addon, add2)==0){
  std::cout << "Sphericity [0.3,0.7] \n";
  CF_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st2.root", foldername.Data()));
}else if(strcmp(addon, add3)==0){
  std::cout << "Sphericity [0.7,1.] \n";
  CF_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st3.root", foldername.Data()));
}else if(strcmp(addon, add4)==0){
  std::cout << "No Sphericity Cuts \n";
  CF_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_full.root", foldername.Data()));
}else if(strcmp(addon, add5)==0){
  std::cout << "Sphericity Cuts [0.9-1]\n";
  CF_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st5.root", foldername.Data()));
}
};


}

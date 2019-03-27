#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void GetCorrectBinCF(const char* filename, const char* filenameMC, const char* prefix,
                     const char* addon = "", const bool epos = false) {

  ReadDreamFile* DreamFiledata = new ReadDreamFile(4, 4);
  ReadDreamFile* DreamFilemc = new ReadDreamFile(4, 4);

  DreamFiledata->SetAnalysisFile(filename, prefix, addon);
  DreamFilemc->SetAnalysisFile(filename, prefix, addon);

  TString add1="1";
  TString add2="2";
  TString add3="3";
  TString add4="4";
  TString add5="5";

  Double_t norm1=0.2;
  Double_t norm2=0.4;

  TFile* _file0data=TFile::Open(filename);
  TFile* _file0mc=TFile::Open(filenameMC);

  TDirectoryFile *dirResultsdata=(TDirectoryFile*)(_file0data->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TDirectoryFile *dirResultsmc=(TDirectoryFile*)(_file0mc->FindObjectAny(Form("%sResults%s", prefix, addon)));

  TList *Resultsdata;
  TList *Resultsmc;
  dirResultsdata->GetObject(Form("%sResults%s", prefix, addon),Resultsdata);
  dirResultsmc->GetObject(Form("%sResults%s", prefix, addon),Resultsmc);
//Proton-AntiProton----------------------------------------------------
  DreamDist* pAp_pairdata = new DreamDist();
  TList* tmpFolderdata=(TList*)Resultsdata->FindObject("Particle0_Particle1");
  TH1F* SEDatapAp_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle0_Particle1");
  TH1F* MEDatapAp_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle0_Particle1");
  TH2F* SEDataMultpAp_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle0_Particle1");
  TH2F* MEDataMultpAp_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle0_Particle1");
   pAp_pairdata->SetSEDist(SEDatapAp_pair , "");
   pAp_pairdata->SetSEMultDist(SEDataMultpAp_pair, "");
   pAp_pairdata->SetMEDist(MEDatapAp_pair, "");
   pAp_pairdata->SetMEMultDist(MEDataMultpAp_pair, "");

  DreamDist* pAp_pairmc = new DreamDist();
  TList* tmpFoldermc=(TList*)Resultsmc->FindObject("Particle0_Particle1");
  TH1F* SEMCpAp_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle0_Particle1");
  TH1F* MEMCpAp_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle0_Particle1");
  TH2F* SEMCMultpAp_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle0_Particle1");
  TH2F* MEMCMultpAp_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle0_Particle1");
   pAp_pairmc->SetSEDist(SEMCpAp_pair , "");
   pAp_pairmc->SetSEMultDist(SEMCMultpAp_pair, "");
   pAp_pairmc->SetMEDist(MEMCpAp_pair, "");
   pAp_pairmc->SetMEMultDist(MEMCMultpAp_pair, "");

   //Lambda-AntiLambda----------------------------------------------------
   DreamDist* LAL_pairdata = new DreamDist();
   tmpFolderdata=(TList*)Resultsdata->FindObject("Particle2_Particle3");
   TH1F* SEDataLAL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle2_Particle3");
   TH1F* MEDataLAL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle2_Particle3");
   TH2F* SEDataMultLAL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle2_Particle3");
   TH2F* MEDataMultLAL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle2_Particle3");
    LAL_pairdata->SetSEDist(SEDataLAL_pair , "_Shifted");
    LAL_pairdata->SetSEMultDist(SEDataMultLAL_pair, "_Shifted");
    LAL_pairdata->SetMEDist(MEDataLAL_pair, "_Shifted");
    LAL_pairdata->SetMEMultDist(MEDataMultLAL_pair, "_Shifted");

   DreamDist* LAL_pairmc = new DreamDist();
   tmpFoldermc=(TList*)Resultsmc->FindObject("Particle2_Particle3");
   TH1F* SEMCLAL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle2_Particle3");
   TH1F* MEMCLAL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle2_Particle3");
   TH2F* SEMCMultLAL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle2_Particle3");
   TH2F* MEMCMultLAL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle2_Particle3");
    LAL_pairmc->SetSEDist(SEMCLAL_pair , "_Shifted");
    LAL_pairmc->SetSEMultDist(SEMCMultLAL_pair, "_Shifted");
    LAL_pairmc->SetMEDist(MEMCLAL_pair, "_Shifted");
    LAL_pairmc->SetMEMultDist(MEMCMultLAL_pair, "_Shifted");

    //Proton-AntiLambda----------------------------------------------------
    DreamDist* pAL_pairdata = new DreamDist();
    tmpFolderdata=(TList*)Resultsdata->FindObject("Particle0_Particle3");
    TH1F* SEDatapAL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle0_Particle3");
    TH1F* MEDatapAL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle0_Particle3");
    TH2F* SEDataMultpAL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle0_Particle3");
    TH2F* MEDataMultpAL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle0_Particle3");
     pAL_pairdata->SetSEDist(SEDatapAL_pair , "_Shifted");
     pAL_pairdata->SetSEMultDist(SEDataMultpAL_pair, "_Shifted");
     pAL_pairdata->SetMEDist(MEDatapAL_pair, "_Shifted");
     pAL_pairdata->SetMEMultDist(MEDataMultpAL_pair, "_Shifted");

     DreamDist* pAL_pairmc = new DreamDist();
     tmpFoldermc=(TList*)Resultsmc->FindObject("Particle0_Particle3");
     TH1F* SEMCpAL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle0_Particle3");
     TH1F* MEMCpAL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle0_Particle3");
     TH2F* SEMCMultpAL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle0_Particle3");
     TH2F* MEMCMultpAL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle0_Particle3");
      pAL_pairmc->SetSEDist(SEMCpAL_pair , "_Shifted");
      pAL_pairmc->SetSEMultDist(SEMCMultpAL_pair, "_Shifted");
      pAL_pairmc->SetMEDist(MEMCpAL_pair, "_Shifted");
      pAL_pairmc->SetMEMultDist(MEMCMultpAL_pair, "_Shifted");

     //AntiProton-Lambda----------------------------------------------------
     DreamDist* ApL_pairdata = new DreamDist();
     tmpFolderdata=(TList*)Resultsdata->FindObject("Particle1_Particle2");
     TH1F* SEDataApL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle1_Particle2");
     TH1F* MEDataApL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle1_Particle2");
     TH2F* SEDataMultApL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle1_Particle2");
     TH2F* MEDataMultApL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle1_Particle2");
      ApL_pairdata->SetSEDist(SEDataApL_pair , "_Shifted");
      ApL_pairdata->SetSEMultDist(SEDataMultApL_pair, "_Shifted");
      ApL_pairdata->SetMEDist(MEDataApL_pair, "_Shifted");
      ApL_pairdata->SetMEMultDist(MEDataMultApL_pair, "_Shifted");

      DreamDist* ApL_pairmc = new DreamDist();
      tmpFoldermc=(TList*)Resultsmc->FindObject("Particle1_Particle2");
      TH1F* SEMCApL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle1_Particle2");
      TH1F* MEMCApL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle1_Particle2");
      TH2F* SEMCMultApL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle1_Particle2");
      TH2F* MEMCMultApL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle1_Particle2");
       ApL_pairmc->SetSEDist(SEMCApL_pair , "_Shifted");
       ApL_pairmc->SetSEMultDist(SEMCMultApL_pair, "_Shifted");
       ApL_pairmc->SetMEDist(MEMCApL_pair, "_Shifted");
       ApL_pairmc->SetMEMultDist(MEMCMultApL_pair, "_Shifted");

   DreamCF* CFData_pAp_App = new DreamCF();
   DreamCF* CFData_ALL_LAL = new DreamCF();
   DreamCF* CFData_pAL_ApL = new DreamCF();


   DreamCF* CFMC_pAp_App = new DreamCF();
   DreamCF* CFMC_ALL_LAL = new DreamCF();
   DreamCF* CFMC_pAL_ApL = new DreamCF();

   //   DreamPair* pAp = new DreamPair("PartAntiPart", norm1, norm2);
   //
   //   DreamCF* CF_pAL_ApL = new DreamCF();
   //   DreamPair* pAL = new DreamPair("PartAntiPart", norm1, norm2);
   //   DreamPair* ApL = new DreamPair("AntiPartPart", norm1, norm2);
   //
   //   DreamCF* CF_ALL_LAL = new DreamCF();
   //   DreamPair* LAL = new DreamPair("PartAntiPart", norm1, norm2);

   DreamPair* pApdata = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* LALdata = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* pALdata = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* ApLdata = new DreamPair("AntiPartPart", norm1, norm2);

   DreamPair* pApmc = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* LALmc = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* pALmc = new DreamPair("PartAntiPart", norm1, norm2);
   DreamPair* ApLmc = new DreamPair("AntiPartPart", norm1, norm2);

     std::cout << "=========================" << std::endl;
     std::cout << "========Pair Set=========" << std::endl;
     std::cout << "=========================" << std::endl;
   pApdata->SetPair(pAp_pairdata);
   pApmc->SetPair(pAp_pairmc);
   LALdata->SetPair(LAL_pairdata);
   LALmc->SetPair(LAL_pairmc);
   pALdata->SetPair(pAL_pairdata);
   pALmc->SetPair(pAL_pairmc);
   ApLdata->SetPair(ApL_pairdata);
   ApLmc->SetPair(ApL_pairmc);

     std::cout << "=========================" << std::endl;
     std::cout << "======Pair Shifted=======" << std::endl;
     std::cout << "=========================" << std::endl;
     pApdata->ShiftForEmpty(pApdata->GetPair());
     LALdata->ShiftForEmpty(LALdata->GetPair());
     pALdata->ShiftForEmpty(pALdata->GetPair());
     ApLdata->ShiftForEmpty(ApLdata->GetPair());

     pApmc->ShiftForEmpty(pApmc->GetPair());
     LALmc->ShiftForEmpty(LALmc->GetPair());
     pALmc->ShiftForEmpty(pALmc->GetPair());
     ApLmc->ShiftForEmpty(ApLmc->GetPair());

     // std::cout<< "Smallest k* for pApdata = "<< pApdata->GetFirstBin()<<std::endl;
     // std::cout<< "Smallest k* for pApmc = "<< pApmc->GetFirstBin()<<std::endl;
     //
     // std::cout<< "Smallest k* for LALdata = "<< LALdata->GetFirstBin()<<std::endl;
     // std::cout<< "Smallest k* for LALmc = "<< LALmc->GetFirstBin()<<std::endl;


      double spApdata = pApdata->GetFirstBin();
      double spApmc = pApmc->GetFirstBin();
      double spApmin;
      double sLALdata = LALdata->GetFirstBin();
      double sLALmc = LALmc->GetFirstBin();
      double sLALmin;

      std::cout<<"Min bin ppdata = "<<spApdata<<std::endl;
      std::cout<<"Min bin ppMC = "<<spApmc<<std::endl;

      double spALdata = pALdata->GetFirstBin();
      double spALmc = pALmc->GetFirstBin();
      double spALmin;

      double sApLdata = ApLdata->GetFirstBin();
      double sApLmc = ApLmc->GetFirstBin();
      double sApLmin;

      double spALApLmin;

      if(spApdata<=spApmc){
        spApmin = spApdata;
      }
      else if(spApmc<=spApdata){
        spApmin = spApmc;
      }
      if(sLALdata<=sLALmc){
        sLALmin = sLALdata;
      }
      else if(sLALmc<=sLALdata){
        sLALmin = sLALmc;
      }

      if(spALdata<=spALmc){
        spALmin = spALdata;
      }
      else if(spALmc<=spALdata){
        spALmin = spALmc;
      }
      if(sApLdata<=sApLmc){
        sApLmin = sApLdata;
      }
      else if(sApLmc<=sApLdata){
        sApLmin = sApLmc;
      }
      if(spALmin<=sApLmin){
      spALApLmin =  spALmin;
      }
      else if(sApLmin<=spALmin){
      spALApLmin = sApLmin;
      }
      std::cout<< "Smallest k* for pAp = "<< spApmin <<std::endl;

      std::cout<< "Smallest k* for LAL = "<< sLALmin<<std::endl;


std::cout << "=========================" << std::endl;
std::cout << "====Pair Fix Shifted=====" << std::endl;
std::cout << "=========================" << std::endl;
//Fix shift singe pair cfs anyway to ensure compatibility!
pApdata->FixShift(pAp_pairdata,pAp_pairdata,
                  spApmin,true);
pApmc->FixShift(pAp_pairmc, pAp_pairdata,
                  spApmin,true);
std::cout<< "FIXSHIFTHING: Smallest k* for pApdata = "<< pApdata->GetFirstBin()<<std::endl;
std::cout<< "FIXSHIFTHING: Smallest k* for pApmc = "<< pApmc->GetFirstBin()<<std::endl;

LALdata->FixShift(LAL_pairdata,LAL_pairdata,
                  sLALmin,true);
LALmc->FixShift(LAL_pairmc, LAL_pairdata,
                  sLALmin,true);
std::cout<< "FIXSHIFTHING: Smallest k* for LALdata = "<< LALdata->GetFirstBin()<<std::endl;
std::cout<< "FIXSHIFTHING: Smallest k* for LALmc = "<< LALmc->GetFirstBin()<<std::endl;

pALdata->FixShift(pAL_pairdata,pAL_pairdata,
                  spALApLmin,true);
pALmc->FixShift(pAL_pairmc, pAL_pairdata,
                  spALApLmin,true);
ApLdata->FixShift(ApL_pairdata,ApL_pairdata,
                  spALApLmin,true);
ApLmc->FixShift(ApL_pairmc, ApL_pairdata,
                  spALApLmin,true);

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
     pALdata->Rebin(pALdata->GetPairFixShifted(0), iReb);
     pALmc->Rebin(pALmc->GetPairFixShifted(0), iReb);
    // pAXi->Rebin(pAXi->GetPairFixShifted(0), iReb);
    LALdata->Rebin(LALdata->GetPairFixShifted(0), iReb);
    LALmc->Rebin(LALmc->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pALdata->ReweightMixedEvent(pALdata->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pALmc->ReweightMixedEvent(pALmc->GetPairRebinned(iReb - 4), 0.2, 0.9);
    LALdata->ReweightMixedEvent(LALdata->GetPairRebinned(iReb - 4), 0.2, 0.9);
    LALmc->ReweightMixedEvent(LALmc->GetPairRebinned(iReb - 4), 0.2, 0.9);
    // pAXi->ReweightMixedEvent(pAXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
     ApLdata->Rebin(ApLdata->GetPairFixShifted(0), iReb);
     ApLmc->Rebin(ApLmc->GetPairFixShifted(0), iReb);
    // ApXi->Rebin(ApXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
     ApLdata->ReweightMixedEvent(ApLdata->GetPairRebinned(iReb - 4), 0.2, 0.9);
     ApLmc->ReweightMixedEvent(ApLmc->GetPairRebinned(iReb - 4), 0.2, 0.9);
    // ApXi->ReweightMixedEvent(ApXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
  }
  pApdata->ReweightMixedEvent(pApdata->GetPairShiftedEmpty(0), 0.2, 0.9);
  pApmc->ReweightMixedEvent(pApmc->GetPairShiftedEmpty(0), 0.2, 0.9);
//
   pALdata->Rebin(pALdata->GetPair(), 4);
   pALdata->Rebin(pALdata->GetPair(), 5);
   pALmc->Rebin(pALmc->GetPair(), 4);
   pALmc->Rebin(pALmc->GetPair(), 5);
   ApLdata->Rebin(ApLdata->GetPair(), 4);
   ApLdata->Rebin(ApLdata->GetPair(), 5);
   ApLmc->Rebin(ApLmc->GetPair(), 4);
   ApLmc->Rebin(ApLmc->GetPair(), 5);
  LALdata->Rebin(LALdata->GetPair(), 4);
  LALmc->Rebin(LALmc->GetPair(), 5);

//
  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;


    TString foldernamedata = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/data/CorrectBin/";
    TString foldernamemc;
    if(!epos){
    foldernamemc = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/";
  }else if(epos){
    foldernamemc = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/mc/CorrectBin/epos/";
  }
    std::cout << "$PWD " << foldernamedata << std::endl;
    std::cout << "pp CF DATA\n";
    std::cout << "Set Pair \n";
    CFData_pAp_App->SetPairs(pApdata, nullptr);
    std::cout << "Get CF \n";
    CFData_pAp_App->GetCorrelations("pApData");
    std::cout << "Write Output \n";

    std::cout << "$PWD " << foldernamemc << std::endl;
    std::cout << "pp CF MC\n";
    std::cout << "Set Pair \n";
    CFMC_pAp_App->SetPairs(pApmc, nullptr);
    std::cout << "Get CF \n";
    CFMC_pAp_App->GetCorrelations("pApMC");
    std::cout << "Write Output \n";

    if(strcmp(addon, add1)==0){
    std::cout << "Sphericity [0.,0.3] \n";
    CFData_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st1.root", foldernamedata.Data()));
    CFMC_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st1.root", foldernamemc.Data()));
  }else if(strcmp(addon, add2)==0){
    std::cout << "Sphericity [0.3,0.7] \n";
    CFData_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st2.root", foldernamedata.Data()));
    CFMC_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st2.root", foldernamemc.Data()));
  }else if(strcmp(addon, add3)==0){
    std::cout << "Sphericity [0.7,1.] \n";
    CFData_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st3.root", foldernamedata.Data()));
    CFMC_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st3.root", foldernamemc.Data()));
    }else if(strcmp(addon, add4)==0){
    std::cout << "No Sphericity Cuts \n";
    CFData_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_full.root", foldernamedata.Data()));
    CFMC_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_full.root", foldernamemc.Data()));
      }
    else if(strcmp(addon, add5)==0){
      std::cout << "Sphericity Cuts [0.9-1]\n";
      CFData_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_st5.root", foldernamedata.Data()));
      CFMC_pAp_App->WriteOutput(Form("%sMCCFOutput_pAp_App_st5.root", foldernamemc.Data()));
          }

    std::cout << "pL CF DATA\n";
    CFData_pAL_ApL->SetPairs(pALdata, ApLdata);
    CFData_pAL_ApL->GetCorrelations("pALData");
    std::cout << "pL CF MC\n";
    CFMC_pAL_ApL->SetPairs(pALmc, ApLmc);
    CFMC_pAL_ApL->GetCorrelations("pALMC");

    if(strcmp(addon, add1)==0){
    std::cout << "Sphericity [0.,0.3] \n";
    CFData_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st1.root", foldernamedata.Data()));
    CFMC_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st1.root", foldernamemc.Data()));
  }else if(strcmp(addon, add2)==0){
    std::cout << "Sphericity [0.3,0.7] \n";
    CFData_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st2.root", foldernamedata.Data()));
    CFMC_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st2.root", foldernamemc.Data()));
  }else if(strcmp(addon, add3)==0){
    std::cout << "Sphericity [0.7,1.] \n";
    CFData_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st3.root", foldernamedata.Data()));
    CFMC_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st3.root", foldernamemc.Data()));
  }else if(strcmp(addon, add4)==0){
    std::cout << "No Sphericity Cuts \n";
    CFData_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_full.root", foldernamedata.Data()));
    CFMC_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_full.root", foldernamemc.Data()));
  }else if(strcmp(addon, add5)==0){
    std::cout << "Sphericity Cuts [0.9-1]\n";
    CFData_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_st5.root", foldernamedata.Data()));
    CFMC_pAL_ApL->WriteOutput(Form("%sMCCFOutput_pAL_ApL_st5.root", foldernamemc.Data()));
  }


    std::cout << "LL CF DATA\n";
    CFData_ALL_LAL->SetPairs(LALdata, nullptr);
    CFData_ALL_LAL->GetCorrelations("LALData");
    std::cout << "LL CF MC\n";
    CFMC_ALL_LAL->SetPairs(LALmc, nullptr);
    CFMC_ALL_LAL->GetCorrelations("LALMC");
    if(strcmp(addon, add1)==0){
    std::cout << "Sphericity [0.,0.3] \n";
    CFData_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st1.root", foldernamedata.Data()));
    CFMC_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st1.root", foldernamemc.Data()));
  }else if(strcmp(addon, add2)==0){
    std::cout << "Sphericity [0.3,0.7] \n";
    CFData_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st2.root", foldernamedata.Data()));
    CFMC_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st2.root", foldernamemc.Data()));
    }else if(strcmp(addon, add3)==0){
    std::cout << "Sphericity [0.7,1.] \n";
    CFData_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st3.root", foldernamedata.Data()));
    CFMC_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st3.root", foldernamemc.Data()));
    }else if(strcmp(addon, add4)==0){
    std::cout << "No Sphericity Cuts \n";
    CFData_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_full.root", foldernamedata.Data()));
    CFMC_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_full.root", foldernamemc.Data()));
    }else if(strcmp(addon, add5)==0){
    std::cout << "Sphericity Cuts [0.9-1]\n";
    CFData_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_st5.root", foldernamedata.Data()));
    CFMC_ALL_LAL->WriteOutput(Form("%sMCCFOutput_LAL_ALL_st5.root", foldernamemc.Data()));
    }


}

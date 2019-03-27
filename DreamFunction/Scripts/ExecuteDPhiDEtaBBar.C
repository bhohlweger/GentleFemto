#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"
#include <iostream>
#include "stdlib.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  bool isMC=false;
//  bool epos=false;

  int decider = (atoi(argv[4])) ? atoi(argv[4]) : 0;
    if (decider > 0) {
       isMC = true;
    } else {
       isMC = false;
    }
    // int deciderep = (atoi(argv[5])) ? atoi(argv[5]) : 0;
    //   if (deciderep > 0) {
    //      epos = true;
    //   } else {
    //      epos = false;
    //   }

  TString add1="1";
  TString add2="2";
  TString add3="3";
  TString add4="4";
  TString add5="5";

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReaddEtadPhiHists(0,filename,prefix,addon);//first entry for mT bins
  //GetdEtadPhiDistribution(iPart1,iPart2,iAPart1,iAPart2,imT)

  DreamdEtadPhi* dEtadPhiProtonAntiProton = DreamFile->GetdEtadPhiDistributionSingle(0,1,-1);
  // DreamdEtadPhi* dEtadPhiProtonProtonmT1 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,0);
  // DreamdEtadPhi* dEtadPhiProtonProtonmT2 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,1);
  // DreamdEtadPhi* dEtadPhiProtonProtonmT3 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,2);
  // DreamdEtadPhi* dEtadPhiProtonProtonmT4 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,3);

  dEtadPhiProtonAntiProton->ShiftAbovePhi();


  dEtadPhiProtonAntiProton->DivideSEandME();


  dEtadPhiProtonAntiProton->ProjectionY();


  // dEtadPhiProtonProtonmT1->ShiftAbovePhi();
  // dEtadPhiProtonProtonmT1->DivideSEandME();
  // dEtadPhiProtonProtonmT1->ProjectionY();
  // dEtadPhiProtonProtonmT2->ShiftAbovePhi();
  // dEtadPhiProtonProtonmT2->DivideSEandME();
  // dEtadPhiProtonProtonmT2->ProjectionY();
  // dEtadPhiProtonProtonmT3->ShiftAbovePhi();
  // dEtadPhiProtonProtonmT3->DivideSEandME();
  // dEtadPhiProtonProtonmT3->ProjectionY();
  // dEtadPhiProtonProtonmT4->ShiftAbovePhi();
  // dEtadPhiProtonProtonmT4->DivideSEandME();
  // dEtadPhiProtonProtonmT4->ProjectionY();
  TString foldername;
  TString foldernamemc;

if(!isMC){
  foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/DEtaDPhi/data/";
} else if(isMC){
  foldernamemc = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/DEtaDPhi/mc/";
}
  // TFile *output = new TFile(foldername+"dEtadPhi.root","RECREATE");
 //  if(!file->FindObjectAny("PAP"))
 //  {
 //  file->mkdir("PAP");
 //  }
 // if(!file->FindObjectAny("PAL+APL"))
 // {
 //  file->mkdir("PAL+APL");
 // }
 // if(!file->FindObjectAny("LAL"))
 // {
 //  file->mkdir("LAL");
 // }

 std::cout << "$PWD " << foldername << std::endl;
  TFile* output = TFile::Open(foldername+"dEtadPhi.root", "RECREATE");
  // output->mkdir("PAP");
  // output->cd("PAP");

  TList* pApList = new TList();
   output->cd();
   if(strcmp(addon, add1)==0){
   pApList->SetName("ProtonAntiProton_Data_st1");
   pApList->SetOwner();
   dEtadPhiProtonAntiProton->WriteOutput(pApList,"ppIntegrated");
   pApList->Write("ProtonAntiProton_Data_st1",1);
 } else if(strcmp(addon, add2)==0){
   pApList->SetName("ProtonAntiProton_Data_st2");
   pApList->SetOwner();
   dEtadPhiProtonAntiProton->WriteOutput(pApList,"ppIntegrated");
   pApList->Write("ProtonAntiProton_Data_st2",1);
 } else if(strcmp(addon, add3)==0){
   pApList->SetName("ProtonAntiProton_Data_st3");
   pApList->SetOwner();
   dEtadPhiProtonAntiProton->WriteOutput(pApList,"ppIntegrated");
   pApList->Write("ProtonAntiProton_Data_st3",1);
 } else if(strcmp(addon, add4)==0){
   pApList->SetName("ProtonAntiProton_Data_full");
   pApList->SetOwner();
   dEtadPhiProtonAntiProton->WriteOutput(pApList,"ppIntegrated");
   pApList->Write("ProtonAntiProton_Data_full",1);
 }

  // // dEtadPhiProtonProtonmT1->WriteOutput(ppList,"ppmT1");
  // // dEtadPhiProtonProtonmT2->WriteOutput(ppList,"ppmT2");
  // // dEtadPhiProtonProtonmT3->WriteOutput(ppList,"ppmT3");
  // // dEtadPhiProtonProtonmT4->WriteOutput(ppList,"ppmT4");

//p-antiLambda + antip-Lambda
   DreamdEtadPhi* dEtadPhiProtonAntiLambda= DreamFile->GetdEtadPhiDistribution(0,3,1,2,-1);
  // // DreamdEtadPhi* dEtadPhiProtonLambdamT1 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,0);
  // // DreamdEtadPhi* dEtadPhiProtonLambdamT2 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,1);
  // // DreamdEtadPhi* dEtadPhiProtonLambdamT3 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,2);
  // // DreamdEtadPhi* dEtadPhiProtonLambdamT4 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,3);
  //
   dEtadPhiProtonAntiLambda->ShiftAbovePhi();
   dEtadPhiProtonAntiLambda->DivideSEandME();
   dEtadPhiProtonAntiLambda->ProjectionY();
  // // dEtadPhiProtonLambdamT1->ShiftAbovePhi();
  // // dEtadPhiProtonLambdamT1->DivideSEandME();
  // // dEtadPhiProtonLambdamT1->ProjectionY();
  // // dEtadPhiProtonLambdamT2->ShiftAbovePhi();
  // // dEtadPhiProtonLambdamT2->DivideSEandME();
  // // dEtadPhiProtonLambdamT2->ProjectionY();
  // // dEtadPhiProtonLambdamT3->ShiftAbovePhi();
  // // dEtadPhiProtonLambdamT3->DivideSEandME();
  // // dEtadPhiProtonLambdamT3->ProjectionY();
  // // dEtadPhiProtonLambdamT4->ShiftAbovePhi();
  // // dEtadPhiProtonLambdamT4->DivideSEandME();
  // // dEtadPhiProtonLambdamT4->ProjectionY();
  //
   output->cd();
   if(strcmp(addon, add1)==0){
   TList* pALList = new TList();
   pALList->SetName("ProtonAntiLambda_st1");
   pALList->SetOwner();
   dEtadPhiProtonAntiLambda->WriteOutput(pALList,"pLIntegrated");
   pALList->Write("ProtonAntiLambda_st1",1);
 } else if(strcmp(addon, add2)==0){
   TList* pALList = new TList();
   pALList->SetName("ProtonAntiLambda_st2");
   pALList->SetOwner();
   dEtadPhiProtonAntiLambda->WriteOutput(pALList,"pLIntegrated");
   pALList->Write("ProtonAntiLambda_st2",1);
 } else if(strcmp(addon, add3)==0){
   TList* pALList = new TList();
   pALList->SetName("ProtonAntiLambda_st3");
   pALList->SetOwner();
   dEtadPhiProtonAntiLambda->WriteOutput(pALList,"pLIntegrated");
   pALList->Write("ProtonAntiLambda_st3",1);
 } else if(strcmp(addon, add4)==0){
   TList* pALList = new TList();
   pALList->SetName("ProtonAntiLambda_full");
   pALList->SetOwner();
   dEtadPhiProtonAntiLambda->WriteOutput(pALList,"pLIntegrated");
   pALList->Write("ProtonAntiLambda_full",1);
 }
  // // dEtadPhiProtonLambdamT1->WriteOutput(pLList,"pLmT1");
  // // dEtadPhiProtonLambdamT2->WriteOutput(pLList,"pLmT2");
  // // dEtadPhiProtonLambdamT3->WriteOutput(pLList,"pLmT3");
  // // dEtadPhiProtonLambdamT4->WriteOutput(pLList,"pLmT4");



   //Lambda-antiLambda
      DreamdEtadPhi* dEtadPhiLambdaAntiLambda= DreamFile->GetdEtadPhiDistributionSingle(2,3,-1);
     // // DreamdEtadPhi* dEtadPhiProtonLambdamT1 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,0);
     // // DreamdEtadPhi* dEtadPhiProtonLambdamT2 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,1);
     // // DreamdEtadPhi* dEtadPhiProtonLambdamT3 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,2);
     // // DreamdEtadPhi* dEtadPhiProtonLambdamT4 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,3);
     //
      dEtadPhiLambdaAntiLambda->ShiftAbovePhi();
      dEtadPhiLambdaAntiLambda->DivideSEandME();
      dEtadPhiLambdaAntiLambda->ProjectionY();


      output->cd();
   if(strcmp(addon, add1)==0){
      TList* pLALList = new TList();
      pLALList->SetName("LambdaAntiLambda_st1");
      pLALList->SetOwner();
      dEtadPhiLambdaAntiLambda->WriteOutput(pLALList,"LALIntegrated");
      pLALList->Write("LambdaAntiLambda_st1",1);
    } else if(strcmp(addon, add2)==0){
      TList* pLALList = new TList();
      pLALList->SetName("LambdaAntiLambda_st2");
      pLALList->SetOwner();
      dEtadPhiLambdaAntiLambda->WriteOutput(pLALList,"LALIntegrated");
      pLALList->Write("LambdaAntiLambda_st2",1);
    } else if(strcmp(addon, add3)==0){
      TList* pLALList = new TList();
      pLALList->SetName("LambdaAntiLambda_st3");
      pLALList->SetOwner();
      dEtadPhiLambdaAntiLambda->WriteOutput(pLALList,"LALIntegrated");
      pLALList->Write("LambdaAntiLambda_st3",1);
    } else if(strcmp(addon, add4)==0){
      TList* pLALList = new TList();
      pLALList->SetName("LambdaAntiLambda_full");
      pLALList->SetOwner();
      dEtadPhiLambdaAntiLambda->WriteOutput(pLALList,"LALIntegrated");
      pLALList->Write("LambdaAntiLambda_full",1);
    }
     // // dEtadPhiProtonLambdamT1->WriteOutput(pLList,"pLmT1");
     // // dEtadPhiProtonLambdamT2->WriteOutput(pLList,"pLmT2");
     // // dEtadPhiProtonLambdamT3->WriteOutput(pLList,"pLmT3");
     // // dEtadPhiProtonLambdamT4->WriteOutput(pLList,"pLmT4");

   output->Close();

  return 1;
}

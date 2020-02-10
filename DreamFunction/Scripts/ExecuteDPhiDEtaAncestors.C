#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"
#include <iostream>
#include "stdlib.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReaddEtadPhiHistsAncestors(0,filename,prefix,addon);//first entry for mT bins

//Proton-Proton
  DreamdEtadPhi* dEtadPhiProtonProtonCommon = DreamFile->GetdEtadPhiDistributionCommon(0,0,1,1,-1);
  dEtadPhiProtonProtonCommon->ShiftAbovePhi();
  dEtadPhiProtonProtonCommon->DivideSEandME();
  dEtadPhiProtonProtonCommon->ProjectionY();

  DreamdEtadPhi* dEtadPhiProtonProtonNonCommon = DreamFile->GetdEtadPhiDistributionNonCommon(0,0,1,1,-1);
  dEtadPhiProtonProtonNonCommon->ShiftAbovePhi();
  dEtadPhiProtonProtonNonCommon->DivideSEandME();
  dEtadPhiProtonProtonNonCommon->ProjectionY();


  TFile* output = TFile::Open("dEtadPhiAncestors.root", "RECREATE");
  TList* ppList = new TList();
  ppList->SetName("ProtonProton");
  ppList->SetOwner();
  dEtadPhiProtonProtonCommon->WriteOutput(ppList,"ppCommon");
  dEtadPhiProtonProtonNonCommon->WriteOutput(ppList,"ppNonCommon");
  ppList->Write("ProtonProton",1);

//Proton-AntiProton
  DreamdEtadPhi* dEtadPhiProtonAntiProtonCommon = DreamFile->GetdEtadPhiDistributionSingleCommon(0,1,-1);
  dEtadPhiProtonAntiProtonCommon->ShiftAbovePhi();
  dEtadPhiProtonAntiProtonCommon->DivideSEandME();
  dEtadPhiProtonAntiProtonCommon->ProjectionY();

  DreamdEtadPhi* dEtadPhiProtonAntiProtonNonCommon = DreamFile->GetdEtadPhiDistributionSingleNonCommon(0,1,-1);
  dEtadPhiProtonAntiProtonNonCommon->ShiftAbovePhi();
  dEtadPhiProtonAntiProtonNonCommon->DivideSEandME();
  dEtadPhiProtonAntiProtonNonCommon->ProjectionY();

  output->cd();
  TList* pApList = new TList();
  pApList->SetName("ProtonAntiProton");
  pApList->SetOwner();
  dEtadPhiProtonAntiProtonCommon->WriteOutput(pApList,"pApCommon");
  dEtadPhiProtonAntiProtonNonCommon->WriteOutput(pApList,"pApNonCommon");
  pApList->Write("ProtonAntiProton",1);

//p-Lambda
  DreamdEtadPhi* dEtadPhiProtonLambdaCommon= DreamFile->GetdEtadPhiDistributionCommon(0,2,1,3,-1);
  dEtadPhiProtonLambdaCommon->ShiftAbovePhi();
  dEtadPhiProtonLambdaCommon->DivideSEandME();
  dEtadPhiProtonLambdaCommon->ProjectionY();

  DreamdEtadPhi* dEtadPhiProtonLambdaNonCommon= DreamFile->GetdEtadPhiDistributionNonCommon(0,2,1,3,-1);
  dEtadPhiProtonLambdaNonCommon->ShiftAbovePhi();
  dEtadPhiProtonLambdaNonCommon->DivideSEandME();
  dEtadPhiProtonLambdaNonCommon->ProjectionY();

  output->cd();
  TList* pLList = new TList();
  pLList->SetName("ProtonLambda");
  pLList->SetOwner();
  dEtadPhiProtonLambdaCommon->WriteOutput(pLList,"pLCommon");
  dEtadPhiProtonLambdaNonCommon->WriteOutput(pLList,"pLNonCommon");
  pLList->Write("ProtonLambda",1);

//p-antiLambda
   DreamdEtadPhi* dEtadPhiProtonAntiLambdaCommon= DreamFile->GetdEtadPhiDistributionCommon(0,3,1,2,-1);
   dEtadPhiProtonAntiLambdaCommon->ShiftAbovePhi();
   dEtadPhiProtonAntiLambdaCommon->DivideSEandME();
   dEtadPhiProtonAntiLambdaCommon->ProjectionY();

   DreamdEtadPhi* dEtadPhiProtonAntiLambdaNonCommon= DreamFile->GetdEtadPhiDistributionNonCommon(0,3,1,2,-1);
   dEtadPhiProtonAntiLambdaNonCommon->ShiftAbovePhi();
   dEtadPhiProtonAntiLambdaNonCommon->DivideSEandME();
   dEtadPhiProtonAntiLambdaNonCommon->ProjectionY();

  output->cd();
  TList* pALList = new TList();
  pALList->SetName("ProtonAntiLambda");
  pALList->SetOwner();
  dEtadPhiProtonAntiLambdaCommon->WriteOutput(pALList,"pALCommon");
  dEtadPhiProtonAntiLambdaNonCommon->WriteOutput(pALList,"pALNonCommon");
  pALList->Write("ProtonAntiLambda",1);

   //Lambda-Lambda
  DreamdEtadPhi* dEtadPhiLambdaLambdaCommon= DreamFile->GetdEtadPhiDistributionCommon(2,2,3,3,-1);
  dEtadPhiLambdaLambdaCommon->ShiftAbovePhi();
  dEtadPhiLambdaLambdaCommon->DivideSEandME();
  dEtadPhiLambdaLambdaCommon->ProjectionY();

  DreamdEtadPhi* dEtadPhiLambdaLambdaNonCommon= DreamFile->GetdEtadPhiDistributionNonCommon(2,2,3,3,-1);
  dEtadPhiLambdaLambdaNonCommon->ShiftAbovePhi();
  dEtadPhiLambdaLambdaNonCommon->DivideSEandME();
  dEtadPhiLambdaLambdaNonCommon->ProjectionY();

  output->cd();
  TList* LLList = new TList();
  LLList->SetName("LambdaLambda");
  LLList->SetOwner();
  dEtadPhiLambdaLambdaCommon->WriteOutput(LLList,"LLCommon");
  dEtadPhiLambdaLambdaNonCommon->WriteOutput(LLList,"LLNonCommon");
  LLList->Write("LambdaLambda",1);

   //Lambda-antiLambda
  DreamdEtadPhi* dEtadPhiLambdaAntiLambdaCommon= DreamFile->GetdEtadPhiDistributionSingleCommon(2,3,-1);
  dEtadPhiLambdaAntiLambdaCommon->ShiftAbovePhi();
  dEtadPhiLambdaAntiLambdaCommon->DivideSEandME();
  dEtadPhiLambdaAntiLambdaCommon->ProjectionY();

  DreamdEtadPhi* dEtadPhiLambdaAntiLambdaNonCommon= DreamFile->GetdEtadPhiDistributionSingleNonCommon(2,3,-1);
  dEtadPhiLambdaAntiLambdaNonCommon->ShiftAbovePhi();
  dEtadPhiLambdaAntiLambdaNonCommon->DivideSEandME();
  dEtadPhiLambdaAntiLambdaNonCommon->ProjectionY();


  output->cd();
  TList* LALList = new TList();
  LALList->SetName("LambdaAntiLambda");
  LALList->SetOwner();
  dEtadPhiLambdaAntiLambdaCommon->WriteOutput(LALList,"LALCommon");
  dEtadPhiLambdaAntiLambdaNonCommon->WriteOutput(LALList,"LALNonCommon");
  LALList->Write("LambdaAntiLambda",1);

  output->Close();

  return 1;
}

#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"
#include <iostream>
#include "stdlib.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReaddEtadPhiHists(0,filename,prefix,addon);//first entry for mT bins

//Proton-Proton
  DreamdEtadPhi* dEtadPhiProtonProton = DreamFile->GetdEtadPhiDistribution(0,0,1,1,-1);
  dEtadPhiProtonProton->ShiftAbovePhi();
  dEtadPhiProtonProton->DivideSEandME();
  dEtadPhiProtonProton->ProjectionY();

  TFile* output = TFile::Open("dEtadPhi.root", "RECREATE");
  TList* ppList = new TList();
  ppList->SetName("ProtonProton");
  ppList->SetOwner();
  dEtadPhiProtonProton->WriteOutput(ppList,"ppIntegrated");
  ppList->Write("ProtonProton",1);

//Proton-AntiProton
  DreamdEtadPhi* dEtadPhiProtonAntiProton = DreamFile->GetdEtadPhiDistributionSingle(0,1,-1);
  dEtadPhiProtonAntiProton->ShiftAbovePhi();
  dEtadPhiProtonAntiProton->DivideSEandME();
  dEtadPhiProtonAntiProton->ProjectionY();

  output->cd();
  TList* pApList = new TList();
  pApList->SetName("ProtonAntiProton");
  pApList->SetOwner();
  dEtadPhiProtonAntiProton->WriteOutput(pApList,"pApIntegrated");
  pApList->Write("ProtonAntiProton",1);

//p-Lambda
  DreamdEtadPhi* dEtadPhiProtonLambda= DreamFile->GetdEtadPhiDistribution(0,2,1,3,-1);
  dEtadPhiProtonLambda->ShiftAbovePhi();
  dEtadPhiProtonLambda->DivideSEandME();
  dEtadPhiProtonLambda->ProjectionY();

  output->cd();
  TList* pLList = new TList();
  pLList->SetName("ProtonLambda");
  pLList->SetOwner();
  dEtadPhiProtonLambda->WriteOutput(pLList,"pLIntegrated");
  pLList->Write("ProtonLambda",1);

//p-antiLambda
   DreamdEtadPhi* dEtadPhiProtonAntiLambda= DreamFile->GetdEtadPhiDistribution(0,3,1,2,-1);
   dEtadPhiProtonAntiLambda->ShiftAbovePhi();
   dEtadPhiProtonAntiLambda->DivideSEandME();
   dEtadPhiProtonAntiLambda->ProjectionY();

  output->cd();
  TList* pALList = new TList();
  pALList->SetName("ProtonAntiLambda");
  pALList->SetOwner();
  dEtadPhiProtonAntiLambda->WriteOutput(pALList,"pALIntegrated");
  pALList->Write("ProtonAntiLambda",1);

   //Lambda-Lambda
  DreamdEtadPhi* dEtadPhiLambdaLambda= DreamFile->GetdEtadPhiDistribution(2,2,3,3,-1);
  dEtadPhiLambdaLambda->ShiftAbovePhi();
  dEtadPhiLambdaLambda->DivideSEandME();
  dEtadPhiLambdaLambda->ProjectionY();

  output->cd();
  TList* LLList = new TList();
  LLList->SetName("LambdaLambda");
  LLList->SetOwner();
  dEtadPhiLambdaLambda->WriteOutput(LLList,"LLIntegrated");
  LLList->Write("LambdaLambda",1);

   //Lambda-antiLambda
  DreamdEtadPhi* dEtadPhiLambdaAntiLambda= DreamFile->GetdEtadPhiDistributionSingle(2,3,-1);
  dEtadPhiLambdaAntiLambda->ShiftAbovePhi();
  dEtadPhiLambdaAntiLambda->DivideSEandME();
  dEtadPhiLambdaAntiLambda->ProjectionY();


  output->cd();
  TList* LALList = new TList();
  LALList->SetName("LambdaAntiLambda");
  LALList->SetOwner();
  dEtadPhiLambdaAntiLambda->WriteOutput(LALList,"LALIntegrated");
  LALList->Write("LambdaAntiLambda",1);

  output->Close();

  return 1;
}

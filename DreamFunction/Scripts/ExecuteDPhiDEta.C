#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->ReaddEtadPhiHists(4,filename,prefix,addon);
  DreamdEtadPhi* dEtadPhiProtonProton = DreamFile->GetdEtadPhiDistribution(0,0,1,1,-1);
  DreamdEtadPhi* dEtadPhiProtonProtonmT1 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,0);
  DreamdEtadPhi* dEtadPhiProtonProtonmT2 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,1);
  DreamdEtadPhi* dEtadPhiProtonProtonmT3 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,2);
  DreamdEtadPhi* dEtadPhiProtonProtonmT4 = DreamFile->GetdEtadPhiDistribution(0,0,1,1,3);

  dEtadPhiProtonProton->ShiftAbovePhi();
  dEtadPhiProtonProton->DivideSEandME();
  dEtadPhiProtonProton->ProjectionY();
  dEtadPhiProtonProtonmT1->ShiftAbovePhi();
  dEtadPhiProtonProtonmT1->DivideSEandME();
  dEtadPhiProtonProtonmT1->ProjectionY();
  dEtadPhiProtonProtonmT2->ShiftAbovePhi();
  dEtadPhiProtonProtonmT2->DivideSEandME();
  dEtadPhiProtonProtonmT2->ProjectionY();
  dEtadPhiProtonProtonmT3->ShiftAbovePhi();
  dEtadPhiProtonProtonmT3->DivideSEandME();
  dEtadPhiProtonProtonmT3->ProjectionY();
  dEtadPhiProtonProtonmT4->ShiftAbovePhi();
  dEtadPhiProtonProtonmT4->DivideSEandME();
  dEtadPhiProtonProtonmT4->ProjectionY();

  TFile* output = TFile::Open("dEtadPhi.root", "RECREATE");
  output->cd();
  TList* ppList = new TList();
  ppList->SetName("ProtonProton");
  ppList->SetOwner();
  dEtadPhiProtonProton->WriteOutput(ppList,"ppIntegrated");
  dEtadPhiProtonProtonmT1->WriteOutput(ppList,"ppmT1");
  dEtadPhiProtonProtonmT2->WriteOutput(ppList,"ppmT2");
  dEtadPhiProtonProtonmT3->WriteOutput(ppList,"ppmT3");
  dEtadPhiProtonProtonmT4->WriteOutput(ppList,"ppmT4");
  ppList->Write("ProtonProton",1);

  DreamdEtadPhi* dEtadPhiProtonLambda= DreamFile->GetdEtadPhiDistribution(0,2,1,3,-1);
  DreamdEtadPhi* dEtadPhiProtonLambdamT1 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,0);
  DreamdEtadPhi* dEtadPhiProtonLambdamT2 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,1);
  DreamdEtadPhi* dEtadPhiProtonLambdamT3 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,2);
  DreamdEtadPhi* dEtadPhiProtonLambdamT4 = DreamFile->GetdEtadPhiDistribution(0,2,1,3,3);

  dEtadPhiProtonLambda->ShiftAbovePhi();
  dEtadPhiProtonLambda->DivideSEandME();
  dEtadPhiProtonLambda->ProjectionY();
  dEtadPhiProtonLambdamT1->ShiftAbovePhi();
  dEtadPhiProtonLambdamT1->DivideSEandME();
  dEtadPhiProtonLambdamT1->ProjectionY();
  dEtadPhiProtonLambdamT2->ShiftAbovePhi();
  dEtadPhiProtonLambdamT2->DivideSEandME();
  dEtadPhiProtonLambdamT2->ProjectionY();
  dEtadPhiProtonLambdamT3->ShiftAbovePhi();
  dEtadPhiProtonLambdamT3->DivideSEandME();
  dEtadPhiProtonLambdamT3->ProjectionY();
  dEtadPhiProtonLambdamT4->ShiftAbovePhi();
  dEtadPhiProtonLambdamT4->DivideSEandME();
  dEtadPhiProtonLambdamT4->ProjectionY();

  output->cd();
  TList* pLList = new TList();
  pLList->SetName("ProtonLambda");
  pLList->SetOwner();
  dEtadPhiProtonLambda->WriteOutput(pLList,"pLIntegrated");
  dEtadPhiProtonLambdamT1->WriteOutput(pLList,"pLmT1");
  dEtadPhiProtonLambdamT2->WriteOutput(pLList,"pLmT2");
  dEtadPhiProtonLambdamT3->WriteOutput(pLList,"pLmT3");
  dEtadPhiProtonLambdamT4->WriteOutput(pLList,"pLmT4");
  pLList->Write("ProtonLambda",1);

  return 1;
}

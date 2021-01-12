#include "ReadDreamFile.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  TString appendix = TString::Format("%s", argv[2]);
  TString suffix = TString::Format("%s", argv[3]);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->ReaddEtadPhiHists(0, filename, appendix.Data(), suffix.Data());
  
  //Proton-D-
  DreamdEtadPhi* dEtadPhiProtonDminus = DreamFile->GetdEtadPhiDistribution(0,3,1,2,-1);
  dEtadPhiProtonDminus->ShiftAbovePhi();
  dEtadPhiProtonDminus->DivideSEandME();
  dEtadPhiProtonDminus->ProjectionY();
  
  //Proton-D+
  DreamdEtadPhi* dEtadPhiProtonDplus = DreamFile->GetdEtadPhiDistribution(0,2,1,3,-1);
  dEtadPhiProtonDplus->ShiftAbovePhi();
  dEtadPhiProtonDplus->DivideSEandME();
  dEtadPhiProtonDplus->ProjectionY();
  
  TString appendixFile =
        (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());
  TFile* output = TFile::Open(TString::Format("dEtadPhi%s.root", appendixFile.Data()), "RECREATE");
  TList* pDminusList = new TList();
  pDminusList->SetName("ProtonDplus");
  pDminusList->SetOwner();
  dEtadPhiProtonDminus->WriteOutput(pDminusList,"pD-");
  pDminusList->Write("ProtonDminus",1);

  TList* pDplusList = new TList();
  pDplusList->SetName("ProtonDplus");
  pDplusList->SetOwner();
  dEtadPhiProtonDplus->WriteOutput(pDplusList,"pD+");
  pDplusList->Write("ProtonDplus",1);
  output->Close();
  
  return 1;
}

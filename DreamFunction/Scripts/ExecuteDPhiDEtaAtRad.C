#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"
#include "TCanvas.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };

  ReadDreamFile* DreamFile = new ReadDreamFile(2, 2);
  DreamFile->ReaddEtadPhiAtRadHists(2, filename, prefix, addon);
  //pp Loop
  TFile* output = TFile::Open("dEtadPhiAtRad.root", "RECREATE");
  auto cpp = new TCanvas("cpp","cpp");
  cpp->Divide(3,3);
  auto cppLowkStar = new TCanvas("cppLowkStar","cppLowkStar");
  cppLowkStar->Divide(3,3);
  auto cppProjY = new TCanvas("cppProjY","cppProjY");
  cppProjY->Divide(3,3);
  auto cppLowkStarProjY = new TCanvas("cppLowkStarProjY","cppLowkStarProjY");
  cppLowkStarProjY->Divide(3,3);
  TList* ppList = new TList();
  ppList->SetOwner();
  ppList->SetName("pp");
  ppList->Add(cpp);
  ppList->Add(cppLowkStar);
  ppList->Add(cppProjY);
  ppList->Add(cppLowkStarProjY);

  for (int iRad = 0; iRad < 9; ++iRad) {
    DreamdEtadPhi *ppRad = DreamFile->GetdEtadPhiAtRadDistribution(0,0,1,1,iRad,0,false);
    ppRad->DivideSEandME();
    ppRad->Draw2D((TPad*)cpp->cd(iRad+1),TPCradii[iRad]);
    ppRad->ProjectionY();
    ppRad->DrawProjectionY((TPad*)cppProjY->cd(iRad+1),TPCradii[iRad]);
    ppRad->WriteOutput(ppList,Form("ppRad%u",iRad));

    DreamdEtadPhi *ppRadLowkStar = DreamFile->GetdEtadPhiAtRadDistribution(0,0,1,1,iRad,0,true);
    ppRadLowkStar->DivideSEandME();
    ppRadLowkStar->Draw2D((TPad*)cppLowkStar->cd(iRad+1),TPCradii[iRad]);
    ppRadLowkStar->ProjectionY();
    ppRadLowkStar->DrawProjectionY((TPad*)cppLowkStarProjY->cd(iRad+1),TPCradii[iRad]);
    ppRadLowkStar->WriteOutput(ppList,Form("ppLowkStarRad%u",iRad));
  }
  output->cd();
  ppList->Write("ppList",1);
  return 0;
}


#include "ReadDreamFile.h"
#include "DreamdEtadPhi.h"
#include "TCanvas.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
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

  auto cpL = new TCanvas("cpL","cpL");
  cpL->Divide(3,3);
  auto cpLLowkStar = new TCanvas("cpLLowkStar","cpLLowkStar");
  cpLLowkStar->Divide(3,3);
  auto cpLProjY = new TCanvas("cpLProjY","cpLProjY");
  cpLProjY->Divide(3,3);
  auto cpLLowkStarProjY = new TCanvas("cpLLowkStarProjY","cpLLowkStarProjY");
  cpLLowkStarProjY->Divide(3,3);
  TList* pLList = new TList();
  pLList->SetOwner();
  pLList->SetName("pL");
  pLList->Add(cpL);
  pLList->Add(cpLLowkStar);
  pLList->Add(cpLProjY);
  pLList->Add(cpLLowkStarProjY);

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

    DreamdEtadPhi *pLRad = DreamFile->GetdEtadPhiAtRadDistribution(0,2,1,3,iRad,0,false);
    pLRad->DivideSEandME();
    pLRad->Draw2D((TPad*)cpL->cd(iRad+1),TPCradii[iRad]);
    pLRad->ProjectionY();
    pLRad->DrawProjectionY((TPad*)cpLProjY->cd(iRad+1),TPCradii[iRad]);
    pLRad->WriteOutput(pLList,Form("pLRad%u",iRad));

    DreamdEtadPhi *pLRadLowkStar = DreamFile->GetdEtadPhiAtRadDistribution(0,2,1,3,iRad,0,true);
    pLRadLowkStar->DivideSEandME();
    pLRadLowkStar->Draw2D((TPad*)cpLLowkStar->cd(iRad+1),TPCradii[iRad]);
    pLRadLowkStar->ProjectionY();
    pLRadLowkStar->DrawProjectionY((TPad*)cpLLowkStarProjY->cd(iRad+1),TPCradii[iRad]);
    pLRadLowkStar->WriteOutput(pLList,Form("pLLowkStarRad%u",iRad));

  }
  output->cd();
  ppList->Write("ppList",1);
  pLList->Write("pLList",1);
  return 0;
}


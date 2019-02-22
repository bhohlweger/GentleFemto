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
  auto cpp = new TCanvas("cpp", "cpp");
  cpp->Divide(3, 3);
  auto cppLowkStar = new TCanvas("cppLowkStar", "cppLowkStar");
  cppLowkStar->Divide(3, 3);
  auto cppProjY = new TCanvas("cppProjY", "cppProjY");
  cppProjY->Divide(3, 3);
  auto cppLowkStarProjY = new TCanvas("cppLowkStarProjY", "cppLowkStarProjY");
  cppLowkStarProjY->Divide(3, 3);
  TList* ppList = new TList();
  ppList->SetOwner();
  ppList->SetName("pp");
  ppList->Add(cpp);
  ppList->Add(cppLowkStar);
  ppList->Add(cppProjY);
  ppList->Add(cppLowkStarProjY);

  auto cpLPPi = new TCanvas("cpLPPi", "cpLPPi");
  cpLPPi->Divide(3, 3);
  auto cpLPPiLowkStar = new TCanvas("cpLPPiLowkStar", "cpLPPiLowkStar");
  cpLPPiLowkStar->Divide(3, 3);
  auto cpLPPiProjY = new TCanvas("cpLPPiProjY", "cpLPPiProjY");
  cpLPPiProjY->Divide(3, 3);
  auto cpLPPiLowkStarProjY = new TCanvas("cpLPPiLowkStarProjY",
                                         "cpLPPiLowkStarProjY");
  cpLPPiLowkStarProjY->Divide(3, 3);
  TList* pLProtonPionList = new TList();
  pLProtonPionList->SetOwner();
  pLProtonPionList->SetName("pLProtonPion");
  pLProtonPionList->Add(cpLPPi);
  pLProtonPionList->Add(cpLPPiLowkStar);
  pLProtonPionList->Add(cpLPPiProjY);
  pLProtonPionList->Add(cpLPPiLowkStarProjY);

  auto cpLPP = new TCanvas("cpLPP", "cpLPP");
  cpLPP->Divide(3, 3);
  auto cpLPPLowkStar = new TCanvas("cpLPPLowkStar", "cpLPPLowkStar");
  cpLPPLowkStar->Divide(3, 3);
  auto cpLPPProjY = new TCanvas("cpLPPProjY", "cpLPPProjY");
  cpLPPProjY->Divide(3, 3);
  auto cpLPPLowkStarProjY = new TCanvas("cpLPPLowkStarProjY",
                                        "cpLPPLowkStarProjY");
  cpLPPLowkStarProjY->Divide(3, 3);
  TList* pLProtonProtonList = new TList();
  pLProtonProtonList->SetOwner();
  pLProtonProtonList->SetName("pLProtonPion");
  pLProtonProtonList->Add(cpLPP);
  pLProtonProtonList->Add(cpLPPLowkStar);
  pLProtonProtonList->Add(cpLPPProjY);
  pLProtonProtonList->Add(cpLPPLowkStarProjY);

  for (int iRad = 0; iRad < 9; ++iRad) {
    DreamdEtadPhi *ppRad = DreamFile->GetdEtadPhiAtRadDistribution(0, 0, 0, 1,
                                                                   1, 0, iRad,
                                                                   false);
    ppRad->DivideSEandME(3);
    ppRad->Draw2D((TPad*) cpp->cd(iRad + 1), TPCradii[iRad]);
    ppRad->ProjectionY();
    ppRad->DrawProjectionY((TPad*) cppProjY->cd(iRad + 1), TPCradii[iRad]);
    ppRad->WriteOutput(ppList, Form("ppRad%u", iRad));


    DreamdEtadPhi *pLRadProtonPion = DreamFile->GetdEtadPhiAtRadDistribution(
        0, 2, 0, 1, 3, 1, iRad, false);
    pLRadProtonPion->DivideSEandME(5);
    pLRadProtonPion->Draw2D((TPad*) cpLPPi->cd(iRad + 1), TPCradii[iRad]);
    pLRadProtonPion->ProjectionY();
    pLRadProtonPion->DrawProjectionY((TPad*) cpLPPiProjY->cd(iRad + 1),
                                     TPCradii[iRad]);
    pLRadProtonPion->WriteOutput(pLProtonPionList, Form("pLPPiRad%u", iRad));


    DreamdEtadPhi *pLRadProtonProton = DreamFile->GetdEtadPhiAtRadDistribution(
        0, 2, 1, 1, 3, 0, iRad, false);
    pLRadProtonProton->DivideSEandME(5);
    pLRadProtonProton->Draw2D((TPad*) cpLPP->cd(iRad + 1), TPCradii[iRad]);
    pLRadProtonProton->ProjectionY();
    pLRadProtonProton->DrawProjectionY((TPad*) cpLPPProjY->cd(iRad + 1),
                                       TPCradii[iRad]);
    pLRadProtonProton->WriteOutput(pLProtonProtonList,
                                   Form("pLPPiRad%u", iRad));

//    DreamdEtadPhi *ppRadLowkStar = DreamFile->GetdEtadPhiAtRadDistribution(
//        0, 0, 0, 1, 1, 0, iRad, true);
//    ppRadLowkStar->DivideSEandME(2);
//    ppRadLowkStar->Draw2D((TPad*) cppLowkStar->cd(iRad + 1), TPCradii[iRad]);
//    ppRadLowkStar->ProjectionY();
//    ppRadLowkStar->DrawProjectionY((TPad*) cppLowkStarProjY->cd(iRad + 1),
//                                   TPCradii[iRad]);
//    ppRadLowkStar->WriteOutput(ppList, Form("ppLowkStarRad%u", iRad));
//
//    DreamdEtadPhi *pLRadProtonProtonLowkStar = DreamFile
//        ->GetdEtadPhiAtRadDistribution(0, 2, 1, 1, 3, 0, iRad, true);
//    pLRadProtonProtonLowkStar->DivideSEandME(4);
//    pLRadProtonProtonLowkStar->Draw2D((TPad*) cpLPPLowkStar->cd(iRad + 1),
//                                      TPCradii[iRad]);
//    pLRadProtonProtonLowkStar->ProjectionY();
//    pLRadProtonProtonLowkStar->DrawProjectionY(
//        (TPad*) cpLPPLowkStarProjY->cd(iRad + 1), TPCradii[iRad]);
//    pLRadProtonProtonLowkStar->WriteOutput(pLProtonProtonList,
//                                           Form("pLPPiLowkStarRad%u", iRad));
//
//    DreamdEtadPhi *pLRadProtonPionLowkStar = DreamFile
//        ->GetdEtadPhiAtRadDistribution(0, 2, 0, 1, 3, 1, iRad, true);
//    pLRadProtonPionLowkStar->DivideSEandME(4);
//    pLRadProtonPionLowkStar->Draw2D((TPad*) cpLPPiLowkStar->cd(iRad + 1),
//                                    TPCradii[iRad]);
//    pLRadProtonPionLowkStar->ProjectionY();
//    pLRadProtonPionLowkStar->DrawProjectionY(
//        (TPad*) cpLPPiLowkStarProjY->cd(iRad + 1), TPCradii[iRad]);
//    pLRadProtonPionLowkStar->WriteOutput(pLProtonPionList,
//                                         Form("pLPPiLowkStarRad%u", iRad));
  }
  output->cd();
  ppList->Write("ppList", 1);
  pLProtonPionList->Write("pLProtonPion", 1);
  pLProtonProtonList->Write("pLProtonProton", 1);
  return 0;
}


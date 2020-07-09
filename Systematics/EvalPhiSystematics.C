#include "CATSInputPhi.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

// #include "SidebandPhi.h"

void PhiEvalSystematics(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int rebin = 2;

  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInputPhi();
  CATSinput->ReadPhiCorrelationFile(InputDir.Data(), trigger.Data(), "0");
  CATSinput->CountPairs(InputDir.Data(), trigger.Data(), "0");
  CATSinput->ObtainCFs(5, 200, 400, rebin, false);
  TString dataHistPhiName = "hCk_ReweightedpPhiMeV_0";
//  TString dataHistSBName = "hCk_ReweightedpPhiSBUpMeV_0";
  auto dataHistPhi = CATSinput->GetCF("pPhi", dataHistPhiName.Data());
  const unsigned int pairCountsDefault = CATSinput->GetFemtoPairs(0, 0.2,
                                                                  "pPhi");
  const int nProtonDefault = CATSinput->GetNProtonTotal();
  const int nPhiDefault = CATSinput->GetNPhi();
  const float purityPhiDefault = CATSinput->GetPhiPurity();

//  auto side = new SidebandPhi();
//  side->SetRebin(rebin * 10);
//  side->SetSideBandFile(InputDir.Data(), trigger.Data(), "0");
//  side->SetNormalizationRange(200, 400);
//  side->SideBandCFs();
//  auto dataHistSB = side->GetSideBands(5);

  DreamSystematics protonphi(DreamSystematics::pPhi);
  protonphi.SetDefaultHist(dataHistPhi);
  protonphi.SetUpperFitRange(500);
  protonphi.SetBarlowUpperRange(500);
  protonphi.SetEstimator(DreamSystematics::Uniform);
//  DreamSystematics protonSB(DreamSystematics::pPhi);
//  protonSB.SetDefaultHist(dataHistSB);
//  protonSB.SetUpperFitRange(500);
//  protonSB.SetBarlowUpperRange(500);
//  protonSB.SetEstimator(DreamSystematics::Uniform);

  for (int i = 1; i <= protonphi.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputPhi();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadPhiCorrelationFile(InputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->CountPairs(InputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    CATSinputVar->ObtainCFs(5, 200, 400, rebin, false);
    protonphi.SetVarHist(
        CATSinputVar->GetCF("pPhi", dataHistPhiName.Data()));
    protonphi.SetPair(pairCountsDefault,
                        CATSinputVar->GetFemtoPairs(0, 0.2, "pPhi"));
    protonphi.SetParticles(nProtonDefault, nPhiDefault,
                             CATSinputVar->GetNProtonTotal(),
                             CATSinputVar->GetNPhi());
//    protonphi.SetPurity(0, purityPhiDefault, 0,
//                          CATSinputVar->GetPhiPurity());
    protonphi.SetPurity(0, purityPhiDefault, 0,
                          CATSinputVar->GetPhiPurity());

//    auto sideVar = new SidebandPhi();
//    sideVar->SetRebin(rebin * 10);
//    sideVar->SetSideBandFile(InputDir.Data(), trigger.Data(), appendixVar.Data());
//    sideVar->SetNormalizationRange(250, 400);
//    sideVar->SideBandCFs();
//    auto dataHistSB = sideVar->GetSideBands(5);

//    protonSB.SetVarHist(sideVar->GetSideBands(5));
//    delete CATSinputVar;
//    delete sideVar;
  }
  protonphi.EvalSystematics();
  protonphi.EvalDifferenceInPairs();
  protonphi.EvalDifferenceInParticles();
  protonphi.EvalDifferenceInPurity();
  protonphi.WriteOutput();

//  protonSB.EvalSystematics();
//  protonSB.WriteOutput("Sidebands");

}

int main(int argc, char* argv[]) {
  PhiEvalSystematics(argv[1], argv[2]);

  return 1;
}

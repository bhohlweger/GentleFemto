#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "SidebandSigma.h"
#include "TCanvas.h"
#include <iostream>

void SigmaEvalSystematics(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int rebin = 2;

  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInputSigma0();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(), "0");
  CATSinput->CountPairs(InputDir.Data(), trigger.Data(), "0");
  CATSinput->ObtainCFs(10, 250, 400, rebin, true);
  TString dataHistSigmaName = "hCk_ReweightedpSigma0MeV_0";
  TString dataHistSBName = "hCk_ReweightedpSigmaSBUpMeV_0";
  auto dataHistSigma = CATSinput->GetCF("pSigma0", dataHistSigmaName.Data());
  const unsigned int pairCountsDefault = CATSinput->GetFemtoPairs(0, 0.2,
                                                                  "pSigma0");
  const int nProtonDefault = CATSinput->GetNProtonTotal();
  const int nSigmaDefault = CATSinput->GetNSigma0();
  const float puritySigmaDefault = CATSinput->GetSigma0Purity();

  auto side = new SidebandSigma();
  side->SetRebin(rebin * 10);
  side->SetSideBandFile(InputDir.Data(), trigger.Data(), "0");
  side->SetNormalizationRange(250, 400);
  side->SideBandCFs();
  auto dataHistSB = side->GetSideBands(5);

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  protonsigma.SetDefaultHist(dataHistSigma);
  protonsigma.SetUpperFitRange(500);
  protonsigma.SetBarlowUpperRange(500);
  protonsigma.SetEstimator(DreamSystematics::Uniform);
  DreamSystematics protonSB(DreamSystematics::pSigma0);
  protonSB.SetDefaultHist(dataHistSB);
  protonSB.SetUpperFitRange(500);
  protonSB.SetBarlowUpperRange(500);
  protonSB.SetEstimator(DreamSystematics::Uniform);

  for (int i = 1; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->CountPairs(InputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    CATSinputVar->ObtainCFs(10, 250, 400, rebin, true);
    protonsigma.SetVarHist(
        CATSinputVar->GetCF("pSigma0", dataHistSigmaName.Data()));
    protonsigma.SetPair(pairCountsDefault,
                        CATSinputVar->GetFemtoPairs(0, 0.2, "pSigma0"));
    protonsigma.SetParticles(nProtonDefault, nSigmaDefault,
                             CATSinputVar->GetNProtonTotal(),
                             CATSinputVar->GetNSigma0());
    protonsigma.SetPurity(0, puritySigmaDefault, 0,
                          CATSinputVar->GetSigma0Purity());

    auto sideVar = new SidebandSigma();
    sideVar->SetRebin(rebin * 10);
    sideVar->SetSideBandFile(InputDir.Data(), trigger.Data(), appendixVar.Data());
    sideVar->SetNormalizationRange(250, 400);
    sideVar->SideBandCFs();
    auto dataHistSB = sideVar->GetSideBands(5);

    protonSB.SetVarHist(sideVar->GetSideBands(5));
    delete CATSinputVar;
    delete sideVar;
  }
  protonsigma.EvalSystematics();
  protonsigma.EvalDifferenceInPairs();
  protonsigma.EvalDifferenceInParticles();
  protonsigma.EvalDifferenceInPurity();
  protonsigma.WriteOutput();

  protonSB.EvalSystematics();
  protonSB.WriteOutput("Sidebands");

}

int main(int argc, char* argv[]) {
  SigmaEvalSystematics(argv[1], argv[2]);

  return 1;
}

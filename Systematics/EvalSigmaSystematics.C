#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

void SigmaEvalSystematics(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int rebin = 1;

  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInputSigma0();
  CATSinput->ReadCorrelationFile(InputDir.Data(), trigger.Data(), "0");
  CATSinput->CountPairs(InputDir.Data(), trigger.Data(), "0");
  CATSinput->ObtainCFs(10, 250, 400, rebin, false);
  TString dataHistSigmaName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHistSigma = CATSinput->GetCF("pSigma0", dataHistSigmaName.Data());
  const unsigned int pairCountsDefault = CATSinput->GetFemtoPairs(0, 0.2,
                                                                  "pSigma0");
  const int nProtonDefault = CATSinput->GetNProtonTotal();
  const int nSigmaDefault = CATSinput->GetNSigma0();
  const float puritySigmaDefault = CATSinput->GetSigma0Purity();

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  protonsigma.SetDefaultHist(dataHistSigma);
  protonsigma.SetUpperFitRange(450);
  protonsigma.SetBarlowUpperRange(400);
  for (int i = 1; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->CountPairs(InputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    CATSinputVar->ObtainCFs(10, 250, 400, rebin, false);
    protonsigma.SetVarHist(
        CATSinputVar->GetCF("pSigma0", dataHistSigmaName.Data()));
    protonsigma.SetPair(pairCountsDefault,
                        CATSinputVar->GetFemtoPairs(0, 0.2, "pSigma0"));
    protonsigma.SetParticles(nProtonDefault, nSigmaDefault,
                             CATSinputVar->GetNProtonTotal(),
                             CATSinputVar->GetNSigma0());
    protonsigma.SetPurity(0, puritySigmaDefault, 0,
                          CATSinputVar->GetSigma0Purity());
    delete CATSinputVar;
  }
  protonsigma.EvalSystematics();
  protonsigma.EvalDifferenceInPairs();
  protonsigma.EvalDifferenceInParticles();
  protonsigma.EvalDifferenceInPurity();
  protonsigma.WriteOutput();
}

int main(int argc, char* argv[]) {
  SigmaEvalSystematics(argv[1], argv[2]);

  return 1;
}

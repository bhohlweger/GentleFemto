#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

void SigmaEvalSystematics(TString InputDir, TString trigger) {
  const int rebin = 3;

  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInputSigma0();
  CATSinput->ReadCorrelationFile(InputDir.Data(), trigger.Data(), "0");
  CATSinput->ObtainCFs(10, 250, 400, rebin);
  TString dataHistSigmaName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHistSigma = CATSinput->GetCF("pSigma0", dataHistSigmaName.Data());

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  protonsigma.SetDefaultHist(dataHistSigma);
  protonsigma.SetUpperFitRange(400);
  for (int i = 1; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->ObtainCFs(10, 250, 400, rebin);
    protonsigma.SetVarHist(
        CATSinputVar->GetCF("pSigma0", dataHistSigmaName.Data()));
    delete CATSinputVar;
  }
  protonsigma.EvalSystematics();
  protonsigma.WriteOutput();
}

int main(int argc, char* argv[]) {
  SigmaEvalSystematics(argv[1], argv[2]);

  return 1;
}

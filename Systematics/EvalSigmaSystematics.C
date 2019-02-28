#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

void SigmaEvalSystematics(TString InputDir, TString appendix) {
  const int rebin = 5;

  DreamPlot::SetStyle();
  auto CATSinput = new CATSInputSigma0();
  auto appendixDefault = TString::Format("%s_0", appendix.Data());
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendixDefault.Data());
  CATSinput->ObtainCFs(10, 340, 440, rebin);
  TString dataHistProtonName = "hCk_ReweightedppMeV_0";
  auto dataHistProton = CATSinput->GetCF("pp", dataHistProtonName.Data());
  TString dataHistSigmaName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHistSigma = CATSinput->GetCF("pSigma0", dataHistSigmaName.Data());

  DreamSystematics protonproton(DreamSystematics::pp);
  protonproton.SetDefaultHist(dataHistProton);
  for (int i = 1; i <= protonproton.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%s_%i", appendix.Data(), i);
    CATSinputVar->ReadSigma0CorrelationFile(InputDir.Data(),
                                            appendixVar.Data());
    CATSinputVar->ObtainCFs(10, 340, 440, rebin);
    protonproton.SetVarHist(
        CATSinputVar->GetCF("pp", dataHistProtonName.Data()));
    delete CATSinputVar;
  }
  protonproton.EvalSystematics();
  protonproton.WriteOutput();

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  protonsigma.SetDefaultHist(dataHistSigma);
  protonsigma.SetUpperFitRange(650);
  for (int i = 1; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%s_%i", appendix.Data(), i);
    CATSinputVar->ReadSigma0CorrelationFile(InputDir.Data(),
                                            appendixVar.Data());
    CATSinputVar->ObtainCFs(10, 340, 440, rebin);
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

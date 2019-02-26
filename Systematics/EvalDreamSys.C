#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix) {
  const int rebin = 5;

  DreamPlot::SetStyle();
  auto CATSinput = new CATSInput();
  CATSinput->ReadCorrelationFile(InputDir.Data(), prefix.Data());
  CATSinput->ObtainCFs(10, 240, 340);
  TString dataHistProtonName = "hCk_RebinnedppMeV_0";
  auto dataHistProton = CATSinput->GetCF("pp", dataHistProtonName.Data());

  DreamSystematics protonproton(DreamSystematics::pp);
  protonproton.SetDefaultHist(dataHistProton);
  const int protonVarStart = 1;
  for (int i = protonVarStart; i < protonVarStart + protonproton.GetNumberOfVars();
      ++i) {
    auto CATSinputVar = new CATSInput( );

    CATSinputVar->ReadCorrelationFile(InputDir.Data(), prefix.Data(), Form("%u",i));

    CATSinputVar->ObtainCFs(10, 240, 340);
    protonproton.SetVarHist(CATSinputVar->GetCF("pp", dataHistProtonName.Data()));
    delete CATSinputVar;
  }
  protonproton.EvalSystematics();
  protonproton.WriteOutput();
  protonproton.DrawAllCF();

}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2]);

  return 1;
}

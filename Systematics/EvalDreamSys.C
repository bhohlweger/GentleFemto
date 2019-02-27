#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix) {
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  const int rebin = 5;

  DreamPlot::SetStyle();
  auto CATSinput = new CATSInput();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix);
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF* CFDef = CATSinput->ObtainCFSyst(10, "ppDef",pp, ApAp);
  auto dataHistProton = CFDef->
  DreamSystematics protonproton(DreamSystematics::pp);
  protonproton.SetDefaultHist(dataHistProton);
  const int protonVarStart = 1;
  for (int i = protonVarStart;
      i < protonVarStart + protonproton.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInput();
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    CATSinputVar->SetNormalization(240, 340);
    DreamCF* CFDef = CATSinputVar->ObtainCFSyst(10, "pp",
                               DreamFile->GetPairDistributions(0, 0, ""),
                               DreamFile->GetPairDistributions(1, 1, ""),
                               pp, ApAp);
    protonproton.SetVarHist(
        CATSinputVar->GetCF("pp", dataHistProtonName.Data()));
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

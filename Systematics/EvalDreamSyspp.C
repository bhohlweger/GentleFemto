#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix) {
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle();
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);

  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix);

  //Proton - Proton
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(10, "ppDef", pp, ApAp);
  DreamSystematics protonproton(DreamSystematics::pp);
  protonproton.SetDefaultHist(CFppDef, "hCk_RebinnedppDefMeV_0");

  const int protonVarStart = 1;
  for (int i = protonVarStart;
      i < protonVarStart + protonproton.GetNumberOfVars(); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", i);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        10, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""), pp, ApAp);
    protonproton.SetVarHist(CFppVar,
                            TString::Format("Rebinned%sMeV", VarName.Data()));
  }
  protonproton.EvalSystematics();
  protonproton.EvalDifferenceInPairs();
  protonproton.WriteOutput();
  protonproton.DrawAllCF();
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2]);

  return 1;
}

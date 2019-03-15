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
const int rebin = 20;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");
  ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(), prefix,
                                                       "0");
  counter->SetNumberOfCandidates(ForgivingFile);
  const int nTracks = counter->GetNumberOfTracks();
  counter->ResetCounter();
  //Proton - Proton
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(rebin, "ppDef", pp, ApAp);
  DreamSystematics protonproton(DreamSystematics::pp);
  protonproton.SetDefaultHist(CFppDef, "hCk_ReweightedppDefMeV_1");

  const int protonVarStart = 1;
  for (int i = protonVarStart;
      i < protonVarStart + protonproton.GetNumberOfVars()-2; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", i);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    protonproton.SetVarHist(
        CFppVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    counter->SetNumberOfCandidates(ForgivingFile);
    protonproton.SetParticles(nTracks, 1, counter->GetNumberOfTracks(), 1);
  }

  for (int i = 43; i < 45; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", i);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    protonproton.SetVarHist(
        CFppVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    counter->SetNumberOfCandidates(ForgivingFile);
    protonproton.SetParticles(nTracks, 1, counter->GetNumberOfTracks(), 1);
  }
  protonproton.EvalSystematics();
  protonproton.EvalDifferenceInPairs();
  protonproton.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonproton.GetPairName().Data()), "update");
  CFppDef->WriteOutput(file, true);
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2]);

  return 1;
}

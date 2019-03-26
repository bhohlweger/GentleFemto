#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix, float upperFitRange) {
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
  const int rebin = 10;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename.Data(), prefix,
                                                       "0");
  counter->SetNumberOfCandidates(ForgivingFileDef);
  const int nTracks = counter->GetNumberOfTracks();
  counter->ResetCounter();
  //Proton - Proton
  CATSinput->SetFixedkStarMinBin(true, 0.004);
  DreamDist* pp = DreamFile->GetPairDistributions(0, 0, "");
  std::cout << "Femto Pairs pp: " << pp->GetFemtoPairs(0, 0.200) << std::endl;
  DreamDist* ApAp = DreamFile->GetPairDistributions(1, 1, "");
  std::cout << "Femto Pairs ApAp: " << ApAp->GetFemtoPairs(0, 0.200)
            << std::endl;
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(rebin, "ppDef", pp, ApAp);
  DreamSystematics protonproton(DreamSystematics::pp);
//  protonproton.SetUpperFitRange(44);
  protonproton.SetDefaultHist(CFppDef, "hCk_ReweightedppDefMeV_1");
  protonproton.SetUpperFitRange(upperFitRange);
  const int protonVarStart = 1;
  for (int i = protonVarStart;
      i < protonVarStart + protonproton.GetNumberOfVars() - 2; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", i);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    protonproton.SetVarHist(
        CFppVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonproton.SetParticles(nTracks, 1, counter->GetNumberOfTracks(), 1);
    counter->ResetCounter();
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
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonproton.SetParticles(nTracks, 1, counter->GetNumberOfTracks(), 1);
    counter->ResetCounter();
  }
  protonproton.EvalSystematics();
  protonproton.EvalDifferenceInParticles();
//  protonproton.EvalDifferenceInPairs();
  protonproton.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonproton.GetPairName().Data()), "update");
  CFppDef->WriteOutput(file, true);
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

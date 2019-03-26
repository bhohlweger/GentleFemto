#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix, int upperFitRange) {
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle();
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
  CATSinput->SetFixedkStarMinBin(true, 0.);
  const int rebin = 10;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

  //Proton - Proton
  DreamSystematics protonproton(DreamSystematics::pp);
  const int protonVarStart = 1;
  const int XiVarStart = 19;

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename.Data(), prefix,
                                                       "0");
  counter->SetNumberOfCandidates(ForgivingFileDef);
  const int nTracks = counter->GetNumberOfTracks();
  const int nCascades = counter->GetNumberOfCascades();
  counter->ResetCounter();
  //Proton - Xi
  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXiDef = CATSinput->ObtainCFSyst(rebin, "ppDef", pXi, ApAXi);
  DreamSystematics protonXi(DreamSystematics::pXi);
  protonXi.SetUpperFitRange(upperFitRange);
  protonXi.SetDefaultHist(CFpXiDef, "hCk_ReweightedppDefMeV_1");
  int iPXICounter = 0;
  for (int i = protonVarStart;
      i < protonVarStart + (protonproton.GetNumberOfVars()-2); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pXiVar%u", i);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""));
    protonXi.SetVarHist(CFpXiVar,
                        TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    std::cout << "nTracks: " << nTracks << '\t' << "nCascades: " << nCascades
              << '\t' << "counter->GetNumberOfTracks(): "
              << counter->GetNumberOfTracks() << '\t'
              << "counter->GetNumberOfCascades(): "
              << counter->GetNumberOfCascades() << std::endl;
    protonXi.SetParticles(nTracks, nCascades, counter->GetNumberOfTracks(),
                          counter->GetNumberOfCascades());
    counter->ResetCounter();
    iPXICounter++;
  }
  for (int i = XiVarStart;
      i
          < XiVarStart + protonXi.GetNumberOfVars()
              - (protonproton.GetNumberOfVars()-2); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pXiVar%u", i);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""));
    protonXi.SetVarHist(CFpXiVar,
                        TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    std::cout << "nTracks: " << nTracks << '\t' << "nCascades: " << nCascades
              << '\t' << "counter->GetNumberOfTracks(): "
              << counter->GetNumberOfTracks() << '\t'
              << "counter->GetNumberOfCascades(): "
              << counter->GetNumberOfCascades() << std::endl;
    protonXi.SetParticles(nTracks, nCascades, counter->GetNumberOfTracks(),
                          counter->GetNumberOfCascades());

    counter->ResetCounter();
    iPXICounter++;
  }
  protonXi.EvalSystematics();
//  protonXi.EvalDifferenceInPairs();
  protonXi.EvalDifferenceInParticles();
  protonXi.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonXi.GetPairName().Data()), "update");
  CFpXiDef->WriteOutput(file, true);
  std::cout << "Worked through " << iPXICounter << " variations" << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atoi(argv[3]));

  return 1;
}

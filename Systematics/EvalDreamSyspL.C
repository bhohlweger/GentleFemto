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
  CATSinput->SetFixedkStarMinBin(true, 0. );
  const int rebin = 10;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

  ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(), prefix,
                                                       "0");
  counter->SetNumberOfCandidates(ForgivingFile);
  const int nTracks = counter->GetNumberOfTracks();
  const int nv0 = counter->GetNumberOfV0s();
  counter->ResetCounter();
  //Proton - L
  DreamDist* pL = DreamFile->GetPairDistributions(0, 2, "");
  DreamDist* ApAL = DreamFile->GetPairDistributions(1, 3, "");
  DreamCF* CFpLDef = CATSinput->ObtainCFSyst(rebin, "pLDef", pL, ApAL);
  DreamSystematics protonL(DreamSystematics::pL);
  protonL.SetDefaultHist(CFpLDef, "hCk_ReweightedpLDefMeV_1");
  protonL.SetUpperFitRange(upperFitRange);
  int iPLCounter = 0;
  for (int i = 1; i <= protonL.GetNumberOfVars(); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pLVar%u", i);
    DreamCF* CFpLVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 2, ""),
        DreamVarFile->GetPairDistributions(1, 3, ""));
    protonL.SetVarHist(CFpLVar,
                        TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    std::cout << "nTracks: " << nTracks << '\t' << "nv0: " << nv0
        << '\t' << "counter->GetNumberOfTracks(): "
        << counter->GetNumberOfTracks() << '\t'
        << "counter->GetNumberOfV0s(): " << counter->GetNumberOfV0s()
        << std::endl;
    protonL.SetParticles(nTracks, nv0, counter->GetNumberOfTracks(),
                          counter->GetNumberOfCascades());
    counter->ResetCounter();
    iPLCounter++;
  }
  protonL.EvalSystematics();
//  protonXi.EvalDifferenceInPairs();
  protonL.EvalDifferenceInParticles();
  protonL.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonL.GetPairName().Data()), "update");
  CFpLDef->WriteOutput(file, true);
  std::cout << "Worked through " << iPLCounter << " variations" << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atoi(argv[3]));

  return 1;
}

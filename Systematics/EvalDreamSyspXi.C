#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix, float upperFitRange) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
  CATSinput->SetFixedkStarMinBin(true, 0.);
  const int rebin = 5;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

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
  const int pairCountsDefault = CFpXiDef->GetFemtoPairs(0, 0.2);
  DreamSystematics protonXi(DreamSystematics::pXi);
  protonXi.SetUpperFitRange(upperFitRange);
  protonXi.SetBarlowUpperRange(400);
  protonXi.SetDefaultHist(CFpXiDef, "hCk_ReweightedppDefMeV_1");
  int iPXICounter = 0;
  for (int i = 1;
      i <= protonXi.GetNumberOfVars(); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pXiVar%u", i);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""));
    int femtoPairVar= CFpXiVar->GetFemtoPairs(0, 0.2);
//    std::cout << "?femtoPairVar: " << femtoPairVar << std::endl;
//    std::cout << ?
    float relDiff = (femtoPairVar-pairCountsDefault)/(float)pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    protonXi.SetVarHist(CFpXiVar,
                        TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonXi.SetParticles(nTracks, nCascades, counter->GetNumberOfTracks(),
                          counter->GetNumberOfCascades());
    protonXi.SetPair(pairCountsDefault, CFpXiVar->GetFemtoPairs(0, 0.2));
    counter->ResetCounter();
    iPXICounter++;
  }
  protonXi.EvalSystematics();
  protonXi.EvalDifferenceInPairs();
  protonXi.EvalDifferenceInParticles();
  protonXi.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonXi.GetPairName().Data()), "update");
  CFpXiDef->WriteOutput(file, true);
  std::cout << "Worked through " << iPXICounter << " variations" << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

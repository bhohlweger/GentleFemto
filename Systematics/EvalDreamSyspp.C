#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TSystem.h"

void EvalDreamSystematics(TString InputDir, TString prefix,
                          float upperFitRange) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
  const int rebin = 5;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename.Data(),
                                                          prefix, "0");
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
  DreamCF* CFppDef = CATSinput->ObtainCFSyst(rebin, "ppVar0", pp, ApAp);
  const int pairCountsDefault = CFppDef->GetFemtoPairs(0, 0.2);
  DreamSystematics protonproton(DreamSystematics::pp);
//  protonproton.SetUpperFitRange(44);
  if (rebin != 1) {
    protonproton.SetDefaultHist(CFppDef, "hCk_ReweightedppVar0MeV_1");
  } else {
    protonproton.SetDefaultHist(CFppDef, "hCk_ReweightedppVar0MeV_0");
  }
  protonproton.SetUpperFitRange(upperFitRange);
  protonproton.SetBarlowUpperRange(400);
  int outCounter = 1;
  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", outCounter);
    DreamCF* CFppVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
        DreamVarFile->GetPairDistributions(1, 1, ""));
    DreamCF* CFppOut = CATSinput->ObtainCFSyst(
            rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 0, ""),
            DreamVarFile->GetPairDistributions(1, 1, ""));
    int femtoPairVar= CFppVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar-pairCountsDefault)/(float)pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    if (rebin != 1) {
      protonproton.SetVarHist(
          CFppVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    } else {
      protonproton.SetVarHist(
          CFppVar, TString::Format("Reweighted%sMeV_0", VarName.Data()));
    }
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonproton.SetPair(pairCountsDefault, CFppVar->GetFemtoPairs(0, 0.2));
    protonproton.SetParticles(nTracks, 1, counter->GetNumberOfTracks(), 1);
    CFppOut->WriteOutput(TString::Format("%s/CF_pp_Var%u.root",gSystem->pwd(),outCounter++).Data());
    counter->ResetCounter();
  }

  protonproton.EvalSystematics();
  protonproton.EvalDifferenceInParticles();
  protonproton.EvalDifferenceInPairs();
  protonproton.WriteOutput();
  auto file = new TFile(
      Form("Systematics_%s.root", protonproton.GetPairName().Data()), "update");
  CFppDef->WriteOutput(TString::Format("%s/CF_pp_Var0.root",gSystem->pwd()).Data());
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

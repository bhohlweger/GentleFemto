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
  CATSinput->SetNormalization(0.2, 0.4);
  CATSinput->SetFixedkStarMinBin(true, 0.);
  const int rebin = 3;
  auto counter = new CandidateCounter();
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);//$ particles covered in the BBar analysis
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");
  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename.Data(),
                                                          prefix, "0");
  counter->SetNumberOfCandidatesBBar(ForgivingFileDef);
  const int nTracks = counter->GetNumberOfTracks();
  printf("Number of Tracks = %i\n", nTracks );
  const int nAntiTracks = counter->GetNumberOfAntiTracks();
  printf("Number of AntiTracks = %i\n", nAntiTracks );
  counter->ResetCounter();
  DreamDist* pAp = DreamFile->GetPairDistributions(0, 1, "");
  DreamCF* CFpApDef = CATSinput->ObtainCFSystBBar(rebin, "pApVar0", pAp, nullptr);
  const int pairCountsDefault = CFpApDef->GetFemtoPairsBBar(0, 0.2);
  printf("pairCountsDefault = %.2i\n", pairCountsDefault);

  DreamSystematics protonAntiproton(DreamSystematics::pAp);
  protonAntiproton.SetUpperFitRange(upperFitRange);
  protonAntiproton.SetBarlowUpperRange(400);
  if (rebin != 1) {
    printf("rebin != 1 -- BE AWARE!!!  \n");
    protonAntiproton.SetDefaultHist(CFpApDef, "hCk_ReweightedpApVar0MeV_1");
//    protonAntiproton.SetDefaultHist(CFpApDef, "hCk_ReweightedpApVar0MeV_2");
  } else {
    protonAntiproton.SetDefaultHist(CFpApDef, "hCk_ReweightedpApVar0MeV_0");
//    protonAntiproton.SetDefaultHist(CFpApDef, "hCk_ReweightedpApVar0MeV_1");
  }

  int outCounter = 1;
  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pApVar%u", outCounter);
    DreamCF* CFpApVar = CATSinput->ObtainCFSystBBar(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 1, ""),nullptr);
    DreamCF* CFpApOut = CATSinput->ObtainCFSystBBar(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 1, ""),nullptr);
    int femtoPairVar = CFpApVar->GetFemtoPairsBBar(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;//exit the loop
    }
    if(TMath::Abs(relDiff)>=0.1 && TMath::Abs(relDiff)<=0.199)printf("--- BIG relDiff (%i) = %.3f ----\n",i,relDiff);

    if (rebin != 1) {
      protonAntiproton.SetVarHist(
          CFpApVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));

    } else {
      protonAntiproton.SetVarHist(
          CFpApVar, TString::Format("Reweighted%sMeV_0", VarName.Data()));
    }
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidatesBBar(ForgivingFile);
    protonAntiproton.SetPair(pairCountsDefault, CFpApVar->GetFemtoPairsBBar(0, 0.2));
    protonAntiproton.SetParticles(nTracks, nAntiTracks, counter->GetNumberOfTracks(),
                      counter->GetNumberOfAntiTracks());
    CFpApOut->WriteOutput(
        TString::Format("%s/CF_pAp_Var%u.root", gSystem->pwd(), outCounter++).Data());
    counter->ResetCounter();
  }
  protonAntiproton.EvalSystematicsBBar(0);
  protonAntiproton.EvalDifferenceInParticles();
  protonAntiproton.EvalDifferenceInPairs();
  protonAntiproton.WriteOutput();
  CFpApDef->WriteOutput(
      TString::Format("%s/CF_pAp_Var0.root", gSystem->pwd()).Data());

  std::cout << "Worked through " << outCounter << " variations" << std::endl;
  std::cout << "Upper fit range " << upperFitRange << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

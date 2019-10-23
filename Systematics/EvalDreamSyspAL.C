#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix, float upperFitRange) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.2, 0.4);
  CATSinput->SetFixedkStarMinBin(true, 0.);
  const int rebin = 5;//default = 4 has binning of 16 MeV
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");

  ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(), prefix,
                                                       "0");
  counter->SetNumberOfCandidatesBBar(ForgivingFile);
  const int nTracks = counter->GetNumberOfTracks();
  const int nv0 = counter->GetNumberOfV0s();
  const int nAntiTracks = counter->GetNumberOfAntiTracks();
  const int nAntiv0 = counter->GetNumberOfAntiV0s();
  counter->ResetCounter();
  //Proton - AntiL
  DreamDist* pAL = DreamFile->GetPairDistributions(0, 3, "");
  DreamDist* ApL = DreamFile->GetPairDistributions(1, 2, "");
  DreamCF* CFpALDef = CATSinput->ObtainCFSyst(rebin, "pALDef", pAL, ApL);
  const int pairCountsDefault = CFpALDef->GetFemtoPairs(0, 0.2);
  printf("pairCountsDefault = %.2i\n", pairCountsDefault);
  DreamSystematics protonAL(DreamSystematics::pAL);
  if (rebin != 1) {
	   printf("rebin != 1 -- BE AWARE!!!  \n");
	  protonAL.SetDefaultHist(CFpALDef, "hCk_ReweightedpALDefMeV_1");
  } else {
	  protonAL.SetDefaultHist(CFpALDef, "hCk_ReweightedpALDefMeV_0");
  }
  protonAL.SetUpperFitRange(upperFitRange);
  int iPLCounter = 0;
  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pALVar%u", i);
    DreamCF* CFpALVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 3, ""),
        DreamVarFile->GetPairDistributions(1, 2, ""));
    DreamCF* CFpALOut = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 3, ""),
        DreamVarFile->GetPairDistributions(1, 2, ""));
    int femtoPairVar= CFpALVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar-pairCountsDefault)/(float)pairCountsDefault;
//    printf("reldiff (%i)=%.2f\n",i,relDiff);
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    if(TMath::Abs(relDiff)>=0.1 && TMath::Abs(relDiff)<=0.199)printf("--- BIG relDiff (%i) = %.3f ----\n",i,relDiff);
    if (rebin != 1) {
    	protonAL.SetVarHist(
    			CFpALVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    } else {
    	protonAL.SetVarHist(
    			CFpALVar, TString::Format("Reweighted%sMeV_0", VarName.Data()));
    }

    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidatesBBar(ForgivingFile);
    protonAL.SetPair(pairCountsDefault, CFpALVar->GetFemtoPairs(0, 0.2));
    protonAL.SetParticles(nTracks, nAntiv0, counter->GetNumberOfTracks(),
                         counter->GetNumberOfAntiV0s());

    CFpALOut->WriteOutput(
        TString::Format("%s/CF_pAL_Var%u.root", gSystem->pwd(), iPLCounter++)
            .Data());
    counter->ResetCounter();
  }

  protonAL.EvalSystematicsBBar(0);
  protonAL.EvalDifferenceInPairs();
  protonAL.EvalDifferenceInParticles();
  protonAL.WriteOutput();
  CFpALDef->WriteOutput(
      TString::Format("%s/CF_pAL_Var0.root", gSystem->pwd()).Data());
  std::cout << "Worked through " << iPLCounter << " variations" << std::endl;
  std::cout << "Upper fit range " << upperFitRange << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

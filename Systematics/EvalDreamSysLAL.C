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
  double norm1 = 0.18;
  double norm2 = 0.28;
  CATSinput->SetNormalization(norm1, norm2);
  CATSinput->SetFixedkStarMinBin(true, 0.);//
  const int rebin = 5;//default = 5 has binning of 20 MeV
  auto counter = new CandidateCounter();
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);//$ particles covered in the BBar analysis
  DreamFile->SetAnalysisFile(filename.Data(), prefix, "0");
  ForgivingReader* ForgivingFileDef = new ForgivingReader(filename.Data(),
                                                          prefix, "0");
  counter->SetNumberOfCandidatesBBar(ForgivingFileDef);
  const int nv0s = counter->GetNumberOfV0s();
//  printf("Number of Lambdas = %i\n", nv0s );
  const int nAntiv0s = counter->GetNumberOfAntiV0s();
//  printf("Number of AntiTracks = %i\n", nAntiv0s );
  counter->ResetCounter();
  DreamDist* LAL = DreamFile->GetPairDistributions(2, 3, "");
  DreamCF* CFLALDef = CATSinput->ObtainCFSystBBar(
    rebin, "LALVar0", DreamFile->GetPairDistributions(2, 3, ""), nullptr);
  const int pairCountsDefault = CFLALDef->GetFemtoPairsBBar(0, 0.2);
  printf("pairCountsDefault = %i\n", pairCountsDefault);

  DreamSystematics LambdaAntiLambda(DreamSystematics::LAL);
  LambdaAntiLambda.SetUpperFitRange(upperFitRange);
  LambdaAntiLambda.SetBarlowUpperRange(400);
  if (rebin != 1) {
	  printf("rebin != 1 -- BE AWARE!!!  \n");
	  LambdaAntiLambda.SetDefaultHist(CFLALDef, "hCk_ReweightedLALVar0MeV_1");
  } else {
	  LambdaAntiLambda.SetDefaultHist(CFLALDef, "hCk_ReweightedLALVar0MeV_0");
  }
  int outCounter = 1;
  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(4, 4);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("LALVar%u", outCounter);
    DreamCF* CFLALVar = CATSinput->ObtainCFSystBBar(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(2, 3, ""),nullptr);
    DreamCF* CFLALOut = CATSinput->ObtainCFSystBBar(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(2, 3, ""),nullptr);
    int femtoPairVar = CFLALVar->GetFemtoPairsBBar(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
//    printf("--- relDiff = %.3f ----\n",relDiff);
    if (TMath::Abs(relDiff) > 0.2) {
      continue;//exit the loop
    }
    if(TMath::Abs(relDiff)>=0.1 && TMath::Abs(relDiff)<=0.199)printf("--- relDiff (%i) = %.3f ----\n",i,relDiff);

    if (rebin != 1) {
    	LambdaAntiLambda.SetVarHist(
          CFLALVar, TString::Format("Reweighted%sMeV_1", VarName.Data()));
    } else {
    	LambdaAntiLambda.SetVarHist(
          CFLALVar, TString::Format("Reweighted%sMeV_0", VarName.Data()));
    }
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidatesBBar(ForgivingFile);
    LambdaAntiLambda.SetPair(pairCountsDefault, CFLALVar->GetFemtoPairsBBar(0, 0.2));
    LambdaAntiLambda.SetParticles(nv0s, nAntiv0s, counter->GetNumberOfV0s(),
                      counter->GetNumberOfAntiV0s());
    CFLALOut->WriteOutput(
        TString::Format("%s/CF_LAL_Var%u.root", gSystem->pwd(), outCounter++)
            .Data());
    counter->ResetCounter();
  }

  LambdaAntiLambda.EvalSystematicsBBar(0);
  LambdaAntiLambda.EvalDifferenceInParticles();
  LambdaAntiLambda.EvalDifferenceInPairs();
  LambdaAntiLambda.WriteOutput();
  CFLALDef->WriteOutput(
      TString::Format("%s/CF_LAL_Var0.root", gSystem->pwd()).Data());
  std::cout << "Worked through " << outCounter << " variations" << std::endl;
  std::cout << "Upper fit range " << upperFitRange << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

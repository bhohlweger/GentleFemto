#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "SidebandSigma.h"
#include "TCanvas.h"
#include <iostream>

TH1F* GetCorrelation(TString filename, TString appendix, TString suffix,
                     TString graphName, const double normLower,
                     const double normUpper, const int rebin, double &nPairs) {
  TH1F *outputGraph;

  ReadDreamFile *DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename.Data(), appendix.Data(), suffix.Data());

  DreamCF *CF_pDminus = new DreamCF();
  DreamPair *pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair *apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  pDminus->ShiftForEmpty(pDminus->GetPair());
  apDplus->ShiftForEmpty(apDplus->GetPair());
  pDminus->FixShift(pDminus->GetPairShiftedEmpty(0),
                    apDplus->GetPairShiftedEmpty(0), apDplus->GetFirstBin());
  apDplus->FixShift(apDplus->GetPairShiftedEmpty(0),
                    pDminus->GetPairShiftedEmpty(0), pDminus->GetFirstBin());
  pDminus->Rebin(pDminus->GetPairFixShifted(0), rebin, true);
  apDplus->Rebin(apDplus->GetPairFixShifted(0), rebin, true);
  pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(0), 0.2, 0.9,
                              pDminus->GetPair());
  apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(0), 0.2, 0.9,
                              apDplus->GetPair());
  CF_pDminus->SetPairs(pDminus, apDplus);
  CF_pDminus->GetCorrelations();
  nPairs = CF_pDminus->GetFemtoPairs(0, 0.2);

  for (auto it : CF_pDminus->GetCorrelationFunctions()) {
    TString itName = it->GetName();
    if (graphName == itName) {
      std::cout << it->GetName() << std::endl;
      outputGraph = it;
    }
  }

  delete apDplus;
  delete pDminus;
  //delete CF_pDminus;
//  delete DreamFile;
  return outputGraph;
}

void EvalSystematics(TString InputDir, int signal) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int rebin = 20;

  DreamPlot::SetStyle(false, true);

  TString name, dirName;
  if (signal == 0) {
    name = "";
    dirName = "HM_CharmFemto_";
  } else if (signal == 1) {
    name = "_SBLeft";
    dirName = "HM_CharmFemto_SBLeft_";
  } else if (signal == 2) {
    name = "_SBRight";
    dirName = "HM_CharmFemto_SBRight_";
  }

  TString dataGrName = "hCk_ReweightedMeV_0";
  TString InputFileName = InputDir;
  InputFileName += "/AnalysisResults.root";
  TString systVar = "0";
  const double normLower = 1.5;
  const double normUpper = 2.;
  double nPairs;
  auto grCF = GetCorrelation(InputFileName, dirName, systVar,
                             dataGrName, normLower, normUpper, rebin, nPairs);
  double nPairsDefault = nPairs;

  DreamSystematics syst(DreamSystematics::pD);
  syst.SetDefaultHist(grCF);
  syst.SetUpperFitRange(1000);
  syst.SetBarlowUpperRange(500);
  syst.SetEstimator(DreamSystematics::Uniform);

  for (int i = 1; i <= 20; ++i) {
    auto grCF = GetCorrelation(InputFileName, dirName, Form("%i", i),
                               dataGrName, normLower, normUpper, rebin, nPairs);
    syst.SetVarHist(grCF);
    syst.SetPair(nPairsDefault,nPairs);

  }
  syst.EvalSystematics();
  syst.EvalDifferenceInPairs();
  syst.WriteOutput(name);
}

int main(int argc, char* argv[]) {
  EvalSystematics(argv[1], atoi(argv[2]));

  return 1;
}

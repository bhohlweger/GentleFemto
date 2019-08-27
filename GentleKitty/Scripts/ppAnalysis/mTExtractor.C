#include "VariationmTAnalysis.h"
#include "TROOT.h"
int main(int argc, char* argv[]) {
//  gROOT->ProcessLine("gErrorIgnoreLevel = 2001");
  const char* DataDir = argv[1];
  const char* FitDirResonance = argv[2];
  const char* FitDirGauss = argv[3];
  const int nMTBins = 7;
  TString avgmTFile = TString::Format("%s/AveragemT.root", DataDir);
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  if (!mTFile) {
    std::cout << "No mT File \n";
    return -1;
  }
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT_ppVar0");
  if (!avgmT) {
    std::cout << "no average mT file " << std::endl;
    return -1;
  }

  VariationmTAnalysis* analyser = new VariationmTAnalysis(2, 26, 81);
  analyser->SetHistName("hCk_RebinnedppVar");
  analyser->SetFileName("OutFileVarpp.root");
  analyser->SetmTAverage(avgmT);
  analyser->SetLegData("p-p #oplus #bar{p}-#bar{p}", "fpe");
  analyser->SetLegModel("Fit Gaussian + Resonance", "f", 2, 3252);
  analyser->SetLegModel("Fit Gaussian", "f", 2, 3225);
  //this implies a folder structure following mTBin_[1,2,3,....]
  for (int imt = 1; imt <= nMTBins; ++imt) {
    std::cout << "========================== \n";
    std::cout << "mTBin: " << imt << std::endl;
    std::cout << "========================== \n";
    TString mTDataDir = Form("%s/mTBin_%u", DataDir, imt);
    analyser->SetSystematic(mTDataDir.Data());
    TString mTFitDirResonance = Form("%s/mTBin_%u", FitDirResonance, imt);
    analyser->SetVariation(mTFitDirResonance.Data(), 0);
    TString mTFitDirGauss = Form("%s/mTBin_%u", FitDirGauss, imt);
    analyser->SetVariation(mTFitDirGauss.Data(), 1);
  }
  std::vector<float> mTBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5 };
  analyser->SetmTBins(mTBins);
  analyser->MakeCFPlotsSingleBand();
  analyser->StoreRadvsmT("RadppvsmTGaussReso.root",0);
  analyser->StoreRadvsmT("RadppvsmTGauss.root",1);
  return 1;
}

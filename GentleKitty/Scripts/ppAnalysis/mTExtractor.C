
#include "VariationmTAnalysis.h"
#include "TROOT.h"
int main(int argc, char* argv[])   {
//  gROOT->ProcessLine("gErrorIgnoreLevel = 2001");
  const char* DataDir = argv[1];
  const char* FitDir = argv[2];
  const int nMTBins = 7;
  TString avgmTFile = TString::Format("%s/AveragemT.root",DataDir);
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  if (!mTFile) {
    std::cout << "No mT File \n";
    return -1;
  }
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT_ppVar0");
  if(! avgmT) {
    std::cout << "no average mT file " << std::endl;
    return -1;
  }

  VariationmTAnalysis* analyser = new VariationmTAnalysis(1,26,81);
  analyser->SetHistName("hCk_RebinnedppVar");
  analyser->SetFileName("OutFileVarpp.root");
  analyser->SetmTAverage(avgmT);
  //this implies a folder structure following mTBin_[1,2,3,....]
  for (int imt = 1; imt <= nMTBins; ++imt) {
    TString mTDataDir = Form("%s/mTBin_%u",DataDir,imt);
    analyser->SetSystematic(mTDataDir.Data());
    TString mTFitDir = Form("%s/mTBin_%u",FitDir,imt);
    analyser->SetVariation(mTFitDir.Data(),0);
  }
  analyser->MakeCFPlotsPP();
  analyser->MakeRadPlotsPP();
  return 1;
}

#include "VariationmTAnalysis.h"
#include "TROOT.h"
int main(int argc, char* argv[])   {
//  gROOT->ProcessLine("gErrorIgnoreLevel = 2001");
  const char* DataDir = argv[1];
  const char* FitDir = argv[2];
  const int nMTBins = 2;
  TString avgmTFile = TString::Format("%s/AveragemT.root",DataDir);
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  if (!mTFile) {
    std::cout << "No mT File \n";
    return -1;
  }
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT_pLVar0");
  if(! avgmT) {
    std::cout << "no average mT file " << std::endl;
    return -1;
  }

  VariationmTAnalysis* analyser = new VariationmTAnalysis(3,42,54);
  analyser->SetHistName("hCk_RebinnedpLVar");
  analyser->SetFileName("OutFileVarpL.root");
  analyser->SetmTAverage(avgmT);
  //this implies a folder structure following mTBin_[1,2,3,....]
  for (int imt = 1; imt <= nMTBins; ++imt) {
    TString mTDataDir = Form("%s/mTBin_%u",DataDir,imt);
    analyser->SetSystematic(mTDataDir.Data());
    TString mTFitDirUsmani = Form("%s/mTBin_%u/1/",FitDir,imt);
    TString mTFitDirNLO = Form("%s/mTBin_%u/2/",FitDir,imt);
    TString mTFitDirLO = Form("%s/mTBin_%u/3/",FitDir,imt);
    analyser->SetVariation(mTFitDirUsmani.Data(),0);
    analyser->SetVariation(mTFitDirNLO.Data(),1);
    analyser->SetVariation(mTFitDirLO.Data(),2);
  }
  analyser->MakeCFPlotsPL();
  analyser->MakeRadPlotsPL();
  return 1;
}

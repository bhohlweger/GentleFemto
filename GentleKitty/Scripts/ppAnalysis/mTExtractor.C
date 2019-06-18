
#include "VariationmTAnalysis.h"

int main(int argc, char* argv[])   {
  const char* DataDir = argv[1];
  const char* FitDir = argv[2];
  const int nMTBins = 6;
  VariationmTAnalysis* analyser = new VariationmTAnalysis();
  analyser->SetHistName("hCk_RebinnedppVar");
  //this implies a folder structure following mTBin_[1,2,3,....]
  for (int imt = 1; imt <= nMTBins; ++imt) {
    TString mTDataDir = Form("%s/mTBin_%u",DataDir,imt);
    analyser->SetSystematic(mTDataDir.Data());
  }
  analyser->MakePlots();
  return 1;
}

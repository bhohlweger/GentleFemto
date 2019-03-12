#include "METoSEReweighting.C"

void METoSEReweightingSigma0(const char* foldername) {
  const char* filenames[4] = { "pp", "pSigma", "pSB_low", "pSB_up" };
  for (int iFile = 0; iFile < 4; ++iFile) {
    TString FileName = Form("%sCFOutput_%s.root", foldername, filenames[iFile]);
    std::cout << FileName.Data() << std::endl;
    TFile* file = TFile::Open(FileName, "update");
    TList* PairDist = (TList*) file->Get("PairDist");
    if (PairDist) {
      ReweightingQA(PairDist);
    } else {
      file->ls();
    }
    TList* AntiPairDist = (TList*) file->Get("AntiPairDist");
    if (AntiPairDist) {
      ReweightingQA(AntiPairDist);
    }
  }
}

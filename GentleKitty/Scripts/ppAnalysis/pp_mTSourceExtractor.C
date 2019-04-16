#include "TFile.h"
#include "TGraphErrors.h"

void mTRadiusExtractor(const char *inputFile, const char* mTFile) {
  TFile* inFile = TFile::Open(inputFile,"READ");
  TFile* inmTFile = TFile::Open(mTFile,"READ");
  TGraphErrors* AveragemT = (TGraphErrors*)inmTFile->Get("AveragemT");
}

int main(int argc, char *argv[]) {
  mTRadiusExtractor(argv[1],argv[2]);
  return 0;
}

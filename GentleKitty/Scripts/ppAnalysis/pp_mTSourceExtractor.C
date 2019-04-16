#include "TFile.h"

void mTRadiusExtractor(const char *inputFile) {
  TFile* inFile = TFile::Open(inputFile);

}

int main(int argc, char *argv[]) {
  mTRadiusExtractor(argv[1]);
  return 0;
}

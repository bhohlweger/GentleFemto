#include "pSigma0.C"
#include "DrawSigma.C"
#include "TROOT.h"

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  FitSigma0(argv);
  DrawSigma(argv);
  return 0;
}

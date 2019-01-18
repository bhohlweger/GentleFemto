#include "GetSigmaCorrelations.C"
#include "METoSEReweightingSigma0.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* trigger = argv[2];
  const char* suffix = argv[3];
  GetSigmaCorrelations(filename, trigger, suffix);

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  METoSEReweightingSigma0(foldername.Data());

  return 1;
}

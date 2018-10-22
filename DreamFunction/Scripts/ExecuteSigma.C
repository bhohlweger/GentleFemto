#include "GetSigmaCorrelations.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* suffix = argv[2];
  GetSigmaCorrelations(filename, suffix);

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
//  METoSEReweighting(foldername.Data());

  return 1;
}

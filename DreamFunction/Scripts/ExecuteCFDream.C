#include "GetCorrelations.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  GetCorrelations(filename, prefix, addon);
  GetCorrelationsBbarB(filename, prefix, addon);
  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  METoSEReweighting(foldername.Data());

  return 1;
}

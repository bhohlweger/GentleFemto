#include "GetCorrelations.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const float fixShift = atof(argv[1]);
  const char* filename = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";
  GetCorrelations(fixShift,filename, prefix, addon);
  GetCorrelationsBbarB(filename, prefix, addon);
  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  METoSEReweighting(foldername.Data());

  return 1;
}

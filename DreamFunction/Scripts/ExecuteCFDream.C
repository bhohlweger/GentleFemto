#include "GetCorrelations.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const float fixShift = atof(argv[1]);
  const char* filename = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";
  GetCorrelations(fixShift,filename, prefix, addon);
  METoSEReweighting(gSystem->pwd());

  return 1;
}

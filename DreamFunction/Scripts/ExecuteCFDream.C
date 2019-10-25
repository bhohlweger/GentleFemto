#include "GetCorrelations.C"
//#include "GetCorrelationsNano.C"
#include "METoSEReweighting.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
//  GetCorrelations(fixShift,filename, prefix, addon);
  GetCorrelations(filename, prefix, addon);
//  METoSEReweighting(gSystem->pwd());

  return 1;
}

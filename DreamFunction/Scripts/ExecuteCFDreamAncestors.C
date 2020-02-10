#include "GetCorrelationsNanoBBAncestors.C"
#include "GetCorrelationsNanoBB.C"
#include "METoSEReweighting.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  GetCorrelationsNanoBBAncestors(filename, prefix, addon, 0.24, 0.34, 0.18, 0.28);
  return 1;
}

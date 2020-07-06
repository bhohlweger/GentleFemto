#include "GetCorrelationsNano.C"
#include "GetCorrelationsNanoBBar.C"
#include "GetCorrelationsBbarB.C"
#include "GetCorrelationsNanoBB.C"
#include "GetCorrelationsNanoSample.C"
#include "METoSEReweightingBBar.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  //  GetCorrelationsNanoBBar(filename, prefix, addon, 0.24, 0.34);
   GetCorrelationsNanoBB(filename, prefix, addon, 0.24, 0.34);
    // GetCorrelationsNanoSample(filename, prefix, addon, 0.24, 0.34,0.18,0.28);
//  METoSEReweightingBBar(gSystem->pwd());

  return 1;
}

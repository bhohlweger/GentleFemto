#include "GetCorrelations.C"
#include "BbarB_QAplots.C"
#include "BbarB_PurityPlots.C"
#include "GetCorrelationsBbarB.C"
#include "GetCorrectBinCF.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filenamedata = argv[1];
  const char* filenamemc = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";
  bool epos = false;
  int decider = (atoi(argv[5])) ? atoi(argv[5]) : 0;
    if (decider > 0) {
       epos = true;
    } else {
       epos = false;
    }

  GetCorrectBinCF(filenamedata, filenamemc,prefix, addon,epos);


  return 1;
}

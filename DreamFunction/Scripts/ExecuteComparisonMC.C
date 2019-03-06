#include "GetCorrelations.C"
#include "BbarB_QAplots.C"
#include "BbarB_PurityPlots.C"
#include "GetCorrelationsBbarB.C"
#include "GetComparisonMC.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename1 = argv[1];
  const char* filename2 = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";

  GetComparisonMC(filename1,filename2,prefix,addon);


  return 1;
}

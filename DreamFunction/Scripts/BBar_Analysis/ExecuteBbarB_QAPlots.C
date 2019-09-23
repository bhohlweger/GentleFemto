#include "BbarB_QAplots.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = argv[3];
  const char* date = (argv[4]) ? argv[4] : "";
  bool isMC=false;

  BbarB_QAplots(filename,prefix,addon,date);

  return 1;
}

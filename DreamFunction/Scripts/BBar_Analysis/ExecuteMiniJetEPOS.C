#include "BbarB_QAplots.C"
#include "BbarB_PurityPlots.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"
#include "BbarB_MiniJetEPOS.C"

int main(int argc, char* argv[]) {
  const char* STflag = argv[1];
bool writefile=false;
// if you do not give an argument isMC is set to false
// atoi cast the integer to a char
int decider = (atoi(argv[2])) ? atoi(argv[2]) : 0;
if (decider > 0) {
   writefile = true;
} else {
   writefile = false;
}
// ./executeMiniJet #st 0 if you do not want to write the file
  BbarB_MiniJetEPOS(STflag,writefile);


  return 1;
}

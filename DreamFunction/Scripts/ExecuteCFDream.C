#include "GetCorrelations.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  bool isMC=false;
  bool epos=false;
// if you do not give an argument isMC is set to false
// atoi cast the integer to a char
int decider = (atoi(argv[4])) ? atoi(argv[4]) : 0;
  if (decider > 0) {
     isMC = true;
  } else {
     isMC = false;
  }
  int deciderep = (atoi(argv[5])) ? atoi(argv[5]) : 0;
    if (deciderep > 0) {
       epos = true;
    } else {
       epos = false;
    }
//  GetCorrelations(filename, prefix, addon);
  GetCorrelationsBbarB(filename, prefix, addon,isMC,epos);


  return 1;
}

#include "GetCorrelations.C"
#include "GetCorrelationsBbarB.C"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {

  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  bool iss0=false;
  bool isMC=false;
  bool epos=false;
// if you do not give an argument isMC is set to false
// atoi cast the integer to a char
int deciderS0 = (atoi(argv[4])) ? atoi(argv[4]) : 0;
  if (deciderS0 > 0) {
     iss0 = true;
  } else {
     iss0 = false;
  }
int decider = (atoi(argv[5])) ? atoi(argv[5]) : 0;
  if (decider > 0) {
     isMC = true;
  } else {
     isMC = false;
  }
  int deciderep = (atoi(argv[6])) ? atoi(argv[6]) : 0;
    if (deciderep > 0) {
       epos = true;
    } else {
       epos = false;
    }
  GetCorrelationsBbarB(filename, prefix, addon,iss0, isMC, epos);


  return 1;
}

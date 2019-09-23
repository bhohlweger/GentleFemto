#include "PlotsCorrelationBBar.C"

int main(int argc, char* argv[]) {
  const char* STflag = argv[1];
//  BbarB_MiniJet(filename,filenameMC);
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
  PlotsCorrelationBBar(STflag);
  return 1;
}

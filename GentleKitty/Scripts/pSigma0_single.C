#include "pSigma0.C"

/// =====================================================================================
int main(int argc, char *argv[]) {
  FitSigma0(atoi(argv[1]), argv[2], argv[3], argv[4], atoi(argv[5]),
            atof(argv[6]), atof(argv[7]), atof(argv[8]));
  return 0;
}

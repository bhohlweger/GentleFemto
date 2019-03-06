#include "pSigma0.C"
#include "DrawSigma.C"

/// =====================================================================================
int main(int argc, char *argv[]) {
  FitSigma0(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], atoi(argv[6]),
            atof(argv[7]), atof(argv[8]), atof(argv[9]));
  DrawSigma(atoi(argv[1]), argv[2], atoi(argv[6]), atof(argv[7]), atof(argv[8]),
            atof(argv[9]));
  return 0;
}

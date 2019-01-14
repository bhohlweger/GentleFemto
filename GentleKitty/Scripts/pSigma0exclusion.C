#include "pSigma0.C"
#include "DrawSigma.C"

/// =====================================================================================
/// \param path/to/data
/// \param Appendix
/// \param path/to/output
/// \param d0
/// \param f0inv
/// \param nSteps
/// \param StepSize

int main(int argc, char *argv[]) {
  const int allSyst = 1;
  const int scattPar = 1;

  const float d0def = atof(argv[4]);
  const float f0invdef = atof(argv[5]);
  const int nSteps = atoi(argv[6]);
  const float stepSize = atof(argv[7]);

  float d0current, f0invcurrent;
  // loop over d0
  for (float d0iter = 0; d0iter < nSteps; ++d0iter) {
    d0current = d0def + d0iter * stepSize;

    for (float f0inviter = 0; f0inviter < nSteps; ++f0inviter) {
      f0invcurrent = f0invdef + f0inviter * stepSize;

      FitSigma0(allSyst, argv[1], argv[2], argv[3], scattPar, d0current,
                f0invcurrent);
      DrawSigma(allSyst, argv[1], argv[2], argv[3], scattPar, d0current,
                f0invcurrent);
    }
  }
  return 0;
}

#include "CATSLambdaParam.h"

void LambdaPP() {
  double PurityProton = 0.9943;
  double protonPrim = 0.873;
  double protonSecLam = 0.089;
  double protonSecSig = 0.037;

  const Particle p1(
      PurityProton,
      protonPrim,
      { protonSecLam, protonSecSig});
  const Particle p2 = p1;
  //
  CATSLambdaParam params(p1, p2, true);
  std::cout << "LAMBDA PP \n";
  params.PrintLambdaParams();
}

void LambdaPXi() {
  double PurityProton = 0.9943;
  double protonPrim = 0.873;
  double protonSecLam = 0.089;
  double protonSecSig = 0.037;

  double PurityXi = 0.915;
  const double Xim1530_to_Xim = 0.32 * (1. / 3.);
  const double Xin1530_to_Xim = 0.32 * (2. / 3.);
  const double Omegam_to_Xim = 0.1;
  const double OmegamXim_BR = 0.086;

  const Particle proton(
      PurityProton,
      protonPrim,
      { protonSecLam, protonSecSig });
  const Particle xi(
      PurityXi,
      1. - Xim1530_to_Xim - Xin1530_to_Xim - Omegam_to_Xim * OmegamXim_BR, { {
          Xim1530_to_Xim, Xin1530_to_Xim, Omegam_to_Xim * OmegamXim_BR } });

  CATSLambdaParam params(proton, xi);
  std::cout << "LAMBDA PXi \n";
  params.PrintLambdaParams();
}

int main(int argc, char* argv[]) {
  LambdaPP();
  LambdaPXi();
  return 0;
}

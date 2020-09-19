#include "CATSLambdaParam.h"

void LambdaPP() {
  double PurityProton = 0.9943;
  double PrimProton = 0.823;
  double SecLamProton = 0.125;  //Fraction of Lambdas
  double SecSigmaProton = 1. - PrimProton - SecLamProton;

  const Particle p1(PurityProton, PrimProton, { SecLamProton, SecSigmaProton });
  const Particle p2 = p1;
  //
  CATSLambdaParam params(p1, p2, true);
  std::cout << "LAMBDA PP \n";
  params.PrintLambdaParams();
}

void LambdaPL() {
  double PurityProton = 0.9943;
  double PrimProton = 0.823;
  double SecLamProton = 0.125;  //Fraction of Lambdas
  double SecSigmaProton = 1. - PrimProton - SecLamProton;

  const Particle proton(PurityProton, PrimProton, { SecLamProton, SecSigmaProton });

  double PurityLambda = 0.96;
  double PrimLambda = 0.54/0.92;
  double SecSigLambda = 0.18/0.92;
  double SecXimLambda = 0.10/0.92;
  double SecXi0Lambda0 = 0.10/0.92;

  const Particle Lambda(PurityLambda, PrimLambda, { SecSigLambda, SecXimLambda, SecXi0Lambda0 });

  CATSLambdaParam params(proton, Lambda);
  std::cout << "LAMBDA PL \n";
  params.PrintLambdaParams();
}
void LambdaPXi() {
  double PurityProton = 0.9943;
  double PrimProton = 0.86;
  double SecLamProton = 0.124;  //Fraction of Lambdas
  double SecSigmaProton = 1. - PrimProton - SecLamProton;

  const Particle proton(PurityProton, PrimProton, { SecLamProton, SecSigmaProton });

  double PurityXi = 0.915;
  /*
  double PurityProton = 0.97;
  double PrimProton = 0.87;
  double SecLamProton = 0.09;  //Fraction of Lambdas
  double SecSigmaProton = 1. - PrimProton - SecLamProton;
  
  const Particle proton(PurityProton, PrimProton, { SecLamProton, SecSigmaProton });
 
  double PurityXi = 0.88;
  */ 
  // Xi01530 Production: dN/dy = 2.6e-3 (https://link.springer.com/content/pdf/10.1140%2Fepjc%2Fs10052-014-3191-x.pdf)
  // Xim1530 Production = Xi01530 Production
  // Xim Production: dN/dy = 7.9e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // -> Production Ratio ~ 1/3

  const double Xi01530XimProdFraction = 1 / 3.;
  const double Xim1530XimProdFraction = 1 / 3.;

  // 2/3 of Xi0(1530) decays via Xi- + pi+ (Isospin considerations)
  const double Xi01530Xim_BR = 2 / 3.;
  // 1/3 of Xi-(1530) decays via Xi- + pi0 (Isospin considerations)
  const double Xim1530Xim_BR = 1 / 3.;

  // Omega production: dN/dy = 0.67e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // -> Production Ratio ~ 1/10
  const double OmegamXimProdFraction = 1 / 9.;
  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %

  // Produce N Xi's -> Produce:
  // 1 ) N* 1/10 Omegas -> See N* 1/10 * 8.6% more Xi's
  // 2)  N* 1/2 Xi0_1530 -> See N*1/3*2/3 = N* 2/9 more Xi's
  // 3)  N* 1/2 Xim_1530 -> See N*1/3*1/3 = N* 1/9 more Xi's
  // Total Sample:  N(1+0.0086+2/9+1/9) ->
  // Primary Fraction = N / N(1+0.0086+2/9+1/9)
  // Secondary Omegas = N*0.0086  / N(1+0.0086+2/9+1/9)
  // etc.
  double XiNormalization = 1 + OmegamXimProdFraction * OmegamXim_BR
      + Xi01530XimProdFraction * Xi01530Xim_BR
      + Xim1530XimProdFraction * Xim1530Xim_BR;
  double SecOmegaXim = OmegamXimProdFraction * OmegamXim_BR
      / (double) XiNormalization;
  double SecXi01530Xim = Xi01530XimProdFraction * Xi01530Xim_BR
      / (double) XiNormalization;
  double SecXim1530Xim = Xim1530XimProdFraction * Xim1530Xim_BR
      / (double) XiNormalization;
  double PrimXim = 1. / (double) XiNormalization;
  const Particle xi(PurityXi, PrimXim, { SecOmegaXim, SecXi01530Xim,
                        SecXim1530Xim });
  CATSLambdaParam params(proton, xi);
  std::cout << "LAMBDA PXi \n";
  params.PrintLambdaParams();
}

int main(int argc, char* argv[]) {
  LambdaPP();
  LambdaPL();
  LambdaPXi();
  return 0;
}

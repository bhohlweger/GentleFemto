#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "TDatabasePDG.h"
#include <iostream>
void testCats() {
  const unsigned momBins = 105;
  const double kMin = 0;
  const double kMax = kMin + 4 * momBins;  //(4 is the bin width)
  const double GaussSourceSize = 1.2;
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  double *pars = new double[6];
  pars[0] = 0;
  pars[1] = 0;
  pars[2] = 0;
  pars[3] = GaussSourceSize * 1.2;
  pars[4] = GaussSourceSize / 1.2;
  pars[5] = 0.5;
//  double pars[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                 GaussSourceSize / 1.2, 0.5 };

  CATS AB_pp;

  const double Weight1S0 = 3. / 12.;
  const double Weight3P0 = 1. / 12.;
  const double Weight3P1 = 3. / 12.;
  const double Weight3P2 = 5. / 12.;
//  double *PotPars1S0 = new double[10];
//  PotPars1S0[0]= 0;
//  PotPars1S0[1]= 0;
//  PotPars1S0[2]= NN_AV18;
//  PotPars1S0[3]= v18_Coupled3P2;
//  PotPars1S0[4]= 1;
//  PotPars1S0[5]= 1;
//  PotPars1S0[6]= 1;
//  PotPars1S0[7]= 0;
//  PotPars1S0[8]= 0;
//  PotPars1S0[9]= 0;
//
//  double *PotPars3P0 = new double[10];
//  PotPars3P0[0]= 0;
//  PotPars3P0[1]= 0;
//  PotPars3P0[2]= NN_AV18;
//  PotPars3P0[3]= v18_Coupled3P2;
//  PotPars3P0[4]= 1;
//  PotPars3P0[5]= 1;
//  PotPars3P0[6]= 1;
//  PotPars3P0[7]= 1;
//  PotPars3P0[8]= 1;
//  PotPars3P0[9]= 0;
//
//  double *PotPars3P1 = new double[10];
//  PotPars3P1[0]= 0;
//  PotPars3P1[1]= 0;
//  PotPars3P1[2]= NN_AV18;
//  PotPars3P1[3]= v18_Coupled3P2;
//  PotPars3P1[4]= 1;
//  PotPars3P1[5]= 1;
//  PotPars3P1[6]= 1;
//  PotPars3P1[7]= 1;
//  PotPars3P1[8]= 1;
//  PotPars3P1[9]= 1;
//
//  double *PotPars3P2 = new double[10];
//  PotPars3P2[0]= 0;
//  PotPars3P2[1]= 0;
//  PotPars3P2[2]= NN_AV18;
//  PotPars3P2[3]= v18_Coupled3P2;
//  PotPars3P2[4]= 1;
//  PotPars3P2[5]= 1;
//  PotPars3P2[6]= 1;
//  PotPars3P2[7]= 1;
//  PotPars3P2[8]= 1;
//  PotPars3P2[9]= 2;
//
//

  double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0,
      0 };
  double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
      0 };
  double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
      1 };
  double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
      2 };
  AB_pp.SetAnaSource(GaussSource, pars);
  AB_pp.SetUseAnalyticSource(true);
  AB_pp.SetThetaDependentSource(false);
  AB_pp.SetExcludeFailedBins(false);
  AB_pp.SetMomBins(momBins, kMin, kMax);
  AB_pp.SetQ1Q2(1);
  AB_pp.SetPdgId(2212, 2212);
  AB_pp.SetRedMass(0.5 * Mass_p);
  AB_pp.SetNumChannels(4);
  AB_pp.SetNumPW(0, 2);
  AB_pp.SetNumPW(1, 2);
  AB_pp.SetNumPW(2, 2);
  AB_pp.SetNumPW(3, 2);
  AB_pp.SetSpin(0, 0);
  AB_pp.SetSpin(1, 1);
  AB_pp.SetSpin(2, 1);
  AB_pp.SetSpin(3, 1);
  AB_pp.SetChannelWeight(0, Weight1S0);
  AB_pp.SetChannelWeight(1, Weight3P0);
  AB_pp.SetChannelWeight(2, Weight3P1);
  AB_pp.SetChannelWeight(3, Weight3P2);
  AB_pp.SetShortRangePotential(0, 0, fDlmPot, PotPars1S0);
  AB_pp.SetShortRangePotential(1, 1, fDlmPot, PotPars3P0);
  AB_pp.SetShortRangePotential(2, 1, fDlmPot, PotPars3P1);
  AB_pp.SetShortRangePotential(3, 1, fDlmPot, PotPars3P2);
  AB_pp.KillTheCat();
  return;
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}

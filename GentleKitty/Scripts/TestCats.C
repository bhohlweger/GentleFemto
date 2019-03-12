#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
void testCats() {
  //  const double Weight1S0 = 3. / 12.;
  //  const double Weight3P0 = 1. / 12.;
  //  const double Weight3P1 = 3. / 12.;
  //  const double Weight3P2 = 5. / 12.;
  //
  //  double PotPars1S0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  //  double PotPars3P0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  //  double PotPars3P1[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
  //  double PotPars3P2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };
  //
  //  CATSparameters *cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,
  //                                                   8, true);
  //  cPotPars1S0->SetParameters(PotPars1S0);
  //  CATSparameters *cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,
  //                                                   8, true);
  //  cPotPars3P0->SetParameters(PotPars3P0);
  //  CATSparameters *cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,
  //                                                   8, true);
  //  cPotPars3P1->SetParameters(PotPars3P1);
  //  CATSparameters *cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,
  //                                                   8, true);
  //  cPotPars3P2->SetParameters(PotPars3P2);

  unsigned int momBins = 50;
  const double kMin = 0;
  const double kMax = 200;

  CATS AB_pp;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, 1.2);
  AB_pp.SetAnaSource(GaussSource, *cPars);
  AB_pp.SetUseAnalyticSource(true);
  AB_pp.SetMomentumDependentSource(true);
  AB_pp.SetThetaDependentSource(false);
  AB_pp.SetExcludeFailedBins(false);
  AB_pp.SetMomBins(momBins, kMin, kMax);
  AB_pp.SetQ1Q2(1);
  AB_pp.SetPdgId(2212, 2212);
  AB_pp.SetNumChannels(1);
  AB_pp.SetNumPW(0, 1);
  AB_pp.SetSpin(0, 0);
  //AB_pp.SetNumChannels(4);
  //  AB_pp->SetNumPW(0, 2);
  //  AB_pp->SetNumPW(1, 2);
  //  AB_pp->SetNumPW(2, 2);
  //  AB_pp->SetNumPW(3, 2);
  //  AB_pp->SetSpin(0, 0);
  //  AB_pp->SetSpin(1, 1);
  //  AB_pp->SetSpin(2, 1);
  //  AB_pp->SetSpin(3, 1);
  //  AB_pp->SetChannelWeight(0, Weight1S0);
  //  AB_pp->SetChannelWeight(1, Weight3P0);
  //  AB_pp->SetChannelWeight(2, Weight3P1);
  //  AB_pp->SetChannelWeight(3, Weight3P2);
  //  AB_pp->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  //  AB_pp->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  //  AB_pp->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  //  AB_pp->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
  AB_pp.KillTheCat();
  return;
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}

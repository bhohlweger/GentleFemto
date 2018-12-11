#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "TDatabasePDG.h"
#include <iostream>
void testCats() {
  const unsigned NumMomBins = 105;
  const double kMin = 0;
  const double kMax = kMin + 4 * NumMomBins;  //(4 is the bin width)
  const double GaussSourceSize = 1.2;
  double Pars_pXim1530[6] = { 0, 0, 0, GaussSourceSize * 1.2,
                 GaussSourceSize / 1.2, 0.5 };
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim1530 =
      TDatabasePDG::Instance()->GetParticle(3314)->Mass() * 1000;
  CATS AB_pXim1530;
  AB_pXim1530.SetAnaSource(GaussSource, Pars_pXim1530);
  AB_pXim1530.SetUseAnalyticSource(true);
  AB_pXim1530.SetThetaDependentSource(false);

  AB_pXim1530.SetExcludeFailedBins(false);
  AB_pXim1530.SetMomBins(NumMomBins, kMin, kMax);

  AB_pXim1530.SetNumChannels(1);
  AB_pXim1530.SetNumPW(0, 1);
  AB_pXim1530.SetSpin(0, 0);
  AB_pXim1530.SetChannelWeight(0, 1.);

  AB_pXim1530.SetQ1Q2(-1);
  AB_pXim1530.SetPdgId(2212, 3122);

  AB_pXim1530.SetRedMass((Mass_p * Mass_Xim1530) / (Mass_p + Mass_Xim1530));
  AB_pXim1530.KillTheCat();
  return;
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}

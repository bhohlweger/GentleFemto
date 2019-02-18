/*
 * TidyCats.cxx
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */

#include <TidyCats.h>
#include <iostream>
#include "TDatabasePDG.h"
#include "CATStools.h"
TidyCats::TidyCats() {
}

TidyCats::~TidyCats() {
  // TODO Auto-generated destructor stub
}

void TidyCats::GetCatsProtonProton(CATS* AB_pp, int momBins, double kMin,
                                   double kMax, bool ResonanceSource) {
  const double Weight1S0 = 3. / 12.;
  const double Weight3P0 = 1. / 12.;
  const double Weight3P1 = 3. / 12.;
  const double Weight3P2 = 5. / 12.;

  double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
  double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };

  CATSparameters cPotPars1S0(CATSparameters::tPotential, 8, true);
  cPotPars1S0.SetParameters(&PotPars1S0[2]);
  CATSparameters cPotPars3P0(CATSparameters::tPotential, 8, true);
  cPotPars3P0.SetParameters(&PotPars3P0[2]);
  CATSparameters cPotPars3P1(CATSparameters::tPotential, 8, true);
  cPotPars3P1.SetParameters(&PotPars3P1[2]);
  CATSparameters cPotPars3P2(CATSparameters::tPotential, 8, true);
  cPotPars3P2.SetParameters(&PotPars3P2[2]);

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  if (ResonanceSource) {
    double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    CATSparameters cPars(CATSparameters::tSource, 6, true);
    double Pars_pp[6] = { 1.4, 1.65, 0.3578, 1361.52, massProton, massPion };
    cPars.SetParameters(Pars_pp);
    AB_pp->SetAnaSource(GaussExpTotIdenticalSimple_2body, cPars);
  } else {
    CATSparameters cPars(CATSparameters::tSource, 1, true);
    cPars.SetParameter(0, 1.2);
    AB_pp->SetAnaSource(GaussSource, cPars);
  }
  AB_pp->SetUseAnalyticSource(true);
  AB_pp->SetThetaDependentSource(false);
  AB_pp->SetExcludeFailedBins(false);
  AB_pp->SetMomBins(momBins, kMin, kMax);
  AB_pp->SetQ1Q2(1);
  AB_pp->SetPdgId(2212, 2212);
  AB_pp->SetRedMass(0.5 * massProton);
  AB_pp->SetNumChannels(4);
  AB_pp->SetNumPW(0, 2);
  AB_pp->SetNumPW(1, 2);
  AB_pp->SetNumPW(2, 2);
  AB_pp->SetNumPW(3, 2);
  AB_pp->SetSpin(0, 0);
  AB_pp->SetSpin(1, 1);
  AB_pp->SetSpin(2, 1);
  AB_pp->SetSpin(3, 1);
  AB_pp->SetChannelWeight(0, Weight1S0);
  AB_pp->SetChannelWeight(1, Weight3P0);
  AB_pp->SetChannelWeight(2, Weight3P1);
  AB_pp->SetChannelWeight(3, Weight3P2);
  AB_pp->SetShortRangePotential(0, 0, fDlmPot, cPotPars1S0);
  AB_pp->SetShortRangePotential(1, 1, fDlmPot, cPotPars3P0);
  AB_pp->SetShortRangePotential(2, 1, fDlmPot, cPotPars3P1);
  AB_pp->SetShortRangePotential(3, 1, fDlmPot, cPotPars3P2);
  return;
}

void TidyCats::GetCatsProtonLambda(CATS* AB_pL, int momBins, double kMin,
                                   double kMax, bool ResonanceSource) {

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass()
      * 1000;
  if (ResonanceSource) {
    double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    CATSparameters cPars(CATSparameters::tSource, 11, true);
    double Pars_pL[11] = { 1.4, 1.65, 0.3578, 1361.52, massProton, massPion,
        4.69, 0.3562, 1462.93, massLambda, massPion };
    cPars.SetParameters(Pars_pL, true);
    AB_pL->SetAnaSource(GaussExpTotSimple_2body, cPars);
  } else {
    CATSparameters cPars(CATSparameters::tSource, 1, true);
    cPars.SetParameter(0, 1.2);
    AB_pL->SetAnaSource(GaussSource, cPars);
  }
  AB_pL->SetUseAnalyticSource(true);
  AB_pL->SetThetaDependentSource(false);
  AB_pL->SetMomBins(momBins, kMin, kMax);
  AB_pL->SetNumChannels(2);
  AB_pL->SetNumPW(0, 1);
  AB_pL->SetNumPW(1, 1);
  AB_pL->SetSpin(0, 0);
  AB_pL->SetSpin(1, 1);
  AB_pL->SetChannelWeight(0, 0.25);
  AB_pL->SetChannelWeight(1, 0.75);
  AB_pL->SetQ1Q2(0);
  AB_pL->SetPdgId(2212, 3122);
  AB_pL->SetRedMass((massProton * massLambda) / (massProton + massLambda));
  //!TEMPORARY: WE USE USMANI, ASK DIMI ABOUT FURTHER CHANGES (LIKE USE THE NLO)
  //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
  double pLamPotPars1S0[8] = { pL_UsmaniOli, 0, 0, 0, 0, 0, 0, 0 };
  double pLamPotPars3S1[8] = { pL_UsmaniOli, 0, 0, 0, 0, 1, 0, 1 };
  CATSparameters cPotPars1S0(CATSparameters::tPotential, 8, true);
  cPotPars1S0.SetParameters(pLamPotPars1S0);
  CATSparameters cPotPars3S1(CATSparameters::tPotential, 8, true);
  cPotPars3S1.SetParameters(pLamPotPars3S1);
  AB_pL->SetShortRangePotential(0, 0, fDlmPot, cPotPars1S0);
  AB_pL->SetShortRangePotential(1, 0, fDlmPot, cPotPars3S1);

  return;
}

void TidyCats::GetCatsProtonXiMinus(CATS* AB_pXim, int momBins, double kMin,
                                    double kMax, bool StrongOn,
                                    double QCDTime) {

  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, 1.2);
  AB_pXim->SetAnaSource(GaussSource, cPars);
  AB_pXim->SetUseAnalyticSource(true);
  AB_pXim->SetThetaDependentSource(false);

  AB_pXim->SetExcludeFailedBins(false);
  AB_pXim->SetMomBins(momBins, kMin, kMax);

  AB_pXim->SetNumChannels(4);
  AB_pXim->SetNumPW(0, 1);
  AB_pXim->SetNumPW(1, 1);
  AB_pXim->SetNumPW(2, 1);
  AB_pXim->SetNumPW(3, 1);
  AB_pXim->SetSpin(0, 0);    //I=0; S=0
  AB_pXim->SetSpin(1, 1);    //I=0; S=1
  AB_pXim->SetSpin(2, 0);    //I=1; S=0
  AB_pXim->SetSpin(3, 1);    //I=1; S=1
  AB_pXim->SetChannelWeight(0, 1. / 8.);
  AB_pXim->SetChannelWeight(1, 3. / 8.);
  AB_pXim->SetChannelWeight(2, 1. / 8.);
  AB_pXim->SetChannelWeight(3, 3. / 8.);
  AB_pXim->SetQ1Q2(-1);
  AB_pXim->SetPdgId(2212, 3122);  //same as Lambda, in case we want to use EPOS pL source

  AB_pXim->SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim->SetMaxRad(64);
  AB_pXim->SetMaxRho(32);
  if (StrongOn) {
    double pXimPotParsI0S0[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1 };  // 9, 10, 11, 12
    double pXimPotParsI1S0[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0 };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1 };  // this value 1-6
    CATSparameters cPotParsI0S0(CATSparameters::tPotential, 8, true);
    cPotParsI0S0.SetParameters(pXimPotParsI0S0);

    CATSparameters cPotParsI0S1(CATSparameters::tPotential, 8, true);
    cPotParsI0S1.SetParameters(pXimPotParsI0S1);

    CATSparameters cPotParsI1S0(CATSparameters::tPotential, 8, true);
    cPotParsI1S0.SetParameters(pXimPotParsI1S0);

    CATSparameters cPotParsI1S1(CATSparameters::tPotential, 8, true);
    cPotParsI1S1.SetParameters(pXimPotParsI1S1);
    AB_pXim->SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, fDlmPot, cPotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, fDlmPot, cPotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, fDlmPot, cPotParsI1S1);
  }

  return;
}

void TidyCats::GetCatsProtonXiMinusCutOff(CATS* AB_pXim, int momBins,
                                          double kMin, double kMax,
                                          bool StrongOn, double QCDTime,
                                          double cutOff) {

  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;

  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, 1.2);
  AB_pXim->SetAnaSource(GaussSource, cPars);
  AB_pXim->SetUseAnalyticSource(true);
  AB_pXim->SetThetaDependentSource(false);

  AB_pXim->SetExcludeFailedBins(false);
  AB_pXim->SetMomBins(momBins, kMin, kMax);

  AB_pXim->SetNumChannels(4);
  AB_pXim->SetNumPW(0, 1);
  AB_pXim->SetNumPW(1, 1);
  AB_pXim->SetNumPW(2, 1);
  AB_pXim->SetNumPW(3, 1);
  AB_pXim->SetSpin(0, 0);    //I=0; S=0
  AB_pXim->SetSpin(1, 1);    //I=0; S=1
  AB_pXim->SetSpin(2, 0);    //I=1; S=0
  AB_pXim->SetSpin(3, 1);    //I=1; S=1
  AB_pXim->SetChannelWeight(0, 1. / 8.);
  AB_pXim->SetChannelWeight(1, 3. / 8.);
  AB_pXim->SetChannelWeight(2, 1. / 8.);
  AB_pXim->SetChannelWeight(3, 3. / 8.);

  AB_pXim->SetQ1Q2(-1);

  AB_pXim->SetPdgId(2212, 3122);  //same as Lambda, in case we want to use EPOS pL source

  AB_pXim->SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim->SetMaxRad(64);
  AB_pXim->SetMaxRho(32);
//
  if (StrongOn) {
    double pXimPotParsI0S0[9] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0,
        cutOff };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[9] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1,
        cutOff };  // 9, 10, 11, 12
    double pXimPotParsI1S0[9] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0,
        cutOff };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[9] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1,
        cutOff };  // this value 1-6
    CATSparameters cPotParsI0S0(CATSparameters::tPotential, 9, true);
    cPotParsI0S0.SetParameters(pXimPotParsI0S0);

    CATSparameters cPotParsI0S1(CATSparameters::tPotential, 9, true);
    cPotParsI0S1.SetParameters(pXimPotParsI0S1);

    CATSparameters cPotParsI1S0(CATSparameters::tPotential, 9, true);
    cPotParsI1S0.SetParameters(pXimPotParsI1S0);

    CATSparameters cPotParsI1S1(CATSparameters::tPotential, 9, true);
    cPotParsI1S1.SetParameters(pXimPotParsI1S1);
    AB_pXim->SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, fDlmPot, cPotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, fDlmPot, cPotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, fDlmPot, cPotParsI1S1);
  }

  return;
}

void TidyCats::GetCatsProtonXiMinus1530(CATS* AB_pXim1530, int momBins,
                                        double kMin, double kMax) {
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim1530 =
      TDatabasePDG::Instance()->GetParticle(3314)->Mass() * 1000;
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, 1.2);
  AB_pXim1530->SetAnaSource(GaussSource, cPars);
  AB_pXim1530->SetUseAnalyticSource(true);
  AB_pXim1530->SetThetaDependentSource(false);

  AB_pXim1530->SetExcludeFailedBins(false);
  AB_pXim1530->SetMomBins(momBins, kMin, kMax);

  AB_pXim1530->SetNumChannels(1);
  AB_pXim1530->SetNumPW(0, 1);
  AB_pXim1530->SetSpin(0, 0);
  AB_pXim1530->SetChannelWeight(0, 1.);

  AB_pXim1530->SetQ1Q2(-1);
  AB_pXim1530->SetPdgId(2212, 3122);

  AB_pXim1530->SetRedMass((Mass_p * Mass_Xim1530) / (Mass_p + Mass_Xim1530));

  return;
}

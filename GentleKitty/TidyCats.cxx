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
#include "DLM_Source.h"
#include "DLM_WfModel.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"

static double DegToRad = 0.01745;

TidyCats::TidyCats()
    : fppCleverLevy(nullptr),
      fppCleverMcLevy(nullptr),
      fpLCleverLevy(nullptr),
      fpLCleverMcLevy(nullptr),
      fpXimCleverLevy(nullptr),
      fpXimCleverMcLevy(nullptr),
      fpXim1530CleverLevy(nullptr),
      fpSigma0CleverMcLevy(nullptr) {
}

TidyCats::~TidyCats() {
  if (fppCleverLevy) {
    delete fppCleverLevy;
  }
  if (fppCleverMcLevy) {
    delete fppCleverMcLevy;
  }
  if (fpLCleverLevy) {
    delete fpLCleverLevy;
  }
  if (fpLCleverMcLevy) {
    delete fpLCleverMcLevy;
  }
  if (fpXimCleverLevy) {
    delete fpXimCleverLevy;
  }
  if (fpXimCleverMcLevy) {
    delete fpXimCleverMcLevy;
  }
  if (fpXim1530CleverLevy) {
    delete fpXim1530CleverLevy;
  }
  if (fpSigma0CleverMcLevy) {
    delete fpSigma0CleverMcLevy;
  }
}

void TidyCats::GetCatsProtonProton(CATS* AB_pp, int momBins, double kMin,
                                   double kMax, TidyCats::Sources source) {
  TString HomeDir = gSystem->GetHomeDirectory().c_str();
  const double Weight1S0 = 3. / 12.;  // also the weight for 1D2
  const double Weight3P0 = 1. / 12.;
  const double Weight3P1 = 3. / 12.;
  const double Weight3P2 = 5. / 12.;

  double PotPars1S0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  double PotPars1D2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 2, 2 };  //the last 3 digits are s,l,j
  double PotPars3P0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  double PotPars3P1[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
//  double PotPars3P2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };
  double PotPars3P2[8] = { NN_AV18, v18_SingleChannelMagic, 1, 1, 1, 1, 1, 2 };  // magically accounts for the coupling to F something

  CATSparameters *cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars1S0->SetParameters(PotPars1S0);
  CATSparameters *cPotPars1D2 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars1D2->SetParameters(PotPars1D2);
  CATSparameters *cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P0->SetParameters(PotPars3P0);
  CATSparameters *cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P1->SetParameters(PotPars3P1);
  CATSparameters *cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P2->SetParameters(PotPars3P2);

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
  double Pars_pp[6] = { 1.4, 1.65, 0.3578, 1361.52, massProton, massPion };
  CATSparameters* cPars;

  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      AB_pp->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance:
      fppCleverMcLevy = new DLM_CleverMcLevyReso();
      fppCleverMcLevy->InitNumMcIter(1000000);
      fppCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fppCleverMcLevy->InitScale(100, 0.2, 2.6);
      fppCleverMcLevy->InitRad(512, 0, 64);
      fppCleverMcLevy->InitType(2);
      fppCleverMcLevy->InitReso(0, 1);  //number of p resonances
      fppCleverMcLevy->InitReso(1, 1);  //number of Xi resonances
      AB_pp->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
      AB_pp->SetAnaSource(0, 1.2);
      AB_pp->SetAnaSource(1, 2.0);
      break;
    case TidyCats::sLevy:
      fppCleverMcLevy = new DLM_CleverMcLevyReso();
      fppCleverMcLevy->InitNumMcIter(1000000);
      fppCleverMcLevy->InitStability(20, 1, 2);
      fppCleverMcLevy->InitScale(100, 0.2, 2.6);
      fppCleverMcLevy->InitRad(512, 0, 64);
      fppCleverMcLevy->InitType(2);
      fppCleverMcLevy->InitReso(0, 1);  //number of p resonances
      fppCleverMcLevy->InitReso(1, 1);  //number of Xi resonances
      fppCleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
                                 massPion);
      fppCleverMcLevy->SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
                                 massPion);
      //Nolan Parameterization.
      AB_pp->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
      AB_pp->SetAnaSource(0, 1.2);  //r0
      AB_pp->SetAnaSource(1, 1.6);  //Stability alpha ( 1= Cauchy, ... 2 = Gauss)
      break;
    default:
      std::cout << "Source not implemented \n";
      break;
  }
  AB_pp->SetEpsilonConv(1e-7);
  AB_pp->SetEpsilonProp(1e-7);
  AB_pp->SetUseAnalyticSource(true);
  AB_pp->SetMomentumDependentSource(false);
  AB_pp->SetThetaDependentSource(false);
  AB_pp->SetExcludeFailedBins(false);
  AB_pp->SetMomBins(momBins, kMin, kMax);
  AB_pp->SetQ1Q2(1);
  AB_pp->SetPdgId(2212, 2212);
  AB_pp->SetRedMass(0.5 * massProton);
  AB_pp->SetNumChannels(4);
  AB_pp->SetNumPW(0, 3);  // the maximum number of partial waves flying around = the number of l (s,p,d -> l = 3)
  AB_pp->SetNumPW(1, 3);
  AB_pp->SetNumPW(2, 3);
  AB_pp->SetNumPW(3, 3);
  AB_pp->SetSpin(0, 0);
  AB_pp->SetSpin(1, 1);
  AB_pp->SetSpin(2, 1);
  AB_pp->SetSpin(3, 1);
  AB_pp->SetChannelWeight(0, Weight1S0);
  AB_pp->SetChannelWeight(1, Weight3P0);
  AB_pp->SetChannelWeight(2, Weight3P1);
  AB_pp->SetChannelWeight(3, Weight3P2);
  AB_pp->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp->SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
  AB_pp->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  AB_pp->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  AB_pp->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);

  return;
}

void TidyCats::GetCatsProtonLambda(CATS* AB_pL, int momBins, double kMin,
                                   double kMax, TidyCats::Sources source,
                                   TidyCats::pLPot pot) {

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
  double Pars_pL[11] = { 1.4, 1.65, 0.3578, 1361.52, massProton, massPion, 4.69,
      0.3562, 1462.93, massLambda, massPion };
  CATSparameters* cPars;
  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      AB_pL->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance:
      fpLCleverMcLevy = new DLM_CleverMcLevyReso();
      fpLCleverMcLevy->InitNumMcIter(1000000);
      fpLCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpLCleverMcLevy->InitScale(100, 0.2, 2.6);
      fpLCleverMcLevy->InitRad(512, 0, 64);
      fpLCleverMcLevy->InitType(2);
      fpLCleverMcLevy->InitReso(0, 1);  //number of p resonances
      fpLCleverMcLevy->InitReso(1, 2);  //number of Xi resonances
      // fpLCleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
      //                            massPion);
      // fpLCleverMcLevy->SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, massLambda,
      //                            massPion);
      AB_pL->SetAnaSource(CatsSourceForwarder, fpLCleverMcLevy, 2);
      AB_pL->SetAnaSource(0, 1.2);
      AB_pL->SetAnaSource(1, 2.0);
      break;
    case TidyCats::sLevy:
      fpLCleverLevy = new DLM_CleverLevy();
      fpLCleverMcLevy = new DLM_CleverMcLevyReso();
      fpLCleverMcLevy->InitNumMcIter(1000000);
      fpLCleverMcLevy->InitStability(20, 1, 2);
      fpLCleverMcLevy->InitScale(100, 0.1, 2.6);
      fpLCleverMcLevy->InitRad(512, 0, 64);
      fpLCleverMcLevy->InitType(2);
      fpLCleverMcLevy->InitReso(0, 1);  //number of p resonances
      fpLCleverMcLevy->InitReso(1, 1);  //number of Xi resonances
      fpLCleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
                                 massPion);
      fpLCleverMcLevy->SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, massLambda,
                                 massPion);
      //Nolan Parameterization.
      AB_pL->SetAnaSource(CatsSourceForwarder, fpLCleverMcLevy, 2);
      AB_pL->SetAnaSource(0, 1.2);  //r0
      AB_pL->SetAnaSource(1, 1.6);  //Stability alpha ( 1= Cauchy, ... 2 = Gauss)
      break;
    default:
      std::cout << "Source not implemented \n";
      break;
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
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;
  if (pot == TidyCats::pUsmani) {
    //!TEMPORARY: WE USE USMANI, ASK DIMI ABOUT FURTHER CHANGES (LIKE USE THE NLO)
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double pLamPotPars1S0[8] = { pL_UsmaniOli, 0, 0, 0, 0, 0, 0, 0 };
    double pLamPotPars3S1[8] = { pL_UsmaniOli, 0, 0, 0, 0, 1, 0, 1 };
    CATSparameters* cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,
                                                     8, true);
    cPotPars1S0->SetParameters(pLamPotPars1S0);
    CATSparameters* cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,
                                                     8, true);
    cPotPars3S1->SetParameters(pLamPotPars3S1);
    AB_pL->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
    AB_pL->SetShortRangePotential(1, 0, fDlmPot, *cPotPars3S1);
  } else if (pot == TidyCats::pNLOWF) {
    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    ExternalWF = Init_pL_Haidenbauer(
        TString::Format("%s/cernbox/WaveFunctions/Haidenbauer/pLambdaNLO/",
                        HomeDir.Data()),
        AB_pL, 10, 600);
    for (unsigned uCh = 0; uCh < AB_pL->GetNumChannels(); uCh++) {
      AB_pL->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                     ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(AB_pL->GetNumChannels(), ExternalWF);
  } else if (pot == TidyCats::pLOWF) {
    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    ExternalWF = Init_pL_Haidenbauer(
        TString::Format("%s/cernbox/WaveFunctions/Haidenbauer/pLambdaLO_600/",
                        HomeDir.Data()),
        AB_pL, 0, 600);
    for (unsigned uCh = 0; uCh < AB_pL->GetNumChannels(); uCh++) {
      AB_pL->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                     ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(AB_pL->GetNumChannels(), ExternalWF);
  } else {
    std::cout << "potential " << pot << "not implemented \n";
  }
  return;

//  if (pot == "LO") {
//    ExternalWF =
//        Init_pL_Haidenbauer(
//            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_600/",
//            Kitty, 0, 600);
//    NumChannels = 2;
//  } else if (pot == "LO_Coupled_S") {
//    ExternalWF =
//        Init_pL_Haidenbauer(
//            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_Coupling/",
//            Kitty, 1, 600);
//    NumChannels = 4;
//  } else if (pot == "NLO") {
//    ExternalWF =
//        Init_pL_Haidenbauer(
//            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",
//            Kitty, 10, 600);
//    NumChannels = 2;
//  } else if (pot == "NLO_Coupled_S") {
//    ExternalWF =
//        Init_pL_Haidenbauer(
//            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO_Coupling/",
//            Kitty, 11, 600);
//    NumChannels = 4;
//  }

}

void TidyCats::GetCatsProtonXiMinus(CATS* AB_pXim, int momBins, double kMin,
                                    double kMax, TidyCats::Sources source,
                                    TidyCats::pXimPot pot, double QCDTime) {
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
  CATSparameters* cPars = nullptr;
  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      AB_pXim->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance:
      fpXimCleverMcLevy = new DLM_CleverMcLevyReso();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpXimCleverMcLevy->InitScale(100, 0.1, 2.6);
      fpXimCleverMcLevy->InitRad(512, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->InitReso(0, 1);    //number of p resonances
//      fpXimCleverMcLevy->InitReso(1,1);//number of Xi resonances
      fpXimCleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p,
                                   massPion);
//      fpXimCleverMcLevy->SetUpReso(1,0,0,1361.52,1.65,Mass_p,massPion);
      AB_pXim->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      AB_pXim->SetAnaSource(0, 1.2);
      AB_pXim->SetAnaSource(1, 2.0);
      break;
    case TidyCats::sLevy:
      fpXimCleverMcLevy = new DLM_CleverMcLevyReso();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(20, 1, 2);
      fpXimCleverMcLevy->InitScale(100, 0.1, 2.6);
      fpXimCleverMcLevy->InitRad(512, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->InitReso(0, 1);    //number of p resonances
//      fpXimCleverMcLevy->InitReso(1,1);//number of Xi resonances
      fpXimCleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p,
                                   massPion);
//      fpXimCleverMcLevy->SetUpReso(1,0,0,1361.52,1.65,Mass_p,massPion);
      AB_pXim->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      AB_pXim->SetAnaSource(0, 1.2);
      AB_pXim->SetAnaSource(1, 2.0);
      break;
    default:
      std::cout << "Source not implemented \n";
      break;
  }
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
//  AB_pXim->SetGamow(true);
  AB_pXim->SetPdgId(2212, 3312);

  AB_pXim->SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim->SetMaxRad(64);
  AB_pXim->SetMaxRho(32);
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;
  if (pot == pHALQCD) {
    double pXimPotParsI0S0[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1 };  // 9, 10, 11, 12
    double pXimPotParsI1S0[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0 };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1 };  // this value 1-6
    CATSparameters* cPotParsI0S0 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI0S0->SetParameters(pXimPotParsI0S0);

    CATSparameters* cPotParsI0S1 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI0S1->SetParameters(pXimPotParsI0S1);

    CATSparameters* cPotParsI1S0 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI1S0->SetParameters(pXimPotParsI1S0);

    CATSparameters* cPotParsI1S1 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI1S1->SetParameters(pXimPotParsI1S1);
    AB_pXim->SetShortRangePotential(0, 0, fDlmPot, *cPotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, fDlmPot, *cPotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, fDlmPot, *cPotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, fDlmPot, *cPotParsI1S1);
  } else if (pot == pHALQCDGamow) {
    double pXimPotParsI0S0[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1 };  // 9, 10, 11, 12
    double pXimPotParsI1S0[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0 };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1 };  // this value 1-6
    CATSparameters* cPotParsI0S0 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI0S0->SetParameters(pXimPotParsI0S0);

    CATSparameters* cPotParsI0S1 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI0S1->SetParameters(pXimPotParsI0S1);

    CATSparameters* cPotParsI1S0 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI1S0->SetParameters(pXimPotParsI1S0);

    CATSparameters* cPotParsI1S1 = new CATSparameters(
        CATSparameters::tPotential, 8, true);
    cPotParsI1S1->SetParameters(pXimPotParsI1S1);
    AB_pXim->SetShortRangePotential(0, 0, fDlmPot, *cPotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, fDlmPot, *cPotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, fDlmPot, *cPotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, fDlmPot, *cPotParsI1S1);
    AB_pXim->SetQ1Q2(0);
    AB_pXim->SetGamow(true);
  } else if (pot == pHaidenbauer) {
    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    ExternalWF = Init_pXi_Haidenbauer(
        TString::Format("%s/cernbox/WaveFunctions/Haidenbauer/pXi_ver1/",
                        HomeDir.Data()),
        AB_pXim);
    for (unsigned uCh = 0; uCh < AB_pXim->GetNumChannels(); uCh++) {
      AB_pXim->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                       ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(AB_pXim->GetNumChannels(), ExternalWF);
  } else if (pot == pRikkenWF) {
    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    ExternalWF = Init_pXi_ESC16_IS(
        TString::Format("%s/cernbox/WaveFunctions/Tom/XN190401/",
                        HomeDir.Data()),
        AB_pXim);
    for (unsigned uCh = 0; uCh < AB_pXim->GetNumChannels(); uCh++) {
      AB_pXim->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                       ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(AB_pXim->GetNumChannels(), ExternalWF);
  } else if (pot == pRikkenPot) {
    CATSparameters* PotParsI0S0 = new CATSparameters(CATSparameters::tPotential,
                                                     2, true);
    PotParsI0S0->SetParameter(0, 0);
    PotParsI0S0->SetParameter(1, 0);

    CATSparameters* PotParsI0S1 = new CATSparameters(CATSparameters::tPotential,
                                                     2, true);
    PotParsI0S1->SetParameter(0, 0);
    PotParsI0S1->SetParameter(1, 1);

    CATSparameters* PotParsI1S0 = new CATSparameters(CATSparameters::tPotential,
                                                     2, true);
    PotParsI1S0->SetParameter(0, 1);
    PotParsI1S0->SetParameter(1, 0);

    CATSparameters* PotParsI1S1 = new CATSparameters(CATSparameters::tPotential,
                                                     2, true);
    PotParsI1S1->SetParameter(0, 1);
    PotParsI1S1->SetParameter(1, 1);

    //define the usual stuff
    AB_pXim->SetShortRangePotential(0, 0, ESC16_pXim_EXAMPLE, *PotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, ESC16_pXim_EXAMPLE, *PotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, ESC16_pXim_EXAMPLE, *PotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, ESC16_pXim_EXAMPLE, *PotParsI1S1);
  } else if (pot == pCoulomb) {
    std::cout << "Coulomb only! \n";
  } else if (pot == pGamow) {
    AB_pXim->SetQ1Q2(0);
    AB_pXim->SetGamow(true);
    std::cout << "Setting no potential at all + Gamow correction \n";
  } else {
    std::cout << "no implemented potential called \n";
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

  CATSparameters* cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, 1.2);
  AB_pXim->SetAnaSource(GaussSource, *cPars);
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
    CATSparameters* cPotParsI0S0 = new CATSparameters(
        CATSparameters::tPotential, 9, true);
    cPotParsI0S0->SetParameters(pXimPotParsI0S0);

    CATSparameters* cPotParsI0S1 = new CATSparameters(
        CATSparameters::tPotential, 9, true);
    cPotParsI0S1->SetParameters(pXimPotParsI0S1);

    CATSparameters* cPotParsI1S0 = new CATSparameters(
        CATSparameters::tPotential, 9, true);
    cPotParsI1S0->SetParameters(pXimPotParsI1S0);

    CATSparameters* cPotParsI1S1 = new CATSparameters(
        CATSparameters::tPotential, 9, true);
    cPotParsI1S1->SetParameters(pXimPotParsI1S1);
    AB_pXim->SetShortRangePotential(0, 0, fDlmPot, *cPotParsI0S0);
    AB_pXim->SetShortRangePotential(1, 0, fDlmPot, *cPotParsI0S1);
    AB_pXim->SetShortRangePotential(2, 0, fDlmPot, *cPotParsI1S0);
    AB_pXim->SetShortRangePotential(3, 0, fDlmPot, *cPotParsI1S1);
  }

  return;
}

void TidyCats::GetCatsProtonXiMinus1530(CATS* AB_pXim1530, int momBins,
                                        double kMin, double kMax,
                                        TidyCats::Sources source) {
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim1530 =
      TDatabasePDG::Instance()->GetParticle(3314)->Mass() * 1000;
  CATSparameters* cPars = nullptr;
  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      AB_pXim1530->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sLevy:
      fpXim1530CleverLevy = new DLM_CleverLevy();
      fpXim1530CleverLevy->InitStability(20, 1, 2);
      fpXim1530CleverLevy->InitScale(35, 0.25, 2.0);
      fpXim1530CleverLevy->InitRad(256, 0, 64);
      fpXim1530CleverLevy->InitType(2);
      //Nolan Parameterization.
      AB_pXim1530->SetAnaSource(CatsSourceForwarder, fpXim1530CleverLevy, 2);
      AB_pXim1530->SetAnaSource(0, 1.2);  //r0
      AB_pXim1530->SetAnaSource(1, 1.6);  //Stability alpha ( 1= Cauchy, ... 2 = Gauss)
      break;
    default:
      std::cout << "Source not implemented \n";
      break;
  }

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

void TidyCats::GetCatsProtonSigma0(CATS* AB_pSigma0, int momBins, double kMin,
                                   double kMax, TidyCats::Sources source,
                                   TidyCats::pSigma0Pot pot) {
  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double massSigma0 = TDatabasePDG::Instance()->GetParticle(3212)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
  CATSparameters* cPars = nullptr;
  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.3);
      AB_pSigma0->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance:
      fpSigma0CleverMcLevy = new DLM_CleverMcLevyReso();
      fpSigma0CleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpSigma0CleverMcLevy->InitScale(35, 0.25, 2.0);
      fpSigma0CleverMcLevy->InitRad(512, 0, 64);
      fpSigma0CleverMcLevy->InitType(2);
      fpSigma0CleverMcLevy->InitReso(0, 1);
      fpSigma0CleverMcLevy->InitReso(1, 1);
      fpSigma0CleverMcLevy->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65,
                                      massProton, massPion);
      fpSigma0CleverMcLevy->SetUpReso(1, 0, 1. - 0.3735, 1581.73, 4.28,
                                      massSigma0, massPion);
      fpSigma0CleverMcLevy->InitNumMcIter(1000000);
      AB_pSigma0->SetAnaSource(CatsSourceForwarder, fpSigma0CleverMcLevy, 2);
      AB_pSigma0->SetAnaSource(0, 0.743);  //r0
      AB_pSigma0->SetAnaSource(1, 2.0);  //Stability alpha ( 1= Cauchy, ... 2 = Gauss)
      break;
    default:
      std::cout << "Source not implemented \n";
      break;
  }

  AB_pSigma0->SetMomBins(momBins, kMin, kMax);
  AB_pSigma0->SetUseAnalyticSource(true);
  AB_pSigma0->SetMomentumDependentSource(false);
  AB_pSigma0->SetThetaDependentSource(false);
  AB_pSigma0->SetExcludeFailedBins(false);
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;
  TString HomeDir = gSystem->GetHomeDirectory().c_str();
  switch (pot) {
    case TidyCats::pSigma0Haidenbauer:
      ExternalWF = Init_pSigma0_Haidenbauer(
          TString::Format("%s/cernbox/SystematicsAndCalib/Sigma0/HaidenbauerWF/", HomeDir.Data()), AB_pSigma0);
      for (unsigned uCh = 0; uCh < AB_pSigma0->GetNumChannels(); uCh++) {
        AB_pSigma0->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                            ExternalWF[1][uCh][0]);
      }
      CleanUpWfHisto(AB_pSigma0->GetNumChannels(), ExternalWF);
      break;
    case TidyCats::pSigma0ESC16:
      ExternalWF = Init_pS0_ESC16(TString::Format("%s/cernbox/SystematicsAndCalib/Sigma0/ESC16/", HomeDir.Data()),
                                  AB_pSigma0);
      AB_pSigma0->SetExternalWaveFunction(0, 0, ExternalWF[0][0][0],
                                          ExternalWF[1][0][0]);
      AB_pSigma0->SetExternalWaveFunction(1, 0, ExternalWF[0][1][0],
                                          ExternalWF[1][1][0]);
      CleanUpWfHisto(AB_pSigma0->GetNumChannels(), ExternalWF);
      break;
    case TidyCats::pSigma0NSC97f:
      ExternalWF = Init_pSigma0_Haidenbauer(
          TString::Format("%s/cernbox/SystematicsAndCalib/Sigma0/NSC97f/", HomeDir.Data()),
          AB_pSigma0);
      for (unsigned uCh = 0; uCh < AB_pSigma0->GetNumChannels(); uCh++) {
        AB_pSigma0->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                            ExternalWF[1][uCh][0]);
      }
      CleanUpWfHisto(AB_pSigma0->GetNumChannels(), ExternalWF);
      break;
    default:
      std::cout << "Potential not implemented \n";
      break;
  }
}

double TidyCats::ESC16_pXim_EXAMPLE(double* Parameters) {
//[0] radius, [1] momentum
//[2] IsoSpin, [3] Spin
  if (Parameters[0] > 2)
    return 0;

  //change the path here
  TString HomeDir = gSystem->GetHomeDirectory().c_str();

  const double& Isospin = Parameters[2];
  const double& Spin = Parameters[3];
  TFile fIn(
      TString::Format(
          "%s/cernbox/WaveFunctions/Tom/Potential/ESC16_190319.root",
          HomeDir.Data()),
      "read");
  TGraph* gPot;

  if (Isospin == 0 && Spin == 0) {
    gPot = (TGraph*) fIn.Get("EffPotI0S0");
  } else if (Isospin == 0 && Spin == 1) {
    gPot = (TGraph*) fIn.Get("EffPotI0S1");
  } else if (Isospin == 1 && Spin == 0) {
    gPot = (TGraph*) fIn.Get("EffPotI1S0");
  } else if (Isospin == 1 && Spin == 1) {
    gPot = (TGraph*) fIn.Get("EffPotI1S1");
  } else {
    std::cout << "Isospin: " << Isospin << " Spin: " << Spin << std::endl;
    printf(" ESC16_pXim_EXAMPLE says ???????????????? \n");
    return 0;
  }
  return gPot->Eval(Parameters[0]);
}

DLM_Histo<double>* TidyCats::ConvertThetaAngleHisto(const TString& FileName, const TString& HistoName, const double kMin, const double kMax, bool convertToRad, int Rebin){
  TFile* InputFile = new TFile(FileName, "read");
  TH2F* InputHisto = NULL;
  if(InputFile){
    InputHisto = (TH2F*)InputFile->Get(HistoName);
  } else {
    std::cout << "no Input file found under the name " << FileName.Data() << " doei.. \n";  
    return nullptr; 
  }
  if (Rebin > 1) {
    InputHisto->Rebin2D(Rebin,Rebin); 
  } 
  TH1D* Projection=NULL;
  if(InputHisto){
    Projection = InputHisto->ProjectionX("ConvertThetaAngleHisto",InputHisto->GetYaxis()->FindBin(kMin),InputHisto->GetYaxis()->FindBin(kMax));
  } else {
    std::cout << "No input Histo found for " << HistoName.Data() << " exiting doei ... \n " ;
    return nullptr; 
  }

  const unsigned NumBins = Projection->GetNbinsX();
  Projection->Scale(1./Projection->Integral(1,NumBins));

  DLM_Histo<double> CummDistr;
  CummDistr.SetUp(1);
  if (convertToRad) {
    CummDistr.SetUp(0,NumBins,Projection->GetBinLowEdge(1)*DegToRad,Projection->GetXaxis()->GetBinUpEdge(NumBins)*DegToRad);
  } else {
    CummDistr.SetUp(0,NumBins,Projection->GetBinLowEdge(1),Projection->GetXaxis()->GetBinUpEdge(NumBins));
  }
  CummDistr.Initialize();

  CummDistr.SetBinContent(unsigned(0),Projection->GetBinContent(1));
  for(unsigned uBin=1; uBin<NumBins; uBin++){
    CummDistr.SetBinContent(uBin,CummDistr.GetBinContent(uBin-1)+Projection->GetBinContent(uBin+1));
  }

  double* MomBins = new double [NumBins+1];
  MomBins[0] = 0;
  //MomBins[1] = CummDistr.GetBinContent(unsigned(0));
  MomBins[NumBins] = 1;
  for(unsigned uBin=1; uBin<NumBins; uBin++){
    MomBins[uBin] = (CummDistr.GetBinContent(uBin)+CummDistr.GetBinContent(uBin-1)+1e-6)*0.5;
    //printf("MomBins[%u] = %f\n",uBin,MomBins[uBin]);
  }

  DLM_Histo<double>* Result = new DLM_Histo<double>();
  Result->SetUp(1);
  Result->SetUp(0,NumBins,MomBins);
  Result->Initialize();

  for(unsigned uBin=0; uBin<NumBins; uBin++){
    if(!Projection) break;
    Result->SetBinCenter(0,uBin,CummDistr.GetBinContent(uBin));//cum. value
    Result->SetBinContent(uBin,CummDistr.GetBinCenter(0,uBin));//momentum value
    //printf("%u: x=%.4f; y=%.4f;\n",uBin,CummDistr.GetBinContent(uBin),CummDistr.GetBinCenter(0,uBin));
  }

  delete [] MomBins;
  delete InputFile;
  return Result;
}

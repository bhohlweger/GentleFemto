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
#include "DLM_Random.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TError.h"

static double DegToRad = 0.01745;

TidyCats::TidyCats()
    : fHomeDir(),
      fkStarCutOff(200.),
      ftauProRes(1.65),
      fmassProRes(1361.52),
      ftauLamRes(4.69),
      fmassLamRes(1462.93),
      fppCleverMcLevy(nullptr),
      fpLCleverMcLevy(nullptr),
      fpXimCleverMcLevy(nullptr),
      fpSigma0CleverMcLevy(nullptr) {
  fHomeDir = gSystem->GetHomeDirectory().c_str();
}

TidyCats::~TidyCats() {
  if (fppCleverMcLevy) {
    delete fppCleverMcLevy;
  }
  if (fpLCleverMcLevy) {
    delete fpLCleverMcLevy;
  }
  if (fpXimCleverMcLevy) {
    delete fpXimCleverMcLevy;
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
  CATSparameters* cPars;

  if (source == TidyCats::sGaussian) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pp->SetAnaSource(GaussSource, *cPars);
  } else if (source == TidyCats::sResonance) {
    fppCleverMcLevy = new DLM_CleverMcLevyResoTM();
    fppCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    fppCleverMcLevy->InitScale(38, 0.15, 2.0);
    fppCleverMcLevy->InitRad(257, 0, 64);
    fppCleverMcLevy->InitType(2);

    fppCleverMcLevy->SetUpReso(0, 0.6422);
    fppCleverMcLevy->SetUpReso(1, 0.6422);
    DLM_Random RanGen(11);
    double RanVal1;
    double RanVal2;

    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;

    TFile* F_EposDisto_p_pReso = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_p_pReso.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_p_pReso = (TNtuple*) F_EposDisto_p_pReso->Get(
        "InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
    T_EposDisto_p_pReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_pReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_pReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_pReso; uEntry++) {
      T_EposDisto_p_pReso->GetEntry(uEntry);
      fM1 = 0;
      fM2 = fmassProRes;
      Tau1 = 0;
      Tau2 = ftauProRes;
      if (k_D > fkStarCutOff)
        continue;
//      RanVal1 = gRandom->Exp(fM2 / (fP2 * Tau2));
      RanVal1 = RanGen.Exponential(fM2 / (fP2 * Tau2));
      fppCleverMcLevy->AddBGT_PR(RanVal1, -cos(AngleRcP2));
      fppCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    TFile* F_EposDisto_pReso_pReso = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_pReso.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*) F_EposDisto_pReso_pReso->Get(
        "InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_pReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_pReso; uEntry++) {
      T_EposDisto_pReso_pReso->GetEntry(uEntry);
      fM1 = fmassProRes;
      fM2 = fmassProRes;
      Tau1 = ftauProRes;
      Tau2 = ftauProRes;
      if (k_D > fkStarCutOff)
        continue;
      RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
      RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
//      RanVal1 = gRandom->Exp(fM1 / (fP1 * Tau1));
//      RanVal2 = gRandom->Exp(fM2 / (fP2 * Tau2));
      fppCleverMcLevy->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2,
                                 cos(AngleRcP2), cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_pReso;

    fppCleverMcLevy->InitNumMcIter(262144);
    AB_pp->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
    AB_pp->SetAnaSource(0, 1.2);
    AB_pp->SetAnaSource(1, 2.0);
  } else if (source == TidyCats::sLevy) {
    fppCleverMcLevy = new DLM_CleverMcLevyResoTM();
    fppCleverMcLevy->InitStability(20, 1.f, 2.f);
    fppCleverMcLevy->InitScale(38, 0.15, 2.0);
    fppCleverMcLevy->InitRad(257, 0, 64);
    fppCleverMcLevy->InitType(2);
    fppCleverMcLevy->InitNumMcIter(262144);
    AB_pp->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
    AB_pp->SetAnaSource(0, 1.2);
    AB_pp->SetAnaSource(1, 2.0);
  } else if (source == TidyCats::sCauchy) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pp->SetAnaSource(CauchySource, *cPars);
  } else {
    std::cout << "Source not implemented \n";
  }
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
  CATSparameters* cPars;
  if (source == TidyCats::sGaussian) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pL->SetAnaSource(GaussSource, *cPars);
  } else if (source == TidyCats::sCauchy) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pL->SetAnaSource(CauchySource, *cPars);
  } else if (source == TidyCats::sResonance) {
    fpLCleverMcLevy = new DLM_CleverMcLevyResoTM();

    fpLCleverMcLevy->SetUpReso(0, 0.6422);
    fpLCleverMcLevy->SetUpReso(1, 0.6438);

    DLM_Random RanGen(11);
    double RanVal1;
    double RanVal2;

    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;

    TFile* F_EposDisto_p_LamReso = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_p_LamReso.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_p_LamReso = (TNtuple*) F_EposDisto_p_LamReso->Get(
        "InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
    T_EposDisto_p_LamReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_LamReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_LamReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_LamReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_LamReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_LamReso; uEntry++) {
      T_EposDisto_p_LamReso->GetEntry(uEntry);
      fM1 = 0;
      fM2 = fmassLamRes;
      Tau1 = 0;
      Tau2 = ftauLamRes;
      if (k_D > fkStarCutOff)
        continue;
      RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
      fpLCleverMcLevy->AddBGT_PR(RanVal2, cos(AngleRcP2));
    }
    delete F_EposDisto_p_LamReso;

    TFile* F_EposDisto_pReso_Lam = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Lam.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_pReso_Lam = (TNtuple*) F_EposDisto_pReso_Lam->Get(
        "InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
    T_EposDisto_pReso_Lam->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Lam->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Lam->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Lam; uEntry++) {
      T_EposDisto_pReso_Lam->GetEntry(uEntry);
      fM1 = fmassProRes;
      fM2 = 0;
      Tau1 = ftauProRes;
      Tau2 = 0;
      if (k_D > fkStarCutOff)
        continue;
      RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
      fpLCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete F_EposDisto_pReso_Lam;

    TFile* F_EposDisto_pReso_LamReso = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_LamReso.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*) F_EposDisto_pReso_LamReso
        ->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_LamReso =
        T_EposDisto_pReso_LamReso->GetEntries();
    T_EposDisto_pReso_LamReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_LamReso; uEntry++) {
      T_EposDisto_pReso_LamReso->GetEntry(uEntry);
      fM1 = fmassProRes;
      fM2 = fmassLamRes;
      Tau1 = ftauProRes;
      Tau2 = ftauLamRes;
      if (k_D > fkStarCutOff)
        continue;
      RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
      RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
      fpLCleverMcLevy->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2,
                                 cos(AngleRcP2), cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_LamReso;
    fpLCleverMcLevy->InitNumMcIter(262144);
    AB_pL->SetAnaSource(CatsSourceForwarder, fpLCleverMcLevy, 2);
    AB_pL->SetAnaSource(0, 1.2);
    AB_pL->SetAnaSource(1, 2.0);
  } else {
    std::cout << "Source not implemented \n";
  }
  AB_pL->SetUseAnalyticSource(true);
  AB_pL->SetMomentumDependentSource(false);
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
  if (source == TidyCats::sGaussian) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pXim->SetAnaSource(GaussSource, *cPars);
  } else if (source == TidyCats::sCauchy) {
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 1.2);
    AB_pXim->SetAnaSource(CauchySource, *cPars);
  } else if (source == TidyCats::sResonance) {
    fpXimCleverMcLevy = new DLM_CleverMcLevyResoTM();
    fpXimCleverMcLevy->InitNumMcIter(1000000);
    fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    fpXimCleverMcLevy->InitScale(38, 0.15, 2.0);
    fpXimCleverMcLevy->InitRad(257, 0, 64);
    fpXimCleverMcLevy->InitType(2);
    fpXimCleverMcLevy->SetUpReso(0, 0.6422);
    DLM_Random RanGen(11);
    double RanVal1;
    double RanVal2;

    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;

    TFile* F_EposDisto_pReso_Xi = new TFile(
        TString::Format(
            "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Xim.root",
            fHomeDir.Data()).Data());
    TNtuple* T_EposDisto_pReso_Xi = (TNtuple*) F_EposDisto_pReso_Xi->Get(
        "InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xi->GetEntries();
    T_EposDisto_pReso_Xi->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Xi->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Xi->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Xi->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Xi->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Xi->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Xi->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Xi->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Xi->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Xi->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Xim; uEntry++) {
      T_EposDisto_pReso_Xi->GetEntry(uEntry);
      fM1 = fmassProRes;
      fM2 = 0;
      Tau1 = ftauProRes;
      Tau2 = 0;
      if (k_D > fkStarCutOff)
        continue;
      RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
      fpXimCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete F_EposDisto_pReso_Xi;
    fpXimCleverMcLevy->InitNumMcIter(262144);
    AB_pXim->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
    AB_pXim->SetAnaSource(0, 1.2);
    AB_pXim->SetAnaSource(1, 2.0);
  } else {
    std::cout << "Source not implemented \n";
  }
  AB_pXim->SetUseAnalyticSource(true);
  AB_pXim->SetThetaDependentSource(false);
  AB_pXim->SetMomentumDependentSource(false);

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
    /*
      double pXimPotParsI0S0[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1 };  // 9, 10, 11, 12
    double pXimPotParsI1S0[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0 };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1 };  // this value 1-6
    */
    
    double pXimPotParsI0S0[8] = { pXim_HALQCDPaper2020, QCDTime, 0, -1, 1, 0, 0, 0 };  //4th argument is the t parameter and can be:
    double pXimPotParsI0S1[8] = { pXim_HALQCDPaper2020, QCDTime, 0, -1, 1, 1, 0, 1 };  // 9, 10, 11, 12
    double pXimPotParsI1S0[8] = { pXim_HALQCDPaper2020, QCDTime, 1,  1, 1, 0, 0, 0 };  //This is shit. Corresponds to 9-14 t
    double pXimPotParsI1S1[8] = { pXim_HALQCDPaper2020, QCDTime, 1,  1, 1, 1, 0, 1 };  // this value 1-6
    
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
    case TidyCats::sCauchy:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      AB_pXim1530->SetAnaSource(CauchySource, *cPars);
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

void TidyCats::GetCatsProtonDplus(CATS* cats, int momBins, double kMin,
                                  double kMax, TidyCats::pDmesonPot pot,
                                  TidyCats::Sources source) {

  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massProton = pdgDatabase->GetParticle(2212)->Mass() * 1000;
  const double massDplus = pdgDatabase->GetParticle(411)->Mass() * 1000;
  CATSparameters* cPars;

  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      cats->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance: {
      fpXimCleverMcLevy = new DLM_CleverMcLevyResoTM();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpXimCleverMcLevy->InitScale(38, 0.15, 2.0);
      fpXimCleverMcLevy->InitRad(257, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->SetUpReso(0, 0.6422);
      DLM_Random RanGen(11);
      double RanVal1, RanVal2;
      float k_D, fP1, fP2, fM1, fM2, Tau1, Tau2, AngleRcP1, AngleRcP2,
          AngleP1P2;

      // TODO For the moment we borrow the kinematics from the Xi which is closest in mass to the D
      TFile* F_EposDisto_pReso_D = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Xim.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_D = (TNtuple*) F_EposDisto_pReso_D->Get(
          "InfoTuple_ClosePairs");
      T_EposDisto_pReso_D->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_D->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_D->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_D->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_D->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_D->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_D->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < T_EposDisto_pReso_D->GetEntries(); uEntry++) {
        T_EposDisto_pReso_D->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = 0;
        Tau1 = ftauProRes;
        Tau2 = 0;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        fpXimCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
      }
      delete F_EposDisto_pReso_D;
      fpXimCleverMcLevy->InitNumMcIter(262144);
      cats->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      cats->SetAnaSource(0, 1.2);
      cats->SetAnaSource(1, 2.0);
      break;
    }
    default:
      std::cout << "Source not implemented \n";
      break;
  }
  cats->SetUseAnalyticSource(true);
  cats->SetMomentumDependentSource(false);
  cats->SetThetaDependentSource(false);
  cats->SetMomBins(momBins, kMin, kMax);

//  TODO To be figured out
//  cats->SetNumChannels(2);
//  cats->SetNumPW(0, 1);
//  cats->SetNumPW(1, 1);
//  cats->SetSpin(0, 0);
//  cats->SetSpin(1, 1);
//  cats->SetChannelWeight(0, 0.25);
//  cats->SetChannelWeight(1, 0.75);
  cats->SetNumChannels(1);
  cats->SetNumPW(0, 1);
  cats->SetSpin(0, 0);
  cats->SetChannelWeight(0, 1.);

  cats->SetQ1Q2(1);
  cats->SetPdgId(2212, 411);
  cats->SetRedMass((massProton * massDplus) / (massProton + massDplus));
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;

  switch (pot) {
    case TidyCats::pCoulombOnly:
      std::cout << "Coulomb only\n";
      break;
    default:
      std::cout << "Potential not implemented\n";
      break;
  }
  return;
}

void TidyCats::GetCatsProtonDminus(CATS* cats, int momBins, double kMin,
                                   double kMax, TidyCats::pDmesonPot pot,
                                   TidyCats::Sources source) {

  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massProton = pdgDatabase->GetParticle(2212)->Mass() * 1000;
  const double massDminus = pdgDatabase->GetParticle(-411)->Mass() * 1000;
  CATSparameters* cPars;

  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      cats->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance: {
      fpXimCleverMcLevy = new DLM_CleverMcLevyResoTM();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpXimCleverMcLevy->InitScale(38, 0.15, 2.0);
      fpXimCleverMcLevy->InitRad(257, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->SetUpReso(0, 0.6422);
      DLM_Random RanGen(11);
      double RanVal1, RanVal2;
      float k_D, fP1, fP2, fM1, fM2, Tau1, Tau2, AngleRcP1, AngleRcP2,
          AngleP1P2;

      // TODO For the moment we borrow the kinematics from the Xi which is closest in mass to the D
      TFile* F_EposDisto_pReso_D = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Xim.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_D = (TNtuple*) F_EposDisto_pReso_D->Get(
          "InfoTuple_ClosePairs");
      T_EposDisto_pReso_D->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_D->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_D->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_D->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_D->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_D->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_D->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < T_EposDisto_pReso_D->GetEntries();
          uEntry++) {
        T_EposDisto_pReso_D->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = 0;
        Tau1 = ftauProRes;
        Tau2 = 0;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        fpXimCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
      }
      delete F_EposDisto_pReso_D;
      fpXimCleverMcLevy->InitNumMcIter(262144);
      cats->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      cats->SetAnaSource(0, 1.2);
      cats->SetAnaSource(1, 2.0);
      break;
    }
    default:
      std::cout << "Source not implemented \n";
      break;
  }
  cats->SetUseAnalyticSource(true);
  cats->SetMomentumDependentSource(false);
  cats->SetThetaDependentSource(false);
  cats->SetMomBins(momBins, kMin, kMax);

//  TODO To be figured out
//  cats->SetNumChannels(2);
//  cats->SetNumPW(0, 1);
//  cats->SetNumPW(1, 1);
//  cats->SetSpin(0, 0);
//  cats->SetSpin(1, 1);
//  cats->SetChannelWeight(0, 0.25);
//  cats->SetChannelWeight(1, 0.75);
  cats->SetNumChannels(1);
  cats->SetNumPW(0, 1);
  cats->SetSpin(0, 0);
  cats->SetChannelWeight(0, 1.);

  cats->SetQ1Q2(-1);
  cats->SetPdgId(2212, -411);
  cats->SetRedMass((massProton * massDminus) / (massProton + massDminus));
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;

  switch (pot) {
    case TidyCats::pCoulombOnly:
      std::cout << "Coulomb only\n";
      break;
    default:
      std::cout << "Potential not implemented\n";
      break;
  }
  return;
}

void TidyCats::GetCatsProtonDstarplus(CATS* cats, int momBins, double kMin,
                                      double kMax, TidyCats::pDmesonPot pot,
                                      TidyCats::Sources source) {

  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massProton = pdgDatabase->GetParticle(2212)->Mass() * 1000;
  const double massDstarplus = pdgDatabase->GetParticle(413)->Mass() * 1000;
  CATSparameters* cPars;

  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      cats->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance: {
      fpXimCleverMcLevy = new DLM_CleverMcLevyResoTM();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpXimCleverMcLevy->InitScale(38, 0.15, 2.0);
      fpXimCleverMcLevy->InitRad(257, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->SetUpReso(0, 0.6422);
      DLM_Random RanGen(11);
      double RanVal1, RanVal2;
      float k_D, fP1, fP2, fM1, fM2, Tau1, Tau2, AngleRcP1, AngleRcP2,
          AngleP1P2;

      // TODO For the moment we borrow the kinematics from the Xi which is closest in mass to the D
      TFile* F_EposDisto_pReso_D = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Xim.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_D = (TNtuple*) F_EposDisto_pReso_D->Get(
          "InfoTuple_ClosePairs");
      T_EposDisto_pReso_D->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_D->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_D->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_D->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_D->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_D->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_D->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < T_EposDisto_pReso_D->GetEntries(); uEntry++) {
        T_EposDisto_pReso_D->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = 0;
        Tau1 = ftauProRes;
        Tau2 = 0;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        fpXimCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
      }
      delete F_EposDisto_pReso_D;
      fpXimCleverMcLevy->InitNumMcIter(262144);
      cats->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      cats->SetAnaSource(0, 1.2);
      cats->SetAnaSource(1, 2.0);
      break;
    }
    default:
      std::cout << "Source not implemented \n";
      break;
  }
  cats->SetUseAnalyticSource(true);
  cats->SetMomentumDependentSource(false);
  cats->SetThetaDependentSource(false);
  cats->SetMomBins(momBins, kMin, kMax);

//  TODO To be figured out
//  cats->SetNumChannels(2);
//  cats->SetNumPW(0, 1);
//  cats->SetNumPW(1, 1);
//  cats->SetSpin(0, 0);
//  cats->SetSpin(1, 1);
//  cats->SetChannelWeight(0, 0.25);
//  cats->SetChannelWeight(1, 0.75);
  cats->SetNumChannels(1);
  cats->SetNumPW(0, 1);
  cats->SetSpin(0, 0);
  cats->SetChannelWeight(0, 1.);

  cats->SetQ1Q2(1);
  cats->SetPdgId(2212, 413);
  cats->SetRedMass((massProton * massDstarplus) / (massProton + massDstarplus));
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;

  switch (pot) {
    case TidyCats::pCoulombOnly:
      std::cout << "Coulomb only\n";
      break;
    default:
      std::cout << "Potential not implemented\n";
      break;
  }
  return;
}

void TidyCats::GetCatsProtonDstarminus(CATS* cats, int momBins, double kMin,
                                       double kMax, TidyCats::pDmesonPot pot,
                                       TidyCats::Sources source) {

  const auto pdgDatabase = TDatabasePDG::Instance();
  const double massProton = pdgDatabase->GetParticle(2212)->Mass() * 1000;
  const double massDstarminus = pdgDatabase->GetParticle(-413)->Mass() * 1000;
  CATSparameters* cPars;

  switch (source) {
    case TidyCats::sGaussian:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.2);
      cats->SetAnaSource(GaussSource, *cPars);
      break;
    case TidyCats::sResonance: {
      fpXimCleverMcLevy = new DLM_CleverMcLevyResoTM();
      fpXimCleverMcLevy->InitNumMcIter(1000000);
      fpXimCleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
      fpXimCleverMcLevy->InitScale(38, 0.15, 2.0);
      fpXimCleverMcLevy->InitRad(257, 0, 64);
      fpXimCleverMcLevy->InitType(2);
      fpXimCleverMcLevy->SetUpReso(0, 0.6422);
      DLM_Random RanGen(11);
      double RanVal1, RanVal2;
      float k_D, fP1, fP2, fM1, fM2, Tau1, Tau2, AngleRcP1, AngleRcP2,
          AngleP1P2;

      // TODO For the moment we borrow the kinematics from the Xi which is closest in mass to the D
      TFile* F_EposDisto_pReso_D = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Xim.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_D = (TNtuple*) F_EposDisto_pReso_D->Get(
          "InfoTuple_ClosePairs");
      T_EposDisto_pReso_D->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_D->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_D->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_D->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_D->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_D->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_D->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_D->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_D->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < T_EposDisto_pReso_D->GetEntries();
          uEntry++) {
        T_EposDisto_pReso_D->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = 0;
        Tau1 = ftauProRes;
        Tau2 = 0;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        fpXimCleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
      }
      delete F_EposDisto_pReso_D;
      fpXimCleverMcLevy->InitNumMcIter(262144);
      cats->SetAnaSource(CatsSourceForwarder, fpXimCleverMcLevy, 2);
      cats->SetAnaSource(0, 1.2);
      cats->SetAnaSource(1, 2.0);
      break;
    }
    default:
      std::cout << "Source not implemented \n";
      break;
  }
  cats->SetUseAnalyticSource(true);
  cats->SetMomentumDependentSource(false);
  cats->SetThetaDependentSource(false);
  cats->SetMomBins(momBins, kMin, kMax);

//  TODO To be figured out
//  cats->SetNumChannels(2);
//  cats->SetNumPW(0, 1);
//  cats->SetNumPW(1, 1);
//  cats->SetSpin(0, 0);
//  cats->SetSpin(1, 1);
//  cats->SetChannelWeight(0, 0.25);
//  cats->SetChannelWeight(1, 0.75);
  cats->SetNumChannels(1);
  cats->SetNumPW(0, 1);
  cats->SetSpin(0, 0);
  cats->SetChannelWeight(0, 1.);

  cats->SetQ1Q2(-1);
  cats->SetPdgId(2212, -411);
  cats->SetRedMass((massProton * massDstarminus) / (massProton + massDstarminus));
  DLM_Histo<complex<double>>*** ExternalWF = nullptr;

  switch (pot) {
    case TidyCats::pCoulombOnly:
      std::cout << "Coulomb only\n";
      break;
    default:
      std::cout << "Potential not implemented\n";
      break;
  }
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
    case TidyCats::sResonance: {
      fmassLamRes = 1581.73;
      ftauLamRes = 4.28;
      fpSigma0CleverMcLevy = new DLM_CleverMcLevyResoTM();
      fpSigma0CleverMcLevy->SetUpReso(0, 0.6422);
      fpSigma0CleverMcLevy->SetUpReso(1, 0.3735);
      DLM_Random RanGen(11);
      double RanVal1;
      double RanVal2;
      Float_t k_D;
      Float_t fP1;
      Float_t fP2;
      Float_t fM1;
      Float_t fM2;
      Float_t Tau1;
      Float_t Tau2;
      Float_t AngleRcP1;
      Float_t AngleRcP2;
      Float_t AngleP1P2;
      // there are no EPOS productions at the moment for p-Sigma0
      // since the masses and lifetimes are close to the ones of the Lambda this is considered
      TFile* F_EposDisto_p_LamReso = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_p_LamReso.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_p_LamReso = (TNtuple*) F_EposDisto_p_LamReso->Get(
          "InfoTuple_ClosePairs");
      unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
      T_EposDisto_p_LamReso->SetBranchAddress("k_D", &k_D);
      T_EposDisto_p_LamReso->SetBranchAddress("P1", &fP1);
      T_EposDisto_p_LamReso->SetBranchAddress("P2", &fP2);
      T_EposDisto_p_LamReso->SetBranchAddress("M1", &fM1);
      T_EposDisto_p_LamReso->SetBranchAddress("M2", &fM2);
      T_EposDisto_p_LamReso->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_p_LamReso->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < N_EposDisto_p_LamReso; uEntry++) {
        T_EposDisto_p_LamReso->GetEntry(uEntry);
        fM1 = 0;
        fM2 = fmassLamRes;
        Tau1 = 0;
        Tau2 = ftauLamRes;
        if (k_D > fkStarCutOff)
          continue;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        fpSigma0CleverMcLevy->AddBGT_PR(RanVal2, cos(AngleRcP2));
      }
      delete F_EposDisto_p_LamReso;

      // there are no EPOS productions at the moment for p-Sigma0
      // since the masses and lifetimes are close to the ones of the Lambda this is considered
      TFile* F_EposDisto_pReso_Lam = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_Lam.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_Lam = (TNtuple*) F_EposDisto_pReso_Lam->Get(
          "InfoTuple_ClosePairs");
      unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
      T_EposDisto_pReso_Lam->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_Lam->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_Lam->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_Lam->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_Lam->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_Lam->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_Lam->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Lam; uEntry++) {
        T_EposDisto_pReso_Lam->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = 0;
        Tau1 = ftauProRes;
        Tau2 = 0;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        fpSigma0CleverMcLevy->AddBGT_RP(RanVal1, cos(AngleRcP1));
      }
      delete F_EposDisto_pReso_Lam;
      // there are no EPOS productions at the moment for p-Sigma0
      // since the masses and lifetimes are close to the ones of the Lambda this is considered
      TFile* F_EposDisto_pReso_LamReso = new TFile(
          TString::Format(
              "%s/cernbox/WaveFunctions/ThetaDist/EposDisto_pReso_LamReso.root",
              fHomeDir.Data()).Data());
      TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*) F_EposDisto_pReso_LamReso
          ->Get("InfoTuple_ClosePairs");
      unsigned N_EposDisto_pReso_LamReso =
          T_EposDisto_pReso_LamReso->GetEntries();
      T_EposDisto_pReso_LamReso->SetBranchAddress("k_D", &k_D);
      T_EposDisto_pReso_LamReso->SetBranchAddress("P1", &fP1);
      T_EposDisto_pReso_LamReso->SetBranchAddress("P2", &fP2);
      T_EposDisto_pReso_LamReso->SetBranchAddress("M1", &fM1);
      T_EposDisto_pReso_LamReso->SetBranchAddress("M2", &fM2);
      T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1", &Tau1);
      T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2", &Tau2);
      T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
      T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
      T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
      for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_LamReso; uEntry++) {
        T_EposDisto_pReso_LamReso->GetEntry(uEntry);
        fM1 = fmassProRes;
        fM2 = fmassLamRes;
        Tau1 = ftauProRes;
        Tau2 = ftauLamRes;
        if (k_D > fkStarCutOff)
          continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        fpSigma0CleverMcLevy->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2,
                                        cos(AngleRcP2), cos(AngleP1P2));
      }
      delete F_EposDisto_pReso_LamReso;
      fpSigma0CleverMcLevy->InitNumMcIter(262144);
      AB_pSigma0->SetAnaSource(CatsSourceForwarder, fpSigma0CleverMcLevy, 2);
      AB_pSigma0->SetAnaSource(0, 1.2);
      AB_pSigma0->SetAnaSource(1, 2.0);
      break;
    }
    case TidyCats::sCauchy:
      cPars = new CATSparameters(CATSparameters::tSource, 1, true);
      cPars->SetParameter(0, 1.);
      AB_pSigma0->SetAnaSource(CauchySource, *cPars);
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
          TString::Format(
              "%s/cernbox/SystematicsAndCalib/Sigma0/HaidenbauerWF/",
              HomeDir.Data()),
          AB_pSigma0);
      for (unsigned uCh = 0; uCh < AB_pSigma0->GetNumChannels(); uCh++) {
        AB_pSigma0->SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                            ExternalWF[1][uCh][0]);
      }
      CleanUpWfHisto(AB_pSigma0->GetNumChannels(), ExternalWF);
      break;
    case TidyCats::pSigma0ESC16:
      ExternalWF = Init_pS0_ESC16(
          TString::Format("%s/cernbox/SystematicsAndCalib/Sigma0/ESC16/",
                          HomeDir.Data()),
          AB_pSigma0);
      AB_pSigma0->SetExternalWaveFunction(0, 0, ExternalWF[0][0][0],
                                          ExternalWF[1][0][0]);
      AB_pSigma0->SetExternalWaveFunction(1, 0, ExternalWF[0][1][0],
                                          ExternalWF[1][1][0]);
      CleanUpWfHisto(AB_pSigma0->GetNumChannels(), ExternalWF);
      break;
    case TidyCats::pSigma0NSC97f:
      ExternalWF = Init_pSigma0_Haidenbauer(
          TString::Format("%s/cernbox/SystematicsAndCalib/Sigma0/NSC97f/",
                          HomeDir.Data()),
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

DLM_Histo<double>* TidyCats::ConvertThetaAngleHisto(const TString& FileName,
                                                    const TString& HistoName,
                                                    const double kMin,
                                                    const double kMax,
                                                    bool convertToRad,
                                                    int Rebin) {
  std::cout << "Angular File Name: " << FileName.Data() << std::endl;
  TFile* InputFile = new TFile(FileName, "read");
  TH2F* InputHisto = NULL;
  if (InputFile) {
    InputHisto = (TH2F*) InputFile->Get(HistoName);
  } else {
    std::cout << "no Input file found under the name " << FileName.Data()
              << " doei.. \n";
    return nullptr;
  }
  if (Rebin > 1) {
    InputHisto->Rebin2D(Rebin, Rebin);
  }
  TH1D* Projection = NULL;
  if (InputHisto) {
    Projection = InputHisto->ProjectionX("ConvertThetaAngleHisto",
                                         InputHisto->GetYaxis()->FindBin(kMin),
                                         InputHisto->GetYaxis()->FindBin(kMax));
  } else {
    std::cout << "No input Histo found for " << HistoName.Data()
              << " exiting doei ... \n ";
    return nullptr;
  }

  const unsigned NumBins = Projection->GetNbinsX();
  Projection->Scale(1. / Projection->Integral(1, NumBins));

  DLM_Histo<double> CummDistr;
  CummDistr.SetUp(1);
  if (convertToRad) {
    CummDistr.SetUp(0, NumBins, Projection->GetBinLowEdge(1) * DegToRad,
                    Projection->GetXaxis()->GetBinUpEdge(NumBins) * DegToRad);
  } else {
    CummDistr.SetUp(0, NumBins, Projection->GetBinLowEdge(1),
                    Projection->GetXaxis()->GetBinUpEdge(NumBins));
  }
  CummDistr.Initialize();

  CummDistr.SetBinContent(unsigned(0), Projection->GetBinContent(1));
  for (unsigned uBin = 1; uBin < NumBins; uBin++) {
    CummDistr.SetBinContent(
        uBin,
        CummDistr.GetBinContent(uBin - 1)
            + Projection->GetBinContent(uBin + 1));
  }

  double* MomBins = new double[NumBins + 1];
  MomBins[0] = 0;
  //MomBins[1] = CummDistr.GetBinContent(unsigned(0));
  //MomBins[NumBins] = 1;
  int theHappening = 1;  //this is done in order to make sure that if two 0 bins follow up on each other shit doesn't go to the toilet.
  for (unsigned uBin = 1; uBin <= NumBins; uBin++) {
    MomBins[uBin] = (CummDistr.GetBinContent(uBin)
        + CummDistr.GetBinContent(uBin - 1)) * 0.5;
    if (MomBins[uBin] <= MomBins[uBin - 1]) {
      MomBins[uBin] = MomBins[uBin - 1] + 1e-6;
    }
    //printf("MomBins[%u] = %e\n",uBin,MomBins[uBin]);
  }
  DLM_Histo<double>* Result = new DLM_Histo<double>();
  Result->SetUp(1);
  Result->SetUp(0, NumBins, MomBins);
  Result->Initialize();

  for (unsigned uBin = 0; uBin < NumBins; uBin++) {
    if (!Projection)
      break;
    Result->SetBinCenter(0, uBin, CummDistr.GetBinContent(uBin));  //cum. value
    Result->SetBinContent(uBin, CummDistr.GetBinCenter(0, uBin));  //momentum value
    //printf("%u: x=%.4f; y=%.4f;\n",uBin,CummDistr.GetBinContent(uBin),CummDistr.GetBinCenter(0,uBin));
  }

  delete[] MomBins;
  delete InputFile;
  return Result;
}

TH2F* TidyCats::ConvertHisto(TH2F* input, int nBins, double kMin, double kMax) { 
  //only converts the histo along the y axis
  
  TString Histname = TString::Format("%sFixShiftRebinned", input->GetName());
  int nBinsX = input->GetXaxis()->GetNbins();
  double xMin = input->GetXaxis()->GetXmin();
  double xMax = input->GetXaxis()->GetXmax();
  //  int nBinsY = input->GetYaxis()->FindBin(kMax) - input->GetYaxis()->FindBin(kMin); 
  std::cout << "Convert Histo ... nBins: " << nBins << " kMin: " << kMin << " kMax: " << kMax << std::endl;
  
  TH2F* out = new TH2F(Histname.Data(),Histname.Data(), nBinsX, xMin, xMax, nBins, kMin, kMax);
  
  for (int iBinX = 1; iBinX <= nBinsX; ++iBinX) {
    for (int iBinY = 1; iBinY <= nBins; ++iBinY) {
      int binY = input->GetYaxis()->FindBin(out->GetYaxis()->GetBinCenter(iBinY));  
      out->SetBinContent(iBinX,iBinY,input->GetBinContent(iBinX, binY));
      out->SetBinError(iBinX,iBinY,input->GetBinError(iBinX, binY)); 
    }
  }
  return out; 
} 

void TidyCats::Smear(CATS& CATS, TH2F* smearing, TH1F* Smeared) {
  int nBinsY = Smeared->GetXaxis()->GetNbins();
  int nBinsX = smearing->GetXaxis()->GetNbins(); 
  for (int iksOut = 1; iksOut <= nBinsY; ++iksOut) {
    double Norm = 0;
    double CkAdded = 0; 
    double ksTrans = Smeared->GetBinCenter(iksOut); 
    int ksMrx = smearing->GetYaxis()->FindBin(ksTrans); 
    for (int iksIn = 1; iksIn <= nBinsX; ++iksIn) {
      double ksInMtrx = smearing->GetXaxis()->GetBinCenter(iksIn);
      double ksInCATS = CATS.GetMomentum(iksIn-1);
      if (TMath::Abs(ksInMtrx-ksInCATS) > 1e-2) {
	std::cout << "Warning in smearing difference in momenta between matrix and CATS object: ksCATS = " << ksInCATS << " and ksMATRIX = " << ksInMtrx << std::endl; 
      } 
      Norm += smearing->GetBinContent(iksIn, ksMrx);
      CkAdded += smearing->GetBinContent(iksIn, ksMrx)*CATS.GetCorrFun(iksIn-1);
      // if (iksOut == 1) {
      // 	std::cout << "iksIn: " << iksIn <<  " Cats: " << CATS.GetCorrFun(iksIn-1) << " Smearer: " << smearing->GetBinContent(iksIn, ksMrx) << std::endl;
      // }
    }
    // normalize
    CkAdded/=Norm; 
    Smeared->SetBinContent(iksOut, CkAdded); 
  }
  return; 
}

DLM_Histo<double>* TidyCats::Convert2LargerOf2Evils(TH1F* CkInput) {
  DLM_Histo<double>* output = new DLM_Histo<double>();
  int nBins = CkInput->GetXaxis()->GetNbins();
  output->SetUp(1);
  output->SetUp(1, nBins, CkInput->GetXaxis()->GetXmin(),
                CkInput->GetXaxis()->GetXmax());
  output->Initialize();
  for (int iBims = 1; iBims < nBins + 1; iBims++) {
    if (TMath::Abs(
        output->GetBinCenter(0, iBims) - CkInput->GetBinCenter(iBims - 1))
        < 1e-3) {
      Error("TidyCats::Convert2LargerOf2Evils", "Bin Centers differ!");
      return nullptr;
    }
    output->SetBinContent(iBims, CkInput->GetBinContent(iBims - 1));
    output->SetBinError(iBims, CkInput->GetBinError(iBims - 1));
  }
  return output;
}

TH1F* TidyCats::Convert2LesserOf2Evils(DLM_Histo<double>* CkInput, TH1F* dim) {
  TH1F* output = nullptr;
  int nBins;
  if (dim) { 
    nBins = dim->GetNbinsX();
    output = new TH1F("TidyCats::RenameMe", "TidyCats::RenameMe", nBins,
                      dim->GetXaxis()->GetXmin(), dim->GetXaxis()->GetXmax());
  } else {
    nBins = CkInput->GetNbins(0);
    output = new TH1F("TidyCats::RenameMe", "TidyCats::RenameMe", nBins,
                      CkInput->GetLowEdge(0), CkInput->GetUpEdge(0));
  }
  for (int iBims = 1; iBims < nBins + 1; iBims++) {
    if (TMath::Abs(
        output->GetBinCenter(iBims) - CkInput->GetBinCenter(0, iBims - 1))
        > 1e-3) {
      Error("TidyCats::Convert2LesserOf2Evils", "Bin Centers differ!");
      std::cout << "iBin: " << iBims << " output->GetBinCenter(iBims): "
                << output->GetBinCenter(iBims)
                << " CkInput->GetBinCenter(0, iBims - 1): "
                << CkInput->GetBinCenter(0, iBims - 1) << std::endl;
      return nullptr;
    }
    output->SetBinContent(iBims, CkInput->GetBinContent(iBims - 1));
    output->SetBinError(iBims, CkInput->GetBinError(iBims - 1));
  }
  return output;
}

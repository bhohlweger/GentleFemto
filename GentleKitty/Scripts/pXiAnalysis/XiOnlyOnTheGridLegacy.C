#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include "stdlib.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "SideBandFit.h"
#include "TMinuit.h"
void GetXiForRadius(const unsigned& NumIter, int system, int iPot, int iSource,
                    TString OutputDir) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  //System: 0 = pPb, 1 = pp MB, 2 = pp HM
  //Potential: 0 = Coulomb, 1 = Gamow, 2 = HAL + Coulomb, 3 = HAL + Gamow, 5 = RikkenWF , 6 = RikkenPot
  //Source: 0 = Gauss, 1 = Resonances, 2 = Levy
  //Iteration number corresponds here to the systematic cut variation being fit.
  bool drawSource = (NumIter == 0);
  bool fitRadius = false;
  double BlRegion[2];
  BlRegion[0] = 500;
  BlRegion[1] = 500;

  double normvarCont[3][2];
  normvarCont[0][0] = 450;
  normvarCont[0][1] = 650;
  normvarCont[1][0] = 400;
  normvarCont[1][1] = 600;
  normvarCont[2][0] = 500;
  normvarCont[2][1] = 700;

  int tQCDVars[3];
  tQCDVars[0] = 11;
  tQCDVars[1] = 12;
  tQCDVars[2] = 13;

  double PurityProton;
  double PurityXi;

  double pp_f0;
  double pp_f1;

  //pPb
  if (system == 0) {
    PurityProton = 0.984266;  //pPb 5 TeV
    PurityXi = 0.88;  //new cuts

    pp_f0 = 0.862814;
    pp_f1 = 0.09603;

  } else {  // pp MB + HM
    PurityProton = 0.991213;
    PurityXi = 0.915;

    pp_f0 = 0.874808;
    pp_f1 = 0.0876342;  //fraction of
  }
  double ProtonPrim = pp_f0;
  double arrayPercLamProton[3] = { pp_f1 / (1. - pp_f0) * 0.8, pp_f1
      / (1. - pp_f0), pp_f1 / (1. - pp_f0) * 1.2 };  //+/- 20%

  const unsigned NumChannels_p = 4;
  double** Purities_p = new double*[3];
  double** Fraction_p = new double*[3];
  for (unsigned uVar = 0; uVar < 3; uVar++) {
    Purities_p[uVar] = new double[NumChannels_p];
    Fraction_p[uVar] = new double[NumChannels_p];

    Purities_p[uVar][0] = PurityProton;
    Purities_p[uVar][1] = PurityProton;
    Purities_p[uVar][2] = PurityProton;
    Purities_p[uVar][3] = 1. - PurityProton;

    Fraction_p[uVar][0] = ProtonPrim;
    Fraction_p[uVar][1] = (1. - ProtonPrim) * (arrayPercLamProton[uVar]);
    Fraction_p[uVar][2] = (1. - ProtonPrim) * (1. - arrayPercLamProton[uVar]);
    Fraction_p[uVar][3] = 1.;
  }

  //ratio Xi-(1530) to Xi-
  const double Xim1530_to_Xim = 0.32 * (1. / 3.);
  //ratio Xi0(1530) to Xi0 (n=neutral)
  const double Xin1530_to_Xim = 0.32 * (2. / 3.);
  const double Omegam_to_Xim = 0.1;
  const double OmegamXim_BR = 0.086;

  const unsigned NumChannels_Xim = 5;
  double** Purities_Xim = new double*[3];
  double** Fraction_Xim = new double*[3];
  for (unsigned uVar = 0; uVar < 3; uVar++) {
    Purities_Xim[uVar] = new double[NumChannels_Xim];
    Fraction_Xim[uVar] = new double[NumChannels_Xim];

    Purities_Xim[uVar][0] = PurityXi;
    Purities_Xim[uVar][1] = PurityXi;
    Purities_Xim[uVar][2] = PurityXi;
    Purities_Xim[uVar][3] = PurityXi;
    Purities_Xim[uVar][4] = 1. - PurityXi;

    //the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
    //hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
    Fraction_Xim[uVar][0] = 1. - Xim1530_to_Xim - Xin1530_to_Xim
        - Omegam_to_Xim * OmegamXim_BR;
    Fraction_Xim[uVar][1] = Xim1530_to_Xim;
    Fraction_Xim[uVar][2] = Xin1530_to_Xim;
    Fraction_Xim[uVar][3] = Omegam_to_Xim * OmegamXim_BR;
    Fraction_Xim[uVar][4] = 1.;
  }
  float ppRadii[3];
  if (system == 0) {
    ppRadii[0] = 1.141;
    ppRadii[1] = 1.427;
    ppRadii[2] = 1.434;
  } else if (system == 2) {
    ppRadii[0] = 0.76;
    ppRadii[1] = 0.79;
    ppRadii[2] = 0.82;
  } else {
    std::cout << "Radii for system " << system
              << " not implemented, exiting \n";
    return;
  }
  std::cout << ppRadii[0] << '\t' << ppRadii[1] << '\t' << ppRadii[2]
            << std::endl;
  TString CalibBaseDir;
  TString DataDir;
  if (system == 0) {
    if (NumIter > 33) {
      std::cout << "The number of variations is not so damn high (" << NumIter
                << "), exiting \n";
      return;
    }
    CalibBaseDir = "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
    DataDir = "~/cernbox/pPb/Systematics_Test/CFs";
  } else if (system == 1) {
    CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  } else if (system == 2) {
    CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  }
  TidyCats::Sources TheSource;
  if (iSource == 0) {
    TheSource = TidyCats::sGaussian;
  } else if (iSource == 1) {
    TheSource = TidyCats::sResonance;
  } else if (iSource == 2) {
    TheSource = TidyCats::sLevy;
  } else {
    std::cout << "Source does not exist! Exiting \n";
    return;
  }

  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  std::cout << "Read Resolution File \n";
//  TidyCats::pXimPot pot = TidyCats::pCoulomb;
  TidyCats::pXimPot pot;
  if (iPot == 0) {
    pot = TidyCats::pCoulomb;
  } else if (iPot == 1) {
    pot = TidyCats::pGamow;
  } else if (iPot == 2) {
    pot = TidyCats::pHALQCD;
  } else if (iPot == 3) {
    pot = TidyCats::pHALQCDGamow;
  } else if (iPot == 5) {
    pot = TidyCats::pRikkenWF;
  } else if (iPot == 6) {
    pot = TidyCats::pRikkenPot;
  } else {
    std::cout << "iPot " << iPot << " not implemented ,exiting \n";
  }

  const int binwidth = 20;
  const int rebin = 5;
  const unsigned NumMomBins_pXim = 30;

  double kMinXiP;

  const char *prefix;
  if (system == 0) {
    CATSinput->SetSigmaFileName("Sample3_MeV_compact.root");
    CATSinput->SetFixedkStarMinBin(true, 0.008);
    kMinXiP = 8.;
    prefix = "MB";
  } else if (system == 1) {
    CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
    prefix = "MB";
    kMinXiP = 0.;
  } else if (system == 2) {
    CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
//    CATSinput->SetFixedkStarMinBin(true, 0.);
    prefix = "HM";
    kMinXiP = 0.;
  }

  const double kMin_pXim = kMinXiP;
  const double kMax_pXim = kMin_pXim + binwidth * NumMomBins_pXim;

  const double NumMomBins_SideBand = 50;
  const double kMin_SideBand = kMinXiP;
  const double kMax_SideBand = 1000 + kMinXiP;

  std::cout << "kMinXiP: " << kMin_pXim << std::endl;
  std::cout << "kMax_pXim: " << kMax_pXim << std::endl;
  std::cout << "binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pXim: " << NumMomBins_pXim << std::endl;

  double FemtoRegion_pXim[3][2];
  FemtoRegion_pXim[0][0] = kMin_pXim;
  FemtoRegion_pXim[0][1] = 360;
  FemtoRegion_pXim[1][0] = kMin_pXim;
  FemtoRegion_pXim[1][1] = 420;
  FemtoRegion_pXim[2][0] = kMin_pXim;
  FemtoRegion_pXim[2][1] = 480;

  CATSinput->ReadSigmaFile();

  std::cout << "Read Sigma File \n";
  TFile* OutFile = new TFile(
      TString::Format("%s/OutFileVarpXi_%u.root", OutputDir.Data(), NumIter),
      "RECREATE");
  TList* CollOut = new TList();
  CollOut->SetOwner();
  CollOut->SetName(TString::Format("Out%u", NumIter));
//  TFile* OutGraphFile = new TFile(
//      TString::Format("%s/OutGraphFileVarpXi_%u.root", OutputDir.Data(), NumIter),
//      "recreate");
  //you save a lot of stuff in an NTuple
  TNtuple* ntResult = new TNtuple(
      "ntResult", "ntResult", "IterID:vFemReg_pXim:vFrac_pXim_pXi1530:"
      "tOut:ppRadius:AlphaLev:AlphaLevErr:varSideNorm:BLSlope:"
      "p_a:p_a_err:p_b:p_b_err:"
      "Chi2NdfGlobal:Chi2NdfLocal:pval:sigma200:sigma100:sigma150");

  Float_t ntBuffer[19];

  int vFemReg_pXim;  //which femto region we use for pXim (1 = default)
  int vFrac_pXim_pXi1530;
  int tOut;
  int ppRadius;
  int varSideNorm;
  bool HaveWeABaseLineSlope = true;

  TidyCats* tidy = new TidyCats();

  SideBandFit* side = new SideBandFit();
  side->SetRebin(rebin);
  if (system == 0) {
    side->SetSideBandFile("~/cernbox/pPb/Sidebands", "MB", "42", "43");
  } else if (system == 2) {
    side->SetSideBandFile("~/cernbox/HM13TeV/AnalysisData/Systematic", "HM",
                          "103", "104");
  }

  int uIter = 1;
//  TString HistpXiDefaultName = "hCk_ReweightedMeV_0";
  TString HistpXiDefaultName = "hCk_RebinnedMeV_1";
  TFile* inFile = TFile::Open(
      TString::Format("%s/CFpXiVariations_%u.root", DataDir.Data(), NumIter),
      "READ");
  TH1F* Prefit = (TH1F*) inFile->Get(HistpXiDefaultName.Data());
  TH1F* StoreHist = (TH1F*) Prefit->Clone("Ck_Input");
//  float p_a_prefit[4] = {0.983778, 1.00007, 1.0160844,0.9840556};
//  float p_b_prefit[4] = {0, -3.56168e-05, -1.3894000e-06, -6.9844200e-05};
  float p_a_prefit[4];
  float p_b_prefit[4];
  TF1* funct_0 = new TF1("myPol0", "pol0", 250, 600);
  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 600);
//  funct_0->SetParameter(0,p_a_prefit[0]);
//  funct_1->SetParameters(p_a_prefit[1],p_b_prefit[1]);

  Prefit->Fit(funct_0, "SNR");
  p_a_prefit[0] = funct_0->GetParameter(0);
  p_b_prefit[0] = 0;

  Prefit->Fit(funct_1, "FNR");

  p_a_prefit[1] = funct_1->GetParameter(0);
  p_b_prefit[1] = funct_1->GetParameter(1);

  gMinuit->SetErrorDef(1);  //note 4 and not 2!
  TGraph *gr12 = (TGraph*) gMinuit->Contour(40, 0, 1);

  p_a_prefit[2] = TMath::MinElement(gr12->GetN(), gr12->GetX());
  p_b_prefit[2] = gr12->Eval(p_a_prefit[2]);

  p_a_prefit[3] = TMath::MaxElement(gr12->GetN(), gr12->GetX());
  p_b_prefit[3] = gr12->Eval(p_a_prefit[3]);

  for (int i = 0; i < 4; ++i) {
    std::cout << i << " pa: " << p_a_prefit[i] << " pb: " << p_b_prefit[i]
              << std::endl;
  }
  delete Prefit;

  CATS AB_pXim;

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins_pXim, kMin_pXim,
                                 kMax_pXim, TidyCats::sGaussian);
  AB_pXim1530.SetAnaSource(1, 1);
  AB_pXim1530.KillTheCat();
  int counter = 1;
  float total = 0;

  int tOutVars = 0;
  if (pot == TidyCats::pHALQCD || pot == TidyCats::pHALQCDGamow) {
    tOutVars = 3;
    total = 324 * 3;
  } else {
    tOutVars = 1;
    total = 324;
  }
  std::vector<bool> SaveSideBands = { true, true, true };
  TCanvas* c1 = new TCanvas(TString::Format("out%u", NumIter));
  c1->SetCanvasSize(1920, 1280);
  StoreHist->SetLineWidth(3);
  StoreHist->SetLineColor(1);
  StoreHist->GetXaxis()->SetRangeUser(0, 400);
  c1->cd();
  StoreHist->DrawCopy();
  StoreHist->GetXaxis()->SetRangeUser(0, 1000);

  for (tOut = 0; tOut < tOutVars; ++tOut) {
    tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins_pXim, kMin_pXim, kMax_pXim,
                               TheSource, pot, tQCDVars[tOut]);
    AB_pXim.KillTheCat();
    for (vFemReg_pXim = 0; vFemReg_pXim < 3; ++vFemReg_pXim) {
      for (vFrac_pXim_pXi1530 = 0; vFrac_pXim_pXi1530 < 3;
          ++vFrac_pXim_pXi1530) {
        for (ppRadius = 0; ppRadius < 3; ++ppRadius) {
          for (varSideNorm = 0; varSideNorm < 3; ++varSideNorm) {
            for (int BaselineSlope = 0; BaselineSlope < 4; ++BaselineSlope) {
              if (BaselineSlope == 0) {
                HaveWeABaseLineSlope = false;  //use baseline
              } else {
                HaveWeABaseLineSlope = true;  // no baseline
              }

              std::cout
                  << "\r Processing progress: "
                  << TString::Format("%.1f %%", counter++ / total * 100.f).Data()
                  << std::flush;

              side->SetNormalizationRange(normvarCont[varSideNorm][0],
                                          normvarCont[varSideNorm][1]);
              double GaussSourceSize = ppRadii[ppRadius];
              //vary this
              side->SideBandCFs(false);
              TH1F* fitme = side->GetSideBands(5);
              double SideBandPars[4];
              side->FitSideBands(fitme, SideBandPars);
              TH1F* OliHisto_pXim = (TH1F*) inFile->Get(
                  HistpXiDefaultName.Data());

              if (!OliHisto_pXim) {
                std::cout << HistpXiDefaultName.Data() << " pXi Missing"
                          << std::endl;
                return;
              }

              const unsigned NumSourcePars = AB_pXim.GetNumSourcePars();
              DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
              DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);
              DLM_Ck* Ck_SideBand = new DLM_Ck(0, 4, NumMomBins_SideBand,
                                               kMin_SideBand, kMax_SideBand,
                                               SideBandFit::Parameterization);

              Ck_SideBand->SetPotPar(0, SideBandPars[0]);
              Ck_SideBand->SetPotPar(1, SideBandPars[1]);
              Ck_SideBand->SetPotPar(2, SideBandPars[2]);
              Ck_SideBand->SetPotPar(3, SideBandPars[3]);
              Ck_pXim->Update();
              Ck_pXim1530->Update();
              Ck_SideBand->Update();
              if (!CATSinput->GetSigmaFile(3)) {
                std::cout << "No Sigma file 3 \n";
                return;
              }
              DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
                                             CATSinput->GetSigmaFile(3));
              DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530,
              NULL);
              DLM_CkDecomposition CkDec_SideBand("pXiSideBand", 0, *Ck_SideBand,
              NULL);
              int vFrac_pp_pL = 1;
              const double lam_pXim = Purities_p[vFrac_pp_pL][0]
                  * Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][0]
                  * Fraction_Xim[0][0];
              const double lam_pXim_pXim1530 = Purities_p[vFrac_pp_pL][0]
                  * Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][1]
                  * Fraction_Xim[0][1];
              const double lam_pXim_fake = Purities_p[vFrac_pp_pL][3]
                  * Purities_Xim[0][0]
                  + Purities_p[vFrac_pp_pL][0] * Purities_Xim[0][4]
                  + Purities_p[vFrac_pp_pL][3] * Purities_Xim[0][4];

//              printf("lam_pXim = %.3f\n", lam_pXim);
//              printf("lam_pXim_pXim1530 = %.3f\n", lam_pXim_pXim1530);
//              printf("lam_pXim_fake = %.3f\n", lam_pXim_fake);
//              printf("\n");

              CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
                                         DLM_CkDecomposition::cFeedDown,
                                         &CkDec_pXim1530,
                                         CATSinput->GetResFile(3));  //from Xi-(1530)
              CkDec_pXim.AddContribution(
                  1, 1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
                  DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
              CkDec_pXim.AddContribution(2, lam_pXim_fake,
                                         DLM_CkDecomposition::cFake,
                                         &CkDec_SideBand);
              CkDec_pXim.Update();
              DLM_Fitter1* fitter = new DLM_Fitter1(1);
              fitter->SetSystem(0, *OliHisto_pXim, 1, CkDec_pXim,
                                FemtoRegion_pXim[vFemReg_pXim][0],
                                FemtoRegion_pXim[vFemReg_pXim][1], BlRegion[0],
                                BlRegion[1]);
              fitter->SetSeparateBL(0, false);  //Simultaneous BL
              if (HaveWeABaseLineSlope) {
                fitter->FixParameter("pXim", DLM_Fitter1::p_a,
                                     p_a_prefit[BaselineSlope]);
                fitter->FixParameter("pXim", DLM_Fitter1::p_b,
                                     p_b_prefit[BaselineSlope]);
              } else {
                fitter->FixParameter("pXim", DLM_Fitter1::p_a, p_a_prefit[0]);
                fitter->FixParameter("pXim", DLM_Fitter1::p_b, p_b_prefit[0]);
              }
              fitter->AddSameSource("pXim1530", "pXim", 1);
              fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);

              fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -1.);

              if (TheSource == TidyCats::sLevy) {
                fitter->SetParameter("pXim", DLM_Fitter1::p_sor1, 1.6, 1.0,
                                     2.0);
              } else if (TheSource == TidyCats::sResonance) {
                fitter->FixParameter("pXim", DLM_Fitter1::p_sor1, 2.0);
              }
              if (fitRadius) {
                fitter->SetParameter("pXim", DLM_Fitter1::p_sor0,
                                     GaussSourceSize, 0.4, 1.0);
              } else {
                fitter->FixParameter("pXim", DLM_Fitter1::p_sor0,
                                     GaussSourceSize);
              }
              fitter->GoBabyGo();
              double p_a_strong = fitter->GetParameter("pXim",
                                                       DLM_Fitter1::p_a);
              double p_a_strong_err = fitter->GetParError("pXim",
                                                          DLM_Fitter1::p_a);
              double p_b_strong = fitter->GetParameter("pXim",
                                                       DLM_Fitter1::p_b);
              double p_b_strong_err = fitter->GetParError("pXim",
                                                          DLM_Fitter1::p_b);
              double Cl_strong = fitter->GetParameter("pXim", DLM_Fitter1::p_c);
              double ChiSqStrongGlobal = fitter->GetChi2Ndf();
              double pValStrongGlobal = fitter->GetPval();

              TGraph FitResult_pXim;

              FitResult_pXim.SetName(TString::Format("pXimGraph"));
              fitter->GetFitGraph(0, FitResult_pXim);
              c1->cd();
              TGraph *pointerFitRes = new TGraph(FitResult_pXim);
              pointerFitRes->SetLineWidth(2);
              pointerFitRes->SetLineColor(kRed);
              pointerFitRes->SetMarkerStyle(24);
              pointerFitRes->SetMarkerColor(kRed);
              pointerFitRes->SetMarkerSize(1);
              pointerFitRes->Draw("CP, same");
              TGraph SideBandStrongWithLambda;
              TGraph SideBandStrongWithOutLambda;
              DLM_Histo<double>* StrongWithLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", true);
              DLM_Histo<double>* StrongWithOutLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", false);
              for (int iBin = 0; iBin < StrongWithLambda->GetNbins(); ++iBin) {
                SideBandStrongWithLambda.SetPoint(
                    iBin,
                    StrongWithLambda->GetBinCenter(0, iBin),
                    (StrongWithLambda->GetBinContent(iBin) + (1 - lam_pXim_fake)));
              }
              for (int iBin = 0; iBin < StrongWithOutLambda->GetNbins();
                  ++iBin) {
                SideBandStrongWithOutLambda.SetPoint(
                    iBin, StrongWithOutLambda->GetBinCenter(0, iBin),
                    StrongWithOutLambda->GetBinContent(iBin));
              }

              double Chi2_pXim = 0;
              unsigned EffNumBins_pXim = 0;

              double Chi2_pXim_kSm100 = 0;
              unsigned EffNumBins_pXim_kSm100 = 0;

              double Chi2_pXim_kSm150 = 0;
              unsigned EffNumBins_pXim_kSm150 = 0;

              int maxkStarBin = OliHisto_pXim->FindBin(200.);
              for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {

                double mom = AB_pXim.GetMomentum(uBin);
                double dataY;
                double dataErr;
                double theoryX;
                double theoryY;

                if (mom > FemtoRegion_pXim[vFemReg_pXim][1])
                  continue;

                FitResult_pXim.GetPoint(uBin, theoryX, theoryY);
                if (mom != theoryX) {
                  std::cout << mom << '\t' << theoryX << std::endl;
                  printf("  PROBLEM pXi!\n");
                }
                dataY = OliHisto_pXim->GetBinContent(uBin + 1);
                dataErr = OliHisto_pXim->GetBinError(uBin + 1);
                Chi2_pXim += (dataY - theoryY) * (dataY - theoryY)
                    / (dataErr * dataErr);
                EffNumBins_pXim++;
                if (mom < 100) {
                  Chi2_pXim_kSm100 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  EffNumBins_pXim_kSm100++;
                }
                if (mom < 150) {
                  Chi2_pXim_kSm150 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  EffNumBins_pXim_kSm150++;
                }
              }

              double pvalXi = TMath::Prob(Chi2_pXim, round(EffNumBins_pXim));
              double nSigmaXi = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXi);

              double pvalXi_kSm100 = TMath::Prob(Chi2_pXim_kSm100,
                                                 round(EffNumBins_pXim_kSm100));
              double nSigmaXi_kSm100 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXi_kSm100);
              double pvalXi_kSm150 = TMath::Prob(Chi2_pXim_kSm150,
                                                 round(EffNumBins_pXim_kSm150));
              double nSigmaXi_kSm150 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXi_kSm150);

              ntBuffer[0] = uIter;
              ntBuffer[1] = vFemReg_pXim;
              ntBuffer[2] = vFrac_pXim_pXi1530;
              ntBuffer[3] = tOut;
              ntBuffer[4] = GaussSourceSize;
              ntBuffer[5] = fitter->GetParameter("pXim", DLM_Fitter1::p_sor1);
              ntBuffer[6] = fitter->GetParError("pXim", DLM_Fitter1::p_sor1);
              ntBuffer[7] = varSideNorm;
              ntBuffer[8] = (float) BaselineSlope;
              ntBuffer[9] = p_a_strong;
              ntBuffer[10] = p_a_strong_err;
              ntBuffer[11] = p_b_strong;
              ntBuffer[12] = p_b_strong_err;
              ntBuffer[13] = ChiSqStrongGlobal;
              ntBuffer[14] = Chi2_pXim / double(EffNumBins_pXim);
              ntBuffer[15] = pvalXi;
              ntBuffer[16] = nSigmaXi;
              ntBuffer[17] = nSigmaXi_kSm100;
              ntBuffer[18] = nSigmaXi_kSm150;
              ntResult->Fill(ntBuffer);
              TList* outList = new TList();
              outList->SetOwner();
              outList->SetName(
                  TString::Format("Graph_Var_%u_iter_%u", NumIter, uIter));
              outList->Add(pointerFitRes);
              CollOut->Add(outList);

              delete fitter;
              delete Ck_pXim;
              delete Ck_pXim1530;
              delete Ck_SideBand;
              delete StrongWithOutLambda;
              delete StrongWithLambda;
              uIter++;
            }
          }
        }
      }
    }
  }
  std::cout << "\n";
  OutFile->cd();
  c1->Write(TString::Format("%s", c1->GetName()));
  ntResult->Write();
  CollOut->Write(CollOut->GetName(), 1);
  OutFile->Close();
  return;
}

int main(int argc, char *argv[]) {
  const char* addon = (argv[7]) ? argv[7] : "";
  GetXiForRadius(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                 argv[5]);
  return 0;
}

//  OutGraphFile->Close();
//  TList* InputList = new TList();
//  InputList->SetOwner();
//  InputList->SetName(TString::Format("Inputs_%u", NumIter));
//
//  gr12->SetName("ContourBaseline");
//  InputList->Add(gr12);
//  InputList->Add(funct_0);
//  InputList->Add(funct_1);
//  InputList->Add(StoreHist);
//  InputList->Write(InputList->GetName(), 1);
//
//              outList->Add(OliHisto_pXim);
//              FitResult_pXim.SetLineWidth(2);
//              FitResult_pXim.SetLineColor(kRed);
//              FitResult_pXim.SetMarkerStyle(24);
//              FitResult_pXim.SetMarkerColor(kRed);
//              FitResult_pXim.SetMarkerSize(1);
//              outList->Write(outList->GetName(),1);
//              outList->Write(outList->GetName(), 1);
//              if (NumIter == 1 && SaveSideBands[varSideNorm]) {
//                TList* SideList;
//                SideList = new TList();
//                SideList->SetOwner();
//                SideList->SetName("SideList");
//                SideBandStrongWithLambda.SetName("SideBandStrongWithLambda");
//                SideList->Add(&SideBandStrongWithLambda);
//                SideBandStrongWithOutLambda.SetName(
//                    "SideBandStrongWithOutLambda");
//                SideList->Add(&SideBandStrongWithOutLambda);
//                SideList->Add(fitme);
//                OutFile->cd();
//                SideList->Write(SideList->GetName(), 1);
//
//                SaveSideBands[varSideNorm] = false;
//              }

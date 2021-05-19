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
#include "CATSLambdaParam.h"
#include "TObject.h"
#include "TDirectoryFile.h"
#include "TROOT.h"

void GetXiForRadius(const unsigned& NumIter, int system, int iPot, int iSource,
                    TString DataFile, TString HistpXiDefaultName,
                    TString OutputDir) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  //System: 0 = pPb, 1 = pp MB, 2 = pp HM
  //Potential: 0 = Coulomb, 1 = Gamow, 2 = HAL + Coulomb, 3 = HAL + Gamow, 5 = RikkenWF , 6 = RikkenPot
  //Source: 0 = Gauss, 1 = Resonances, 2 = Levy
  //Iteration number corresponds here to the systematic cut variation being fit.
//  bool drawSource = (NumIter == 0);
  bool fitRadius = false;
//  TString HistpXiDefaultName = "hCk_ReweightedMeV_0";
  TFile* inFile = TFile::Open(DataFile.Data(), "READ");
  TH1F* StoreHist = (TH1F*) inFile->Get(HistpXiDefaultName.Data());
  if (!StoreHist) {
    Error("GetXiForRadius", "No Histogram Loaded!\n ");
    return;
  }
  TH1F* Prefit = (TH1F*) StoreHist->Clone("CkPrefit");

  const int binwidth = 20;
  const int rebin = 5;
  const unsigned NumMomBins_pXim = 30;

//  double kMinXiP;

  const double kMin = Prefit->GetXaxis()->GetXmin();
  const double kMax = kMin + binwidth * NumMomBins_pXim;

  const double NumMomBins_SideBand = 50;
  const double kMax_SideBand = 1000 + kMin;

  std::cout << "kMinXiP: " << kMin << std::endl;
  std::cout << "kMax: " << kMax << std::endl;
  std::cout << "binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pXim: " << NumMomBins_pXim << std::endl;

  double FemtoRegion_pXim[3];
  FemtoRegion_pXim[0] = 360;
  FemtoRegion_pXim[1] = 420;
  FemtoRegion_pXim[2] = 480;

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

  float ppRadii[3];

  TString CalibBaseDir = "";
  TString SigmaFileName = "";
  double PurityProton, PrimProton, SecLamProton;
  double PurityXi;

  std::cout << "SYSTEM: " << system << std::endl;
  if (system == 0) {
    if (NumIter > 34) {
      std::cout << "The number of variations is not so damn high (" << NumIter
                << "), exiting \n";
      return;
    }
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
    SigmaFileName += "Sample3_MeV_compact.root";

    ppRadii[0] = 1.141;
    ppRadii[1] = 1.427;
    ppRadii[2] = 1.434;

    PurityProton = 0.984266;  //pPb 5 TeV
    PrimProton = 0.862814;
    SecLamProton = 0.09603;

    PurityXi = 0.88;  //new cuts
  } else if (system == 2) {  // pp HM
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
    SigmaFileName += "Sample6_MeV_compact.root";
//    mT integrated
    ppRadii[0] = 0.76;
    ppRadii[1] = 0.79;
    ppRadii[2] = 0.82;
    //mT 0-1.74
//    ppRadii[0] = 0.87;
//    ppRadii[1] = 0.9;
//    ppRadii[2] = 0.93;
    //mT 1.74-4.5
//    ppRadii[0] = 0.69;
//    ppRadii[1] = 0.72;
//    ppRadii[2] = 0.76;

    PurityProton = 0.9943;
    PrimProton = 0.873;
    SecLamProton = 0.089;  //Fraction of Lambdas

    PurityXi = 0.915;
  } else {
    std::cout << "System " << system << " not supported, extiting \n";
    return;
  }

  // Xi01530 Production: dN/dy = 2.6e-3 (https://link.springer.com/content/pdf/10.1140%2Fepjc%2Fs10052-014-3191-x.pdf)
  // Xim1530 Production = Xi01530 Production
  // Xim Production: dN/dy = 5.3e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // -> Production Ratio ~ 1/2

  const double Xi01530XimProdFraction = 1 / 2.;
  const double Xim1530XimProdFraction = 1 / 2.;

  // 2/3 of Xi0(1530) decays via Xi- + pi+ (Isospin considerations)
  const double Xi01530Xim_BR = 2 / 3.;
  // 1/3 of Xi-(1530) decays via Xi- + pi0 (Isospin considerations)
  const double Xim1530Xim_BR = 1 / 3.;

  // Omega production: dN/dy = 0.67e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // Xim Production: dN/dy = 5.3e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // -> Production Ratio ~ 1/10
  const double OmegamXimProdFraction = 1 / 10.;
  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %

  // Produce N Xi's -> Produce:
  // 1 ) N* 1/10 Omegas -> See N* 1/10 * 8.6% more Xi's
  // 2)  N* 1/2 Xi0_1530 -> See N*1/2*2/3 = N* 1/3 more Xi's
  // 3)  N* 1/2 Xim_1530 -> See N*1/2*1/3 = N* 1/6 more Xi's
  // Total Sample:  N(1+0.0086+1/3+1/6) ->
  // Primary Fraction = N / N(1+0.0086+1/3+1/6)
  // Secondary Omegas = N*0.0086  / N(1+0.0086+1/3+1/6)
  // etc.

  std::vector<double> Variation = { 0.8, 1.0, 1.2 };
  Particle Proton[3];  // 1) variation of the Secondary Comp.
  Particle Xi[3][3];  //1) variation of dN/dy Omega 2) variation of dN/dy Xi1530

  int iVar1 = 0;
  for (auto it : Variation) {
    double SecFracSigma = 1. - PrimProton - it * SecLamProton;
    Proton[iVar1] = Particle(PurityProton, PrimProton, { it * SecLamProton,
                                 SecFracSigma });
    int iVar2 = 0;
    for (auto itXim : Variation) {
      double XiNormalization = 1 + it * OmegamXimProdFraction * OmegamXim_BR
          + itXim * Xi01530XimProdFraction * Xi01530Xim_BR
          + itXim * Xim1530XimProdFraction * Xim1530Xim_BR;
      double SecOmegaXim = it * OmegamXimProdFraction * OmegamXim_BR
          / (double) XiNormalization;
      double SecXi01530Xim = itXim * Xi01530XimProdFraction * Xi01530Xim_BR
          / (double) XiNormalization;
      double SecXim1530Xim = itXim * Xim1530XimProdFraction * Xim1530Xim_BR
          / (double) XiNormalization;
      double PrimXim = 1. / (double) XiNormalization;
      Xi[iVar1][iVar2] = Particle(PurityXi, PrimXim, { SecOmegaXim,
                                      SecXi01530Xim, SecXim1530Xim });
      iVar2++;
    }
    iVar1++;
  }

  std::cout << ppRadii[0] << '\t' << ppRadii[1] << '\t' << ppRadii[2]
            << std::endl;
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
  CATSinput->SetSigmaFileName(SigmaFileName.Data());
  CATSinput->ReadSigmaFile();

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
    return;
  }

  TFile* OutFile = new TFile(
      TString::Format("%s/OutFileVarpXi_%u.root", OutputDir.Data(), NumIter),
      "RECREATE");
  TDirectoryFile* CollOut = new TDirectoryFile(
      TString::Format("Out%u", NumIter), TString::Format("Out%u", NumIter));
//  CollOut->SetOwner();
//  CollOut->SetName();
  //you save a lot of stuff in an NTuple
  TNtuple* ntResult = new TNtuple(
      "ntResult", "ntResult",
      "NumIter:IterID:FemtoRegion:lampXi:lampXi1530:lampXiFake:"
      "tOut:ppRadius:AlphaLev:AlphaLevErr:varSideNorm:BLSlope:"
      "p_a:p_a_err:p_b:p_b_err:p_c:p_c_err:"
      "Chi2NdfGlobal:chisqPerndf:pval:sigma200:sigma100:sigma150");

  Float_t ntBuffer[23];

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
    side->SetSideBandFile("~/cernbox/pPb/Sidebands_redo", "MB", "42", "43");
  } else if (system == 2) {
    side->SetSideBandFile("~/cernbox/HM13TeV/AnalysisData/latestSystematic",
                          "HM", "103", "104");
  }

  int uIter = 1;

  float p_a_prefit[4];
  float p_b_prefit[4];
  TF1* funct_0 = new TF1("myPol0", "pol0", 250, 450);
  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 450);

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
  funct_0->SetRange(0, 600);
  funct_1->SetRange(0, 600);
//  delete Prefit;

  CATS AB_pXim;

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins_pXim, kMin, kMax,
                                 TidyCats::sGaussian);
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

  CollOut->Add(c1);
  CollOut->Add(funct_0);
  CollOut->Add(funct_1);
  CollOut->Add(Prefit);
  c1->cd();
  StoreHist->DrawCopy();
  StoreHist->GetXaxis()->SetRangeUser(0, 1000);

  for (tOut = 0; tOut < tOutVars; ++tOut) {
    tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins_pXim, kMin, kMax, TheSource,
                               pot, tQCDVars[tOut]);
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
              DLM_Ck* Ck_SideBand = new DLM_Ck(0, 4, NumMomBins_SideBand, kMin,
                                               kMax_SideBand,
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

              CATSLambdaParam pXiLam(Proton[1], Xi[1][vFrac_pXim_pXi1530]);
              const double lam_pXim = pXiLam.GetLambdaParam(
                  CATSLambdaParam::Primary, CATSLambdaParam::Primary);
              ;
              const double lam_pXim_pXim1530 = pXiLam.GetLambdaParam(
                  CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
              const double lam_pXim_fake = pXiLam.GetLambdaParam(
                  CATSLambdaParam::Primary, CATSLambdaParam::Fake)
                  + pXiLam.GetLambdaParam(CATSLambdaParam::Fake,
                                          CATSLambdaParam::Primary)
                  + pXiLam.GetLambdaParam(CATSLambdaParam::Fake,
                                          CATSLambdaParam::Fake);
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
              fitter->SetSystem(0, *OliHisto_pXim, 1, CkDec_pXim, kMin,
                                FemtoRegion_pXim[vFemReg_pXim],
                                FemtoRegion_pXim[vFemReg_pXim],
                                FemtoRegion_pXim[vFemReg_pXim]);
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

              FitResult_pXim.SetName(TString::Format("Graph_Var_%u_Iter_%u", NumIter, uIter));
              fitter->GetFitGraph(0, FitResult_pXim);
              c1->cd();
              TGraph *pointerFitRes = new TGraph(FitResult_pXim);
              pointerFitRes->SetName(
                  TString::Format("Graph_Var_%u_Iter_%u", NumIter, uIter));
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

                if (mom > FemtoRegion_pXim[vFemReg_pXim])
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

              ntBuffer[0] = NumIter;
              ntBuffer[1] = uIter;
              ntBuffer[2] = FemtoRegion_pXim[vFemReg_pXim];
              ntBuffer[3] = lam_pXim;
              ntBuffer[4] = lam_pXim_pXim1530;
              ntBuffer[5] = lam_pXim_fake;
              ntBuffer[6] = tOut;
              ntBuffer[7] = GaussSourceSize;
              ntBuffer[8] = fitter->GetParameter("pXim", DLM_Fitter1::p_sor1);
              ntBuffer[9] = fitter->GetParError("pXim", DLM_Fitter1::p_sor1);
              ntBuffer[10] = varSideNorm;
              ntBuffer[11] = (float) BaselineSlope;
              ntBuffer[12] = p_a_strong;
              ntBuffer[13] = p_a_strong_err;
              ntBuffer[14] = p_b_strong;
              ntBuffer[15] = p_b_strong_err;
              ntBuffer[16] = 0;
              ntBuffer[17] = 0;
              ntBuffer[18] = ChiSqStrongGlobal;
              ntBuffer[19] = Chi2_pXim / double(EffNumBins_pXim);
              ntBuffer[20] = pvalXi;
              ntBuffer[21] = nSigmaXi;
              ntBuffer[22] = nSigmaXi_kSm100;
              ntBuffer[23] = nSigmaXi_kSm150;
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
  ntResult->Write();
  CollOut->Write(TString::Format("Out%u", NumIter), TObject::kSingleKey);
  OutFile->Close();
  delete tidy;
  return;
}

int main(int argc, char *argv[]) {
  GetXiForRadius(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                 argv[5], argv[6], argv[7]);
  return 0;
}

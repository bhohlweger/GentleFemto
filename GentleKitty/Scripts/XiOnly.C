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
void GetXiForRadius(const unsigned& NumIter, TString InputDir, TString ppFile,
                    TString OutputDir) {
  bool FAST_PLOT = true;
  TRandom3 rangen(0);
  const int binwidth = 20;
  const unsigned NumMomBins_pXim = 30;

  double kMinXiP = 8.;
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

  double PurityProton = 0.984265;  //pPb 5 TeV
  double PurityXi = 0.88;  //new cuts

  double pp_f0 = 0.862814;
  double pp_f1 = 0.09603;

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
  TFile *ppVars = TFile::Open(ppFile, "READ");
  TNtuple* ppTuple = (TNtuple*) ppVars->Get("outTuple");
  float ppRadii[3];
  ppTuple->SetBranchAddress("Rad_pp", &ppRadii[0]);
  ppTuple->SetBranchAddress("RadLow_pXi", &ppRadii[1]);
  ppTuple->SetBranchAddress("RadUp_pXi", &ppRadii[2]);
  ppTuple->GetEntry(0);
  std::cout << ppRadii[0] << '\t' << ppRadii[1] << '\t' << ppRadii[2]
            << std::endl;
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
  TString ResMatrixFileName = TString::Format("%s/run2_decay_matrices_old.root",
                                              CalibBaseDir.Data());
  TString SigmaMatrixFileName = TString::Format("%s/Sample3_MeV_compact.root",
                                                CalibBaseDir.Data());

  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample3_MeV_compact.root");
  CATSinput->ReadSigmaFile();

  CATSinput->ReadCorrelationFile(InputDir.Data());

  TFile* OutFile = new TFile(
      TString::Format("%s/OutFilepXi_%s_Iter%u.root", OutputDir.Data(),
                      "CutVarAdd", NumIter),
      "recreate");
  //you save a lot of stuff in an NTuple
  TNtuple* ntResult = new TNtuple(
      "ntResult", "ntResult", "IterID:vFemReg_pXim:vFrac_pXim_pXi1530:"
      "tOut:ppRadius:varSideNorm:BLSlope:"
      "p_a:p_a_err:p_b:p_b_err:"
      "Chi2NdfGlobal:Chi2NdfLocal:pval:sigma:"
      "p_a_coulomb:p_a_err_coulomb:p_b_coulomb:p_b_err_coulomb"
      ":Chi2NdfGlobal_coulomb:Chi2NdfLocal_coulomb:"
      "pval_coulomb:sigma_coulomb");

  Float_t ntBuffer[23];

  int vFemReg_pXim;  //which femto region we use for pXim (1 = default)
  int vFrac_pXim_pXi1530;
  int tOut;
  int ppRadius;
  int varSideNorm;
  bool HaveWeABaseLineSlope = true;

  TidyCats* tidy = new TidyCats();

  SideBandFit* side = new SideBandFit();
  side->SetRebin(5);
  side->SetSideBandFile("/home/hohlweger/cernbox/pPb/Sidebands", "42", "43");
  int uIter = 1;

  CATSinput->ObtainCFs(5, 240, 340);
  TString HistpXiName = "hCk_ReweightedpXiMeV_0";

  float p_a_prefit[4];
  float p_b_prefit[4];
  TH1F* Prefit = CATSinput->GetCF("pXi", HistpXiName.Data());
  TF1* funct_0 = new TF1("myPol0", "pol0", 250, 600);
  Prefit->Fit(funct_0, "FSNRM");
  p_a_prefit[0] = funct_0->GetParameter(0);
  p_b_prefit[0] = 0;

  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 600);
  Prefit->Fit(funct_1, "FSNRM");
//  TVirtualFitter::SetDefaultFitter("Minuit");
  gMinuit->SetErrorDef(4); //note 4 and not 2!

  p_a_prefit[1] = funct_1->GetParameter(0);

  p_b_prefit[1] = funct_1->GetParameter(1);

  TGraph *gr12 = (TGraph*)gMinuit->Contour(40,0,1);

  p_a_prefit[2] = TMath::MinElement(gr12->GetN(),gr12->GetX());
  p_b_prefit[2] = gr12->Eval(p_a_prefit[2]);

  p_a_prefit[3] = TMath::MaxElement(gr12->GetN(),gr12->GetX());
  p_b_prefit[3] = gr12->Eval(p_a_prefit[3]);

  for (int i = 0; i<4; ++i) {
    std::cout << i << " pa: " << p_a_prefit[i] << " pb: " << p_b_prefit[i] << std::endl;
  }

  delete Prefit;
//  delete funct_0;
//  delete funct_1;
//
  for (vFemReg_pXim = 0; vFemReg_pXim < 3; ++vFemReg_pXim) {
    for (vFrac_pXim_pXi1530 = 0; vFrac_pXim_pXi1530 < 3; ++vFrac_pXim_pXi1530) {
      for (tOut = 0; tOut < 3; ++tOut) {
        for (ppRadius = 0; ppRadius < 3; ++ppRadius) {
          for (varSideNorm = 0; varSideNorm < 3; ++varSideNorm) {
            for (int BaselineSlope = 0; BaselineSlope < 4; ++BaselineSlope) {
              if (BaselineSlope == 0) {
                HaveWeABaseLineSlope = false;  //use baseline
              } else {
                HaveWeABaseLineSlope = true;  // no baseline
              }
              if (NumIter == 0) {  //in that case run the default
                vFemReg_pXim = 1;
                vFrac_pXim_pXi1530 = 1;
                tOut = 1;
                ppRadius = 0;
                varSideNorm = 1;
                HaveWeABaseLineSlope = false;
              }
              side->SetNormalizationRange(normvarCont[varSideNorm][0],
                                          normvarCont[varSideNorm][1]);
              double GaussSourceSize = ppRadii[ppRadius];
              //vary this
              side->SideBandCFs(false);
              TH1F* fitme = side->GetSideBands(5);
              double SideBandPars[4];
              side->FitSideBands(fitme, SideBandPars);
              double Pars_pXi[6] = { 0, 0, 0, GaussSourceSize * 1.2,
                  GaussSourceSize / 1.2, 0.5 };
              CATS AB_pXim;
              tidy->GetCatsProtonXiMinus(&AB_pXim, GaussSourceSize, Pars_pXi,
                                         NumMomBins_pXim, kMin_pXim, kMax_pXim,
                                         true, tQCDVars[tOut]);
              AB_pXim.KillTheCat();

              double Pars_pXim1530[6] = { 0, 0, 0, GaussSourceSize * 1.2,
                  GaussSourceSize / 1.2, 0.5 };
              CATS AB_pXim1530;
              tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, GaussSourceSize,
                                             Pars_pXim1530, NumMomBins_pXim,
                                             kMin_pXim, kMax_pXim);
              AB_pXim1530.KillTheCat();

//              TString HistpXiName = "hCk_ReweightedpXiMeV_0";

              TH1F* OliHisto_pXim = CATSinput->GetCF("pXi", HistpXiName.Data());
              if (!OliHisto_pXim)
                std::cout << HistpXiName.Data() << " pXi Missing" << std::endl;

              TH1F *OliHisto_pXimFornSigma = nullptr;
              OliHisto_pXimFornSigma = (TH1F*) OliHisto_pXim->Clone(
                  "pXiForNSigma");
              CATSinput->AddSystematics("C2totalsysPXi.root", OliHisto_pXim);

              const unsigned NumSourcePars = 1;
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
//              const double lam_pXim = Purities_p[vFrac_pXim_pXi1530][0]
//                  * Fraction_p[vFrac_pXim_pXi1530][0] * Purities_Xim[0][0]
//                  * Fraction_Xim[0][0];
//              const double lam_pXim_pXim1530 = Purities_p[vFrac_pXim_pXi1530][0]
//                  * Fraction_p[vFrac_pXim_pXi1530][0] * Purities_Xim[0][1]
//                  * Fraction_Xim[0][1];
//              const double lam_pXim_fake = Purities_p[vFrac_pXim_pXi1530][3]
//                  * Purities_Xim[0][0]
//                  + Purities_p[vFrac_pXim_pXi1530][0] * Purities_Xim[0][4]
//                  + Purities_p[vFrac_pXim_pXi1530][3] * Purities_Xim[0][4];
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

              printf("lam_pXim = %.3f\n", lam_pXim);
              printf("lam_pXim_pXim1530 = %.3f\n", lam_pXim_pXim1530);
              printf("lam_pXim_fake = %.3f\n", lam_pXim_fake);
              printf("\n");

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
                fitter->FixParameter("pXim", DLM_Fitter1::p_a, p_a_prefit[BaselineSlope]);
                fitter->FixParameter("pXim", DLM_Fitter1::p_b, p_b_prefit[BaselineSlope]);
                std::cout << "Fitting ranges for BL set \n";
              } else {
                fitter->FixParameter("pXim", DLM_Fitter1::p_a, p_a_prefit[0]);
                fitter->FixParameter("pXim", DLM_Fitter1::p_b, p_b_prefit[0]);;
              }
              fitter->AddSameSource("pXim1530", "pXim", 1);
              fitter->AddSameSource("pXiSideBand", "pXim", 1);
              //Global Fit default
              //Fit BL & Normalization
              fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);

              fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -1.);

              fitter->FixParameter("pXim", DLM_Fitter1::p_sor0,
                                   GaussSourceSize);

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

              std::cout << "paStrongCoulomb = " << p_a_strong << " pm " << p_a_strong_err <<std::endl;
              std::cout << "pbStrongCoulomb = " << p_b_strong << " pm " << p_b_strong_err <<std::endl;

              TGraph FitResult_pXim;

              FitResult_pXim.SetName(TString::Format("pXimGraph"));
              fitter->GetFitGraph(0, FitResult_pXim);

              TGraph SideBandStrongWithLambda;
              TGraph SideBandStrongWithOutLambda;
              CATShisto<double>* StrongWithLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", true);
              CATShisto<double>* StrongWithOutLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", false);
              for (int iBin = 0; iBin < StrongWithLambda->GetNbins(); ++iBin) {
                SideBandStrongWithLambda.SetPoint(
                    iBin,
                    StrongWithLambda->GetBinCenter(iBin),
                    (StrongWithLambda->GetBinContent(iBin) + (1 - lam_pXim_fake)));
              }
              for (int iBin = 0; iBin < StrongWithOutLambda->GetNbins();
                  ++iBin) {
                SideBandStrongWithOutLambda.SetPoint(
                    iBin, StrongWithOutLambda->GetBinCenter(iBin),
                    StrongWithOutLambda->GetBinContent(iBin));
              }

              double Chi2 = fitter->GetChi2();
              unsigned NDF = fitter->GetNdf();

              double Chi2_pXim = 0;
              unsigned EffNumBins_pXim = 0;
              double Chi2_pXim_exclPeak = 0;
              unsigned EffNumBins_pXim_exclPeak = 0;
              double Chi2_pXim_kSm60 = 0;
              unsigned EffNumBins_pXim_kSm60 = 0;
              int maxkStarBin = OliHisto_pXimFornSigma->FindBin(200);
              for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {

                double mom = AB_pXim.GetMomentum(uBin);
                //double dataX;
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
                dataY = OliHisto_pXimFornSigma->GetBinContent(uBin + 1);
                dataErr = OliHisto_pXimFornSigma->GetBinError(uBin + 1);
                Chi2_pXim += (dataY - theoryY) * (dataY - theoryY)
                    / (dataErr * dataErr);
                EffNumBins_pXim++;
                if (mom < 60) {
                  Chi2_pXim_kSm60 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  EffNumBins_pXim_kSm60++;
                }
                if (!(mom > 60 && mom < 100)) {
                  Chi2_pXim_exclPeak += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  EffNumBins_pXim_exclPeak++;
                }
              }
              EffNumBins_pXim--;
              EffNumBins_pXim_kSm60--;
              EffNumBins_pXim_exclPeak--;
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //++++++++++++++++Colomb Only++++++++++++++++++++++++++++++++
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

              printf("AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));
              AB_pXim.RemoveShortRangePotential(0, 0);
              AB_pXim.RemoveShortRangePotential(1, 0);
              AB_pXim.RemoveShortRangePotential(2, 0);
              AB_pXim.RemoveShortRangePotential(3, 0);
              AB_pXim.KillTheCat(CATS::kPotentialChanged);
              printf("NEW AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));

              CkDec_pXim.Update(true);
              fitter->GoBabyGo();

              TGraph FitResult_pXim_COULOMB;
              FitResult_pXim_COULOMB.SetName(
                  TString::Format("pXimGraph_COULOMB"));
              fitter->GetFitGraph(0, FitResult_pXim_COULOMB);
              double p_a_coulomb = fitter->GetParameter("pXim",
                                                        DLM_Fitter1::p_a);
              double p_a_coulomb_err = fitter->GetParError("pXim",
                                                           DLM_Fitter1::p_a);
              double p_b_coulomb = fitter->GetParameter("pXim",
                                                        DLM_Fitter1::p_b);
              double p_b_coulomb_err = fitter->GetParError("pXim",
                                                           DLM_Fitter1::p_b);
              double Cl_coulomb = fitter->GetParameter("pXim",
                                                       DLM_Fitter1::p_c);
              double ChiSqCoulombGlobal = fitter->GetChi2Ndf();
              double pValCoulombGlobal = fitter->GetPval();
              std::cout << "paCoulomb = " << p_a_coulomb << " pm " << p_a_coulomb_err <<std::endl;
              std::cout << "pbCoulomb = " << p_b_coulomb << " pm " << p_b_coulomb_err <<std::endl;
              TGraph SideBandCoulombWithLambda;
              TGraph SideBandCoulombWithOutLambda;
              CATShisto<double>* CoulombWithLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", true);
              CATShisto<double>* CoulombWithOutLambda = CkDec_pXim
                  .GetChildContribution("pXiSideBand", false);
              for (int iBin = 0; iBin < CoulombWithLambda->GetNbins(); ++iBin) {
                SideBandCoulombWithLambda.SetPoint(
                    iBin,
                    CoulombWithLambda->GetBinCenter(iBin),
                    CoulombWithLambda->GetBinContent(iBin)
                        + (1 - lam_pXim_fake));
              }
              for (int iBin = 0; iBin < CoulombWithOutLambda->GetNbins();
                  ++iBin) {
                SideBandCoulombWithOutLambda.SetPoint(
                    iBin, CoulombWithOutLambda->GetBinCenter(iBin),
                    CoulombWithOutLambda->GetBinContent(iBin));
              }

              double Chi2_pXim_COULOMB = 0;
              unsigned EffNumBins_pXim_COULOMB = 0;
              double Chi2_pXim_COULOMB_exclPeak = 0;
              unsigned EffNumBins_pXim_COULOMB_exclPeak = 0;
              double Chi2_pXim_COULOMB_kSm60 = 0;
              unsigned EffNumBins_pXim_COULOMB_kSm60 = 0;

              for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {

                double mom = AB_pXim.GetMomentum(uBin);
                double dataY;
                double dataErr;
                double theoryX;
                double theoryY;

                if (mom > FemtoRegion_pXim[vFemReg_pXim][1])
                  continue;
                FitResult_pXim_COULOMB.GetPoint(uBin, theoryX, theoryY);
                if (mom != theoryX) {
                  std::cout << mom << '\t' << theoryX << std::endl;
                  printf("  PROBLEM pXi!\n");
                }
                dataY = OliHisto_pXimFornSigma->GetBinContent(uBin + 1);
                dataErr = OliHisto_pXimFornSigma->GetBinError(uBin + 1);

                Chi2_pXim_COULOMB += (dataY - theoryY) * (dataY - theoryY)
                    / (dataErr * dataErr);
                EffNumBins_pXim_COULOMB++;
                if (mom < 60) {
                  Chi2_pXim_COULOMB_kSm60 += (dataY - theoryY)
                      * (dataY - theoryY) / (dataErr * dataErr);
                  EffNumBins_pXim_COULOMB_kSm60++;
                }
                if (!(mom > 60 && mom < 100)) {
                  Chi2_pXim_COULOMB_exclPeak += (dataY - theoryY)
                      * (dataY - theoryY) / (dataErr * dataErr);
                  EffNumBins_pXim_COULOMB_exclPeak++;
                }
              }
              EffNumBins_pXim_COULOMB--;
              EffNumBins_pXim_COULOMB_kSm60--;
              EffNumBins_pXim_COULOMB_exclPeak--;

              double pvalXi = TMath::Prob(Chi2_pXim, round(EffNumBins_pXim));
              double nSigmaXi = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXi);

              double pvalXi_kSm60 = TMath::Prob(Chi2_pXim_kSm60,
                                                round(EffNumBins_pXim_kSm60));
              double nSigmaXi_kSm60 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXi_kSm60);

              double pvalXi_exclPeak = TMath::Prob(
                  Chi2_pXim_exclPeak, round(EffNumBins_pXim_exclPeak));
              double nSigmaXi_exclPeak = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXi_exclPeak);

              double pvalXiCoulomb = TMath::Prob(
                  Chi2_pXim_COULOMB, round(EffNumBins_pXim_COULOMB));
              double nSigmaXiCoulomb = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXiCoulomb);

              double pvalXiCoulomb_kSm60 = TMath::Prob(
                  Chi2_pXim_COULOMB_kSm60,
                  round(EffNumBins_pXim_COULOMB_kSm60));
              double nSigmaXiCoulomb_kSm60 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXiCoulomb_kSm60);

              double pvalXiCoulomb_exclPeak = TMath::Prob(
                  Chi2_pXim_COULOMB_exclPeak,
                  round(EffNumBins_pXim_COULOMB_exclPeak));
              double nSigmaXiCoulomb_exclPeak = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalXiCoulomb_exclPeak);

              ntBuffer[0] = uIter;
              ntBuffer[1] = vFemReg_pXim;
              ntBuffer[2] = vFrac_pXim_pXi1530;
              ntBuffer[3] = tOut;
              ntBuffer[4] = GaussSourceSize;
              ntBuffer[5] = varSideNorm;
              ntBuffer[6] = (float) BaselineSlope;
              ntBuffer[7] = p_a_strong;
              ntBuffer[8] = p_a_strong_err;
              ntBuffer[9] = p_b_strong;
              ntBuffer[10] = p_b_strong_err;
              ntBuffer[11] = ChiSqStrongGlobal;
              ntBuffer[12] = Chi2_pXim / double(EffNumBins_pXim);
              ntBuffer[13] = pvalXi;
              ntBuffer[14] = nSigmaXi;
              ntBuffer[15] = p_a_coulomb;
              ntBuffer[16] = p_a_coulomb_err;
              ntBuffer[17] = p_b_coulomb;
              ntBuffer[18] = p_b_coulomb_err;
              ntBuffer[19] = ChiSqCoulombGlobal;
              ntBuffer[20] = Chi2_pXim_COULOMB / double(EffNumBins_pXim);
              ntBuffer[21] = pvalXiCoulomb;
              ntBuffer[22] = nSigmaXiCoulomb;
              ntResult->Fill(ntBuffer);

              TString outFileName = TString::Format(
                  "%s/GraphFile_pXi_iter%u_Var%u.root", OutputDir.Data(),
                  NumIter, uIter);
              TFile* file = TFile::Open(outFileName.Data(), "RECREATE");
              file->cd();
              FitResult_pXim_COULOMB.Write();
              FitResult_pXim.Write();
              //we need to add the side band stuff here
              SideBandStrongWithLambda.SetName("SideBandStrongWithLambda");
              SideBandStrongWithLambda.Write();
              SideBandStrongWithOutLambda.SetName(
                  "SideBandStrongWithOutLambda");
              SideBandStrongWithOutLambda.Write();
              SideBandCoulombWithLambda.SetName("SideBandCoulombWithLambda");
              SideBandCoulombWithLambda.Write();
              SideBandCoulombWithOutLambda.SetName(
                  "SideBandCoulombWithOutLambda");
              SideBandCoulombWithOutLambda.Write();
              fitme->Write();
              file->Close();

              if (FAST_PLOT) {
                TPaveText* info4 = new TPaveText(0.2, 0.5, 0.9, 0.95, "blNDC");  //lbrt
                info4->SetName("info4");
                info4->SetBorderSize(1);
                info4->SetTextSize(0.04);
                info4->SetFillColor(kWhite);
                info4->SetTextFont(22);
                TString SOURCE_NAME = "Gauss";
                double Yoffset = 1.2;

                info4->AddText(
                    TString::Format(
                        "R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
                        fitter->GetParameter("pXim", DLM_Fitter1::p_sor0),
                        fitter->GetParError("pXim", DLM_Fitter1::p_sor0)));
                info4->AddText(
                    TString::Format("p_a (COU + STRONG) = %.3f #pm %.5f",
                                    p_a_strong, p_a_strong_err));
                info4->AddText(
                    TString::Format("p_a (COU) = %.3f #pm %.5f", p_a_coulomb,
                                    p_a_coulomb_err));

                if (HaveWeABaseLineSlope) {
                  info4->AddText(
                      TString::Format(
                          "p_b (COU + STRONG) = (%.2f #pm %.2f )1e-4",
                          p_b_strong * 1e4, p_b_strong_err * 1e4));
                  info4->AddText(
                      TString::Format("p_b (COU) = (%.2f #pm %.2f )1e-4",
                                      p_b_coulomb * 1e4,
                                      p_b_coulomb_err * 1e4));
                }
                info4->AddText(
                    TString::Format(
                        "Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.5f, n#sigma=%.3f",
                        Chi2_pXim, EffNumBins_pXim,
                        Chi2_pXim / double(EffNumBins_pXim), pvalXi, nSigmaXi));
                info4->AddText(
                    TString::Format(
                        "Coulomb #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.5f, n#sigma=%.3f",
                        Chi2_pXim_COULOMB, EffNumBins_pXim,
                        Chi2_pXim_COULOMB / double(EffNumBins_pXim),
                        pvalXiCoulomb, nSigmaXiCoulomb));

                TH1F* hAxis_pXim = new TH1F("hAxis_pXim", "hAxis_pXim", 600, 0,
                                            600);
                hAxis_pXim->SetStats(false);
                hAxis_pXim->SetTitle("");
                hAxis_pXim->GetXaxis()->SetLabelSize(0.065);
                hAxis_pXim->GetXaxis()->CenterTitle();
                hAxis_pXim->GetXaxis()->SetTitleOffset(1.35);
                hAxis_pXim->GetXaxis()->SetLabelOffset(0.02);
                hAxis_pXim->GetXaxis()->SetTitleSize(0.075);
                hAxis_pXim->GetYaxis()->SetLabelSize(0.065);
                hAxis_pXim->GetYaxis()->CenterTitle();
                hAxis_pXim->GetYaxis()->SetTitleOffset(Yoffset);
                hAxis_pXim->GetYaxis()->SetTitleSize(0.075);
                hAxis_pXim->GetXaxis()->SetRangeUser(0, 500);
                hAxis_pXim->GetYaxis()->SetRangeUser(0.7, 4.5);    //pPb

                TCanvas* cfast = new TCanvas(
                    TString::Format("cfast_r%.2f", GaussSourceSize),
                    TString::Format("cfast_r%.2f", GaussSourceSize), 1);
//                TPad* padFast = (TPad*)cfast->cd(0);
                cfast->SetCanvasSize(1920, 1280);
                cfast->SetMargin(0.15, 0.05, 0.2, 0.05);    //lrbt
                cfast->Update();
                OliHisto_pXim->SetTitle("p#Xi^{#minus}");
                OliHisto_pXim->SetLineWidth(2);
                OliHisto_pXim->SetLineColor(kBlack);
                FitResult_pXim.SetLineWidth(2);
                FitResult_pXim.SetLineColor(kRed);
                FitResult_pXim.SetMarkerStyle(24);
                FitResult_pXim.SetMarkerColor(kRed);
                FitResult_pXim.SetMarkerSize(1);

                FitResult_pXim_COULOMB.SetLineWidth(2);
                FitResult_pXim_COULOMB.SetLineColor(kGreen);
                FitResult_pXim_COULOMB.SetMarkerStyle(24);
                FitResult_pXim_COULOMB.SetMarkerColor(kGreen);
                FitResult_pXim_COULOMB.SetMarkerSize(1);

                hAxis_pXim->Draw("axis");
                OliHisto_pXim->Draw("same");
                FitResult_pXim.Draw("CP,same");
                FitResult_pXim_COULOMB.Draw("CP,same");
                SideBandStrongWithLambda.Draw("CP,SAME");
                SideBandStrongWithOutLambda.Draw("CP,SAME");
                SideBandCoulombWithLambda.Draw("CP,SAME");
                SideBandCoulombWithOutLambda.Draw("CP,SAME");
                info4->Draw("same");
                TF1* blPXiStrong = new TF1("blPP", "pol1", 0, 500);
                blPXiStrong->SetParameters(p_a_strong, p_b_strong);
                blPXiStrong->SetLineColor(2);
                blPXiStrong->SetLineStyle(8);

                TF1* blPXiCoulomb = new TF1("blPP", "pol1", 0, 500);
                blPXiCoulomb->SetParameters(p_a_coulomb, p_b_coulomb);
                blPXiCoulomb->SetLineColor(3);
                blPXiCoulomb->SetLineStyle(4);

                TH1F* hAxis_pXimILOVEROOT = (TH1F*) hAxis_pXim->Clone(
                    "ILOVEROOT");

                hAxis_pXim->Draw("axis");
                hAxis_pXimILOVEROOT->Draw("axig same");
                OliHisto_pXim->Draw("same");
                FitResult_pXim.Draw("CP,same");
                FitResult_pXim_COULOMB.Draw("CP,same");
                SideBandStrongWithLambda.SetLineColor(4);
                SideBandStrongWithLambda.Draw("CP,SAME");
                blPXiStrong->Draw("Same");
                blPXiCoulomb->Draw("Same");
                info4->Draw("same");
                cfast->SetGridx();
                cfast->SetGridy();
                cfast->Update();
                cfast->SaveAs(
                    TString::Format("%s/cfast_%uNumIter_%uVar.png",
                                    OutputDir.Data(), NumIter, uIter));
                delete info4;
                delete hAxis_pXim;
                delete cfast;
                delete blPXiStrong;
                delete blPXiCoulomb;
              }

              delete fitter;
              delete Ck_pXim;
              delete Ck_pXim1530;
              delete Ck_SideBand;
              delete StrongWithOutLambda;
              delete StrongWithLambda;
              delete CoulombWithOutLambda;
              delete CoulombWithLambda;
              uIter++;
              if (NumIter == 0) {
                goto outofloop;
              }
            }
          }
        }
      }
    }
  }
  outofloop: OutFile->cd();
  ntResult->Write();
  gr12->SetName("ContourBaseline");
  gr12->Write();
  funct_0->Write();
  funct_1->Write();
  OutFile->Close();

//  delete ntResult;
  return;
}

int main(int argc, char *argv[]) {
  GetXiForRadius(atoi(argv[1]), argv[2], argv[3], argv[4]);
  return 0;
}
//
//
////Global Fit default
//	//baseline
////	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.06826);
////	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 1.2624e-13);
////	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);
////	//gaussian radius
////	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);
////	//normalization
////	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -0.933815);
//
//	//Standalone Fit default
//	//baseline
////	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.05447);
////	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 2.10185e-07);
////	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);
////	//gaussian radius
////	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);
////	//normalization
////	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -0.933603);
//
//	//Fit BL & Normalization
//	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);
//
//	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 0);
////	fitter->SetParameter("pXim", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
////	fitter->SetParameter("pXim", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
//
////	fitter->SetParameter("pXim", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -1.);
//
//	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);

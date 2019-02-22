#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "ReadDreamFile.h"
#include "DreamCF.h"
#include "DreamPair.h"
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
void FitPPVariations(const unsigned& NumIter, int system, TString InputDir,
                     TString OutputDir);

int main(int argc, char *argv[]) {
  FitPPVariations(atoi(argv[1]), atoi(argv[2]), argv[3], argv[4]);
  return 0;
}

void FitPPVariations(const unsigned& NumIter, int system, TString InputDir,
                     TString OutputDir) {
//	TRandom3 rangen(1 + JobID);
  TRandom3 rangen(0);
  TString HistppName = "hCk_ReweightedppMeV_0";

  int Rebin = 5;
  //!SETTING THAT YOU MIGHT NEED
  //What source to use: 0 = Gauss; 1=Cauchy; 2=DoubleGauss
  //int vSource = rangen.Integer(2); if(vSource==1) vSource=3; //chose Gauss/EPOS at random

  const bool FAST_PLOT = true;
  const bool FULL_PLOT = false;
  TString CalibBaseDir = "";
  std::cout << "SYSTEM: " << system << std::endl;
  if (system == 0) {  // pPb MB
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
  } else if (system == 1) {  // pp MB
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  } else if (system == 2) {  // pp HM
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  } else if (system == 11) {  // pPb at the Grid
    CalibBaseDir +=
        "/home/gu74req/Analysis/CATS_Input/SystematicsAndCalib/pPbRun2_MB/";
//    system = 0;
  }

  //just some temp folder for temp files. Really not important, I will try to get rid of this soon
  //  const TString OutputDir = "/Users/bernhardhohlweger/CATSOutput/";

  //perform the systematics by using the different cut variations as further permutations, do NOT add the systematic errors in quadrature
  //TString SystematicsType = "CutVarIterate";
  //perform the systematics by adding the systematics errors of the bin quadratically to the data points
  //  TString SystematicsType = "CutVarAdd";
  //  const bool ExcludeXi = false; //(fit only pp/pL/LL just as in run1)

  //!Binning
  //This is for the CATS objects, make sure it covers the full femto range
  const int binwidth = 4 * Rebin;
  const unsigned NumMomBins = 105;
  const double kMin = 0;
  const double kMax = kMin + 4 * NumMomBins;  //(4 is the bin width)
  unsigned int momBins = 0;

  //if you modify you may need to change the CATS ranges somewhere below
  double FemtoRegion_pp[3][2];
  FemtoRegion_pp[0][0] = kMin;
  FemtoRegion_pp[0][1] = 350;
  FemtoRegion_pp[1][0] = kMin;
  FemtoRegion_pp[1][1] = 375;
  FemtoRegion_pp[2][0] = kMin;
  FemtoRegion_pp[2][1] = 400;

  //	220-320 / 260-360/ 240-340
  //!The baseline region (the same for all other systems)
  double BlRegion[2];
  BlRegion[0] = 420;
  BlRegion[1] = 420;

  double PurityProton;
  double PurityLambda;
  double PurityXi;
  double pp_f0;
  double pp_f1;
  double pL_f0;
  double pL_f1;
  double pL_f2;

  //(single particle quantities)
  //pPb
  if (system == 0) {
    PurityProton = 0.984266;  //pPb 5 TeV
    PurityLambda = 0.937761;
    PurityXi = 0.88;  //new cuts

    pp_f0 = 0.862814;
    pp_f1 = 0.09603;

    pL_f0 = 0.521433;  //fraction of primary Lambdas
    pL_f1 = 0.173811;  //fraction of Sigma0
    pL_f2 = 0.152378;  //fractions of Xi0/m
  } else {  // pp MB + HM
    PurityProton = 0.991213;
    PurityLambda = 0.965964;
    PurityXi = 0.956;

    pp_f0 = 0.874808;
    pp_f1 = 0.0876342;  //fraction of

    pL_f0 = 0.619493;  //fraction of primary Lambdas
    pL_f1 = 0.206498;  //fraction of Sigma0
    pL_f2 = 0.0870044;  //fractions of Xi0/m
  }
  //ADVANCED***
  double ProtonPrim = pp_f0;
  double arrayPercLamProton[3] = { pp_f1 / (1. - pp_f0) * 0.8, pp_f1
      / (1. - pp_f0), pp_f1 / (1. - pp_f0) * 1.2 };  //+/- 20%

  const double SigLambdaPrimDir = pL_f0 + pL_f1;
  double arrayPercSigLambda[3] = { 0.8 * pL_f1 / pL_f0, pL_f1 / pL_f0, 1.2
      * pL_f1 / pL_f0 };  //1/3 +/- 20%
  double arrayPercXiLambda[3] = { pL_f2 / (1. - pL_f0 - pL_f1) * 0.8, pL_f2
      / (1. - pL_f0 - pL_f1), pL_f2 / (1. - pL_f0 - pL_f1) * 1.2 };  //+/- 20%

  //ratio Xi-(1530) to Xi-
  const double Xim1530_to_Xim = 0.32 * (1. / 3.);
  //ratio Xi0(1530) to Xi0 (n=neutral)
  const double Xin1530_to_Xim = 0.32 * (2. / 3.);
  const double Omegam_to_Xim = 0.1;
  const double OmegamXim_BR = 0.086;

  //following my lambda pars with the 3 possible modifications
  //for the proton:
  //0 = primary
  //1 = from Lambda
  //2 = other feeddown (flat)
  //3 = missidentified
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

  //for the Lambda:
  //0 = primary
  //1 = from Sigma0
  //2 = from Xim
  //3 = from Xi0
  //4 = missidentified
  const unsigned NumChannels_L = 5;
  double*** Purities_L = new double**[3];
  double*** Fraction_L = new double**[3];
  for (unsigned uVarSL = 0; uVarSL < 3; uVarSL++) {
    Purities_L[uVarSL] = new double*[3];
    Fraction_L[uVarSL] = new double*[3];
    for (unsigned uVarXi = 0; uVarXi < 3; uVarXi++) {
      Purities_L[uVarSL][uVarXi] = new double[NumChannels_L];
      Fraction_L[uVarSL][uVarXi] = new double[NumChannels_L];

      Purities_L[uVarSL][uVarXi][0] = PurityLambda;
      Purities_L[uVarSL][uVarXi][1] = PurityLambda;
      Purities_L[uVarSL][uVarXi][2] = PurityLambda;
      Purities_L[uVarSL][uVarXi][3] = PurityLambda;
      Purities_L[uVarSL][uVarXi][4] = 1. - PurityLambda;

      //the array is r = S/L, and S+L=1 are the fractions of Sigmas and Lambdas
      double FracOfLambda = 1. / (1. + arrayPercSigLambda[uVarSL]);
      Fraction_L[uVarSL][uVarXi][0] = SigLambdaPrimDir * FracOfLambda;
      Fraction_L[uVarSL][uVarXi][1] = SigLambdaPrimDir * (1. - FracOfLambda);
      Fraction_L[uVarSL][uVarXi][2] = (1. - SigLambdaPrimDir)
          * (arrayPercXiLambda[uVarXi]);
      Fraction_L[uVarSL][uVarXi][3] = (1. - SigLambdaPrimDir)
          * (1. - arrayPercXiLambda[uVarXi]);
      Fraction_L[uVarSL][uVarXi][4] = 1.;
    }
  }

  //for the Xi:
  //0 = primary
  //1 = from Xi-(1530)
  //2 = from Xi0(1530)
  //3 = from Omega
  //4 = missidentified
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

  //starting value, do not worry about it too much
  const double GaussSourceSize = 1.2;

  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  const char *prefix;
  if (system == 0) {
    CATSinput->SetSigmaFileName("Sample3_MeV_compact.root");
    prefix = "MB";
  } else if (system == 1) {
    CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
    prefix = "MB";
  } else if (system == 2) {
    CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
    prefix = "HM";
  }
  CATSinput->ReadSigmaFile();
  CATSinput->ReadCorrelationFile(InputDir.Data(), prefix);

  //each JOB produces a separate output file
  TFile* OutFile = new TFile(
      TString::Format("%s/OutFile_%s_Iter%u.root", OutputDir.Data(),
                      "CutVarAdd", NumIter),
      "recreate");
  //you save a lot of stuff in an NTuple
  TNtuple* ntResult = new TNtuple(
      "ntResult", "ntResult",
      "IterID:vFemReg_pp:vFrac_pp_pL:vModpL:vFrac_pL_pSigma0:"
      "vFrac_pL_pXim:BLSlope:Radius_pp:RadiusErr_pp:"
      "pa_pp:paErr_pp:pb_pp:pbErr_pp:pCl_pp:pClErr_pp:"
      "Chi2NdfGlobal:pval:Chi2NdfLocal");

  Float_t ntBuffer[18];

  int vFemReg_pp;  //which femto region we use for pp (1 = default)
  int vMod_pL;  //which pL function to use: //0=exact NLO (at the moment temporary it is Usmani); 1=Ledni NLO; 2=Ledni LO; 3=ESC08
  int vFrac_pp_pL;  //fraction of protons coming from Lambda variation (1 = default)
  int vFrac_pL_pSigma0;  //fraction of Lambdas coming from Sigma0 variation (1 = default)
  int vFrac_pL_pXim;  //fraction of Lambdas coming from Xim (1 = default)
  bool HaveWeABaseLine = true;

  TidyCats* tidy = new TidyCats();
  int uIter = 1;

  for (vFemReg_pp = 0; vFemReg_pp < 3; ++vFemReg_pp) {
    for (vMod_pL = 1; vMod_pL < 3; ++vMod_pL) {
      for (vFrac_pL_pSigma0 = 0; vFrac_pL_pSigma0 < 3; ++vFrac_pL_pSigma0) {
        for (vFrac_pL_pXim = 0; vFrac_pL_pXim < 3; ++vFrac_pL_pXim) {
          for (int BaselineSlope = 0; BaselineSlope < 2; ++BaselineSlope) {
            if (BaselineSlope == 0) {
              HaveWeABaseLine = true;  //use baseline
            } else {
              HaveWeABaseLine = false;  // no baseline
            }
            if (NumIter == 0) {
              vFemReg_pp = 1;
              vFrac_pp_pL = 1;
              vFrac_pL_pSigma0 = 1;
              vFrac_pL_pXim = 1;
              vMod_pL = 1;
              HaveWeABaseLine = false;  // no base line in default
            }
//            double Pars_pp[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                GaussSourceSize / 1.2, 0.5 };
            CATS AB_pp;
            tidy->GetCatsProtonProton(&AB_pp, NumMomBins, kMin, kMax, false);
            AB_pp.KillTheCat();

//            double Pars_pL[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                GaussSourceSize / 1.2, 0.5 };
            CATS AB_pL;
            tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax);
            AB_pL.KillTheCat();
//            double Pars_pXi[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                GaussSourceSize / 1.2, 0.5 };
            CATS AB_pXim;
            tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins, kMin, kMax, true,
                                       13);
            AB_pXim.KillTheCat();

//            double Pars_pXim1530[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                GaussSourceSize / 1.2, 0.5 };
            CATS AB_pXim1530;
            tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins, kMin,
                                           kMax);
            AB_pXim1530.KillTheCat();

            std::cout << "Reading Data \n";

            CATSinput->ObtainCFs(Rebin, 240, 340);

            TH1F* OliHisto_pp = CATSinput->GetCF("pp", HistppName);
            if (!OliHisto_pp)
              std::cout << HistppName.Data() << " Missing" << std::endl;
            CATSinput->AddSystematics("C2totalsysPP.root", OliHisto_pp);
            std::cout << "pp Done\n";

            //!CHANGE PATH HERE

            const unsigned NumSourcePars = 1;

            //this way you define a correlation function using a CATS object.
            //needed inputs: num source/pot pars, CATS obj
            DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp);
            //    DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL);
            DLM_Ck* Ck_pL = new DLM_Ck(1, 4, NumMomBins, kMin, kMax,
                                       Lednicky_SingletTriplet);
            //this way you define a correlation function using Lednicky.
            //needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
            DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins, kMin, kMax,
                                            Lednicky_gauss_Sigma0);
            Ck_pSigma0->SetSourcePar(0, GaussSourceSize);
            DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
            DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);

            if (vMod_pL == 1) {
              Ck_pL->SetPotPar(0, 2.91);
              Ck_pL->SetPotPar(1, 2.78);
              Ck_pL->SetPotPar(2, 1.54);
              Ck_pL->SetPotPar(3, 2.72);
            } else if (vMod_pL == 2) {
              Ck_pL->SetPotPar(0, 1.91);
              Ck_pL->SetPotPar(1, 1.4);
              Ck_pL->SetPotPar(2, 1.23);
              Ck_pL->SetPotPar(3, 2.13);
            }

            Ck_pp->Update();
            Ck_pL->Update();
            Ck_pSigma0->Update();
            Ck_pXim->Update();
            Ck_pXim1530->Update();
            if (!CATSinput->GetSigmaFile(0)) {
              std::cout << "No Sigma file 0 \n";
              return;
            }
            if (!CATSinput->GetSigmaFile(1)) {
              std::cout << "No Sigma file 1 \n";
              return;
            }
            if (!CATSinput->GetSigmaFile(2)) {
              std::cout << "No Sigma file 2 \n";
              return;
            }
            if (!CATSinput->GetSigmaFile(3)) {
              std::cout << "No Sigma file 3 \n";
              return;
            }
            DLM_CkDecomposition CkDec_pp("pp", 3, *Ck_pp,
                                         CATSinput->GetSigmaFile(0));
            DLM_CkDecomposition CkDec_pL("pLambda", 4, *Ck_pL,
                                         CATSinput->GetSigmaFile(1));
            DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0,
            NULL);
            DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
                                           CATSinput->GetSigmaFile(3));
            DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530,
            NULL);

            double lam_pp = Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0];
            double lam_pp_pL = Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][1]
                * Fraction_p[vFrac_pp_pL][1] * 2;
            double lam_pp_fake = Purities_p[vFrac_pp_pL][3]
                * Purities_p[vFrac_pp_pL][0]
                + Purities_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][3]
                + Purities_p[vFrac_pp_pL][3] * Purities_p[vFrac_pp_pL][3];

            printf("lam_pp = %.3f\n", lam_pp);
            printf("lam_pp_pL = %.3f\n", lam_pp_pL);
            printf("lam_pp_fake = %.3f\n", lam_pp_fake);
            printf("\n");
            if (!CATSinput->GetResFile(0)) {
              std::cout << "No Calib 0 \n";
              return;
            }
            CkDec_pp.AddContribution(0, lam_pp_pL,
                                     DLM_CkDecomposition::cFeedDown, &CkDec_pL,
                                     CATSinput->GetResFile(0));
            CkDec_pp.AddContribution(1, 1. - lam_pp - lam_pp_pL - lam_pp_fake,
                                     DLM_CkDecomposition::cFeedDown);
            CkDec_pp.AddContribution(2, lam_pp_fake,
                                     DLM_CkDecomposition::cFake);  //0.02

            double lam_pL = Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0]
                * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
                * Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
            double lam_pL_pS0 = Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0]
                * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1]
                * Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1];
            double lam_pL_pXm = Purities_p[vFrac_pp_pL][0]
                * Fraction_p[vFrac_pp_pL][0]
                * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2]
                * Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2];
            double lam_pL_fake = Purities_p[vFrac_pp_pL][3]
                * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
                + Purities_p[vFrac_pp_pL][0]
                    * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]
                + Purities_p[vFrac_pp_pL][3]
                    * Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];

            printf("lam_pL=%.3f\n", lam_pL);
            printf("lam_pL_pS0=%.3f\n", lam_pL_pS0);
            printf("lam_pL_pXm=%.3f\n", lam_pL_pXm);
            printf("lam_pL_fake=%.3f\n", lam_pL_fake);
            printf("\n");
            if (!CATSinput->GetResFile(1)) {
              std::cout << "No Calib 1 \n";
              return;
            }
            if (!CATSinput->GetResFile(2)) {
              std::cout << "No Calib 2 \n";
              return;
            }

            CkDec_pL.AddContribution(0, lam_pL_pS0,
                                     DLM_CkDecomposition::cFeedDown,
                                     &CkDec_pSigma0, CATSinput->GetResFile(1));
            CkDec_pL.AddContribution(1, lam_pL_pXm,
                                     DLM_CkDecomposition::cFeedDown,
                                     &CkDec_pXim, CATSinput->GetResFile(2));
            CkDec_pL.AddContribution(
                2, 1. - lam_pL - lam_pL_pS0 - lam_pL_pXm - lam_pL_fake,
                DLM_CkDecomposition::cFeedDown);
            CkDec_pL.AddContribution(3, lam_pL_fake,
                                     DLM_CkDecomposition::cFake);  //0.03

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
            if (!CATSinput->GetResFile(3)) {
              std::cout << "No Calib 3 \n";
              return;
            }
            CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
                                       DLM_CkDecomposition::cFeedDown,
                                       &CkDec_pXim1530,
                                       CATSinput->GetResFile(3));  //from Xi-(1530)
            CkDec_pXim.AddContribution(
                1, 1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
                DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
            CkDec_pXim.AddContribution(2, lam_pXim_fake,
                                       DLM_CkDecomposition::cFake);

            DLM_Fitter1* fitter = new DLM_Fitter1(1);
            fitter->SetOutputDir(OutputDir.Data());

            std::cout << "pp" << std::endl;
            fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp,
                              FemtoRegion_pp[vFemReg_pp][0],
                              FemtoRegion_pp[vFemReg_pp][1], BlRegion[0],
                              BlRegion[1]);
            fitter->SetSeparateBL(0, false);
            fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
            if (HaveWeABaseLine) {
              fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, -2e-3, 2e-3);
              std::cout << "Fitting ranges for BL set \n";
            } else {
              fitter->FixParameter("pp", DLM_Fitter1::p_b, 0);
            }

            fitter->AddSameSource("pLambda", "pp", 1);
            fitter->AddSameSource("pSigma0", "pp", 1);
            fitter->AddSameSource("pXim", "pp", 1);
            fitter->AddSameSource("pXim1530", "pp", 1);

            fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);

            fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.2, 0.8, 1.8);
            fitter->FixParameter("pp", DLM_Fitter1::p_Cl, -1);
            std::cout << "CL Fixed \n";

            //    fitter->FixParameter("pp", DLM_Fitter1::p_sor0, 1.383);
            CkDec_pp.Update();
            CkDec_pL.Update();
            CkDec_pXim.Update();
            fitter->GoBabyGo();

            TFile* GraphFile = new TFile(
                TString::Format("%s/GraphFile_%s_Iter%u_uIter%u.root",
                                OutputDir.Data(), "CutVarAdd", NumIter, uIter),
                "recreate");

            GraphFile->cd();

            TGraph FitResult_pp;
            FitResult_pp.SetName(TString::Format("FitResult_pp_%u", uIter));
            fitter->GetFitGraph(0, FitResult_pp);

            double Chi2 = fitter->GetChi2();
            unsigned NDF = fitter->GetNdf();
            double RadiusResult = fitter->GetParameter("pp",
                                                       DLM_Fitter1::p_sor0);
            double RadiusError = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
            if (Chi2 / double(NDF) != fitter->GetChi2Ndf()) {
              printf("Oh boy...\n");
            }

            double Chi2_pp = 0;
            unsigned EffNumBins_pp = 0;

            for (unsigned uBin = 0; uBin < 50; uBin++) {

              double mom = AB_pp.GetMomentum(uBin);
              double dataY;
              double dataErr;
              double theoryX;
              double theoryY;

              if (mom > FemtoRegion_pp[vFemReg_pp][1])
                continue;

              FitResult_pp.GetPoint(uBin, theoryX, theoryY);
              if (mom != theoryX) {
                std::cout << mom << '\t' << theoryX << std::endl;
                printf("  PROBLEM pp!\n");
              }
              dataY = OliHisto_pp->GetBinContent(uBin + 1);
              dataErr = OliHisto_pp->GetBinError(uBin + 1);
              if (dataErr < 1e-5) {
                std::cout << dataErr << '\t'
                          << "WARNING POINT NOT CONSIDERED \n";
                continue;
              }
              Chi2_pp += (dataY - theoryY) * (dataY - theoryY)
                  / (dataErr * dataErr);
              EffNumBins_pp++;
            }
            ntBuffer[0] = uIter;
            ntBuffer[1] = vFemReg_pp;
            ntBuffer[2] = vFrac_pp_pL;
            ntBuffer[3] = vMod_pL;
            ntBuffer[4] = vFrac_pL_pSigma0;
            ntBuffer[5] = vFrac_pL_pXim;
            ntBuffer[6] = (int) HaveWeABaseLine;
            ntBuffer[7] = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
            ntBuffer[8] = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
            ntBuffer[9] = fitter->GetParameter("pp", DLM_Fitter1::p_a);
            ntBuffer[10] = fitter->GetParError("pp", DLM_Fitter1::p_a);
            ntBuffer[11] = fitter->GetParameter("pp", DLM_Fitter1::p_b);
            ntBuffer[12] = fitter->GetParError("pp", DLM_Fitter1::p_b);
            ntBuffer[13] = fitter->GetParameter("pp", DLM_Fitter1::p_Cl);
            ntBuffer[14] = fitter->GetParError("pp", DLM_Fitter1::p_Cl);
            ntBuffer[15] = fitter->GetChi2Ndf();
            ntBuffer[16] = fitter->GetPval();
            ntBuffer[17] = Chi2_pp / double(EffNumBins_pp);

            //    ntBuffer[23] = (int) FIX_CL;
            ntResult->Fill(ntBuffer);
            printf("Chi2_pp/bins = %.2f/%u = %.2f\n", Chi2_pp, EffNumBins_pp,
                   Chi2_pp / double(EffNumBins_pp));

            std::cout << "float radius = "
                      << fitter->GetParameter("pp", DLM_Fitter1::p_sor0) << ";"
                      << std::endl;
            std::cout << "float ppBL0 ="
                      << fitter->GetParameter("pp", DLM_Fitter1::p_a) << ";"
                      << std::endl;
            std::cout << "err ppBL0 = "
                      << fitter->GetParError("pp", DLM_Fitter1::p_a) << ";"
                      << std::endl;

            std::cout << "float ppBL1 ="
                      << fitter->GetParameter("pp", DLM_Fitter1::p_b) << ";"
                      << std::endl;
            std::cout << "err ppBL1 = "
                      << fitter->GetParError("pp", DLM_Fitter1::p_b) << ";"
                      << std::endl;

            std::cout << "float ppNorm ="
                      << fitter->GetParameter("pp", DLM_Fitter1::p_Cl) << ";"
                      << std::endl;
            std::cout << "err ppNorm = "
                      << fitter->GetParError("pp", DLM_Fitter1::p_Cl) << ";"
                      << std::endl;

            GraphFile->cd();
            FitResult_pp.Write();
            fitter->GetFit()->SetName(TString::Format("GlobalFit_%u", uIter));
            fitter->GetFit()->Write();

            if (FAST_PLOT) {
              TPaveText* info1 = new TPaveText(0.45, 0.65, 0.9, 0.95, "blNDC");  //lbrt
              info1->SetName("info1");
              info1->SetBorderSize(1);
              info1->SetTextSize(0.04);
              info1->SetFillColor(kWhite);
              info1->SetTextFont(22);
              TString SOURCE_NAME;
              SOURCE_NAME = "Gauss";
              info1->AddText(
                  TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
                                  RadiusResult, RadiusError));
              info1->AddText(
                  TString::Format("p_a = %.3f #pm %.5f",
                                  fitter->GetParameter("pp", DLM_Fitter1::p_a),
                                  fitter->GetParError("pp", DLM_Fitter1::p_a)));

              if (HaveWeABaseLine) {
                info1->AddText(
                    TString::Format(
                        "p_b = (%.2f #pm %.2f )1e-4",
                        fitter->GetParameter("pp", DLM_Fitter1::p_b) * 1e4,
                        fitter->GetParError("pp", DLM_Fitter1::p_b) * 1e4));
              }
              info1->AddText(
                  TString::Format(
                      "Global #chi^{2}_{ndf}=%.0f/%u=%.2f, p_{val}=%.3f", Chi2,
                      NDF, Chi2 / double(NDF), TMath::Prob(Chi2, round(NDF))));
              info1->AddText(
                  TString::Format(
                      "Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f", Chi2_pp,
                      EffNumBins_pp, Chi2_pp / double(EffNumBins_pp),
                      TMath::Prob(Chi2_pp, round(EffNumBins_pp))));

              double Yoffset = 1.2;
              TH1F* hAxis_pp = new TH1F("hAxis_pp", "hAxis_pp", 600, 0, 600);
              hAxis_pp->SetStats(false);
              hAxis_pp->SetTitle("");
              hAxis_pp->GetXaxis()->SetLabelSize(0.065);
              hAxis_pp->GetXaxis()->CenterTitle();
              hAxis_pp->GetXaxis()->SetTitleOffset(1.35);
              hAxis_pp->GetXaxis()->SetLabelOffset(0.02);
              hAxis_pp->GetXaxis()->SetTitleSize(0.075);
              hAxis_pp->GetYaxis()->SetLabelSize(0.065);
              hAxis_pp->GetYaxis()->CenterTitle();
              hAxis_pp->GetYaxis()->SetTitleOffset(Yoffset);
              hAxis_pp->GetYaxis()->SetTitleSize(0.075);
              hAxis_pp->GetXaxis()->SetRangeUser(0, 600);
              hAxis_pp->GetYaxis()->SetRangeUser(0.5, 4.5);
              TF1* blPP = new TF1("blPP", "pol1", 0, 500);
              blPP->SetParameters(fitter->GetParameter("pp", DLM_Fitter1::p_a),
                                  fitter->GetParameter("pp", DLM_Fitter1::p_b));
              blPP->SetLineColor(7);
              blPP->SetLineWidth(4);

              TString CanName = Form("cfast");
              //          TString CanName = Form("cFastBL_%u_%u.root",
              //              (int) BlRegion_pp[1][0],
              //              (int) BlRegion_pp[1][1]);
              TCanvas* cfast = new TCanvas(CanName.Data(), CanName.Data(), 1);
              cfast->cd(0);
              cfast->SetCanvasSize(1920, 1280);
              cfast->SetMargin(0.15, 0.05, 0.2, 0.05);    //lrbt

              cfast->cd(1);
              OliHisto_pp->SetStats(false);
              OliHisto_pp->SetTitle("pp");
              OliHisto_pp->SetLineWidth(2);
              OliHisto_pp->SetLineColor(kBlack);
              FitResult_pp.SetLineWidth(2);
              FitResult_pp.SetLineColor(kRed);
              FitResult_pp.SetMarkerStyle(24);
              FitResult_pp.SetMarkerColor(kRed);
              FitResult_pp.SetMarkerSize(1);

              hAxis_pp->Draw("axis");
              OliHisto_pp->Draw("same");
              blPP->Draw("same");
              FitResult_pp.Draw("CP,same");
              info1->Draw("same");

              cfast->Write();
              cfast->SaveAs(
                  TString::Format("%s/Iter%u_uIter%u.png", OutputDir.Data(),
                                  NumIter, uIter));
              delete hAxis_pp;
              delete info1;
              delete cfast;

            }      //FAST_PLOT
            GraphFile->Close();
            delete OliHisto_pp;
            delete fitter;
            delete Ck_pp;
            delete Ck_pL;
            delete Ck_pSigma0;
            delete Ck_pXim;
            delete Ck_pXim1530;

            delete GraphFile;
            uIter++;

            if (NumIter == 0) {
              goto outofloop;
            }
          }
        }
      }
    }
  }
  outofloop:
  OutFile->cd();
  ntResult->Write();

  delete ntResult;
  delete OutFile;

  for (unsigned uVar = 0; uVar < 3; uVar++) {
    delete[] Purities_p[uVar];
    delete[] Fraction_p[uVar];
    delete[] Purities_L[uVar];
    delete[] Fraction_L[uVar];
  }
  delete[] Purities_p;
  delete[] Fraction_p;
  delete[] Purities_L;
  delete[] Fraction_L;

}

//if (uIter == 0) {
//			vBlReg = 0;
//			HaveWeABaseLine = false; // no base line in default
//		} else if (uIter == 1) {
//			vBlReg = 1;
//			HaveWeABaseLine = false; // no base line in default
//		} else if (uIter == 2) {
//			vBlReg = 2;
//			HaveWeABaseLine = false; // no base line in default
//		} else if (uIter == 3) {
//			vBlReg = 0;
//			HaveWeABaseLine = true; // no base line in default
//		} else if (uIter == 4) {
//			vBlReg = 1;
//			HaveWeABaseLine = true; // no base line in default
//		} else if (uIter == 5) {
//			vBlReg = 2;
//			HaveWeABaseLine = true; // no base line in default
//		}

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
void ppmTBins(TString InputDir, TString OutputDir, int system, int numkTBins);

int main(int argc, char *argv[]) {
  ppmTBins(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
  return 0;
}

void ppmTBins(TString InputDir, TString OutputDir, int system, int numkTBins) {
//  TRandom3 rangen(1 + JobID);
  TString HistppName = "hCkTotNormWeightMeV_mTBin_";

  //!SETTING THAT YOU MIGHT NEED
  //What source to use: 0 = Gauss; 1=Cauchy; 2=DoubleGauss
  TFile* inFile = TFile::Open(Form("%s/CFOutputALL_mT_pp_HM.root",InputDir.Data()));
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

  //perform the systematics by using the different cut variations as further permutations, do NOT add the systematic errors in quadrature
  //TString SystematicsType = "CutVarIterate";
  //perform the systematics by adding the systematics errors of the bin quadratically to the data points
  //  TString SystematicsType = "CutVarAdd";
  //  const bool ExcludeXi = false; //(fit only pp/pL/LL just as in run1)

  //!Binning
  //This is for the CATS objects, make sure it covers the full femto range
  const int binwidth = 4;
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

  //  220-320 / 260-360/ 240-340
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

  //starting value, do not worry about it too much
  const double GaussSourceSize = 0.8;
  const double GaussSourceSize_pp = 1.2;

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
//  CATSinput->ReadCorrelationFile(InputDir.Data(), prefix);

  //each JOB produces a separate output file
  TFile* OutFile = new TFile(
      TString::Format("%s/OutFile.root", OutputDir.Data()), "recreate");

  TFile* GraphFile = new TFile(
      TString::Format("%s/GraphFile.root", OutputDir.Data()), "recreate");

  //you save a lot of stuff in an NTuple
  TNtuple* ntResult = new TNtuple(
      "ntResult", "ntResult",
      "IterID:vFemReg_pp:vFrac_pp_pL:vModpL:vFrac_pL_pSigma0:"
      "vFrac_pL_pXim:BLSlope:Radius_pp:RadiusErr_pp:"
      "pa_pp:paErr_pp:pb_pp:pbErr_pp:pCl_pp:pClErr_pp:"
      "Chi2NdfGlobal:pval:Chi2NdfLocal");
  Float_t ntBuffer[18];
  for (int ikT = 0; ikT < numkTBins; ++ikT) {

    TidyCats* tidy = new TidyCats();
//    double Pars_pp[6] = { 0, 0, 0, GaussSourceSize * 1.2,
//                  GaussSourceSize / 1.2, 0.5 };
    double massProton = 938.272;
    double massPion =   139.570;
    double Pars_pp[9] =
        { 0, 0, 0, GaussSourceSize_pp * 1.2, 1.65, 0.3578, 1361.52, massProton, massPion };
    CATS AB_pp;
    tidy->GetCatsProtonProton(&AB_pp, Pars_pp, NumMomBins,
                              kMin, kMax, true);
    AB_pp.KillTheCat();

    double Pars_pL[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize / 1.2,
        0.5 };
    CATS AB_pL;
    tidy->GetCatsProtonLambda(&AB_pL, Pars_pL, NumMomBins,
                              kMin, kMax);
    AB_pL.KillTheCat();

    std::cout << "Reading Data \n";

//    CATSinput->ObtainCFs(5, 240, 340);
    TString HistNamemT=HistppName;
    HistNamemT+=ikT;
    TH1F* OliHisto_pp = (TH1F*)inFile->Get(HistNamemT.Data());
    if (!OliHisto_pp)
      std::cout << HistNamemT.Data() << " Missing" << std::endl;
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

    Ck_pp->Update();
    Ck_pL->Update();
    if (!CATSinput->GetSigmaFile(0)) {
      std::cout << "No Sigma file 0 \n";
      return;
    }
    if (!CATSinput->GetSigmaFile(1)) {
      std::cout << "No Sigma file 1 \n";
      return;
    }

    DLM_CkDecomposition CkDec_pp("pp", 3, *Ck_pp, CATSinput->GetSigmaFile(0));
    DLM_CkDecomposition CkDec_pL("pLambda", 4, *Ck_pL,
                                 CATSinput->GetSigmaFile(1));

    double lam_pp = Purities_p[1][0] * Fraction_p[1][0] * Purities_p[1][0]
        * Fraction_p[1][0];
    double lam_pp_pL = Purities_p[1][0] * Fraction_p[1][0] * Purities_p[1][1]
        * Fraction_p[1][1] * 2;
    double lam_pp_fake = Purities_p[1][3] * Purities_p[1][0]
        + Purities_p[1][0] * Purities_p[1][3]
        + Purities_p[1][3] * Purities_p[1][3];

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
//    CkDec_pp.AddContribution(1, lam_pp_fake+1-lam_pp, DLM_CkDecomposition::cFake);  //0.02

    double lam_pL = Purities_p[1][0] * Fraction_p[1][0] * Purities_L[1][1][0]
        * Fraction_L[1][1][0];
    double lam_pL_pS0 = Purities_p[1][0] * Fraction_p[1][0]
        * Purities_L[1][1][1] * Fraction_L[1][1][1];
    double lam_pL_pXm = Purities_p[1][0] * Fraction_p[1][0]
        * Purities_L[1][1][2] * Fraction_L[1][1][2];
    double lam_pL_fake = Purities_p[1][3] * Purities_L[1][1][0]
        + Purities_p[1][0] * Purities_L[1][1][4]
        + Purities_p[1][3] * Purities_L[1][1][4];

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

    CkDec_pL.AddContribution(0, lam_pL_pS0, DLM_CkDecomposition::cFeedDown,
                             &CkDec_pSigma0, CATSinput->GetResFile(1));
    CkDec_pL.AddContribution(1, lam_pL_pXm, DLM_CkDecomposition::cFeedDown,
                             &CkDec_pXim, CATSinput->GetResFile(2));
    CkDec_pL.AddContribution(
        2, 1. - lam_pL - lam_pL_pS0 - lam_pL_pXm - lam_pL_fake,
        DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3, lam_pL_fake, DLM_CkDecomposition::cFake);  //0.03

    DLM_Fitter1* fitter = new DLM_Fitter1(1);
    fitter->SetOutputDir(OutputDir.Data());

    std::cout << "pp" << std::endl;
    fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp, FemtoRegion_pp[1][0],
                      FemtoRegion_pp[1][1], BlRegion[0], BlRegion[1]);
    fitter->SetSeparateBL(0, false);
    fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
    fitter->FixParameter("pp", DLM_Fitter1::p_b, 0);

//    fitter->FixParameter("pLambda", DLM_Fitter1::p_sor0, 1.2);
//    fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0, 1.2);
//    fitter->FixParameter("pXim",    DLM_Fitter1::p_sor0, 1.2);
//    fitter->FixParameter("pXim1530",DLM_Fitter1::p_sor0, 1.2);

    fitter->AddSameSource("pLambda", "pp", 1);

    fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);

    fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.2, 0.7, 1.8);
    fitter->FixParameter("pp", DLM_Fitter1::p_Cl, -1);
    std::cout << "CL Fixed \n";

    //    fitter->FixParameter("pp", DLM_Fitter1::p_sor0, 1.383);
    CkDec_pp.Update();
    CkDec_pL.Update();
    fitter->GoBabyGo();


    TGraph FitResult_pp;
    FitResult_pp.SetName(TString::Format("FitResult_pp_%u",ikT));
    fitter->GetFitGraph(0, FitResult_pp);

    double Chi2 = fitter->GetChi2();
    unsigned NDF = fitter->GetNdf();
    double RadiusResult = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
    double RadiusError = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
    if (Chi2 / double(NDF) != fitter->GetChi2Ndf()) {
      printf("Oh boy...\n");
    }

    double Chi2_pp = 0;
    unsigned EffNumBins_pp = 0;
    for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {

      double mom = AB_pp.GetMomentum(uBin);
      double dataY;
      double dataErr;
      double theoryX;
      double theoryY;

      if (mom > FemtoRegion_pp[1][1])
        continue;

      FitResult_pp.GetPoint(uBin, theoryX, theoryY);
      if (mom != theoryX) {
        std::cout << mom << '\t' << theoryX << std::endl;
        printf("  PROBLEM pp!\n");
      }
      dataY = OliHisto_pp->GetBinContent(uBin + 1);
      dataErr = OliHisto_pp->GetBinError(uBin + 1);

      Chi2_pp += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
      EffNumBins_pp++;
    }
    ntBuffer[0] = ikT;
    ntBuffer[1] = 1;
    ntBuffer[2] = 1;
    ntBuffer[3] = 1;
    ntBuffer[4] = 1;
    ntBuffer[5] = 1;
    ntBuffer[6] = (int) false;
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
    std::cout << "float ppBL0 =" << fitter->GetParameter("pp", DLM_Fitter1::p_a)
              << ";" << std::endl;
    std::cout << "err ppBL0 = " << fitter->GetParError("pp", DLM_Fitter1::p_a)
              << ";" << std::endl;

    std::cout << "float ppBL1 =" << fitter->GetParameter("pp", DLM_Fitter1::p_b)
              << ";" << std::endl;
    std::cout << "err ppBL1 = " << fitter->GetParError("pp", DLM_Fitter1::p_b)
              << ";" << std::endl;

    std::cout << "float ppNorm ="
              << fitter->GetParameter("pp", DLM_Fitter1::p_Cl) << ";"
              << std::endl;
    std::cout << "err ppNorm = " << fitter->GetParError("pp", DLM_Fitter1::p_Cl)
              << ";" << std::endl;

    GraphFile->cd();
    FitResult_pp.Write();

    TH1F* copyOliHist = (TH1F*)OliHisto_pp->Clone(Form("CFkT_%u",ikT));
    copyOliHist->Write(copyOliHist->GetName());

    fitter->GetFit()->SetName(TString::Format("GlobalFit"));
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
          TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(), RadiusResult,
                          RadiusError));
      info1->AddText(
          TString::Format("p_a = %.3f #pm %.5f",
                          fitter->GetParameter("pp", DLM_Fitter1::p_a),
                          fitter->GetParError("pp", DLM_Fitter1::p_a)));

      info1->AddText(
          TString::Format("Global #chi^{2}_{ndf}=%.0f/%u=%.2f, p_{val}=%.3f",
                          Chi2, NDF, Chi2 / double(NDF),
                          TMath::Prob(Chi2, round(NDF))));
      info1->AddText(
          TString::Format("Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f",
                          Chi2_pp, EffNumBins_pp,
                          Chi2_pp / double(EffNumBins_pp),
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
      OutFile->cd();
      cfast->Write();
      cfast->SaveAs(TString::Format("%s/FitppmT_%u.png", OutputDir.Data(),ikT));
      delete hAxis_pp;
      delete info1;
      delete cfast;

    }      //FAST_PLOT
    delete OliHisto_pp;
    delete copyOliHist;
    delete fitter;
    delete Ck_pp;
    delete Ck_pL;

  }
  OutFile->cd();
  ntResult->Write();
  OutFile->Close();
  GraphFile->Close();
//  delete ntResult;
//  delete OutFile;
//  delete GraphFile;
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

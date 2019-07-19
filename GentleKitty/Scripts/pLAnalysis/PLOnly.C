#include "CATS.h"
#include "CATStools.h"
#include "DLM_WfModel.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"

#include "TidyCats.h"
#include "CATSLambdaParam.h"
#include "CATSInput.h"
#include "ReadDreamFile.h"
#include "DreamCF.h"
#include "DreamPair.h"

#include "TROOT.h"
#include "TNtuple.h"
#include "TCanvas.h"

#include <iostream>
#include "stdlib.h"
#include <chrono>
#include <ctime>

void FitPPVariations(const unsigned& NumIter, int imTBin, int system,
                     int source, int iPotential, TString InputFile,
                     TString HistoName, TString OutputDir) {
  auto start = std::chrono::system_clock::now();

  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  //What source to use: 0 = Gauss; 1=Resonance; 2=Levy
//  What potential to use: 0 = Scattering para; 1 = Usmani; 2 = NLO; 3 = LO;
  TString HistppName = HistoName.Data();
  TFile* inFile = TFile::Open(TString::Format("%s", InputFile.Data()), "READ");
  if (!inFile) {
    std::cout << "No input file found exiting \n";
    return;
  }
  TH1F* StoreHist = (TH1F*) inFile->Get(HistppName.Data());
  if (!StoreHist) {
    std::cout << "Histogram " << HistppName.Data() << " not found, exiting \n";
    inFile->ls();
    return;
  }
  //This is for the CATS objects, make sure it covers the full femto range

  const unsigned NumMomBins = 20;
  const double kMin = StoreHist->GetXaxis()->GetXmin();
  const double kMax = kMin + StoreHist->GetXaxis()->GetBinWidth(1) * NumMomBins;  //(4 is the bin width)

  //if you modify you may need to change the CATS ranges somewhere below
  double FemtoRegion[3];
  FemtoRegion[0] = 180;
  FemtoRegion[1] = 204;
  FemtoRegion[2] = 228;

  double BaseLineRegion[3][2];
  BaseLineRegion[0][0] = 336;
  BaseLineRegion[0][1] = 552;
  BaseLineRegion[1][0] = 360;
  BaseLineRegion[1][1] = 576;
  BaseLineRegion[2][0] = 284;
  BaseLineRegion[2][1] = 600;

  TidyCats::Sources TheSource;
  TidyCats::Sources FeeddownSource;
  if (source == 0) {
    TheSource = TidyCats::sGaussian;
    FeeddownSource = TheSource;
  } else if (source == 1) {
    TheSource = TidyCats::sResonance;
    FeeddownSource = TidyCats::sGaussian;
  } else if (source == 2) {
    TheSource = TidyCats::sLevy;
    FeeddownSource = TidyCats::sGaussian;
  } else {
    std::cout << "Source does not exist! Exiting \n";
    return;
  }

  TString CalibBaseDir = "";
  TString SigmaFileName = "";
  double PurityProton, PrimProton, SecLamProton;
  double PurityLambda, PrimLambdaAndSigma, SecLambda;
  double PurityXi;

  std::cout << "SYSTEM: " << system << std::endl;
  if (system == 0) {
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
    SigmaFileName += "Sample3_MeV_compact.root";
    PurityProton = 0.984266;  //pPb 5 TeV
    PrimProton = 0.862814;
    SecLamProton = 0.09603;

    PurityLambda = 0.937761;
    PrimLambdaAndSigma = 0.79;  //fraction of primary Lambdas + Sigma 0
    SecLambda = 0.30;  //fraction of weak decay Lambdas

    PurityXi = 0.88;  //new cuts
  } else if (system == 1) {  // pp HM
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
    SigmaFileName += "Sample6_MeV_compact.root";
    PurityProton = 0.991213;
    PrimProton = 0.874808;
    SecLamProton = 0.0876342;

    PurityLambda = 0.965964;
    PrimLambdaAndSigma = 0.806;  //fraction of primary Lambdas + Sigma 0
    SecLambda = 0.194;  //fraction of weak decay Lambdas

    PurityXi = 0.915;
  } else if (system == 2) {  // pp HM
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
    SigmaFileName += "Sample6_MeV_compact.root";
    PurityProton = 0.9943;
    PrimProton = 0.873;
    SecLamProton = 0.089;  //Fraction of Lambdas

    PurityLambda = 0.961;
    PrimLambdaAndSigma = 0.785;  //fraction of primary Lambdas + Sigma 0
    SecLambda = 0.215;  //fraction of weak decay Lambdas

    PurityXi = 0.915;
  } else {
    std::cout << "System " << system << " not implmented, extiting \n";
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
  Particle Lambda[3][3];  // 1) variation of Lambda/Sigma Ratio, 2) variation of Xi0/Xim Ratio
  Particle Xi[3][3];  //1) variation of dN/dy Omega 2) variation of dN/dy Xi1530

  int iVar1 = 0;
  for (auto it : Variation) {
    double SecFracSigma = 1. - PrimProton - it * SecLamProton;
    Proton[iVar1] = Particle(PurityProton, PrimProton, { it * SecLamProton,
                                 SecFracSigma });
    std::cout << "it: " << it << " PurityProton: " << PurityProton
              << " it * SecLamProton: " << it * SecLamProton << " SecFracSigma:"
              << SecFracSigma << std::endl;

    double LamSigProdFraction = 3 * it / 4. < 1 ? 3 * it / 4. : 1;
    double PrimLambda = LamSigProdFraction * PrimLambdaAndSigma;
    double SecSigLambda = (1. - LamSigProdFraction) * PrimLambdaAndSigma;  // decay probability = 100%!
    int iVar2 = 0;
    for (auto itXim : Variation) {
      double SecXimLambda = itXim * SecLambda / 2.;
      double SecXi0Lambda0 = (1 - itXim / 2.) * SecLambda;
      Lambda[iVar1][iVar2] = Particle(PurityLambda, PrimLambda, { SecSigLambda,
                                          SecXimLambda, SecXi0Lambda0 });

      double XiNormalization = 1 + it * OmegamXimProdFraction * OmegamXim_BR
          + itXim
              * (Xi01530XimProdFraction * Xi01530Xim_BR
                  + Xim1530XimProdFraction * Xim1530Xim_BR);
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
  //Sigma/Lam variations
  CATSLambdaParam pLLam0(Proton[0], Lambda[0][1]);
  CATSLambdaParam pLLam1(Proton[1], Lambda[1][1]);
  CATSLambdaParam pLLam2(Proton[2], Lambda[2][1]);

  const std::vector<double> lam_pL =
      { pLLam0.GetLambdaParam(CATSLambdaParam::Primary,
                              CATSLambdaParam::Primary), pLLam1.GetLambdaParam(
          CATSLambdaParam::Primary, CATSLambdaParam::Primary), pLLam2
          .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary) };

  const std::vector<double> lam_pL_pS0 =
      { pLLam0.GetLambdaParam(CATSLambdaParam::Primary,
                              CATSLambdaParam::FeedDown, 0, 0), pLLam1
          .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown,
                          0, 0), pLLam2.GetLambdaParam(
          CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0) };

  const std::vector<double> lam_pL_pXm =
      { pLLam0.GetLambdaParam(CATSLambdaParam::Primary,
                              CATSLambdaParam::FeedDown, 0, 1), pLLam1
          .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown,
                          0, 1), pLLam2.GetLambdaParam(
          CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1) };

  const std::vector<double> lam_pL_fake =
      { pLLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
          + pLLam0.GetLambdaParam(CATSLambdaParam::Fake,
                                  CATSLambdaParam::Primary)
          + pLLam0.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake),
          pLLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
              + pLLam1.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
              + pLLam1.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake), pLLam2
              .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
              + pLLam2.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
              + pLLam2.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake) };

  CATSLambdaParam pXiLam(Proton[1], Xi[1][1]);
  const double lam_pXim = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                CATSLambdaParam::Primary);
  const double lam_pXim_pXim1530 = pXiLam.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
  const double lam_pXim_fake = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                     CATSLambdaParam::Fake)
      + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
      + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);

  for (int vFrac_pL = 0; vFrac_pL < 3; vFrac_pL++) {
    std::cout << "=================== \n";
    std::cout << "vFrac_pL: " << vFrac_pL << std::endl;
    std::cout << " lam_pL: " << lam_pL[vFrac_pL] << " lam_pL_fake: "
              << lam_pL_fake[vFrac_pL] << " lam_pL_pS0: "
              << lam_pL_pS0[vFrac_pL] << " lam_pL_pXm: " << lam_pL_pXm[vFrac_pL]
              << std::endl;
    std::cout << "=================== \n";
  }
  std::cout << "lam_pXim: " << lam_pXim << " lam_pXim_pXim1530: "
            << lam_pXim_pXim1530 << " lam_pXim_fake:" << lam_pXim_fake
            << std::endl;

  const double GaussSourceSize = 1.2;
  //insert p-Sigma0 radius for different mT bins from r_core p-p &
  //effective gaussian fit of the  p-Sigma0 source including resonancess.
  std::vector<float> pSigma0Radii = { 1.473, 1.421, 1.368, 1.295, 1.220, 1.124 };
  const double pSigma0Radius = pSigma0Radii[imTBin];
  std::cout << "===========================\n";
  std::cout << "==pSigma0Radius: " << pSigma0Radius << "fm ==\n";
  std::cout << "===========================";
  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName(SigmaFileName.Data());
  CATSinput->ReadSigmaFile();

  TFile* OutFile = new TFile(
      TString::Format("%s/OutFileVarpL_%u.root", OutputDir.Data(), NumIter),
      "RECREATE");
  if (!OutFile) {
    return;
  }
  TDirectoryFile* CollOut = new TDirectoryFile(
      TString::Format("Out%u", NumIter), TString::Format("Out%u", NumIter));
  CollOut->SetName(TString::Format("Out%u", NumIter));
  TNtuple* ntResult = new TNtuple("ntResult", "ntResult", "NumIter:"
                                  "IterID:"
                                  "FemtoRegion:"
                                  "vModpL:"
                                  "PolBaseLine:"
                                  "lam_pp:"
                                  "lam_ppL:"
                                  "lam_pL:"
                                  "lam_pLS0:"
                                  "lam_pLXim:"
                                  "Source:"
                                  "Radius_pp:RadiusErr_pp:"
                                  "Stab_pp:StabErr_pp:"
                                  "p_a:p_aErr:"
                                  "p_b:p_bErr:"
                                  "p_c:p_cErr:"
                                  "chisqPerndf:"
                                  "BaseLineMin:"
                                  "BaseLineMax");

  Float_t ntBuffer[24];

  int uIter = 1;
  float total = 162;
  int counter = 1;
  int vFemReg;  //which femto region we use for pp (1 = default)
  int vMod_pL = iPotential;  //which pL function to use: //0=exact NLO (at the moment temporary it is Usmani); 1=Ledni NLO; 2=Ledni LO; 3=ESC08
  int vFrac_pL;  //fraction of protons coming from Lambda variation (1 = default)
  int iNorm = 1;

  TidyCats* tidy = new TidyCats();
  TCanvas* c1 = new TCanvas(TString::Format("out%u", NumIter));
  c1->SetCanvasSize(1920, 1280);
  StoreHist->SetLineWidth(3);
  StoreHist->SetLineColor(1);
  StoreHist->GetXaxis()->SetRangeUser(0, 600);
  c1->cd();
  StoreHist->DrawCopy();

  CollOut->Add(c1);
  StoreHist->GetXaxis()->SetRangeUser(0, 1000);
  CollOut->Add(StoreHist);

  CATS AB_pXim;
  tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins, kMin, kMax, FeeddownSource,
                             TidyCats::pHALQCD, 12);
  AB_pXim.KillTheCat();

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins, kMin, kMax,
                                 FeeddownSource);
  AB_pXim1530.KillTheCat();

  CATS AB_pL;
  for (vMod_pL = 1; vMod_pL < 4; ++vMod_pL) {
    if (vMod_pL == 1) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pUsmani);
      AB_pL.KillTheCat();
    } else if (vMod_pL == 2) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pNLOWF);
      AB_pL.KillTheCat();
    } else if (vMod_pL == 3) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pLOWF);
      AB_pL.KillTheCat();
    }
    for (vFemReg = 0; vFemReg < 3; ++vFemReg) {
      for (vFrac_pL = 0; vFrac_pL < 3; ++vFrac_pL) {
        for (int iBL = 0; iBL < 3; iBL++) {
          for (int BaselineSlope = 0; BaselineSlope < 3; ++BaselineSlope) {
            if (BaselineSlope == 1) {
              //no pol1 baseline.
              continue;
            }

            // Some computation here
            auto end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - start;

            std::cout
                << "\r Processing progress: "
                << TString::Format("%.1f %%", counter++ / total * 100.f).Data()
                << " elapsed time: " << elapsed_seconds.count()/60. << std::flush;
            TH1F* OliHisto_pp = (TH1F*) inFile->Get(HistppName.Data());
            if (!OliHisto_pp) {
              std::cout << HistppName.Data() << " Missing" << std::endl;
              return;
            }
            //!CHANGE PATH HERE

            const unsigned NumSourcePars =
                (TheSource == TidyCats::sLevy ? 2 : 1);

            //this way you define a correlation function using a CATS object.
            //needed inputs: num source/pot pars, CATS obj
            DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL);
            //this way you define a correlation function using Lednicky.
            //needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
            DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins, kMin, kMax,
                                            Lednicky_gauss_Sigma0);
            Ck_pSigma0->SetSourcePar(0, pSigma0Radius);
            DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
            Ck_pXim->SetSourcePar(0, 0.92);
            DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);
            Ck_pXim1530->SetSourcePar(0, 0.92);
            Ck_pL->Update();
            Ck_pSigma0->Update();
            Ck_pXim->Update();
            Ck_pXim1530->Update();
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
            DLM_CkDecomposition CkDec_pL("pLambda",
                                         TheSource == TidyCats::sLevy ? 3 : 4,
                                         *Ck_pL, CATSinput->GetSigmaFile(1));
            DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0,
            NULL);
            DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
                                           CATSinput->GetSigmaFile(3));
            DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530,
            NULL);
            if (!CATSinput->GetResFile(1)) {
              std::cout << "No Calib 1 \n";
              return;
            }
            if (!CATSinput->GetResFile(2)) {
              std::cout << "No Calib 2 \n";
              return;
            }
//
            if (TheSource == TidyCats::sLevy) {
              CkDec_pL.AddContribution(0, lam_pL_pXm.at(vFrac_pL),
                                       DLM_CkDecomposition::cFeedDown,
                                       &CkDec_pXim, CATSinput->GetResFile(2));
              CkDec_pL.AddContribution(
                  1,
                  1. - lam_pL.at(vFrac_pL) - lam_pL_pXm.at(vFrac_pL)
                      - lam_pL_fake.at(vFrac_pL),
                  DLM_CkDecomposition::cFeedDown);
              CkDec_pL.AddContribution(2, lam_pL_fake.at(vFrac_pL),
                                       DLM_CkDecomposition::cFake);  //0.03
            } else {
              CkDec_pL.AddContribution(0, lam_pL_pS0.at(vFrac_pL),
                                       DLM_CkDecomposition::cFeedDown,
                                       &CkDec_pSigma0,
                                       CATSinput->GetResFile(1));
              CkDec_pL.AddContribution(1, lam_pL_pXm.at(vFrac_pL),
                                       DLM_CkDecomposition::cFeedDown,
                                       &CkDec_pXim, CATSinput->GetResFile(2));
              CkDec_pL.AddContribution(
                  2,
                  1. - lam_pL.at(vFrac_pL) - lam_pL_pS0.at(vFrac_pL)
                      - lam_pL_pXm.at(vFrac_pL) - lam_pL_fake.at(vFrac_pL),
                  DLM_CkDecomposition::cFeedDown);
              CkDec_pL.AddContribution(3, lam_pL_fake.at(vFrac_pL),
                                       DLM_CkDecomposition::cFake);  //0.03
            }

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

            DLM_Fitter1* fitter;
            if (TheSource == TidyCats::sLevy) {
              fitter = new DLM_Fitter1(1);
              fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pL, kMin,
                                FemtoRegion[vFemReg], BaseLineRegion[iBL][0],
                                BaseLineRegion[iBL][1]);
//            fitter->AddSameSource("pXim", "pLambda", 2);
//            fitter->AddSameSource("pXim1530", "pLambda", 2);

              fitter->SetParameter("pLambda", DLM_Fitter1::p_sor0, 1.4, 0.5,
                                   2.5);
              fitter->SetParameter("pLambda", DLM_Fitter1::p_sor1, 1.7, 1., 2.);
            } else {
              fitter = new DLM_Fitter1(1);

              fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pL, kMin,
                                FemtoRegion[vFemReg], BaseLineRegion[iBL][0],
                                BaseLineRegion[iBL][1]);

//            fitter->AddSameSource("pLambda", "pLambda", 1);
//            fitter->AddSameSource("pSigma0", "pLambda", 1);
//            fitter->AddSameSource("pXim", "pLambda", 1);
//            fitter->AddSameSource("pXim1530", "pLambda", 1);

              fitter->SetParameter("pLambda", DLM_Fitter1::p_sor0, 1.4, 0.5,
                                   2.5);
            }
            fitter->SetOutputDir(OutputDir.Data());

            fitter->SetSeparateBL(0, false);
            fitter->SetParameter("pLambda", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
            if (BaselineSlope == 1) {
              fitter->SetParameter("pLambda", DLM_Fitter1::p_b, 1e-4, -2e-3,
                                   2e-3);
              fitter->FixParameter("pLambda", DLM_Fitter1::p_c, 0);
            } else if (BaselineSlope == 2) {
              fitter->SetParameter("pLambda", DLM_Fitter1::p_b, 1e-4, -2e-3,
                                   2e-3);
              fitter->SetParameter("pLambda", DLM_Fitter1::p_c, 1e-5, -1e-4,
                                   1e-4);
            } else {
              fitter->FixParameter("pLambda", DLM_Fitter1::p_b, 0);
              fitter->FixParameter("pLambda", DLM_Fitter1::p_c, 0);
            }

            fitter->FixParameter("pLambda", DLM_Fitter1::p_Cl, -1);

            CkDec_pL.Update();
            CkDec_pXim.Update();
            fitter->GoBabyGo();

            TGraph FitResult;
            FitResult.SetName(
                TString::Format("Graph_Var_%u_Iter_%u", NumIter, uIter));
            fitter->GetFitGraph(0, FitResult);
            c1->cd();
            TGraph* pointerFitResult = new TGraph(FitResult);
            pointerFitResult->SetName(
                TString::Format("Graph_Var_%u_Iter_%u", NumIter, uIter));
            pointerFitResult->SetLineWidth(2);
            pointerFitResult->SetLineColor(kRed);
            pointerFitResult->SetMarkerStyle(24);
            pointerFitResult->SetMarkerColor(kRed);
            pointerFitResult->SetMarkerSize(1);
            pointerFitResult->Draw("CP, same");

            double Chi2 = 0;
            unsigned EffNumBins = 0;
            if (BaselineSlope == 0) {
              EffNumBins = -2;  // radius and normalization
            } else if (BaselineSlope == 0) {
              EffNumBins = -3;  // radius, normalization and slope
            } else if (BaselineSlope == 0) {
              EffNumBins = -4;  // radius, normalization, slope and skewness
            }
            for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {

              double mom = AB_pL.GetMomentum(uBin);
              double dataY;
              double dataErr;
              double theoryX;
              double theoryY;

              if (mom > FemtoRegion[vFemReg])
                continue;

              FitResult.GetPoint(uBin, theoryX, theoryY);
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
              Chi2 += (dataY - theoryY) * (dataY - theoryY)
                  / (dataErr * dataErr);
              EffNumBins++;
            }

            ntBuffer[0] = NumIter;
            ntBuffer[1] = uIter;
            ntBuffer[2] = FemtoRegion[vFemReg];
            ntBuffer[3] = vMod_pL;
            ntBuffer[4] = BaselineSlope;
            ntBuffer[5] = 0.;
            ntBuffer[6] = 0.;
            ntBuffer[7] = lam_pL.at(vFrac_pL);
            ntBuffer[8] = lam_pL_pS0.at(vFrac_pL);
            ntBuffer[9] = lam_pL_pXm.at(vFrac_pL);
            ntBuffer[10] = TheSource;
            ntBuffer[11] = fitter->GetParameter("pLambda", DLM_Fitter1::p_sor0);
            ntBuffer[12] = fitter->GetParError("pLambda", DLM_Fitter1::p_sor0);
            ntBuffer[13] = fitter->GetParameter("pLambda", DLM_Fitter1::p_sor1);
            ntBuffer[14] = fitter->GetParError("pLambda", DLM_Fitter1::p_sor1);
            ntBuffer[15] = fitter->GetParameter("pLambda", DLM_Fitter1::p_a);
            ntBuffer[16] = fitter->GetParError("pLambda", DLM_Fitter1::p_a);
            ntBuffer[17] = fitter->GetParameter("pLambda", DLM_Fitter1::p_b);
            ntBuffer[18] = fitter->GetParError("pLambda", DLM_Fitter1::p_b);
            ntBuffer[19] = fitter->GetParameter("pLambda", DLM_Fitter1::p_c);
            ntBuffer[20] = fitter->GetParError("pLambda", DLM_Fitter1::p_c);
            ntBuffer[21] = Chi2 / EffNumBins;
            ntBuffer[22] = BaseLineRegion[iBL][0];
            ntBuffer[23] = BaseLineRegion[iBL][1];
            ntResult->Fill(ntBuffer);

            DLM_Histo<double>* CkpL_pS0 = CkDec_pL.GetChildContribution(
                (const unsigned int) 0, false);
            DLM_Histo<double>* CkpL_pXim = CkDec_pL.GetChildContribution(
                (const unsigned int) 1, false);
            TGraph* grCkpL_pS0 = new TGraph();
            grCkpL_pS0->SetName(
                TString::Format("grCkpL_pS0_Var_%u_Iter_%u", NumIter, uIter));

            TGraph* grCkpL_pS0Scaled = new TGraph();
            grCkpL_pS0Scaled->SetName(
                TString::Format("grCkpL_pS0Scaled_Var_%u_Iter_%u", NumIter,
                                uIter));

            TGraph* grCkpS0 = new TGraph();
            grCkpS0->SetName(
                TString::Format("grCkpS0_Var_%u_Iter_%u", NumIter, uIter));

            TGraph* grCkpL_pXi = new TGraph();
            grCkpL_pXi->SetName(
                TString::Format("grCkpL_pXi_Var_%u_Iter_%u", NumIter, uIter));

            TGraph* grCkpL_pXiScaled = new TGraph();
            grCkpL_pXiScaled->SetName(
                TString::Format("grCkpL_pXiScaled_Var_%u_Iter_%u", NumIter,
                                uIter));

            TGraph* grCkpXi = new TGraph();
            grCkpXi->SetName(
                TString::Format("grCkpXi_Var_%u_Iter_%u", NumIter, uIter));

            for (int iBin = 0; iBin < FitResult.GetN(); iBin++) {
              grCkpL_pS0->SetPoint(iBin, CkpL_pS0->GetBinCenter(0, iBin),
                                   CkpL_pS0->GetBinContent(iBin));
              grCkpL_pS0Scaled->SetPoint(
                  iBin,
                  CkpL_pS0->GetBinCenter(0, iBin),
                  1
                      + (CkpL_pS0->GetBinContent(iBin) - 1)
                          * lam_pL_pS0.at(vFrac_pL));
              grCkpS0->SetPoint(
                  iBin, CkpL_pS0->GetBinCenter(0, iBin),
                  CkDec_pSigma0.EvalCk(CkpL_pS0->GetBinCenter(0, iBin)));

              grCkpL_pXi->SetPoint(iBin, CkpL_pXim->GetBinCenter(0, iBin),
                                   CkpL_pXim->GetBinContent(iBin));
              grCkpL_pXiScaled->SetPoint(
                  iBin,
                  CkpL_pXim->GetBinCenter(0, iBin),
                  1
                      + (CkpL_pXim->GetBinContent(iBin) - 1)
                          * lam_pL_pXm.at(vFrac_pL));
              grCkpXi->SetPoint(
                  iBin, CkpL_pXim->GetBinCenter(0, iBin),
                  CkDec_pXim.EvalCk(CkpL_pXim->GetBinCenter(0, iBin)));
            }

            TList* outList = new TList();
            outList->SetOwner();
            outList->SetName(
                TString::Format("Graph_Var_%u_iter_%u", NumIter, uIter));
            outList->Add(pointerFitResult);
            outList->Add(grCkpL_pS0);
            outList->Add(grCkpL_pS0Scaled);
            outList->Add(grCkpS0);

            outList->Add(grCkpL_pXi);
            outList->Add(grCkpL_pXiScaled);
            outList->Add(grCkpXi);
            CollOut->Add(outList);

            uIter++;

            delete Ck_pL;
            delete Ck_pSigma0;
            delete Ck_pXim;
            delete Ck_pXim1530;
            delete fitter;
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
}

int main(int argc, char *argv[]) {
  FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                  atoi(argv[5]), argv[6], argv[7], argv[8]);
  return 0;
}


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

#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TCanvas.h"

#include <iostream>
#include "stdlib.h"
#include <chrono>
#include <ctime>

void FitPPVariations(const unsigned& NumIter, int system, int source,
                     TString InputFile, TString HistoName, TString OutputDir) {
  auto start = std::chrono::system_clock::now();
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  //What source to use: 0 = Gauss; 1=Resonance; 2=Levy
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

  const unsigned NumMomBins = 60;
  const double kMin = StoreHist->GetXaxis()->GetXmin();
  const double kMax = kMin + StoreHist->GetXaxis()->GetBinWidth(1) * NumMomBins;  //(4 is the bin width)

  //if you modify you may need to change the CATS ranges somewhere below
  double FemtoRegion[3];
  FemtoRegion[0] = 350;
  FemtoRegion[1] = 375;
  FemtoRegion[2] = 400;

  TidyCats::Sources TheSource;
  TidyCats::Sources FeeddownSource;
  bool RandomizedEmission = false; 
  if (source == 0) {
    TheSource = TidyCats::sGaussian;
    FeeddownSource = TheSource;
  } else if (source == 1) {
    TheSource = TidyCats::sResonance;
    FeeddownSource = TidyCats::sGaussian;
    RandomizedEmission = true; //in case we have a resonance source this steers the way the emission of resonances is handeled. 
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
  CATSLambdaParam ppLam0(Proton[0], Proton[0], true);
  CATSLambdaParam ppLam1(Proton[1], Proton[1], true);
  CATSLambdaParam ppLam2(Proton[2], Proton[2], true);

  const std::vector<double> lam_pp =
      { ppLam0.GetLambdaParam(CATSLambdaParam::Primary,
                              CATSLambdaParam::Primary), ppLam1.GetLambdaParam(
          CATSLambdaParam::Primary, CATSLambdaParam::Primary), ppLam2
          .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary) };
  const std::vector<double> lam_pp_pL =
      { ppLam0.GetLambdaParam(CATSLambdaParam::Primary,
                              CATSLambdaParam::FeedDown, 0, 0), ppLam1
          .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown,
                          0, 0), ppLam2.GetLambdaParam(
          CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0) };
  const std::vector<double> lam_pp_fake =
      {
          (ppLam0.GetLambdaParam(CATSLambdaParam::Primary,
                                 CATSLambdaParam::Fake)
              + ppLam0.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
              + ppLam0.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake)), (ppLam1
              .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
              + ppLam1.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
              + ppLam1.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake)), (ppLam2
              .GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
              + ppLam2.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
              + ppLam2.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake)) };
  CATSLambdaParam pLLam(Proton[1], Lambda[1][1]);
  const double lam_pL = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
                                             CATSLambdaParam::Primary);
  const double lam_pL_pS0 = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                 CATSLambdaParam::FeedDown, 0,
                                                 0);
  const double lam_pL_pXm = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                 CATSLambdaParam::FeedDown, 0,
                                                 1);
  const double lam_pL_fake = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                  CATSLambdaParam::Fake)
      + pLLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
      + pLLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);
  CATSLambdaParam pXiLam(Proton[1], Xi[1][1]);
  const double lam_pXim = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                CATSLambdaParam::Primary);
  const double lam_pXim_pXim1530 = pXiLam.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
  const double lam_pXim_fake = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
                                                     CATSLambdaParam::Fake)
      + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
      + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);

  for (int vFrac_pp_pL = 0; vFrac_pp_pL < 3; vFrac_pp_pL++) {
    std::cout << "lam_pp: " << lam_pp.at(vFrac_pp_pL) << " lam_pp_pL: "
              << lam_pp_pL.at(vFrac_pp_pL) << " lam_pp_fake: "
              << lam_pp_fake.at(vFrac_pp_pL) << std::endl;
  }
  std::cout << "lam_pL_pXm: " << lam_pL_pXm << " lam_pL: " << lam_pL
            << " lam_pL_fake: " << lam_pL_fake << std::endl;
  std::cout << "lam_pXim: " << lam_pXim << " lam_pXim_pXim1530: "
            << lam_pXim_pXim1530 << " lam_pXim_fake:" << lam_pXim_fake
            << std::endl;
  const double GaussSourceSize = 1.2;

  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName(SigmaFileName.Data());
  CATSinput->ReadSigmaFile();

  TFile* OutFile = new TFile(
      TString::Format("%s/OutFileVarpp_%u.root", OutputDir.Data(), NumIter),
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
                                  "chisqPerndf");

  Float_t ntBuffer[22];

  int uIter = 1;
  float total = TheSource == TidyCats::sLevy ? 54 : 81;
  int counter = 1;
  int vFemReg;  //which femto region we use for pp (1 = default)
  int vMod_pL = 1;  //which pL function to use: //0=exact NLO (at the moment temporary it is Usmani); 1=Ledni NLO; 2=Ledni LO; 3=ESC08
  int vFrac_pp_pL;  //fraction of protons coming from Lambda variation (1 = default)
  int iNorm = 1;

  TidyCats* tidy = new TidyCats();
  TCanvas* c1 = new TCanvas(TString::Format("out%u", NumIter));
  c1->SetCanvasSize(1920, 1280);
  StoreHist->SetLineWidth(3);
  StoreHist->SetLineColor(1);
  StoreHist->GetXaxis()->SetRangeUser(0, 400);
  c1->cd();
  StoreHist->DrawCopy();
  CollOut->Add(c1);
  StoreHist->GetXaxis()->SetRangeUser(0, 1000);
  CollOut->Add(StoreHist);
  
  CATS AB_pp;
  tidy->GetCatsProtonProton(&AB_pp, NumMomBins, kMin, kMax, TheSource);
  if (TheSource == TidyCats::sResonance) {
    const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass() * 1000;
    const double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    DLM_CleverMcLevyReso* source = tidy->GetSourceProtonProton(); 
    source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
		      massPion,false,false,RandomizedEmission?DLM_CleverMcLevyReso::rdtRandom:DLM_CleverMcLevyReso::rdtBackwards);
    source->SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
		      massPion,false,false,RandomizedEmission?DLM_CleverMcLevyReso::rdtRandom:DLM_CleverMcLevyReso::rdtBackwards);
    if (RandomizedEmission) {
      std::cout << "Sir, Emission will be fully randomized, commencing countdown ... 3 ....\n"; 
      const char* PhiFile = "DimiPhi_pp_HM.root"; 
      DLM_Histo<double>* HISTO = tidy->ConvertThetaAngleHisto(TString::Format("~/cernbox/WaveFunctions/ThetaDist/%s",PhiFile).Data(),"h_rkAngle_Mom2",400,600);
      source->SetUpResoEmission(0,0,HISTO);
      source->SetUpResoEmission(1,0,HISTO);
    }  
  }
  AB_pp.KillTheCat();
  if (tidy->GetSourceProtonProton()) {
    TGraph* SourceDist = new TGraph();
    SourceDist->SetName(TString::Format("SourceDist_NumIter_%i", NumIter));
    for (int iRad = 0; iRad < 200; ++iRad) {
      std::cout << "\r Source progress: "
          << TString::Format("%.1f %%", (iRad + 1) / 200 * 100.f).Data()
          << std::flush;
      double rad = 0.04 * iRad;
      double pars[5] = { 0, rad, 0, 1.2, 1.7 };
      SourceDist->SetPoint(iRad, rad, tidy->GetSourceProtonProton()->Eval(pars));
    }
    std::cout << std::endl;
    CollOut->Add(SourceDist);
  }
  CATS AB_pXim;
  tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins, kMin, kMax, FeeddownSource,
                             TidyCats::pHALQCD, 12);
  AB_pXim.KillTheCat();

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins, kMin, kMax,
                                 FeeddownSource);
  AB_pXim1530.KillTheCat();

  for (vMod_pL = TheSource == TidyCats::sLevy ? 1 : 0; vMod_pL < 3; ++vMod_pL) {
    TidyCats::pLPot PLpot;
    CATS AB_pL;
    if (vMod_pL == 1) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, FeeddownSource,
                                TidyCats::pUsmani);
      AB_pL.KillTheCat();
    } else if (vMod_pL == 2) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, FeeddownSource,
                                TidyCats::pNLOWF);
      AB_pL.KillTheCat();
    }
    for (vFemReg = 0; vFemReg < 3; ++vFemReg) {
      for (vFrac_pp_pL = 0; vFrac_pp_pL < 3; ++vFrac_pp_pL) {
        for (int BaselineSlope = 0; BaselineSlope < 3; ++BaselineSlope) {
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
//          CATSinput->AddSystematics("C2totalsysPP.root", OliHisto_pp);

          //!CHANGE PATH HERE

          const unsigned NumSourcePars = (TheSource == TidyCats::sLevy ? 2 : 1);

          //this way you define a correlation function using a CATS object.
          //needed inputs: num source/pot pars, CATS obj
          DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp);
          DLM_Ck* Ck_pL =
              vMod_pL == 0 ?
                  new DLM_Ck(1, 4, NumMomBins, kMin, kMax,
                             Lednicky_SingletTriplet) :
                  new DLM_Ck(NumSourcePars, 0, AB_pL);
          //this way you define a correlation function using Lednicky.
          //needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
          DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins, kMin, kMax,
                                          Lednicky_gauss_Sigma0);
          Ck_pSigma0->SetSourcePar(0, GaussSourceSize);
          DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
          DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);
          if (vMod_pL == 0) {
            Ck_pL->SetPotPar(0, 2.91);
            Ck_pL->SetPotPar(1, 2.78);
            Ck_pL->SetPotPar(2, 1.54);
            Ck_pL->SetPotPar(3, 2.72);
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
          DLM_CkDecomposition CkDec_pL("pLambda",
                                       TheSource == TidyCats::sLevy ? 3 : 4,
                                       *Ck_pL, CATSinput->GetSigmaFile(1));
          DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0,
          NULL);
          DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
                                         CATSinput->GetSigmaFile(3));
          DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530,
          NULL);
          if (!CATSinput->GetResFile(0)) {
            std::cout << "No Calib 0 \n";
            return;
          }

          CkDec_pp.AddContribution(0, lam_pp_pL.at(vFrac_pp_pL),
                                   DLM_CkDecomposition::cFeedDown, &CkDec_pL,
                                   CATSinput->GetResFile(0));
          CkDec_pp.AddContribution(
              1,
              1. - lam_pp.at(vFrac_pp_pL) - lam_pp_pL.at(vFrac_pp_pL)
                  - lam_pp_fake.at(vFrac_pp_pL),
              DLM_CkDecomposition::cFeedDown);
          CkDec_pp.AddContribution(2, lam_pp_fake.at(vFrac_pp_pL),
                                   DLM_CkDecomposition::cFake);  //0.02

          if (!CATSinput->GetResFile(1)) {
            std::cout << "No Calib 1 \n";
            return;
          }
          if (!CATSinput->GetResFile(2)) {
            std::cout << "No Calib 2 \n";
            return;
          }

          if (TheSource == TidyCats::sLevy) {
            CkDec_pL.AddContribution(0, lam_pL_pXm,
                                     DLM_CkDecomposition::cFeedDown,
                                     &CkDec_pXim, CATSinput->GetResFile(2));
            CkDec_pL.AddContribution(1, 1. - lam_pL - lam_pL_pXm - lam_pL_fake,
                                     DLM_CkDecomposition::cFeedDown);
            CkDec_pL.AddContribution(2, lam_pL_fake,
                                     DLM_CkDecomposition::cFake);  //0.03
          } else {
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
          }

          if (!CATSinput->GetResFile(3)) {
            std::cout << "No Calib 3 \n";
            return;
          }
          CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
                                     DLM_CkDecomposition::cFeedDown,
                                     &CkDec_pXim1530, CATSinput->GetResFile(3));  //from Xi-(1530)
          CkDec_pXim.AddContribution(
              1, 1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
              DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
          CkDec_pXim.AddContribution(2, lam_pXim_fake,
                                     DLM_CkDecomposition::cFake);

          DLM_Fitter1* fitter;
          if (TheSource == TidyCats::sLevy) {
            fitter = new DLM_Fitter1(1);
            fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp, kMin,
                              FemtoRegion[vFemReg], FemtoRegion[vFemReg],
                              FemtoRegion[vFemReg]);
            fitter->AddSameSource("pLambda", "pp", 2);
            fitter->AddSameSource("pXim", "pp", 2);
            fitter->AddSameSource("pXim1530", "pp", 2);

            fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.4, 0.5, 2.5);
            fitter->SetParameter("pp", DLM_Fitter1::p_sor1, 1.7, 1., 2.);
          } else {
            fitter = new DLM_Fitter1(1);

            fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp, kMin,
                              FemtoRegion[vFemReg], FemtoRegion[vFemReg],
                              FemtoRegion[vFemReg]);

            fitter->AddSameSource("pLambda", "pp", 1);
            fitter->AddSameSource("pSigma0", "pp", 1);
            fitter->AddSameSource("pXim", "pp", 1);
            fitter->AddSameSource("pXim1530", "pp", 1);

            fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.4, 0.5, 2.5);
          }
          fitter->SetOutputDir(OutputDir.Data());

          fitter->SetSeparateBL(0, false);
          fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
          if (BaselineSlope == 1) {
            fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, -2e-3, 2e-3);
          } else if (BaselineSlope == 2) {
            fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, -2e-3, 2e-3);
            fitter->SetParameter("pp", DLM_Fitter1::p_c, 1e-5, -1e-4, 1e-4);
          } else {
            fitter->FixParameter("pp", DLM_Fitter1::p_b, 0);
            fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);
          }

          fitter->FixParameter("pp", DLM_Fitter1::p_Cl, -1);

          CkDec_pp.Update();
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
          for (unsigned uBin = 0; uBin < 50; uBin++) {

            double mom = AB_pp.GetMomentum(uBin);
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
              std::cout << dataErr << '\t' << "WARNING POINT NOT CONSIDERED \n";
              continue;
            }
            Chi2 += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
            EffNumBins++;
          }

          ntBuffer[0] = NumIter;
          ntBuffer[1] = uIter;
          ntBuffer[2] = FemtoRegion[vFemReg];
          ntBuffer[3] = vMod_pL;
          ntBuffer[4] = BaselineSlope;
          ntBuffer[5] = lam_pp.at(vFrac_pp_pL);
          ntBuffer[6] = lam_pp_pL.at(vFrac_pp_pL);
          ntBuffer[7] = lam_pL;
          ntBuffer[8] = lam_pL_pS0;
          ntBuffer[9] = lam_pL_pXm;
          ntBuffer[10] = TheSource;
          ntBuffer[11] = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
          ntBuffer[12] = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
          ntBuffer[13] = fitter->GetParameter("pp", DLM_Fitter1::p_sor1);
          ntBuffer[14] = fitter->GetParError("pp", DLM_Fitter1::p_sor1);
          ntBuffer[15] = fitter->GetParameter("pp", DLM_Fitter1::p_a);
          ntBuffer[16] = fitter->GetParError("pp", DLM_Fitter1::p_a);
          ntBuffer[17] = fitter->GetParameter("pp", DLM_Fitter1::p_b);
          ntBuffer[18] = fitter->GetParError("pp", DLM_Fitter1::p_b);
          ntBuffer[19] = fitter->GetParameter("pp", DLM_Fitter1::p_c);
          ntBuffer[20] = fitter->GetParError("pp", DLM_Fitter1::p_c);
          ntBuffer[21] = Chi2 / EffNumBins;
          ntResult->Fill(ntBuffer);

          TList* outList = new TList();
          outList->SetOwner();
          outList->SetName(
              TString::Format("Graph_Var_%u_iter_%u", NumIter, uIter));
          outList->Add(pointerFitResult);
          CollOut->Add(outList);

          uIter++;

          delete Ck_pp;
          delete Ck_pL;
          delete Ck_pSigma0;
          delete Ck_pXim;
          delete Ck_pXim1530;
          delete fitter;
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
  FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), argv[4], argv[5],
                  argv[6]);
  return 0;
}


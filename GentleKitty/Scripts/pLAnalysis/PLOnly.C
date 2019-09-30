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

void FitPPVariations(const unsigned& NumIter, int imTBin, int system,
                     int source, int iPotential, unsigned int iAngDist, int iRange,
                     TString InputFile, TString HistoName, TString OutputDir) {
  //What source to use: 0 = Gauss; 1=Resonance; 2=Levy
  auto start = std::chrono::system_clock::now();
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  
  //Setup Various Inputs
  
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
  std::vector<double> ProSecondary = {0.78, 0.8,  0.8,  0.8,  0.82, 0.84, 0.87};
  std::vector<double> LamPurity = {0.92, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95 };
  std::vector<double> LamSecondary = { 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78};
  std::vector<float> mTValues = { 1.21236, 1.28964, 1.37596, 1.54074, 1.75601, 2.25937 }; 

  std::cout << "SYSTEM: " << system << std::endl;
  
  //Setup the input vales for the Lambda parameters & Specific to the system  
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

  //Setup the input resolution and smearing
  
  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName(SigmaFileName.Data());
  CATSinput->ReadSigmaFile();

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
  
  //Calculate the Lambda Parameters and their variations 

  std::vector<double> Variation = { 0.95, 1.0, 1.05 };

  std::vector<Particle> Proton;  // 1) variation of the Prim Fraction 2) variation of the Secondary Comp.
  std::vector<Particle> Lambda;  // 1) variation of Lambda/Sigma Ratio, 2) variation of Xi0/Xim Ratio
  std::vector<Particle> Xi;  //1) variation of dN/dy Omega 2) variation of dN/dy Xi1530

  for (auto itFrac : Variation) {
    const double varPrimProton = itFrac* PrimProton; 
    const double varSecProton = (1-(double)varPrimProton);
    for (auto itComp : Variation) {
      double Composition = 0.7*itComp; 
      double varSecFracLamb = Composition*varSecProton; 
      double varSecFracSigma = 1. - varPrimProton - varSecFracLamb;
      Proton.push_back(Particle(PurityProton , varPrimProton, { varSecFracLamb,
	      varSecFracSigma }));
      double sum = varPrimProton+varSecFracLamb+varSecFracSigma; 
    }
  }

  for (auto itSLRatio : Variation) {
    double LamSigProdFraction = 3 * itSLRatio / 4. < 1 ? 3 * itSLRatio / 4. : 1;
    double PrimLambda = LamSigProdFraction * PrimLambdaAndSigma; 
    double SecSigLambda = (1. - LamSigProdFraction) * PrimLambdaAndSigma;  // decay probability = 100%!
    for (auto itXi02mRatio : Variation) {
      double SecXimLambda = itXi02mRatio * SecLambda / 2.;
      double SecXi0Lambda0 = (1 - itXi02mRatio / 2.) * SecLambda;
      Lambda.push_back(Particle(PurityLambda, PrimLambda, { SecSigLambda,
	      SecXimLambda, SecXi0Lambda0 }));
    }
  }
  
  for (auto itXimOmega : Variation) {
    for (auto itXi1530 : Variation) {
      double XiNormalization = 1 + itXimOmega * OmegamXimProdFraction * OmegamXim_BR
	+ itXi1530 * (Xi01530XimProdFraction * Xi01530Xim_BR
		      + Xim1530XimProdFraction * Xim1530Xim_BR);
      double SecOmegaXim = itXimOmega * OmegamXimProdFraction * OmegamXim_BR
	/ (double) XiNormalization;
      double SecXi01530Xim = itXi1530 * Xi01530XimProdFraction * Xi01530Xim_BR
	/ (double) XiNormalization;
      double SecXim1530Xim = itXi1530 * Xim1530XimProdFraction * Xim1530Xim_BR
	/ (double) XiNormalization;
      double PrimXim = 1. / (double) XiNormalization;
      Xi.push_back(Particle(PurityXi, PrimXim, { SecOmegaXim,
	      SecXi01530Xim, SecXim1530Xim }));
    }
  }
  //Sigma/Lam variations
  CATSLambdaParam pLLam0(Proton[0], Lambda[0]);
  CATSLambdaParam pLLam1(Proton[0], Lambda[4]);
  CATSLambdaParam pLLam2(Proton[0], Lambda[8]);
  CATSLambdaParam pLLam3(Proton[4], Lambda[0]);
  CATSLambdaParam pLLam4(Proton[4], Lambda[4]);
  CATSLambdaParam pLLam5(Proton[4], Lambda[8]);
  CATSLambdaParam pLLam6(Proton[8], Lambda[0]);
  CATSLambdaParam pLLam7(Proton[8], Lambda[4]);
  CATSLambdaParam pLLam8(Proton[8], Lambda[8]);
  
  const std::vector<double> lam_pL = {
    pLLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    pLLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary)
  };
  const std::vector<double> lam_pL_pS0 = {
    pLLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    pLLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
  };
  
  const std::vector<double> lam_pL_pXm = {
    pLLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
    pLLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1),
  };
  
  const std::vector<double> lam_pL_fake = {
    pLLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam0.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam0.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam1.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam1.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam2.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam2.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam3.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam3.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam4.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam4.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam5.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam5.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam6.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam6.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam7.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam7.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake),
    pLLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
    pLLam8.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Primary) +
    pLLam8.GetLambdaParam(CATSLambdaParam::Fake,    CATSLambdaParam::Fake)
  };

  CATSLambdaParam pXiLam(Proton[4], Xi[4]);
  const double lam_pXim =
    pXiLam.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary);
  const double lam_pXim_pXim1530 =
    pXiLam.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
  const double lam_pXim_fake =
    pXiLam.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
    + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
    + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);
  
  // for (int vFrac_pL = 0; vFrac_pL < 9; vFrac_pL++) {
  //   std::cout << "=================== \n";
  //   std::cout << "vFrac_pL: " << vFrac_pL << std::endl;
  //   std::cout << " lam_pL: " << lam_pL[vFrac_pL] << " lam_pL_fake: "
  //             << lam_pL_fake[vFrac_pL] << " lam_pL_pS0: "
  //             << lam_pL_pS0[vFrac_pL] << " lam_pL_pXm: " << lam_pL_pXm[vFrac_pL]
  //             << std::endl;
  //   std::cout << "=================== \n";
  // }
  // std::cout << "lam_pXim: " << lam_pXim << " lam_pXim_pXim1530: "
  //           << lam_pXim_pXim1530 << " lam_pXim_fake:" << lam_pXim_fake
  //           << std::endl;

  //insert p-Sigma0 radius for different mT bins from r_core p-p &
  //effective gaussian fit of the  p-Sigma0 source including resonances.
  std::vector<float> pSigma0Radii = { 1.473, 1.421, 1.368, 1.295, 1.220, 1.124 };
  const double pSigma0Radius = pSigma0Radii[imTBin];
  std::cout << "===========================\n";
  std::cout << "==pSigma0Radius: " << pSigma0Radius << "fm ==\n";
  std::cout << "===========================\n";

  TFile* OutFile = new TFile(
			     TString::Format("%s/OutFileVarpL_%u.root", OutputDir.Data(), NumIter),
			     "RECREATE");
  if (!OutFile) {
    return;
  }
  double AvgRadius = 0.;
  
  float total = 324;
  int numIter = NumIter; 
  int uIter = 1; 
  int vFemReg;  //which femto region we use for pp (1 = default)
  float thismT = mTValues[imTBin]; 
  float FemtoFitMax = 0;
  float BaseLineMin = 0;
  float BaseLineMax = 0; 
  unsigned int vMod_pL;
  unsigned int BaseLine;
  float pa; 
  float pb;  
  float pc; 
  int vFrac_pL;  //fraction of protons coming from Lambda variation (1 = default)
  float lmb_pp = 0; 
  float lmb_ppL = 0;
  float lmb_pL;
  float lmb_pLS0;
  float lmb_pLXim;
  float ResultRadius;
  float ResultRadiusErr;
  float FeedDownRadius = pSigma0Radius; 
  float Stability;
  float outChiSqNDF;
  TGraph* pointerFitResult = nullptr;

  TTree* outTree = new TTree("ppTree","ppTree");
  outTree->Branch("DataVarID", &numIter, "DataVarID/i"); 
  outTree->Branch("FitVarID", &uIter, "FitVarID/i"); 
  outTree->Branch("imT",&imTBin, "imT/i");
  outTree->Branch("mTValue",&thismT, "mTValue/F");
  outTree->Branch("iAng",&iAngDist,"iAng/i");
  outTree->Branch("iRange",&iRange,"iRange/i");
  outTree->Branch("FemtoFitMax",&FemtoFitMax,"FemtoFitMax/F");
  outTree->Branch("BaseLineMin",&BaseLineMin,"BaseLineMin/F");
  outTree->Branch("BaseLineMax",&BaseLineMax,"BaseLineMax/F");
  outTree->Branch("ModPL",&vMod_pL,"ModPL/i");
  outTree->Branch("PolBL",&BaseLine,"PolBL/i");
  outTree->Branch("pa",&pa,"pa/F");
  outTree->Branch("pb",&pb,"pb/F");
  outTree->Branch("pc",&pc,"pc/F");
  outTree->Branch("lam_pp",&lmb_pp,"lam_pp/F");
  outTree->Branch("lam_ppL",&lmb_ppL,"lam_ppL/F");
  outTree->Branch("lam_pL",&lmb_pL,"lam_pL/F");
  outTree->Branch("lam_pLS0",&lmb_pLS0,"lam_pLS0/F");
  outTree->Branch("lam_pLXim",&lmb_pLXim,"lam_pLXim/F");
  outTree->Branch("Source",&TheSource,"Source/i");
  outTree->Branch("Radius",&ResultRadius,"Radius/F");
  outTree->Branch("RadiusErr",&ResultRadiusErr,"RadiusErr/F");
  outTree->Branch("FeedDownRadius",&FeedDownRadius,"FeedDownRadius/F");
  outTree->Branch("Stab",&Stability,"Stab/F");
  outTree->Branch("chiSqNDF",&outChiSqNDF,"chiSqNDF/F");
  outTree->Branch("CorrHist","TH1F",&StoreHist,sizeof(TH1F));
  outTree->Branch("FitResult","TGraph",&pointerFitResult,sizeof(TGraph));
  
  TidyCats* tidy = new TidyCats();
  
  CATS AB_pXim;
  tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins, kMin, kMax, FeeddownSource,
                             TidyCats::pHALQCD, 12);
  AB_pXim.SetAnaSource(0, pSigma0Radius);
  AB_pXim.KillTheCat();
  std::cout << "pXim CaT\n"; 
  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins, kMin, kMax,
                                 FeeddownSource);
  AB_pXim1530.SetAnaSource(0, pSigma0Radius);
  AB_pXim1530.KillTheCat();
  std::cout << "pXim1530 CaT\n"; 
  CATS AB_pL;
  for (vMod_pL = 1; vMod_pL < 4; ++vMod_pL) {
    if (vMod_pL == 1) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pUsmani);
    } else if (vMod_pL == 2) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pNLOWF);
    } else if (vMod_pL == 3) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, TheSource,
                                TidyCats::pLOWF);
    }
    if (TheSource == TidyCats::sResonance) {
      const double massProton = TDatabasePDG::Instance()->GetParticle(2212)
	->Mass() * 1000;
      const double massLambda = TDatabasePDG::Instance()->GetParticle(3122)
	->Mass() * 1000;
      const double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass()
	* 1000;
      DLM_CleverMcLevyReso* source = tidy->GetSourceProtonLambda();
      if (iAngDist == 0) {
	source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
			  massPion);
	source->SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, massLambda,
			  massPion);
    
	std::cout << "Emission source is boring and back to back \n"; 
      } else if (iAngDist == 1) {
	source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
			  massPion, false, false,
			  DLM_CleverMcLevyReso::rdtRandomBackwards);
	source->SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, massLambda,
			  massPion, false, false,
			  DLM_CleverMcLevyReso::rdtRandomBackwards);
 	std::cout
	  << "Sir, Emission will be fully randomized with a tad bit above 90 degrees, commencing countdown ... 3 ....\n";
      } else if(iAngDist > 1) { 
	source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
			  massPion, false, false,
			  DLM_CleverMcLevyReso::rdtRandom);
	source->SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, massLambda,
			  massPion, false, false,
			  DLM_CleverMcLevyReso::rdtRandom);
	std::cout
	  << "Sir, Emission will be fully randomized, commencing countdown ... 3 ..2 .., but wait, how random is random ...? \n";
	if (iAngDist > 2) {
	  const char* PhiFile;
	  const char* PhiHistoName;
	  int PhiRebin = 1;
	  int PhiConversion = 1; 
	  double RangeProtonMin, RangeProtonMax;
	  double RangeLambdaMin, RangeLambdaMax;
	  if (iAngDist == 3) {
	    PhiFile = "DimiPhi_pp_HM.root";
	    PhiHistoName = "h_rkAngle_Mom2"; 
	  } else if (iAngDist == 4){
	    PhiFile = "DimiPhi_pLambda_HM.root";
	    PhiHistoName = "h_rkAngle_Mom2"; 
	  } else if (iAngDist == 5) {
	    PhiFile = "DimiPhi_LambdaLambda_HM.root";
	    PhiHistoName = "h_rkAngle_Mom2"; 
	  } else if (iAngDist == 6) {
	    PhiFile = "DimiPhi_pXim_HM.root";
	    PhiHistoName = "h_rkAngle_Mom2"; 
	  } else if (iAngDist == 7) {
	    PhiFile = "Output60.root";//60 degree cut off
	    PhiHistoName = "fAngleDistTwo_boosted";
	    PhiRebin = 2;
	    PhiConversion = 1000; 
	  } else if (iAngDist == 8) {
	    PhiFile = "Output75.root";//75 degree cut off 
	    PhiHistoName = "fAngleDistTwo_boosted";
	    PhiRebin = 2;
	    PhiConversion = 1000; 
	  } else if (iAngDist == 9) {
	    PhiFile = "Output90.root";//90 degree cut off
	    PhiHistoName = "fAngleDistTwo_boosted";
	    PhiRebin = 2;
	    PhiConversion = 1000; 
	  } else if (iAngDist == 10) {
	    PhiFile = "Output180.root";//no cut off 
	    PhiHistoName = "fAngleDistTwo_boosted";
	    PhiRebin = 2;
	    PhiConversion = 1000; 
	  } else {
	    std::cout << "Option iAngDist == " << iAngDist << " not viable \n";
	    return;
	  }
	  if (iRange == 0) {
	    RangeProtonMin = 400/(double)PhiConversion;
	    RangeProtonMax = 600/(double)PhiConversion;

	    RangeLambdaMin = 250/(double)PhiConversion;
	    RangeLambdaMax = 450/(double)PhiConversion;
	  } else if (iRange == 1) {
	    RangeProtonMin = 450/(double)PhiConversion;
	    RangeProtonMax = 550/(double)PhiConversion;

	    RangeLambdaMin = 300/(double)PhiConversion;
	    RangeLambdaMax = 400/(double)PhiConversion;
	  } else if (iRange == 2) {
	    RangeProtonMin = 250/(double)PhiConversion;
	    RangeProtonMax = 650/(double)PhiConversion;

	    RangeLambdaMin = 200/(double)PhiConversion;
	    RangeLambdaMax = 500/(double)PhiConversion;
	  } else {
	    std::cout << "Option iRange == " << iRange
		      << " is not a viable option \n";
	    return;
	  }
	  std::cout << "Using file: " << PhiFile << " in the range ProtonMin: "
		    << RangeProtonMin << " to Proton Max: " << RangeProtonMax
		    << std::endl;
	  DLM_Histo<double>* HISTO_PROTON = tidy->ConvertThetaAngleHisto(TString::Format("~/cernbox/WaveFunctions/ThetaDist/%s", PhiFile).Data(),TString::Format("%s",PhiHistoName).Data(), RangeProtonMin, RangeProtonMax, false, PhiRebin);
	  DLM_Histo<double>* HISTO_LAMBDA = tidy->ConvertThetaAngleHisto(TString::Format("~/cernbox/WaveFunctions/ThetaDist/%s", PhiFile).Data(),TString::Format("%s",PhiHistoName).Data(), RangeLambdaMin, RangeLambdaMax, false, PhiRebin);
	    
	  source->SetUpResoEmission(0, 0, HISTO_PROTON);
	  source->SetUpResoEmission(1, 0, HISTO_LAMBDA);
      	}
      } else {
	std::cout << "Case for iAng == :" << iAngDist << " not implemented, exiting. \n";
	return; 
      }
    }
    AB_pL.SetNotifications(CATS::nError);
    std::cout << "Killing pL Cat \n";
    AB_pL.KillTheCat();
    for (vFemReg = 0; vFemReg < 3; ++vFemReg) {
      FemtoFitMax = FemtoRegion[vFemReg];
      for (vFrac_pL = 0; vFrac_pL < 9; ++vFrac_pL) {
	lmb_pL = lam_pL.at(vFrac_pL); 
	lmb_pLS0 = lam_pL_pS0.at(vFrac_pL);
	lmb_pLXim = lam_pL_pXm.at(vFrac_pL);
        for (int iBL = 0; iBL < 3; iBL++) {
	  BaseLineMin = FemtoRegion[iBL];
	  BaseLineMax = FemtoRegion[iBL];
          for (BaseLine = 0; BaseLine < 3; ++BaseLine) {
            if (BaseLine == 1) {
              //no pol1 baseline.
              continue;
            }
            // Some computation here
            auto end = std::chrono::system_clock::now();
	    double runningAvg = uIter == 1?1:AvgRadius/(double)(uIter-1);
            std::chrono::duration<double> elapsed_seconds = end - start;
            std::cout
	      << "\r Processing progress: "
	      << TString::Format("%.1f %%", uIter / total * 100.f).Data()
	      << " elapsed time: " << elapsed_seconds.count() / 60. << " avg. Radius: " <<   runningAvg
	      << std::flush;
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
	    //Ck_pL->SetSourcePar(0,1.3); 
            DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins, kMin, kMax,
                                            Lednicky_gauss_Sigma0);
            Ck_pSigma0->SetSourcePar(0, pSigma0Radius);
            DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
            Ck_pXim->SetSourcePar(0, pSigma0Radius);
            DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);
            Ck_pXim1530->SetSourcePar(0, pSigma0Radius);
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
            if (BaseLine == 1) {
              fitter->SetParameter("pLambda", DLM_Fitter1::p_b, 1e-4, -2e-3,
                                   2e-3);
              fitter->FixParameter("pLambda", DLM_Fitter1::p_c, 0);
            } else if (BaseLine == 2) {
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
            
	    pointerFitResult = new TGraph(FitResult.GetN(), FitResult.GetX(), FitResult.GetY());

            pointerFitResult->SetLineWidth(2);
            pointerFitResult->SetLineColor(kRed);
            pointerFitResult->SetMarkerStyle(24);
            pointerFitResult->SetMarkerColor(kRed);
            pointerFitResult->SetMarkerSize(1);

            double Chi2 = 0;
            unsigned EffNumBins = 0;
            if (BaseLine == 0) {
              EffNumBins = -2;  // radius and normalization
            } else if (BaseLine == 0) {
              EffNumBins = -3;  // radius, normalization and slope
            } else if (BaseLine == 0) {
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
	    AvgRadius += fitter->GetParameter("pLambda", DLM_Fitter1::p_sor0);
	    
	    pa = fitter->GetParameter("pLambda", DLM_Fitter1::p_a);
	    pb = fitter->GetParameter("pLambda", DLM_Fitter1::p_b);
	    pc = fitter->GetParameter("pLambda", DLM_Fitter1::p_c);
	    ResultRadius  = fitter->GetParameter("pLambda", DLM_Fitter1::p_sor0);
	    ResultRadiusErr = fitter->GetParError("pLambda", DLM_Fitter1::p_sor0);
	    Stability = fitter->GetParameter("pLambda", DLM_Fitter1::p_sor1);
	    outChiSqNDF = Chi2 / EffNumBins;
	    outTree->Fill();
	    
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
  outTree->Write();   
  OutFile->Close();
  delete tidy;
}

int main(int argc, char *argv[]) {
  FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                  atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), argv[8], argv[9],
                  argv[10]);
  return 0;
}


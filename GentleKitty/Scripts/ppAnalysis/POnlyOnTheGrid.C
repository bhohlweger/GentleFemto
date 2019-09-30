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

void FitPPVariations(const unsigned& NumIter, int imTBin, int system, int source,
                     unsigned int iAngDist, int iRange, TString InputFile,
                     TString HistoName, TString OutputDir) {
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
  std::vector<double> ProSecondary = {0.79, 0.8,  0.81,  0.82,  0.83, 0.84, 0.87};
  std::vector<double> LamPurity = {0.92, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95 };
  std::vector<double> LamSecondary = { 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78};
  std::vector<float> mTValues = { 1.10765, 1.16831, 1.22838, 1.31558, 1.46282, 1.68716, 2.21156 }; 
  
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
    // PrimProton = 0.873; 
    // SecLamProton = 0.089;  //Fraction of Lambdas

    // PurityLambda = 0.961;
    // PrimLambdaAndSigma = 0.785;  //fraction of primary Lambdas + Sigma 0
    // SecLambda = 0.215;  //fraction of weak decay Lambdas

    PrimProton = ProSecondary[imTBin]; 
    SecLamProton = 0.7*(1-(double)PrimProton);  //Fraction of Lambdas

    PurityLambda = LamPurity[imTBin];
    PrimLambdaAndSigma = LamSecondary[imTBin];  //fraction of primary Lambdas + Sigma 0
    SecLambda = 1-PrimLambdaAndSigma;  //fraction of weak decay Lambdas

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
  CATSLambdaParam ppLam0(Proton[0], Proton[0], true);
  CATSLambdaParam ppLam1(Proton[1], Proton[1], true);
  CATSLambdaParam ppLam2(Proton[2], Proton[2], true);
  CATSLambdaParam ppLam3(Proton[3], Proton[3], true);
  CATSLambdaParam ppLam4(Proton[4], Proton[4], true);
  CATSLambdaParam ppLam5(Proton[5], Proton[5], true);
  CATSLambdaParam ppLam6(Proton[6], Proton[6], true);
  CATSLambdaParam ppLam7(Proton[7], Proton[7], true);
  CATSLambdaParam ppLam8(Proton[8], Proton[8], true);
  
  const std::vector<double> lam_pp = {
    ppLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
    ppLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Primary),
  };
  const std::vector<double> lam_pp_pL = {
    ppLam0.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0), 
    ppLam3.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam6.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0),
    ppLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0)
  };
  const std::vector<double> lam_pp_fake = {
    (ppLam0.GetLambdaParam(CATSLambdaParam::Primary,CATSLambdaParam::Fake) +
     ppLam0.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam0.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),    
    (ppLam1.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam1.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam1.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam2.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam2.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam2.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam3.GetLambdaParam(CATSLambdaParam::Primary,CATSLambdaParam::Fake) +
     ppLam3.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam3.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam4.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam4.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam4.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam5.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam5.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam5.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam6.GetLambdaParam(CATSLambdaParam::Primary,CATSLambdaParam::Fake) +
     ppLam6.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam6.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam7.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam7.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam7.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake)),
    (ppLam8.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake) +
     ppLam8.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary) +
     ppLam8.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake))
  };
  CATSLambdaParam pLLam(Proton[1], Lambda[1]);
  const float lam_pL = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
					    CATSLambdaParam::Primary);
  const float lam_pL_pS0 = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
						CATSLambdaParam::FeedDown, 0,
						0);
  const float lam_pL_pXm = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
						CATSLambdaParam::FeedDown, 0,
						1);
  const float lam_pL_fake = pLLam.GetLambdaParam(CATSLambdaParam::Primary,
						 CATSLambdaParam::Fake)
    + pLLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
    + pLLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);
  CATSLambdaParam pXiLam(Proton[1], Xi[1]);
  const float lam_pXim = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
					       CATSLambdaParam::Primary);
  const float lam_pXim_pXim1530 = pXiLam.GetLambdaParam(
							CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
  const float lam_pXim_fake = pXiLam.GetLambdaParam(CATSLambdaParam::Primary,
						    CATSLambdaParam::Fake)
    + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Primary)
    + pXiLam.GetLambdaParam(CATSLambdaParam::Fake, CATSLambdaParam::Fake);

  // for (int vFrac_pp_pL = 0; vFrac_pp_pL < 9; vFrac_pp_pL++) {
  //   std::cout << "lam_pp: " << lam_pp.at(vFrac_pp_pL) << " lam_pp_pL: "
  //             << lam_pp_pL.at(vFrac_pp_pL) << " lam_pp_fake: "
  //             << lam_pp_fake.at(vFrac_pp_pL) << std::endl;
  // }
  // std::cout << " lam_pL: " << lam_pL << "lam_pL_pXm: " << lam_pL_pXm
  // 	    << " lam_pL_fake: " << lam_pL_fake << std::endl;
  // std::cout << "lam_pXim: " << lam_pXim << " lam_pXim_pXim1530: "
  //           << lam_pXim_pXim1530 << " lam_pXim_fake:" << lam_pXim_fake
  //           << std::endl;
   
  std::vector<float> pSigma0Radii = { 1.55, 1.473, 1.421, 1.368, 1.295, 1.220, 1.124 };
  const float pFeeddownRadius = pSigma0Radii[imTBin];
  std::cout << "===========================\n";
  std::cout << "==pFeeddownRadius: " << pFeeddownRadius << "fm ==\n";
  std::cout << "===========================\n";

  //Link the output 
  TFile* OutFile = new TFile(TString::Format("%s/OutFileVarpp_%u.root", OutputDir.Data(), NumIter), "RECREATE");
  if (!OutFile) {
    return;
  }
  
  float total = TheSource == TidyCats::sLevy ? 162 : 243;
  int numIter = NumIter; 
  int uIter = 1; 
  int vFemReg;  //which femto region we use for pp (1 = default)
  float thismT = mTValues[imTBin]; 
  float FemtoFitMax = 0;
  float BaseLineMin = 0;
  float BaseLineMax = 0; 
  unsigned int vMod_pL = 1;  //which pL function to use: //0=exact NLO (at the moment temporary it is Usmani); 1=Ledni NLO; 2=Ledni LO; 3=ESC08
  unsigned int BaseLine = 0;
  float pa; 
  float pb;  
  float pc; 
  int vFrac_pp_pL;  //fraction of protons coming from Lambda variation (1 = default)
  float lmb_pp; 
  float lmb_ppL;
  float lmb_pL = lam_pL;
  float lmb_pLS0 = lam_pL_pS0;
  float lmb_pLXim = lam_pL_pXm; 
  float ResultRadius;
  float ResultRadiusErr;
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
  outTree->Branch("Stab",&Stability,"Stab/F");
  outTree->Branch("chiSqNDF",&outChiSqNDF,"chiSqNDF/F");
  outTree->Branch("CorrHist","TH1F",&StoreHist,sizeof(TH1F));
  outTree->Branch("FitResult","TGraph",&pointerFitResult,sizeof(TGraph));
  
  TidyCats* tidy = new TidyCats();


  //Setup the source stuff

  
  CATS AB_pp;
  tidy->GetCatsProtonProton(&AB_pp, NumMomBins, kMin, kMax, TheSource);
  if (TheSource == TidyCats::sResonance) {
    const double massProton =
      TDatabasePDG::Instance()->GetParticle(2212)->Mass() * 1000;
    const double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass()
      * 1000;
    DLM_CleverMcLevyReso* source = tidy->GetSourceProtonProton();
    if (iAngDist == 0) {
      source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion);
      source->SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion);
      std::cout << "Emission source is boring and back to back \n"; 
    } else if (iAngDist == 1) {
      //backwards > 90 Degree case
      source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion,
			false, false, DLM_CleverMcLevyReso::rdtRandomBackwards);
      source->SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion,
			false, false, DLM_CleverMcLevyReso::rdtRandomBackwards);
      std::cout
	<< "Sir, Emission will be fully randomized with a tad bit above 90 degrees, commencing countdown ... 3 ....\n";
    } else if (iAngDist > 1) {
      source->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion,
			false, false, DLM_CleverMcLevyReso::rdtRandom);
      source->SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, massProton, massPion,
			false, false, DLM_CleverMcLevyReso::rdtRandom);
      std::cout
	<< "Sir, Emission will be fully randomized, commencing countdown ... 3 ..2 .., but wait, how random is random ...? \n";
      if (iAngDist > 2) { 
	const char* PhiFile;
	const char* PhiHistoName; 
	int PhiRebin = 1;
	int PhiConversion = 1; 
	double RangeProtonMin, RangeProtonMax;
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
	} else if (iRange == 1) {
	  RangeProtonMin = 450/(double)PhiConversion;
	  RangeProtonMax = 550/(double)PhiConversion;
	} else if (iRange == 2) {
	  RangeProtonMin = 250/(double)PhiConversion;
	  RangeProtonMax = 650/(double)PhiConversion;
	} else {
	  std::cout << "Option iRange == " << iRange
		    << " is not a viable option \n";
	  return;
	}
	std::cout << "Using file: " << PhiFile << " in the range ProtonMin: "
		  << RangeProtonMin << " to Proton Max: " << RangeProtonMax
		  << std::endl;
	DLM_Histo<double>* HISTO = tidy->ConvertThetaAngleHisto(TString::Format("~/cernbox/WaveFunctions/ThetaDist/%s", PhiFile).Data(),TString::Format("%s",PhiHistoName).Data(), RangeProtonMin, RangeProtonMax, false, PhiRebin);
	source->SetUpResoEmission(0, 0, HISTO);
	source->SetUpResoEmission(1, 0, HISTO);
      }
    } else {
      std::cout << "Case for iAng == :" << iAngDist << " not implemented, exiting. \n";
      return;
    }
  }  
  AB_pp.KillTheCat();

  //Setup the feed down objects
  
  CATS AB_pXim;
  tidy->GetCatsProtonXiMinus(&AB_pXim, NumMomBins, kMin, kMax, FeeddownSource,
                             TidyCats::pHALQCD, 12);
  AB_pXim.SetAnaSource(0, pFeeddownRadius);
  AB_pXim.KillTheCat();

  CATS AB_pXim1530;
  AB_pXim1530.SetAnaSource(0, pFeeddownRadius);
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, NumMomBins, kMin, kMax,
                                 FeeddownSource);
  AB_pXim1530.KillTheCat();

  for (vMod_pL = TheSource == TidyCats::sLevy ? 1 : 0; vMod_pL < 3; ++vMod_pL) {
    TidyCats::pLPot PLpot;
    CATS AB_pL;
    if (vMod_pL == 1) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, FeeddownSource,
                                TidyCats::pUsmani);
      AB_pL.SetAnaSource(0, pFeeddownRadius);
      AB_pL.KillTheCat();
    } else if (vMod_pL == 2) {
      tidy->GetCatsProtonLambda(&AB_pL, NumMomBins, kMin, kMax, FeeddownSource,
                                TidyCats::pNLOWF);
      AB_pL.SetAnaSource(0, pFeeddownRadius);
      AB_pL.KillTheCat();
    }
    for (vFemReg = 0; vFemReg < 3; ++vFemReg) {
      FemtoFitMax = FemtoRegion[vFemReg];
      BaseLineMin = FemtoRegion[vFemReg];
      BaseLineMax = FemtoRegion[vFemReg];
      for (vFrac_pp_pL = 0; vFrac_pp_pL < 9; ++vFrac_pp_pL) {
	lmb_pp = lam_pp.at(vFrac_pp_pL);
	lmb_ppL = lam_pp_pL.at(vFrac_pp_pL);
        for (BaseLine = 0; BaseLine < 3; ++BaseLine) {
	  // Some computation here
          auto end = std::chrono::system_clock::now();

          std::chrono::duration<double> elapsed_seconds = end - start;
          std::cout
              << "\r Processing progress: "
              << TString::Format("%.1f %%", uIter / total * 100.f).Data()
              << " elapsed time: " << elapsed_seconds.count() / 60.
              << std::flush;

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
          Ck_pL->SetSourcePar(0, pFeeddownRadius);
          //this way you define a correlation function using Lednicky.
          //needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
          DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins, kMin, kMax,
                                          Lednicky_gauss_Sigma0);
          Ck_pSigma0->SetSourcePar(0, pFeeddownRadius);
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
            fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.4, 0.5, 2.5);
          }
          fitter->SetOutputDir(OutputDir.Data());

          fitter->SetSeparateBL(0, false);
          fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
          if (BaseLine == 1) {
            fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, -2e-3, 2e-3);
          } else if (BaseLine == 2) {
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
	  
	  pa = fitter->GetParameter("pp", DLM_Fitter1::p_a);
	  pb = fitter->GetParameter("pp", DLM_Fitter1::p_b);
	  pc = fitter->GetParameter("pp", DLM_Fitter1::p_c);
	  ResultRadius  = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
	  ResultRadiusErr = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
	  Stability = fitter->GetParameter("pp", DLM_Fitter1::p_sor1);
	  outChiSqNDF = Chi2 / EffNumBins;
	  outTree->Fill();
	  
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
  outTree->Write();
  OutFile->Close();
  delete tidy;
}

int main(int argc, char *argv[]) {
  FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                  atoi(argv[5]), atoi(argv[6]), argv[7], argv[8], argv[9]);
  return 0;
}


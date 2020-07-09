#include "ReadDreamFile.h"
#include "DreamPlot.h"
#include "CATSLambdaParam.h"
#include "TROOT.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

/// Ingredients for the fitter
/// 1. The p-phi correlation function that will be fitted
/// 2. The p-sideband correlation function from which we get the parametrization
/// 3. The momentum resolution matrix, from MC

/// 4. The lambda parameters (which need purity and primary fractions)
/// 5. The femtoscopic radius (for now we assume 1.3 fm)

/// 6. The baseline fit - taking care of the long-range correlations. Maybe we can put this to the sidebands
/// 7. The fit of the sideband correlation function

/// 8. The total fit of the correlation function, including the momentum resolution, lambda params, etc.
/// 9. The fitter itself, with a couple of start parameters etc.

///// Number of parameters for the sideband fit
//const int nSidebandPars = 0;

///// =====================================================================================
///// Fit for the sidebands - now a pol1
//auto sidebandFit = [ ] (double *x, double *p) {
//  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
//};

///// =====================================================================================
///// Function to cast the nice lambda from above to CATS...
//double sidebandFitCATS(const double &Momentum, const double *SourcePar,
//                       const double *PotPar) {
//  double *x = new double[1];
//  x[0] = Momentum;
//  double *p = const_cast<double*>(PotPar);
//  return sidebandFit(x, p);
//}

/// =====================================================================================
/// Get Correlationfunction

TH1F* GetCF2( const char* prefix,  const char* a) {
  const char* name=Form("/home/emma/FemtoPhiHM_BG2/CFOutput_pPhi_%s_%s.root",prefix, a);
  auto file = TFile::Open(name);
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;
}

TH1F* GetCF( const char* filename) {
  auto file = TFile::Open(filename);
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;
}

/// =====================================================================================
/// Get Minijets

TH1F* GetMJ( const char* prefix, int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/FemtoPhiHM_ROTMC/phitruth10001200/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_BG2/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
  TH1F* MJ = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return MJ;
}

/// =====================================================================================
/// Fit for the MJ - now a pol1
///
const int nMJPars = 6;

auto MJFit = [ ] (double *x, double *p) {
  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] * x[0] + p[4] * x[0] * x[0] * x[0] * x[0] + p[5] * x[0] * x[0] * x[0] * x[0] * x[0];
};

/// =====================================================================================
/// Function to cast the nice lambda from above to CATS...
double MJFitCATS(const double &Momentum, const double *SourcePar,
                       const double *PotPar) {
  double *x = new double[1];
  x[0] = Momentum;
  double *p = const_cast<double*>(PotPar);
  return MJFit(x, p);
}


/// =====================================================================================
/// Get Sidebands

TH1F* GetSB( const char* prefix, int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/FemtoPhiHM_BG2/240340/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_BG2/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
  TH1F* SB = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return SB;
}

TGraph* GetSBGraph( const char* prefix) {
  const char* name=Form("/home/emma/FemtoPhiHM_allSB/parasmer/%s",prefix);
//  const char* name=Form("/home/emma/FemtoPhiHM_BG2/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
  TGraph* SB = (TGraph*)file->FindObjectAny("cfSMEAR2");
  return SB;
}


TH1F* MJH =GetMJ("HMPhi",0);

TGraph* SB=GetSBGraph("SidebandParallelPOL2POL1.root");



/// Number of parameters for the sideband fit
const int nSidebandPars = 0;

/// =====================================================================================
/// Fit for the sidebands - now a pol1
auto sidebandFit = [ ] (double *x, double *p) {
    return SB->Eval(x[0]);
};

/// =====================================================================================
/// Function to cast the nice lambda from above to CATS...
double sidebandFitCATS(const double &Momentum, const double *SourcePar,
                       const double *PotPar) {
  double *x = new double[1];
  x[0] = Momentum;
  double *p = const_cast<double*>(PotPar);
  return sidebandFit(x, p);
}

/// Number of parameters for the real sideband fit (sb without MC truth minijet bg)
TGraph* GetSBrealGraph( const char* prefix) {
  const char* name=Form("/home/emma/FemtoPhiHM_allSB/parasmer/%s",prefix);
//  const char* name=Form("/home/emma/FemtoPhiHM_BG2/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
  TGraph* SB = (TGraph*)file->FindObjectAny("Ratio");
  return SB;
}

TGraph* SBreal=GetSBrealGraph("SidebandParallelRatio2.root");

auto sidebandrealFit = [ ] (double *x, double *p) {
    return SBreal->Eval(x[0]);
};

double sidebandrealFitCATS(const double &Momentum, const double *SourcePar,
                       const double *PotPar) {
  double *x = new double[1];
  x[0] = Momentum;
  double *p = const_cast<double*>(PotPar);
  return sidebandrealFit(x, p);
}




/// =====================================================================================
/// Get  Momentum smearing matrix

TH2F* GetMS( const char* prefix, int a, int b,  int c) {
  const char* p1=Form("%d",b);
  const char* p2=Form("%d",c);
  const char* addon=Form("%d",a);
 cout<<"M1"<<endl;
  const char* name=Form("/home/emma/FemtoPhiHM_AODMC2/AnalysisResults.root");
  auto file = TFile::Open(name);
  //file->ls();
  TDirectoryFile *dir=(TDirectoryFile*)(file->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *dirResults;
  dir->GetObject(Form("%sResults%s", prefix, addon), dirResults);
  TList* QAList = (TList*)dirResults->FindObject("PairQA");
 // QAList->ls();
  TList* resList = (TList*)QAList->FindObject(Form("QA_Particle%s_Particle%s", p1,p2));
  TH2F* MSgev = (TH2F*)resList->FindObject(Form("MomentumResolutionSE_Particle%s_Particle%s", p1,p2));

  cout<<"M3"<<endl;

  int nbinsx=MSgev->GetNbinsX();
  int nbinsy=MSgev->GetNbinsY();

  TH2F* MSmev = new TH2F ("MR MeV","Momentum Resolution in MeV",500,0,1000,500,0,1000);
  int binx;
  int biny;
  for (binx=0; binx<=500; binx++)
  {
      for (biny=0; biny<=500; biny++)
      {
          MSmev->SetBinContent(binx,biny,MSgev->GetBinContent(binx,biny));
      }
  }

  return MSmev;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");  // make ROOT shut up...
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
//  float p1=atof(argv[4]);
//  float p2=atof(argv[5]);
//  std::cout << "p1: " << p1 << " p2: " << p2 << std::endl;
  const char* filename = Form("%s/CFOutput_pPhi_%s_%s.root",argv[1],prefix,addon);
  DreamPlot::SetStyle();

  /// Set up the output file
  auto outputFile = new TFile("PhiOutput.root", "RECREATE");
  /// -----------------------------------------------------------------------------------
  /// 1. Get the correlation function

  /*
  ReadDreamFile* DreamFile = new ReadDreamFile(3, 3);
  DreamFile->SetAnalysisFile(filename, "Results", prefix, addon);
  DreamCF* CF_pPhi = new DreamCF();
  DreamPair* pPhi = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApPhi = new DreamPair("AntiPart", 0.24, 0.34);

  pPhi->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApPhi->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  pPhi->ShiftForEmpty(pPhi->GetPair());
  ApPhi->ShiftForEmpty(ApPhi->GetPair());
  pPhi->FixShift(pPhi->GetPairShiftedEmpty(0), ApPhi->GetPairShiftedEmpty(0),
                 ApPhi->GetFirstBin());
  ApPhi->FixShift(ApPhi->GetPairShiftedEmpty(0), pPhi->GetPairShiftedEmpty(0),
                  pPhi->GetFirstBin());
  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    pPhi->Rebin(pPhi->GetPairFixShifted(0), rebinVec[iReb]);
    pPhi->ReweightMixedEvent(pPhi->GetPairRebinned(iReb), 0.2, 0.9);
    ApPhi->Rebin(ApPhi->GetPairFixShifted(0), rebinVec[iReb]);
    ApPhi->ReweightMixedEvent(ApPhi->GetPairRebinned(iReb), 0.2, 0.9);
  }
  CF_pPhi->SetPairs(pPhi, ApPhi);
  CF_pPhi->GetCorrelations();
  std::vector<TH1F*> histVec = CF_pPhi->GetCorrelationFunctions();  // this vector contains all histograms in the order you find the CFOutput file
  TH1F* dataHist = histVec.at(13);  // TODO for now I randomly pick one of the histograms...

  outputFile->cd();
  dataHist->Write();
  */

  TH1F* dataHist=GetCF(filename);
  outputFile->cd();
  dataHist->Write();


  cout<<"P1"<<endl;
  /// -----------------------------------------------------------------------------------
  /// 2. Get the sideband correlation function
  // TODO For now I just copy the p-Phi CF
  // auto sidebandHist = (TH1F*) dataHist->Clone("sidebandFake");

 //  TH1F* sidebandHist =GetSB("HMPhi",2);

  TGraph* SBGr=GetSBGraph("SidebandParallelPOL2POL1.root");
  cout<<"P2"<<endl;

  TGraph* SBGrREAL=GetSBrealGraph("SidebandParallelRatio2.root");


  TH1F* MINIJETHist =GetMJ("HMPhi",0);
  cout<<"p222:"<<dataHist->GetNbinsX()<<endl;
  //TH1F* pPHI=nullptr;
  //divide minijtbg
//  TH1F *pPHI = new TH1F("h1", "pphi", dataHist->GetNbinsX(), 0, 3000);
//  for (int i=1; i<dataHist->GetNbinsX()+1;i++){
//      double B1=dataHist->GetBinContent(i);
//      double B2=MINIJETHist->GetBinContent(i);
//      cout<<"b1:"<<B1<<" B2:"<<B2<<endl;
//      double A=(B1/B2);
//      cout<<"A:"<<A<<endl;

//      pPHI->SetBinContent(i,A);
//  }


  dataHist->Divide(MINIJETHist);



  /// -----------------------------------------------------------------------------------
  /// 3. Momentum smearing matrix for the finite detector resolution
  /// TODO should come from MC

  /// Fake momentum smearing matrix with perfect resolution
//  TH2F *momentumSmearing = new TH2F("momSmear", "", 1000, 0, 3000, 1000, 0,
//                                    3000);
//  for (int i = 0; i < 3000; ++i) {
//    momentumSmearing->Fill(i, i);
//  }
  cout<<"P3"<<endl;

   TH2F *momentumSmearing02 = GetMS("HM", 0, 0, 2);

      cout<<"P4"<<endl;
   TH2F *momentumSmearing12 = GetMS("HM", 0, 1, 2);
   TH2F *momentumSmearing = new TH2F ("MS","Momsmearing",500,0,1000,500,0,1000);
   int binx;
   int biny;
   for (binx=0; binx<=500; binx++)
   {
       for (biny=0; biny<=500; biny++)
       {
           momentumSmearing->SetBinContent(binx,biny,(momentumSmearing02->GetBinContent(binx,biny)+momentumSmearing12->GetBinContent(binx,biny)));
       }
   }


   cout<<"P4"<<endl;

  /// -----------------------------------------------------------------------------------
  /// 4. Lambda parameters

  // Proton
  const double protonPurity = 0.9943;
  const double protonPrimary = 0.877;
  const double protonLambda = 0.089;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  // Phi
  const double phiPurity = 0.673;  // TODO
  const double phiPrimary = 1.;  // TODO
  const Particle phi(phiPurity, phiPrimary, { { 0 } });  // TODO for now we assume no secondary contributions

  const CATSLambdaParam lambdaParamDefault(phi, proton);
  const float PrimaryFraction = lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Primary);
  float SidebandFraction = lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::Primary, 0, 0);
  SidebandFraction += lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 0);
  SidebandFraction += lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 1);
  SidebandFraction += lambdaParamDefault.GetLambdaParam(CATSLambdaParam::Fake,
                                                        CATSLambdaParam::Fake);
  const float flatFraction = 1.f - PrimaryFraction - SidebandFraction;

  std::cout << "Lambda parameters\n";
  std::cout << " Primary  " << PrimaryFraction << "\n";
  std::cout << " Sideband " << SidebandFraction << "\n";
  std::cout << " Flat     " << flatFraction << "\n";

  /// -----------------------------------------------------------------------------------
  /// 5. Femtoscopic radius
  const double sourceRadius = 1.3;  // TODO we have to figure out what to put here

  /// -----------------------------------------------------------------------------------
  /// 6. Baseline fit - let's make sure that correlations present at larger k* are considered in the fit
//  TF1 *baselinePol1 = new TF1("baseline", "pol1", 600, 900);  // TODO check that the ranges make sense
//  baselinePol1->FixParameter(1,0);
//  baselinePol1->FixParameter(0,0.97);


//  TF1 *baselinePol1 = new TF1("baseline", "pol3", 600, 1300);  // TODO check that the ranges make sense
//  baselinePol1->FixParameter(1,0);
//  baselinePol1->FixParameter(0,0.97);

//  dataHist->Fit(baselinePol1, "FSNRMQ");

  //TF1 *baselinePol1 = new TF1("baseline", "pol4", 0, 1200);  // TODO check that the ranges make sense
//  baselinePol1->FixParameter(1,0);
//  baselinePol1->FixParameter(0,0.97);
  TF1 *baselinePol1 = new TF1("baseline", "pol2", 800, 1200);  // TODO check that the ranges make sense


 // MINIJETHist->Fit(baselinePol1, "FSNRMQ");
  dataHist->Fit(baselinePol1, "FSNRMQ");


//cout<<"MJ PAR 1:"<< baselinePol1->GetParameter(0)<<" P2:"<<baselinePol1->GetParameter(1)<<" P3:"<<baselinePol1->GetParameter(2)<<" P4:"<<baselinePol1->GetParameter(3)<<" P5:"<<baselinePol1->GetParameter(4)<<endl;
  /// -----------------------------------------------------------------------------------
  /// 7. Sideband fit
//  TF1 *sidebandFitFct = new TF1("sidebandFit", sidebandFit, 0, 600,
//                                nSidebandPars);  // TODO check that the ranges make sense
 // sidebandHist->Fit(sidebandFitFct, "NRMQ");


    TF1 *sidebandFitFct = new TF1("sidebandFit", sidebandFit, 0, 900,
                                  nSidebandPars);
    TF1 *sidebandREALFitFct = new TF1("sidebandREALFit", sidebandrealFit, 0, 1500,
                                  nSidebandPars);

  /// -----------------------------------------------------------------------------------
  /// 7. MJ fit
    TF1 *MJFitFct = new TF1("MJFit", MJFit, 0, 1500,
                                  nMJPars);  // TODO check that the ranges make sense
    MINIJETHist->Fit(MJFitFct, "NRMQ");
//---------------------------------------------------------------------------------------------
   ///Get SB without minijet cotribution!


  /// -----------------------------------------------------------------------------------
  /// 8. The full correlation function model
  const int binwidth = dataHist->GetBinWidth(1);
  const int NumMomBins = int(1000 / binwidth);
  //const double kMin = dataHist->GetBinCenter(1) - binwidth / 2.;
  const double kMin = 0;
  const double kMax = kMin + binwidth * NumMomBins;

  /// Lednicky model for one channel only, takes the inverse of the scattering length and the effective range as parameters
  /// These we want to obtain from the fit itself
  DLM_Ck *pPhiModel = new DLM_Ck(1, 2, NumMomBins, kMin, kMax,
                                 Lednicky_Singlet_InvScatLen);
  pPhiModel->SetPotPar(1, 0.5);
  pPhiModel->SetPotPar(0, 6);
  pPhiModel->SetSourcePar(0, sourceRadius);
  pPhiModel->Update();

  DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins, kMin, kMax,
                                   sidebandFitCATS);  // we make a correlation function out of this...
  for (unsigned i = 0; i < nSidebandPars; ++i) {
    Ck_SideBand->SetPotPar(i, sidebandFitFct->GetParameter(i));
  }
  Ck_SideBand->Update();

  DLM_Ck* Ck_SideBandREAL = new DLM_Ck(0, nSidebandPars, NumMomBins, kMin, kMax,
                                   sidebandrealFitCATS);  // we make a correlation function out of this...
  for (unsigned i = 0; i < nSidebandPars; ++i) {
    Ck_SideBand->SetPotPar(i, sidebandREALFitFct->GetParameter(i));
  }
  Ck_SideBand->Update();


  DLM_Ck* Ck_MJ = new DLM_Ck(0, nMJPars, NumMomBins, kMin, kMax,
                                   MJFitCATS);  // we make a correlation function out of this...
  for (unsigned i = 0; i < nMJPars; ++i) {
    Ck_MJ->SetPotPar(i, MJFitFct->GetParameter(i));
  }
  Ck_MJ->Update();




  DLM_CkDecomposition pPhiFullCF("pPhi", 2, *pPhiModel, momentumSmearing);

  //DLM_CkDecomposition sidebandFullCF("SideBand", 0, *Ck_SideBand, nullptr);

  DLM_CkDecomposition sidebandREALFullCF("SideBand", 0, *Ck_SideBandREAL, nullptr);

  //DLM_CkDecomposition MJFullCF("Minijets", 0, *Ck_MJ, nullptr); //CHECKKKKK if pars=0


//  pPhiFullCF.AddContribution(0, SidebandFraction, DLM_CkDecomposition::cFake,
//                             &sidebandFullCF);

  pPhiFullCF.AddContribution(0, SidebandFraction, DLM_CkDecomposition::cFake,
                             &sidebandREALFullCF);

  pPhiFullCF.AddContribution(1, flatFraction, DLM_CkDecomposition::cFeedDown);

//  pPhiFullCF.AddContribution(2, 1.f, DLM_CkDecomposition::cFeedDown,
//                             &MJFullCF);


  pPhiFullCF.Update();

  /// -----------------------------------------------------------------------------------
  /// 9. The fitter
  DLM_Fitter1* fitter = new DLM_Fitter1(1);
  fitter->SetSystem(0, *dataHist, 1, pPhiFullCF, kMin, 300, 800, 1200);
//  fitter->SetSystem(0, *dataHist, 1, pPhiFullCF, kMin, 300,400, 1000);
//  fitter->SetSystem(0, *dataHist, 1, pPhiFullCF, kMin, 150,350, 1000);

  fitter->SetSeparateBL(0, false);        // no simultaneous fit of the baseline
  fitter->FixParameter("pPhi", DLM_Fitter1::p_a, baselinePol1->GetParameter(0));
  fitter->FixParameter("pPhi", DLM_Fitter1::p_b, baselinePol1->GetParameter(1));
  fitter->FixParameter("pPhi", DLM_Fitter1::p_c, baselinePol1->GetParameter(2));
//  fitter->FixParameter("pPhi", DLM_Fitter1::p_3, baselinePol1->GetParameter(3));
//  fitter->FixParameter("pPhi", DLM_Fitter1::p_4, baselinePol1->GetParameter(4));

//  fitter->AddSameSource("SideBand", "pPhi", 1);
//  fitter->FixParameter("pPhi", DLM_Fitter1::p_c, 0);
//  fitter->FixParameter("pPhi", DLM_Fitter1::p_Cl, -1.);
//  fitter->FixParameter("pPhi", DLM_Fitter1::p_sor0, sourceRadius);

  // Here you can tune the start values for the scattering parameters
  //fitter->SetParameter("pPhi", DLM_Fitter1::p_pot0, p1);
  //  fitter->SetParameter("pPhi", DLM_Fitter1::p_pot1, p2);
//     fitter->SetParameter("pPhi", DLM_Fitter1::p_pot0, p1, -10., 10.);
//     fitter->SetParameter("pPhi", DLM_Fitter1::p_pot1, p2, 0, 20.);
//     fitter->SetParameter("pPhi", DLM_Fitter1::p_pot0, p1, p1-5, p1+5);
//     fitter->SetParameter("pPhi", DLM_Fitter1::p_pot1, p2, p2-5, p2+5);

  // In case you want to check what some specific correlation function looks like
//   fitter->FixParameter("pPhi", DLM_Fitter1::p_pot0, p1);
//   fitter->FixParameter("pPhi", DLM_Fitter1::p_pot1, p2);

  //std::cout <<"fit start \n";

  /// Run the fit
  fitter->GoBabyGo();

  /// -----------------------------------------------------------------------------------
  /// Post-processing
//  std::cout <<"fit finished \n";

  /// Get the parameters from the fit
  const double bl_a = fitter->GetParameter("pPhi", DLM_Fitter1::p_a);
  const double bl_a_err = fitter->GetParError("pPhi", DLM_Fitter1::p_a);
  const double bl_b = fitter->GetParameter("pPhi", DLM_Fitter1::p_b);
  const double bl_b_err = fitter->GetParError("pPhi", DLM_Fitter1::p_b);
  const double bl_c = fitter->GetParameter("pPhi", DLM_Fitter1::p_c);
  const double bl_c_err = fitter->GetParError("pPhi", DLM_Fitter1::p_c);
//  const double bl_3 = fitter->GetParameter("pPhi", DLM_Fitter1::p_3);
//  const double bl_3_err = fitter->GetParError("pPhi", DLM_Fitter1::p_3);
//  const double bl_4 = fitter->GetParameter("pPhi", DLM_Fitter1::p_4);
//  const double bl_4_err = fitter->GetParError("pPhi", DLM_Fitter1::p_4);


  const double invScatteringLength = fitter->GetParameter("pPhi",
                                                          DLM_Fitter1::p_pot0);
  const double invScatteringLengthErr = fitter->GetParError(
      "pPhi", DLM_Fitter1::p_pot0);
 // std::cout <<"fit parameter 1 finished \n";
  const double effectiveRange = fitter->GetParameter("pPhi",
                                                     DLM_Fitter1::p_pot1);
  const double effectiveRangeErr = fitter->GetParError("pPhi",
                                                       DLM_Fitter1::p_pot1);
  //std::cout <<"fit parameter 2 finished \n";






  std::cout << "++++++++++++++++++++++++++\n";
  std::cout << "Result of the fit \n";
  std::cout << "Inv. scattering length: " << invScatteringLength << " +/- "
            << invScatteringLengthErr << " fm^-1 \n";
  std::cout << "Effective range: " << effectiveRange << " +/- "
            << effectiveRangeErr << " fm \n";


//  std::cout << "Inv. scattering length: " << invScatteringLength << " +/- "
//            << invScatteringLengthErr << "  ";
//  std::cout << "Effective range: " << effectiveRange << " +/- "
//            << effectiveRangeErr << "   ";


  TGraph grFitResult;
  grFitResult.SetName("LednickyFit");
  grFitResult.SetTitle("LednickyFit");
  fitter->GetFitGraph(0, grFitResult);
  // Globales chi2:

  double Chi2=0.0;
  Chi2 = fitter->GetChi2();

  std::cout << "global chi2: " << Chi2<< "\n";
  std::cout << "Baseline p0: " << baselinePol1->GetParameter(0)<< "\n";
  std::cout << "Baseline p1: " << baselinePol1->GetParameter(1)<< "\n";
  std::cout << "Baseline p2: " << baselinePol1->GetParameter(2)<< "\n";
//  std::cout << "Baseline p3: " << baselinePol1->GetParameter(3)<< "\n";
//  std::cout << "Baseline p4: " << baselinePol1->GetParameter(4)<< "\n";

//  std::cout << "Sideband p0: " <<  MJFitFct->GetParameter(0)<< "\n";
//  std::cout << "Sideband p1: " <<  MJFitFct->GetParameter(1)<< "\n";
//  std::cout << "Sideband p2: " <<  MJFitFct->GetParameter(2)<< "\n";
//  std::cout << "Sideband p3: " <<  MJFitFct->GetParameter(3)<< "\n";
//  std::cout << "Sideband p4: " <<  MJFitFct->GetParameter(4)<< "\n";
//  std::cout << "Sideband p5: " <<  MJFitFct->GetParameter(5)<< "\n";


  Chi2=0.0;

  // Lokales chi2:
  int counter=0;

  for (unsigned uBin = 1; uBin < dataHist->GetNbinsX()-1; uBin++) {

    double mom = dataHist->GetBinCenter(uBin);

    double dataY;

    double dataErr;

    double theoryX;

    double theoryY;

    if (mom > 200) continue;

    grFitResult.GetPoint(uBin-1, theoryX, theoryY);

    if (mom != theoryX) {

      std::cout << mom << '\t' << theoryX << std::endl;

      printf("  PROBLEM!\n");

    }

    dataY = dataHist->GetBinContent(uBin);

    dataErr = dataHist->GetBinError(uBin);

    if (dataErr < 1e-5) {

      std::cout << dataErr << '\t' << "WARNING POINT NOT CONSIDERED \n";

                continue;

    }

    Chi2 += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
  //  std::cout <<uBin << "chi2: " << Chi2<< "\n";
    counter++;
  //  std::cout << "counter: " << counter << " mom:"<< mom<<"\n";


  }

  std::cout << "local chi2: " << Chi2<< "  ";
  std::cout << "chi2 over points: " << Chi2/(counter-2)<< "\n";
//  int aa= a*10;
//  int bb= b*10;
//  const char* l=Form("%d",aa);;
//  const char* m=Form("%d",bb);

 /// write to file
//  const char* txtname = Form("cout_%s_%s.txt"),l,m);
//  file.open ("txtname");
//  file << "++++++++++++++++++++++++++\n";
//  file << "Result of the fit \n";
//  file << "Inv. scattering length: " << invScatteringLength << " +/- "
//            << invScatteringLengthErr << " fm^-1 \n";
//  file << "Effective range: " << effectiveRange << " +/- "
//            << effectiveRangeErr << " fm \n";
//  file <<"local Chi2: "<<Chi2<<endl;
//  file.close();

//  TF1 *SBScaled= new TF1("SBscaled", "1+(([0]+[1]*x+[2]*x*x-1)*[3])",0,900);
//  SBScaled->SetParameter(0,sidebandFitFct->GetParameter(0));
//  SBScaled->SetParameter(1,sidebandFitFct->GetParameter(1));
//  SBScaled->SetParameter(2,sidebandFitFct->GetParameter(2));
//  SBScaled->SetParameter(3,SidebandFraction);

  int nP=1000;

  TGraph *SBScaled=new TGraph(nP);
   for (int i=0; i<nP; i++){
         SBScaled->SetPoint(i, SBGr->GetX()[i], 1+((SBGr->GetY()[i]-1)*SidebandFraction));

}

   TGraph *SBrealScaled=new TGraph(nP);
    for (int i=0; i<nP; i++){
          SBrealScaled->SetPoint(i, SBGrREAL->GetX()[i], 1+((SBGrREAL->GetY()[i]-1)*SidebandFraction));

 }

  std::fstream output;

  output.open ("output.dat", std::fstream::in | std::fstream::out | std::fstream::app);

  output << Form("%.2f ",invScatteringLength) << Form("%.2f ",effectiveRange)<< Form("%.3f ",Chi2) << Form("%.3f",Chi2/(counter-2)) << "\n";

  output.close();

  TGraph grGenuine;
  grGenuine.SetName("GenuineCF");
  TGraph grFeedCF;
  grFeedCF.SetName("GenuineCFSmearedLambda");
  TGraph grSidebandCF;
  grSidebandCF.SetName("SidebandCF");
  for (int i = 0; i < NumMomBins; ++i) {
    const float mom = pPhiModel->GetBinCenter(0, i);
    grGenuine.SetPoint(i, mom, pPhiFullCF.EvalMain(mom));  // genuine p-Phi CF with the parameters obtained in the fit
    grFeedCF.SetPoint(i, mom, pPhiFullCF.EvalMainFeed(mom));  // same as above, scaled by lambda params and momentum smearing
 //   grSidebandCF.SetPoint(i, mom, sidebandFullCF.EvalMain(mom));
  }

  /// Write all the relevant stuff out
  outputFile->cd();
  dataHist->Write("dicide");
 // pPHI->Write("pphiw/oMJ");
  baselinePol1->Write("Baseline");
  sidebandFitFct->Write("SB+MJ");
  sidebandREALFitFct->Write("SB only");
  SBGr->Write("SB+MJGraph");
  SBGrREAL->Write("SBonlyGraph");
  MINIJETHist->Write("Minijet(MCTRUTH)");
  //MJFitFct->Write();
  momentumSmearing->Write();
  grFitResult.Write();
  grGenuine.Write();
  grFeedCF.Write();
  //grSidebandCF.Write();
  SBScaled->Write("SB+MJscaled");
  SBrealScaled->Write("SBonlyscaled");
  outputFile->Close();

  return 0;
}

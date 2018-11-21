#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "CATSInputSigma0.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include "CATSLambdaParam.h"
#include "SideBandFit.h"
#include "TMinuit.h"
#include "TMath.h"

const int nSidebandPars = 5;

// Fit for the sidebands
auto sidebandFit = [ ] (double *x, double *p) {
  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + std::exp(p[3] + p[4] * x[0]);
};

// Function to cast the nice lambda from above to CATS...
double sidebandFitCATS(const double &Momentum, const double *SourcePar,
                       const double *PotPar) {
  double *x = new double[1];
  x[0] = Momentum;
  double *p = const_cast<double*>(PotPar);
  return sidebandFit(x, p);
}

void FitSigma0(const unsigned& NumIter, TString InputDir, TString appendix,
               TString ppFile, TString OutputDir) {
  bool fastPlot = true;
  TRandom3 rangen(0);
  const int binwidth = 10;
  const unsigned NumMomBins_pSigma = 100;

  double kMin_pSigma = 4.;
  const double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;

  std::cout << "kMin_pSigma: " << kMin_pSigma << std::endl;
  std::cout << "kMax_pSigma: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pSigma: " << NumMomBins_pSigma << std::endl;

  // TO BE DONE - Get the radius
  const double sourceSize = 1.36;

  double FemtoRegion_pSigma[3][2];
  FemtoRegion_pSigma[0][0] = kMin_pSigma;
  FemtoRegion_pSigma[0][1] = 650;
  FemtoRegion_pSigma[1][0] = kMin_pSigma;
  FemtoRegion_pSigma[1][1] = 1000;
  FemtoRegion_pSigma[2][0] = kMin_pSigma;
  FemtoRegion_pSigma[2][1] = 1000;

  // Lambda parameters
  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;

  const double sigmaPurity = 0.199;
  const double sigmaPrimary = 1.;

  const Particle proton(
      protonPurity,
      protonPrimary,
      { { protonLambda, (1. - protonPrimary)
          * (1. - protonLambda / (1. - protonPrimary)) } });

  const Particle sigma0(0.2, sigmaPrimary, { { 0 } });

  const CATSLambdaParam lambdaParams(proton, sigma0);
  lambdaParams.PrintLambdaParams();

  const float primFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Fake);

  // CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendix.Data());
  CATSinput->ObtainCFs(10, 350, 550);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHist = CATSinput->GetCF("pSigma0", dataHistName.Data());
  if (!dataHist) {
    std::cerr << "ERROR pSigma0 fitter: p-Sigma0 histogram missing\n";
    return;
  }
  auto sidebandHistUp = CATSinput->GetCF("pSigmaSBUp",
                                         "hCk_ReweightedpSigmaSBUpMeV_0");

  auto sidebandHistLow = CATSinput->GetCF("pSigmaSBLow",
                                          "hCk_ReweightedpSigmaSBLowMeV_0");
  if (!sidebandHistLow || !sidebandHistUp) {
    std::cerr << "ERROR pSigma0 fitter: p-pSigmaSB histogram missing\n";
    return;
  }

  // Fit the sideband
  SideBandFit* side = new SideBandFit();
  auto SBmerge = side->AddCF(sidebandHistUp, sidebandHistLow, "SBmerge");
  auto sideband = new TF1("sideband", sidebandFit, 0, 1000, nSidebandPars);
  SBmerge->Fit(sideband, "R");

  // Prefit for the baseline
  float p_a_prefit[4];
  float p_b_prefit[4];
  auto Prefit = (TH1F*) dataHist->Clone(Form("%s_prefit", dataHist->GetName()));
  auto funct_0 = new TF1("myPol0", "pol0", 250, 650);
  Prefit->Fit(funct_0, "FSNRMQ");
  p_a_prefit[0] = funct_0->GetParameter(0);
  p_b_prefit[0] = 0;

  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 650);
  Prefit->Fit(funct_1, "FSNRMQ");
  gMinuit->SetErrorDef(4);  //note 4 and not 2!
  p_a_prefit[1] = funct_1->GetParameter(0);
  p_b_prefit[1] = funct_1->GetParameter(1);

  auto grPrefitContour = (TGraph*) gMinuit->Contour(40, 0, 1);

  p_a_prefit[2] = TMath::MinElement(grPrefitContour->GetN(),
                                    grPrefitContour->GetX());
  p_b_prefit[2] = grPrefitContour->Eval(p_a_prefit[2]);
  p_a_prefit[3] = TMath::MaxElement(grPrefitContour->GetN(),
                                    grPrefitContour->GetX());
  p_b_prefit[3] = grPrefitContour->Eval(p_a_prefit[3]);

  std::cout << "Result of the prefit to constrain the baseline\n";
  for (int i = 0; i < 4; ++i) {
    std::cout << i << " a: " << p_a_prefit[i] << " b: " << p_b_prefit[i]
              << "\n";
  }
  delete Prefit;

  // Correlation function
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile
  DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins_pSigma, kMin_pSigma,
                                  kMax_pSigma, Lednicky_gauss_Sigma0);
  Ck_pSigma0->SetSourcePar(0, sourceSize);
  Ck_pSigma0->Update();

  // TO BE FIXED - MOMENTUM RESOLUTION MATRIX!
  DLM_CkDecomposition CkDec_pSigma0("pSigma0", 1, *Ck_pSigma0,
                                    CATSinput->GetSigmaFile(0));

  // Sidebands
  DLM_Ck* Ck_SideBand = new DLM_Ck(0, 5, NumMomBins_pSigma, kMin_pSigma,
                                   kMax_pSigma, sidebandFitCATS);
  std::cout << "Result of the sideband fit\n";
  for (unsigned i = 0; i < sideband->GetNumberFreeParameters(); ++i) {
    std::cout << i << " " << sideband->GetParameter(i) << std::endl;
    Ck_SideBand->SetPotPar(i, sideband->GetParameter(i));
  }
  Ck_SideBand->Update();

  DLM_CkDecomposition CkDec_SideBand("pSigma0SideBand", 0, *Ck_SideBand,
                                     nullptr);
  CkDec_pSigma0.AddContribution(
      0, lambdaParams.GetLambdaParam(CATSLambdaParam::Fake),
      DLM_CkDecomposition::cFake, &CkDec_SideBand);
  CkDec_pSigma0.Update();

  // Fitter
  DLM_Fitter1* fitter = new DLM_Fitter1(1);
  fitter->SetSystem(0, *dataHist, 1, CkDec_pSigma0, FemtoRegion_pSigma[0][0],
                    FemtoRegion_pSigma[0][1], 1000, 1000);
  fitter->SetSeparateBL(0, false);  //Simultaneous BL
//  if (HaveWeABaseLineSlope) {
//    fitter->FixParameter("pSigma0", DLM_Fitter1::p_a, p_a_prefit[BaselineSlope]);
//    fitter->FixParameter("pSigma0", DLM_Fitter1::p_b, p_b_prefit[BaselineSlope]);
//    std::cout << "Fitting ranges for BL set \n";
//  } else {
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_a, p_a_prefit[0]);
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_b, p_b_prefit[0]);
//  }
  fitter->AddSameSource("pSigma0SideBand", "pSigma0", 1);
  //Fit BL & Normalization
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_c, 0);
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_Cl, -1.);
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0, sourceSize);

  fitter->GoBabyGo();

  // Get the parameters from the fit
  const double bl_a = fitter->GetParameter("pSigma0", DLM_Fitter1::p_a);
  const double bl_a_err = fitter->GetParError("pSigma0", DLM_Fitter1::p_a);
  const double bl_b = fitter->GetParameter("pSigma0", DLM_Fitter1::p_b);
  const double bl_b_err = fitter->GetParError("pSigma0", DLM_Fitter1::p_b);
  const double Cl = fitter->GetParameter("pSigma0", DLM_Fitter1::p_c);
  const double chi2 = fitter->GetChi2Ndf();
  const double pval = fitter->GetPval();

  std::cout << "=============\n";
  std::cout << "Fitter output\n";
  std::cout << "BL a " << bl_a << " " << bl_a_err << "\n";
  std::cout << "BL b " << bl_b << " " << bl_b_err << "\n";
  std::cout << "Cl   " << Cl << "\n";
  std::cout << "Chi2 " << chi2 << "\n";
  std::cout << "=============\n";

  TGraph FitResult_pSigma0;
  FitResult_pSigma0.SetName(TString::Format("pSigma0Graph"));
  fitter->GetFitGraph(0, FitResult_pSigma0);

// Write to output
  TGraph grCFSigmaRaw;
  grCFSigmaRaw.SetName("Sigma0Raw");
  TGraph grCFSigmaMain;
  grCFSigmaMain.SetName("Sigma0Main");
  TGraph grCFSigmaFeed;
  grCFSigmaFeed.SetName("Sigma0Feed");
  TGraph grCFSigmaSideband;
  grCFSigmaSideband.SetName("Sigma0Sideband");
  for (unsigned int i = 0; i < Ck_pSigma0->GetNbins(); ++i) {
    grCFSigmaRaw.SetPoint(i, Ck_pSigma0->GetBinCenter(i),
                          CkDec_pSigma0.EvalCk(Ck_pSigma0->GetBinCenter(i)));
    grCFSigmaMain.SetPoint(
        i,
        Ck_pSigma0->GetBinCenter(i),
        ((CkDec_pSigma0.EvalMain(Ck_pSigma0->GetBinCenter(i)) - 1.)
            * lambdaParams.GetLambdaParam(CATSLambdaParam::Primary)) + 1);
    grCFSigmaFeed.SetPoint(
        i, Ck_pSigma0->GetBinCenter(i),
        CkDec_pSigma0.EvalMainFeed(Ck_pSigma0->GetBinCenter(i)));
    grCFSigmaSideband.SetPoint(
        i,
        Ck_pSigma0->GetBinCenter(i),
        ((Ck_SideBand->Eval(Ck_pSigma0->GetBinCenter(i)) - 1.)
            * lambdaParams.GetLambdaParam(CATSLambdaParam::Fake)) + 1);
  }

  TString outfilename = OutputDir;
  outfilename += "/CF_Sigma0.root";
  auto out = new TFile(outfilename, "RECREATE");
  grPrefitContour->Write("fitContour");
  grCFSigmaRaw.Write();
  grCFSigmaMain.Write();
  grCFSigmaFeed.Write();
  grCFSigmaSideband.Write();
  sidebandHistLow->Write("SidebandLow");
  sidebandHistUp->Write("SidebandUp");
  SBmerge->Write("SidebandMerged");
  dataHist->Write("CF");
  FitResult_pSigma0.Write("Fit");
  out->Close();
}

int main(int argc, char *argv[]) {
  FitSigma0(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5]);
  return 0;
}

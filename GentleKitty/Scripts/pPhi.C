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

/// Number of parameters for the sideband fit
const int nSidebandPars = 3;

/// =====================================================================================
/// Fit for the sidebands - now a pol1
auto sidebandFit = [ ] (double *x, double *p) {
  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
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

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");  // make ROOT shut up...
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  DreamPlot::SetStyle();

  /// Set up the output file
  auto outputFile = new TFile("PhiOutput.root", "RECREATE");

  /// -----------------------------------------------------------------------------------
  /// 1. Get the correlation function
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

  /// -----------------------------------------------------------------------------------
  /// 2. Get the sideband correlation function
  // TODO For now I just copy the p-Phi CF
  auto sidebandHist = (TH1F*) dataHist->Clone("sidebandFake");

  /// -----------------------------------------------------------------------------------
  /// 3. Momentum smearing matrix for the finite detector resolution
  /// TODO should come from MC

  /// Fake momentum smearing matrix with perfect resolution
  TH2F *momentumSmearing = new TH2F("momSmear", "", 1000, 0, 3000, 1000, 0,
                                    3000);
  for (int i = 0; i < 3000; ++i) {
    momentumSmearing->Fill(i, i);
  }

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
  const double phiPurity = 0.7;  // TODO
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
  TF1 *baselinePol1 = new TF1("baseline", "pol1", 250, 600);  // TODO check that the ranges make sense
  dataHist->Fit(baselinePol1, "FSNRMQ");

  /// -----------------------------------------------------------------------------------
  /// 7. Sideband fit
  TF1 *sidebandFitFct = new TF1("sidebandFit", sidebandFit, 0, 600,
                                nSidebandPars);  // TODO check that the ranges make sense
  sidebandHist->Fit(sidebandFitFct, "NRMQ");

  /// -----------------------------------------------------------------------------------
  /// 8. The full correlation function model
  const int binwidth = dataHist->GetBinWidth(1);
  const int NumMomBins = int(600 / binwidth);
  const double kMin = dataHist->GetBinCenter(1) - binwidth / 2.;
  const double kMax = kMin + binwidth * NumMomBins;

  /// Lednicky model for one channel only, takes the inverse of the scattering length and the effective range as parameters
  /// These we want to obtain from the fit itself
  DLM_Ck *pPhiModel = new DLM_Ck(1, 2, NumMomBins, kMin, kMax,
                                 Lednicky_Singlet_InvScatLen);
  pPhiModel->SetSourcePar(0, sourceRadius);
  pPhiModel->Update();

  DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins, kMin, kMax,
                                   sidebandFitCATS);  // we make a correlation function out of this...
  for (unsigned i = 0; i < nSidebandPars; ++i) {
    Ck_SideBand->SetPotPar(i, sidebandFitFct->GetParameter(i));
  }
  Ck_SideBand->Update();

  DLM_CkDecomposition pPhiFullCF("pPhi", 2, *pPhiModel, momentumSmearing);

  DLM_CkDecomposition sidebandFullCF("SideBand", 0, *Ck_SideBand, nullptr);

  pPhiFullCF.AddContribution(0, SidebandFraction, DLM_CkDecomposition::cFake,
                             &sidebandFullCF);
  pPhiFullCF.AddContribution(1, flatFraction, DLM_CkDecomposition::cFeedDown);
  pPhiFullCF.Update();

  /// -----------------------------------------------------------------------------------
  /// 9. The fitter
  DLM_Fitter1* fitter = new DLM_Fitter1(1);
  fitter->SetSystem(0, *dataHist, 1, pPhiFullCF, kMin, 600, 900, 900);
  fitter->SetSeparateBL(0, false);        // no simultaneous fit of the baseline
  fitter->FixParameter("pPhi", DLM_Fitter1::p_a, baselinePol1->GetParameter(0));
  fitter->FixParameter("pPhi", DLM_Fitter1::p_b, baselinePol1->GetParameter(1));
  fitter->AddSameSource("SideBand", "pPhi", 1);
  fitter->FixParameter("pPhi", DLM_Fitter1::p_c, 0);
  fitter->FixParameter("pPhi", DLM_Fitter1::p_Cl, -1.);
  fitter->FixParameter("pPhi", DLM_Fitter1::p_sor0, sourceRadius);

  /// Here you can tune the start values for the scattering parameters
  fitter->SetParameter("pPhi", DLM_Fitter1::p_pot0, 2);
  fitter->SetParameter("pPhi", DLM_Fitter1::p_pot1, 0.5);

  // In case you want to check what some specific correlation function looks like
  // fitter->FixParameter("pPhi", DLM_Fitter1::p_pot0, 2);
  // fitter->FixParameter("pPhi", DLM_Fitter1::p_pot1, 0.5);

  /// Run the fit
  fitter->GoBabyGo();

  /// -----------------------------------------------------------------------------------
  /// Post-processing

  /// Get the parameters from the fit
  const double bl_a = fitter->GetParameter("pPhi", DLM_Fitter1::p_a);
  const double bl_a_err = fitter->GetParError("pPhi", DLM_Fitter1::p_a);
  const double bl_b = fitter->GetParameter("pPhi", DLM_Fitter1::p_b);
  const double bl_b_err = fitter->GetParError("pPhi", DLM_Fitter1::p_b);

  const double invScatteringLength = fitter->GetParameter("pPhi",
                                                          DLM_Fitter1::p_pot0);
  const double invScatteringLengthErr = fitter->GetParError(
      "pPhi", DLM_Fitter1::p_pot0);
  const double effectiveRange = fitter->GetParameter("pPhi",
                                                     DLM_Fitter1::p_pot1);
  const double effectiveRangeErr = fitter->GetParError("pPhi",
                                                       DLM_Fitter1::p_pot1);

  std::cout << "++++++++++++++++++++++++++\n";
  std::cout << "Result of the fit \n";
  std::cout << "Inv. scattering length: " << invScatteringLength << " +/- "
            << invScatteringLengthErr << " fm^-1 \n";
  std::cout << "Effective range: " << effectiveRange << " +/- "
            << effectiveRangeErr << " fm \n";

  TGraph grFitResult;
  grFitResult.SetName("LednickyFit");
  grFitResult.SetTitle("LednickyFit");
  fitter->GetFitGraph(0, grFitResult);

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
    grSidebandCF.SetPoint(i, mom, sidebandFullCF.EvalMain(mom));
  }

  /// Write all the relevant stuff out
  outputFile->cd();
  baselinePol1->Write();
  sidebandFitFct->Write();
  momentumSmearing->Write();
  grFitResult.Write();
  grGenuine.Write();
  grFeedCF.Write();
  grSidebandCF.Write();
  outputFile->Close();
  return 0;
}

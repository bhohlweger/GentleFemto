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
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include "CATSLambdaParam.h"
#include "SideBandFit.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"

/// Number of parameters for the sideband fit
const int nSidebandPars = 6;

/// =====================================================================================
/// Fit for the sidebands
auto sidebandFit =
    [ ] (double *x, double *p) {
      return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] *x[0] + std::exp(p[4] + p[5] * x[0]);
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
void FitSigma0exclusion(const unsigned& NumIter, TString InputDir,
                        TString appendix, TString ppFile, TString OutputDir,
                        const float rebin) {
  DreamPlot::SetStyle();
  bool fastPlot = (NumIter == 0) ? true : false;
  TRandom3 rangen(0);
  int iterID = 0;
  bool useBaselineSlope = true;
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  TString exclusionfilename = TString::Format("%s/Exclusion_pSigma0_%i.root",
                                              OutputDir.Data(), NumIter);
  auto exclusionFile = new TFile(exclusionfilename, "RECREATE");

  int NumBins_f0Inv = 1000;
  double Min_f0Inv = -5;
  double Max_f0Inv = 5;
  int NumBins_d0 = 1000;
  double Min_d0 = 0.;
  double Max_d0 = 15.;

  NumBins_f0Inv /= rebin;
  NumBins_d0 /= rebin;

  auto hChi2Local = new TH2F(
      "hChi2Local",
      "; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #chi^{2} _{local}",
      NumBins_f0Inv, Min_f0Inv, Max_f0Inv, NumBins_d0, Min_d0, Max_d0);
  auto hChi2Global = new TH2F(
      "hChi2Global",
      "; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #chi^{2} _{global}",
      NumBins_f0Inv, Min_f0Inv, Max_f0Inv, NumBins_d0, Min_d0, Max_d0);
  auto hDeltaChi2Plot = new TH2F(
      "hDeltaChi2Plot",
      "; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); #Delta#chi^{2}",
      NumBins_f0Inv, Min_f0Inv, Max_f0Inv, NumBins_d0, Min_d0, Max_d0);
  auto hnSigmaPlot = new TH2F(
      "hnSigmaPlot", "; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm); n_{#sigma}",
      NumBins_f0Inv, Min_f0Inv, Max_f0Inv, NumBins_d0, Min_d0, Max_d0);
  auto hNegCk = new TH2F("hNegCk", "; 1/#it{f}_{0} (fm^{-1}); #it{d}_{0} (fm)",
                         NumBins_f0Inv, Min_f0Inv, Max_f0Inv, NumBins_d0,
                         Min_d0, Max_d0);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendix.Data());
  CATSinput->ObtainCFs(10, 340, 440);
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

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.
  const int binwidth = 40;
  const unsigned NumMomBins_pSigma = 20;
  double kMin_pSigma = dataHist->GetBinCenter(1)
      - dataHist->GetBinWidth(1) / 2.f;
  const double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;

  std::cout << "kMin_pSigma: " << kMin_pSigma << std::endl;
  std::cout << "kMax_pSigma: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pSigma: " << NumMomBins_pSigma << std::endl;

  // TODO TO BE DONE - Get the radius from the pp file
  double sourceSize[3];
  sourceSize[0] = 1.36;
  sourceSize[1] = 1.36 * 0.8;
  sourceSize[2] = 1.36 * 1.2;

  double FemtoRegion_pSigma[3][2];
  FemtoRegion_pSigma[0][0] = kMin_pSigma;
  FemtoRegion_pSigma[0][1] = 500;
  FemtoRegion_pSigma[1][0] = kMin_pSigma;
  FemtoRegion_pSigma[1][1] = 550;
  FemtoRegion_pSigma[2][0] = kMin_pSigma;
  FemtoRegion_pSigma[2][1] = 600;

  /// Lambda parameters
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

  const Particle sigma0(sigmaPurity, sigmaPrimary, { { 0 } });

  const CATSLambdaParam lambdaParams(proton, sigma0);
  lambdaParams.PrintLambdaParams();

  const float primFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Fake);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Fit the sideband
  SideBandFit* side = new SideBandFit();
  auto SBmerge = side->AddCF(sidebandHistUp, sidebandHistLow, "SBmerge");
  auto sideband = new TF1("sideband", sidebandFit, kMin_pSigma, 650,
                          nSidebandPars);
  SBmerge->Fit(sideband, "FSNRMQ");

  /// Prefit for the baseline
  float p_a_prefit[4];
  float p_b_prefit[4];
  auto Prefit = (TH1F*) dataHist->Clone(Form("%s_prefit", dataHist->GetName()));
  auto funct_0 = new TF1("myPol0", "pol0", 250, 750);
  Prefit->Fit(funct_0, "FSNRMQ");
  p_a_prefit[0] = funct_0->GetParameter(0);
  p_b_prefit[0] = 0;

  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 750);
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

  std::cout << "\n\nStarting the systematic variations\n\n";

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// "Systematic" variations
  unsigned int FitReg_pSigma = 0;
  unsigned int BaselineSlope = 0;
  unsigned int SourceSizeIter = 0;

//  for (unsigned int FitReg_pSigma = 0; FitReg_pSigma < 3; ++FitReg_pSigma) {
//    for (unsigned int BaselineSlope = 0; BaselineSlope < 4; ++BaselineSlope) {
//      for (unsigned int SourceSizeIter = 0; SourceSizeIter < 3;
//          ++SourceSizeIter) {
//
  if (BaselineSlope == 0) {
    useBaselineSlope = false;  //use baseline
  } else {
    useBaselineSlope = true;  // no baseline
  }

  std::cout << "\n\n_________________________\n";
  std::cout << "Processing iteration " << iterID << "\n";
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Correlation function
  /// Here we use the Ledicky model with inverse effective scattering length and effective effective range
  DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 2, NumMomBins_pSigma, kMin_pSigma,
                                  kMax_pSigma, Lednicky_Singlet_InvScatLen);

  Ck_pSigma0->SetPotPar(1, 1);
  Ck_pSigma0->SetSourcePar(0, sourceSize[SourceSizeIter]);
  Ck_pSigma0->Update();

  // TODO TO BE FIXED - MOMENTUM RESOLUTION MATRIX!
  DLM_CkDecomposition CkDec_pSigma0("pSigma0", 1, *Ck_pSigma0,
                                    CATSinput->GetSigmaFile(0));

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Sidebands
  DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins_pSigma,
                                   kMin_pSigma, kMax_pSigma, sidebandFitCATS);
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

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Fitter
  DLM_Fitter1* fitter = new DLM_Fitter1(1);
  fitter->SetSystem(0, *dataHist, 1, CkDec_pSigma0,
                    FemtoRegion_pSigma[FitReg_pSigma][0],
                    FemtoRegion_pSigma[FitReg_pSigma][1], 1000, 1000);
  fitter->SetSeparateBL(0, false);  //Simultaneous BL
  if (useBaselineSlope) {
    fitter->FixParameter("pSigma0", DLM_Fitter1::p_a,
                         p_a_prefit[BaselineSlope]);
    fitter->FixParameter("pSigma0", DLM_Fitter1::p_b,
                         p_b_prefit[BaselineSlope]);
    std::cout << "Use slope for baseline\n";
  } else {
    fitter->FixParameter("pSigma0", DLM_Fitter1::p_a, p_a_prefit[0]);
    fitter->FixParameter("pSigma0", DLM_Fitter1::p_b, p_b_prefit[0]);
  }
  fitter->AddSameSource("pSigma0SideBand", "pSigma0", 1);

  //Fit BL & Normalization
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_c, 0);
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_Cl, -1.);
  fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0,
                       sourceSize[SourceSizeIter]);

  double Val_f0Inv;
  double Val_d0;
  double lowestChi2 = 1e30;
  double GlobChi2 = 0;
  double LocChi2 = 0;
  double EffNumBins_pSigma0 = 0;
  int maxkStarBinChi2 = dataHist->FindBin(250);
  double dataY;
  double dataErr;
  double theoryX;
  double theoryY;
  TGraph FitResult_pSigma0;
  for (unsigned uBin_f0 = 0; uBin_f0 < NumBins_f0Inv; uBin_f0++) {
    Val_f0Inv = hChi2Global->GetXaxis()->GetBinCenter(uBin_f0 + 1);
    for (unsigned uBin_d0 = 0; uBin_d0 < NumBins_d0; uBin_d0++) {
      Val_d0 = hChi2Global->GetYaxis()->GetBinCenter(uBin_d0 + 1);
//      std::cout << "Processing scattering parameters 1/f0 " << Val_f0Inv
//                << " d0 " << Val_d0 << std::endl;
      hNegCk->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, 0);
      fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot0, Val_f0Inv);
      fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot1, Val_d0);
      fitter->GoBabyGo();

      // Compute local chi2
      LocChi2 = 0;
      for (unsigned uBin = 1; uBin <= maxkStarBinChi2; uBin++) {
        double mom = dataHist->GetBinCenter(uBin);
        if (mom > FemtoRegion_pSigma[FitReg_pSigma][1]) {
          continue;
        }
        FitResult_pSigma0.SetName(TString::Format("pSigma0Graph"));
        fitter->GetFitGraph(0, FitResult_pSigma0);
        FitResult_pSigma0.GetPoint(uBin - 1, theoryX, theoryY);

        if (mom != theoryX) {
          std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX << std::endl;
        }
        dataY = dataHist->GetBinContent(uBin);
        dataErr = dataHist->GetBinError(uBin);
        LocChi2 += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
      }

      GlobChi2 = fitter->GetChi2();
      if (hNegCk->GetBinContent(uBin_f0 + 1, uBin_d0 + 1) == 0
          && fitter->CheckNegativeCk()) {
        hNegCk->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, 1);
      }
      if (LocChi2 < lowestChi2) {
        lowestChi2 = LocChi2;
      }
      hChi2Global->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, GlobChi2);
      hChi2Local->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, LocChi2);
      hnSigmaPlot->SetBinContent(
          uBin_f0 + 1, uBin_d0 + 1,
          std::sqrt(2) * TMath::ErfcInverse(fitter->GetPval()));
    }
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Get the parameters from the fit
  const double bl_a = fitter->GetParameter("pSigma0", DLM_Fitter1::p_a);
  const double bl_a_err = fitter->GetParError("pSigma0", DLM_Fitter1::p_a);
  const double bl_b = fitter->GetParameter("pSigma0", DLM_Fitter1::p_b);
  const double bl_b_err = fitter->GetParError("pSigma0", DLM_Fitter1::p_b);
  const double Cl = fitter->GetParameter("pSigma0", DLM_Fitter1::p_c);
  const double chi2 = fitter->GetChi2Ndf();
  const double pval = fitter->GetPval();

  FitResult_pSigma0.SetName(TString::Format("pSigma0Graph"));
  fitter->GetFitGraph(0, FitResult_pSigma0);

  double Chi2_pSigma0 = 0;
  int maxkStarBin = dataHist->FindBin(200);
  for (unsigned uBin = 1; uBin <= maxkStarBin; uBin++) {

    double mom = dataHist->GetBinCenter(uBin);
    //double dataX;
    double dataY;
    double dataErr;
    double theoryX;
    double theoryY;

    if (mom > FemtoRegion_pSigma[FitReg_pSigma][1]) {
      continue;
    }

    FitResult_pSigma0.GetPoint(uBin - 1, theoryX, theoryY);
    if (mom != theoryX) {
      std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX << std::endl;
    }
    dataY = dataHist->GetBinContent(uBin + 1);
    dataErr = dataHist->GetBinError(uBin + 1);
    Chi2_pSigma0 += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
    ++EffNumBins_pSigma0;
  }
  double pvalpSigma0 = TMath::Prob(Chi2_pSigma0, round(EffNumBins_pSigma0));
  double nSigmapSigma0 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalpSigma0);

  std::cout << "=============\n";
  std::cout << "Fitter output\n";
  std::cout << "BL a  " << bl_a << " " << bl_a_err << "\n";
  std::cout << "BL b  " << bl_b << " " << bl_b_err << "\n";
  std::cout << "Cl    " << Cl << "\n";
  std::cout << "Chi2\n";
  std::cout << " glob " << chi2 << "\n";
  std::cout << " loc  " << Chi2_pSigma0 / round(EffNumBins_pSigma0) << "\n";
  std::cout << "p-val\n";
  std::cout << " glob " << pval << "\n";
  std::cout << " loc  " << pvalpSigma0 << "\n";
  std::cout << "=============\n";

  delete Ck_SideBand;
  delete Ck_pSigma0;
  delete fitter;

//        if (NumIter == 0) {
//          goto exitThroughTheGiftShop;
//        }
//      }
//    }
//  }
//  exitThroughTheGiftShop:

  int nBinsX = hNegCk->GetXaxis()->GetNbins();
  float xMin = hNegCk->GetXaxis()->GetBinLowEdge(1);
  float xMax = hNegCk->GetXaxis()->GetBinUpEdge(nBinsX);
  int nBinsY = hNegCk->GetYaxis()->GetNbins();
  float yMin = hNegCk->GetYaxis()->GetBinLowEdge(1);
  float yMax = hNegCk->GetYaxis()->GetBinUpEdge(nBinsY);
  auto lednickyReallySucks = new TGraphAsymmErrors();
  int iPoint = 0;
  for (int iXbin = 1; iXbin <= nBinsX; ++iXbin) {
    float yBinMin = 0;
    float yBinMax = yMax;
    float binWidth = 0;
    for (int iYbin = 1; iYbin <= nBinsY; ++iYbin) {
      binWidth = hNegCk->GetYaxis()->GetBinWidth(iYbin);
      if (hNegCk->GetBinContent(iXbin, iYbin) != 0) {
        yBinMin = hNegCk->GetYaxis()->GetBinLowEdge(iYbin);
        break;
      }
    }
    if (yBinMin == 0) {
      continue;
    }
    float err = yBinMax - yBinMin;
    float xValue = hNegCk->GetXaxis()->GetBinLowEdge(iXbin);
    lednickyReallySucks->SetPoint(iPoint, xValue, yBinMin);
    lednickyReallySucks->SetPointError(iPoint++, 0, 0, 0, 20);
  }

  for (int iXbin = 1; iXbin <= hChi2Local->GetXaxis()->GetNbins(); ++iXbin) {
    for (int iYbin = 1; iYbin <= hChi2Local->GetYaxis()->GetNbins(); ++iYbin) {
      hDeltaChi2Plot->SetBinContent(
          iXbin, iYbin, hChi2Local->GetBinContent(iXbin, iYbin) - lowestChi2);
    }
  }

  exclusionFile->cd();
  hChi2Local->Write("chi2local");
  hChi2Global->Write("chi2global");
  hnSigmaPlot->Write("nSigma");
  hNegCk->Write("negCk");
  lednickyReallySucks->Write("negCkGraph");

  double contours[3];
  contours[0] = 2.3;
  contours[1] = 6.18;
  contours[2] = 11.83;
  hDeltaChi2Plot->SetContour(3, contours);
  hDeltaChi2Plot->Write("DeltaChi2");

  const int NRGBs = 5;
  Double_t stops[NRGBs];
  for (int i = 0; i < NRGBs; ++i)
    stops[i] = float(i) / (NRGBs - 1);
  Double_t red[NRGBs] = { 1., 29. / 255., 25. / 255., 27. / 255., 32. / 255. };
  Double_t green[NRGBs] = { 1., 221. / 255., 160. / 255., 113. / 255., 74.
      / 255. };
  Double_t blue[NRGBs] = { 1., 221. / 255., 184. / 255., 154. / 255., 129.
      / 255. };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);

  auto c = new TCanvas();
  c->SetRightMargin(0.16);
  hDeltaChi2Plot->SetMaximum(25);
  hDeltaChi2Plot->Draw("colz");

  lednickyReallySucks->SetFillStyle(3004);
//  lednickyReallySucks->SetLineStyle(3004);
  lednickyReallySucks->SetLineWidth(0);
  lednickyReallySucks->SetFillColor(kBlack);
  lednickyReallySucks->SetLineColorAlpha(kBlack, 0.7);
  lednickyReallySucks->Draw("F SAME");

  c->Print(Form("%s/pSigma0_exclusion.pdf", OutputDir.Data()));
  c->Print(Form("%s/pSigma0_exclusion.png", OutputDir.Data()));
  c->Write("exclusion");
  exclusionFile->Close();

  return;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  FitSigma0exclusion(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5],
                     atof(argv[6]));
  return 0;
}

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
#include "SidebandSigma.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"

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
  bool useBaseline = true;
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  TString exclusionfilename = TString::Format("%s/Exclusion_pSigma0_%i.root",
                                              OutputDir.Data(), NumIter);
  auto exclusionFile = new TFile(exclusionfilename, "RECREATE");

  int NumBins_f0Inv = 1000;
  double Min_f0Inv = -5;
  double Max_f0Inv = 5;
  int NumBins_d0 = 1000;
  double Min_d0 = 0.;
  double Max_d0 = 10.;

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
  const int binwidth = dataHist->GetBinWidth(1);
  const unsigned NumMomBins_pSigma = int(800 / binwidth);
  double kMin_pSigma = dataHist->GetBinCenter(1) - binwidth / 2.f;
  const double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;

  std::cout << "kMin_pSigma: " << kMin_pSigma << std::endl;
  std::cout << "kMax_pSigma: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pSigma: " << NumMomBins_pSigma << std::endl;

  // pp radius systematic variations
  const double ppRadius = 1.36;
  const std::vector<double> sourceSize = { { ppRadius, ppRadius * 0.8, ppRadius
      * 1.2 } };

  // femto fit region systematic variations
  const std::vector<double> femtoFitRegionUp = { { 550, 500, 600 } };

  /// Lambda parameters
  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;

  // proton secondary contribution systematic variations
  const std::vector<double> protonSecondary = { { protonLambda
      / (1. - protonPrimary), protonLambda / (1. - protonPrimary) * 0.8,
      protonLambda / (1. - protonPrimary) * 1.2 } };
  std::vector<CATSLambdaParam> lambdaParams;

  const double sigmaPurity = 0.199;
  const double sigmaPrimary = 1.;
  const Particle sigma0(sigmaPurity, sigmaPrimary, { { 0 } });

  for (size_t lambdaIter = 0; lambdaIter < protonSecondary.size();
      ++lambdaIter) {

    std::cout << "\n\n_________________________\n";
    std::cout << "Processing iteration " << iterID << "\n\n";

    const Particle proton(
        protonPurity,
        protonPrimary,
        { { (1. - protonPrimary) * protonSecondary[lambdaIter], (1.
            - protonPrimary) * (1 - protonSecondary[lambdaIter]) } });

    lambdaParams.push_back( { proton, sigma0 });
    lambdaParams[lambdaIter].PrintLambdaParams();

  }

  // sideband fit normalization range systematic variation
  const std::vector<double> sidebandNormDown = { { 340, 400, 460 } };
  const std::vector<double> sidebandNormUp = { { 440, 500, 560 } };

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Fit the sideband
  auto side = new SidebandSigma();
  side->SetRebin(10);
  side->SetSideBandFile(InputDir.Data(), appendix.Data());

  /// Prefit for the baseline
  auto Prefit = (TH1F*) dataHist->Clone(Form("%s_prefit", dataHist->GetName()));
  auto funct_0 = new TF1("myPol0", "pol0", 250, 750);
  Prefit->Fit(funct_0, "FSNRMQ");

  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 750);
  Prefit->Fit(funct_1, "FSNRMQ");
  gMinuit->SetErrorDef(1);  // 1 corresponds to 1 sigma contour, 4 to 2 sigma
  auto grPrefitContour = (TGraph*) gMinuit->Contour(40, 0, 1);

  std::vector<double> prefit_a;
  prefit_a.emplace_back(funct_0->GetParameter(0));
  prefit_a.emplace_back(funct_1->GetParameter(0));
  prefit_a.emplace_back(
      TMath::MinElement(grPrefitContour->GetN(), grPrefitContour->GetX()));
  prefit_a.emplace_back(
      TMath::MaxElement(grPrefitContour->GetN(), grPrefitContour->GetX()));

  std::vector<double> prefit_b;
  prefit_b.emplace_back(0);
  prefit_b.emplace_back(funct_1->GetParameter(1));
  prefit_b.emplace_back(grPrefitContour->Eval(prefit_a[2]));
  prefit_b.emplace_back(grPrefitContour->Eval(prefit_a[3]));

  std::cout << "Result of the prefit to constrain the baseline\n";
  for (size_t i = 0; i < prefit_a.size(); ++i) {
    std::cout << i << " a: " << prefit_a[i] << " b: " << prefit_b[i] << "\n";
  }

  delete funct_0;
  delete funct_1;
  delete Prefit;

  if (NumIter != 0) {
    std::cout << "\n\nStarting the systematic variations\n";
    std::cout << "Number of variations of the fit region: "
              << femtoFitRegionUp.size() << "\n";
    std::cout << "Number of variations of the baseline:   " << prefit_a.size()
              << "\n";
    std::cout << "Number of variations of the source size : "
              << sourceSize.size() << "\n";
    std::cout << "Number of variations of the sideband normalization: "
              << sidebandNormDown.size() << "\n";
    std::cout << "Number of variations of the lambda param: "
              << protonSecondary.size() << "\n\n";
  }

  // Set up the model, fitter, etc.
  DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 2, NumMomBins_pSigma, kMin_pSigma,
                                  kMax_pSigma, Lednicky_Singlet_InvScatLen);
  Ck_pSigma0->SetPotPar(1, 1);

  double highestChi2 = -1.;
  double Val_f0Inv;
  double Val_d0;
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
      highestChi2 = -1.;
      std::cout << "Processing scattering parameters 1/f0 " << Val_f0Inv
                << " d0 " << Val_d0 << std::endl;
      hNegCk->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, 0);

      /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      /// "Systematic" variations

      // 1. Femto fit range
      for (size_t femtoFitIter = 0; femtoFitIter < femtoFitRegionUp.size();
          ++femtoFitIter) {

        // 2. Baseline treatment
        for (size_t blIter = 0; blIter < prefit_a.size(); ++blIter) {

          if (blIter == 0) {
            useBaseline = false;  //use baseline
          } else {
            useBaseline = true;  // no baseline
          }

          // 3. Sideband normalization
          for (size_t sbNormIter = 0; sbNormIter < sidebandNormDown.size();
              ++sbNormIter) {

            side->SetNormalizationRange(sidebandNormDown[sbNormIter],
                                        sidebandNormUp[sbNormIter]);

            side->SideBandCFs();
            auto SBmerge = side->GetSideBands(5);
            const int firstBin = SBmerge->GetBinCenter(1);
            auto sideband = new TF1("sideband", sidebandFit, firstBin, 650,
                                    nSidebandPars);
            sideband->SetParameter(0, -1.5);
            sideband->SetParameter(1, 0.);
            sideband->SetParameter(2, 0);
            sideband->SetParameter(3, 0);
            sideband->SetParameter(4, 1);
            sideband->SetParameter(5, 0);
            SBmerge->Fit(sideband, "FSNRMQ");

            DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars,
                                             NumMomBins_pSigma, kMin_pSigma,
                                             kMax_pSigma, sidebandFitCATS);
            for (unsigned i = 0; i < sideband->GetNumberFreeParameters(); ++i) {
              Ck_SideBand->SetPotPar(i, sideband->GetParameter(i));
            }
            Ck_SideBand->Update();

            // 4. Source size
            for (size_t sizeIter = 0; sizeIter < sourceSize.size();
                ++sizeIter) {

              // 5. Lambda parameters
              for (size_t lambdaIter = 0; lambdaIter < lambdaParams.size();
                  ++lambdaIter) {

                /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                /// Correlation function
                Ck_pSigma0->SetSourcePar(0, sourceSize[sizeIter]);
                Ck_pSigma0->Update();

                // TODO TO BE FIXED - MOMENTUM RESOLUTION MATRIX!
                DLM_CkDecomposition CkDec_pSigma0("pSigma0", 1, *Ck_pSigma0,
                                                  CATSinput->GetSigmaFile(0));

                /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                DLM_CkDecomposition CkDec_SideBand("pSigma0SideBand", 0,
                                                   *Ck_SideBand, nullptr);
                CkDec_pSigma0.AddContribution(
                    0,
                    lambdaParams[lambdaIter].GetLambdaParam(
                        CATSLambdaParam::Fake),
                    DLM_CkDecomposition::cFake, &CkDec_SideBand);
                CkDec_pSigma0.Update();

                /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                /// Fitter
                DLM_Fitter1* fitter = new DLM_Fitter1(1);
                fitter->SetSystem(0, *dataHist, 1, CkDec_pSigma0, kMin_pSigma,
                                  femtoFitRegionUp[femtoFitIter], 1000, 1000);
                fitter->SetSeparateBL(0, false);              //Simultaneous BL
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_a,
                                     prefit_a[blIter]);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_b,
                                     prefit_b[blIter]);
                fitter->AddSameSource("pSigma0SideBand", "pSigma0", 1);

                //Fit BL & Normalization
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_c, 0);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_Cl, -1.);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0,
                                     sourceSize[sizeIter]);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot0, Val_f0Inv);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot1, Val_d0);

                fitter->GoBabyGo();
                // Compute local chi2
                LocChi2 = 0;
                for (unsigned uBin = 1; uBin <= maxkStarBinChi2; uBin++) {
                  double mom = dataHist->GetBinCenter(uBin);
                  if (mom > femtoFitRegionUp[femtoFitIter]) {
                    continue;
                  }
                  FitResult_pSigma0.SetName(TString::Format("pSigma0Graph"));
                  fitter->GetFitGraph(0, FitResult_pSigma0);
                  FitResult_pSigma0.GetPoint(uBin - 1, theoryX, theoryY);

                  if (mom != theoryX) {
                    std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX
                              << std::endl;
                  }
                  dataY = dataHist->GetBinContent(uBin);
                  dataErr = dataHist->GetBinError(uBin);
                  LocChi2 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                }

                GlobChi2 = fitter->GetChi2();
                if (hNegCk->GetBinContent(uBin_f0 + 1, uBin_d0 + 1) == 0
                    && fitter->CheckNegativeCk()) {
                  hNegCk->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, 1);
                }
                if (LocChi2 > highestChi2) {
                  highestChi2 = LocChi2;
                }
                hChi2Global->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, GlobChi2);
                hnSigmaPlot->SetBinContent(
                    uBin_f0 + 1, uBin_d0 + 1,
                    std::sqrt(2) * TMath::ErfcInverse(fitter->GetPval()));

                delete fitter;

                if (NumIter == 0) {
                  std::cout << "Skipping all systematic variations \n";
                  goto exitThroughTheGiftShop;
                }
              }
            }
            delete sideband;
            delete Ck_SideBand;
          }
        }
      }
      exitThroughTheGiftShop:

      hChi2Local->SetBinContent(uBin_f0 + 1, uBin_d0 + 1, LocChi2);
    }
  }

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
          iXbin, iYbin, hChi2Local->GetBinContent(iXbin, iYbin) - highestChi2);
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
  lednickyReallySucks->Draw("3 SAME");

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
